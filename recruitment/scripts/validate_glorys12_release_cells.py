"""
Acceptance tests of the release-cell extraction
(preprocess_glorys12_release_cells.py).

Usage
-----
    # 1. Coverage: is the grid point nearest to every release position among
    #    the extracted cells? Needs no extraction file.
    python3 validate_glorys12_release_cells.py --releases ../data

    # 2. Day 0: does the extraction reproduce the sea ice and temperature a
    #    trajectory recorded at its release? Needs the extraction file of
    #    the cohorts' spawning year.
    python3 validate_glorys12_release_cells.py \\
        $KRICO_RUNS/KRICO_0001/1993_11_15.nc $KRICO_RUNS/KRICO_0003/1994_01_15.nc

Coverage
--------
--releases DIR reads release_lon and release_lat from every classification
file YYYY_MM_DD.nc in DIR (the v2.0.0 files hold them for every particle,
exactly as the trajectories store them), finds the nearest grid point by the
rule of M1 (sea_ice._sample_nearest) and reports how many particles fall in
a release cell, in a neighbour, or outside the extracted cells. Any particle
outside means the cell set is too small: the extraction must not be run.

Day 0
-----
Parcels samples a tracer on the C-grid without spatial interpolation: a
particle reads the value of the cell that holds it, bounded by the f-points
(glamf, gphif of glorys12_coordinates.nc), at the level whose bounds
(depthw) hold its depth. Particles are released at 12:00, the time stamp of
the daily means, so the values recorded at day 0 are those of the release
day. The extraction must therefore return them exactly.

For each trajectory file the test compares the trajectory's
sea_ice_area_fraction and temperature at obs 0 with the extraction

  - at the cell holding the release position (the Parcels rule) and at the
    grid point nearest to it (the rule of the classification); the two
    differ for the particles described in krico_recruitment/release_cells.py;
  - on the release day and, as controls of the indexing, one day earlier and
    later, and one level higher and lower. The controls must be worse.

Verdict per file: PASS when no particle is outside the extracted cells, the
unshifted comparison is the best one, and at least 99.9% of the particles
are reproduced exactly at the cell holding them.

Environment
-----------
    KRICO_GLORYS12                 required; grid, level bounds, f-points
    KRICO_GLORYS12_RELEASE_CELLS   required for the day-0 test
"""

import argparse
import os
import sys
import time
from datetime import date
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from krico_recruitment import release_cells  # noqa: E402

EXACT_SHARE_REQUIRED = 0.999
TOLERANCE = 0.01          # sea-ice fraction or degrees Celsius, for the counts


def _env(name):
    v = os.environ.get(name)
    if not v:
        sys.exit(f"ERROR: {name} environment variable not set.")
    return Path(v)


# ----------------------------------------------------------------- coverage

def coverage(release_dir, grid, cells):
    files = sorted(Path(release_dir).glob("[0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9].nc"))
    if not files:
        sys.exit(f"ERROR: no YYYY_MM_DD.nc file in {release_dir}")
    shape = grid["bathymetry"].shape
    ncell = len(cells["cell_y"])
    lookup = release_cells.cell_lookup(cells["cell_y"], cells["cell_x"], shape)
    per_cell = np.zeros(ncell, dtype=np.int64)
    n_total = n_outside = 0
    worst = []
    started = time.time()
    print(f"Coverage of {len(files)} files in {release_dir}", flush=True)
    for n, path in enumerate(files, 1):
        with netCDF4.Dataset(path) as nc:
            lon = np.asarray(nc["release_lon"][:], dtype=np.float64)
            lat = np.asarray(nc["release_lat"][:], dtype=np.float64)
        idx = release_cells.locate_cells(lon, lat, grid["grid_lon"],
                                         grid["grid_lat"], lookup)
        inside = idx >= 0
        per_cell += np.bincount(idx[inside], minlength=ncell)
        n_total += idx.size
        if not inside.all():
            n_outside += int((~inside).sum())
            if len(worst) < 5:
                k = int(np.nonzero(~inside)[0][0])
                worst.append(f"{path.name}: lon {lon[k]:.5f}, lat {lat[k]:.5f}")
        if n % 500 == 0:
            print(f"  {n} files, {time.time() - started:.0f} s", flush=True)

    # Where the neighbours that received particles lie, relative to the
    # release cells: the cell south of a release cell has one to its north.
    release = np.zeros(shape, dtype=bool)
    release[cells["cell_y"][cells["is_release"] == 1],
            cells["cell_x"][cells["is_release"] == 1]] = True
    cy, cx = cells["cell_y"], cells["cell_x"]
    is_release = cells["is_release"] == 1
    used = per_cell > 0
    padded = np.pad(release, 1)                  # no release cell off the grid
    north = padded[cy + 2, cx + 1]
    east = padded[cy + 1, cx + 2]
    north_east = padded[cy + 2, cx + 2]
    classes = (
        ("south of a release cell", ~is_release & north),
        ("west of a release cell", ~is_release & ~north & east),
        ("south-west of a release cell", ~is_release & ~north & ~east & north_east),
        ("other neighbours", ~is_release & ~north & ~east & ~north_east),
    )
    print(f"  particles: {n_total} in {len(files)} files, "
          f"{time.time() - started:.0f} s")
    print(f"  in a release cell:  {int(per_cell[is_release].sum()):>13d} particles, "
          f"{int((used & is_release).sum())} of {int(is_release.sum())} cells used")
    for label, sel in classes:
        print(f"  {label + ':':<29s}{int(per_cell[sel].sum()):>4d} particles, "
              f"{int((used & sel).sum())} of {int(sel.sum())} cells used")
    print(f"  outside the extracted cells: {n_outside} particles")
    for line in worst:
        print(f"    {line}")
    ok = n_outside == 0
    print("  verdict: " + ("PASS, every release position has its nearest grid "
                           "point among the extracted cells" if ok else
                           "FAIL, the cell set is too small; do not run the "
                           "extraction"))
    return ok


# -------------------------------------------------------------------- day 0

def _f_axes(glorys_dir):
    """Ascending 1-D axes of the f-points (cell corners), or None if absent."""
    path = Path(glorys_dir) / "glorys12_coordinates.nc"
    if not path.is_file():
        return None
    with xr.open_dataset(path, decode_times=False) as ds:
        glamf = np.asarray(ds["glamf"].values, dtype=np.float64)
        gphif = np.asarray(ds["gphif"].values, dtype=np.float64)
    glamf = glamf.reshape(glamf.shape[-2:])
    gphif = gphif.reshape(gphif.shape[-2:])
    if not (np.allclose(glamf, glamf[0:1, :]) and np.allclose(gphif, gphif[:, 0:1])):
        raise ValueError("glamf/gphif are not separable axes over this domain")
    return glamf[0, :], gphif[:, 0]


def _holding_cells(lon, lat, f_lon, f_lat, lookup):
    """
    Cell index of the T-cell that holds each position: T-cell (j, i) lies
    between f-points i-1 and i in longitude, j-1 and j in latitude. -1 where
    that cell is not extracted or the position is off the grid.
    """
    i = np.searchsorted(f_lon, lon, side="right")
    j = np.searchsorted(f_lat, lat, side="right")
    ny, nx = lookup.shape
    ok = (i >= 1) & (i < nx) & (j >= 1) & (j < ny)
    out = np.full(lon.shape, -1, dtype=np.int32)
    out[ok] = lookup[j[ok], i[ok]]
    return out


def _stats(traj, ref):
    """Trajectory against extraction where both hold a value."""
    both = np.isfinite(traj) & np.isfinite(ref)
    n = int(both.sum())
    if n == 0:
        return {"n": 0, "exact": np.nan, "mean": np.nan, "max": np.nan, "above": 0}
    d = np.abs(traj[both].astype(np.float64) - ref[both].astype(np.float64))
    return {"n": n, "exact": float((d == 0).mean()), "mean": float(d.mean()),
            "max": float(d.max()), "above": int((d > TOLERANCE).sum())}


def _row(label, s, n_particles):
    lost = n_particles - s["n"]
    return (f"    {label:<34s}{100 * s['exact']:9.4f}%  {s['mean']:10.3e}  "
            f"{s['max']:10.3e}  {s['above']:9d}  {lost:9d}")


def _read(var):
    values = np.asarray(var, dtype=np.float32)
    values[~(np.abs(values) < 1e30)] = np.nan
    return values


def day0(path, extraction_dir, f_axes):
    started = time.time()
    with xr.open_dataset(path, decode_timedelta=False) as ds:
        first = ds.isel(obs=0)
        lon = np.asarray(first["lon"].values, dtype=np.float64)
        lat = np.asarray(first["lat"].values, dtype=np.float64)
        depth = np.asarray(first["z"].values, dtype=np.float64)
        traj_temp = np.asarray(first["temperature"].values, dtype=np.float32)
        traj_ice = np.asarray(first["sea_ice_area_fraction"].values, dtype=np.float32)
        stamps = pd.DatetimeIndex(np.asarray(first["time"].values).ravel())
    n = lon.size
    release_times = stamps.unique()
    release_day = release_times[0].date()
    print(f"{Path(path).name}: {n} particles, released "
          + ", ".join(str(t) for t in release_times[:3]))
    if len(release_times) != 1:
        print("  verdict: FAIL, more than one release time in the file")
        return False
    year = release_cells.spawning_year(release_day)
    extraction = Path(extraction_dir) / release_cells.FILE_PATTERN.format(year=year)
    if not extraction.is_file():
        sys.exit(f"ERROR: {extraction} is missing; run the extraction of "
                 f"spawning year {year} first")

    with netCDF4.Dataset(extraction) as nc:
        first_day = date.fromisoformat(nc.first_day)
        t = (release_day - first_day).days
        if not 1 <= t < nc.dimensions["time"].size - 1:
            sys.exit(f"ERROR: {release_day} is not inside {extraction.name}")
        stamp = netCDF4.num2date(nc["time"][t], nc["time"].units,
                                 calendar=nc["time"].calendar,
                                 only_use_cftime_datetimes=False,
                                 only_use_python_datetimes=True)
        if stamp.date() != release_day:
            sys.exit(f"ERROR: {extraction.name}: time index {t} is {stamp}, "
                     f"expected {release_day}")
        cell_y = np.asarray(nc["cell_y"][:])
        cell_x = np.asarray(nc["cell_x"][:])
        is_release = np.asarray(nc["is_release_cell"][:]) == 1
        grid_lon = np.asarray(nc["grid_lon"][:], dtype=np.float64)
        grid_lat = np.asarray(nc["grid_lat"][:], dtype=np.float64)
        bounds = np.asarray(nc["depth_bounds"][:], dtype=np.float64)
        for name in ("sea_ice_area_fraction", "temperature"):
            nc[name].set_auto_maskandscale(False)
        ice = _read(nc["sea_ice_area_fraction"][t - 1:t + 2, :])      # (3, cell)
        temp = _read(nc["temperature"][t - 1:t + 2, :, :])            # (3, depth, cell)
    nlev = temp.shape[1]
    lookup = release_cells.cell_lookup(cell_y, cell_x, (grid_lat.size, grid_lon.size))

    nearest = release_cells.locate_cells(lon, lat, grid_lon, grid_lat, lookup)
    outside = int((nearest < 0).sum())
    in_release = int(is_release[nearest[nearest >= 0]].sum())
    print(f"  nearest grid point: {in_release} particles in a release cell, "
          f"{n - outside - in_release} in a neighbour, {outside} outside the "
          "extracted cells")
    rules = [("nearest", nearest)]
    if f_axes is not None:
        f_lon, f_lat = f_axes
        if (f_lon.size != grid_lon.size or f_lat.size != grid_lat.size
                or not (np.all(f_lon[:-1] < grid_lon[1:]) and np.all(grid_lon[1:] < f_lon[1:])
                        and np.all(f_lat[:-1] < grid_lat[1:]) and np.all(grid_lat[1:] < f_lat[1:]))):
            sys.exit("ERROR: the f-points of glorys12_coordinates.nc do not "
                     "bracket the cell centres as expected (T-point i between "
                     "f-points i-1 and i)")
        holding = _holding_cells(lon, lat, f_axes[0], f_axes[1], lookup)
        print(f"  cell holding the position: {int((holding < 0).sum())} outside "
              f"the extracted cells; differs from the nearest grid point for "
              f"{int((holding != nearest).sum())} particles")
        rules.insert(0, ("holding", holding))
    else:
        print("  glorys12_coordinates.nc not found: nearest grid point only")
    level = release_cells.locate_levels(depth, bounds)
    print(f"  release depths {np.nanmin(depth):.1f} to {np.nanmax(depth):.1f} m, "
          f"levels {int(level.min()) + 1} to {int(level.max()) + 1}")
    if level.min() < 1 or level.max() > nlev - 2:
        sys.exit("ERROR: release depths outside the levels that can be tested")

    def sample_ice(cell, dt):
        ref = np.full(n, np.nan, dtype=np.float32)
        ok = cell >= 0
        ref[ok] = ice[1 + dt, cell[ok]]
        return ref

    def sample_temp(cell, dt, dk):
        ref = np.full(n, np.nan, dtype=np.float32)
        ok = cell >= 0
        ref[ok] = temp[1 + dt, level[ok] + dk, cell[ok]]
        return ref

    header = (f"    {'':<34s}{'exact':>10s}  {'mean |d|':>10s}  {'max |d|':>10s}  "
              f"{'|d|>' + str(TOLERANCE):>9s}  {'no value':>9s}")
    main_rule, main_cell = rules[0]
    results = {}
    print("  sea ice at day 0, trajectory minus extraction")
    print(header)
    for label, cell, dt in ([(f"{main_rule} cell, release day", main_cell, 0)]
                            + [(f"{r} cell, release day", c, 0) for r, c in rules[1:]]
                            + [(f"{main_rule} cell, one day earlier", main_cell, -1),
                               (f"{main_rule} cell, one day later", main_cell, 1)]):
        results["ice", label] = s = _stats(traj_ice, sample_ice(cell, dt))
        print(_row(label, s, n))
    print("  temperature at day 0, trajectory minus extraction")
    print(header)
    for label, cell, dt, dk in (
            [(f"{main_rule} cell, release day", main_cell, 0, 0)]
            + [(f"{r} cell, release day", c, 0, 0) for r, c in rules[1:]]
            + [(f"{main_rule} cell, one day earlier", main_cell, -1, 0),
               (f"{main_rule} cell, one day later", main_cell, 1, 0),
               (f"{main_rule} cell, one level higher", main_cell, 0, -1),
               (f"{main_rule} cell, one level lower", main_cell, 0, 1)]):
        results["temp", label] = s = _stats(traj_temp, sample_temp(cell, dt, dk))
        print(_row(label, s, n))

    key = f"{main_rule} cell, release day"
    reasons = []
    if outside:
        reasons.append(f"{outside} particles outside the extracted cells")
    for var, name in (("ice", "sea ice"), ("temp", "temperature")):
        best = results[var, key]
        if best["n"] < n:
            reasons.append(f"{name}: no value for {n - best['n']} particles")
        shifted = [s["mean"] for (v, label), s in results.items()
                   if v == var and label.startswith(main_rule) and label != key]
        if not all(best["mean"] < m for m in shifted):
            reasons.append(f"{name}: a shifted comparison is as good as the "
                           "unshifted one (indexing)")
        if not best["exact"] >= EXACT_SHARE_REQUIRED:
            reasons.append(f"{name}: {100 * best['exact']:.4f}% reproduced "
                           f"exactly, below {100 * EXACT_SHARE_REQUIRED}%")
    ok = not reasons
    print("  verdict: " + ("PASS" if ok else "FAIL, " + "; ".join(reasons))
          + f"  ({time.time() - started:.0f} s)")
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("trajectories", nargs="*", type=Path,
                        help="trajectory files for the day-0 test")
    parser.add_argument("--releases", type=Path, metavar="DIR",
                        help="folder of classification files for the coverage test")
    parser.add_argument("--cells", type=Path, default=release_cells.TABLE,
                        help="release-cell table (default: the committed one)")
    args = parser.parse_args()
    if not args.trajectories and args.releases is None:
        parser.error("give --releases DIR, trajectory files, or both")
    glorys_dir = _env("KRICO_GLORYS12")
    ok = True
    try:
        if args.releases is not None:
            grid = release_cells.read_grid(glorys_dir)
            cells = release_cells.build_cells(grid, release_cells.read_table(args.cells))
            ok &= coverage(args.releases, grid, cells)
        if args.trajectories:
            extraction_dir = _env("KRICO_GLORYS12_RELEASE_CELLS")
            f_axes = _f_axes(glorys_dir)
            for path in args.trajectories:
                ok &= day0(path, extraction_dir, f_axes)
    except (ValueError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: {exc}")
    print("ALL PASS" if ok else "AT LEAST ONE FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
