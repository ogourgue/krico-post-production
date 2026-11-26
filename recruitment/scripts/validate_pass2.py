"""
Tests of pass 2 (process_cohort_chl.py) on files it has written.

Usage
-----
    python3 validate_pass2.py ../data_v3.0_chl_test/1994_01_15.nc
    python3 validate_pass2.py FILE... --pass1 ../data_v3.0_glorys12 --sample 300

For each file:

  1. Layout: the variables are those of the schema.
  2. Consistency, on every particle: no window longer than its bounds
     allow; the days with a value never exceed the days of the window, a
     minimum count of 10 never keeps more days than a minimum of 1, the
     matched days never exceed either product's; a mean exists where, and
     only where, its number of days is not zero; every mean is positive
     and already rounded.
  3. Windows against pass 1: calyptopis_end_day and advance_day equal
     those of the pass 1 file of the cohort for every particle, and the
     particles are the same, in the same order. The pass 1 file is looked
     for in --pass1 (default: recruitment/data_v3.0_glorys12); the test is
     skipped, and said to be, if it is not there.
  4. Recomputation: for a sample of particles, every index is computed
     again here from the trajectory file and the two products, by another
     route than the package: the files are opened and decoded by xarray,
     the cell is found by nearest-coordinate selection, the level by the
     nearest depth, and the sums are plain loops. Stored and recomputed
     values must be equal after rounding, and the numbers of days equal.
     The sample holds --sample random particles (default 300) and
     particles picked for the cases that matter (no advance, empty
     furcilia window, deleted early, model without a value at depth,
     satellite without any value). A particle that lies within EDGE
     degree of a cell edge of either grid on a day of its windows is left
     out and counted: the coordinates of the files are float32, so which
     of two cells is the nearest is then a matter of rounding.
  5. Report: per window, the share of particles with a value and the
     typical values; how the two satellite statistics and the two products
     compare; the bytes each variable takes on disk (needs h5py).

Environment
-----------
    KRICO_RUNS, KRICO_GLOBCOLOUR, KRICO_PISCES   required for step 4
"""

import argparse
import os
import sys
from datetime import date, timedelta
from pathlib import Path

import netCDF4
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from krico_recruitment import chlorophyll, classify, inputs, io as krico_io  # noqa: E402

RECRUITMENT = Path(__file__).resolve().parents[1]
DEFAULT_PASS1 = RECRUITMENT / "data_v3.0_glorys12"
NONE = classify.DAY_NONE
# Distance (degree) to a cell edge under which the cell holding a position
# depends on the rounding of the float32 coordinates; about half a metre.
EDGE = 5e-6
PAIRS = (("chl_satellite_mean", "chl_satellite_days"),
         ("chl_satellite_median", "chl_satellite_days"),
         ("chl_model_surface", "chl_model_days"),
         ("chl_model_mean_200m", "chl_model_days"),
         ("chl_model_at_depth", "chl_model_at_depth_days"),
         ("chl_model_surface_matched", "chl_matched_days"))


def _rounded(x):
    """Positive values rounded to the stored precision, by another route."""
    x = np.asarray(x, dtype=np.float32).astype(np.float64)
    scale = 2.0 ** chlorophyll.SIGNIFICANT_BITS
    mantissa, exponent = np.frexp(x)
    return np.ldexp(np.floor(mantissa * scale + 0.5) / scale, exponent)


def _same(stored, reference):
    """Stored value against a recomputed one, both possibly missing."""
    if np.isnan(reference):
        return bool(np.isnan(stored))
    if np.isnan(stored):
        return False
    return any(float(stored) == float(_rounded(reference * f))
               for f in (1.0, 1.0 + 1e-6, 1.0 - 1e-6))


def _read(path):
    d = {}
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_maskandscale(False)
        for name, var in nc.variables.items():
            d[name] = np.asarray(var[:])
        attrs = {k: nc.getncattr(k) for k in nc.ncattrs()}
    return d, attrs


def _storage(path):
    try:
        import h5py
    except ImportError:
        return None
    with h5py.File(path, "r") as f:
        return {name: f[name].id.get_storage_size() for name in f.keys()}


# -------------------------------------------------------------- consistency

def consistency(d, n_obs):
    reasons = []
    end = np.where(d["calyptopis_end_day"] == NONE, n_obs,
                   d["calyptopis_end_day"]).astype(np.int64)
    advance = np.where(d["advance_day"] == NONE, n_obs,
                       d["advance_day"]).astype(np.int64)
    longest = np.stack([end, np.maximum(advance - end, 0)], axis=1)
    window = d["window_days"].astype(np.int64)
    bad = int((window > longest).sum())
    print(f"  window_days above what the bounds allow: {bad}")
    if bad:
        reasons.append(f"{bad} windows longer than their bounds")
    tests = (
        ("satellite days above window days",
         d["chl_satellite_days"] > d["window_days"][:, :, None]),
        ("more satellite days with the larger minimum count",
         np.diff(d["chl_satellite_days"].astype(np.int64), axis=2) > 0),
        ("model days above window days", d["chl_model_days"] > d["window_days"]),
        ("model days at depth above model days",
         d["chl_model_at_depth_days"] > d["chl_model_days"]),
        ("matched days above satellite days",
         d["chl_matched_days"] > d["chl_satellite_days"]),
        ("matched days above model days",
         d["chl_matched_days"] > d["chl_model_days"][:, :, None]),
    )
    for label, wrong in tests:
        bad = int(wrong.sum())
        print(f"  {label}: {bad}")
        if bad:
            reasons.append(f"{label}: {bad}")
    for value, days in PAIRS:
        v, c = d[value], d[days]
        mismatch = int((np.isfinite(v) != (c > 0)).sum())
        with np.errstate(invalid="ignore"):
            out = int((v <= 0).sum() + (v > 1000).sum())
            unrounded = int((np.isfinite(v) & (_rounded(v) != v)).sum())
        print(f"  {value}: {mismatch} values without days or days without a "
              f"value, {out} not positive or above 1000, {unrounded} not rounded")
        if mismatch or out or unrounded:
            reasons.append(f"{value}: {mismatch + out + unrounded} inconsistent values")
    return reasons


# ---------------------------------------------------- windows against pass 1

def against_pass1(d, path, pass1_dir):
    other = Path(pass1_dir) / path.name
    if not other.is_file():
        print(f"  pass 1 file {other} not found: windows NOT compared with pass 1")
        return []
    with netCDF4.Dataset(other) as nc:
        nc.set_auto_maskandscale(False)
        ice = getattr(nc, "sic_source", "?")
        version = getattr(nc, "code_version", "?")
        p1 = {name: np.asarray(nc[name][:])
              for name in ("particle", "calyptopis_end_day", "advance_day")}
    print(f"  pass 1 file: {other} ({version}, sea ice {ice})")
    reasons = []
    if p1["particle"].shape != d["particle"].shape \
            or not np.array_equal(p1["particle"], d["particle"]):
        print("  particles: NOT the same as in the pass 1 file")
        return ["particles differ from the pass 1 file"]
    for name in ("calyptopis_end_day", "advance_day"):
        bad = int((p1[name] != d[name]).sum())
        print(f"  {name} against pass 1: {bad} of {d[name].size} differ")
        if bad:
            reasons.append(f"{name} differs from pass 1 for {bad} particles")
    return reasons


# ------------------------------------------------------------- recomputation

def _nearest(field, lon, lat):
    """
    Values of an xarray field at positions, by nearest-coordinate
    selection; NaN where the position is more than half a cell outside.
    Also returns which positions lie within EDGE of a cell edge.
    """
    step_lat = float(field.latitude[1] - field.latitude[0])
    step_lon = float(field.longitude[1] - field.longitude[0])
    picked = field.sel(latitude=xr.DataArray(lat, dims="p"),
                       longitude=xr.DataArray(lon, dims="p"), method="nearest")
    off_lat = np.abs(picked.latitude.values.astype(np.float64) - lat)
    off_lon = np.abs(picked.longitude.values.astype(np.float64) - lon)
    outside = (off_lat > 0.5 * step_lat) | (off_lon > 0.5 * step_lon)
    edge = ((np.abs(off_lat - 0.5 * step_lat) < EDGE)
            | (np.abs(off_lon - 0.5 * step_lon) < EDGE))
    values = np.array(picked.values, dtype=np.float64)
    values[..., outside] = np.nan
    return values, edge


def recompute(d, attrs, rows, runs_dir, satellite_dir, model_dir):
    """
    Every index of the particles `rows`, by plain loops. Returns a dict,
    and which particles lie on a cell edge on a day of their windows.
    """
    release = date(*map(int, str(attrs["release_date"])[:10].split("-")))
    with xr.open_dataset(Path(runs_dir) / attrs["trajectory_source"],
                         decode_timedelta=False) as ds:
        sub = ds[["lon", "lat", "z"]].isel(trajectory=rows).load()
        n_obs = ds.sizes["obs"]
    lon, lat, z = (sub[v].values.astype(np.float64) for v in ("lon", "lat", "z"))
    satellite = xr.open_dataset(Path(satellite_dir) / chlorophyll.SATELLITE_FILE)
    month_day = list(zip(satellite.time.dt.month.values, satellite.time.dt.day.values))
    m = len(rows)
    n_w, n_c = len(chlorophyll.WINDOWS), len(chlorophyll.MIN_COUNTS)
    names2 = ("surface", "layer", "depth")
    names3 = ("mean", "median", "matched")
    total = {k: np.zeros((m, n_w)) for k in names2}
    total.update({k: np.zeros((m, n_w, n_c)) for k in names3})
    days = {k: np.zeros((m, n_w), dtype=int) for k in ("window", "surface", "depth")}
    days.update({k: np.zeros((m, n_w, n_c), dtype=int) for k in ("satellite", "matched")})
    model, model_month = None, None
    bounds = None
    on_edge = np.zeros(m, dtype=bool)
    for k in range(n_obs):
        today = release + timedelta(days=k)
        window = np.full(m, -1)
        for q, p in enumerate(rows):
            if not (np.isfinite(lon[q, k]) and np.isfinite(lat[q, k])
                    and np.isfinite(z[q, k])):
                continue
            end, advance = int(d["calyptopis_end_day"][p]), int(d["advance_day"][p])
            if end == NONE or k < end:
                window[q] = 0
            elif advance == NONE or k < advance:
                window[q] = 1
        active = np.flatnonzero(window >= 0)
        if active.size == 0:
            continue
        x, y = lon[active, k], lat[active, k]
        # Satellite: the day of the year.
        i = month_day.index((today.month, today.day))
        sat = satellite.isel(time=i).load()
        s_mean, edge = _nearest(sat["CHL_mean"], x, y)
        s_median, _ = _nearest(sat["CHL_median"], x, y)
        s_count, _ = _nearest(sat["CHL_count"], x, y)
        on_edge[active] |= edge
        # Model: the actual date.
        if model_month != (today.year, today.month):
            if model is not None:
                model.close()
            model = xr.open_dataset(Path(model_dir) / chlorophyll.MODEL_PATTERN.format(
                year=today.year, month=today.month))
            model_month = (today.year, today.month)
            if bounds is None:
                centre = model.depth.values.astype(np.float64)
                bounds = [0.0] + [0.5 * (centre[j] + centre[j + 1])
                                  for j in range(centre.size - 1)]
                bounds.append(2 * centre[-1] - bounds[-1])
        field = model["chl"].sel(time=np.datetime64(today.isoformat())).load()
        column, edge = _nearest(field, x, y)                 # (depth, p)
        on_edge[active] |= edge
        on_edge[active] |= np.abs(z[active, k][:, None]
                                  - np.asarray(bounds)[None, :]).min(axis=1) < 1e-6
        for a, q in enumerate(active):
            w = window[q]
            days["window"][q, w] += 1
            col = column[:, a]
            has_model = np.isfinite(col[0])
            if has_model:
                top = thickness = 0.0
                for j in range(col.size):
                    part = min(bounds[j + 1], chlorophyll.LAYER_DEPTH) \
                        - min(bounds[j], chlorophyll.LAYER_DEPTH)
                    if part > 0 and np.isfinite(col[j]):
                        top += col[j] * part
                        thickness += part
                total["surface"][q, w] += col[0]
                total["layer"][q, w] += top / thickness
                days["surface"][q, w] += 1
                level = int(np.argmin(np.abs(centre - z[q, k])))
                if z[q, k] >= 0 and z[q, k] < bounds[-1] and np.isfinite(col[level]):
                    total["depth"][q, w] += col[level]
                    days["depth"][q, w] += 1
            for c, minimum in enumerate(chlorophyll.MIN_COUNTS):
                if (np.isfinite(s_mean[a]) and np.isfinite(s_median[a])
                        and s_mean[a] >= 0 and s_median[a] >= 0
                        and np.isfinite(s_count[a]) and s_count[a] >= minimum):
                    total["mean"][q, w, c] += s_mean[a]
                    total["median"][q, w, c] += s_median[a]
                    days["satellite"][q, w, c] += 1
                    if has_model:
                        total["matched"][q, w, c] += col[0]
                        days["matched"][q, w, c] += 1
    satellite.close()
    if model is not None:
        model.close()

    def mean(name, count):
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(count > 0, total[name] / count, np.nan)
    return on_edge, {
        "window_days": days["window"],
        "chl_satellite_mean": mean("mean", days["satellite"]),
        "chl_satellite_median": mean("median", days["satellite"]),
        "chl_satellite_days": days["satellite"],
        "chl_model_surface": mean("surface", days["surface"]),
        "chl_model_mean_200m": mean("layer", days["surface"]),
        "chl_model_days": days["surface"],
        "chl_model_at_depth": mean("depth", days["depth"]),
        "chl_model_at_depth_days": days["depth"],
        "chl_model_surface_matched": mean("matched", days["matched"]),
        "chl_matched_days": days["matched"],
    }


def pick_sample(d, n_random, n_obs, rng):
    n = d["particle"].size
    rows = set(rng.choice(n, min(n_random, n), replace=False).tolist())
    cases = {
        "no advance": d["advance_day"] == NONE,
        "furcilia I not reached": d["calyptopis_end_day"] == NONE,
        "advance during the calyptopis stages":
            (d["advance_day"] != NONE) & (d["advance_day"] <= d["calyptopis_end_day"]),
        "deleted during the calyptopis stages":
            d["window_days"][:, 0].astype(int)
            < np.where(d["calyptopis_end_day"] == NONE, n_obs, d["calyptopis_end_day"]),
        "model without a value at depth on some days":
            (d["chl_model_at_depth_days"] < d["chl_model_days"]).any(axis=1),
        "model without a value on some days":
            (d["chl_model_days"] < d["window_days"]).any(axis=1),
        "satellite without any value in a window":
            ((d["chl_satellite_days"][:, :, 0] == 0) & (d["window_days"] > 0)).any(axis=1),
        "minimum count removes days":
            (d["chl_satellite_days"][:, :, 1] < d["chl_satellite_days"][:, :, 0]).any(axis=1),
    }
    found = {}
    for label, mask in cases.items():
        some = np.flatnonzero(mask)[:15]
        found[label] = int(mask.sum())
        rows |= set(some.tolist())
    return np.array(sorted(rows)), found


def recomputation(d, attrs, n_random, n_obs, rng):
    directories = [os.environ.get(k) for k in
                   ("KRICO_RUNS", "KRICO_GLOBCOLOUR", "KRICO_PISCES")]
    if not all(directories):
        print("  KRICO_RUNS, KRICO_GLOBCOLOUR or KRICO_PISCES not set: "
              "indices NOT recomputed")
        return ["the indices could not be recomputed (environment)"]
    rows, cases = pick_sample(d, n_random, n_obs, rng)
    print(f"  recomputing {rows.size} particles; cases in the whole file: "
          + "; ".join(f"{label} {count}" for label, count in cases.items()))
    on_edge, reference = recompute(d, attrs, rows, *directories)
    keep = ~on_edge
    print(f"  {int(on_edge.sum())} of them lie within {EDGE} degree of a cell "
          f"edge on a day of their windows and are left out; "
          f"{int(keep.sum())} compared")
    reasons = []
    if keep.sum() < 0.5 * rows.size:
        reasons.append("more than half of the sample lies on a cell edge")
    for name, ref in reference.items():
        stored, ref = d[name][rows][keep], ref[keep]
        if name.endswith("days"):
            bad = int((stored.astype(int) != ref).sum())
            total = ref.size
        else:
            flat_s, flat_r = stored.reshape(-1), ref.reshape(-1)
            bad = sum(not _same(s, r) for s, r in zip(flat_s, flat_r))
            total = flat_r.size
        print(f"  {name:28s} {bad} of {total} values differ")
        if bad:
            reasons.append(f"{name}: {bad} values differ from the recomputation")
    return reasons


# -------------------------------------------------------------------- report

def describe(d):
    n = d["particle"].size

    def line(label, value, days, length):
        has = days > 0
        if not has.any():
            print(f"    {label:44s} no particle with a value")
            return
        q = np.percentile(value[has].astype(np.float64), (5, 50, 95))
        print(f"    {label:44s} {100 * has.sum() / n:5.1f}% of particles, "
              f"{100 * np.median(days[has] / length[has]):5.1f}% of their days "
              f"(median); mg m-3: 5% {q[0]:.3f}, median {q[1]:.3f}, 95% {q[2]:.3f}")

    for w, window in enumerate(chlorophyll.WINDOWS):
        length = d["window_days"][:, w].astype(np.float64)
        some = length > 0
        print(f"  window {window}: {100 * some.mean():.1f}% of particles have "
              f"at least one day; days 5% {np.percentile(length[some], 5):.0f}, "
              f"median {np.median(length[some]):.0f}, 95% "
              f"{np.percentile(length[some], 95):.0f}" if some.any() else
              f"  window {window}: no particle has a day")
        line("model, surface", d["chl_model_surface"][:, w], d["chl_model_days"][:, w], length)
        line("model, top 200 m", d["chl_model_mean_200m"][:, w], d["chl_model_days"][:, w], length)
        line("model, at depth", d["chl_model_at_depth"][:, w],
             d["chl_model_at_depth_days"][:, w], length)
        for c, minimum in enumerate(chlorophyll.MIN_COUNTS):
            line(f"satellite mean, count >= {minimum}",
                 d["chl_satellite_mean"][:, w, c], d["chl_satellite_days"][:, w, c], length)
            line(f"satellite median, count >= {minimum}",
                 d["chl_satellite_median"][:, w, c], d["chl_satellite_days"][:, w, c], length)
        for c, minimum in enumerate(chlorophyll.MIN_COUNTS):
            mean = d["chl_satellite_mean"][:, w, c].astype(np.float64)
            median = d["chl_satellite_median"][:, w, c].astype(np.float64)
            model = d["chl_model_surface_matched"][:, w, c].astype(np.float64)
            same_days = (d["chl_matched_days"][:, w, c] == d["chl_satellite_days"][:, w, c]) \
                & (d["chl_matched_days"][:, w, c] > 0)
            has = np.isfinite(mean)
            if has.any():
                r = np.percentile(mean[has] / median[has], (5, 50, 95))
                print(f"    count >= {minimum:2d}: satellite mean / median: 5% "
                      f"{r[0]:.2f}, median {r[1]:.2f}, 95% {r[2]:.2f}")
            if same_days.sum() > 2:
                r = np.percentile(model[same_days] / mean[same_days], (5, 50, 95))
                rho = np.corrcoef(np.log(model[same_days]), np.log(mean[same_days]))[0, 1]
                print(f"    count >= {minimum:2d}: model surface / satellite mean "
                      f"on the same days ({100 * same_days.mean():.1f}% of "
                      f"particles): 5% {r[0]:.2f}, median {r[1]:.2f}, 95% "
                      f"{r[2]:.2f}; correlation of the logarithms {rho:.2f}")


def sizes(path, d):
    n = d["particle"].size
    stored = _storage(path)
    total = path.stat().st_size
    if stored is None:
        print(f"  h5py not available: {total} bytes in all, {total / n:.2f} per particle")
        return
    print("  bytes on disk per variable:")
    for name in d:
        print(f"    {name:28s} {stored[name]:11d}  {stored[name] / n:6.2f} per "
              f"particle  (uncompressed {d[name].nbytes / n:5.1f})")
    print(f"    {'file':28s} {total:11d}  {total / n:6.2f} per particle")


def report(path, pass1_dir, n_random, rng):
    path = Path(path)
    print(f"== {path}")
    d, attrs = _read(path)
    n_obs = inputs.N_OBS
    print(f"  {d['particle'].size} particles, release {attrs.get('release_date')}, "
          f"code {attrs.get('code_version')} ({str(attrs.get('code_commit'))[:7]}), "
          f"dataset {attrs.get('dataset_version')}")
    print(f"  satellite: {attrs.get('satellite_source')}")
    print(f"  model: {attrs.get('model_source')}")
    reasons = []
    names = [v for v in d if v not in ("particle", "window", "min_count")]
    if names != krico_io.chlorophyll_variable_names():
        reasons.append("variables differ from the schema: "
                       f"{sorted(set(names) ^ set(krico_io.chlorophyll_variable_names()))}")
    if tuple(d["min_count"].tolist()) != chlorophyll.MIN_COUNTS:
        reasons.append(f"min_count is {d['min_count'].tolist()}")
    print(" consistency, every particle:")
    reasons += consistency(d, n_obs)
    print(" windows against pass 1:")
    reasons += against_pass1(d, path, pass1_dir)
    print(" recomputation by another route:")
    reasons += recomputation(d, attrs, n_random, n_obs, rng)
    print(" what the file holds:")
    describe(d)
    sizes(path, d)
    ok = not reasons
    print("  verdict: " + ("PASS" if ok else "FAIL, " + "; ".join(reasons)))
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("files", nargs="+", type=Path, help="pass 2 files")
    parser.add_argument("--pass1", type=Path, default=DEFAULT_PASS1, metavar="DIR",
                        help="folder of the pass 1 files with GLORYS12 sea ice "
                             "(default: recruitment/data_v3.0_glorys12)")
    parser.add_argument("--sample", type=int, default=300,
                        help="random particles recomputed per file (default 300)")
    args = parser.parse_args()
    rng = np.random.default_rng(2026)
    ok = True
    try:
        for path in args.files:
            ok &= report(path, args.pass1, args.sample, rng)
    except (RuntimeError, ValueError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: {exc}")
    print("ALL PASS" if ok else "AT LEAST ONE FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
