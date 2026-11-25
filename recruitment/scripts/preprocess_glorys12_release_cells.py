"""
Extract the daily GLORYS12 sea-ice concentration and temperature column at
the release cells of the KRICO simulations, one file per spawning year.

This prepares the inputs of the rules that act before release (M1 at
spawning, M3 during the descent and the ascent), of the thermal limits at
spawning, and of the no-transport control. It evaluates no rule and reads no
trajectory. The values are copied unchanged (float32, lossless compression),
so that a trajectory's own sea ice and temperature at day 0 can be compared
with them exactly (validate_glorys12_release_cells.py).

Usage
-----
    python3 preprocess_glorys12_release_cells.py --check        # inputs of the 32 years
    python3 preprocess_glorys12_release_cells.py YEAR --check   # inputs of one year
    python3 preprocess_glorys12_release_cells.py YEAR           # extract one spawning year

Input:   $KRICO_GLORYS12/glorys12_ice_YYYY_MM.nc   ileadfra(time_counter, y, x)
         $KRICO_GLORYS12/glorys12_temp_YYYY_MM.nc  votemper(time_counter, deptht, y, x)
         $KRICO_GLORYS12/glorys12_w_YYYY_MM.nc     depthw, the level bounds (any month)
         $KRICO_GLORYS12/glorys12_bathymetry.nc    hdept, nav_lon, nav_lat
         glorys12_release_cells.csv                the 35,303 release cells
Output:  $KRICO_GLORYS12_RELEASE_CELLS/glorys12_release_cells_YYYY.nc

Cells. Every release cell and its eight neighbours, 55,034 cells in all, so
that the grid point nearest to any release position is among them
(krico_recruitment/release_cells.py explains why the neighbours are needed).
The cell dimension is sorted by (j, i) and is the same in every file; the
classification links a particle to a cell by its index along it.

Period. Spawning year Y runs from 22 October Y-1 (first spawning date) to
30 September Y (last tracked day): 344 days, 345 when Y is a leap year. Days
are found by their date in the monthly files, not by their position:
glorys12_*_1993_10.nc holds 22 to 31 October only.

Levels. From the surface to the sea floor of the deepest release cell
(40 levels, the last one centred on 1941.9 m). A neighbour deeper than the
bottom of level 40 (2080 m) is cut there. Values below the sea floor and on
land are missing.

Output variables
----------------
    sea_ice_area_fraction(time, cell)   ileadfra, fraction 0-1
    temperature(time, depth, cell)      votemper, degrees Celsius
    time(time)                          as in the source, daily means at 12:00
    depth(depth), depth_bounds(depth_edge)   deptht, and depthw: level k lies
                                        between depth_bounds[k] and [k + 1]
    cell_y, cell_x, lon, lat, bathymetry, is_release_cell, n_release  (cell)
    grid_lon(grid_x), grid_lat(grid_y)  axes of the whole grid, needed to find
                                        the grid point nearest to a position

A year whose file exists is skipped. The file is written under a temporary
name and renamed when complete, so an interrupted task leaves no partial
file under the final name.

Environment
-----------
    KRICO_GLORYS12                 required; GLORYS12 monthly files
    KRICO_GLORYS12_RELEASE_CELLS   required; output,
                                   e.g. /scratch/cvan/KRICO/Pre/GLORYS12_RELEASE_CELLS
"""

import argparse
import os
import sys
import time
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from krico_recruitment import provenance, release_cells  # noqa: E402

ICE_PATTERN = "glorys12_ice_{year:04d}_{month:02d}.nc"
TEMP_PATTERN = "glorys12_temp_{year:04d}_{month:02d}.nc"
ICE_VAR = "ileadfra"
TEMP_VAR = "votemper"
FILL = np.float32(netCDF4.default_fillvals["f4"])
MISSING_ABOVE = 1e30       # the source marks land and sub-floor values so
COMPLEVEL = 4


def _env(name):
    v = os.environ.get(name)
    if not v:
        sys.exit(f"ERROR: {name} environment variable not set.")
    return Path(v)


def _months(days):
    """Days grouped by calendar month, in order: [((year, month), [days])]."""
    groups = []
    for d in days:
        if not groups or groups[-1][0] != (d.year, d.month):
            groups.append(((d.year, d.month), []))
        groups[-1][1].append(d)
    return groups


def _records_by_day(nc, name):
    """
    {date: (record index, raw time value)} of a monthly file, with the units
    and calendar of its time axis. Each record must be a daily mean stamped
    at 12:00, one per day.
    """
    var = nc["time_counter"]
    raw = np.asarray(var[:], dtype=np.float64)
    calendar = getattr(var, "calendar", "standard")
    stamps = netCDF4.num2date(raw, var.units, calendar=calendar,
                              only_use_cftime_datetimes=False,
                              only_use_python_datetimes=True)
    records = {}
    for k, stamp in enumerate(np.atleast_1d(stamps)):
        if (stamp.hour, stamp.minute, stamp.second) != (12, 0, 0):
            sys.exit(f"ERROR: {name}: record {k} is stamped {stamp}, not 12:00")
        if stamp.date() in records:
            sys.exit(f"ERROR: {name}: two records for {stamp.date()}")
        records[stamp.date()] = (k, raw[k])
    return records, var.units, calendar


def _check_variable(var, name, dims, shape_tail):
    if var.dimensions != dims:
        sys.exit(f"ERROR: {name}: dimensions {var.dimensions}, expected {dims}")
    if var.shape[-len(shape_tail):] != shape_tail:
        sys.exit(f"ERROR: {name}: shape {var.shape}, expected (..., "
                 f"{', '.join(map(str, shape_tail))})")
    if var.dtype.kind != "f":
        sys.exit(f"ERROR: {name}: stored as {var.dtype}, expected floating point")
    # The Mercator files carry scale_factor = 1 and add_offset = 0 on float32
    # values: the stored values are the physical ones. Any other packing
    # would make a raw copy wrong.
    for attr, neutral in (("scale_factor", 1.0), ("add_offset", 0.0)):
        if attr in var.ncattrs() and float(var.getncattr(attr)) != neutral:
            sys.exit(f"ERROR: {name}: packed values ({attr} = "
                     f"{var.getncattr(attr)}); the extraction copies raw "
                     "values and cannot read them")
    var.set_auto_maskandscale(False)


def _static(glorys_dir, table_path):
    grid = release_cells.read_grid(glorys_dir)
    cells = release_cells.build_cells(grid, release_cells.read_table(table_path))
    nlev = release_cells.number_of_levels(grid, cells)
    if nlev >= len(grid["depthw"]):
        sys.exit("ERROR: the deepest release cell reaches the last level; "
                 "its lower bound is not in depthw")
    return grid, cells, nlev


def _take(block, rows, cols):
    """Values of a (..., rows, x) block at the cells, missing set to FILL."""
    values = block[..., rows, cols].astype(np.float32)
    values[~(np.abs(values) < MISSING_ABOVE)] = FILL      # also catches NaN
    return values


def extract_year(year, glorys_dir, out_dir, table_path, allow_dirty=False):
    out = out_dir / release_cells.FILE_PATTERN.format(year=year)
    if out.is_file():
        print(f"{out.name} exists, skipped")
        return out
    identity = provenance.code_identity(script=__file__, allow_dirty=allow_dirty)
    started = time.time()
    grid, cells, nlev = _static(glorys_dir, table_path)
    ny, nx = grid["bathymetry"].shape
    ncell = len(cells["cell_y"])
    days = release_cells.days(year)
    # Rows are read over the full width of the grid, between the first and
    # the last row holding a cell: one contiguous block per level.
    j0, j1 = int(cells["cell_y"].min()), int(cells["cell_y"].max()) + 1
    rows = cells["cell_y"] - j0
    cols = cells["cell_x"]
    print(f"Spawning year {year}: {days[0]} to {days[-1]}, {len(days)} days")
    print(f"  cells: {ncell} ({int(cells['is_release'].sum())} release cells "
          f"and their neighbours), levels: {nlev} "
          f"(to {grid['depthw'][nlev]:.1f} m), grid rows {j0} to {j1 - 1}")
    print(f"  code: {identity['code_version']} ({identity['code_commit'][:7]})")

    tmp = out.with_suffix(".nc.part")
    sources = []
    with netCDF4.Dataset(tmp, "w", format="NETCDF4") as ds:
        ds.createDimension("time", len(days))
        ds.createDimension("depth", nlev)
        ds.createDimension("depth_edge", nlev + 1)
        ds.createDimension("cell", ncell)
        ds.createDimension("grid_y", ny)
        ds.createDimension("grid_x", nx)

        vtime = ds.createVariable("time", "f8", ("time",))
        vtime.standard_name = "time"
        vtime.long_name = "time of the daily mean (12:00)"
        vtime.axis = "T"

        v = ds.createVariable("depth", "f4", ("depth",))
        v.standard_name = "depth"
        v.long_name = "depth of the temperature levels (deptht)"
        v.units = "m"
        v.positive = "down"
        v.axis = "Z"
        v[:] = grid["deptht"][:nlev]
        v = ds.createVariable("depth_bounds", "f4", ("depth_edge",))
        v.long_name = ("bounds of the levels (depthw): level k lies between "
                       "depth_bounds[k] and depth_bounds[k + 1]; the last "
                       "wet level of a column ends at the sea floor")
        v.units = "m"
        v.positive = "down"
        v[:] = grid["depthw"][:nlev + 1]

        v = ds.createVariable("cell", "i4", ("cell",))
        v.long_name = ("cell index; the classification links a particle to "
                       "a cell by this index")
        v[:] = np.arange(ncell, dtype=np.int32)
        for name, dtype, long_name, units in (
            ("cell_y", "i4", "row of the cell in the GLORYS12 files (y, 0-based)", None),
            ("cell_x", "i4", "column of the cell in the GLORYS12 files (x, 0-based)", None),
            ("lon", "f4", "longitude of the cell centre (nav_lon)", "degrees_east"),
            ("lat", "f4", "latitude of the cell centre (nav_lat)", "degrees_north"),
            ("bathymetry", "f4", "sea-floor depth of the cell (hdept)", "m"),
        ):
            v = ds.createVariable(name, dtype, ("cell",), zlib=True, complevel=COMPLEVEL)
            v.long_name = long_name
            if units:
                v.units = units
            v[:] = cells[name]
        v = ds.createVariable("is_release_cell", "i1", ("cell",), zlib=True,
                              complevel=COMPLEVEL)
        v.long_name = ("1 where particles are released (sea floor of 1000 to "
                       "2000 m inside the CCAMLR subareas), 0 for a "
                       "neighbour of such a cell")
        v.flag_values = np.array([0, 1], dtype=np.int8)
        v.flag_meanings = "neighbour release_cell"
        v[:] = cells["is_release"]
        v = ds.createVariable("n_release", "i2", ("cell",), zlib=True,
                              complevel=COMPLEVEL)
        v.long_name = ("particles generated in the cell on each release day "
                       "(run.py of krico-templates, before the removal of "
                       "positions south of 76S); 0 for a neighbour")
        v[:] = cells["n_release"]
        v = ds.createVariable("grid_lon", "f8", ("grid_x",))
        v.long_name = "longitude axis of the whole GLORYS12 grid (nav_lon)"
        v.units = "degrees_east"
        v[:] = grid["grid_lon"]
        v = ds.createVariable("grid_lat", "f8", ("grid_y",))
        v.long_name = "latitude axis of the whole GLORYS12 grid (nav_lat)"
        v.units = "degrees_north"
        v[:] = grid["grid_lat"]

        vice = ds.createVariable("sea_ice_area_fraction", "f4", ("time", "cell"),
                                 zlib=True, shuffle=True, complevel=COMPLEVEL,
                                 chunksizes=(1, ncell), fill_value=FILL)
        vice.set_auto_maskandscale(False)
        vice.standard_name = "sea_ice_area_fraction"
        vice.long_name = "GLORYS12v1 sea-ice concentration (ileadfra), daily mean"
        vice.units = "1"
        vice.coordinates = "lon lat"
        vtemp = ds.createVariable("temperature", "f4", ("time", "depth", "cell"),
                                  zlib=True, shuffle=True, complevel=COMPLEVEL,
                                  chunksizes=(1, 1, ncell), fill_value=FILL)
        vtemp.set_auto_maskandscale(False)
        vtemp.standard_name = "sea_water_potential_temperature"
        vtemp.long_name = "GLORYS12v1 temperature (votemper), daily mean"
        vtemp.units = "degC"
        vtemp.coordinates = "lon lat"

        k = 0
        time_units = None
        for (yy, mm), month_days in _months(days):
            t0 = time.time()
            ice_path = glorys_dir / ICE_PATTERN.format(year=yy, month=mm)
            temp_path = glorys_dir / TEMP_PATTERN.format(year=yy, month=mm)
            for path in (ice_path, temp_path):
                if not path.is_file():
                    sys.exit(f"ERROR: {path} is missing")
            n_missing = np.zeros(2, dtype=np.int64)
            with netCDF4.Dataset(ice_path) as fi, netCDF4.Dataset(temp_path) as ft:
                vi, vt = fi[ICE_VAR], ft[TEMP_VAR]
                _check_variable(vi, ice_path.name, ("time_counter", "y", "x"), (ny, nx))
                _check_variable(vt, temp_path.name,
                                ("time_counter", "deptht", "y", "x"),
                                (len(grid["deptht"]), ny, nx))
                for nc, path in ((fi, ice_path), (ft, temp_path)):
                    difference = release_cells.coordinate_difference(
                        nc["nav_lon"][:], nc["nav_lat"][:], grid)
                    if not difference <= release_cells.COORDINATE_TOLERANCE:
                        sys.exit(f"ERROR: {path.name}: coordinates differ from "
                                 f"{grid['files']['ice']} by {difference:.2e} degrees")
                if not np.array_equal(np.asarray(ft["deptht"][:], dtype=np.float64),
                                      grid["deptht"]):
                    sys.exit(f"ERROR: {temp_path.name}: deptht differs from "
                             f"{grid['files']['temp']}")
                rec_i, units_i, cal_i = _records_by_day(fi, ice_path.name)
                rec_t, units_t, cal_t = _records_by_day(ft, temp_path.name)
                if (units_i, cal_i) != (units_t, cal_t):
                    sys.exit(f"ERROR: {ice_path.name} and {temp_path.name} "
                             "have different time units")
                if time_units is None:
                    time_units = (units_i, cal_i)
                    vtime.units, vtime.calendar = time_units
                elif time_units != (units_i, cal_i):
                    sys.exit(f"ERROR: {ice_path.name}: time units "
                             f"{units_i!r} differ from the previous months")
                for d in month_days:
                    for rec, path in ((rec_i, ice_path), (rec_t, temp_path)):
                        if d not in rec:
                            sys.exit(f"ERROR: {path.name} holds no record for {d}")
                    (ki, stamp_i), (kt, stamp_t) = rec_i[d], rec_t[d]
                    if stamp_i != stamp_t:
                        sys.exit(f"ERROR: {d}: sea-ice and temperature records "
                                 "have different times")
                    ice = _take(vi[ki, j0:j1, :], rows, cols)
                    temp = _take(vt[kt, :nlev, j0:j1, :], rows, cols)
                    n_missing[0] += int((ice == FILL).sum())
                    n_missing[1] += int((temp == FILL).sum())
                    vice[k, :] = ice
                    vtemp[k, :, :] = temp
                    vtime[k] = stamp_i
                    k += 1
            sources += [ice_path.name, temp_path.name]
            n = len(month_days)
            print(f"  {yy}-{mm:02d}: {n:2d} days in {time.time() - t0:6.1f} s, "
                  f"missing: sea ice {n_missing[0] / (n * ncell):.4f}, "
                  f"temperature {n_missing[1] / (n * ncell * nlev):.4f}",
                  flush=True)
        if k != len(days):
            sys.exit(f"ERROR: wrote {k} days, expected {len(days)}")

        ds.Conventions = "CF-1.7"
        ds.title = ("GLORYS12v1 daily sea-ice concentration and temperature "
                    "column at the release cells of the KRICO simulations")
        ds.spawning_year = np.int32(year)
        ds.first_day = str(days[0])
        ds.last_day = str(days[-1])
        ds.source_product = ("GLORYS12v1 (Mercator Ocean), daily means on the "
                             "native ORCA12 grid, KRICO domain subset; values "
                             "copied unchanged")
        ds.source_variables = f"{ICE_VAR}, {TEMP_VAR}"
        ds.source_files = " ".join(sources)
        ds.cell_definition = (
            "release cells of run.py in krico-templates (sea floor of 1000 "
            "to 2000 m, centre inside CCAMLR subareas 48.1 to 48.6 and 88.3) "
            "and their eight neighbours; a particle is linked to the grid "
            "point nearest to its release position")
        ds.release_cell_table = Path(table_path).name
        for key, value in identity.items():
            ds.setncattr(key, value)
        ds.history = provenance.history_line(
            f"preprocess_glorys12_release_cells.py {year}")
    os.replace(tmp, out)
    elapsed = time.time() - started
    size = out.stat().st_size
    raw = len(days) * ncell * (nlev + 1) * 4
    print(f"  {out.name}: {size / 1e9:.2f} GB on disk "
          f"({raw / 1e9:.2f} GB uncompressed), {elapsed / 60:.1f} min")
    return out


def check(glorys_dir, table_path, years):
    """Static inputs and monthly files of the spawning years; no output."""
    grid, cells, nlev = _static(glorys_dir, table_path)
    ny, nx = grid["bathymetry"].shape
    ncell = len(cells["cell_y"])
    release = cells["is_release"] == 1
    print(f"Grid: {ny} x {nx}, {len(grid['deptht'])} levels; coordinates of the "
          "bathymetry, sea-ice and temperature files "
          + ("are identical" if grid["coordinate_difference"] == 0 else
             f"differ by at most {grid['coordinate_difference']:.1e} degrees"))
    print(f"Release cells: {int(release.sum())}, table matches the bathymetry "
          f"file ({cells['bathymetry'][release].min():.1f} to "
          f"{cells['bathymetry'][release].max():.1f} m, "
          f"{int(cells['n_release'].sum())} particles per release day)")
    print(f"Cells to extract, with neighbours: {ncell}")
    print(f"Levels to extract: {nlev}, level {nlev} is centred on "
          f"{grid['deptht'][nlev - 1]:.1f} m and ends at {grid['depthw'][nlev]:.1f} m")
    deeper = int((cells["bathymetry"][~release] > grid["depthw"][nlev]).sum())
    print(f"Neighbours deeper than that, cut at level {nlev}: {deeper}")

    # Wet levels of the temperature field against the bathymetry: the level
    # holding the sea floor must be the last one with a value.
    temp_path = glorys_dir / grid["files"]["temp"]
    with netCDF4.Dataset(temp_path) as ft:
        vt = ft[TEMP_VAR]
        _check_variable(vt, temp_path.name, ("time_counter", "deptht", "y", "x"),
                        (len(grid["deptht"]), ny, nx))
        j0, j1 = int(cells["cell_y"].min()), int(cells["cell_y"].max()) + 1
        column = vt[0, :, j0:j1, :][:, cells["cell_y"] - j0, cells["cell_x"]]
        storage = "contiguous" if vt.chunking() == "contiguous" else f"chunks {vt.chunking()}"
        compressed = [name for name, on in vt.filters().items() if on and name != "complevel"]
        print(f"{temp_path.name}: {vt.dtype}, {storage}, "
              + ("compressed (" + ", ".join(compressed) + ")" if compressed
                 else "not compressed"))
    wet = np.abs(column) < MISSING_ABOVE
    n_wet = wet.sum(axis=0)
    expected = np.searchsorted(grid["depthw"], cells["bathymetry"], side="left")
    gaps = int((wet[:-1] < wet[1:]).any(axis=0).sum())
    same = n_wet[release] == expected[release]
    print(f"Wet levels against the bathymetry, release cells: {int(same.sum())} "
          f"of {int(release.sum())} agree; columns with a gap: {gaps}")
    if n_wet[release].max() > nlev:
        sys.exit(f"ERROR: a release cell has {n_wet[release].max()} wet levels, "
                 f"more than the {nlev} extracted")
    land = int((n_wet[~release] == 0).sum())
    print(f"Neighbours on land (no wet level): {land}")

    # Monthly files and their days, for every spawning year.
    print(f"Monthly files, spawning years {years[0]} to {years[-1]}:", flush=True)
    problems = []
    seen = {}
    for year in years:
        for (yy, mm), month_days in _months(release_cells.days(year)):
            for pattern in (ICE_PATTERN, TEMP_PATTERN):
                path = glorys_dir / pattern.format(year=yy, month=mm)
                if path not in seen:
                    if not path.is_file():
                        seen[path] = None
                    else:
                        with netCDF4.Dataset(path) as nc:
                            if pattern == ICE_PATTERN:
                                _check_variable(nc[ICE_VAR], path.name,
                                                ("time_counter", "y", "x"), (ny, nx))
                            else:
                                _check_variable(nc[TEMP_VAR], path.name,
                                                ("time_counter", "deptht", "y", "x"),
                                                (len(grid["deptht"]), ny, nx))
                            seen[path] = set(_records_by_day(nc, path.name)[0])
                if seen[path] is None:
                    problems.append(f"{path.name} is missing (spawning year {year})")
                else:
                    absent = [d for d in month_days if d not in seen[path]]
                    if absent:
                        problems.append(f"{path.name} lacks {absent[0]} to "
                                        f"{absent[-1]} (spawning year {year})")
    n_days = sum(len(release_cells.days(y)) for y in years)
    total = n_days * ncell * (nlev + 1) * 4
    if problems:
        print("\n".join(f"  {p}" for p in sorted(set(problems))))
        sys.exit(f"ERROR: {len(set(problems))} problem(s) in the monthly files")
    print(f"  {len(seen)} files present and readable as raw values, every "
          f"day of the {len(years)} period(s) found ({n_days} days)")
    print(f"Output, {len(years)} file(s): {total / 1e9:.1f} GB before compression "
          f"({total / len(years) / 1e9:.2f} GB per year)")
    print("Check passed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("year", nargs="?", type=int, help="spawning year, 1994-2025")
    parser.add_argument("--check", action="store_true",
                        help="check the inputs (of YEAR, or of the 32 years) "
                             "and write nothing")
    parser.add_argument("--cells", type=Path, default=release_cells.TABLE,
                        help="release-cell table (default: the committed one)")
    parser.add_argument("--allow-dirty", action="store_true",
                        help="for tests: accept modified or unpushed code")
    args = parser.parse_args()
    glorys_dir = _env("KRICO_GLORYS12")
    if args.year is not None and not (release_cells.FIRST_YEAR <= args.year
                                      <= release_cells.LAST_YEAR):
        sys.exit(f"year must be in {release_cells.FIRST_YEAR}-{release_cells.LAST_YEAR}")
    try:
        if args.check:
            years = ([args.year] if args.year is not None else
                     list(range(release_cells.FIRST_YEAR, release_cells.LAST_YEAR + 1)))
            check(glorys_dir, args.cells, years)
            return
        if args.year is None:
            parser.error("give a spawning year or --check")
        out_dir = _env("KRICO_GLORYS12_RELEASE_CELLS")
        out_dir.mkdir(parents=True, exist_ok=True)
        extract_year(args.year, glorys_dir, out_dir, args.cells,
                     allow_dirty=args.allow_dirty)
    except (RuntimeError, ValueError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: {exc}")


if __name__ == "__main__":
    main()
