"""
Stack the daily NSIDC G02202 v6 files into monthly files on the native grid,
one per calendar month, so that the recruitment classification can read ice
along trajectories the way it reads GLORYS12 ice: one small file per month,
indexed by day.

Usage
-----
    python3 preprocess_nsidc.py YEAR          # months of YEAR within the period
    python3 preprocess_nsidc.py --check-grid  # verify nsidc_grid against the
                                              # ancillary lat/lon, print masks

Input:   $KRICO_NSIDC/daily/YYYY/sic_pss25_YYYYMMDD_<sat>_v06r00.nc
Output:  $KRICO_NSIDC/nsidc_ice_YYYY_MM.nc

The monthly file keeps cdr_seaice_conc exactly as in the source: ubyte 0-100,
scale_factor 0.01 (so xarray decodes it to a 0-1 fraction), _FillValue 255,
dimensions (time, y, x) on the 332 x 316 grid, one time step per day of the
month. Nothing is masked, filled or regridded here; land, coast and lake
handling is the sampler's job, using the ancillary surface_type mask. The
satellite tag of every source file is recorded in the output so that the
sensor changes within the record (F11, F13, F17, am2) remain traceable.

The script refuses to write a month with a day missing: the provenance step
of the download established that the archive has one file per day of the
period, and a silent gap here would propagate into the along-trajectory
sampling as a NaN, which the filters treat as "below threshold".

Period: 1993-02-01 to 2026-02-28, as in download_nsidc.py.

Environment
-----------
    KRICO_NSIDC   required; e.g. /scratch/cvan/KRICO/Pre/NSIDC
"""

import argparse
import calendar
import os
import re
import sys
from datetime import date, datetime, timezone
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from krico_recruitment import nsidc_grid  # noqa: E402

DATE_MIN = date(1993, 2, 1)
DATE_MAX = date(2026, 2, 28)
ANCILLARY = "G02202-ancillary-pss25-v06r00.nc"
FILE_RE = re.compile(r"sic_pss25_(\d{8})_([A-Za-z0-9]+)_v06r00\.nc")
VAR = "cdr_seaice_conc"
EPOCH = date(1970, 1, 1)


def _root():
    root = os.environ.get("KRICO_NSIDC")
    if not root:
        sys.exit("ERROR: KRICO_NSIDC environment variable not set.")
    return Path(root)


def _daily_files(root, year, month):
    """{day: (path, satellite)} for one month, from the file names."""
    found = {}
    for path in (root / "daily" / str(year)).glob(f"sic_pss25_{year}{month:02d}*_v06r00.nc"):
        m = FILE_RE.fullmatch(path.name)
        if not m:
            continue
        d = date(int(m.group(1)[:4]), int(m.group(1)[4:6]), int(m.group(1)[6:8]))
        if d in found:
            sys.exit(f"ERROR: two files for {d}: {found[d][0].name}, {path.name}")
        found[d] = (path, m.group(2))
    return found


def process_month(root, year, month):
    ndays = calendar.monthrange(year, month)[1]
    days = [date(year, month, d) for d in range(1, ndays + 1)]
    days = [d for d in days if DATE_MIN <= d <= DATE_MAX]
    if not days:
        return None

    files = _daily_files(root, year, month)
    missing = [d for d in days if d not in files]
    if missing:
        sys.exit(f"ERROR: {year}-{month:02d}: no daily file for "
                 + ", ".join(str(d) for d in missing))

    out = root / f"nsidc_ice_{year}_{month:02d}.nc"
    if out.is_file():
        print(f"  {out.name} exists, skipped")
        return out

    conc = np.full((len(days), nsidc_grid.NROW, nsidc_grid.NCOL), 255, dtype=np.uint8)
    sats = []
    for k, d in enumerate(days):
        path, sat = files[d]
        with netCDF4.Dataset(path) as ds:
            ds.set_auto_maskandscale(False)      # keep raw ubyte + 255 fill
            v = ds[VAR]
            if v.shape != (1, nsidc_grid.NROW, nsidc_grid.NCOL):
                sys.exit(f"ERROR: {path.name}: unexpected shape {v.shape}")
            t = int(ds["time"][0])
            if t != (d - EPOCH).days:
                sys.exit(f"ERROR: {path.name}: time coordinate is "
                         f"{EPOCH.fromordinal(EPOCH.toordinal() + t)}, expected {d}")
            if abs(ds[VAR].scale_factor - 0.01) > 1e-12 or ds[VAR]._FillValue != 255:
                sys.exit(f"ERROR: {path.name}: unexpected packing of {VAR}")
            conc[k] = v[0]
            if k == 0:
                x = ds["x"][:]
                y = ds["y"][:]
                crs_attrs = {a: ds["crs"].getncattr(a) for a in ds["crs"].ncattrs()}
                product_version = ds.getncattr("product_version")
        sats.append(sat)

    if not (np.allclose(x, nsidc_grid.X_CENTRES) and np.allclose(y, nsidc_grid.Y_CENTRES)):
        sys.exit("ERROR: x/y coordinates differ from nsidc_grid cell centres")

    tmp = out.with_suffix(".nc.part")
    with netCDF4.Dataset(tmp, "w", format="NETCDF4") as ds:
        ds.createDimension("time", len(days))
        ds.createDimension("y", nsidc_grid.NROW)
        ds.createDimension("x", nsidc_grid.NCOL)

        vt = ds.createVariable("time", "i8", ("time",))
        vt.units = "days since 1970-01-01"
        vt.calendar = "standard"
        vt.standard_name = "time"
        vt.axis = "T"
        vt[:] = [(d - EPOCH).days for d in days]

        vx = ds.createVariable("x", "f8", ("x",))
        vx.standard_name = "projection_x_coordinate"
        vx.units = "meters"
        vx.axis = "X"
        vx[:] = x
        vy = ds.createVariable("y", "f8", ("y",))
        vy.standard_name = "projection_y_coordinate"
        vy.units = "meters"
        vy.axis = "Y"
        vy[:] = y

        vc = ds.createVariable("crs", "i4")
        vc.setncatts(crs_attrs)

        vs = ds.createVariable(VAR, "u1", ("time", "y", "x"),
                               zlib=True, complevel=4, fill_value=255,
                               chunksizes=(1, nsidc_grid.NROW, nsidc_grid.NCOL))
        vs.set_auto_maskandscale(False)
        vs.standard_name = "sea_ice_area_fraction"
        vs.long_name = "NOAA/NSIDC CDR of Passive Microwave Sea Ice Concentration"
        vs.units = "1"
        vs.scale_factor = 0.01
        vs.add_offset = 0.0
        vs.valid_range = np.array([0, 100], dtype=np.uint8)
        vs.grid_mapping = "crs"
        vs.coordinates = "time y x"
        vs[:] = conc

        vsat = ds.createVariable("satellite", str, ("time",))
        vsat.long_name = "source satellite tag of the daily file"
        for k, s in enumerate(sats):
            vsat[k] = s

        ds.title = ("NOAA/NSIDC CDR of Passive Microwave Sea Ice Concentration v6, "
                    "southern hemisphere, daily, stacked by month for KRICO")
        ds.source_product = f"G02202 {product_version}, doi:10.7265/b18j-z797"
        ds.source_files = " ".join(files[d][0].name for d in days)
        ds.created = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        ds.created_by = "krico-post-production/recruitment/scripts/preprocess_nsidc.py"
    os.replace(tmp, out)

    nfill = int((conc == 255).sum()) / conc.size
    print(f"  {out.name}: {len(days)} days, satellites {sorted(set(sats))}, "
          f"fill fraction {nfill:.3f}")
    return out


def process_year(year):
    root = _root()
    for month in range(1, 13):
        if date(year, month, 1) > DATE_MAX or date(year, month,
                calendar.monthrange(year, month)[1]) < DATE_MIN:
            continue
        process_month(root, year, month)


def check_grid():
    root = _root()
    with netCDF4.Dataset(root / "ancillary" / ANCILLARY) as ds:
        ds.set_auto_maskandscale(False)      # raw values; adj123 has no _FillValue
        lat = np.asarray(ds["latitude"][:])
        lon = np.asarray(ds["longitude"][:])
        surface = np.asarray(ds["surface_type"][:])
        adj = np.asarray(ds["adj123"][:])
        x = np.asarray(ds["x"][:])
        y = np.asarray(ds["y"][:])
    assert np.allclose(x, nsidc_grid.X_CENTRES) and np.allclose(y, nsidc_grid.Y_CENTRES), \
        "x/y coordinates of the ancillary file differ from nsidc_grid"

    rows, cols = np.meshgrid(np.arange(nsidc_grid.NROW), np.arange(nsidc_grid.NCOL),
                             indexing="ij")
    r, c = nsidc_grid.lonlat_to_index(lon, lat)
    ok = (r == rows) & (c == cols)
    xx, yy = nsidc_grid.lonlat_to_xy(lon, lat)
    dist = np.hypot(xx - nsidc_grid.X_CENTRES[cols], yy - nsidc_grid.Y_CENTRES[rows])
    print(f"cell centres mapping back to their own cell: {ok.sum()} / {ok.size}")
    print(f"max distance between projected centre and nominal centre: {np.nanmax(dist):.3f} m")
    if not ok.all():
        bad = np.argwhere(~ok)[:5]
        sys.exit(f"ERROR: grid mismatch, first failures (row, col): {bad.tolist()}")

    names = {50: "ocean", 75: "lake", 200: "coast", 250: "land"}
    print("surface_type cells:", {names.get(int(v), int(v)): int((surface == v).sum())
                                  for v in np.unique(surface)})
    print("adj123 cells:", {int(v): int((adj == v).sum()) for v in np.unique(adj)})


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("year", nargs="?", type=int)
    parser.add_argument("--check-grid", action="store_true")
    args = parser.parse_args()
    if args.check_grid:
        check_grid()
    elif args.year is not None:
        if not DATE_MIN.year <= args.year <= DATE_MAX.year:
            sys.exit(f"year must be in {DATE_MIN.year}-{DATE_MAX.year}")
        print(f"{args.year}:")
        process_year(args.year)
    else:
        parser.error("give a year or --check-grid")


if __name__ == "__main__":
    main()
