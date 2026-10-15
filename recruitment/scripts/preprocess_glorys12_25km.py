"""
Bin the GLORYS12 monthly sea-ice concentration onto the NSIDC 25 km
polar-stereographic grid, and write it in the format of the NSIDC monthly
stacks, so that the classification can read GLORYS12 ice exactly as it
reads NSIDC ice: by containing 25 km cell, with the same land/coast fill.

This is the third ice source of the forcing-change test. With it,
    native GLORYS12  minus  GLORYS12 on 25 km   = resolution and sampling
    GLORYS12 on 25 km  minus  NSIDC             = product difference at
                                                   equal resolution

Usage
-----
    python3 preprocess_glorys12_25km.py YEAR      # months of YEAR within the period
    python3 preprocess_glorys12_25km.py --check   # static map diagnostics only

Input:   $KRICO_GLORYS12/glorys12_ice_YYYY_MM.nc  (ileadfra(time_counter, y, x),
         nav_lon/nav_lat (y, x), daily at 12:00)
         $KRICO_NSIDC/ancillary/G02202-ancillary-pss25-v06r00.nc  (surface_type)
Output:  $KRICO_GLORYS12_25KM/glorys12_25km_ice_YYYY_MM.nc
         $KRICO_GLORYS12_25KM/ancillary/  (copy of the NSIDC ancillary file, so
         the directory is self-contained for the sampler and the validation)

Method. Every GLORYS12 cell centre is assigned once to the NSIDC cell that
contains it (krico_recruitment.nsidc_grid); each day, the concentration of
a 25 km cell is the plain mean of the GLORYS12 ocean cells it contains
(land cells, NaN in GLORYS12, are left out; latitude varies too little
within 25 km for area weighting to matter). Cells that are land, coast or
lake in the NSIDC surface_type mask are set to fill (255) whatever GLORYS12
holds there, so that a trajectory position reads NaN in exactly the cells
where it reads NaN from NSIDC. Cells outside the GLORYS12 domain (115W-40E,
south of 40S) are fill too; no trajectory is there. Values are packed as
ubyte 0-100 (scale 0.01) like the CDR.

Period: 1993-10-01 to 2025-10-31, as for NSIDC. A month without a GLORYS12
file is skipped with a message; the sampler then fails loudly for any cohort
that needs it, so check that glorys12_ice_1993_10.nc exists (the first
cohorts' M1 reads 22 October 1993) before submitting the run.

Environment
-----------
    KRICO_GLORYS12        required; GLORYS12 monthly ice
    KRICO_NSIDC           required; for the ancillary mask
    KRICO_GLORYS12_25KM   required; output, e.g. /scratch/cvan/KRICO/Pre/GLORYS12_25KM
"""

import argparse
import calendar
import os
import shutil
import sys
from datetime import date, datetime, timezone
from pathlib import Path

import netCDF4
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from krico_recruitment import nsidc_grid  # noqa: E402

DATE_MIN = date(1993, 10, 1)
DATE_MAX = date(2025, 10, 31)
GLORYS_PATTERN = "glorys12_ice_{year:04d}_{month:02d}.nc"
OUT_PATTERN = "glorys12_25km_ice_{year:04d}_{month:02d}.nc"
ANCILLARY = "G02202-ancillary-pss25-v06r00.nc"
VAR = "cdr_seaice_conc"
FILL = 255
SURFACE_OCEAN = 50
EPOCH = date(1970, 1, 1)
NCELL = nsidc_grid.NROW * nsidc_grid.NCOL


def _env(name):
    v = os.environ.get(name)
    if not v:
        sys.exit(f"ERROR: {name} environment variable not set.")
    return Path(v)


def _static_map(glorys_dir, nsidc_dir):
    """
    Flat NSIDC cell index of every GLORYS12 cell (-1 off the NSIDC grid),
    the NSIDC ocean mask (flat, bool), and the GLORYS12 grid shape.
    """
    sample = sorted(glorys_dir.glob("glorys12_ice_*.nc"))
    if not sample:
        sys.exit(f"ERROR: no glorys12_ice_*.nc in {glorys_dir}")
    with xr.open_dataset(sample[0]) as ds:
        lon = np.asarray(ds["nav_lon"].values, dtype="float64")
        lat = np.asarray(ds["nav_lat"].values, dtype="float64")
    if lon.ndim != 2:
        sys.exit("ERROR: expected 2-D nav_lon/nav_lat")
    row, col = nsidc_grid.lonlat_to_index(lon, lat)
    flat = np.where((row >= 0) & (col >= 0), row * nsidc_grid.NCOL + col, -1)
    flat = flat.astype(np.int64).reshape(-1)

    with xr.open_dataset(nsidc_dir / "ancillary" / ANCILLARY, mask_and_scale=False) as ds:
        surface = np.asarray(ds["surface_type"].values)
    ocean = (surface == SURFACE_OCEAN).reshape(-1)
    return flat, ocean, lon.shape


def bin_day(field, flat, ocean):
    """
    One daily GLORYS12 field (y, x), float with NaN over land -> packed
    (NROW, NCOL) ubyte on the NSIDC grid.
    """
    valid = np.isfinite(field.reshape(-1)) & (flat >= 0)
    idx = flat[valid]
    s = np.bincount(idx, weights=field.reshape(-1)[valid], minlength=NCELL)
    n = np.bincount(idx, minlength=NCELL)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s / n
    packed = np.full(NCELL, FILL, dtype=np.uint8)
    good = (n > 0) & ocean
    packed[good] = np.clip(np.rint(100.0 * mean[good]), 0, 100).astype(np.uint8)
    return packed.reshape(nsidc_grid.NROW, nsidc_grid.NCOL), good


def process_month(glorys_dir, out_dir, flat, ocean, year, month):
    ndays = calendar.monthrange(year, month)[1]
    days = [date(year, month, d) for d in range(1, ndays + 1)]
    days = [d for d in days if DATE_MIN <= d <= DATE_MAX]
    if not days:
        return None
    src = glorys_dir / GLORYS_PATTERN.format(year=year, month=month)
    if not src.is_file():
        print(f"  {src.name} absent, skipped (no GLORYS12 ice for this month)")
        return None
    out = out_dir / OUT_PATTERN.format(year=year, month=month)
    if out.is_file():
        print(f"  {out.name} exists, skipped")
        return out

    with xr.open_dataset(src) as ds:
        times = ds["time_counter"].values
        got = sorted({np.datetime64(t, "D").astype(object) for t in times})
        if got != days:
            sys.exit(f"ERROR: {src.name} covers {got[0]}..{got[-1]} ({len(got)} days), "
                     f"expected {days[0]}..{days[-1]} ({len(days)})")
        conc = np.empty((len(days), nsidc_grid.NROW, nsidc_grid.NCOL), dtype=np.uint8)
        for k, d in enumerate(days):
            target = np.datetime64(f"{d:%Y-%m-%d}T12:00:00")
            field = np.asarray(ds["ileadfra"].sel(time_counter=target, method="nearest").values,
                               dtype="float64").squeeze()
            conc[k], good = bin_day(field, flat, ocean)

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
        vx[:] = nsidc_grid.X_CENTRES
        vy = ds.createVariable("y", "f8", ("y",))
        vy.standard_name = "projection_y_coordinate"
        vy.units = "meters"
        vy.axis = "Y"
        vy[:] = nsidc_grid.Y_CENTRES
        vc = ds.createVariable("crs", "i4")
        vc.grid_mapping_name = "polar_stereographic"
        vc.latitude_of_projection_origin = -90.0
        vc.standard_parallel = -70.0
        vc.straight_vertical_longitude_from_pole = 0.0
        vc.semi_major_axis = nsidc_grid.A
        vc.inverse_flattening = nsidc_grid.INV_F
        vc.srid = "urn:ogc:def:crs:EPSG::3412"
        vc.GeoTransform = "-3950000 25000 0 4350000 0 -25000"
        vs = ds.createVariable(VAR, "u1", ("time", "y", "x"), zlib=True, complevel=4,
                               fill_value=FILL,
                               chunksizes=(1, nsidc_grid.NROW, nsidc_grid.NCOL))
        vs.set_auto_maskandscale(False)
        vs.standard_name = "sea_ice_area_fraction"
        vs.long_name = ("GLORYS12v1 sea-ice concentration (ileadfra), mean of the 1/12 deg "
                        "cells inside each 25 km NSIDC cell, NSIDC land/coast/lake as fill")
        vs.units = "1"
        vs.scale_factor = 0.01
        vs.add_offset = 0.0
        vs.valid_range = np.array([0, 100], dtype=np.uint8)
        vs.grid_mapping = "crs"
        vs.coordinates = "time y x"
        vs[:] = conc
        vsat = ds.createVariable("satellite", str, ("time",))
        vsat.long_name = "source tag; 'glorys12' for every day of this product"
        for k in range(len(days)):
            vsat[k] = "glorys12"
        ds.title = ("GLORYS12v1 sea-ice concentration binned onto the NSIDC 25 km "
                    "polar-stereographic south grid, daily, stacked by month for KRICO")
        ds.source_product = "GLORYS12v1 ileadfra (Mercator Ocean), KRICO domain subset"
        ds.source_file = src.name
        ds.setncattr("mask", f"NSIDC G02202 v6 ancillary surface_type != {SURFACE_OCEAN} -> fill")
        ds.created = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        ds.created_by = "krico-post-production/recruitment/scripts/preprocess_glorys12_25km.py"
    os.replace(tmp, out)
    nfill = float((conc == FILL).mean())
    print(f"  {out.name}: {len(days)} days, fill fraction {nfill:.3f}, "
          f"ocean cells with data {int(good.sum())}")
    return out


def check(glorys_dir, nsidc_dir):
    flat, ocean, shape = _static_map(glorys_dir, nsidc_dir)
    on_grid = flat >= 0
    n_per = np.bincount(flat[on_grid], minlength=NCELL)
    covered = (n_per > 0) & ocean
    print(f"GLORYS12 grid {shape}: {on_grid.sum()} of {flat.size} cells fall on the NSIDC grid")
    print(f"NSIDC ocean cells receiving GLORYS12 cells: {covered.sum()} of {ocean.sum()}")
    print(f"GLORYS12 cells per covered NSIDC cell: min {n_per[covered].min()}, "
          f"median {int(np.median(n_per[covered]))}, max {n_per[covered].max()}")
    # expected ~ (25 km / (1/12 deg * 111 km * cos lat))^2: 20-45 depending on latitude


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("year", nargs="?", type=int)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    glorys_dir = _env("KRICO_GLORYS12")
    nsidc_dir = _env("KRICO_NSIDC")
    if args.check:
        check(glorys_dir, nsidc_dir)
        return
    if args.year is None:
        parser.error("give a year or --check")
    if not DATE_MIN.year <= args.year <= DATE_MAX.year:
        sys.exit(f"year must be in {DATE_MIN.year}-{DATE_MAX.year}")
    out_dir = _env("KRICO_GLORYS12_25KM")
    out_dir.mkdir(parents=True, exist_ok=True)
    anc_dst = out_dir / "ancillary" / ANCILLARY
    if not anc_dst.is_file():
        anc_dst.parent.mkdir(exist_ok=True)
        shutil.copy2(nsidc_dir / "ancillary" / ANCILLARY, anc_dst)
    flat, ocean, _ = _static_map(glorys_dir, nsidc_dir)
    print(f"{args.year}:")
    for month in range(1, 13):
        first = date(args.year, month, 1)
        last = date(args.year, month, calendar.monthrange(args.year, month)[1])
        if first > DATE_MAX or last < DATE_MIN:
            continue
        process_month(glorys_dir, out_dir, flat, ocean, args.year, month)


if __name__ == "__main__":
    main()
