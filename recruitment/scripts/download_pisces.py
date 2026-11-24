"""
Download the daily 3D chlorophyll of the Copernicus Marine global
biogeochemical hindcast (GLOBAL_MULTIYEAR_BGC_001_029, PISCES, 1/4 degree,
75 levels) for the period and domain of the KRICO trajectories.

Usage
-----
    python3 download_pisces.py YEAR             # the months of one year, into $KRICO_PISCES
    python3 download_pisces.py YEAR --month M   # a single month
    python3 download_pisces.py YEAR --dry-run   # sizes only, nothing is downloaded
    python3 download_pisces.py --provenance     # after all years: what is on disk

Period: 1993-11 to 2025-10, the months covered by the GLORYS12 forcing of the
trajectories. Chlorophyll is sampled along the trajectories only, so no month
is needed before the first release (15 November 1993).

Domain: the GLORYS12 box of the trajectories (115W to 40E, south of 40S) with
a margin of 0.5 degree, so that every particle position has grid points on
both sides. All 75 depth levels are kept, because particles move in 3D.

Layout, mirroring the GLORYS12 and NSIDC preprocessing directories:
    $KRICO_PISCES/pisces_chl_YYYY_MM.nc   one file per month, daily means

Files are about 0.9 GB each, 384 in total (about 360 GB), uncompressed like
the GLORYS12 monthly files.

The download is resumable: a month already present with the expected number
of days and levels is skipped, and the toolbox writes to a temporary name that
is renamed only on success, so an interrupted job never leaves a truncated
file that a re-run would skip.

Unlike the NSIDC daily files, these files are subsets written on this side by
the toolbox, so a second download is not byte-identical and checksums would
verify nothing. --provenance records instead the product, DOI, dataset version,
toolbox version, subset and access date, with the days and size of each file,
in pisces_chl_provenance.txt (next to this script, tracked in git).

Access needs the Copernicus Marine toolbox (pip install copernicusmarine,
Python >= 3.10) and a free Copernicus Marine account, stored once on the
machine with:
    copernicusmarine login

Environment
-----------
    KRICO_PISCES     required; e.g. /scratch/cvan/KRICO/Pre/PISCES
"""

import argparse
import calendar
import datetime as dt
import os
import sys
from pathlib import Path

PRODUCT = "GLOBAL_MULTIYEAR_BGC_001_029"
DOI = "10.48670/moi-00019"
DATASET = "cmems_mod_glo_bgc_my_0.25deg_P1D-m"
VARIABLE = "chl"
N_LEVELS = 75
FIRST = (1993, 11)
LAST = (2025, 10)
WEST, EAST = -115.5, 40.5
SOUTH, NORTH = -80.0, -39.5


def _root():
    root = os.environ.get("KRICO_PISCES")
    if not root:
        sys.exit(
            "ERROR: KRICO_PISCES environment variable not set.\n"
            "Set it to the PISCES preprocessing directory, e.g.:\n"
            "  export KRICO_PISCES=/scratch/cvan/KRICO/Pre/PISCES"
        )
    return Path(root)


def months_of(year):
    """Months of one calendar year that fall inside the period."""
    return [m for m in range(1, 13) if FIRST <= (year, m) <= LAST]


def filename(year, month):
    return f"pisces_chl_{year}_{month:02d}.nc"


def dataset_version():
    """Label of the dataset version the toolbox selects by default."""
    import copernicusmarine
    catalogue = copernicusmarine.describe(dataset_id=DATASET,
                                          disable_progress_bar=True)
    return catalogue.products[0].datasets[0].versions[0].label


def inspect(path):
    """(n_days, n_levels) of a monthly file, or None if it cannot be read."""
    import xarray as xr
    try:
        with xr.open_dataset(path) as ds:
            return ds.sizes["time"], ds.sizes["depth"]
    except Exception:
        return None


def is_complete(path, year, month):
    expected = (calendar.monthrange(year, month)[1], N_LEVELS)
    return Path(path).is_file() and inspect(path) == expected


def subset(year, month, version, dry_run):
    import copernicusmarine
    n_days = calendar.monthrange(year, month)[1]
    return copernicusmarine.subset(
        dataset_id=DATASET,
        dataset_version=version,
        variables=[VARIABLE],
        minimum_longitude=WEST,
        maximum_longitude=EAST,
        minimum_latitude=SOUTH,
        maximum_latitude=NORTH,
        start_datetime=f"{year}-{month:02d}-01T00:00:00",
        end_datetime=f"{year}-{month:02d}-{n_days:02d}T23:59:59",
        output_directory=_root(),
        output_filename=filename(year, month),
        overwrite=True,
        dry_run=dry_run,
        disable_progress_bar=True,
    )


def download_year(year, months, dry_run):
    root = _root()
    root.mkdir(parents=True, exist_ok=True)
    version = dataset_version()
    print(f"{DATASET}, version {version} -> {root}")
    n_ok = n_skip = n_fail = 0
    total_mb = 0.0
    for month in months:
        name = filename(year, month)
        if dry_run:
            response = subset(year, month, version, dry_run=True)
            total_mb += response.file_size or 0.0
            print(f"  {name}: file {response.file_size:.0f} MB, "
                  f"transfer {response.data_transfer_size:.0f} MB, "
                  f"variables {response.variables}")
            if month == months[0]:
                for extent in response.coordinates_extent:
                    print(f"    {extent.model_dump()}")
            continue
        if is_complete(root / name, year, month):
            n_skip += 1
            continue
        try:
            subset(year, month, version, dry_run=False)
        except Exception as exc:
            print(f"  FAILED {name}: {exc}", file=sys.stderr)
            n_fail += 1
            continue
        if is_complete(root / name, year, month):
            n_ok += 1
            print(f"  got  {name}  {(root / name).stat().st_size / 1e9:.2f} GB")
        else:
            n_fail += 1
            print(f"  INCOMPLETE {name}: (days, levels) = {inspect(root / name)}",
                  file=sys.stderr)
    if dry_run:
        print(f"{year}: {len(months)} months, {total_mb / 1e3:.2f} GB estimated")
        return
    print(f"{year}: downloaded {n_ok}, already present {n_skip}, failed {n_fail}")
    if n_fail:
        sys.exit(1)


def write_provenance():
    import copernicusmarine
    root = _root()
    out = Path(__file__).resolve().parent / "pisces_chl_provenance.txt"
    expected = [(y, m) for y in range(FIRST[0], LAST[0] + 1) for m in months_of(y)]
    rows, missing = [], []
    for year, month in expected:
        path = root / filename(year, month)
        if is_complete(path, year, month):
            rows.append(f"{calendar.monthrange(year, month)[1]:>4}  "
                        f"{path.stat().st_size:>11}  {path.name}")
        else:
            missing.append(path.name)
    print(f"{len(rows)} complete files of {len(expected)} expected, "
          f"{len(missing)} missing or incomplete")
    if missing:
        print("  re-run the download for the years concerned before committing "
              "the provenance file", file=sys.stderr)
    accessed = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        "# Daily 3D chlorophyll inputs for the chlorophyll indices of the",
        "# recruitment dataset",
        "#",
        "# product    : Global Ocean Biogeochemistry Hindcast",
        f"#              {PRODUCT} (PISCES, no data assimilation)",
        f"# doi        : {DOI}",
        f"# dataset    : {DATASET}",
        f"# version    : {dataset_version()}",
        f"# toolbox    : copernicusmarine {copernicusmarine.__version__}",
        f"# variable   : {VARIABLE}, all {N_LEVELS} depth levels",
        f"# subset     : longitude {WEST} to {EAST}, latitude {SOUTH} to {NORTH},",
        f"#              daily, {FIRST[0]}-{FIRST[1]:02d} to {LAST[0]}-{LAST[1]:02d}",
        f"# accessed   : {accessed}",
        f"# n_files    : {len(rows)} (expected: {len(expected)})",
        "#",
        "# The files themselves live under $KRICO_PISCES and are not in the",
        "# repository. They are subsets written by the toolbox, so a second",
        "# download is not byte-identical; no checksums are given.",
        "#",
        "# Missing or incomplete:",
    ]
    lines += [f"#   {name}" for name in missing] or ["#   (none)"]
    lines += ["#", "# days        bytes  file (relative to $KRICO_PISCES)"]
    out.write_text("\n".join(lines + rows) + "\n")
    print(f"Provenance written to {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("year", nargs="?", type=int,
                        help="calendar year to download (monthly files)")
    parser.add_argument("--month", type=int,
                        help="download this month only")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the estimated sizes, download nothing")
    parser.add_argument("--provenance", action="store_true",
                        help="write pisces_chl_provenance.txt from what is on disk")
    args = parser.parse_args()

    if args.provenance:
        write_provenance()
    elif args.year is not None:
        months = months_of(args.year)
        if args.month is not None:
            months = [m for m in months if m == args.month]
        if not months:
            sys.exit(f"nothing to download: the period is {FIRST[0]}-{FIRST[1]:02d} "
                     f"to {LAST[0]}-{LAST[1]:02d}")
        download_year(args.year, months, args.dry_run)
    else:
        parser.error("give a year or --provenance")


if __name__ == "__main__":
    main()
