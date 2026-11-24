"""
Download the daily surface chlorophyll climatology of the Copernicus-GlobColour
satellite product (OCEANCOLOUR_GLO_BGC_L4_MY_009_104, multi-sensor, 4 km) for
the domain of the KRICO trajectories.

Usage
-----
    python3 download_globcolour.py              # one file, into $KRICO_GLOBCOLOUR
    python3 download_globcolour.py --dry-run    # sizes only, nothing is downloaded
    python3 download_globcolour.py --provenance # after the download: what is on disk

The dataset is a day-of-year climatology: 366 daily fields, carried on the
calendar of the leap year 2008 (2008-01-01 to 2008-12-31). It is sampled by
day of year, so there is no period to choose.

Variables: of the eleven the dataset offers, the mean and the median of the
climatology and the number of values behind them (CHL_mean, CHL_median,
CHL_count). Edit VARIABLES to take others (CHL_standard_deviation, CHL_minimum,
CHL_maximum, CHL_percentile_3, CHL_percentile_97,
CHL_percentile_standard_deviation, flags); each is about 5 GB.

Domain: the GLORYS12 box of the trajectories (115W to 40E, south of 40S) with
a margin of 0.5 degree, as in download_pisces.py.

Layout, mirroring the other preprocessing directories:
    $KRICO_GLOBCOLOUR/globcolour_chl_climatology.nc   one file, 366 days

The download is resumable in the sense that a complete file is not downloaded
again; the toolbox writes to a temporary name that is renamed only on success,
so an interrupted job never leaves a truncated file under the final name.

As for PISCES, the file is a subset written on this side by the toolbox, so
checksums would verify nothing. --provenance records the product, DOI, dataset
version, toolbox version, subset and access date, with the size of the file,
in globcolour_chl_provenance.txt (next to this script, tracked in git).

Access needs the Copernicus Marine toolbox, version 2 (provided by the python3
module of the ECMWF HPC; elsewhere pip install copernicusmarine, Python >=
3.10), and a free Copernicus Marine account, stored once on the machine with:
    copernicusmarine login

Environment
-----------
    KRICO_GLOBCOLOUR   required; e.g. /scratch/cvan/KRICO/Pre/GLOBCOLOUR
"""

import argparse
import datetime as dt
import os
import sys
from pathlib import Path

PRODUCT = "OCEANCOLOUR_GLO_BGC_L4_MY_009_104"
DOI = "10.48670/moi-00281"
DATASET = "cmems_obs-oc_glo_bgc-plankton_my_l4-multi-climatology-4km_P1D"
VARIABLES = ["CHL_mean", "CHL_median", "CHL_count"]
N_DAYS = 366
START, END = "2008-01-01T00:00:00", "2008-12-31T23:59:59"
WEST, EAST = -115.5, 40.5
SOUTH, NORTH = -80.0, -39.5
FILENAME = "globcolour_chl_climatology.nc"


def _root():
    root = os.environ.get("KRICO_GLOBCOLOUR")
    if not root:
        sys.exit(
            "ERROR: KRICO_GLOBCOLOUR environment variable not set.\n"
            "Set it to the GlobColour preprocessing directory, e.g.:\n"
            "  export KRICO_GLOBCOLOUR=/scratch/cvan/KRICO/Pre/GLOBCOLOUR"
        )
    return Path(root)


def dataset_version():
    """Label of the dataset version the toolbox selects by default."""
    import copernicusmarine
    catalogue = copernicusmarine.describe(dataset_id=DATASET,
                                          disable_progress_bar=True)
    return catalogue.products[0].datasets[0].versions[0].label


def inspect(path):
    """(n_days, variables present) of the file, or None if it cannot be read."""
    import xarray as xr
    try:
        with xr.open_dataset(path) as ds:
            return ds.sizes["time"], sorted(v for v in VARIABLES if v in ds)
    except Exception:
        return None


def is_complete(path):
    return Path(path).is_file() and inspect(path) == (N_DAYS, sorted(VARIABLES))


def subset(version, dry_run):
    import copernicusmarine
    return copernicusmarine.subset(
        dataset_id=DATASET,
        dataset_version=version,
        variables=VARIABLES,
        minimum_longitude=WEST,
        maximum_longitude=EAST,
        minimum_latitude=SOUTH,
        maximum_latitude=NORTH,
        start_datetime=START,
        end_datetime=END,
        output_directory=_root(),
        output_filename=FILENAME,
        overwrite=True,
        dry_run=dry_run,
        disable_progress_bar=True,
    )


def download(dry_run):
    root = _root()
    root.mkdir(parents=True, exist_ok=True)
    path = root / FILENAME
    version = dataset_version()
    print(f"{DATASET}, version {version} -> {root}")
    if dry_run:
        response = subset(version, dry_run=True)
        print(f"  {FILENAME}: file {response.file_size / 1e3:.2f} GB, "
              f"transfer {response.data_transfer_size / 1e3:.2f} GB, "
              f"variables {response.variables}")
        for extent in response.coordinates_extent:
            print(f"    {extent.model_dump()}")
        return
    if is_complete(path):
        print(f"  already present: {FILENAME}")
        return
    subset(version, dry_run=False)
    if not is_complete(path):
        sys.exit(f"ERROR: {FILENAME} is incomplete: "
                 f"(days, variables) = {inspect(path)}")
    print(f"  got  {FILENAME}  {path.stat().st_size / 1e9:.2f} GB")


def write_provenance():
    import copernicusmarine
    import xarray as xr
    root = _root()
    path = root / FILENAME
    out = Path(__file__).resolve().parent / "globcolour_chl_provenance.txt"
    if not is_complete(path):
        sys.exit(f"ERROR: {path} is missing or incomplete; run the download first")
    with xr.open_dataset(path) as ds:
        grid = (f"{ds.sizes['latitude']} x {ds.sizes['longitude']} cells, "
                f"longitude {float(ds.longitude.min()):.4f} to "
                f"{float(ds.longitude.max()):.4f}, latitude "
                f"{float(ds.latitude.min()):.4f} to {float(ds.latitude.max()):.4f}")
    accessed = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        "# Daily surface chlorophyll climatology for the chlorophyll indices of",
        "# the recruitment dataset",
        "#",
        "# product    : Global Ocean Colour (Copernicus-GlobColour),",
        "#              Bio-Geo-Chemical, L4, from Satellite Observations",
        f"#              {PRODUCT}",
        f"# doi        : {DOI}",
        f"# dataset    : {DATASET}",
        f"# version    : {dataset_version()}",
        f"# toolbox    : copernicusmarine {copernicusmarine.__version__}",
        f"# variables  : {', '.join(VARIABLES)}",
        f"# subset     : requested longitude {WEST} to {EAST}, latitude {SOUTH}",
        f"#              to {NORTH}; day-of-year climatology of {N_DAYS} days,",
        "#              carried on the calendar of 2008",
        f"# grid       : {grid}",
        f"# accessed   : {accessed}",
        "#",
        "# The file itself lives under $KRICO_GLOBCOLOUR and is not in the",
        "# repository. It is a subset written by the toolbox, so a second",
        "# download is not byte-identical; no checksum is given.",
        "#",
        "# days        bytes  file (relative to $KRICO_GLOBCOLOUR)",
        f"{N_DAYS:>4}  {path.stat().st_size:>11}  {FILENAME}",
    ]
    out.write_text("\n".join(lines) + "\n")
    print(f"Provenance written to {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--dry-run", action="store_true",
                        help="print the estimated sizes, download nothing")
    parser.add_argument("--provenance", action="store_true",
                        help="write globcolour_chl_provenance.txt from what is on disk")
    args = parser.parse_args()

    if args.provenance:
        write_provenance()
    else:
        download(args.dry_run)


if __name__ == "__main__":
    main()
