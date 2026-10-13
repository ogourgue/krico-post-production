"""
Download the daily southern-hemisphere files of the NOAA/NSIDC Climate Data
Record of Passive Microwave Sea Ice Concentration, Version 6 (G02202 v6),
for the period the KRICO recruitment classification needs.

Usage
-----
    python3 download_nsidc.py YEAR          # one year, into $KRICO_NSIDC/daily/YEAR/
    python3 download_nsidc.py --provenance  # after all years: checksums + missing dates

Period: 1993-10-01 to 2025-10-31. The first cohort is released on
15 November 1993 and M1 reads back 24 days to 22 October 1993; the last
cohort is released on 14 March 2025 and tracked 200 days, to 30 September
2025. One month of margin on either side.

Layout, mirroring the GLORYS12 preprocessing directory:
    $KRICO_NSIDC/ancillary/G02202-ancillary-pss25-v06r00.nc   grid + masks
    $KRICO_NSIDC/daily/YYYY/sic_pss25_YYYYMMDD_<sat>_v06r00.nc raw daily files

The daily archive is organised by year, and the satellite tag <sat> (F11,
F13, F17, am2, ...) changes through the record, so the server listing is
matched on the date stamp rather than a fixed file name. Files are ~0.35 MB
each, ~11,700 in total (~4 GB).

The download is resumable: a file already present and non-empty is skipped,
and a transfer in progress writes to a temporary name that is renamed only on
success, so an interrupted job never leaves a truncated file that a re-run
would skip. Re-run until --provenance reports every listed file present.

The CDR has days with no file at all (sensor outages). --provenance lists
them in nsidc_daily_provenance.txt (next to this script, tracked in git)
together with the DOI, source URL, access date and sha256 of every file, so
the inputs stay verifiable after NSIDC retires the version 6 distribution
path. Filling those days is the job of the preprocessing step, not of this
script. Model: krico-paper1/F1_domain_map/download_sic.sh.

Access is anonymous HTTPS; no Earthdata account is required.

Environment
-----------
    KRICO_NSIDC      required; e.g. /scratch/cvan/KRICO/Pre/NSIDC
    NSIDC_BASE_URL   optional; alternative archive root, used to test
                     against a mirror
"""

import argparse
import hashlib
from concurrent.futures import ThreadPoolExecutor
import os
import re
import ssl
import sys
import urllib.error
import urllib.request
from datetime import date, timedelta
from pathlib import Path

BASE_URL = os.environ.get(
    "NSIDC_BASE_URL", "https://noaadata.apps.nsidc.org/NOAA/G02202_V6"
)
DOI = "10.7265/b18j-z797"
DATE_MIN = date(1993, 10, 1)
DATE_MAX = date(2025, 10, 31)
ANCILLARY = "G02202-ancillary-pss25-v06r00.nc"
FILE_RE = re.compile(r"sic_pss25_(\d{8})_[A-Za-z0-9]+_v06r00\.nc")
RETRIES = 5

# NSIDC documents --no-check-certificate for this server. Verification is
# attempted first and only relaxed if it fails, so the fallback is visible.
_SSL_CONTEXTS = [ssl.create_default_context(), ssl._create_unverified_context()]


def _root():
    root = os.environ.get("KRICO_NSIDC")
    if not root:
        sys.exit(
            "ERROR: KRICO_NSIDC environment variable not set.\n"
            "Set it to the NSIDC preprocessing directory, e.g.:\n"
            "  export KRICO_NSIDC=/scratch/cvan/KRICO/Pre/NSIDC"
        )
    return Path(root)


def _urlopen(url):
    """Open a URL, trying a verified SSL context before an unverified one."""
    last = None
    for ctx in _SSL_CONTEXTS:
        try:
            return urllib.request.urlopen(url, timeout=60, context=ctx)
        except urllib.error.URLError as exc:
            last = exc
            if not isinstance(exc.reason, ssl.SSLError):
                raise
    raise last


def list_year(year):
    """Names of the daily files for one year that fall inside the period."""
    with _urlopen(f"{BASE_URL}/south/daily/{year}/") as resp:
        index = resp.read().decode("utf-8", errors="replace")
    names = set()
    for match in FILE_RE.finditer(index):
        stamp = match.group(1)
        d = date(int(stamp[:4]), int(stamp[4:6]), int(stamp[6:8]))
        if DATE_MIN <= d <= DATE_MAX:
            names.add(match.group(0))
    return sorted(names)


def fetch(url, dest):
    """Download url to dest unless dest already exists. Returns True on success."""
    dest = Path(dest)
    if dest.is_file() and dest.stat().st_size > 0:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f"{dest.name}.part.{os.getpid()}")
    for attempt in range(1, RETRIES + 1):
        try:
            with _urlopen(url) as resp, open(tmp, "wb") as fh:
                while chunk := resp.read(1 << 20):
                    fh.write(chunk)
            if tmp.stat().st_size == 0:
                raise OSError("empty response")
            os.replace(tmp, dest)
            return True
        except (urllib.error.URLError, OSError) as exc:
            tmp.unlink(missing_ok=True)
            print(f"  attempt {attempt}/{RETRIES} failed for {dest.name}: {exc}",
                  file=sys.stderr)
    return False


def download_year(year):
    root = _root()
    print(f"Ancillary file -> {root / 'ancillary'}")
    if not fetch(f"{BASE_URL}/ancillary/{ANCILLARY}", root / "ancillary" / ANCILLARY):
        print("  WARNING: ancillary file not retrieved", file=sys.stderr)

    try:
        names = list_year(year)
    except (urllib.error.URLError, OSError) as exc:
        sys.exit(f"ERROR: could not list {BASE_URL}/south/daily/{year}/: {exc}")
    print(f"{year}: {len(names)} files listed on server within "
          f"{DATE_MIN} to {DATE_MAX}")
    n_ok = n_skip = n_fail = 0
    for name in names:
        dest = root / "daily" / str(year) / name
        if dest.is_file() and dest.stat().st_size > 0:
            n_skip += 1
            continue
        if fetch(f"{BASE_URL}/south/daily/{year}/{name}", dest):
            n_ok += 1
            print(f"  got  {name}")
        else:
            n_fail += 1
    print(f"{year}: downloaded {n_ok}, already present {n_skip}, failed {n_fail}")
    if n_fail:
        sys.exit(1)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def write_provenance():
    root = _root()
    out = Path(__file__).resolve().parent / "nsidc_daily_provenance.txt"

    # Dates present on disk, from file names.
    present = {}
    for path in sorted((root / "daily").rglob("sic_pss25_*_v06r00.nc")):
        m = FILE_RE.fullmatch(path.name)
        if m:
            stamp = m.group(1)
            present[date(int(stamp[:4]), int(stamp[4:6]), int(stamp[6:8]))] = path

    # Dates the server offers, so that "missing on disk" and "missing in the
    # archive" are told apart.
    def _listed_dates(year):
        try:
            return [FILE_RE.fullmatch(n).group(1) for n in list_year(year)]
        except (urllib.error.URLError, OSError) as exc:
            print(f"  WARNING: could not list {year}: {exc}", file=sys.stderr)
            return []

    years = range(DATE_MIN.year, DATE_MAX.year + 1)
    listed = set()
    with ThreadPoolExecutor(max_workers=8) as pool:
        for stamps in pool.map(_listed_dates, years):
            for stamp in stamps:
                listed.add(date(int(stamp[:4]), int(stamp[4:6]), int(stamp[6:8])))

    n_days = (DATE_MAX - DATE_MIN).days + 1
    all_days = [DATE_MIN + timedelta(days=i) for i in range(n_days)]
    no_file = [d for d in all_days if d not in listed]
    not_downloaded = sorted(listed - set(present))

    print(f"{len(present)} files on disk, {len(listed)} listed on server, "
          f"{len(no_file)} dates without a file in the archive, "
          f"{len(not_downloaded)} listed files not downloaded")
    if not_downloaded:
        print("  re-run the download for the years concerned before committing "
              "the provenance file", file=sys.stderr)

    import datetime as _dt
    accessed = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    lines = [
        "# Daily sea-ice concentration inputs for the NSIDC recruitment classification",
        "#",
        "# dataset    : NOAA/NSIDC Climate Data Record of Passive Microwave",
        "#              Sea Ice Concentration, G02202 Version 6",
        f"# doi        : {DOI}",
        f"# source_url : {BASE_URL}/south/daily/YYYY/",
        f"# ancillary  : {BASE_URL}/ancillary/{ANCILLARY}",
        f"# subset     : southern hemisphere, daily, {DATE_MIN} to {DATE_MAX}",
        f"# accessed   : {accessed}",
        f"# n_files    : {len(present)} (listed on server: {len(listed)})",
        f"# n_days_without_file : {len(no_file)}",
        f"# n_listed_not_downloaded : {len(not_downloaded)}",
        "#",
        "# The files themselves live under $KRICO_NSIDC/daily/YYYY/ and are not",
        "# in the repository. Re-download with download_nsidc.py and verify",
        "# against the checksums below.",
        "#",
        "# Dates in the period for which the archive offers no file:",
    ]
    lines += [f"#   {d}" for d in no_file] or ["#   (none)"]
    if not_downloaded:
        lines += ["#", "# Listed on the server but not on disk:"]
        lines += [f"#   {d}" for d in not_downloaded]
    lines += ["#", "# sha256  path (relative to $KRICO_NSIDC)"]

    # Hashing is bound by per-file open latency on a parallel file system, not
    # by CPU, so it is spread over threads.
    paths = []
    anc = root / "ancillary" / ANCILLARY
    if anc.is_file():
        paths.append(anc)
    paths += [present[d] for d in sorted(present)]
    print(f"hashing {len(paths)} files ...")
    with ThreadPoolExecutor(max_workers=16) as pool:
        for path, digest in zip(paths, pool.map(_sha256, paths)):
            lines.append(f"{digest}  {path.relative_to(root)}")

    out.write_text("\n".join(lines) + "\n")
    print(f"Provenance written to {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("year", nargs="?", type=int,
                        help="year to download (daily files)")
    parser.add_argument("--provenance", action="store_true",
                        help="write nsidc_daily_provenance.txt from what is on disk")
    args = parser.parse_args()

    if args.provenance:
        write_provenance()
    elif args.year is not None:
        if not DATE_MIN.year <= args.year <= DATE_MAX.year:
            sys.exit(f"year must be in {DATE_MIN.year}-{DATE_MAX.year}")
        download_year(args.year)
    else:
        parser.error("give a year or --provenance")


if __name__ == "__main__":
    main()
