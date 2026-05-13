"""
download_from_zenodo.py
========================

Download the KRICO recruitment outcomes dataset from Zenodo into
recruitment/archives/.

Fetches the file list via the Zenodo REST API, then downloads each
recruitment_YYYY.tar.gz archive with MD5 verification against the
Zenodo metadata.

Behavior:
  - Idempotent. Re-running skips archives that already exist with the
    correct size and MD5 checksum.
  - Resumable. Partial downloads are picked up using HTTP Range requests.
  - Retries. Network errors and MD5 mismatches trigger up to MAX_RETRIES
    attempts per archive before failing.
  - Selective. The --years argument restricts the download to a subset.

The default DOI points to v1.0.0 of the dataset (10.5281/zenodo.20101159).
Override with --doi for a different version.

Usage
-----
    cd recruitment/scripts
    python download_from_zenodo.py                          # All 32 archives
    python download_from_zenodo.py --years 1994             # Single year
    python download_from_zenodo.py --years 1994-1998        # Year range
    python download_from_zenodo.py --years 1994,1996,2000   # Specific years
    python download_from_zenodo.py --doi 10.5281/zenodo.X   # Override DOI

Requires only the Python standard library (Python >= 3.11).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Default Zenodo DOI for the KRICO recruitment dataset (v1.0.0).
# Update when releasing a new dataset version, then re-tag this repo so the
# new code-data pairing is captured. Override at runtime with --doi.
DEFAULT_DOI = "10.5281/zenodo.20101159"

# Filename pattern for recruitment archives: recruitment_YYYY.tar.gz.
ARCHIVE_PATTERN = re.compile(r"^recruitment_(\d{4})\.tar\.gz$")

# Download buffer (1 MB) and progress refresh rate (1 s).
CHUNK_SIZE = 1024 * 1024
PROGRESS_INTERVAL_S = 1.0

# Retry policy: per-archive, applies to both network errors and MD5 failures.
MAX_RETRIES = 3
RETRY_BACKOFF_S = 5

# Timeouts.
METADATA_TIMEOUT_S = 30
DOWNLOAD_TIMEOUT_S = 60


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_record_id(doi_or_url: str) -> str:
    """Extract the Zenodo record ID from a DOI string or Zenodo URL."""
    m = re.search(r"zenodo\.(\d+)", doi_or_url)
    if not m:
        raise ValueError(
            f"Could not parse Zenodo record ID from: {doi_or_url!r}\n"
            f"Expected format like '10.5281/zenodo.XXXXXXX' "
            f"or 'https://zenodo.org/records/XXXXXXX'."
        )
    return m.group(1)


def parse_years(years_arg: str) -> set[int]:
    """
    Parse --years into a set of integers.

    Accepts a comma-separated list of single years and inclusive ranges,
    e.g. '1994', '1994-1998', '1994,1996,2000', '1994-1996,2000'.
    """
    years: set[int] = set()
    for chunk in years_arg.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" in chunk:
            start, end = chunk.split("-", 1)
            years.update(range(int(start), int(end) + 1))
        else:
            years.add(int(chunk))
    return years


def format_bytes(n: float) -> str:
    """Format a byte count as a human-readable string (decimal units, matching Zenodo)."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1000:
            return f"{n:.2f} {unit}"
        n /= 1000
    return f"{n:.2f} PB"


def file_md5(path: Path) -> str:
    """Compute MD5 hash of a file."""
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(CHUNK_SIZE), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Zenodo metadata
# ---------------------------------------------------------------------------

def fetch_record_metadata(record_id: str) -> dict:
    """Fetch Zenodo record metadata via the REST API."""
    url = f"https://zenodo.org/api/records/{record_id}"
    print(f"Fetching record metadata: {url}")
    try:
        with urllib.request.urlopen(url, timeout=METADATA_TIMEOUT_S) as response:
            return json.load(response)
    except urllib.error.HTTPError as e:
        raise RuntimeError(
            f"HTTP {e.code} fetching Zenodo record {record_id}: {e.reason}"
        ) from e
    except urllib.error.URLError as e:
        raise RuntimeError(
            f"Network error fetching Zenodo record {record_id}: {e.reason}\n"
            f"Check internet connection."
        ) from e


def collect_archive_files(
    metadata: dict, record_id: str, years_filter: set[int] | None,
) -> list[dict]:
    """
    Filter the record's files to recruitment_YYYY.tar.gz archives matching
    years_filter (or all years if filter is None).

    Returns a list of dicts: {year, key, size, md5, url}, sorted by year.
    """
    files = metadata.get("files", [])
    if not files:
        raise RuntimeError("Zenodo record contains no files.")

    archives = []
    for f in files:
        key = f["key"]
        m = ARCHIVE_PATTERN.match(key)
        if not m:
            print(f"  [skip] {key} (not a recruitment archive)")
            continue
        year = int(m.group(1))
        if years_filter is not None and year not in years_filter:
            continue
        archives.append({
            "year": year,
            "key": key,
            "size": f["size"],
            "md5": f["checksum"].replace("md5:", ""),
            # Use the public file URL pattern. It works without API auth and
            # is stable across Zenodo's API revisions.
            "url": f"https://zenodo.org/records/{record_id}/files/{key}?download=1",
        })

    archives.sort(key=lambda a: a["year"])
    return archives


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def _stream_download(url: str, dest: Path, expected_size: int) -> None:
    """
    Start or resume a download to dest. No checksum or retry logic — the
    caller handles those.
    """
    start_byte = dest.stat().st_size if dest.exists() else 0
    if start_byte > expected_size:
        # Partial download larger than expected — restart cleanly.
        dest.unlink()
        start_byte = 0
    if start_byte == expected_size:
        # Bytes already on disk; caller verifies MD5.
        return

    req = urllib.request.Request(url)
    if start_byte > 0:
        req.add_header("Range", f"bytes={start_byte}-")
        action = "resuming"
    else:
        action = "downloading"

    print(f"  [{action}] {dest.name} "
          f"({format_bytes(start_byte)} / {format_bytes(expected_size)})")

    mode = "ab" if start_byte > 0 else "wb"
    bytes_so_far = start_byte
    last_print = time.monotonic()

    with urllib.request.urlopen(req, timeout=DOWNLOAD_TIMEOUT_S) as response, \
         dest.open(mode) as out:
        while True:
            chunk = response.read(CHUNK_SIZE)
            if not chunk:
                break
            out.write(chunk)
            bytes_so_far += len(chunk)

            now = time.monotonic()
            if now - last_print >= PROGRESS_INTERVAL_S:
                pct = 100.0 * bytes_so_far / expected_size
                print(f"    {format_bytes(bytes_so_far)} / "
                      f"{format_bytes(expected_size)} ({pct:5.1f}%)",
                      end="\r", flush=True)
                last_print = now

    # Newline so the next message doesn't overwrite the last progress line.
    print()


def download_archive(archive: dict, dest: Path) -> None:
    """
    Download one archive with resume, retry, and MD5 verification.

    Skips immediately if dest already exists with matching size and MD5.
    Raises RuntimeError if the file cannot be downloaded successfully
    within MAX_RETRIES attempts.
    """
    expected_size = archive["size"]
    expected_md5 = archive["md5"]

    # Fast path: already complete and verified.
    if dest.exists() and dest.stat().st_size == expected_size:
        print(f"  [check] {dest.name} present, verifying MD5...")
        if file_md5(dest) == expected_md5:
            print(f"  [skip]  {dest.name} (complete, MD5 OK)")
            return
        print(f"  [warn]  {dest.name} size matches but MD5 differs; redownloading")
        dest.unlink()

    last_error: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            _stream_download(archive["url"], dest, expected_size)

            if dest.stat().st_size != expected_size:
                raise ValueError(
                    f"size mismatch: expected {expected_size:,}, "
                    f"got {dest.stat().st_size:,}"
                )

            print(f"  [verify] computing MD5 for {dest.name}...")
            actual_md5 = file_md5(dest)
            if actual_md5 != expected_md5:
                raise ValueError(
                    f"MD5 mismatch: expected {expected_md5}, got {actual_md5}"
                )

            print(f"  [OK]    {dest.name} downloaded and verified")
            return

        except (urllib.error.URLError, ConnectionResetError, TimeoutError,
                ValueError, OSError) as e:
            last_error = e
            if attempt < MAX_RETRIES:
                wait = RETRY_BACKOFF_S * attempt
                print(f"  [retry] attempt {attempt} failed ({e}); "
                      f"waiting {wait}s before retry {attempt + 1}/{MAX_RETRIES}")
                # Remove the failed partial file before retrying when the
                # failure was a content-level mismatch (size or MD5).
                if isinstance(e, ValueError) and dest.exists():
                    dest.unlink()
                time.sleep(wait)
            else:
                print(f"  [FAIL]  {dest.name} after {MAX_RETRIES} attempts: {e}")

    raise RuntimeError(
        f"Failed to download {dest.name} after {MAX_RETRIES} attempts: {last_error}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download the KRICO recruitment outcomes dataset from Zenodo "
            "into recruitment/archives/."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python download_from_zenodo.py                          # All 32 archives\n"
            "  python download_from_zenodo.py --years 1994             # Single year\n"
            "  python download_from_zenodo.py --years 1994-1998        # Year range\n"
            "  python download_from_zenodo.py --years 1994,1996,2000   # Specific years\n"
            "  python download_from_zenodo.py --doi 10.5281/zenodo.X   # Override DOI\n"
        ),
    )
    parser.add_argument(
        "--doi",
        default=DEFAULT_DOI,
        help=f"Zenodo DOI to download from (default: {DEFAULT_DOI})",
    )
    parser.add_argument(
        "--years",
        default=None,
        help="Subset of years to download. e.g. '1994', '1994-1998', '1994,1996'.",
    )
    args = parser.parse_args()

    # Resolve output directory: recruitment/archives/ relative to this script.
    # Mirrors the layout used by archive_by_year.sh and extract_archives.sh.
    script_dir = Path(__file__).resolve().parent
    recruitment_dir = script_dir.parent
    archives_dir = recruitment_dir / "archives"
    archives_dir.mkdir(parents=True, exist_ok=True)

    # Parse year filter.
    years_filter = parse_years(args.years) if args.years else None

    # Fetch and summarize record metadata.
    record_id = parse_record_id(args.doi)
    metadata = fetch_record_metadata(record_id)
    record_meta = metadata.get("metadata", {})

    print()
    print("=" * 72)
    print(f"Zenodo record: {record_meta.get('title', '(no title)')}")
    print(f"Version:       {record_meta.get('version', '?')}")
    print(f"Record ID:     {record_id}")
    print(f"DOI:           {record_meta.get('doi', args.doi)}")
    print("=" * 72)
    print()

    # Filter to the archives we want.
    archives = collect_archive_files(metadata, record_id, years_filter)
    if not archives:
        print("No matching archives to download.", file=sys.stderr)
        sys.exit(1)

    total_size = sum(a["size"] for a in archives)
    print(f"To process: {len(archives)} archives, {format_bytes(total_size)} total")
    print(f"Destination: {archives_dir}")
    print()

    # Download.
    n_ok = 0
    n_fail = 0
    failed_keys: list[str] = []

    for archive in archives:
        dest = archives_dir / archive["key"]
        print(f"--- {archive['key']} ({format_bytes(archive['size'])}) ---")
        try:
            download_archive(archive, dest)
            n_ok += 1
        except RuntimeError as e:
            print(f"  [ERROR] {e}", file=sys.stderr)
            n_fail += 1
            failed_keys.append(archive["key"])

    # Summary.
    print()
    print("=" * 72)
    print(f"Done: {n_ok} archives in {archives_dir}, {n_fail} failed")
    print("=" * 72)
    print()

    if n_fail > 0:
        print("Failed archives (re-run the script to retry):", file=sys.stderr)
        for k in failed_keys:
            print(f"  {k}", file=sys.stderr)
        sys.exit(1)

    print("Next step: extract archives into data/ with extract_archives.sh:")
    print()
    print(f"  cd {script_dir}")
    print(f"  ./extract_archives.sh")


if __name__ == "__main__":
    main()