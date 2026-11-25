"""
Check the inputs of pass 1 (process_cohort.py) for whole spawning years,
before any job is submitted.

Usage
-----
    python3 check_inputs.py            # the 32 spawning years
    python3 check_inputs.py 1994       # one spawning year
    python3 check_inputs.py 1994 quick # headers and release times only
    python3 check_inputs.py 1994 nsidc # for a run with another ice source

Every cohort goes through krico_recruitment.inputs.check_cohort, the
function process_cohort.py calls before it computes anything, so a cohort
accepted here is accepted by the job:

  - the file is named YYYY_MM_DD.nc and holds the six trajectory variables
    and time on (trajectory, obs), with 201 daily values;
  - every particle is released at 12:00 on the day of the file name;
  - the release-cell file of the spawning year exists, has one record per
    day at 12:00, levels down to the deepest hatching depth, and holds the
    spawning day (release minus 24 days) and the release day;
  - every particle has a position at release, and the grid point nearest
    to it is a cell of the release-cell file (skipped with "quick");
  - with nsidc or glorys12_25km: the monthly files of that ice source hold
    every day from the spawning day to the last tracked day.

The check does not stop at the first problem: every cohort is tested and
every problem is listed. It also compares the cohorts found with the
release days expected for the year (15 November to 14 March); a missing or
unexpected day is reported as a warning, since the job classifies the
cohorts it finds.

Reading the release positions takes a few seconds per cohort: run whole
years through check_inputs.sh rather than on a login node.

Environment
-----------
    KRICO_RUNS                     required; KRICO_NNNN folders of trajectories
    KRICO_GLORYS12_RELEASE_CELLS   required
    KRICO_NSIDC, KRICO_GLORYS12_25KM   monthly files of the other ice sources
"""

import os
import sys
import time
from datetime import date, timedelta
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from krico_recruitment import inputs, nsidc_ice, release_cells  # noqa: E402


def folders(year):
    """The four simulation folders of a spawning year (as run_year.sh)."""
    first = 4 * (year - release_cells.FIRST_YEAR) + 1
    return [f"KRICO_{first + k:04d}" for k in range(4)]


def expected_days(year):
    first = date(year - 1, *release_cells.FIRST_RELEASE)
    last = date(year, *release_cells.LAST_RELEASE)
    return {first + timedelta(days=k) for k in range((last - first).days + 1)}


def check_year(year, runs_dir, cells_dir, read_positions,
               ice_source="glorys12", ice_dir=None):
    started = time.time()
    problems, warnings = [], []
    files = []
    for folder in folders(year):
        path = runs_dir / folder
        if not path.is_dir():
            problems.append(f"folder {path} is missing")
            continue
        files += sorted(path.glob("[0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9].nc"))
    print(f"Spawning year {year}: {len(files)} cohorts in {', '.join(folders(year))}"
          + ("" if ice_source == "glorys12" else
             f"; sea ice from {ice_source} in {ice_dir}"), flush=True)

    found = {}
    for path in files:
        day = inputs.cohort_day(path)
        if day in found:
            problems.append(f"{path.name} is in {found[day].parent.name} and "
                            f"in {path.parent.name}")
        found[day] = path
    expected = expected_days(year)
    absent = sorted(expected - set(found))
    if absent:
        warnings.append(f"{len(absent)} of the {len(expected)} release days of "
                        f"the year have no cohort ({absent[0]} to {absent[-1]}"
                        + (", not all of them" if len(absent) <
                           (absent[-1] - absent[0]).days + 1 else "") + ")")
    for day in sorted(set(found) - expected):
        warnings.append(f"{found[day].parent.name}/{found[day].name} is outside "
                        "the release season of the year")

    statics = {}
    n_particles = {}
    in_release = in_neighbour = accepted = 0
    is_release = None
    for k, path in enumerate(files, 1):
        try:
            checked = inputs.check_cohort(path, cells_dir, statics=statics,
                                          read_positions=read_positions,
                                          ice_source=ice_source, ice_dir=ice_dir)
        except (ValueError, FileNotFoundError, OSError) as exc:
            problems.append(str(exc))
            continue
        accepted += 1
        n = checked["n_particles"]
        n_particles[n] = n_particles.get(n, 0) + 1
        if read_positions:
            if is_release is None:
                with netCDF4.Dataset(checked["static"]["path"]) as nc:
                    is_release = np.asarray(nc["is_release_cell"][:]) == 1
            inside = int(is_release[checked["cell"]].sum())
            in_release += inside
            in_neighbour += n - inside
        if k % 20 == 0:
            print(f"  {k} cohorts, {time.time() - started:.0f} s", flush=True)

    print(f"  accepted: {accepted} of {len(files)} cohorts; particles per "
          f"cohort: " + ", ".join(f"{n} ({c} cohorts)"
                                  for n, c in sorted(n_particles.items())))
    if read_positions and accepted:
        print(f"  nearest grid point: {in_release} particles in a release "
              f"cell, {in_neighbour} in a neighbour, none outside")
    for static in statics.values():
        print(f"  {static['path'].name}: {static['n_time']} days from "
              f"{static['first_day']}, {static['n_cell']} cells, "
              f"code {static['code_version']}")
    for line in warnings:
        print(f"  WARNING: {line}")
    for line in problems:
        print(f"  PROBLEM: {line}")
    print(f"  verdict {year}: " + ("PASS" if not problems else
                                   f"FAIL, {len(problems)} problem(s)")
          + (f", {len(warnings)} warning(s)" if warnings else "")
          + f"  ({time.time() - started:.0f} s)", flush=True)
    return not problems


def main():
    args = sys.argv[1:]
    if any(a in ("-h", "--help") for a in args):
        print(__doc__)
        return
    quick = "quick" in args
    sources = [a for a in args if a in nsidc_ice.SOURCES or a == "glorys12"]
    if len(sources) > 1:
        sys.exit("give one ice source at most")
    ice_source = sources[0] if sources else "glorys12"
    try:
        years = [int(a) for a in args if a != "quick" and a not in sources]
    except ValueError:
        sys.exit("arguments: years, 'quick', and an ice source among glorys12, "
                 + ", ".join(nsidc_ice.SOURCES))
    for year in years:
        if not release_cells.FIRST_YEAR <= year <= release_cells.LAST_YEAR:
            sys.exit(f"year must be in {release_cells.FIRST_YEAR}-"
                     f"{release_cells.LAST_YEAR}")
    if not years:
        years = list(range(release_cells.FIRST_YEAR, release_cells.LAST_YEAR + 1))
    directories = []
    names = ["KRICO_RUNS", "KRICO_GLORYS12_RELEASE_CELLS"]
    if ice_source != "glorys12":
        names.append(nsidc_ice.SOURCES[ice_source]["env"])
    for name in names:
        value = os.environ.get(name)
        if not value:
            sys.exit(f"ERROR: {name} environment variable not set.")
        directories.append(Path(value))
    ice_dir = directories[2] if ice_source != "glorys12" else None
    ok = True
    for year in years:
        ok &= check_year(year, directories[0], directories[1], not quick,
                         ice_source=ice_source, ice_dir=ice_dir)
    print("Check passed." if ok else "CHECK FAILED.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
