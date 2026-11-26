"""
Check the inputs of pass 2 (process_cohort_chl.py) for whole spawning
years, before any job is submitted.

Usage
-----
    python3 check_inputs_chl.py            # the 32 spawning years
    python3 check_inputs_chl.py 1994       # one spawning year

Every cohort goes through the checks process_cohort_chl.py itself makes
before it computes anything (its function check), so a cohort accepted
here is accepted by the job:

  - the file is named YYYY_MM_DD.nc and holds the six trajectory variables
    and time on (trajectory, obs), with 201 daily values;
  - every particle is released at 12:00 on the day of the file name;
  - the satellite climatology holds its three variables on a regular grid
    and one record for each of the 366 days of the year;
  - the monthly files of the model hold chl on a regular grid, the same
    in every file, with 75 levels and one record for each of the 201 days
    of the cohort.

The check does not stop at the first problem: every cohort is tested and
every problem is listed. It also compares the cohorts found with the
release days expected for the year (15 November to 14 March), as
check_inputs.py does, and reports the share of release positions that
fall on a cell of each grid, as information: a position outside a grid
reads no value and stops nothing.

Only headers, time axes and release positions are read: a year takes well
under a minute.

Environment
-----------
    KRICO_RUNS, KRICO_GLOBCOLOUR, KRICO_PISCES   required
"""

import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from krico_recruitment import chlorophyll, inputs, release_cells  # noqa: E402
import process_cohort_chl  # noqa: E402
from check_inputs import expected_days, folders  # noqa: E402


def check_year(year, runs_dir, satellite, model):
    started = time.time()
    problems, warnings = [], []
    files = []
    for folder in folders(year):
        path = runs_dir / folder
        if not path.is_dir():
            problems.append(f"folder {path} is missing")
            continue
        files += sorted(path.glob("[0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9].nc"))
    print(f"Spawning year {year}: {len(files)} cohorts in "
          f"{', '.join(folders(year))}", flush=True)

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
                        f"the year have no cohort ({absent[0]} to {absent[-1]})")
    for day in sorted(set(found) - expected):
        warnings.append(f"{found[day].parent.name}/{found[day].name} is outside "
                        "the release season of the year")

    accepted = particles = on_satellite = on_model = 0
    for path in files:
        try:
            checked = process_cohort_chl.check(path, satellite, model)
        except (ValueError, FileNotFoundError, OSError) as exc:
            problems.append(str(exc))
            continue
        accepted += 1
        particles += checked["n_particles"]
        for grid, name in ((satellite, "satellite"), (model, "model")):
            inside = grid.lat.index(checked["lat"])[1] & grid.lon.index(checked["lon"])[1]
            if name == "satellite":
                on_satellite += int(inside.sum())
            else:
                on_model += int(inside.sum())

    print(f"  accepted: {accepted} of {len(files)} cohorts, {particles} particles")
    if particles:
        print(f"  release positions on the satellite grid: "
              f"{100 * on_satellite / particles:.3f}%; on the model grid: "
              f"{100 * on_model / particles:.3f}%")
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
    try:
        years = [int(a) for a in args]
    except ValueError:
        sys.exit("arguments: spawning years")
    for year in years:
        if not release_cells.FIRST_YEAR <= year <= release_cells.LAST_YEAR:
            sys.exit(f"year must be in {release_cells.FIRST_YEAR}-"
                     f"{release_cells.LAST_YEAR}")
    if not years:
        years = list(range(release_cells.FIRST_YEAR, release_cells.LAST_YEAR + 1))
    runs = os.environ.get("KRICO_RUNS")
    if not runs:
        sys.exit("ERROR: KRICO_RUNS environment variable not set.")
    satellite_dir, model_dir = process_cohort_chl.directories()
    try:
        satellite, model = process_cohort_chl.open_products(satellite_dir, model_dir)
    except (ValueError, FileNotFoundError, OSError) as exc:
        sys.exit(f"ERROR: {exc}")
    print(f"Satellite: {satellite.path} ({chlorophyll.source('satellite')})")
    print(f"Model: {model_dir} ({chlorophyll.source('model')})")
    ok = True
    for year in years:
        ok &= check_year(year, Path(runs), satellite, model)
    print("Check passed." if ok else "CHECK FAILED.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
