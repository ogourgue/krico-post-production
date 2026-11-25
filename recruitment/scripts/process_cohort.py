"""
Pass 1 of dataset v3.0.0: classify one release-day cohort.

Reads one trajectory file and the release-cell file of its spawning year,
computes the diagnostics of every particle, derives the outcome from them
at the baseline parameters and writes one file per cohort.

Usage
-----
    python3 process_cohort.py TRAJECTORY.nc OUTPUT.nc
    python3 process_cohort.py TRAJECTORY.nc --check     # inputs only, no output
    python3 process_cohort.py TRAJECTORY.nc OUTPUT.nc --ice-source nsidc
    python3 process_cohort.py TRAJECTORY.nc OUTPUT.nc --ice-source glorys12_25km

Input:   $KRICO_RUNS/KRICO_NNNN/YYYY_MM_DD.nc    one cohort, 201 daily values
         $KRICO_GLORYS12_RELEASE_CELLS/glorys12_release_cells_YYYY.nc
                                                 sea ice and temperature at the
                                                 release cells (pass 0)
Output:  the file named on the command line, e.g.
         recruitment/data_v3.0_glorys12/YYYY_MM_DD.nc

No GLORYS12 forcing file is read: what the rules need before release (M1
at spawning, M3 during the descent and the ascent, the temperature at
spawning) comes from the release-cell file, and the rest from the
trajectory.

Steps
-----
1. Inputs are checked (krico_recruitment.inputs.check_cohort, the function
   check_inputs.py runs on every cohort beforehand) and every particle is
   linked to its release cell, the grid point nearest to its release
   position.
2. Before release: sea ice at spawning, M3 temperatures, temperature at
   spawning (krico_recruitment.pre_release).
3. Along the trajectory: development, sea-ice runs, advance, positions
   (krico_recruitment.classify.trajectory_diagnostics).
4. The outcome is computed from the diagnostics in their stored form
   (classify.outcome), so that any reader of the file recomputes it
   exactly. The position, day and path length at which fate is recorded
   follow from it, as in v2.0.0; a particle removed by M3 is recorded like
   one removed by M1, at its release position on day 0.
5. Day-0 check: the trajectory's own sea ice and temperature at day 0
   against the release-cell file. The result is written in the global
   attributes; it does not stop the run.

Ice source
----------
By default every sea-ice concentration is that of GLORYS12: at spawning
from the release-cell file, along the trajectory from the trajectory file.
With --ice-source nsidc, every concentration the classification reads (M1
at spawning, the runs behind M4, the advance, censoring, the daily bits)
comes from the NSIDC G02202 v6 record instead, sampled at the same
positions and days; the trajectories, and the temperatures behind M3 and
development, are unchanged. This is the ice-source test for the
ERA-Interim/ERA5 forcing change in GLORYS12 (krico_recruitment.nsidc_ice).
--ice-source glorys12_25km reads GLORYS12 ice binned onto the NSIDC grid
through the same sampler, so that the difference between the two is the
ice values alone. The monthly files are found in --ice-dir, or in
$KRICO_NSIDC and $KRICO_GLORYS12_25KM. The day-0 check always compares the
GLORYS12 values of the trajectory file with the release-cell file.

The output records the commit that produced it: the script stops if
tracked files are modified or if the commit is not pushed (--allow-dirty
for tests).

Environment
-----------
    KRICO_GLORYS12_RELEASE_CELLS   required
    KRICO_NSIDC, KRICO_GLORYS12_25KM   monthly files of the other ice sources
"""

import argparse
import os
import resource
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

# Allow running the script directly from the `scripts/` folder: add the
# parent directory (which contains the krico_recruitment package) to sys.path.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from krico_recruitment import (classify, inputs, nsidc_ice, pre_release,  # noqa: E402
                               provenance, sea_ice, trajectory, io as krico_io)
from krico_recruitment.outcome import OUTCOME_MEANINGS  # noqa: E402

# Share of particles whose day-0 values must lie within
# pre_release.DAY0_TOLERANCE of the release-cell file for the check to be
# reported as expected. The share exactly equal is reported too, but drops
# for sea ice from 19 December 2018, whose records are stamped 11:52:30.
DAY0_EXPECTED_SHARE = 0.999


def classify_cohort(input_path, release_cells_dir,
                    parameters=classify.BASELINE, log=print,
                    ice_source="glorys12", ice_dir=None):
    """
    Diagnostics and outcome of one cohort, in memory.

    Returns (particle_id, variables, info): the variables of a v3.0.0 file
    as a dict of arrays, and info with release_day, static, n_obs and the
    day-0 check. parameters other than the baseline are for tests (the
    regression against v2.0.0 switches M3 off); a file is only ever
    written at the baseline. ice_source and ice_dir: see the module
    docstring.
    """
    log(f"Checking inputs of {input_path}")
    checked = inputs.check_cohort(input_path, release_cells_dir,
                                  ice_source=ice_source, ice_dir=ice_dir)
    gridded = nsidc_ice.SOURCES.get(ice_source)
    day, static, cell = checked["release_day"], checked["static"], checked["cell"]
    release_date = pd.Timestamp(datetime.combine(day, datetime.min.time())
                                + timedelta(hours=12))
    log(f"  release {day}, {checked['n_particles']} particles, "
        f"{checked['n_obs']} daily values; release cells from "
        f"{static['path'].name} ({static['code_version']})")

    log("Before release: sea ice at spawning, M3 temperatures...")
    before = pre_release.evaluate(static, day, cell)

    log("Reading trajectories...")
    ds = krico_io.read_trajectory_file(input_path)
    lon = ds["lon"].values
    lat = ds["lat"].values
    depth = ds["z"].values
    temperature = ds["temperature"].values
    bathymetry = ds["bathymetry"].values
    if gridded:
        # GLORYS12 sea ice of the trajectory at day 0, for the day-0 check.
        sic_day0 = ds["sea_ice_area_fraction"].isel(obs=0).values
    else:
        sic = ds["sea_ice_area_fraction"].values
        sic_day0 = sic[:, 0]
    if "trajectory" in ds.coords:
        particle_id = ds["trajectory"].values.astype(np.int64)
    else:
        particle_id = np.arange(lon.shape[0], dtype=np.int64)
    ds.close()

    day0 = pre_release.day0_check(static, day, cell, depth[:, 0], sic_day0,
                                  temperature[:, 0])

    if gridded:
        # The other ice source replaces GLORYS12 everywhere the rules read
        # sea ice: at spawning, and along the trajectory.
        log(f"Sampling {gridded['id']} from {ice_dir}...")
        before["sic_spawning"] = nsidc_ice.spawning_sic(
            lon[:, 0], lat[:, 0], release_date, ice_dir,
            sea_ice.SPAWNING_OFFSET_DAYS, pattern=gridded["pattern"])
        sic = nsidc_ice.along_trajectory_sic(lon, lat, release_date, ice_dir,
                                             pattern=gridded["pattern"])
    d = classify.pre_release_diagnostics(before, cell)

    log("Along the trajectory: development, sea-ice runs, advance...")
    diagnostics, last_valid = classify.trajectory_diagnostics(
        lon, lat, depth, temperature, bathymetry, sic, release_date)
    d.update(diagnostics)
    del temperature, bathymetry, sic

    log("Outcome from the diagnostics...")
    code = classify.outcome(d, parameters)
    fate = classify.fate_day(d, code, last_valid)

    variables = {
        "release_lon": lon[:, 0].astype(np.float32),
        "release_lat": lat[:, 0].astype(np.float32),
        "release_depth": depth[:, 0].astype(np.float32),
        "outcome": code,
        "final_lon": trajectory.extract_at_day(lon, fate).astype(np.float32),
        "final_lat": trajectory.extract_at_day(lat, fate).astype(np.float32),
        "final_depth": trajectory.extract_at_day(depth, fate).astype(np.float32),
        "travel_time": fate.astype(np.int32),
    }
    del depth
    cum_path = trajectory.cumulative_path_length_m(lon, lat)
    variables["trajectory_path_length"] = trajectory.extract_at_day(
        cum_path, fate).astype(np.float32)
    del cum_path, lon, lat
    variables.update(d)

    info = {"release_day": day, "release_date": release_date, "static": static,
            "n_obs": checked["n_obs"], "day0": day0,
            "sic_source": gridded["id"] if gridded else "GLORYS12"}
    return particle_id, variables, info


def summary(variables, info, log=print):
    """Outcome shares and day-0 check, for the log."""
    code = variables["outcome"]
    n = code.size
    log("Outcome summary:")
    for value, label in enumerate(OUTCOME_MEANINGS):
        count = int((code == value).sum())
        log(f"  {value} {label:26s} {count:8d}  ({100.0 * count / n:5.2f}%)")
    log("Day-0 check, trajectory against the release-cell file:")
    for name, s in info["day0"].items():
        note = ""
        if not s["within_tolerance_share"] >= DAY0_EXPECTED_SHARE:
            note = f"   WARNING: below {100 * DAY0_EXPECTED_SHARE}% within tolerance"
        log(f"  {name:12s} {100 * s['exact_share']:8.4f}% exactly equal, "
            f"{100 * s['within_tolerance_share']:8.4f}% within "
            f"{pre_release.DAY0_TOLERANCE}, mean difference "
            f"{s['mean_abs_difference']:.3e}, largest "
            f"{s['max_abs_difference']:.3e}, {s['compared']} particles "
            f"compared{note}")


def process_cohort(input_path, output_path, release_cells_dir,
                   allow_dirty=False, ice_source="glorys12", ice_dir=None):
    """Classify one cohort and write its v3.0.0 file."""
    started = time.time()
    input_path, output_path = Path(input_path), Path(output_path)
    if not output_path.parent.is_dir():
        raise FileNotFoundError(f"output folder {output_path.parent} does not "
                                "exist; create it first")
    identity = provenance.code_identity(script=__file__, allow_dirty=allow_dirty)
    print(f"Code: {identity['code_version']} ({identity['code_commit'][:7]})")

    particle_id, variables, info = classify_cohort(
        input_path, release_cells_dir, ice_source=ice_source, ice_dir=ice_dir)
    static, day = info["static"], info["release_day"]
    command = f"process_cohort.py {input_path.name} {output_path.name}"
    if ice_source != "glorys12":
        command += f" --ice-source {ice_source}"

    attributes = {
        "release_date": str(info["release_date"]),
        "source": "krico_recruitment post-processing pipeline",
        "dataset_version": krico_io.DATASET_VERSION,
        **identity,
        "history": provenance.history_line(command),
        "trajectory_source": f"{input_path.parent.name}/{input_path.name}",
        "release_cells_source": (f"{static['path'].name} "
                                 f"({static['code_version']}, "
                                 f"{static['code_commit'][:7]})"),
        "m1_spawning_offset_days": np.int64(sea_ice.SPAWNING_OFFSET_DAYS),
        "m1_spawning_date": str(day - timedelta(days=sea_ice.SPAWNING_OFFSET_DAYS)),
        "sic_source": info["sic_source"],
    }
    if ice_source != "glorys12":
        attributes["sic_sampling"] = nsidc_ice.SAMPLING
        attributes["sic_comment"] = (
            "every sea-ice concentration read by the classification (at "
            "spawning and along the trajectory) comes from sic_source; the "
            "trajectories and the temperatures are those of GLORYS12, and "
            "the day-0 check compares the GLORYS12 sea ice of the "
            "trajectory file with the release-cell file")
    for name, s in info["day0"].items():
        attributes[f"day0_{name}_exact_share"] = np.float64(s["exact_share"])
        attributes[f"day0_{name}_within_tolerance_share"] = np.float64(
            s["within_tolerance_share"])
        attributes[f"day0_{name}_mean_abs_difference"] = np.float64(
            s["mean_abs_difference"])
        attributes[f"day0_{name}_max_abs_difference"] = np.float64(s["max_abs_difference"])
        attributes[f"day0_{name}_particles_compared"] = np.int64(s["compared"])
    attributes["day0_tolerance"] = np.float64(pre_release.DAY0_TOLERANCE)
    attributes["day0_comment"] = (
        "sea ice and temperature of the trajectory at day 0 against the "
        "release-cell file on the release day, at the release cell and at "
        "the level holding the release depth: share of particles exactly "
        "equal, share within day0_tolerance, mean and largest absolute "
        "difference. From 19 December 2018 the sea-ice records of the "
        "forcing are stamped 11:52:30, so a trajectory's sea ice at day 0 "
        "holds 0.5% of the next day and is no longer exactly equal")

    print(f"Writing {output_path}")
    krico_io.write_classification_file(output_path, particle_id, variables,
                                       info["n_obs"], attributes)
    print()
    summary(variables, info)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6   # kB on Linux
    print(f"Done: {output_path.name}, {output_path.stat().st_size / 1e6:.1f} MB, "
          f"{time.time() - started:.0f} s, peak memory {peak:.1f} GB")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("input", type=Path,
                        help="trajectory file, YYYY_MM_DD.nc")
    parser.add_argument("output", type=Path, nargs="?",
                        help="classification file to write")
    parser.add_argument("--check", action="store_true",
                        help="check the inputs of this cohort and write nothing")
    parser.add_argument("--allow-dirty", action="store_true",
                        help="for tests: accept modified or unpushed code")
    parser.add_argument("--ice-source", default="glorys12",
                        choices=("glorys12", *nsidc_ice.SOURCES),
                        help="where sea-ice concentration is read from "
                             "(default glorys12)")
    parser.add_argument("--ice-dir", type=Path, default=None,
                        help="folder of the monthly files of the ice source "
                             "(default: $KRICO_NSIDC or $KRICO_GLORYS12_25KM)")
    args = parser.parse_args()
    directory = os.environ.get("KRICO_GLORYS12_RELEASE_CELLS")
    if not directory:
        sys.exit("ERROR: KRICO_GLORYS12_RELEASE_CELLS environment variable not set.")
    ice_dir = args.ice_dir
    if args.ice_source != "glorys12" and ice_dir is None:
        variable = nsidc_ice.SOURCES[args.ice_source]["env"]
        ice_dir = os.environ.get(variable)
        if not ice_dir:
            sys.exit(f"ERROR: {variable} environment variable not set "
                     "(or give --ice-dir).")
    try:
        if args.check:
            checked = inputs.check_cohort(args.input, directory,
                                          ice_source=args.ice_source,
                                          ice_dir=ice_dir)
            print(f"{args.input.name}: inputs accepted, "
                  f"{checked['n_particles']} particles")
            return
        if args.output is None:
            parser.error("give the output file, or --check")
        process_cohort(args.input, args.output, directory,
                       allow_dirty=args.allow_dirty,
                       ice_source=args.ice_source, ice_dir=ice_dir)
    except (RuntimeError, ValueError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: {exc}")


if __name__ == "__main__":
    main()
