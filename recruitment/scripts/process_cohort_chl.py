"""
Pass 2 of dataset v3.0.0: chlorophyll along the trajectories of one
release-day cohort.

Reads one trajectory file and the two chlorophyll products, and writes one
file per cohort with, for every particle, the mean chlorophyll at its
positions over the calyptopis stages and over the furcilia stages up to
its sea-ice advance. Chlorophyll is a diagnostic: no outcome depends on
it. The method is described in krico_recruitment.chlorophyll.

Usage
-----
    python3 process_cohort_chl.py TRAJECTORY.nc OUTPUT.nc
    python3 process_cohort_chl.py TRAJECTORY.nc --check     # inputs only, no output

Input:   $KRICO_RUNS/KRICO_NNNN/YYYY_MM_DD.nc    one cohort, 201 daily values
         $KRICO_GLOBCOLOUR/globcolour_chl_climatology.nc
                                                 satellite climatology,
                                                 366 days of the year
         $KRICO_PISCES/pisces_chl_YYYY_MM.nc     model hindcast, daily, the
                                                 months of the 201 days
Output:  the file named on the command line, e.g.
         recruitment/data_v3.0_chl/YYYY_MM_DD.nc

Neither the pass 1 file of the cohort nor the release-cell file is read.
The two days that bound the windows (calyptopis_end_day, advance_day) are
computed from the trajectory by classify.trajectory_diagnostics, the
function that computes them in pass 1, and written in the output, so that
they can be compared with the pass 1 file (validate_pass2.py).

The output records the commit that produced it: the script stops if
tracked files are modified or if the commit is not pushed (--allow-dirty
for tests).

Environment
-----------
    KRICO_GLOBCOLOUR, KRICO_PISCES   required
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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from krico_recruitment import (chlorophyll, classify, inputs,  # noqa: E402
                               provenance, io as krico_io)


def open_products(satellite_dir, model_dir):
    """The two chlorophyll products, opened and checked."""
    return chlorophyll.Satellite(satellite_dir), chlorophyll.Model(model_dir)


def check(input_path, satellite, model, read_positions=True):
    """
    The checks made before anything is computed, shared with
    check_inputs_chl.py: the trajectory file, and both products over the
    201 days of the cohort.
    """
    checked = inputs.check_trajectory(input_path, read_positions=read_positions)
    chlorophyll.check_inputs(satellite, model, checked["release_day"],
                             checked["n_obs"])
    return checked


def chlorophyll_of_cohort(input_path, satellite, model, log=print):
    """
    The variables of a pass 2 file for one cohort, in memory.

    Returns (particle_id, variables, info).
    """
    log(f"Checking inputs of {input_path}")
    checked = check(input_path, satellite, model)
    day, n_obs = checked["release_day"], checked["n_obs"]
    release_date = pd.Timestamp(datetime.combine(day, datetime.min.time())
                                + timedelta(hours=12))
    log(f"  release {day}, {checked['n_particles']} particles, "
        f"{n_obs} daily values")

    log("Reading trajectories...")
    ds = krico_io.read_trajectory_file(input_path)
    lon = ds["lon"].values
    lat = ds["lat"].values
    depth = ds["z"].values
    temperature = ds["temperature"].values
    bathymetry = ds["bathymetry"].values
    sic = ds["sea_ice_area_fraction"].values
    if "trajectory" in ds.coords:
        particle_id = ds["trajectory"].values.astype(np.int64)
    else:
        particle_id = np.arange(lon.shape[0], dtype=np.int64)
    ds.close()

    log("Windows: development and sea-ice advance, as in pass 1...")
    diagnostics, _ = classify.trajectory_diagnostics(
        lon, lat, depth, temperature, bathymetry, sic, release_date)
    del temperature, bathymetry, sic
    variables = {name: diagnostics[name]
                 for name in ("calyptopis_end_day", "advance_day")}
    del diagnostics
    start, stop = chlorophyll.window_bounds(
        variables["calyptopis_end_day"], variables["advance_day"], n_obs,
        classify.DAY_NONE)

    log("Chlorophyll along the trajectories...")
    variables.update(chlorophyll.along_trajectory(
        lon, lat, depth, start, stop, day, satellite, model, log=log))
    info = {"release_day": day, "release_date": release_date, "n_obs": n_obs}
    return particle_id, variables, info


def summary(variables, log=print):
    """Shares of particles with a value and typical values, for the log."""
    n = variables["window_days"].shape[0]
    for w, window in enumerate(chlorophyll.WINDOWS):
        length = variables["window_days"][:, w]
        some = length > 0
        typical = np.median(length[some]) if some.any() else 0
        log(f"Window {window}: {100 * some.mean():.1f}% of particles have at "
            f"least one day, median {typical:.0f} days")
        rows = [("model, surface", variables["chl_model_surface"][:, w],
                 variables["chl_model_days"][:, w]),
                ("model, top 200 m", variables["chl_model_mean_200m"][:, w],
                 variables["chl_model_days"][:, w]),
                ("model, at depth", variables["chl_model_at_depth"][:, w],
                 variables["chl_model_at_depth_days"][:, w])]
        for j, minimum in enumerate(chlorophyll.MIN_COUNTS):
            rows += [(f"satellite mean, count >= {minimum}",
                      variables["chl_satellite_mean"][:, w, j],
                      variables["chl_satellite_days"][:, w, j]),
                     (f"satellite median, count >= {minimum}",
                      variables["chl_satellite_median"][:, w, j],
                      variables["chl_satellite_days"][:, w, j]),
                     (f"model surface on those days, count >= {minimum}",
                      variables["chl_model_surface_matched"][:, w, j],
                      variables["chl_matched_days"][:, w, j])]
        for label, value, days in rows:
            has = days > 0
            if not has.any():
                log(f"  {label:42s} no particle with a value")
                continue
            covered = days[has] / length[has]
            q = np.percentile(value[has].astype(np.float64), (5, 50, 95))
            log(f"  {label:42s} {100 * has.sum() / n:5.1f}% of particles, "
                f"{100 * np.median(covered):5.1f}% of their days (median); "
                f"mg m-3: 5% {q[0]:.3f}, median {q[1]:.3f}, 95% {q[2]:.3f}")


def process_cohort(input_path, output_path, satellite_dir, model_dir,
                   allow_dirty=False):
    """Compute the chlorophyll indices of one cohort and write its file."""
    started = time.time()
    input_path, output_path = Path(input_path), Path(output_path)
    if not output_path.parent.is_dir():
        raise FileNotFoundError(f"output folder {output_path.parent} does not "
                                "exist; create it first")
    identity = provenance.code_identity(script=__file__, allow_dirty=allow_dirty)
    print(f"Code: {identity['code_version']} ({identity['code_commit'][:7]})")

    satellite, model = open_products(satellite_dir, model_dir)
    try:
        particle_id, variables, info = chlorophyll_of_cohort(
            input_path, satellite, model)
    finally:
        satellite.close()
        model.close()

    attributes = {
        "release_date": str(info["release_date"]),
        "source": "krico_recruitment post-processing pipeline",
        "dataset_version": krico_io.DATASET_VERSION,
        **identity,
        "history": provenance.history_line(
            f"process_cohort_chl.py {input_path.name} {output_path.name}"),
        "trajectory_source": f"{input_path.parent.name}/{input_path.name}",
        "satellite_source": chlorophyll.source("satellite"),
        "satellite_comment": (
            "daily climatology of surface chlorophyll, read by day of the "
            "year at each position; it describes the years of the product, "
            "not the year of the trajectory, and holds no value under sea "
            "ice or in the dark"),
        "model_source": chlorophyll.source("model"),
        "model_comment": (
            "daily means read on the actual date of each position; the "
            "hindcast has its own ocean and sea ice, which are not those "
            "of the trajectories"),
        "sampling": (
            "value of the cell holding the position (nearest grid point), "
            "no interpolation, no filling from neighbours; day k of a "
            "trajectory is the calendar day of release plus k"),
        "windows_comment": (
            "the windows are bounded by calyptopis_end_day and advance_day, "
            "computed from the trajectory as in pass 1 with the sea ice of "
            "the trajectory file (GLORYS12); every particle has its "
            "indices, whatever its outcome"),
    }

    print(f"Writing {output_path}")
    krico_io.write_chlorophyll_file(output_path, particle_id, variables,
                                    attributes)
    print()
    summary(variables)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6   # kB on Linux
    print(f"Done: {output_path.name}, {output_path.stat().st_size / 1e6:.1f} MB, "
          f"{time.time() - started:.0f} s, peak memory {peak:.1f} GB")


def directories():
    """The folders of the two products, from the environment."""
    found = []
    for name in ("KRICO_GLOBCOLOUR", "KRICO_PISCES"):
        value = os.environ.get(name)
        if not value:
            sys.exit(f"ERROR: {name} environment variable not set.")
        found.append(Path(value))
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("input", type=Path,
                        help="trajectory file, YYYY_MM_DD.nc")
    parser.add_argument("output", type=Path, nargs="?",
                        help="chlorophyll file to write")
    parser.add_argument("--check", action="store_true",
                        help="check the inputs of this cohort and write nothing")
    parser.add_argument("--allow-dirty", action="store_true",
                        help="for tests: accept modified or unpushed code")
    args = parser.parse_args()
    satellite_dir, model_dir = directories()
    try:
        if args.check:
            satellite, model = open_products(satellite_dir, model_dir)
            checked = check(args.input, satellite, model)
            print(f"{args.input.name}: inputs accepted, "
                  f"{checked['n_particles']} particles")
            return
        if args.output is None:
            parser.error("give the output file, or --check")
        process_cohort(args.input, args.output, satellite_dir, model_dir,
                       allow_dirty=args.allow_dirty)
    except (RuntimeError, ValueError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: {exc}")


if __name__ == "__main__":
    main()
