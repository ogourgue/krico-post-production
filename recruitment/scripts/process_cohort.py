"""
Process one release-day cohort: read trajectory NetCDF, evaluate all
filters, determine per-particle outcomes, write per-particle recruitment
NetCDF.

Usage
-----
    python process_cohort.py <input_trajectory.nc> <output_recruitment.nc>

The input file is one of the KRICO simulation outputs (e.g. 1993_11_15.nc).
The output file contains one row per particle with the 9-column schema
described in the Paper 1 methodology document.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Allow running the script directly from the `scripts/` folder: add the
# parent directory (which contains the krico_recruitment package) to sys.path.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from krico_recruitment import development, filters, sea_ice, trajectory, io as krico_io
from krico_recruitment.outcome import OUTCOME_CODES


# Bathymetry threshold for the shelf / shelf-slope recruitment zone (meters).
SHELF_DEPTH_THRESHOLD = 2000.0


def process_cohort(input_path, output_path):
    """
    Run the full pipeline on one cohort and write the output file.
    """
    print(f"Reading trajectory file: {input_path}")
    ds = krico_io.read_trajectory_file(input_path)

    n_particles = ds.sizes["trajectory"]
    n_obs = ds.sizes["obs"]
    print(f"  particles: {n_particles}, obs per particle: {n_obs}")

    # Extract arrays. Shape (n_particles, n_obs).
    lon = ds["lon"].values
    lat = ds["lat"].values
    depth = ds["z"].values
    temperature = ds["temperature"].values
    sic = ds["sea_ice_area_fraction"].values
    bathymetry = ds["bathymetry"].values
    time_ns = ds["time"].values  # nanoseconds since 1970-01-01

    # Release date from day 0 time of particle 0.
    # (All particles in a cohort share the release date by design.)
    release_date = pd.Timestamp(time_ns[0, 0])
    print(f"  release date: {release_date.date()}")

    # ------------------------------------------------------------------
    # Step 2: T-dependent development
    # ------------------------------------------------------------------
    print("Integrating T-dependent development...")
    cumulative_dev = development.integrate_development(temperature)
    calyptope_end_day = development.calyptope_window_end(cumulative_dev)
    fiv_day = development.days_to_reach_stage(
        cumulative_dev, development.FIV_INDEX
    )

    # ------------------------------------------------------------------
    # Step 3: M1 (spawning SIC >= 80%)
    # ------------------------------------------------------------------
    print("Evaluating M1...")
    sic_at_release = sic[:, 0]
    killed_M1 = filters.evaluate_M1(sic_at_release)

    # ------------------------------------------------------------------
    # Step 4: M4 (calyptope starvation)
    # ------------------------------------------------------------------
    print("Evaluating M4...")
    killed_M4, M4_kill_day = filters.evaluate_M4(sic, calyptope_end_day)
    # A particle killed by M1 cannot subsequently be killed by M4.
    killed_M4 = killed_M4 & ~killed_M1

    # ------------------------------------------------------------------
    # Step 5: sea-ice advance detection
    # ------------------------------------------------------------------
    print("Detecting sea-ice advance...")
    april1_day = sea_ice.find_april1_day_index(release_date, n_obs)
    print(f"  April 1 day index: {april1_day}")
    advance_day = sea_ice.detect_sea_ice_advance(sic, april1_day)

    # Censoring: no advance but SIC rising at end of tracking.
    censored = sea_ice.check_censoring(sic, n_obs, advance_day)
    # Censoring does not apply to particles already killed by M1 or M4.
    censored = censored & ~killed_M1 & ~killed_M4

    # ------------------------------------------------------------------
    # Step 6: evaluate M5 / M6 at the advance event (for still-live particles)
    # ------------------------------------------------------------------
    print("Evaluating M5 / M6 at advance event...")

    # Candidates: alive after M1 and M4.
    alive = ~killed_M1 & ~killed_M4

    # Particles with no advance detected AND not censored → killed_M6_no_advance.
    no_advance = alive & (advance_day == -1) & ~censored

    # For particles with an advance event, check M5 conditions at that day.
    has_advance = alive & (advance_day != -1)

    # Reached FIV and on shelf at advance event.
    # "On shelf" means water depth (bathymetry) < 2000 m, i.e. the particle
    # is over the continental shelf or shelf-slope. The particle's own
    # z-coordinate is always shallow (upper 200 m).
    dev_at_advance = trajectory.extract_at_day(cumulative_dev, advance_day)
    bathy_at_advance = trajectory.extract_at_day(bathymetry, advance_day)

    reached_FIV_at_advance = dev_at_advance >= development.FIV_INDEX
    on_shelf_at_advance = bathy_at_advance < SHELF_DEPTH_THRESHOLD

    killed_M5_no_FIV = has_advance & ~reached_FIV_at_advance
    killed_M5_not_on_shelf = (
        has_advance & reached_FIV_at_advance & ~on_shelf_at_advance
    )
    success = has_advance & reached_FIV_at_advance & on_shelf_at_advance

    # ------------------------------------------------------------------
    # Step 7: assemble outcome, final position, travel time, path length
    # ------------------------------------------------------------------
    print("Assembling outcomes...")

    # "Moment fate is determined" day index per particle.
    # Default to -1, then fill in per outcome.
    fate_day = np.full(n_particles, -1, dtype=np.int64)

    # killed_M1: day 0 (release).
    fate_day[killed_M1] = 0

    # killed_M4: day threshold was crossed.
    fate_day[killed_M4] = M4_kill_day[killed_M4]

    # success, killed_M5_no_FIV, killed_M5_not_on_shelf: advance day.
    advance_mask = success | killed_M5_no_FIV | killed_M5_not_on_shelf
    fate_day[advance_mask] = advance_day[advance_mask]

    # censored and killed_M6_no_advance: last day of tracking.
    end_mask = censored | no_advance
    fate_day[end_mask] = n_obs - 1

    # Outcome code.
    outcome_code = np.full(n_particles, -1, dtype=np.int8)
    outcome_code[success] = OUTCOME_CODES["success"]
    outcome_code[censored] = OUTCOME_CODES["censored"]
    outcome_code[killed_M1] = OUTCOME_CODES["killed_M1"]
    outcome_code[killed_M4] = OUTCOME_CODES["killed_M4"]
    outcome_code[killed_M5_no_FIV] = OUTCOME_CODES["killed_M5_no_FIV"]
    outcome_code[killed_M5_not_on_shelf] = OUTCOME_CODES["killed_M5_not_on_shelf"]
    outcome_code[no_advance] = OUTCOME_CODES["killed_M6_no_advance"]

    # Sanity check: every particle should have a single outcome.
    assigned = (outcome_code >= 0)
    if not assigned.all():
        n_unassigned = (~assigned).sum()
        raise RuntimeError(
            f"{n_unassigned} particles not assigned to any outcome. Bug in logic."
        )

    # Final positions at fate_day.
    # For particles with fate_day valid, pick position at that day.
    release_lon_arr = lon[:, 0]
    release_lat_arr = lat[:, 0]
    release_depth_arr = depth[:, 0]

    final_lon_arr = trajectory.extract_at_day(lon, fate_day)
    final_lat_arr = trajectory.extract_at_day(lat, fate_day)
    final_depth_arr = trajectory.extract_at_day(depth, fate_day)

    # Travel time = fate_day in days (daily output → day index == days).
    travel_time_arr = fate_day.astype(np.int32)
    # For killed_M1, fate_day == 0, travel_time == 0. Consistent with spec.

    # Trajectory path length up to fate_day.
    print("Computing trajectory path length...")
    cum_path = trajectory.cumulative_path_length_m(lon, lat)
    path_length_arr = trajectory.extract_at_day(cum_path, fate_day)

    # Particle IDs from the trajectory coord (if present) or 0..n-1.
    if "trajectory" in ds.coords:
        particle_id = ds["trajectory"].values.astype(np.int64)
    else:
        particle_id = np.arange(n_particles, dtype=np.int64)

    # ------------------------------------------------------------------
    # Step 8: write output file
    # ------------------------------------------------------------------
    print(f"Writing output: {output_path}")
    krico_io.write_recruitment_file(
        path=output_path,
        particle_id=particle_id,
        release_lon=release_lon_arr,
        release_lat=release_lat_arr,
        release_depth=release_depth_arr,
        outcome=outcome_code,
        final_lon=final_lon_arr,
        final_lat=final_lat_arr,
        final_depth=final_depth_arr,
        travel_time=travel_time_arr,
        trajectory_path_length=path_length_arr,
        release_date=release_date,
    )

    # Quick summary to stdout for validation.
    print()
    print("Outcome summary:")
    unique, counts = np.unique(outcome_code, return_counts=True)
    for code, count in zip(unique, counts):
        label = list(OUTCOME_CODES.keys())[list(OUTCOME_CODES.values()).index(code)]
        pct = 100.0 * count / n_particles
        print(f"  {label:28s} {count:8d}  ({pct:5.2f}%)")


def main():
    parser = argparse.ArgumentParser(
        description="Process one KRICO cohort into a per-particle recruitment file."
    )
    parser.add_argument("input", type=str, help="Input trajectory NetCDF (e.g. 1993_11_15.nc)")
    parser.add_argument("output", type=str, help="Output recruitment NetCDF")
    args = parser.parse_args()

    process_cohort(args.input, args.output)


if __name__ == "__main__":
    main()
