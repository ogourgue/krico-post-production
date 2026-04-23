"""
Synthetic smoke test for the pipeline: generates a tiny fake trajectory
file with a handful of particles designed to exercise each outcome state,
runs process_cohort.py, and inspects the output.
"""

import sys
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))

from scripts.process_cohort import process_cohort
from krico_recruitment.outcome import OUTCOME_MEANINGS


def make_fake_cohort():
    """
    Six particles, each designed to exercise one outcome state.

    Particle 0: success         — ice advances on day 150, at FIV, on shelf
    Particle 1: censored        — no advance, but SIC rising at end of tracking
    Particle 2: killed_M1       — SIC at release >= 80%
    Particle 3: killed_M4       — 15 consecutive days under SIC >40% early on
    Particle 4: killed_M5_no_FIV — ice advances day 130, in cold water (FIV not yet reached)
    Particle 5: killed_M5_not_on_shelf — ice advances at FIV but bathymetry >2000
    Particle 6: killed_M6       — no advance, SIC low/flat
    """
    n_particles = 7
    n_obs = 201

    lon = np.zeros((n_particles, n_obs), dtype=np.float32)
    lat = np.full((n_particles, n_obs), -65.0, dtype=np.float32)  # warm-ish
    z = np.full((n_particles, n_obs), 100.0, dtype=np.float32)
    temperature = np.full((n_particles, n_obs), 0.5, dtype=np.float32)  # moderate
    sic = np.zeros((n_particles, n_obs), dtype=np.float32)
    bathymetry = np.full((n_particles, n_obs), 1500.0, dtype=np.float32)  # on-shelf

    # Release date: 15 Nov 2015 (so April 1 → day 138).
    release = pd.Timestamp("2015-11-15")
    times = np.array(
        [release + pd.Timedelta(days=t) for t in range(n_obs)],
        dtype="datetime64[ns]",
    )
    time_arr = np.broadcast_to(times, (n_particles, n_obs)).astype("int64")

    # Particle 0: success. Ice advance triggered from day 150.
    sic[0, 150:] = 0.6

    # Particle 1: censored. SIC rising in the last 30 days but not crossing
    # the 5-consecutive-day >15% threshold.
    sic[1, -30:] = np.linspace(0.05, 0.14, 30)

    # Particle 2: killed_M1. SIC at release >= 80%.
    sic[2, 0] = 0.85

    # Particle 3: killed_M4. 15 days >40% SIC early on (before CIII finishes,
    # which at T=0.5°C happens around day 35-40 based on Thorpe Table 1).
    sic[3, 5:25] = 0.5

    # Particle 4: killed_M5_no_FIV. Very cold water so development is slow,
    # then ice advance happens before FIV is reached.
    temperature[4, :] = -1.0  # cold; FIV reached around day ~110-120
    sic[4, 70:] = 0.6  # ice advance triggered on day 70 (before FIV)

    # Particle 5: killed_M5_not_on_shelf. Ice advance at FIV but in deep water.
    bathymetry[5, :] = 3500.0  # off-shelf
    sic[5, 150:] = 0.6

    # Particle 6: killed_M6. No advance event, no rising trend.
    # SIC stays at 0 throughout — default.

    return xr.Dataset(
        data_vars={
            "lon": (("trajectory", "obs"), lon),
            "lat": (("trajectory", "obs"), lat),
            "z": (("trajectory", "obs"), z),
            "temperature": (("trajectory", "obs"), temperature),
            "sea_ice_area_fraction": (("trajectory", "obs"), sic),
            "bathymetry": (("trajectory", "obs"), bathymetry),
            "time": (("trajectory", "obs"), time_arr),
        },
        coords={
            "trajectory": np.arange(n_particles, dtype=np.int64),
            "obs": np.arange(n_obs, dtype=np.int64),
        },
    )


def main():
    ds_fake = make_fake_cohort()

    with tempfile.TemporaryDirectory() as tmp:
        in_path = Path(tmp) / "fake_trajectory.nc"
        out_path = Path(tmp) / "fake_recruitment.nc"

        ds_fake.to_netcdf(in_path)

        process_cohort(str(in_path), str(out_path))

        print()
        print("Reading back and checking each particle:")
        ds_out = xr.open_dataset(out_path)
        for i in range(ds_out.sizes["particle"]):
            code = int(ds_out["outcome"].values[i])
            label = OUTCOME_MEANINGS[code]
            tt = int(ds_out["travel_time"].values[i])
            pl = float(ds_out["trajectory_path_length"].values[i])
            print(f"  particle {i}: {label:28s}  travel_time={tt:3d}  path_length={pl:.1f} m")


if __name__ == "__main__":
    main()
