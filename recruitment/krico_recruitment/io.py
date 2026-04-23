"""
NetCDF I/O for the recruitment pipeline.

- read_trajectory_file: load one Parcels trajectory NetCDF into xarray.
- write_recruitment_file: save per-particle recruitment outcomes as
  a CF-compliant NetCDF file, with integer-coded outcome.
"""

import numpy as np
import xarray as xr

from krico_recruitment.outcome import OUTCOME_MEANINGS


# Variable names as written by the production simulation script.
# Matches the ncdump of 1993_11_15.nc supplied during pipeline design.
TRAJ_VARIABLES = (
    "lon",
    "lat",
    "z",
    "time",
    "temperature",
    "bathymetry",
    "sea_ice_area_fraction",
)


def read_trajectory_file(path):
    """
    Load a Parcels trajectory file.

    Parameters
    ----------
    path : str or pathlib.Path

    Returns
    -------
    xr.Dataset
        Dimensions: (trajectory, obs). Variables as listed in TRAJ_VARIABLES.
    """
    ds = xr.open_dataset(path, decode_timedelta=False)
    return ds


def write_recruitment_file(
    path,
    particle_id,
    release_lon,
    release_lat,
    release_depth,
    outcome,
    final_lon,
    final_lat,
    final_depth,
    travel_time,
    trajectory_path_length,
    release_date,
):
    """
    Write per-particle recruitment outcomes to a CF-compliant NetCDF file.

    All per-particle arrays must share the same shape (n_particles,).

    Parameters
    ----------
    path : str or pathlib.Path
        Output file path.
    particle_id : ndarray, int64
    release_lon, release_lat : ndarray, float
        Degrees.
    release_depth : ndarray, float
        Meters, positive down.
    outcome : ndarray, int8
        Integer-coded outcome (see OUTCOME_MEANINGS).
    final_lon, final_lat : ndarray, float
        Degrees.
    final_depth : ndarray, float
        Meters, positive down.
    travel_time : ndarray, int32
        Days from release to moment fate is determined.
    trajectory_path_length : ndarray, float
        Meters.
    release_date : pandas.Timestamp or str
        Release date, stored as a global attribute.
    """
    n_particles = particle_id.size

    ds = xr.Dataset(
        data_vars={
            "release_lon": (
                "particle",
                release_lon.astype(np.float32),
                {
                    "standard_name": "longitude",
                    "long_name": "release longitude",
                    "units": "degrees_east",
                },
            ),
            "release_lat": (
                "particle",
                release_lat.astype(np.float32),
                {
                    "standard_name": "latitude",
                    "long_name": "release latitude",
                    "units": "degrees_north",
                },
            ),
            "release_depth": (
                "particle",
                release_depth.astype(np.float32),
                {
                    "standard_name": "depth",
                    "long_name": "release depth",
                    "units": "m",
                    "positive": "down",
                },
            ),
            "outcome": (
                "particle",
                outcome.astype(np.int8),
                {
                    "long_name": "recruitment outcome",
                    "flag_values": np.arange(len(OUTCOME_MEANINGS), dtype=np.int8),
                    "flag_meanings": " ".join(OUTCOME_MEANINGS),
                },
            ),
            "final_lon": (
                "particle",
                final_lon.astype(np.float32),
                {
                    "standard_name": "longitude",
                    "long_name": "longitude at the moment fate is determined",
                    "units": "degrees_east",
                },
            ),
            "final_lat": (
                "particle",
                final_lat.astype(np.float32),
                {
                    "standard_name": "latitude",
                    "long_name": "latitude at the moment fate is determined",
                    "units": "degrees_north",
                },
            ),
            "final_depth": (
                "particle",
                final_depth.astype(np.float32),
                {
                    "standard_name": "depth",
                    "long_name": "depth at the moment fate is determined",
                    "units": "m",
                    "positive": "down",
                },
            ),
            "travel_time": (
                "particle",
                travel_time.astype(np.int32),
                {
                    "long_name": "days from release to moment fate is determined",
                    "units": "days",
                },
            ),
            "trajectory_path_length": (
                "particle",
                trajectory_path_length.astype(np.float32),
                {
                    "long_name": "cumulative along-track path length up to moment fate is determined",
                    "units": "m",
                },
            ),
        },
        coords={
            "particle": ("particle", particle_id.astype(np.int64)),
        },
        attrs={
            "Conventions": "CF-1.7",
            "title": "KRICO per-particle recruitment outcomes",
            "release_date": str(release_date),
            "source": "krico_recruitment post-processing pipeline",
        },
    )

    # Compression settings matching upstream KRICO convention.
    encoding = {
        var: {"zlib": True, "complevel": 5}
        for var in ds.data_vars
    }

    ds.to_netcdf(path, encoding=encoding)
