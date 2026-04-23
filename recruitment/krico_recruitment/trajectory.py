"""
Trajectory-level utilities: along-track path length, final position
extraction at arbitrary per-particle day indices.
"""

import numpy as np


# Earth radius for great-circle distance (meters).
EARTH_RADIUS_M = 6_371_000.0


def haversine_distance_m(lon1, lat1, lon2, lat2):
    """
    Great-circle distance between two sets of points, in meters.

    All inputs in degrees; broadcasting allowed.
    """
    lon1_r = np.radians(lon1)
    lat1_r = np.radians(lat1)
    lon2_r = np.radians(lon2)
    lat2_r = np.radians(lat2)

    dlon = lon2_r - lon1_r
    dlat = lat2_r - lat1_r

    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1_r) * np.cos(lat2_r) * np.sin(dlon / 2.0) ** 2
    c = 2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))
    return EARTH_RADIUS_M * c


def cumulative_path_length_m(lon, lat):
    """
    Cumulative along-track path length (meters) for each particle,
    computed segment by segment with haversine distance.

    Parameters
    ----------
    lon, lat : ndarray of shape (n_particles, n_obs)
        Longitude and latitude in degrees. NaN segments are skipped
        (cumulative distance does not increase when either endpoint is NaN).

    Returns
    -------
    cum_length : ndarray of shape (n_particles, n_obs), dtype float32
        Cumulative distance from day 0 to day t, in meters. NaN where
        the position itself is NaN.
    """
    n_particles, n_obs = lon.shape

    # Segment distances from day t to day t+1.
    seg = haversine_distance_m(
        lon[:, :-1], lat[:, :-1], lon[:, 1:], lat[:, 1:]
    )
    # Zero out NaN segments (e.g. post-deletion).
    seg = np.nan_to_num(seg, nan=0.0)

    cum = np.zeros((n_particles, n_obs), dtype=np.float32)
    cum[:, 1:] = np.cumsum(seg, axis=1)

    # Propagate NaN where the position itself is NaN.
    cum = np.where(np.isnan(lon) | np.isnan(lat), np.nan, cum)
    return cum


def extract_at_day(array, day_indices):
    """
    Pick array[i, day_indices[i]] for each particle i, with safe handling
    of -1 / out-of-bounds values (returns NaN for those).

    Parameters
    ----------
    array : ndarray of shape (n_particles, n_obs)
    day_indices : ndarray of shape (n_particles,), dtype int

    Returns
    -------
    values : ndarray of shape (n_particles,)
        array[i, day_indices[i]] where day_indices[i] is a valid index,
        NaN otherwise.
    """
    n_particles, n_obs = array.shape
    valid = (day_indices >= 0) & (day_indices < n_obs)
    safe_idx = np.where(valid, day_indices, 0)
    values = array[np.arange(n_particles), safe_idx]

    # Ensure float output for NaN handling; caller must pass floatable data.
    if values.dtype.kind not in ("f", "c"):
        values = values.astype(np.float32)
    values = np.where(valid, values, np.nan)
    return values
