"""
Sea-ice advance detection (per-particle winter onset), and sea-ice
concentration at spawning.

Advance
-------
For each particle, winter onset is defined as the first date on or after
April 1 of the relevant austral autumn where sea ice concentration at the
particle's location exceeds 15% for 5 consecutive days.

The April 1 cutoff avoids spurious detection from lingering previous-winter
ice. The austral autumn is defined relative to the release date:
  - release in Jan-Mar → April 1 of same calendar year
  - release in Nov-Dec → April 1 of following calendar year
Equivalently: the first April 1 that falls within the tracking window.

Spawning
--------
M1 represents a constraint acting on the spawning adult, but the 23-26 day
descent-ascent cycle from spawning to calyptopis I is not simulated, so
particles enter the model as calyptopis I and their release date postdates
spawning. Reading sea-ice concentration from the trajectory at day 0
therefore samples the field at the wrong time. spawning_sic samples the
GLORYS12 field at the release position on the spawning date instead.

The release position stands as the spawning position, which is consistent
with the model's neglect of transport during the descent-ascent interval.

Sampling is nearest-neighbour on the reanalysis grid, whereas Parcels
interpolates when sampling along trajectories. The two agree to a mean
absolute difference of 2.4e-06, with disagreement above 0.01 confined to
0.005% of particles in coastal cells where the field has a sharp gradient.
See krico-paper1/S1_m1_offset_sensitivity/.
"""

import numpy as np
import pandas as pd
import xarray as xr


# Sea-ice advance parameters (Thorpe 2019, Stammerjohn et al. 2008).
ADVANCE_SIC_THRESHOLD = 0.15       # fraction, 0-1; strictly greater than this
ADVANCE_MIN_CONSECUTIVE_DAYS = 5   # SIC above threshold for this many days

# Tolerance for "alive at end of tracking". Parcels can produce a NaN
# position on the very last output step due to time-interpolation precision
# at the simulation end, even for particles that survived. We therefore
# treat last_valid in [n_obs - 1 - END_OF_TRACKING_TOLERANCE, n_obs - 1]
# as "alive at end".
END_OF_TRACKING_TOLERANCE = 1

# Descent-ascent duration from spawning to calyptopis I. Thorpe et al. (2019)
# give 23-26 days; the midpoint is applied as a constant. Varying it across
# that range changes the domain-wide M1 fraction by less than one percentage
# point — see krico-paper1/S1_m1_offset_sensitivity/.
SPAWNING_OFFSET_DAYS = 24

GLORYS_ICE_PATTERN = "glorys12_ice_{year:04d}_{month:02d}.nc"

# Monthly ice files are opened lazily and cached: consecutive cohorts almost
# always fall in the same month.
_ICE_CACHE = {}


def _open_ice_month(glorys_dir, year, month):
    """Open a monthly GLORYS12 sea-ice file, with a small cache."""
    key = (str(glorys_dir), year, month)
    if key not in _ICE_CACHE:
        path = f"{glorys_dir}/" + GLORYS_ICE_PATTERN.format(year=year,
                                                            month=month)
        _ICE_CACHE[key] = xr.open_dataset(path)
        # Keep the cache small; a cohort never needs more than two months.
        if len(_ICE_CACHE) > 3:
            for stale in list(_ICE_CACHE):
                if stale != key:
                    _ICE_CACHE.pop(stale).close()
                    break
    return _ICE_CACHE[key]


def _grid_axes(ds):
    """
    Return ascending 1-D longitude and latitude axes from the 2-D NEMO
    coordinate arrays.

    The monthly files carry nav_lon and nav_lat on the ORCA grid. Over the
    KRICO domain that grid is regular, so a row and a column are valid axes,
    but this is checked rather than assumed.
    """
    lon = np.asarray(ds["nav_lon"].values)
    lat = np.asarray(ds["nav_lat"].values)

    if lon.ndim == 2:
        if not np.allclose(lon, lon[0:1, :], equal_nan=True):
            raise ValueError(
                "nav_lon varies along y: the grid is not regular over this "
                "domain and cannot be indexed as separable axes."
            )
        lon = lon[0, :]
    if lat.ndim == 2:
        if not np.allclose(lat, lat[:, 0:1], equal_nan=True):
            raise ValueError(
                "nav_lat varies along x: the grid is not regular over this "
                "domain and cannot be indexed as separable axes."
            )
        lat = lat[:, 0]

    if lon.size < 2 or lat.size < 2:
        raise ValueError("degenerate coordinate axes in the sea-ice file")
    if lon[1] < lon[0] or lat[1] < lat[0]:
        raise ValueError("coordinate axes are not ascending")
    return lon, lat


def _sample_nearest(field, lon_axis, lat_axis, lon, lat):
    """Nearest-neighbour sample of a regular grid at scattered points."""
    lon_mid = 0.5 * (lon_axis[:-1] + lon_axis[1:])
    lat_mid = 0.5 * (lat_axis[:-1] + lat_axis[1:])
    i = np.searchsorted(lon_mid, lon)
    j = np.searchsorted(lat_mid, lat)
    np.clip(i, 0, lon_axis.size - 1, out=i)
    np.clip(j, 0, lat_axis.size - 1, out=j)
    return field[j, i]


def spawning_sic(lon, lat, release_date, glorys_dir,
                 offset_days=SPAWNING_OFFSET_DAYS):
    """
    Sea-ice concentration at the spawning position and date.

    Parameters
    ----------
    lon, lat : ndarray of shape (n_particles,)
        Release positions, taken as the spawning positions.
    release_date : pandas.Timestamp or datetime-like
        Release date of the cohort (day 0 of the trajectory).
    glorys_dir : str or Path
        Directory holding glorys12_ice_YYYY_MM.nc.
    offset_days : int
        Days from spawning to calyptopis I. The field is sampled
        offset_days before release_date. Zero samples the release date
        itself, which reproduces the trajectory value at day 0 and is
        used as a regression test.

    Returns
    -------
    ndarray of shape (n_particles,)
        Concentration as a fraction (0-1), NaN over land and where the
        field carries its fill value.
    """
    release_date = pd.Timestamp(release_date).normalize()
    spawning_date = release_date - pd.Timedelta(days=int(offset_days))

    ds = _open_ice_month(glorys_dir, spawning_date.year, spawning_date.month)

    # Daily means are stamped at 12:00, so a midnight target would be
    # equidistant from two days. Ask for midday, then verify the day that
    # came back is the one requested.
    target = np.datetime64(spawning_date.strftime("%Y-%m-%dT12:00:00"))
    da = ds["ileadfra"].sel(time_counter=target, method="nearest")
    selected = pd.Timestamp(da["time_counter"].values).normalize()
    if selected != spawning_date:
        raise ValueError(
            f"requested sea ice for {spawning_date.date()}, field returned "
            f"{selected.date()}; check monthly file coverage in {glorys_dir}"
        )

    field = np.squeeze(np.asarray(da.values, dtype="float64"))
    if field.ndim != 2:
        raise ValueError(
            f"expected a 2-D sea-ice field, got shape {field.shape}"
        )

    lon_axis, lat_axis = _grid_axes(ds)
    return _sample_nearest(field, lon_axis, lat_axis,
                           np.asarray(lon, dtype="float64"),
                           np.asarray(lat, dtype="float64"))


def find_april1_day_index(release_date, n_obs):
    """
    Find the day index (0-based, into the trajectory) corresponding to
    April 1 of the relevant austral autumn following the release date.

    Parameters
    ----------
    release_date : pandas.Timestamp or datetime-like
        Release date (day 0 of the trajectory).
    n_obs : int
        Number of observation days in the trajectory (typically 201).

    Returns
    -------
    day_index : int
        Day index corresponding to April 1. If the relevant April 1
        falls outside the tracking window (should not happen with
        200+ day tracking and Nov-Mar releases), returns n_obs.
    """
    release_date = pd.Timestamp(release_date)

    # Austral-autumn April 1: same year if release month >= January
    # (i.e. release in Jan, Feb, Mar); else following year (Nov, Dec).
    if release_date.month >= 1 and release_date.month <= 3:
        april1_year = release_date.year
    else:  # Nov, Dec
        april1_year = release_date.year + 1

    april1 = pd.Timestamp(year=april1_year, month=4, day=1)
    day_index = (april1 - release_date.normalize()).days

    if day_index < 0 or day_index >= n_obs:
        return n_obs

    return int(day_index)


def detect_sea_ice_advance(sic, april1_day):
    """
    Detect sea-ice advance event per particle.

    Advance is the first day index (>= april1_day) such that sic exceeds
    ADVANCE_SIC_THRESHOLD for ADVANCE_MIN_CONSECUTIVE_DAYS consecutive days
    starting at that day.

    Parameters
    ----------
    sic : ndarray of shape (n_particles, n_obs)
        Sea ice concentration along trajectory (fraction, 0-1).
        NaN treated as "below threshold" (no advance detected for
        deleted particles).
    april1_day : int
        First day index from which detection is allowed.

    Returns
    -------
    advance_day : ndarray of shape (n_particles,), dtype int64
        Day index of the first day of the qualifying 5-day run.
        -1 where no advance event is detected within the trajectory.
    """
    n_particles, n_obs = sic.shape

    # Boolean: SIC above threshold (NaN -> False).
    with np.errstate(invalid="ignore"):
        above = np.nan_to_num(sic, nan=0.0) > ADVANCE_SIC_THRESHOLD

    advance_day = np.full(n_particles, -1, dtype=np.int64)

    # For efficiency, compute the length of the run starting at each day
    # by running backward once. run_len[:, t] = number of consecutive True
    # starting at day t.
    run_len = np.zeros((n_particles, n_obs), dtype=np.int32)
    run_len[:, -1] = above[:, -1].astype(np.int32)
    for t in range(n_obs - 2, -1, -1):
        run_len[:, t] = np.where(above[:, t], run_len[:, t + 1] + 1, 0)

    # We want the first t >= april1_day where run_len[:, t] >= min_days.
    # Because the run_len array already accounts for continuation, a single
    # day t with run_len[t] >= k means "starting at t, next k days all True".
    qualifies = run_len >= ADVANCE_MIN_CONSECUTIVE_DAYS

    # Mask out days before april1_day.
    if april1_day > 0:
        qualifies[:, :april1_day] = False

    any_qualifies = qualifies.any(axis=1)
    advance_day = np.where(any_qualifies, qualifies.argmax(axis=1), -1)

    return advance_day.astype(np.int64)


def check_censoring(sic, n_obs, advance_day, last_valid):
    """
    For particles with no detected sea-ice advance, identify which are
    "censored" (SIC rising at end of tracking, advance may be pending)
    versus "no advance" (SIC flat/low, advance not expected).

    A particle is considered censored if it is alive at end of tracking
    (last_valid within END_OF_TRACKING_TOLERANCE of the final day) AND its
    SIC at the last valid day is non-negligible AND trending upward over
    the last 30 days.

    Parameters
    ----------
    sic : ndarray of shape (n_particles, n_obs)
        Sea ice concentration along trajectory.
    n_obs : int
        Number of observation days.
    advance_day : ndarray of shape (n_particles,)
        Output of detect_sea_ice_advance. Censoring is only evaluated for
        particles with advance_day == -1.
    last_valid : ndarray of shape (n_particles,)
        Output of trajectory.last_valid_day. Used to identify particles
        that are still alive (or essentially alive, modulo a 1-day Parcels
        precision artifact) at the end of tracking.

    Returns
    -------
    censored : ndarray of shape (n_particles,), dtype bool
        True for particles flagged as censored.
    """
    n_particles = sic.shape[0]

    # Only candidates: those without a detected advance event.
    candidates = advance_day == -1

    # "Alive at end of tracking" — tolerant to a single NaN day at the very
    # end caused by Parcels time-interpolation precision.
    alive_at_end = last_valid >= n_obs - 1 - END_OF_TRACKING_TOLERANCE

    # SIC at the last valid day for each particle.
    safe_idx = np.where(last_valid >= 0, last_valid, 0)
    sic_last = sic[np.arange(n_particles), safe_idx]

    # Censoring heuristic: SIC rising in the last 30 days.
    # Compare mean of last 10 days to mean of 20-30 days before end.
    # nanmean tolerates the precision NaN at the very last day; particles
    # that exited the domain mid-trajectory have all-NaN tails and yield
    # NaN means, which fail the rising comparison below (NaN comparisons
    # are False).
    import warnings

    if n_obs >= 30:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            recent = np.nanmean(sic[:, -10:], axis=1)
            earlier = np.nanmean(sic[:, -30:-10], axis=1)
    else:
        # Trajectory too short to assess trend; fall back to a simple
        # "non-zero at end" flag.
        recent = sic_last
        earlier = np.zeros(n_particles)

    # Require recent SIC to be non-negligible AND higher than earlier.
    rising = (recent > 0.05) & (recent > earlier)

    censored = candidates & alive_at_end & rising
    return censored
