"""
Sea-ice advance detection (per-particle winter onset).

For each particle, winter onset is defined as the first date on or after
April 1 of the relevant austral autumn where sea ice concentration at the
particle's location exceeds 15% for 5 consecutive days.

The April 1 cutoff avoids spurious detection from lingering previous-winter
ice. The austral autumn is defined relative to the release date:
  - release in Jan-Mar → April 1 of same calendar year
  - release in Nov-Dec → April 1 of following calendar year
Equivalently: the first April 1 that falls within the tracking window.
"""

import numpy as np
import pandas as pd


# Sea-ice advance parameters (Thorpe 2019, Stammerjohn et al. 2008).
ADVANCE_SIC_THRESHOLD = 0.15       # fraction, 0-1; strictly greater than this
ADVANCE_MIN_CONSECUTIVE_DAYS = 5   # SIC above threshold for this many days

# Tolerance for "alive at end of tracking". Parcels can produce a NaN
# position on the very last output step due to time-interpolation precision
# at the simulation end, even for particles that survived. We therefore
# treat last_valid in [n_obs - 1 - END_OF_TRACKING_TOLERANCE, n_obs - 1]
# as "alive at end".
END_OF_TRACKING_TOLERANCE = 1


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