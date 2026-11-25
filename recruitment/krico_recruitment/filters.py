"""
Mortality rules M1, M3 and M4 of Thorpe et al. (2019), and the run-length
routines behind M4.

M1: Spawning cannot occur where sea ice concentration >= 80% at spawning.
M3: Embryos, nauplii and the metanauplius die if the mean water temperature
    is -1 degC or colder during either the descent or the ascent.
M4: Calyptopis stages (CI-CIII) starve after >10 consecutive days under
    sea ice concentration > 40%. The duration of CI-CIII is
    temperature-dependent and trajectory-specific (development module).

M1 and M3 act during the 24 days between spawning and release, which are
not simulated. Their inputs are read at the release cell, taken as the
spawning position (pre_release module); the functions here only apply the
thresholds.

M2 (hatching on the sea floor) holds by construction: embryos hatch at
1000 m or shallower and the release cells are 1000 m or deeper.

M5 and M6 are evaluated at each particle's own sea-ice advance, in
classify.outcome.
"""

import numpy as np


# M1 threshold: sea ice concentration at spawning (fraction, 0-1).
M1_SIC_THRESHOLD = 0.80

# M3 threshold: mean temperature (degC) during the descent or the ascent.
M3_TEMPERATURE_THRESHOLD = -1.0

# M4 parameters.
M4_SIC_THRESHOLD = 0.40           # fraction, 0-1
M4_MAX_CONSECUTIVE_DAYS = 10      # strictly greater than this triggers starvation


def evaluate_M1(sic_at_spawning, threshold=M1_SIC_THRESHOLD):
    """
    Evaluate M1 (spawning under compact sea ice) for all particles.

    Parameters
    ----------
    sic_at_spawning : ndarray of shape (n_particles,)
        Sea ice concentration (fraction, 0-1) at the release cell on the
        spawning date. NaN (land, or outside the field) is treated as below
        threshold.
    threshold : float

    Returns
    -------
    killed : ndarray of shape (n_particles,), dtype bool
        True where the particle is removed by M1.
    """
    # Compared in double precision, as since v2.0.0, whatever the type the
    # concentration was stored with.
    sic = np.nan_to_num(np.asarray(sic_at_spawning, dtype=np.float64), nan=0.0)
    return sic >= threshold


def evaluate_M3(descent_temperature, ascent_temperature,
                threshold=M3_TEMPERATURE_THRESHOLD):
    """
    Evaluate M3 (cold sensitivity before calyptopis I) for all particles.

    Parameters
    ----------
    descent_temperature, ascent_temperature : ndarray of shape (n_particles,)
        Mean temperature (degC) experienced during the descent and during
        the ascent (pre_release.evaluate). NaN (a cell with no value)
        removes nothing, as for M1.
    threshold : float

    Returns
    -------
    killed : ndarray of shape (n_particles,), dtype bool
        True where the colder of the two means is at or below the threshold.
    """
    coldest = np.fmin(np.asarray(descent_temperature, dtype=np.float64),
                      np.asarray(ascent_temperature, dtype=np.float64))
    with np.errstate(invalid="ignore"):
        return np.nan_to_num(coldest, nan=np.inf) <= threshold


def under_ice(sic, threshold):
    """
    Boolean (n_particles, n_obs): concentration strictly above threshold.
    NaN (after deletion) counts as not under ice. The comparison is made in
    the precision of sic, as since v2.0.0.
    """
    # float() keeps the threshold a Python scalar, which numpy compares in
    # the precision of the array; a numpy float64 would promote it.
    with np.errstate(invalid="ignore"):
        return np.nan_to_num(sic, nan=0.0) > float(threshold)


def run_lengths(flags):
    """
    Length of the run of consecutive True ending at each day.

    Parameters
    ----------
    flags : ndarray of shape (n_particles, n_obs), dtype bool

    Returns
    -------
    ndarray of shape (n_particles, n_obs), dtype int16
        0 where the flag is False; otherwise the number of consecutive True
        days up to and including this one.
    """
    n_particles, n_obs = flags.shape
    out = np.zeros((n_particles, n_obs), dtype=np.int16)
    running = np.zeros(n_particles, dtype=np.int16)
    for t in range(n_obs):
        running = np.where(flags[:, t], running + 1, 0).astype(np.int16)
        out[:, t] = running
    return out


def longest_run(runs, end_day):
    """
    Longest run within the window of days 0 to end_day (exclusive).

    Parameters
    ----------
    runs : ndarray of shape (n_particles, n_obs)
        Output of run_lengths.
    end_day : ndarray of shape (n_particles,)
        First day outside each particle's window; n_obs or more for the
        whole track.

    Returns
    -------
    ndarray of shape (n_particles,), dtype int16
        0 where the window is empty or holds no True day.
    """
    n_obs = runs.shape[1]
    inside = np.arange(n_obs)[None, :] < np.asarray(end_day)[:, None]
    return np.where(inside, runs, 0).max(axis=1).astype(np.int16)


def first_long_run_day(runs, max_days):
    """
    First day on which a run becomes longer than max_days, over the whole
    track: the day a particle would starve if it were still a calyptopis.
    -1 where no run is that long.
    """
    crossed = runs > max_days
    found = crossed.any(axis=1)
    return np.where(found, crossed.argmax(axis=1), -1).astype(np.int64)


def evaluate_M4(sic, calyptopis_end_day, threshold=M4_SIC_THRESHOLD,
                max_days=M4_MAX_CONSECUTIVE_DAYS):
    """
    Evaluate M4 (calyptopis starvation under sea ice) for all particles.

    A particle is removed if, within its calyptopis window (day 0 to
    calyptopis_end_day exclusive), it spends more than max_days consecutive
    days with concentration above threshold.

    Parameters
    ----------
    sic : ndarray of shape (n_particles, n_obs)
        Sea ice concentration along trajectory (fraction, 0-1).
    calyptopis_end_day : ndarray of shape (n_particles,), dtype int64
        Day index at which CIII completes
        (development.calyptopis_window_end).

    Returns
    -------
    killed : ndarray of shape (n_particles,), dtype bool
    kill_day : ndarray of shape (n_particles,), dtype int64
        Day on which the run became longer than max_days; -1 where M4 did
        not act.
    """
    day = first_long_run_day(run_lengths(under_ice(sic, threshold)), max_days)
    killed = (day >= 0) & (day < np.asarray(calyptopis_end_day))
    return killed, np.where(killed, day, -1).astype(np.int64)
