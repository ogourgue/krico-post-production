"""
Evaluation of mortality filters M1 and M4 from Thorpe et al. (2019).

M1: Spawning cannot occur where sea ice concentration >= 80% at release.
M4: Calyptope stages (CI-CIII) starve after >10 consecutive days under
    sea ice concentration > 40%. The duration of CI-CIII is
    temperature-dependent and trajectory-specific (from development module).

M2 (seafloor hatching) is handled implicitly via the 1000-2000 m bathymetric
release zone and is not evaluated here.

M3 (cold sensitivity during descent/ascent) is omitted because particles
begin post-ascent at 50-200 m depth.
"""

import numpy as np


# M1 threshold: sea ice concentration at release (fraction, 0-1).
M1_SIC_THRESHOLD = 0.80

# M4 parameters.
M4_SIC_THRESHOLD = 0.40           # fraction, 0-1
M4_MAX_CONSECUTIVE_DAYS = 10      # strictly greater than this triggers starvation


def evaluate_M1(sic_at_release):
    """
    Evaluate M1 (spawning under dense sea ice) for all particles.

    Parameters
    ----------
    sic_at_release : ndarray of shape (n_particles,)
        Sea ice concentration (fraction, 0-1) at the particle's release
        position at release time (day 0).

    Returns
    -------
    killed : ndarray of shape (n_particles,), dtype bool
        True where the particle is killed by M1.
    """
    return sic_at_release >= M1_SIC_THRESHOLD


def evaluate_M4(sic, calyptope_end_day):
    """
    Evaluate M4 (calyptope starvation under sea ice) for all particles.

    For each particle, examine sea ice concentration along the trajectory
    during its calyptope window (day 0 to calyptope_end_day exclusive).
    If at any point in that window the particle experiences more than
    M4_MAX_CONSECUTIVE_DAYS consecutive days with SIC > M4_SIC_THRESHOLD,
    it is killed by M4.

    Parameters
    ----------
    sic : ndarray of shape (n_particles, n_obs)
        Sea ice concentration along trajectory (fraction, 0-1).
    calyptope_end_day : ndarray of shape (n_particles,), dtype int64
        For each particle, day index at which CIII completes
        (from development.calyptope_window_end).

    Returns
    -------
    killed : ndarray of shape (n_particles,), dtype bool
        True where the particle is killed by M4.
    kill_day : ndarray of shape (n_particles,), dtype int64
        Day index at which the M4 threshold was crossed. -1 where M4
        did not trigger.
    """
    n_particles, n_obs = sic.shape
    killed = np.zeros(n_particles, dtype=bool)
    kill_day = np.full(n_particles, -1, dtype=np.int64)

    # Boolean array: particle under heavy ice on this day.
    # NaN in SIC (after deletion) treated as False here — doesn't matter
    # because those particles will be handled at outcome-determination time.
    with np.errstate(invalid="ignore"):
        under_ice = np.nan_to_num(sic, nan=0.0) > M4_SIC_THRESHOLD

    # Vectorized consecutive-run detection per particle would be elegant
    # but tricky. Loop over particles — n_particles is ~500k, but each
    # particle's inner operation is cheap (numpy on 1D array of length n_obs).
    # For now this is O(n_particles * n_obs) which is fine (~1e8 ops).
    #
    # TODO: if profiling shows this is a bottleneck, consider a cumulative-
    # reset trick for fully vectorized run-length on 2D arrays.

    for i in range(n_particles):
        end = calyptope_end_day[i]
        if end <= 0:
            continue  # calyptope window is empty (already past CIII at t=0)

        window = under_ice[i, :end]
        if not window.any():
            continue

        # Count consecutive True runs. Any run > M4_MAX_CONSECUTIVE_DAYS
        # triggers M4.
        run_length = 0
        for t in range(end):
            if window[t]:
                run_length += 1
                if run_length > M4_MAX_CONSECUTIVE_DAYS:
                    killed[i] = True
                    kill_day[i] = t
                    break
            else:
                run_length = 0

    return killed, kill_day
