"""
Temperature-dependent development model for Antarctic krill larvae.

Implements Thorpe et al. (2019) Eq. 2 with stage-specific coefficients
from Table 1. Integrates fractional development along trajectory using
daily-mean temperature at the particle's actual depth.

Particles begin at stage CI (calyptopis I), consistent with the PoC
assumption that release at 50-200 m represents post-ascent larvae.
Naupliar and metanaupliar stages (NI, NII, MN) are skipped.

References
----------
Thorpe, S.E., Tarling, G.A., Murphy, E.J. (2019). Circumpolar patterns in
    Antarctic krill larval recruitment: an environmentally driven model.
    Marine Ecology Progress Series 613, 77-96.
"""

import numpy as np


# Exponential decay coefficient (shared across stages).
# From Thorpe (2019) Section 2.1.2: f = y0 + a * exp(-b * T), with b = 1.0479.
B = 1.0479

# Stage-specific coefficients from Thorpe (2019) Table 1.
# Only stages from CI onward are needed (particles start post-ascent).
# Order: CI, CII, CIII, FI, FII, FIII, FIV, FV, FVI.
STAGES = ("CI", "CII", "CIII", "FI", "FII", "FIII", "FIV", "FV", "FVI")

Y0 = np.array([
    11.0449,  # CI
    7.1003,   # CII
    7.8892,   # CIII
    9.4671,   # FI
    8.6871,   # FII
    10.2560,  # FIII
    10.2560,  # FIV
    11.8338,  # FV
    9.4671,   # FVI
])

A = np.array([
    2.9568,   # CI
    1.9008,   # CII
    2.1120,   # CIII
    2.5344,   # FI
    2.3232,   # FII
    2.7456,   # FIII
    2.7456,   # FIV
    3.1680,   # FV
    2.5344,   # FVI
])

# Index of FIV within STAGES (for the "reached FIV" threshold).
# A particle has reached FIV when cumulative fractional development
# has crossed the sum of durations for CI, CII, CIII, FI, FII, FIII.
FIV_INDEX = STAGES.index("FIV")


def stage_duration(temperature, stage_index):
    """
    Compute development duration (days) for a single stage at given temperature.

    Uses Thorpe (2019) Eq. 2:  duration = y0 + a * exp(-b * T)

    Parameters
    ----------
    temperature : float or ndarray
        Temperature in degrees Celsius.
    stage_index : int
        Index into STAGES (0 = CI, 1 = CII, ..., 8 = FVI).

    Returns
    -------
    float or ndarray
        Stage duration in days (same shape as temperature).
    """
    return Y0[stage_index] + A[stage_index] * np.exp(-B * temperature)


def integrate_development(temperature):
    """
    Integrate fractional development along trajectory for many particles.

    At each daily timestep, each particle's current stage is identified from
    its cumulative fractional development, and the stage-specific duration is
    computed from that day's temperature. The increment 1/duration is added
    to the cumulative fractional development.

    Parameters
    ----------
    temperature : ndarray of shape (n_particles, n_obs)
        Daily-mean temperature (°C) at the particle's actual depth.
        May contain NaN for deleted particles (post-deletion positions).

    Returns
    -------
    cumulative_dev : ndarray of shape (n_particles, n_obs)
        Cumulative fractional development, expressed in "stage units".
        A value of k means the particle has completed k stages
        (starting count from CI = stage 0). Particle has reached stage
        STAGES[i] when cumulative_dev >= i.
        NaN where temperature is NaN (propagated).
    """
    n_particles, n_obs = temperature.shape
    cumulative_dev = np.zeros((n_particles, n_obs), dtype=np.float32)

    # Running cumulative development (1D, one value per particle).
    running = np.zeros(n_particles, dtype=np.float32)

    # We need stage-duration-at-current-temperature for the stage each
    # particle is currently in. The current stage is int(running), clipped
    # to valid range.

    # Pre-compute durations for all stages at all (particle, obs) pairs?
    # That would allocate (n_stages, n_particles, n_obs) which could be large.
    # Instead, loop over days and use fancy indexing: per-particle current
    # stage picks the right y0/a.

    n_stages = len(STAGES)

    for t in range(n_obs):
        T_t = temperature[:, t]  # shape (n_particles,)

        # Current stage index per particle, clipped to [0, n_stages-1].
        # Particles past FVI stop developing (their stage index saturates).
        current_stage = np.clip(
            np.floor(running).astype(np.int64), 0, n_stages - 1
        )

        # Pick y0 and a for each particle's current stage.
        y0_t = Y0[current_stage]
        a_t = A[current_stage]

        # Stage duration at today's temperature for each particle.
        duration_t = y0_t + a_t * np.exp(-B * T_t)

        # Daily increment: 1 day / stage duration.
        # Handle NaN temperatures: increment becomes NaN, so running becomes
        # NaN and stays NaN. We don't want that — use masked increment.
        with np.errstate(invalid="ignore"):
            increment = np.where(
                np.isnan(T_t),
                0.0,  # no development on missing days (particle deleted)
                1.0 / duration_t,
            )

        # Saturate at n_stages: once past FVI, no more development tracked.
        running = np.where(
            running >= n_stages,
            running,
            running + increment,
        )

        cumulative_dev[:, t] = running

    # Propagate NaN where temperature was NaN (post-deletion).
    cumulative_dev = np.where(np.isnan(temperature), np.nan, cumulative_dev)

    return cumulative_dev


def days_to_reach_stage(cumulative_dev, stage_index):
    """
    For each particle, find the first day index at which cumulative fractional
    development reaches or exceeds stage_index.

    Parameters
    ----------
    cumulative_dev : ndarray of shape (n_particles, n_obs)
        Cumulative fractional development from `integrate_development`.
    stage_index : int
        Target stage index (e.g. FIV_INDEX for "reached FIV").

    Returns
    -------
    days : ndarray of shape (n_particles,), dtype int64
        Day index at which the stage was reached. -1 if never reached
        within the trajectory.
    """
    reached = cumulative_dev >= stage_index
    # For each particle, first True along axis=1.
    any_reached = reached.any(axis=1)
    first_day = np.where(any_reached, reached.argmax(axis=1), -1)
    return first_day.astype(np.int64)


def calyptope_window_end(cumulative_dev):
    """
    For each particle, find the day index at which it completes CIII
    (cumulative_dev reaches 3, i.e. transitions from CIII to FI).

    This defines the end of the CI-CIII calyptope window used by M4.

    Parameters
    ----------
    cumulative_dev : ndarray of shape (n_particles, n_obs)

    Returns
    -------
    day_end : ndarray of shape (n_particles,), dtype int64
        Day index at which CIII completes. If the particle never finishes
        CIII within the tracking period, returns n_obs (i.e., calyptope
        window extends to end of tracking).
    """
    # Stage index 3 = FI, so reaching 3 means CIII just completed.
    n_obs = cumulative_dev.shape[1]
    reached = cumulative_dev >= 3.0
    any_reached = reached.any(axis=1)
    day_end = np.where(any_reached, reached.argmax(axis=1), n_obs)
    return day_end.astype(np.int64)
