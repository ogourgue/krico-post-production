"""
The classification, in two steps that never touch a file:

    diagnostics   what each particle experienced, from its trajectory
                  (trajectory_diagnostics) and from its release cell
                  (pre_release.evaluate)
    outcome       the state of each particle, from the diagnostics alone
                  and a set of parameters (outcome)

The diagnostics are what the dataset stores, for every particle, including
those removed by M1 or M3. The outcome of the dataset is outcome() at the
baseline parameters, computed from the diagnostics in their stored form
(to_stored), so that reading a file and calling outcome() returns the
stored outcome exactly, and any other threshold can be evaluated the same
way without the trajectories.

Order of the rules: M1, M3, M4, then the state at the particle's own
sea-ice advance (recruit on or off the shelf/slope, or M5), then M6,
censoring and domain exit for particles with no advance.

What the stored diagnostics can vary, one parameter at a time:
    m1_sic            any value
    m3_temperature    any value; None switches M3 off
    m3_hatching_depth one of the stored depths (530 or 1000 m)
    m4_sic            one of the stored thresholds (0.20 to 0.60)
    m4_days           any value
    stage             any stage index (FIV = 6), exactly
    shelf_depth       any value
    rate              any factor from 0.633, with m4_sic and m4_days at
                      their baseline (see development.apply_rate)
The definition of the advance (0.15, 5 days, from 1 April) is baked into
advance_day; ice15_daily and shelf2000_daily allow another definition to
be evaluated outside this function.
"""

from dataclasses import dataclass

import numpy as np

from krico_recruitment import development, filters, sea_ice, trajectory
from krico_recruitment.outcome import OUTCOME_CODES

# Stored form of a day index: uint8, DAY_NONE where there is none.
DAY_NONE = 255

# Sea-ice thresholds at which the longest calyptopis run is stored.
M4_ICE_THRESHOLDS = (0.20, 0.30, 0.40, 0.50, 0.60)

# Sea-floor depth (m) separating the two classes of recruits.
SHELF_DEPTH = 2000.0

# Lower limit (degC) of the range over which the development function was
# fitted (-1 to +4 degC).
FIT_LOWER_LIMIT = -1.0


@dataclass(frozen=True)
class Parameters:
    """Parameters of outcome(). The defaults are the baseline of v3.0.0."""
    m1_sic: float = filters.M1_SIC_THRESHOLD
    m3_temperature: float | None = filters.M3_TEMPERATURE_THRESHOLD
    m3_hatching_depth: float = 1000.0
    m4_sic: float = filters.M4_SIC_THRESHOLD
    m4_days: int = filters.M4_MAX_CONSECUTIVE_DAYS
    stage: int = development.FIV_INDEX
    shelf_depth: float = SHELF_DEPTH
    rate: float = 1.0


BASELINE = Parameters()


def _day(values, none):
    """Day indices as uint8, with DAY_NONE where values == none."""
    values = np.asarray(values)
    if values.size and (values.max() >= DAY_NONE or values.min() < -1):
        raise ValueError("day index outside the range of the stored type")
    return np.where(values == none, DAY_NONE, values).astype(np.uint8)


def trajectory_diagnostics(lon, lat, depth, temperature, bathymetry, sic,
                           release_date):
    """
    Diagnostics of every particle from its trajectory, in stored form.

    Parameters
    ----------
    lon, lat, depth, temperature, bathymetry, sic :
        ndarray of shape (n_particles, n_obs), daily, NaN after deletion.
    release_date : pandas.Timestamp or datetime-like

    Returns
    -------
    diagnostics : dict of arrays over the particle dimension
    last_valid : ndarray, int64
        Last day with a valid position (trajectory.last_valid_day); needed
        for the position at which fate is recorded, not stored.
    """
    n_particles, n_obs = lon.shape
    rows = np.arange(n_particles)
    last_valid = trajectory.last_valid_day(lon, lat, depth)
    alive_at_end = last_valid >= n_obs - 1 - sea_ice.END_OF_TRACKING_TOLERANCE

    # Development.
    dev = development.integrate_development(temperature)
    calyptopis_end = development.calyptopis_window_end(dev)      # n_obs: never
    fiv_day = development.days_to_reach_stage(dev, development.FIV_INDEX)

    # Sea-ice advance. An advance detected after the particle was deleted
    # was never experienced.
    opening = sea_ice.find_april1_day_index(release_date, n_obs)
    advance = sea_ice.detect_sea_ice_advance(sic, opening)
    advance = np.where((advance != -1) & (advance > last_valid), -1, advance)
    has_advance = advance != -1

    # Calyptopis starvation: run lengths under ice.
    longest = np.zeros((n_particles, len(M4_ICE_THRESHOLDS)), dtype=np.uint8)
    for k, threshold in enumerate(M4_ICE_THRESHOLDS):
        runs = filters.run_lengths(filters.under_ice(sic, threshold))
        longest[:, k] = filters.longest_run(runs, calyptopis_end)
        if threshold == filters.M4_SIC_THRESHOLD:
            first_run = filters.first_long_run_day(
                runs, filters.M4_MAX_CONSECUTIVE_DAYS)
    del runs

    # Days below the fitted range, up to the advance or the end of tracking.
    with np.errstate(invalid="ignore"):
        below = np.cumsum(temperature < FIT_LOWER_LIMIT, axis=1, dtype=np.int16)
    days_below = below[rows, np.where(has_advance, advance, n_obs - 1)]
    del below

    with np.errstate(invalid="ignore"):
        ice15 = np.nan_to_num(sic, nan=0.0) > sea_ice.ADVANCE_SIC_THRESHOLD
        shelf = np.nan_to_num(bathymetry, nan=np.inf) < SHELF_DEPTH

    diagnostics = {
        "m4_longest_run": longest,
        "m4_first_run_day": _day(first_run, -1),
        "m4_first_run_dev": development.quantise(
            trajectory.extract_at_day(dev, first_run)),
        "calyptopis_end_day": _day(calyptopis_end, n_obs),
        "fiv_day": _day(fiv_day, -1),
        "advance_day": _day(advance, -1),
        "advance_dev": development.quantise(
            trajectory.extract_at_day(dev, advance)),
        "advance_bathymetry": trajectory.extract_at_day(
            bathymetry, advance).astype(np.float32),
        "advance_lon": trajectory.extract_at_day(lon, advance).astype(np.float32),
        "advance_lat": trajectory.extract_at_day(lat, advance).astype(np.float32),
        "alive_at_end": alive_at_end.astype(np.int8),
        "sic_rising_at_end": sea_ice.rising_at_end(
            sic, n_obs, last_valid).astype(np.int8),
        "release_bathymetry": bathymetry[:, 0].astype(np.float32),
        "days_below_fit": days_below.astype(np.uint8),
        "ice15_daily": np.packbits(ice15, axis=1),
        "shelf2000_daily": np.packbits(shelf, axis=1),
    }
    return diagnostics, last_valid


def pre_release_diagnostics(values, cell):
    """
    Stored form of the output of pre_release.evaluate, with the cell index.
    """
    mean = values["spawning_temperature_mean"]
    return {
        "release_cell": np.asarray(cell, dtype=np.int32),
        "sic_spawning": values["sic_spawning"].astype(np.float32),
        "m3_descent_temperature": values["m3_descent_temperature"].astype(np.float32),
        "m3_ascent_temperature": values["m3_ascent_temperature"].astype(np.float32),
        "spawning_temperature_surface":
            values["spawning_temperature_surface"].astype(np.float32),
        "spawning_temperature_mean_100m": mean[:, 0].astype(np.float32),
        "spawning_temperature_mean_200m": mean[:, 1].astype(np.float32),
    }


def _index(value, stored, what):
    match = np.flatnonzero(np.isclose(np.asarray(stored, dtype=np.float64), value))
    if match.size != 1:
        raise ValueError(f"{what} {value} is not stored; choose among "
                         f"{[float(v) for v in stored]}")
    return int(match[0])


def masks(d, parameters=BASELINE, hatching_depths=(530.0, 1000.0),
          ice_thresholds=M4_ICE_THRESHOLDS):
    """
    Boolean mask of each outcome, from the stored diagnostics d.

    d holds the variables of a v3.0.0 file as raw arrays (not masked, not
    scaled). hatching_depths and ice_thresholds are the coordinates of the
    two-dimensional diagnostics. Returns a dict keyed by outcome label;
    the masks are mutually exclusive and cover every particle.
    """
    p = parameters
    killed_M1 = filters.evaluate_M1(d["sic_spawning"], p.m1_sic)

    if p.m3_temperature is None:
        killed_M3 = np.zeros(killed_M1.shape, dtype=bool)
    else:
        h = _index(p.m3_hatching_depth, hatching_depths, "hatching depth")
        killed_M3 = ~killed_M1 & filters.evaluate_M3(
            d["m3_descent_temperature"][:, h],
            d["m3_ascent_temperature"][:, h], p.m3_temperature)

    advance_day = np.asarray(d["advance_day"])
    advance_dev = np.asarray(d["advance_dev"])
    has_advance = advance_day != DAY_NONE
    has_dev = advance_dev != development.DEV_MISSING

    if p.rate == 1.0:
        k = _index(p.m4_sic, ice_thresholds, "M4 sea-ice threshold")
        starved = np.asarray(d["m4_longest_run"])[:, k] > p.m4_days
        # A development missing at the advance (no temperature that day)
        # counts as below the stage, as since v2.0.0.
        reached = has_dev & (advance_dev.astype(np.int64)
                             >= p.stage * development.DEV_SCALE)
    else:
        if (not np.isclose(p.m4_sic, filters.M4_SIC_THRESHOLD)
                or p.m4_days != filters.M4_MAX_CONSECUTIVE_DAYS):
            raise ValueError("a development rate other than 1 can only be "
                             "combined with the baseline M4 threshold and "
                             "duration: the files store the first run of "
                             "more than 10 days above 0.40")
        run_day = np.asarray(d["m4_first_run_day"]).astype(np.int64)
        run_dev = np.asarray(d["m4_first_run_dev"])
        end = np.asarray(d["calyptopis_end_day"]).astype(np.int64)
        scaled = development.apply_rate(development.dequantise(run_dev), p.rate)
        with np.errstate(invalid="ignore"):
            # Still a calyptopis when the run becomes too long. Where the
            # development of that day is missing, the baseline answer.
            starved = np.where(run_dev != development.DEV_MISSING,
                               scaled < 3.0, run_day < end)
        starved &= run_day != DAY_NONE
        scaled = development.apply_rate(development.dequantise(advance_dev), p.rate)
        with np.errstate(invalid="ignore"):
            reached = has_dev & (scaled >= p.stage)
    killed_M4 = ~killed_M1 & ~killed_M3 & starved

    alive = ~killed_M1 & ~killed_M3 & ~killed_M4
    alive_at_end = np.asarray(d["alive_at_end"]) == 1
    rising = np.asarray(d["sic_rising_at_end"]) == 1
    with np.errstate(invalid="ignore"):
        on_shelf = np.nan_to_num(np.asarray(d["advance_bathymetry"],
                                            dtype=np.float64),
                                 nan=np.inf) < p.shelf_depth

    at_advance = alive & has_advance
    no_advance = alive & ~has_advance
    return {
        "recruit_on_shelf_slope": at_advance & reached & on_shelf,
        "recruit_off_shelf_slope": at_advance & reached & ~on_shelf,
        "M1_spawning_under_ice": killed_M1,
        "M3_cold_sensitivity": killed_M3,
        "M4_calyptopis_starvation": killed_M4,
        "M5_overwintering_stage": at_advance & ~reached,
        "M6_overwintering_habitat": no_advance & alive_at_end & ~rising,
        "censored": no_advance & alive_at_end & rising,
        "exited_domain": no_advance & ~alive_at_end,
    }


def outcome(d, parameters=BASELINE, **coordinates):
    """
    Outcome code (int8, see krico_recruitment.outcome) of every particle
    from the stored diagnostics d. coordinates are passed to masks().
    """
    selected = masks(d, parameters, **coordinates)
    code = np.full(next(iter(selected.values())).shape, -1, dtype=np.int8)
    count = np.zeros(code.shape, dtype=np.int8)
    for label, mask in selected.items():
        code[mask] = OUTCOME_CODES[label]
        count += mask
    if not (count == 1).all():
        raise RuntimeError(f"{int((count != 1).sum())} particles are not in "
                           "exactly one outcome; error in the rules")
    return code


def fate_day(d, code, last_valid):
    """
    Day at which the fate of each particle is recorded, at the baseline:
    day 0 for M1 and M3 (the rules act before release, at the release
    position); the day the starvation run becomes too long for M4; the
    advance for the recruits and M5; the last valid day otherwise.
    """
    day = np.asarray(last_valid, dtype=np.int64).copy()
    before = np.isin(code, (OUTCOME_CODES["M1_spawning_under_ice"],
                            OUTCOME_CODES["M3_cold_sensitivity"]))
    day[before] = 0
    starved = code == OUTCOME_CODES["M4_calyptopis_starvation"]
    day[starved] = np.asarray(d["m4_first_run_day"])[starved]
    at_advance = np.isin(code, (OUTCOME_CODES["recruit_on_shelf_slope"],
                                OUTCOME_CODES["recruit_off_shelf_slope"],
                                OUTCOME_CODES["M5_overwintering_stage"]))
    day[at_advance] = np.asarray(d["advance_day"])[at_advance]
    return day
