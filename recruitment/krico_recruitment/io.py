"""
NetCDF I/O for the recruitment pipeline.

- read_trajectory_file: load one Parcels trajectory NetCDF into xarray.
- write_classification_file: save the per-particle outcome and diagnostics
  of one cohort (dataset v3.0.0).
- read_diagnostics: read them back in the form classify.outcome expects.

Layout of a v3.0.0 file
-----------------------
Every variable is on the particle dimension; three diagnostics have a
second one (hatching_depth, ice_threshold, packed_day).

  The nine variables of v2.0.0, with the same names and types:
    release_lon, release_lat, release_depth, outcome, final_lon, final_lat,
    final_depth, travel_time, trajectory_path_length

  Before release (pre_release module, at the release cell):
    release_cell, sic_spawning, m3_descent_temperature,
    m3_ascent_temperature, spawning_temperature_surface,
    spawning_temperature_mean_100m, spawning_temperature_mean_200m

  Along the trajectory (classify.trajectory_diagnostics):
    m4_longest_run, m4_first_run_day, m4_first_run_dev, calyptopis_end_day,
    fiv_day, advance_day, advance_dev, advance_bathymetry, advance_lon,
    advance_lat, alive_at_end, sic_rising_at_end, release_bathymetry,
    days_below_fit, ice15_daily, shelf2000_daily

Days are counted from release (day 0) and stored as uint8, 255 where there
is none. Development is stored as uint16 thousandths of a stage counted
from the start of calyptopis I, rounded down, 65535 where there is none.
The two daily flags hold one bit per day, 8 days per byte, most significant
bit first: numpy.unpackbits(x, axis=1)[:, :n_days].

Each diagnostic carries as attributes the parameters that are baked into
it and cannot be varied from the file; the outcome variable carries the
baseline parameters.
"""

import os
from pathlib import Path

import netCDF4
import numpy as np
import xarray as xr

from krico_recruitment import (classify, development, filters, pre_release,
                               sea_ice)
from krico_recruitment.outcome import OUTCOME_MEANINGS

DATASET_VERSION = "3.0.0"
COMPLEVEL = 5


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



_LON = {"standard_name": "longitude", "units": "degrees_east"}
_LAT = {"standard_name": "latitude", "units": "degrees_north"}
_DEPTH = {"standard_name": "depth", "units": "m", "positive": "down"}
# "day", not "days" or "days since release": xarray decodes the first as a
# time interval and fails to open a file holding the second. travel_time
# keeps the "days" it has had since v2.0.0.
_DAY = {"units": "day", "_FillValue": np.uint8(classify.DAY_NONE),
        "comment": "day index along the trajectory: 0 is the release day"}
_DEV = {"units": "0.001 stage",
        "_FillValue": np.uint16(development.DEV_MISSING),
        "comment": "cumulative development from the start of calyptopis I, "
                   "rounded down: 3000 is the start of furcilia I, 6000 the "
                   "start of furcilia IV, 9000 the end of furcilia VI, where "
                   "it saturates"}
_NAN = {"_FillValue": np.float32(np.nan)}
_FLAG = {"flag_values": np.array([0, 1], dtype=np.int8)}


def _schema(n_days):
    """Variables of a v3.0.0 file, in order: name -> (type, dims, attributes)."""
    offset = sea_ice.SPAWNING_OFFSET_DAYS
    descent = pre_release.DESCENT_DAYS
    column = ("thickness-weighted mean over the model levels from the "
              "surface to the hatching depth, cut at the sea floor of the "
              "cell, then over the days")
    packed = ("one bit per day from release, 8 days per byte, most "
              f"significant bit first: numpy.unpackbits(x, axis=1)[:, :{n_days}]")
    p = ("particle",)
    return {
        # ------------------------------------------------ as in v2.0.0
        "release_lon": ("f4", p, {**_LON, **_NAN, "long_name": "release longitude"}),
        "release_lat": ("f4", p, {**_LAT, **_NAN, "long_name": "release latitude"}),
        "release_depth": ("f4", p, {**_DEPTH, **_NAN, "long_name": "release depth"}),
        "outcome": ("i1", p, {
            "long_name": "recruitment outcome",
            "flag_values": np.arange(len(OUTCOME_MEANINGS), dtype=np.int8),
            "flag_meanings": " ".join(OUTCOME_MEANINGS),
            "comment": "krico_recruitment.classify.outcome at the baseline "
                       "parameters below, from the diagnostics of this file. "
                       "Codes 0 and 1 are the recruits; the codes differ "
                       "from those of v2.0.0",
            "m1_sic_threshold": filters.M1_SIC_THRESHOLD,
            "m3_temperature_threshold": filters.M3_TEMPERATURE_THRESHOLD,
            "m3_hatching_depth": pre_release.BASELINE_HATCHING_DEPTH,
            "m4_sic_threshold": filters.M4_SIC_THRESHOLD,
            "m4_consecutive_days": np.int32(filters.M4_MAX_CONSECUTIVE_DAYS),
            "overwintering_stage": development.STAGES[development.FIV_INDEX],
            "shelf_slope_depth": classify.SHELF_DEPTH,
            "development_rate": 1.0,
        }),
        "final_lon": ("f4", p, {**_LON, **_NAN,
                                "long_name": "longitude at the moment fate is determined"}),
        "final_lat": ("f4", p, {**_LAT, **_NAN,
                                "long_name": "latitude at the moment fate is determined"}),
        "final_depth": ("f4", p, {**_DEPTH, **_NAN,
                                  "long_name": "depth at the moment fate is determined"}),
        "travel_time": ("i4", p, {
            "long_name": "days from release to moment fate is determined",
            "units": "days"}),
        "trajectory_path_length": ("f4", p, {
            **_NAN, "units": "m",
            "long_name": "cumulative along-track path length up to moment "
                         "fate is determined"}),
        # ------------------------------------------------ before release
        "release_cell": ("i4", p, {
            "long_name": "release cell: index along the cell dimension of "
                         "the release-cell file of the spawning year "
                         "(glorys12_release_cells_YYYY.nc), the grid point "
                         "nearest to the release position"}),
        "sic_spawning": ("f4", p, {
            **_NAN, "units": "1",
            "long_name": "sea-ice concentration at the release cell on the "
                         "spawning date",
            "spawning_offset_days": np.int32(offset)}),
        "m3_descent_temperature": ("f4", ("particle", "hatching_depth"), {
            **_NAN, "units": "degC",
            "long_name": "mean temperature of the embryo during the descent, "
                         "at the release cell, under steady sinking",
            "cell_methods_comment": column,
            "first_day_after_spawning": np.int32(0),
            "last_day_after_spawning": np.int32(descent - 1)}),
        "m3_ascent_temperature": ("f4", ("particle", "hatching_depth"), {
            **_NAN, "units": "degC",
            "long_name": "mean temperature of the water column between the "
                         "hatching depth and the surface during the ascent, "
                         "at the release cell",
            "cell_methods_comment": column,
            "first_day_after_spawning": np.int32(descent),
            "last_day_after_spawning": np.int32(offset - 1)}),
        "spawning_temperature_surface": ("f4", p, {
            **_NAN, "units": "degC",
            "long_name": "temperature of the first model level at the "
                         "release cell on the spawning date",
            "spawning_offset_days": np.int32(offset)}),
        **{f"spawning_temperature_mean_{int(z)}m": ("f4", p, {
            **_NAN, "units": "degC",
            "long_name": f"mean temperature of the top {int(z)} m at the "
                         "release cell on the spawning date",
            "cell_methods_comment": "thickness-weighted mean over the model "
                                    "levels, cut at the sea floor of the cell",
            "spawning_offset_days": np.int32(offset)})
           for z in pre_release.SPAWNING_LAYERS},
        # ------------------------------------------------ along the trajectory
        "m4_longest_run": ("u1", ("particle", "ice_threshold"), {
            "units": "day",
            "long_name": "longest run of consecutive days with sea ice above "
                         "the threshold during the calyptopis stages (day 0 "
                         "to calyptopis_end_day, exclusive)"}),
        "m4_first_run_day": ("u1", p, {
            **_DAY,
            "long_name": "day on which a run of consecutive days under sea "
                         "ice first becomes longer than the limit, over the "
                         "whole track",
            "sic_threshold": filters.M4_SIC_THRESHOLD,
            "consecutive_days": np.int32(filters.M4_MAX_CONSECUTIVE_DAYS)}),
        "m4_first_run_dev": ("u2", p, {
            **_DEV, "long_name": "development on m4_first_run_day"}),
        "calyptopis_end_day": ("u1", p, {
            **_DAY, "long_name": "first day at furcilia I or beyond; none if "
                                 "not reached within tracking"}),
        "fiv_day": ("u1", p, {
            **_DAY, "long_name": "first day at furcilia IV or beyond; none "
                                 "if not reached within tracking"}),
        "advance_day": ("u1", p, {
            **_DAY,
            "long_name": "sea-ice advance: first day of the first run of "
                         "consecutive days with sea ice above the threshold, "
                         "from the opening date; none if no advance while "
                         "the particle is tracked",
            "sic_threshold": sea_ice.ADVANCE_SIC_THRESHOLD,
            "consecutive_days": np.int32(sea_ice.ADVANCE_MIN_CONSECUTIVE_DAYS),
            "opening_date": "1 April"}),
        "advance_dev": ("u2", p, {**_DEV, "long_name": "development on advance_day"}),
        "advance_bathymetry": ("f4", p, {
            **_NAN, "units": "m",
            "long_name": "sea-floor depth at the position on advance_day"}),
        "advance_lon": ("f4", p, {**_LON, **_NAN, "long_name": "longitude on advance_day"}),
        "advance_lat": ("f4", p, {**_LAT, **_NAN, "long_name": "latitude on advance_day"}),
        "alive_at_end": ("i1", p, {
            **_FLAG, "flag_meanings": "left_the_domain tracked_to_the_end",
            "long_name": "particle still tracked at the end, whatever its "
                         "outcome",
            "tolerance_days": np.int32(sea_ice.END_OF_TRACKING_TOLERANCE)}),
        "sic_rising_at_end": ("i1", p, {
            **_FLAG, "flag_meanings": "not_rising rising",
            "long_name": "sea ice rising at the end of tracking, whatever "
                         "the outcome: mean of the last 10 days above 0.05 "
                         "and above the mean of the 20 days before"}),
        "release_bathymetry": ("f4", p, {
            **_NAN, "units": "m",
            "long_name": "sea-floor depth at the release position"}),
        "days_below_fit": ("u1", p, {
            "units": "day",
            "long_name": "days with temperature below the fitted range of "
                         "the development function, from release to "
                         "advance_day, or to the end of tracking if none",
            "temperature_threshold": classify.FIT_LOWER_LIMIT}),
        "ice15_daily": ("u1", ("particle", "packed_day"), {
            "long_name": "daily flag: sea ice above the threshold",
            "sic_threshold": sea_ice.ADVANCE_SIC_THRESHOLD,
            "n_days": np.int32(n_days), "comment": packed}),
        "shelf2000_daily": ("u1", ("particle", "packed_day"), {
            "long_name": "daily flag: sea floor shallower than the threshold "
                         "at the position of the day",
            "depth_threshold": classify.SHELF_DEPTH,
            "n_days": np.int32(n_days), "comment": packed}),
    }


def variable_names(n_days=201):
    """Names of the variables of a v3.0.0 file, in order."""
    return list(_schema(n_days))


def write_classification_file(path, particle_id, variables, n_days,
                              attributes):
    """
    Write the outcome and diagnostics of one cohort (dataset v3.0.0).

    Parameters
    ----------
    path : str or pathlib.Path
        Output file. It is written under a temporary name and renamed when
        complete, so an interrupted run leaves no partial file under it.
    particle_id : ndarray, int64
    variables : dict of ndarray
        Every variable of the schema, over the particle dimension.
    n_days : int
        Daily values per trajectory (201).
    attributes : dict
        Global attributes, written after Conventions and title.
    """
    path = Path(path)
    schema = _schema(n_days)
    missing = set(schema) - set(variables)
    extra = set(variables) - set(schema)
    if missing or extra:
        raise ValueError(f"variables do not match the schema: missing "
                         f"{sorted(missing)}, unexpected {sorted(extra)}")
    n = particle_id.size
    sizes = {"particle": n,
             "hatching_depth": len(pre_release.HATCHING_DEPTHS),
             "ice_threshold": len(classify.M4_ICE_THRESHOLDS),
             "packed_day": (n_days + 7) // 8}
    tmp = path.with_name(path.name + ".part")
    with netCDF4.Dataset(tmp, "w", format="NETCDF4") as nc:
        for name, size in sizes.items():
            nc.createDimension(name, size)
        v = nc.createVariable("particle", "i8", ("particle",), zlib=True,
                              shuffle=True, complevel=COMPLEVEL)
        v.long_name = "particle identifier, as in the trajectory file"
        v[:] = np.asarray(particle_id, dtype=np.int64)
        v = nc.createVariable("hatching_depth", "f4", ("hatching_depth",))
        v.units = "m"
        v.positive = "down"
        v.long_name = "depth at which the embryo hatches"
        v[:] = pre_release.HATCHING_DEPTHS
        v = nc.createVariable("ice_threshold", "f4", ("ice_threshold",))
        v.units = "1"
        v.long_name = "sea-ice concentration threshold"
        v[:] = classify.M4_ICE_THRESHOLDS
        for name, (dtype, dims, attrs) in schema.items():
            values = np.asarray(variables[name])
            shape = tuple(sizes[d] for d in dims)
            if values.shape != shape:
                raise ValueError(f"{name}: shape {values.shape}, expected {shape}")
            if values.dtype != np.dtype(dtype):
                raise ValueError(f"{name}: type {values.dtype}, expected {dtype}")
            attrs = dict(attrs)
            v = nc.createVariable(name, dtype, dims, zlib=True, shuffle=True,
                                  complevel=COMPLEVEL, chunksizes=shape,
                                  fill_value=attrs.pop("_FillValue", None))
            v.set_auto_maskandscale(False)
            v.setncatts(attrs)
            v[:] = values
        nc.Conventions = "CF-1.7"
        nc.title = "KRICO per-particle recruitment outcomes and diagnostics"
        nc.setncatts(attributes)
    os.replace(tmp, path)


def read_diagnostics(path):
    """
    Variables of a v3.0.0 file as raw arrays (not masked, not scaled), the
    form classify.outcome expects, with the coordinates it needs:

        d, coordinates = read_diagnostics(path)
        code = classify.outcome(d, classify.Parameters(m1_sic=0.7), **coordinates)
    """
    d = {}
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_maskandscale(False)
        for name, var in nc.variables.items():
            d[name] = np.asarray(var[:])
    coordinates = {
        "hatching_depths": tuple(float(h) for h in d.pop("hatching_depth")),
        "ice_thresholds": tuple(float(t) for t in d.pop("ice_threshold")),
    }
    return d, coordinates
