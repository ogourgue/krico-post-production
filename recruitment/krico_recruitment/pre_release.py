"""
What acts before release: the inputs of M1 and M3, and the temperature at
spawning, read at the release cell from the release-cell files
(scripts/preprocess_glorys12_release_cells.py). No trajectory is needed.

The 24 days from spawning to calyptopis I (sea_ice.SPAWNING_OFFSET_DAYS)
are not simulated. The release position stands as the spawning position,
and the days are counted from the spawning day, day 0:

    day 0            spawning; M1 reads the sea-ice concentration
    days 0 to 4      descent of the embryo, 5 daily means
    days 5 to 23     ascent of the nauplii and metanauplius, 19 daily means
    day 24           release as calyptopis I; first day of the trajectory

M3, as computed here
--------------------
Thorpe et al. (2019) remove embryos and larvae whose mean temperature is
-1 degC or colder during either the descent or the ascent. Their embryo
sinks at a density-dependent rate and hatches after about 5 days at 530 to
1000 m; the ascent temperature is the depth-weighted mean from the hatching
depth to the surface. The forcing held for KRICO has no salinity, so the
sinking rate cannot be computed. Two approximations follow:

  - the hatching depth is bracketed by the two ends of their range, 1000 m
    (baseline) and 530 m, and both are stored;
  - sinking is steady, so the embryo spends equal time at every depth and
    the descent mean is the mean of the water column from the surface to
    the hatching depth, over the descent days. The ascent mean is the same
    column mean over the ascent days.

Column means are thickness-weighted over the model levels, cut at the
hatching depth and at the sea floor of the cell. Where the cell is
shallower than the hatching depth the mean is taken over the levels that
exist; a cell with no value (land) gives NaN, which removes nothing.

Temperature at spawning
-----------------------
On the spawning day, at the release cell: the surface value (first level)
and the means over the top 100 m and 200 m, for the warm limit at spawning
(a sensitivity, not part of the outcome).

The same function serves the classification (one value per particle,
through its cell index), the no-transport control and any later change of
an M3 parameter, which needs only the release-cell files.
"""

from datetime import date, datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np

from krico_recruitment import release_cells, sea_ice

# Days of descent, counted from the spawning day (Thorpe et al. 2019: the
# embryo hatches after 5.0 to 5.6 days).
DESCENT_DAYS = 5

# Hatching depths (m), the two ends of the range of Thorpe et al. (2019).
HATCHING_DEPTHS = (530.0, 1000.0)
BASELINE_HATCHING_DEPTH = 1000.0

# Layers (m, from the surface) of the mean temperature at spawning.
SPAWNING_LAYERS = (100.0, 200.0)

# The release-cell files mark missing values with the netCDF default fill.
MISSING_ABOVE = 1e30


def release_cells_path(directory, release_day):
    """Release-cell file of the spawning year that holds a release day."""
    year = release_cells.spawning_year(release_day)
    return Path(directory) / release_cells.FILE_PATTERN.format(year=year)


def open_static(path):
    """
    Static content of a release-cell file, with the checks the
    classification relies on. Raises ValueError when one fails.

    Returns a dict: path, first_day (datetime.date), n_time, n_cell,
    depth_bounds, bathymetry (cell), lookup ((ny, nx) map to the cell
    index), grid_lon, grid_lat, code_commit, code_version.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"{path} is missing")
    with netCDF4.Dataset(path) as nc:
        for name in ("sea_ice_area_fraction", "temperature", "time",
                     "depth_bounds", "cell_y", "cell_x", "bathymetry",
                     "grid_lon", "grid_lat"):
            if name not in nc.variables:
                raise ValueError(f"{path.name}: variable {name} is missing")
        if nc["sea_ice_area_fraction"].dimensions != ("time", "cell"):
            raise ValueError(f"{path.name}: sea_ice_area_fraction is not (time, cell)")
        if nc["temperature"].dimensions != ("time", "depth", "cell"):
            raise ValueError(f"{path.name}: temperature is not (time, depth, cell)")
        first_day = date.fromisoformat(nc.first_day)
        var = nc["time"]
        stamps = netCDF4.num2date(var[:], var.units, calendar=var.calendar,
                                  only_use_cftime_datetimes=False,
                                  only_use_python_datetimes=True)
        expected = [datetime.combine(first_day + timedelta(days=k),
                                     datetime.min.time()) + timedelta(hours=12)
                    for k in range(len(stamps))]
        if list(np.atleast_1d(stamps)) != expected:
            raise ValueError(f"{path.name}: the time axis is not one day per "
                             f"record at 12:00 from {first_day}")
        bounds = np.asarray(nc["depth_bounds"][:], dtype=np.float64)
        n_lev = nc.dimensions["depth"].size
        if bounds.size != n_lev + 1 or not np.all(np.diff(bounds) > 0):
            raise ValueError(f"{path.name}: depth_bounds does not bound the "
                             f"{n_lev} levels")
        if bounds[-1] < max(HATCHING_DEPTHS):
            raise ValueError(f"{path.name}: levels stop at {bounds[-1]:.0f} m, "
                             "above the deepest hatching depth")
        cell_y = np.asarray(nc["cell_y"][:])
        cell_x = np.asarray(nc["cell_x"][:])
        grid_lon = np.asarray(nc["grid_lon"][:], dtype=np.float64)
        grid_lat = np.asarray(nc["grid_lat"][:], dtype=np.float64)
        static = {
            "path": path,
            "first_day": first_day,
            "n_time": len(stamps),
            "n_cell": cell_y.size,
            "depth_bounds": bounds,
            "bathymetry": np.asarray(nc["bathymetry"][:], dtype=np.float64),
            "grid_lon": grid_lon,
            "grid_lat": grid_lat,
            "lookup": release_cells.cell_lookup(
                cell_y, cell_x, (grid_lat.size, grid_lon.size)),
            "code_commit": getattr(nc, "code_commit", "unknown"),
            "code_version": getattr(nc, "code_version", "unknown"),
        }
    return static


def day_index(static, day):
    """Index of a day along the time axis; ValueError outside the file."""
    t = (day - static["first_day"]).days
    if not 0 <= t < static["n_time"]:
        raise ValueError(f"{day} is not inside {static['path'].name} "
                         f"(from {static['first_day']}, {static['n_time']} days)")
    return t


def _clean(values):
    values = np.asarray(values, dtype=np.float32)
    values[~(np.abs(values) < MISSING_ABOVE)] = np.nan
    return values


def read_ice(static, t):
    """Sea-ice concentration of time index t, (cell,) float32, NaN missing."""
    with netCDF4.Dataset(static["path"]) as nc:
        var = nc["sea_ice_area_fraction"]
        var.set_auto_maskandscale(False)
        return _clean(var[t, :])


def read_temperature(static, t):
    """Temperature of time index t, (depth, cell) float32, NaN missing."""
    with netCDF4.Dataset(static["path"]) as nc:
        var = nc["temperature"]
        var.set_auto_maskandscale(False)
        return _clean(var[t, :, :])


def layer_thickness(depth_bounds, bathymetry, limit):
    """
    Thickness (m) of each level between the surface and limit, cut at the
    sea floor: (depth, cell), zero for a level below either.
    """
    top = depth_bounds[:-1][:, None]
    bottom = np.minimum(depth_bounds[1:][:, None],
                        np.minimum(limit, bathymetry)[None, :])
    return np.clip(bottom - top, 0.0, None)


class _ColumnMean:
    """Thickness-weighted mean of a column, accumulated over days."""

    def __init__(self, thickness):
        self.thickness = thickness
        self.total = np.zeros(thickness.shape[1], dtype=np.float64)
        self.weight = np.zeros(thickness.shape[1], dtype=np.float64)

    def add(self, temperature):
        valid = np.isfinite(temperature)
        w = np.where(valid, self.thickness, 0.0)
        self.total += (np.where(valid, temperature, 0.0) * w).sum(axis=0)
        self.weight += w.sum(axis=0)

    def mean(self):
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(self.weight > 0, self.total / self.weight, np.nan)


def evaluate_cells(static, release_day,
                   offset_days=sea_ice.SPAWNING_OFFSET_DAYS,
                   descent_days=DESCENT_DAYS,
                   hatching_depths=HATCHING_DEPTHS,
                   layers=SPAWNING_LAYERS):
    """
    Inputs of M1 and M3 and the temperature at spawning, for every cell of
    a release-cell file and one release day.

    Returns a dict of float32 arrays over the cell dimension:
      sic_spawning                      (cell,)
      m3_descent_temperature            (cell, hatching depth)
      m3_ascent_temperature             (cell, hatching depth)
      spawning_temperature_surface      (cell,)
      spawning_temperature_mean         (cell, layer)
    """
    if not 0 < descent_days < offset_days:
        raise ValueError("the descent must end before release")
    spawning_day = release_day - timedelta(days=int(offset_days))
    t0 = day_index(static, spawning_day)
    day_index(static, release_day - timedelta(days=1))   # last day read
    bounds, bathymetry = static["depth_bounds"], static["bathymetry"]

    descent = [_ColumnMean(layer_thickness(bounds, bathymetry, h))
               for h in hatching_depths]
    ascent = [_ColumnMean(m.thickness) for m in descent]
    at_spawning = [_ColumnMean(layer_thickness(bounds, bathymetry, z))
                   for z in layers]
    surface = None
    for k in range(int(offset_days)):
        temperature = read_temperature(static, t0 + k)
        if k == 0:
            surface = temperature[0].copy()
            for m in at_spawning:
                m.add(temperature)
        for m in (descent if k < descent_days else ascent):
            m.add(temperature)
    return {
        "sic_spawning": read_ice(static, t0),
        "m3_descent_temperature": np.stack(
            [m.mean() for m in descent], axis=1).astype(np.float32),
        "m3_ascent_temperature": np.stack(
            [m.mean() for m in ascent], axis=1).astype(np.float32),
        "spawning_temperature_surface": surface,
        "spawning_temperature_mean": np.stack(
            [m.mean() for m in at_spawning], axis=1).astype(np.float32),
    }


def locate(static, lon, lat):
    """
    Cell index of the grid point nearest to each release position, the rule
    of M1 since v2.0.0 (release_cells.locate_cells); -1 where that point is
    not in the file.
    """
    return release_cells.locate_cells(lon, lat, static["grid_lon"],
                                      static["grid_lat"], static["lookup"])


def evaluate(static, release_day, cell, **parameters):
    """
    evaluate_cells, taken at the cell of each particle: the same dict, over
    the particle dimension. cell must hold valid indices.
    """
    return {name: values[cell]
            for name, values in evaluate_cells(static, release_day,
                                               **parameters).items()}


# Difference (sea-ice fraction or degC) below which a day-0 value counts as
# reproduced; the tolerance of validate_glorys12_release_cells.py.
DAY0_TOLERANCE = 0.01


def _compare(trajectory, reference):
    both = np.isfinite(trajectory) & np.isfinite(reference)
    n = int(both.sum())
    if n == 0:
        nan = float("nan")
        return {"compared": 0, "exact_share": nan, "within_tolerance_share": nan,
                "mean_abs_difference": nan, "max_abs_difference": nan}
    d = np.abs(trajectory[both].astype(np.float64)
               - reference[both].astype(np.float64))
    return {"compared": n, "exact_share": float((d == 0).mean()),
            "within_tolerance_share": float((d <= DAY0_TOLERANCE).mean()),
            "mean_abs_difference": float(d.mean()),
            "max_abs_difference": float(d.max())}


def day0_check(static, release_day, cell, depth, sic_day0, temperature_day0):
    """
    The trajectory's own sea ice and temperature at day 0 against the
    release-cell file on the release day, at the particle's cell and at the
    level holding its release depth.

    Parcels reads a tracer as the value of the cell holding the particle,
    so the two are equal except for particles within rounding of a cell
    edge or a level bound, and for sea ice from 19 December 2018, whose
    records are stamped 11:52:30 (Parcels then mixes in 0.5% of the next
    day). Expect 99.96 to 99.99% exactly equal, except for sea ice from
    that date, where the share within DAY0_TOLERANCE and the mean
    difference are the meaningful numbers.

    Returns {"sea_ice": {...}, "temperature": {...}}, each with compared,
    exact_share, within_tolerance_share, mean_abs_difference and
    max_abs_difference.
    """
    t = day_index(static, release_day)
    reference_ice = read_ice(static, t)[cell]
    level = release_cells.locate_levels(depth, static["depth_bounds"])
    n_lev = static["depth_bounds"].size - 1
    ok = (level >= 0) & (level < n_lev)
    reference_temperature = np.full(cell.shape, np.nan, dtype=np.float32)
    reference_temperature[ok] = read_temperature(static, t)[level[ok], cell[ok]]
    return {
        "sea_ice": _compare(np.asarray(sic_day0, dtype=np.float32), reference_ice),
        "temperature": _compare(np.asarray(temperature_day0, dtype=np.float32),
                                reference_temperature),
    }
