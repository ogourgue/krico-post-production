"""
Chlorophyll along the trajectories: pass 2 of dataset v3.0.0.

Chlorophyll is a diagnostic, not a rule: it changes no outcome. For every
particle, including those removed by M1 or M3, the chlorophyll at its
position is read day by day from two products and averaged over two
periods of the larva's life.

Products
--------
    satellite   Copernicus-GlobColour daily climatology of surface
                chlorophyll (scripts/download_globcolour.py): 366 days of
                the year, 1/24 degree. Read by day of year. Each value
                comes with the number of observations behind it
                (CHL_count); it has none under ice or in the dark.
    model       PISCES hindcast (scripts/download_pisces.py): daily means,
                1/4 degree, 75 levels. Read on the actual date.

Windows
-------
    calyptopis  from release (day 0) to the day before the first day at
                furcilia I: the days M4 looks at
    furcilia    from the first day at furcilia I to the day before the
                sea-ice advance, or to the end of tracking if there is no
                advance; empty if the advance comes during the calyptopis
                stages or if furcilia I is never reached

The two bounds, calyptopis_end_day and advance_day, are those of pass 1:
they come from classify.trajectory_diagnostics, the function that writes
them in the pass 1 files, with the sea ice of the trajectory file
(GLORYS12). Only days on which the particle has a position count.

Indices, each the arithmetic mean over the days of the window that hold a
value, stored with the number of those days
---------------------------------------------------------------------
    satellite mean, satellite median
                the climatological mean and median, over the days on which
                both exist and CHL_count reaches a minimum; given for each
                minimum of MIN_COUNTS (1 is no minimum)
    model surface
                first level (0.5 m)
    model at depth
                level holding the particle's depth
    model upper layer
                thickness-weighted mean from the surface to LAYER_DEPTH,
                cut at the sea floor of the model
    model surface, matched
                first level over the days on which the satellite index
                holds a value too, for each minimum count: the comparison
                of the two products on the same days

Sampling
--------
The value of the cell that holds the position: the nearest grid point of a
regular grid, with no interpolation and no filling from neighbours, as
Parcels reads a tracer and as the other ice sources are read. A position
outside the grid, on a land cell, or below the sea floor of the model,
reads no value. Model level bounds are taken midway between the level
depths, since the files carry no bounds; the level holding a depth is
therefore the nearest one.

Storage
-------
Means are stored as float32 rounded to SIGNIFICANT_BITS significant binary
digits (round_bits): a relative error of at most 2**-SIGNIFICANT_BITS,
under 1%, far below the uncertainty of either product, for about a third
of the size.
"""

import re
from datetime import timedelta
from pathlib import Path

import netCDF4
import numpy as np

# Minimum number of observations behind a satellite value (CHL_count) for
# the day to count. 1 is no minimum. Of 220 possible (20 years, 11 days).
MIN_COUNTS = (1, 10)

WINDOWS = ("calyptopis", "furcilia")

# Lower limit (m) of the upper layer of the model index.
LAYER_DEPTH = 200.0

# Significant binary digits kept in the stored means.
SIGNIFICANT_BITS = 7

SATELLITE_FILE = "globcolour_chl_climatology.nc"
SATELLITE_VARIABLES = ("CHL_mean", "CHL_median", "CHL_count")
MODEL_PATTERN = "pisces_chl_{year:04d}_{month:02d}.nc"
MODEL_VARIABLE = "chl"
MODEL_LEVELS = 75

_SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
PROVENANCE = {"satellite": _SCRIPTS / "globcolour_chl_provenance.txt",
              "model": _SCRIPTS / "pisces_chl_provenance.txt"}

# Above this, a float is a fill value (the model files use 9.96921e36).
_HUGE = 1e30


def round_bits(values, bits=SIGNIFICANT_BITS):
    """
    float32 rounded to `bits` significant binary digits, half away from
    zero; NaN stays NaN. The relative error is at most 2**-bits.
    """
    x = np.asarray(values, dtype=np.float32)
    drop = 24 - bits                       # of the 24 significant bits of a float32
    u = x.view(np.uint32).astype(np.uint64)
    u = (u + (1 << (drop - 1))) & ~np.uint64((1 << drop) - 1)
    out = u.astype(np.uint32).view(np.float32)
    return np.where(np.isfinite(x), out, x).astype(np.float32)


def source(product):
    """
    One line naming a product as it was downloaded (product identifier,
    DOI, dataset and version), from the provenance file its download script
    wrote.
    """
    text = PROVENANCE[product].read_text()
    found = {}
    for key in ("doi", "dataset", "version"):
        match = re.search(rf"^# {key}\s*:\s*(.+)$", text, flags=re.MULTILINE)
        if not match:
            raise ValueError(f"{PROVENANCE[product].name}: no line '{key}'")
        found[key] = match.group(1).strip()
    code = re.search(r"^#\s+([A-Z0-9_]{12,})", text, flags=re.MULTILINE)
    if not code:
        raise ValueError(f"{PROVENANCE[product].name}: no product identifier")
    return (f"{code.group(1)}, doi:{found['doi']}, dataset "
            f"{found['dataset']}, version {found['version']}")


class _Axis:
    """A regular, ascending coordinate axis and the cell holding a value."""

    def __init__(self, values, name, where):
        values = np.asarray(values, dtype=np.float64)
        if values.ndim != 1 or values.size < 2:
            raise ValueError(f"{where}: {name} is not an axis")
        self.n = values.size
        self.first = float(values[0])
        self.step = float(values[-1] - values[0]) / (self.n - 1)
        if self.step <= 0:
            raise ValueError(f"{where}: {name} is not ascending")
        even = self.first + self.step * np.arange(self.n)
        if np.abs(values - even).max() > 0.01 * self.step:
            raise ValueError(f"{where}: {name} is not evenly spaced")

    def index(self, x):
        """Index of the cell holding each x, and whether x is on the grid."""
        with np.errstate(invalid="ignore"):
            i = np.floor((np.asarray(x, dtype=np.float64) - self.first)
                         / self.step + 0.5)
            inside = (i >= 0) & (i < self.n)
        return np.where(inside, i, 0).astype(np.int64), inside

    def same(self, other):
        return (self.n == other.n and abs(self.first - other.first) < 1e-6
                and abs(self.step - other.step) < 1e-9)


def _dates(nc):
    """Calendar days of the time axis of an open file."""
    time = nc["time"]
    stamps = netCDF4.num2date(np.asarray(time[:], dtype=np.float64), time.units,
                              getattr(time, "calendar", "standard"))
    return [(int(s.year), int(s.month), int(s.day)) for s in np.atleast_1d(stamps)]


class Satellite:
    """The daily climatology of surface chlorophyll, read by day of year."""

    def __init__(self, directory):
        self.path = Path(directory) / SATELLITE_FILE
        name = self.path.name
        if not self.path.is_file():
            raise FileNotFoundError(f"{self.path} is missing")
        self.nc = netCDF4.Dataset(self.path)
        self.nc.set_auto_maskandscale(False)
        for var in SATELLITE_VARIABLES:
            if var not in self.nc.variables:
                raise ValueError(f"{name}: variable {var} is missing")
            if self.nc[var].dimensions != ("time", "latitude", "longitude"):
                raise ValueError(f"{name}: {var} has dimensions "
                                 f"{self.nc[var].dimensions}, expected "
                                 "('time', 'latitude', 'longitude')")
            for attribute in ("scale_factor", "add_offset"):
                if attribute in self.nc[var].ncattrs():
                    raise ValueError(f"{name}: {var} is packed ({attribute}); "
                                     "the reader expects plain values")
        self.lat = _Axis(self.nc["latitude"][:], "latitude", name)
        self.lon = _Axis(self.nc["longitude"][:], "longitude", name)
        days = [(m, d) for _, m, d in _dates(self.nc)]
        if len(days) != 366 or len(set(days)) != 366:
            raise ValueError(f"{name}: {len(set(days))} distinct days of the "
                             f"year in {len(days)} records, expected 366")
        self.day_index = {day: i for i, day in enumerate(days)}
        self.fill = {var: getattr(self.nc[var], "_FillValue", None)
                     for var in SATELLITE_VARIABLES}
        self._loaded = (None, None)

    def close(self):
        self.nc.close()

    def fields(self, day):
        """
        (mean, median, count) of a calendar day (datetime.date): float32
        with NaN where there is no value, and the count as int32 with 0.
        """
        i = self.day_index[(day.month, day.day)]
        if self._loaded[0] != i:
            out = []
            for var in SATELLITE_VARIABLES[:2]:
                a = np.asarray(self.nc[var][i], dtype=np.float32)
                with np.errstate(invalid="ignore"):
                    missing = ~np.isfinite(a) | (a < 0) | (np.abs(a) > _HUGE)
                if self.fill[var] is not None:
                    missing |= a == self.fill[var]
                a[missing] = np.nan
                out.append(a)
            count = np.asarray(self.nc["CHL_count"][i]).astype(np.int32)
            count[count < 0] = 0
            out.append(count)
            self._loaded = (i, tuple(out))
        return self._loaded[1]

    def sample(self, day, lon, lat):
        """
        Mean, median and count at positions on a calendar day. A day counts
        when the mean and the median both hold a value; elsewhere the two
        are NaN and the count 0.
        """
        mean, median, count = self.fields(day)
        row, in_lat = self.lat.index(lat)
        col, in_lon = self.lon.index(lon)
        a, b, c = mean[row, col], median[row, col], count[row, col]
        ok = in_lat & in_lon & np.isfinite(a) & np.isfinite(b) & (c >= 1)
        return (np.where(ok, a, np.nan), np.where(ok, b, np.nan),
                np.where(ok, c, 0))


class Model:
    """The daily 3D chlorophyll of the hindcast, read on the actual date."""

    def __init__(self, directory):
        self.directory = Path(directory)
        self._month = (None, None, None)      # (year, month), dataset, days
        self._loaded = (None, 0, None)        # day, levels read, fields
        self._known = {}                      # (year, month) -> days it holds
        self.lat = self.lon = self.depth = None

    def close(self):
        if self._month[1] is not None:
            self._month[1].close()
        self._month = (None, None, None)

    def _open(self, year, month):
        if self._month[0] != (year, month):
            self.close()
            path = self.directory / MODEL_PATTERN.format(year=year, month=month)
            name = path.name
            if not path.is_file():
                raise FileNotFoundError(f"{path} is missing")
            nc = netCDF4.Dataset(path)
            nc.set_auto_maskandscale(False)
            try:
                if MODEL_VARIABLE not in nc.variables:
                    raise ValueError(f"{name}: variable {MODEL_VARIABLE} is missing")
                var = nc[MODEL_VARIABLE]
                if var.dimensions != ("time", "depth", "latitude", "longitude"):
                    raise ValueError(f"{name}: {MODEL_VARIABLE} has dimensions "
                                     f"{var.dimensions}, expected ('time', "
                                     "'depth', 'latitude', 'longitude')")
                for attribute in ("scale_factor", "add_offset"):
                    if attribute in var.ncattrs():
                        raise ValueError(f"{name}: {MODEL_VARIABLE} is packed "
                                         f"({attribute}); the reader expects "
                                         "plain values")
                lat = _Axis(nc["latitude"][:], "latitude", name)
                lon = _Axis(nc["longitude"][:], "longitude", name)
                depth = np.asarray(nc["depth"][:], dtype=np.float64)
                if depth.size != MODEL_LEVELS or (np.diff(depth) <= 0).any():
                    raise ValueError(f"{name}: {depth.size} depth levels, "
                                     f"expected {MODEL_LEVELS} increasing ones")
                if self.depth is None:
                    self.lat, self.lon, self.depth = lat, lon, depth
                    middle = 0.5 * (depth[:-1] + depth[1:])
                    self.bounds = np.concatenate(
                        ([0.0], middle, [2 * depth[-1] - middle[-1]]))
                    top = np.minimum(self.bounds[:-1], LAYER_DEPTH)
                    bottom = np.minimum(self.bounds[1:], LAYER_DEPTH)
                    self.layer_thickness = bottom - top
                    self.layer_levels = int(np.count_nonzero(self.layer_thickness))
                elif not (lat.same(self.lat) and lon.same(self.lon)
                          and np.array_equal(depth, self.depth)):
                    raise ValueError(f"{name}: not on the grid of the other "
                                     "monthly files")
                days = _dates(nc)
                if any((y, m) != (year, month) for y, m, _ in days) \
                        or len(set(days)) != len(days):
                    raise ValueError(f"{name}: records outside {year}-"
                                     f"{month:02d}, or repeated days")
            except Exception:
                nc.close()
                raise
            self._month = ((year, month), nc, {d: k for k, (_, _, d) in enumerate(days)})
            self._known[(year, month)] = frozenset(self._month[2])
        return self._month[1], self._month[2]

    def check_day(self, day):
        """
        Raise unless the monthly file of `day` exists, is on the grid of
        the others and holds that day. A month is opened and checked once.
        """
        key = (day.year, day.month)
        if key not in self._known:
            self._open(*key)
        if day.day not in self._known[key]:
            raise ValueError(
                f"{MODEL_PATTERN.format(year=day.year, month=day.month)}: "
                f"no record for {day}")

    def fields(self, day, levels):
        """The first `levels` levels of a calendar day, float32, NaN where none."""
        if self._loaded[0] != day or self._loaded[1] < levels:
            self.check_day(day)
            nc, days = self._open(day.year, day.month)
            a = np.asarray(nc[MODEL_VARIABLE][days[day.day], :levels],
                           dtype=np.float32)
            with np.errstate(invalid="ignore"):
                a[~np.isfinite(a) | (np.abs(a) > _HUGE) | (a < 0)] = np.nan
            self._loaded = (day, levels, a)
        return self._loaded[2]

    def sample(self, day, lon, lat, depth):
        """
        Chlorophyll at positions on a calendar day: at the first level, at
        the level holding each depth, and the mean from the surface to
        LAYER_DEPTH. NaN where the model has no value.
        """
        self.check_day(day)                   # also sets the grid, once
        row, in_lat = self.lat.index(lat)
        col, in_lon = self.lon.index(lon)
        inside = in_lat & in_lon
        depth = np.asarray(depth, dtype=np.float64)
        level = np.searchsorted(self.bounds, depth, side="right") - 1
        with np.errstate(invalid="ignore"):
            in_depth = inside & (depth >= 0) & (level < MODEL_LEVELS)
        level = np.where(in_depth, level, 0)
        needed = max(self.layer_levels, int(level.max()) + 1 if level.size else 0)
        a = self.fields(day, needed)

        surface = np.where(inside, a[0][row, col], np.nan)
        at_depth = np.where(in_depth, a[level, row, col], np.nan)
        column = a[:self.layer_levels, row, col].astype(np.float64)   # (level, particle)
        weight = np.where(np.isfinite(column),
                          self.layer_thickness[:self.layer_levels, None], 0.0)
        total = weight.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            layer = np.where(inside & (total > 0),
                             (np.nan_to_num(column) * weight).sum(axis=0) / total,
                             np.nan)
        return surface, at_depth, layer


def window_bounds(calyptopis_end_day, advance_day, n_obs, none):
    """
    First day and day after the last of each window, per particle, from
    the two bounds in their stored form (`none` where there is none).
    Returns (start, stop), each of shape (n_particles, 2).
    """
    end = np.asarray(calyptopis_end_day).astype(np.int64)
    end = np.where(end == none, n_obs, end)
    advance = np.asarray(advance_day).astype(np.int64)
    advance = np.where(advance == none, n_obs, advance)
    start = np.stack([np.zeros_like(end), end], axis=1)
    stop = np.stack([end, advance], axis=1)
    return start, stop


def check_inputs(satellite, model, release_day, n_obs):
    """
    Raise ValueError or FileNotFoundError unless both products hold every
    day from release_day to the last tracked day. The satellite object has
    been checked when it was opened.
    """
    for k in range(n_obs):
        day = release_day + timedelta(days=k)
        if (day.month, day.day) not in satellite.day_index:
            raise ValueError(f"{satellite.path.name}: no record for "
                             f"{day:%d %B}")
        model.check_day(day)


def along_trajectory(lon, lat, depth, start, stop, release_day, satellite,
                     model, log=None):
    """
    The chlorophyll indices of every particle, in stored form.

    Parameters
    ----------
    lon, lat, depth : ndarray (n_particles, n_obs), daily, NaN after deletion
    start, stop : ndarray (n_particles, 2), from window_bounds
    release_day : datetime.date; day k of a trajectory is release_day + k
    satellite, model : Satellite, Model

    Returns
    -------
    dict of arrays over (particle, window) or (particle, window, min_count)
    """
    n, n_obs = lon.shape
    n_w, n_c = len(WINDOWS), len(MIN_COUNTS)
    total = {name: np.zeros((n, n_w), dtype=np.float64)
             for name in ("surface", "depth", "layer")}
    total.update({name: np.zeros((n, n_w, n_c), dtype=np.float64)
                  for name in ("mean", "median", "matched")})
    days = {name: np.zeros((n, n_w), dtype=np.int32)
            for name in ("window", "surface", "depth")}
    days.update({name: np.zeros((n, n_w, n_c), dtype=np.int32)
                 for name in ("satellite", "matched")})

    for k in range(n_obs):
        x, y, z = lon[:, k], lat[:, k], depth[:, k]
        here = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        member = here[:, None] & (start <= k) & (k < stop)
        rows = np.flatnonzero(member.any(axis=1))
        if rows.size == 0:
            continue
        # The windows do not overlap: one window per particle and day.
        w = np.argmax(member[rows], axis=1)
        day = release_day + timedelta(days=k)
        days["window"][rows, w] += 1

        surface, at_depth, layer = model.sample(day, x[rows], y[rows], z[rows])
        has_model = np.isfinite(surface) & np.isfinite(layer)
        r, c = rows[has_model], w[has_model]
        total["surface"][r, c] += surface[has_model]
        total["layer"][r, c] += layer[has_model]
        days["surface"][r, c] += 1
        ok = has_model & np.isfinite(at_depth)
        total["depth"][rows[ok], w[ok]] += at_depth[ok]
        days["depth"][rows[ok], w[ok]] += 1

        mean, median, count = satellite.sample(day, x[rows], y[rows])
        for j, minimum in enumerate(MIN_COUNTS):
            ok = count >= minimum
            r, c = rows[ok], w[ok]
            total["mean"][r, c, j] += mean[ok]
            total["median"][r, c, j] += median[ok]
            days["satellite"][r, c, j] += 1
            ok &= has_model
            r, c = rows[ok], w[ok]
            total["matched"][r, c, j] += surface[ok]
            days["matched"][r, c, j] += 1
        if log and (k + 1) % 50 == 0:
            log(f"  day {k + 1} of {n_obs}")

    def mean_of(name, count):
        with np.errstate(invalid="ignore", divide="ignore"):
            return round_bits(np.where(count > 0, total[name] / count, np.nan))

    return {
        "window_days": days["window"].astype(np.uint8),
        "chl_satellite_mean": mean_of("mean", days["satellite"]),
        "chl_satellite_median": mean_of("median", days["satellite"]),
        "chl_satellite_days": days["satellite"].astype(np.uint8),
        "chl_model_surface": mean_of("surface", days["surface"]),
        "chl_model_mean_200m": mean_of("layer", days["surface"]),
        "chl_model_days": days["surface"].astype(np.uint8),
        "chl_model_at_depth": mean_of("depth", days["depth"]),
        "chl_model_at_depth_days": days["depth"].astype(np.uint8),
        "chl_model_surface_matched": mean_of("matched", days["matched"]),
        "chl_matched_days": days["matched"].astype(np.uint8),
    }
