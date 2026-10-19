"""
Sea-ice concentration from the NSIDC G02202 v6 CDR, sampled where the
recruitment classification samples GLORYS12 ice: at the release position on
the spawning date (M1) and along the stored trajectory (M4, advance,
censoring, hence M5a/M5b/M6).

This is the ice-source swap for the ERA-Interim/ERA5 forcing test: the
trajectories, and therefore the circulation, stay GLORYS12; only the
concentration read at each position and day comes from the observational
record instead of the reanalysis. The two functions return arrays with the
same shape and meaning as their GLORYS12 counterparts, so everything
downstream in the pipeline is unchanged:

    spawning_sic(lon, lat, release_date, nsidc_dir, offset_days)
        -> (n_particles,)         like sea_ice.spawning_sic
    along_trajectory_sic(lon, lat, release_date, nsidc_dir)
        -> (n_particles, n_obs)   like ds["sea_ice_area_fraction"].values

Input files are the monthly stacks written by scripts/preprocess_nsidc.py,
nsidc_ice_YYYY_MM.nc in nsidc_dir: cdr_seaice_conc as ubyte 0-100 with fill
255, one time step per day, on the 25 km EPSG:3412 grid. Sampling is by the
cell containing the position (krico_recruitment.nsidc_grid), with no
interpolation: a 25 km cell average is read as such.

Differences from the GLORYS12 reading that the validation must quantify:

  * Resolution. 25 km cells against ~4 km, so the shelf-break spawning band
    and the coastal ice edge are coarser.
  * Coast cells carry no value in the CDR (fill, like land), whereas
    GLORYS12 has ice right up to its coastline. A position on a coast cell
    therefore reads NaN, which the filters treat as "below threshold".
    With nearest_ocean=True the value of the nearest ocean cell is read
    instead (nearest in grid distance, from the ancillary surface_type
    mask). Off by default; switch on only if the validation shows the
    coast cells matter.
  * NaN positions (particle deleted, or the Parcels end-of-tracking
    artefact) read NaN, as in the trajectory array.

The day of obs k is release_date + k days, the convention the rest of the
pipeline uses (sea_ice.find_april1_day_index); the validation checks it
against the trajectory's own time variable.
"""

import numpy as np
import pandas as pd
import xarray as xr

from krico_recruitment import nsidc_grid
from krico_recruitment.sea_ice import SPAWNING_OFFSET_DAYS

NSIDC_ICE_PATTERN = "nsidc_ice_{year:04d}_{month:02d}.nc"
ANCILLARY = "ancillary/G02202-ancillary-pss25-v06r00.nc"
VAR = "cdr_seaice_conc"
FILL = 255
SURFACE_OCEAN = 50

SOURCE_ID = "NSIDC G02202 v6 (doi:10.7265/b18j-z797)"

# A second product in the same file format: GLORYS12 ice binned onto the
# NSIDC grid by scripts/preprocess_glorys12_25km.py, read with the same code
# so that the only difference between the two runs is the ice values.
SOURCES = {
    "nsidc": {
        "pattern": NSIDC_ICE_PATTERN,
        "id": SOURCE_ID,
    },
    "glorys12_25km": {
        "pattern": "glorys12_25km_ice_{year:04d}_{month:02d}.nc",
        "id": "GLORYS12v1 ice binned onto the NSIDC G02202 25 km grid",
    },
}

# Raw monthly arrays (uint8, (ndays, NROW, NCOL), ~3 MB each) cached by
# (dir, year, month). A cohort needs at most 9 consecutive months (spawning
# date to release + 200 days); consecutive cohorts overlap almost entirely.
_MONTH_CACHE = {}
_MONTH_CACHE_MAX = 16

_NEAREST_OCEAN = {}

# Particles per chunk when sampling along trajectories. A cohort has ~545k
# particles x 201 obs; indexing it in one go allocates several full-size
# temporaries (several GB). Chunking keeps the peak to the output array plus
# a few chunk-size temporaries.
CHUNK = 20_000


def _load_month(nsidc_dir, year, month, pattern=NSIDC_ICE_PATTERN):
    """Raw uint8 concentration of one month, (ndays, NROW, NCOL), and its days."""
    key = (str(nsidc_dir), pattern, year, month)
    if key not in _MONTH_CACHE:
        path = f"{nsidc_dir}/" + pattern.format(year=year, month=month)
        with xr.open_dataset(path, mask_and_scale=False) as ds:
            raw = np.asarray(ds[VAR].values, dtype=np.uint8)
            days = pd.DatetimeIndex(ds["time"].values).normalize()
        if raw.shape[1:] != (nsidc_grid.NROW, nsidc_grid.NCOL):
            raise ValueError(f"{path}: unexpected grid shape {raw.shape[1:]}")
        if len(_MONTH_CACHE) >= _MONTH_CACHE_MAX:
            _MONTH_CACHE.pop(next(iter(_MONTH_CACHE)))
        _MONTH_CACHE[key] = (raw, days)
    return _MONTH_CACHE[key]


def _load_range(nsidc_dir, first_day, last_day, pattern=NSIDC_ICE_PATTERN):
    """
    Raw concentration for every day from first_day to last_day inclusive,
    as one (n_days, NROW, NCOL) uint8 array. Raises if a day is missing.
    """
    first_day = pd.Timestamp(first_day).normalize()
    last_day = pd.Timestamp(last_day).normalize()
    wanted = pd.date_range(first_day, last_day, freq="D")

    blocks, got = [], []
    for period in pd.period_range(first_day, last_day, freq="M"):
        raw, days = _load_month(nsidc_dir, period.year, period.month, pattern)
        keep = (days >= first_day) & (days <= last_day)
        blocks.append(raw[keep])
        got.append(days[keep])
    got = got[0].append(got[1:]) if len(got) > 1 else got[0]
    if not got.equals(wanted):
        missing = wanted.difference(got)
        raise ValueError(
            f"NSIDC monthly files in {nsidc_dir} do not cover "
            f"{first_day.date()} to {last_day.date()}; missing "
            f"{len(missing)} days, first {missing[0].date()}"
        )
    return np.concatenate(blocks, axis=0)


def _nearest_ocean_map(nsidc_dir):
    """
    (2, NROW, NCOL) int array: row and column of the nearest ocean cell to
    every cell (itself where the cell is ocean), from the ancillary
    surface_type mask.
    """
    key = str(nsidc_dir)
    if key not in _NEAREST_OCEAN:
        from scipy.ndimage import distance_transform_edt
        with xr.open_dataset(f"{nsidc_dir}/{ANCILLARY}", mask_and_scale=False) as ds:
            surface = np.asarray(ds["surface_type"].values)
        if surface.shape != (nsidc_grid.NROW, nsidc_grid.NCOL):
            raise ValueError("ancillary surface_type has an unexpected shape")
        idx = distance_transform_edt(surface != SURFACE_OCEAN,
                                     return_distances=False, return_indices=True)
        _NEAREST_OCEAN[key] = np.asarray(idx)
    return _NEAREST_OCEAN[key]


def cell_index(lon, lat):
    """
    Row and column (int32) of the grid cell containing each position; -1
    where the position is NaN or off the grid. Any shape.
    """
    lon = np.asarray(lon, dtype="float64")
    lat = np.asarray(lat, dtype="float64")
    finite = np.isfinite(lon) & np.isfinite(lat)
    row, col = nsidc_grid.lonlat_to_index(np.where(finite, lon, 0.0),
                                          np.where(finite, lat, -70.0))
    row = np.where(finite, row, -1).astype(np.int32)
    col = np.where(finite, col, -1).astype(np.int32)
    return row, col


def _read(raw, t, row, col, nsidc_dir, nearest_ocean):
    """
    Read raw[t, row, col] as a 0-1 float32 fraction; NaN where the index is
    -1 or the cell is fill. With nearest_ocean, fill cells are redirected to
    the nearest ocean cell before reading.
    """
    valid = (row >= 0) & (col >= 0)
    r = np.where(valid, row, 0)
    c = np.where(valid, col, 0)
    if nearest_ocean:
        near = _nearest_ocean_map(nsidc_dir)
        r, c = near[0][r, c], near[1][r, c]
    packed = raw[t, r, c]
    value = packed.astype("float32")
    value /= 100.0
    value[(packed == FILL) | ~valid] = np.nan
    return value


def spawning_sic(lon, lat, release_date, nsidc_dir,
                 offset_days=SPAWNING_OFFSET_DAYS, nearest_ocean=False,
                 pattern=NSIDC_ICE_PATTERN):
    """
    Concentration at the release position on the spawning date, from the
    NSIDC record (or another product in the same format, via pattern).
    Same contract as sea_ice.spawning_sic.
    """
    release_date = pd.Timestamp(release_date).normalize()
    spawning_date = release_date - pd.Timedelta(days=int(offset_days))
    raw = _load_range(nsidc_dir, spawning_date, spawning_date, pattern)
    row, col = cell_index(lon, lat)
    t = np.zeros(row.shape, dtype=np.int64)
    return _read(raw, t, row, col, nsidc_dir, nearest_ocean)


def along_trajectory_sic(lon, lat, release_date, nsidc_dir, nearest_ocean=False,
                         pattern=NSIDC_ICE_PATTERN):
    """
    Concentration along the trajectory from the NSIDC record (or another
    product in the same format, via pattern): obs k of each particle is read
    on release_date + k days at (lon[:, k], lat[:, k]).

    Returns (n_particles, n_obs) float64, fraction 0-1, NaN where the
    position is NaN, off the grid, or on a fill (land/coast/lake) cell.
    """
    lon = np.asarray(lon, dtype="float64")
    lat = np.asarray(lat, dtype="float64")
    if lon.ndim != 2 or lon.shape != lat.shape:
        raise ValueError("lon and lat must be (n_particles, n_obs) arrays")
    n_particles, n_obs = lon.shape

    release_date = pd.Timestamp(release_date).normalize()
    raw = _load_range(nsidc_dir, release_date,
                      release_date + pd.Timedelta(days=n_obs - 1), pattern)

    t = np.broadcast_to(np.arange(n_obs, dtype=np.int32), (CHUNK, n_obs))
    out = np.empty((n_particles, n_obs), dtype="float32")
    for i0 in range(0, n_particles, CHUNK):
        i1 = min(i0 + CHUNK, n_particles)
        row, col = cell_index(lon[i0:i1], lat[i0:i1])
        out[i0:i1] = _read(raw, t[: i1 - i0], row, col, nsidc_dir, nearest_ocean)
    return out


def cell_surface_type(lon, lat, nsidc_dir):
    """
    surface_type of the cell under each position (50 ocean, 75 lake,
    200 coast, 250 land), 0 where the position is NaN or off the grid, and
    adj123 (0 non-ocean, 1-3 ocean cells from land, 255 open ocean).
    For validation and reporting.
    """
    with xr.open_dataset(f"{nsidc_dir}/{ANCILLARY}", mask_and_scale=False) as ds:
        surface = np.asarray(ds["surface_type"].values)
        adj = np.asarray(ds["adj123"].values)
    lon = np.asarray(lon)
    lat = np.asarray(lat)
    out_s = np.zeros(lon.shape, dtype=np.uint8)
    out_a = np.zeros(lon.shape, dtype=np.uint8)
    flat_lon, flat_lat = lon.reshape(-1), lat.reshape(-1)
    fs, fa = out_s.reshape(-1), out_a.reshape(-1)
    step = CHUNK * 64
    for i0 in range(0, flat_lon.size, step):
        i1 = min(i0 + step, flat_lon.size)
        row, col = cell_index(flat_lon[i0:i1], flat_lat[i0:i1])
        valid = (row >= 0) & (col >= 0)
        r = np.where(valid, row, 0)
        c = np.where(valid, col, 0)
        fs[i0:i1] = np.where(valid, surface[r, c], 0)
        fa[i0:i1] = np.where(valid, adj[r, c], 0)
    return out_s, out_a
