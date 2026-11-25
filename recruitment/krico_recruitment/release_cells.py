"""
Release cells of the KRICO simulations on the GLORYS12 grid, and the daily
series extracted at those cells (scripts/preprocess_glorys12_release_cells.py).

The rules that act before release (M1 at spawning, M3 during the descent
and the ascent) depend only on the date and on the grid cell of the release
position, not on the trajectory. Their inputs are therefore extracted once
per spawning year at the cells where particles are released, and every
particle is linked to its cell.

Cells
-----
The release cells are the 35,303 GLORYS12 cells in which run.py of
krico-templates generates particles: sea floor between 1000 and 2000 m,
centre inside CCAMLR subareas 48.1 to 48.6 and 88.3. They are listed in
scripts/glorys12_release_cells.csv (0-based j, i in the 968 x 1858 KRICO
subset of the grid).

A particle is linked to the grid point nearest to its release position,
the rule M1 has used since v2.0.0 (sea_ice._sample_nearest). That point is
not always the cell the particle was generated in:

  - the latitude spacing shrinks southward, and run.py jitters a particle by
    up to half the northward spacing, so about 0.06% of particles lie
    nearest to the cell south of their own;
  - the float32 longitudes of the grid are not evenly spaced to the last
    bit, so a few particles per cohort lie nearest to the cell west of
    their own.

The extraction therefore holds every release cell and its eight neighbours
(NEIGHBOUR_OFFSETS). A jitter never exceeds half a cell, so no particle can
lie nearest to a cell outside that set.

Checked on the 2,100,265,336 release positions of the v2.0.0 dataset
(validate_glorys12_release_cells.py, 25 November 2026): none has its nearest
grid point outside the set; 242,288 have it in a cell south of a release
cell, 1,955 in a cell west of one, 4 in another neighbour and none
south-west. The southern neighbours alone would not have been enough.

Period
------
One file per spawning year Y, from 22 October Y-1 to 30 September Y: the
first spawning date (the release of 15 November minus the 24 days of
sea_ice.SPAWNING_OFFSET_DAYS) to the last tracked day (the release of
14 March plus 200 days).
"""

from datetime import date, timedelta
from pathlib import Path

import numpy as np
import xarray as xr

from krico_recruitment import sea_ice

FILE_PATTERN = "glorys12_release_cells_{year:04d}.nc"
TABLE = Path(__file__).resolve().parents[1] / "scripts" / "glorys12_release_cells.csv"

FIRST_YEAR = 1994
LAST_YEAR = 2025

# Release window of a spawning year and tracking duration (krico-templates).
FIRST_RELEASE = (11, 15)   # of the previous calendar year
LAST_RELEASE = (3, 14)
TRACKING_DAYS = 200

# Bathymetry band of the release cells (run.py), for checking the table.
DEPTH_MIN = 1000.0
DEPTH_MAX = 2000.0

# Files of the same grid may store their coordinates with different
# rounding; beyond this difference (degrees, about 1 m) they are another grid.
COORDINATE_TOLERANCE = 1e-5

# A release cell and its eight neighbours, as (dj, di).
NEIGHBOUR_OFFSETS = tuple((dj, di) for dj in (-1, 0, 1) for di in (-1, 0, 1))


def period(year):
    """First and last day (datetime.date) held by the file of spawning year."""
    first = (date(year - 1, *FIRST_RELEASE)
             - timedelta(days=sea_ice.SPAWNING_OFFSET_DAYS))
    last = date(year, *LAST_RELEASE) + timedelta(days=TRACKING_DAYS)
    return first, last


def days(year):
    """Every day of the period of a spawning year, as datetime.date."""
    first, last = period(year)
    return [first + timedelta(days=k) for k in range((last - first).days + 1)]


def spawning_year(day):
    """Spawning year whose period holds a given day (October onward: next)."""
    return day.year + 1 if day.month >= 10 else day.year


def read_table(path=TABLE):
    """Release-cell table as a dict of arrays: j, i, lon, lat, depth_m, n_release."""
    raw = np.genfromtxt(path, delimiter=",", names=True)
    table = {
        "j": raw["j"].astype(np.int64),
        "i": raw["i"].astype(np.int64),
        "lon": raw["lon"].astype(np.float64),
        "lat": raw["lat"].astype(np.float64),
        "depth_m": raw["depth_m"].astype(np.float64),
        "n_release": raw["n_release"].astype(np.int64),
    }
    pairs = table["j"] * 100000 + table["i"]
    if np.unique(pairs).size != pairs.size:
        raise ValueError(f"duplicated cells in {path}")
    return table


def _first(glorys_dir, variable):
    files = sorted(Path(glorys_dir).glob(f"glorys12_{variable}_????_??.nc"))
    if not files:
        raise FileNotFoundError(f"no glorys12_{variable}_YYYY_MM.nc in {glorys_dir}")
    return files[0]


def coordinate_difference(lon, lat, grid):
    """Largest difference (degrees) between coordinates and those of the grid."""
    if lon.shape != grid["nav_lon"].shape:
        raise ValueError(f"grid {lon.shape} differs from the sea-ice grid "
                         f"{grid['nav_lon'].shape}")
    return max(float(np.abs(np.asarray(lon, dtype=np.float64) - grid["nav_lon"]).max()),
               float(np.abs(np.asarray(lat, dtype=np.float64) - grid["nav_lat"]).max()))


def read_grid(glorys_dir):
    """
    Static description of the GLORYS12 subset in glorys_dir.

    Returns a dict with:
      nav_lon, nav_lat   (ny, nx), as stored in the sea-ice files, the
                         coordinates M1 has read since v2.0.0
      grid_lon, grid_lat their ascending 1-D axes (float64), exactly as M1
                         derives them (sea_ice._grid_axes)
      bathymetry         (ny, nx) float32, hdept
      deptht             depth of the temperature levels (m)
      depthw             depth of the top of each level (m); the bottom of
                         level k is depthw[k + 1]
      coordinate_difference   largest difference (degrees) between those
                         coordinates and the ones of the bathymetry and
                         temperature files; above COORDINATE_TOLERANCE the
                         files do not describe the same grid (ValueError)
    """
    glorys_dir = Path(glorys_dir)
    ice_file = _first(glorys_dir, "ice")
    with xr.open_dataset(ice_file, decode_times=False) as ds:
        grid_lon, grid_lat = sea_ice._grid_axes(ds)
        nav_lon = np.asarray(ds["nav_lon"].values)
        nav_lat = np.asarray(ds["nav_lat"].values)
    grid = {"nav_lon": nav_lon, "nav_lat": nav_lat,
            "grid_lon": np.asarray(grid_lon, dtype=np.float64),
            "grid_lat": np.asarray(grid_lat, dtype=np.float64)}
    with xr.open_dataset(glorys_dir / "glorys12_bathymetry.nc",
                         decode_times=False) as ds:
        hdept = np.asarray(ds["hdept"].values, dtype=np.float32)
        hdept = hdept.reshape(hdept.shape[-2:])
        diff_bathy = coordinate_difference(ds["nav_lon"].values,
                                           ds["nav_lat"].values, grid)
    temp_file = _first(glorys_dir, "temp")
    with xr.open_dataset(temp_file, decode_times=False) as ds:
        deptht = np.asarray(ds["deptht"].values, dtype=np.float64)
        diff_temp = coordinate_difference(ds["nav_lon"].values,
                                          ds["nav_lat"].values, grid)
    with xr.open_dataset(_first(glorys_dir, "w"), decode_times=False) as ds:
        depthw = np.asarray(ds["depthw"].values, dtype=np.float64)

    difference = max(diff_bathy, diff_temp)
    if not difference <= COORDINATE_TOLERANCE:
        raise ValueError(
            "the bathymetry, sea-ice and temperature files are not on the "
            f"same grid: coordinates differ by up to {difference:.2e} degrees")
    if depthw.shape != deptht.shape:
        raise ValueError("depthw and deptht have different lengths")
    if not (np.all(depthw[:-1] < deptht[:-1]) and np.all(deptht[:-1] < depthw[1:])):
        raise ValueError("deptht levels do not lie between the depthw bounds")
    grid.update({"bathymetry": hdept, "deptht": deptht, "depthw": depthw,
                 "coordinate_difference": difference,
                 "files": {"ice": ice_file.name, "temp": temp_file.name}})
    return grid


def build_cells(grid, table):
    """
    Check the release-cell table against the grid and add the neighbours.

    Returns a dict of arrays over the cell dimension, sorted by (j, i):
      cell_y, cell_x   0-based grid indices
      lon, lat         cell centre (nav_lon, nav_lat)
      bathymetry       hdept (m)
      is_release       1 for a release cell, 0 for a neighbour
      n_release        particles generated per release day (0 for a neighbour)
    """
    ny, nx = grid["bathymetry"].shape
    j, i = table["j"], table["i"]
    if j.min() < 0 or j.max() >= ny or i.min() < 0 or i.max() >= nx:
        raise ValueError(f"table indices fall outside the {ny} x {nx} grid")
    dlon = np.abs(grid["nav_lon"][j, i] - table["lon"]).max()
    dlat = np.abs(grid["nav_lat"][j, i] - table["lat"]).max()
    ddep = np.abs(grid["bathymetry"][j, i] - table["depth_m"]).max()
    # The table is rounded to 5 decimals (degrees) and 1 decimal (metres).
    if dlon > 2e-5 or dlat > 2e-5 or ddep > 0.051:
        raise ValueError(
            "the release-cell table does not match the grid: largest "
            f"difference {dlon:.2e} deg (lon), {dlat:.2e} deg (lat), "
            f"{ddep:.3f} m (depth)")
    depth = grid["bathymetry"][j, i]
    if depth.min() < DEPTH_MIN or depth.max() > DEPTH_MAX:
        raise ValueError("release cells outside the 1000-2000 m band")

    n_release = np.zeros((ny, nx), dtype=np.int64)
    n_release[j, i] = table["n_release"]
    is_release = np.zeros((ny, nx), dtype=bool)
    is_release[j, i] = True
    keep = np.zeros((ny, nx), dtype=bool)
    for dj, di in NEIGHBOUR_OFFSETS:
        jj, ii = j + dj, i + di
        ok = (jj >= 0) & (jj < ny) & (ii >= 0) & (ii < nx)
        if not ok.all():
            raise ValueError("a release cell lies on the edge of the grid; "
                             "its neighbours cannot be extracted")
        keep[jj, ii] = True
    cell_y, cell_x = np.nonzero(keep)          # sorted by (j, i)
    return {
        "cell_y": cell_y.astype(np.int32),
        "cell_x": cell_x.astype(np.int32),
        "lon": grid["nav_lon"][cell_y, cell_x].astype(np.float32),
        "lat": grid["nav_lat"][cell_y, cell_x].astype(np.float32),
        "bathymetry": grid["bathymetry"][cell_y, cell_x].astype(np.float32),
        "is_release": is_release[cell_y, cell_x].astype(np.int8),
        "n_release": n_release[cell_y, cell_x].astype(np.int16),
    }


def number_of_levels(grid, cells):
    """
    Number of levels to extract: every level whose top lies above the sea
    floor of the deepest release cell, so that the column of a release cell
    reaches its sea floor.
    """
    deepest = float(cells["bathymetry"][cells["is_release"] == 1].max())
    return int(np.count_nonzero(grid["depthw"] < deepest))


def cell_lookup(cell_y, cell_x, shape):
    """(ny, nx) int32 map: index along the cell dimension, -1 elsewhere."""
    lookup = np.full(shape, -1, dtype=np.int32)
    lookup[cell_y, cell_x] = np.arange(len(cell_y), dtype=np.int32)
    return lookup


def locate_cells(lon, lat, grid_lon, grid_lat, lookup):
    """
    Index along the cell dimension of the grid point nearest to each
    position, -1 where that point is not an extracted cell. This is the
    sampling rule of M1 (sea_ice._sample_nearest), applied to the lookup map.
    """
    return sea_ice._sample_nearest(lookup, grid_lon, grid_lat,
                                   np.asarray(lon, dtype=np.float64),
                                   np.asarray(lat, dtype=np.float64))


def locate_levels(depth, depth_bounds):
    """
    Level holding each depth: k such that depth_bounds[k] <= depth <
    depth_bounds[k + 1]; -1 above the first bound or at and below the last.
    """
    depth = np.asarray(depth, dtype=np.float64)
    k = np.searchsorted(depth_bounds, depth, side="right") - 1
    k[(k < 0) | (k >= len(depth_bounds) - 1) | ~np.isfinite(depth)] = -1
    return k
