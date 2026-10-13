"""
Geometry of the NSIDC 25 km polar-stereographic south grid (EPSG:3412), as
used by the G02202 sea-ice concentration CDR.

The grid is 316 columns by 332 rows of 25 km cells. The projection is polar
stereographic on the Hughes 1980 ellipsoid, true at 70S, with 0E pointing up
(+y) and 90E pointing right (+x). The cell-edge origin and spacing come from
the GeoTransform attribute of the product files:
    "-3950000 25000 0 4350000 0 -25000"
i.e. the top-left corner of cell (row 0, col 0) is at x = -3950000,
y = +4350000 metres, columns step +25 km in x, rows step -25 km in y.

Two functions:
    lonlat_to_xy(lon, lat)      -> x, y in metres (Snyder 1987, eq. 21-33/34)
    lonlat_to_index(lon, lat)   -> row, col of the containing cell, or -1
                                   for points off the grid

Both take and return numpy arrays of any shape. Nearest-cell lookup is then
``field[row, col]``, which is how the sampler reads concentration at a
trajectory position: the 25 km product is a cell average and is read as such,
without interpolation.

The formulas are verified against the ancillary file's own latitude/longitude
arrays by ``scripts/preprocess_nsidc.py --check-grid`` (every cell centre must
map back to its own row/col) and, offline, against pyproj.
"""

import numpy as np

# Hughes 1980 ellipsoid, from the product's crs attributes.
A = 6378273.0
INV_F = 298.279411123064
E = np.sqrt(2.0 / INV_F - 1.0 / INV_F**2)

LAT_TS = -70.0      # standard parallel (degrees)
LON_0 = 0.0         # central meridian (degrees)

# Cell-edge origin and spacing (metres), from the GeoTransform.
X_EDGE_0 = -3950000.0
Y_EDGE_0 = 4350000.0
DX = 25000.0
NCOL = 316
NROW = 332

# Cell centres, matching the x and y coordinate variables of the files.
X_CENTRES = X_EDGE_0 + DX * (np.arange(NCOL) + 0.5)
Y_CENTRES = Y_EDGE_0 - DX * (np.arange(NROW) + 0.5)


def _t(phi):
    """Snyder's t(phi) for the ellipsoidal polar stereographic."""
    sin_phi = np.sin(phi)
    return np.tan(np.pi / 4 - phi / 2) / (
        (1 - E * sin_phi) / (1 + E * sin_phi)) ** (E / 2)


def lonlat_to_xy(lon, lat):
    """Project geographic coordinates (degrees) to EPSG:3412 x, y (metres).

    South polar aspect: Snyder's north-polar formulas with the signs of
    latitude, longitude and the resulting x, y reversed.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    phi = np.radians(-lat)
    lam = np.radians(-lon)
    lam0 = np.radians(-LON_0)
    phi_c = np.radians(-LAT_TS)

    m_c = np.cos(phi_c) / np.sqrt(1 - E**2 * np.sin(phi_c) ** 2)
    rho = A * m_c * _t(phi) / _t(phi_c)
    x = rho * np.sin(lam - lam0)
    y = -rho * np.cos(lam - lam0)
    return -x, -y


def xy_to_index(x, y):
    """Row and column of the cell containing (x, y); -1 where off the grid."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    col = np.floor((x - X_EDGE_0) / DX).astype(np.int64)
    row = np.floor((Y_EDGE_0 - y) / DX).astype(np.int64)
    inside = (col >= 0) & (col < NCOL) & (row >= 0) & (row < NROW)
    col = np.where(inside, col, -1)
    row = np.where(inside, row, -1)
    return row, col


def lonlat_to_index(lon, lat):
    """Row and column of the grid cell containing each (lon, lat); -1 off-grid."""
    x, y = lonlat_to_xy(lon, lat)
    return xy_to_index(x, y)
