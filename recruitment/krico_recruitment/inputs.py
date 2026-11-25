"""
Checks of the inputs of the classification of one cohort: the trajectory
file and the release-cell file of its spawning year.

check_cohort is called by process_cohort.py before anything is computed,
and by check_inputs.py on every cohort before any job is submitted, so
that the two can never disagree on what an acceptable input is. Every
failure raises ValueError (or FileNotFoundError) with the file name.
"""

import re
from datetime import date, datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr

from krico_recruitment import nsidc_ice, pre_release, release_cells, sea_ice

# Float variables read along each trajectory, as written by run.py of
# krico-templates: (trajectory, obs), NaN after deletion.
TRAJECTORY_VARIABLES = ("lon", "lat", "z", "temperature", "bathymetry",
                        "sea_ice_area_fraction")

# Daily values per trajectory: the release day and TRACKING_DAYS more.
N_OBS = release_cells.TRACKING_DAYS + 1

_NAME = re.compile(r"^(\d{4})_(\d{2})_(\d{2})\.nc$")


def cohort_day(path):
    """Release day (datetime.date) from a file name YYYY_MM_DD.nc."""
    match = _NAME.match(Path(path).name)
    if not match:
        raise ValueError(f"{Path(path).name}: not named YYYY_MM_DD.nc")
    return date(*map(int, match.groups()))


def check_ice_source(ice_source, ice_dir):
    """The ice source is one the classification knows, with its folder."""
    if ice_source != "glorys12":
        if ice_source not in nsidc_ice.SOURCES:
            raise ValueError(f"unknown ice source {ice_source!r}: glorys12, "
                             + ", ".join(nsidc_ice.SOURCES))
        if not ice_dir:
            raise ValueError(f"ice source {ice_source} needs the folder of "
                             "its monthly files ($"
                             f"{nsidc_ice.SOURCES[ice_source]['env']})")


def check_cohort(trajectory_path, release_cells_dir, statics=None,
                 read_positions=True, ice_source="glorys12", ice_dir=None):
    """
    Check the inputs of one cohort and return what the classification
    needs from them.

    Parameters
    ----------
    trajectory_path : path of YYYY_MM_DD.nc
    release_cells_dir : folder of glorys12_release_cells_YYYY.nc
    statics : dict, optional
        Cache of pre_release.open_static results by path, for a caller that
        checks many cohorts of the same spawning year.
    read_positions : bool
        False skips the checks that read the release positions (headers
        and release time only).
    ice_source, ice_dir : str
        "glorys12" (the sea ice of the trajectory file and of the
        release-cell file), or a product of nsidc_ice.SOURCES with the
        folder of its monthly files, which must then hold every day from
        the spawning day to the last tracked day.

    Returns
    -------
    dict with release_day, n_particles, n_obs, static (pre_release
    .open_static) and, if read_positions, cell (index of the release cell
    of every particle) and depth (release depths).
    """
    check_ice_source(ice_source, ice_dir)
    path = Path(trajectory_path)
    name = path.name
    day = cohort_day(path)
    if not path.is_file():
        raise FileNotFoundError(f"{path} is missing")

    # Opened as the classification opens it (io.read_trajectory_file).
    with xr.open_dataset(path, decode_timedelta=False) as ds:
        for dim in ("trajectory", "obs"):
            if dim not in ds.sizes:
                raise ValueError(f"{name}: no dimension {dim}")
        n = ds.sizes["trajectory"]
        n_obs = ds.sizes["obs"]
        if n == 0:
            raise ValueError(f"{name}: no particle")
        if n_obs != N_OBS:
            raise ValueError(f"{name}: {n_obs} daily values per trajectory, "
                             f"expected {N_OBS}")
        for var in TRAJECTORY_VARIABLES + ("time",):
            if var not in ds.variables:
                raise ValueError(f"{name}: variable {var} is missing")
            if ds[var].dims != ("trajectory", "obs"):
                raise ValueError(f"{name}: {var} has dimensions "
                                 f"{ds[var].dims}, expected "
                                 "('trajectory', 'obs')")
        for var in TRAJECTORY_VARIABLES:
            if ds[var].dtype.kind != "f":
                raise ValueError(f"{name}: {var} is read as {ds[var].dtype}, "
                                 "expected floating point")
        if ds["time"].dtype.kind != "M":
            raise ValueError(f"{name}: time is read as {ds['time'].dtype}, "
                             "expected dates")

        # Release time: the same for every particle, at 12:00 on the day
        # of the file name. Daily means are stamped 12:00, so this is what
        # makes day k of a trajectory the daily mean of release day + k.
        stamps = np.unique(ds["time"].isel(obs=0).values)
        expected = np.datetime64(datetime.combine(day, datetime.min.time())
                                 + timedelta(hours=12), "ns")
        if stamps.size != 1 or stamps[0] != expected:
            shown = ", ".join(str(s) for s in stamps[:3])
            raise ValueError(f"{name}: released at {shown}"
                             + (" ..." if stamps.size > 3 else "")
                             + f", expected {expected} for every particle")

        if read_positions:
            first = ds.isel(obs=0)
            lon = np.asarray(first["lon"].values, dtype=np.float64)
            lat = np.asarray(first["lat"].values, dtype=np.float64)
            depth = np.asarray(first["z"].values, dtype=np.float64)

    # Release-cell file of the spawning year: the 24 days before release
    # (M1, M3) and the release day (day-0 check) must be inside it.
    cells_path = pre_release.release_cells_path(release_cells_dir, day)
    if statics is not None and cells_path in statics:
        static = statics[cells_path]
    else:
        static = pre_release.open_static(cells_path)
        if statics is not None:
            statics[cells_path] = static
    pre_release.day_index(static, day - timedelta(days=sea_ice.SPAWNING_OFFSET_DAYS))
    pre_release.day_index(static, day)

    if ice_source != "glorys12":
        nsidc_ice.check_coverage(
            ice_dir, day - timedelta(days=sea_ice.SPAWNING_OFFSET_DAYS),
            day + timedelta(days=n_obs - 1),
            nsidc_ice.SOURCES[ice_source]["pattern"])

    result = {"release_day": day, "n_particles": n, "n_obs": n_obs,
              "static": static}
    if read_positions:
        bad = ~(np.isfinite(lon) & np.isfinite(lat))
        if bad.any():
            raise ValueError(f"{name}: {int(bad.sum())} particles have no "
                             "position at release")
        cell = pre_release.locate(static, lon, lat)
        outside = cell < 0
        if outside.any():
            k = int(np.flatnonzero(outside)[0])
            raise ValueError(
                f"{name}: {int(outside.sum())} particles have their nearest "
                f"grid point outside the cells of {cells_path.name} (first: "
                f"lon {lon[k]:.5f}, lat {lat[k]:.5f})")
        result.update({"cell": cell, "depth": depth})
    return result
