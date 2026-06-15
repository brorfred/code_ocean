"""OSTIA daily and monthly SST data source.

Provides access to the Copernicus OSTIA (Operational Sea Surface
Temperature and Ice Analysis) product.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm

from . import config
from . import copernicus
from .utils.grids import setup_plate_carree_grid

settings = config.settings.from_env("ostia")

__seagrid_data_source__ = True

_DATADIR = Path(__file__).parent.parent / "data" / "ostia"
_DATADIR.mkdir(parents=True, exist_ok=True)


def datadir():
    return str(_DATADIR)

grid = setup_plate_carree_grid(-180, 180, -90, 90, 7200, 3600)


def vprint(string, level=30):
    """Print text to stdout.

    Parameters
    ----------
    text : str
    """
    print(string)


def setup_grid(**kw):
    """Create a pyresample grid from the OSTIA dataset.

    Returns
    -------
    pyresample.geometry.GridDefinition
    """
    return setup_latlon_grid(open_dataset(), **kw)


def filename(dtm="2010-01-01", timetype="day", version=None):
    """Generate an OSTIA filename.

    Parameters
    ----------
    dtm : str or datetime-like
        Date.
    timetype : str, optional
        'day' for daily, or any string containing 'm' for monthly.
    version : optional
        Unused, reserved.

    Returns
    -------
    str
    """
    dtm = pd.to_datetime(dtm)
    if "m" in timetype.lower():
        return settings["filestamp"].format(dtm=dtm.strftime("%Y%m") + "_MON_")
    else:
        return settings["filestamp"].format(dtm=dtm.strftime("%Y%m%d"))


def dtm_vec():
    """Return a date range from settings.

    Returns
    -------
    pandas.DatetimeIndex
    """
    return pd.date_range(
        settings["dtm_min"], settings["dtm_max"], freq=settings["dtm_freq"]
    )


def download(dtm="2010-01-01"):
    """Download an OSTIA file via Copernicus Marine.

    Parameters
    ----------
    dtm : str or datetime-like
        Date to download.
    """
    dtm = pd.to_datetime(dtm)
    copernicus.get(dtm=dtm, settings=settings)


def _harmonize_coordinates(ds):
    """Standardize OSTIA dataset coordinates to uniform spacing.

    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset.

    Returns
    -------
    xarray.Dataset
        Dataset with harmonized lat/lon coordinates.
    """
    return ds.assign_coords(
        lat=np.linspace(-89.975, 89.975, 3600), lon=np.linspace(-180, 180, 7200)
    )


def open_dataset(
    dtm="2010-01-01", data_var="sst", version=2.0, timetype="daily", nrt=False, **xrargs
):
    """Open an OSTIA dataset, downloading if necessary.

    SST is converted from Kelvin to Celsius and the variable is
    renamed from 'analysed_sst' to 'sst'.

    Parameters
    ----------
    dtm : str or datetime-like
        Date.
    data_var : str, optional
        Variable name.
    version : float, optional
        Dataset version.
    timetype : str, optional
        'daily' for daily data, or string containing 'm' for monthly.
    nrt : bool, optional
        If True, use near-real-time data source.
    **xrargs : dict
        Additional arguments passed to ``xarray.open_dataset``.

    Returns
    -------
    xarray.Dataset
    """
    if nrt and (pd.to_datetime(dtm) >= pd.to_datetime("2007-01-01")):
        settings.setenv("ostia-nrt")
    else:
        settings.setenv("ostia")
    if "m" in timetype.lower():
        return _open_monthly_dataset(
            dtm=dtm, version=version, data_var=data_var, nrt=nrt
        )
    local_filename = os.path.join(datadir(), filename(dtm=dtm, version=version))
    # vprint(local_filename)
    if not os.path.isfile(local_filename):
        download(dtm=dtm)
    ds = xr.open_dataset(local_filename, **xrargs)
    ds = ds.rename({"analysed_sst": "sst"})
    ds["sst"] = ds["sst"] - 273.15
    ds.attrs["seagrid_note"] = (
        "analysed_sst is renamed to sst and converted from Kelvin to Celsius"
    )
    return _harmonize_coordinates(ds)


def _open_monthly_dataset(dtm="2010-01-01", version=2.0, data_var="sst", nrt=False):
    """Aggregate daily OSTIA data into a monthly mean.

    Parameters
    ----------
    dtm : str or datetime-like
        Any date within the target month.
    version : float, optional
        Dataset version.
    data_var : str, optional
        Variable to average.
    nrt : bool, optional
        Use near-real-time source.

    Returns
    -------
    xarray.Dataset
        Monthly-averaged dataset.
    """
    dtm = pd.to_datetime(dtm) if isinstance(dtm, str) else dtm
    fn = os.path.join(datadir(), filename(dtm=dtm, timetype="month"))
    if os.path.isfile(fn):
        return xr.open_dataset(fn)
    dtm1 = dtm.replace(day=1)
    dtm2 = dtm1 + pd.DateOffset(months=1) - pd.DateOffset(days=1)
    dslist = []
    for dtm in tqdm(pd.date_range(dtm1, dtm2)):
        dslist.append(
            open_dataset(dtm=dtm, version=version, nrt=nrt, data_var=data_var)
        )
        dslist.append(
            open_dataset(dtm=dtm, version=version, nrt=nrt, data_var=data_var)
        )
    ds = xr.concat(dslist, dim="time")
    ds = ds.mean(dim="time", skipna=True)
    ds.to_netcdf(fn)
    return _harmonize_coordinates(ds)
