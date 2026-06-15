import warnings
from pathlib import Path

import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

import pandas as pd
import xarray as xr
import haversine as hs

from abstemp.seagrid import ostia
from abstemp.seagrid import gridtools
from abstemp.seagrid import cmip6

def open_dataset(
    dtm1: str = "1985-01-01",
    dtm2: str = "1995-01-01",
    data_var: str = "sst",
    nrt: bool = False,
) -> xr.Dataset:
    """Load monthly OSTIA SST fields and return them as a time-stacked dataset.

    Parameters
    ----------
    dtm1 : str, optional
        Start date (inclusive) in ``"YYYY-MM-DD"`` format.
        Default is ``"1985-01-01"``.
    dtm2 : str, optional
        End date (inclusive) in ``"YYYY-MM-DD"`` format.
        Default is ``"1995-01-01"``.
    data_var : str, optional
        Variable name to extract from each monthly OSTIA file.
        Default is ``"sst"``.
    nrt : bool, optional
        If True, load from the near-real-time OSTIA product instead of
        the reanalysis.  Default is False.

    Returns
    -------
    xarray.Dataset
        Dataset with dimension ``time`` (one step per month-start in
        ``[dtm1, dtm2]``) and spatial dimensions ``lat`` / ``lon``
        at 0.05° resolution.
    """
    dtmvec = pd.date_range(dtm1, dtm2, freq="MS")
    dalist = []
    for dtm in dtmvec:
        ds = ostia.open_dataset(dtm=dtm, timetype="M", nrt=nrt)
        ds = ds.assign_coords(lat=np.linspace(-89.975, 89.975,3600),
                              lon=np.linspace(-180,180,7200))
        ds = ds.assign_coords({"time":dtm})
        dalist.append(ds[data_var])
    return xr.concat(dalist, dim="time").to_dataset()


def monthly_climatology(ds: xr.Dataset, data_var: str = "sst") -> xr.Dataset:
    """Compute a 12-month SST climatology by averaging over all available years.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with a ``time`` dimension and a variable named *data_var*.
    data_var : str, optional
        Name of the SST variable.  Default is ``"sst"``.

    Returns
    -------
    xarray.Dataset
        Dataset with dimension ``month`` (1–12) containing the mean SST
        for each calendar month across all years in *ds*.
    """
    if hasattr(ds.indexes["time"], "to_datetimeindex"):
        with warnings.catch_warnings(action="ignore"):
            months = ds.indexes["time"].to_datetimeindex().month
    else:
        months = pd.to_datetime(ds.time).month
    mnlist = []
    for month in range(1,13):
        #mask = pd.to_datetime(ds.time).month == month
        with warnings.catch_warnings(action="ignore"):
            mask = months == month
        clim = ds[data_var][mask,:,:].mean(dim="time")
        clim = clim.assign_coords({"month":month})
        mnlist.append(clim)
    return xr.concat(mnlist, dim="month").to_dataset()

def max_month(ds: xr.Dataset, data_var: str = "sst") -> np.ndarray:
    """Return the calendar month (1–12) of the climatological maximum SST at each pixel.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset passed to :func:`monthly_climatology`.
    data_var : str, optional
        SST variable name.  Default is ``"sst"``.

    Returns
    -------
    numpy.ndarray
        2-D integer array (lat × lon) with values 1–12 indicating the
        warmest calendar month at each grid cell.
    """
    cl = monthly_climatology(ds, data_var=data_var)
    return np.argmax(cl[data_var].data, axis=0) + 1

def area(ds: xr.Dataset) -> np.ndarray:
    """Compute grid-cell area (km²) for each pixel in a dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with either 1-D ``lat``/``lon`` coordinates (regular grid)
        or 2-D ``lats``/``lons`` data variables (curvilinear grid).

    Returns
    -------
    numpy.ndarray
        2-D array of grid-cell areas in km².
    """
    if "lat" in ds:
        lons,lats = np.meshgrid(ds.lon,ds.lat)
    else:
        lons,lats = ds.lons.data, ds.lats.data
    return gridtools.area(lats=lats, lons=lons)

def warmest_monthly_sst(ds: xr.Dataset, data_var: str = "sst") -> np.ndarray:
    """Return the climatological mean SST at the warmest month for each pixel.

    Computes a 12-month climatology via :func:`monthly_climatology` and
    returns the maximum across months at each grid cell.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with a ``time`` dimension and an ``sst`` variable spanning
        at least one full year.
    data_var : str, optional
        SST variable name.  Default is ``"sst"``.

    Returns
    -------
    numpy.ndarray
        2-D array (lat × lon) of climatological mean SST (°C) at the
        warmest calendar month for each pixel.  Ocean-free cells are NaN.
    """
    cl = monthly_climatology(ds, data_var=data_var)
    return cl[data_var].max(dim="month").values


def haversine_area(ds: xr.Dataset) -> np.ndarray:
    """Compute grid-cell area (km²) using haversine distances between adjacent cell centres.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with 1-D ``lat``/``lon`` coordinates or 2-D ``lats``/``lons``
        data variables.

    Returns
    -------
    numpy.ndarray
        2-D array of approximate cell areas in km².  Edge cells are estimated
        by averaging the neighbouring cell's dimension.
    """
    if "lat" in ds:
        lons,lats = np.meshgrid(ds.lon,ds.lat)
    else:
        lons,lats = ds.lons.data.copy(), ds.lats.data.copy()
    lons[lons>180] = lons[lons>180]-360

    pos1 = np.array([lats[:,:-1].flatten(), lons[:,:-1].flatten()]).T
    pos2 = np.array([lats[:,1: ].flatten(), lons[:,1: ].flatten()]).T
    xdist = hs.haversine_vector(pos1, pos2).reshape(lats[:,1:].shape)

    pos1 = np.array([lats[:-1,:].flatten(), lons[:-1,:].flatten()]).T
    pos2 = np.array([lats[1:, :].flatten(), lons[1:, :].flatten()]).T
    ydist = hs.haversine_vector(pos1, pos2).reshape(lats[1:,:].shape)

    area = np.zeros(lats.shape)
    area[:-1,:] = ydist
    area[:,:-1] = area[:,:-1] * xdist
    area[:,-1] = (area[:,0] + area[:,-2]) / 2
    area[-1,:] = (area[0,:] + area[-2,:]) / 2

    return area


def save_ostia_files() -> None:
    """Generate and save OSTIA warmest-month SST files to the package data directory.

    Writes ``ostia_sst_1985-1990.nc`` and ``ostia_sst_2019-2023.nc`` by
    loading monthly OSTIA fields via :func:`open_dataset`.  Requires access
    to the OSTIA archive through njord.

    For CMIP6 2095–2100 files use
    :func:`abstemp.data.download.maxmonsst_fields` instead.

    Returns
    -------
    None
    """
    datadir = Path(__file__).parent / "data"
    datadir.mkdir(exist_ok=True)

    ds = open_dataset("1985-01-01", "1990-12-31", "sst", nrt=False)
    ds.to_netcdf(datadir / "ostia_sst_1985-1990.nc")

    ds = open_dataset("2019-01-01", "2023-12-31", "sst", nrt=True)
    ds.to_netcdf(datadir / "ostia_sst_2019-2023.nc")
