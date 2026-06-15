"""Retrieve and store daily ocean reanalysis data from Mercator Ocean."""

import pathlib

import numpy as np
import xarray as xr
import pandas as pd
import requests

from xarray.backends import PydapDataStore
from dask.diagnostics import ProgressBar


HOST = "dap2://tds.mercator-ocean.fr"
#PRODUCT = "psy4v3r1"
PRODUCT = "glorys12v1"

vardict = dict(
    uvel=["U", "vozocrtx"],
    vvel=["V", "vomecrty"],
    wvel=["W", "vovecrtz"],
    temp=["T", "votemper"],
    salt=["S", "vosaline"]
)

def datadir(dtm="2020-01-01"):
    """Return the local data directory for a given dat, creating it if needed.

    Parameters
    ----------
    dtm : str or datetime-like
        Date string or datetime-like object.

    Returns
    -------
    pathlib.Path
        Path to the year-specific data directory.
    """
    dtm = pd.to_datetime(dtm)
    path = f"/data/raid/mercator/{PRODUCT}/{dtm.year}/"
    path = pathlib.Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path

def filename(dtm="2020-01-01", data_var="uvel"):
    """Build the Zarr filename for a given dat and variable.

    Parameters
    ----------
    dtm : str or datetime-like
        Date string or datetime-like object.
    data_var : str
        Variable key (one of 'uvel', 'vvel', 'wvel', 'temp', 'salt').

    Returns
    -------
    str
        Filename string in the form ``GLORYS12V1_<date>_grid<X>.zarr``.
    """
    def dtmstr(dtm):
        return f"{dtm.year}{dtm.month:02}{dtm.day:02}"
    dtm = pd.to_datetime(dtm)
    chr = vardict[data_var][0]
    if PRODUCT == "psy4v3r1":
        #ext-PSY4V3R1_1dAV_20121123_20121124_gridW_R20121128.nc
        dtmplus = dtmstr(dtm + pd.DateOffset(1, "D"))
        fn = (
           "ext-PSY4V3R1_1dAV_" +
           f"{dtmstr(dtm)}_{dtmplus}_grid{chr}_" +
           f"R{dtmstr(pd.Timestamp.now())}.zarr")
    else:
        fn = f"{PRODUCT}_{dtm.date()}_grid{chr}.zarr"
    return fn

def retrieve(dtm="2020-01-01", data_var="uvel", force=False):
    """Download a single daily field from the Mercator THREDDS server and save it as Zarr.

    Opens the remote dataset via OPeNDAP, selects the time step closest to `dtm`,
    and writes the result to a local Zarr store.

    Parameters
    ----------
    dtm : str or datetime-like
        Date string or datetime-like object for the desired time step.
    data_var : str
        Variable key (one of 'uvel', 'vvel', 'wvel', 'temp', 'salt').
        force : bool
            If True, re-download the file even if it already exists locally.
    """
    dtm = pd.to_datetime(dtm)
    fname = filename(dtm=dtm, data_var=data_var)
    ddir  = datadir(dtm=dtm)
    if (ddir / fname).exists() and not force:
       return
    chunks = {'time_counter':1, 'deptht':-1, 'y':1000, 'x':1000}
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(
        max_retries=requests.adapters.Retry(total=5, backoff_factor=1),
        pool_connections=8,
        pool_maxsize=16
    )
    session.mount("http://", adapter)
    path = f"/thredds/dodsC/{PRODUCT}-daily-grid"
    url = HOST + path + vardict[data_var][0]
    store = PydapDataStore.open(url, session=session)
    ds = xr.open_dataset(store, chunks=chunks).astype('float32')

    dtm_index = pd.to_datetime(ds.time_counter).normalize()
    iloc = np.argmin(np.abs(dtm_index-dtm))
    da = ds[vardict[data_var][1]].isel(time_counter=iloc).to_dataset()

    write_job = da.to_zarr(str(ddir / fname), mode="w", compute=False)
    with ProgressBar():
        write_job.compute()

def open_dataset(dtm="2020-01-01", data_var="uvel",  force=False):
    """Open a mercator dataset, downloading it first if necessary.

    Parameters
    ----------
    dtm : str or datetime-like
        Date string or datetime-like object for the desired time step.
    data_var : str
        Variable key (one of 'uvel', 'vvel', 'wvel', 'temp', 'salt').
    force : bool
        If True, re-download the file even if it already exists locally.

    Returns
    -------
    xarray.Dataset
        Dataset containing the requested variable for the given date.
    """
    dtm = pd.to_datetime(dtm)
    fname = (datadir(dtm=dtm) /
        filename(dtm=dtm, data_var=data_var))
    if not fname.exists() or force:
        retrieve(dtm=dtm, data_var=data_var, force=force)
    ds = xr.open_dataset(fname)
    return ds

def open_dataarray(dtm="2020-01-01", data_var="uvel",  force=False):
    """Open a mercator dataarray, downloading it first if necessary.

    Parameters
    ----------
    dtm : str or datetime-like
        Date string or datetime-like object for the desired time step.
    data_var : str
        Variable key (one of 'uvel', 'vvel', 'wvel', 'temp', 'salt').
    force : bool
        If True, re-download the file even if it already exists locally.

    Returns
    -------
    xarray.Dataset
        Dataset containing the requested variable for the given date.
    """
    ds = open_dataset(dtm=dtm, data_var=data_var, force=force)
    return ds[data_var]
