import numpy as np
import xarray as xr
import pandas as pd
import dask
from tqdm import tqdm
from pyresample import bilinear
from pyresample.bucket import BucketResampler
from sklearn.neighbors import BallTree

import abstemp
from abstemp.seagrid import cmip6

MODEL = "ec_earth3_cc"


def max_month():
    """Find the per-year maximum monthly SST and the month it occurs.

    For each year in the EC-Earth3-CC dataset, iterates over monthly SST
    fields and records the highest SST value and the corresponding month at
    every grid cell.

    Returns
    -------
    maxarr : numpy.ndarray
        3-D array of shape (n_years, 3600, 7200) with the maximum monthly SST
        (°C) for each year.  Cells that were never updated are NaN.
    monarr : numpy.ndarray of uint16
        3-D array of shape (n_years, 3600, 7200) with the calendar month
        (1–12) at which the maximum SST was reached.
    """
    ds = cmip6.center_on_gmt(cmip6.open_dataset("ec_earth3_cc"))

    shape = np.array(ds.sst.shape)
    shape[0] = shape[0]//12 + 1
    maxarr = np.zeros(shape) - 6
    monarr = np.zeros(shape, np.uint16)

    for mn, sst in enumerate(ds.sst):
        yr = mn // 12
        mask = sst.values > maxarr[yr]
        maxarr[yr][mask] = sst.values[mask]
        monarr[yr][mask] = np.mod(mn, 12) + 1
    maxarr[maxarr==-6] = np.nan
    return maxarr, monarr

def min_month():
    """Find the per-year minimum monthly SST and the month it occurs.

    For each year in the EC-Earth3-CC dataset, iterates over monthly SST
    fields and records the lowest SST value and the corresponding month at
    every grid cell.

    Returns
    -------
    minarr : numpy.ndarray
        3-D array of shape (n_years, 3600, 7200) with the minimum monthly SST
        (°C) for each year.  Cells that were never updated are NaN.
    monarr : numpy.ndarray of uint16
        3-D array of shape (n_years, 3600, 7200) with the calendar month
        (1–12) at which the minimum SST was reached.
    """
    ds = cmip6.center_on_gmt(cmip6.open_dataset("ec_earth3_cc"))

    shape = np.array(ds.sst.shape)
    shape[0] = shape[0]//12 + 1
    minarr = np.zeros(shape) - 6
    monarr = np.zeros(shape, np.uint16)

    for mn, sst in enumerate(ds.sst):
        yr = mn // 12
        mask = sst.values < minarr[yr]
        minarr[yr][mask] = sst.values[mask]
        monarr[yr][mask] = np.mod(mn, 12) + 1
    minarr[minarr==-6] = np.nan
    return minarr, monarr



def clim_max_month():
    """Compute the climatological maximum monthly SST and the warmest month.

    Averages EC-Earth3-CC SST for each calendar month across all available
    years, then takes the maximum across months at every grid cell to produce
    a climatological peak SST and the corresponding month index.

    Returns
    -------
    maxarr : numpy.ndarray
        2-D array of shape (3600, 7200) with the climatological maximum SST
        (°C) across all calendar months.
    monarr : numpy.ndarray
        2-D array of shape (3600, 7200) with the calendar month (1–12) of the
        climatological maximum.
    """
    ds = cmip6.center_on_gmt(cmip6.open_dataset("ec_earth3_cc"))

    shape = np.array(ds.sst.shape)
    shape[0] = 12
    sumarr = np.zeros(shape)
    cntarr = np.zeros(shape, np.uint16)

    for mn, sst in enumerate(ds.sst):
        month_idx = mn % 12
        sumarr[month_idx] += sst.values
        cntarr[month_idx] += np.isfinite(sst.values)
    maxarr = np.max(sumarr / cntarr, axis=0)
    monarr = np.argmax(sumarr / cntarr, axis=0) + 1
    return maxarr, monarr


def nearest(ds, return_dist=False):
    """Find the nearest mintmat region for each OSTIA full-resolution pixel.

    Builds a haversine BallTree on the mintmat region centroids (from
    ``data/mintmat_2001-2009.nc``) and queries it with every pixel in the
    monthly OSTIA SST dataset to obtain a region assignment for each pixel.

    Parameters
    ----------
    return_dist : bool, optional
        If True, also return the haversine distances (in radians) to the
        nearest region centroid.  Default is False.

    Returns
    -------
    ij : numpy.ndarray
        1-D array of region indices for each OSTIA pixel (flattened).
    dist : numpy.ndarray
        Haversine distances in radians (only returned when
        ``return_dist=True``).
    """
    gl = abstemp.open_mintmat_ds()
    latlon = np.deg2rad(np.array([gl.reglat,gl.reglon]).T)
    tree = BallTree(latlon[1:,:], metric="haversine")

    if "lons" in ds:
        lons,lats = ds.lons.values, ds.lats.values
    else:
        lons,lats = np.meshgrid(ds.lon,ds.lat)
    latlon2 = np.deg2rad(np.array([lats.flatten(),lons.flatten()]).T)
    dist,ij = tree.query(latlon2)
    if return_dist:
        return dist,ij
    return ij


def to_reg(ds, ij=None):
    """Average CMIP6 SST values into mintmat regions.

    Maps each model grid pixel to the nearest mintmat region centroid and
    computes the mean SST within each region.  Accepts both a full monthly
    dataset (``sst`` + ``time``) and a pre-computed max-month dataset
    (``maxarr``, no ``time`` dimension).

    Parameters
    ----------
    ds : xarray.Dataset
        Either a monthly dataset with a ``sst`` variable and ``time``
        dimension, or a pre-computed statistics dataset with a ``maxarr``
        variable (as returned by :func:`abstemp.data.max_min_month`).
    ij : numpy.ndarray, optional
        Precomputed region indices from :func:`nearest`.  Computed
        on-the-fly if not provided.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by region number (starting at 0) with a single
        column ``sst`` containing the mean SST (°C) for each region.
        Region 0 is set to NaN.
    """
    ij = nearest(ds) if ij is None else ij
    arr = ds.maxarr.values if "maxarr" in ds else ds.sst.max(dim="time").values
    df = pd.DataFrame({"ij":np.squeeze(ij+1), "sst":arr.flatten()})
    svec = df.groupby("ij").mean().reset_index()
    svec = pd.concat([pd.DataFrame({"ij":0, "sst":[np.nan]}),svec], axis=0)
    svec.set_index("ij", inplace=True)
    svec = svec.reindex(range(svec.index.max() + 1))
    return svec
