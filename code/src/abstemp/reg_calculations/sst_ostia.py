import numpy as np
import xarray as xr
import pandas as pd
import dask
from tqdm import tqdm
from pyresample import bilinear
from pyresample.bucket import BucketResampler
from sklearn.neighbors import BallTree

from abstemp.seagrid import ostia
from abstemp.seagrid.utils import grids

import abstemp

def mean_month():
    """Compute the temporal mean SST over all months from 2001 to 2009.

    Loads monthly OSTIA SST fields and accumulates a sum and a finite-value
    count, then returns their ratio as a mean field.

    Returns
    -------
    numpy.ndarray
        2-D array of shape (3600, 7200) with mean SST values (°C).
        Pixels with no valid observations are NaN.
    """
    sumarr = np.zeros((3600,7200))
    cntarr = np.zeros((3600,7200))

    dtmvec = pd.date_range("2001-01-01", "2009-12-31", freq="MS")
    for dtm in tqdm(dtmvec):
        ds = ostia.open_dataset(dtm)
        sumarr += ds.sst.values[0,:,:]
        cntarr += np.isfinite(ds.sst.values[0,:,:])
    return sumarr / cntarr

def max_month(year1=2001, year2=2009):
    """Find the per-year maximum monthly SST and the month it occurs.

    For each year in ``range(year1, year2)`` (``year2`` exclusive), iterates
    over monthly OSTIA SST fields and records the highest SST value and the
    corresponding month at every grid cell.

    Returns
    -------
    maxarr : numpy.ndarray
        3-D array of shape (n_years, 3600, 7200) with the maximum monthly SST
        (°C) for each year.  Cells that were never updated are NaN.
    monarr : numpy.ndarray of uint16
        3-D array of shape (n_years, 3600, 7200) with the calendar month
        (1–12) at which the maximum SST was reached.
    """
    years = range(year1, year2)

    maxarr = np.zeros((len(years), 3600, 7200)) - 6
    monarr = np.zeros((len(years), 3600, 7200), np.uint16)

    pos = 0
    for year in tqdm(years):
        dtmvec = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="MS")
        for dtm in tqdm(dtmvec, leave=False):
            ds = ostia.open_dataset(dtm)
            mask = ds.sst.values[0,:,:] > maxarr[pos]
            maxarr[pos][mask] = ds.sst.values[0,:,:][mask]
            monarr[pos][mask] = dtm.month
        pos += 1
    maxarr[maxarr==-6] = np.nan
    return maxarr, monarr

def min_month():
    """Find the per-year minimum monthly SST and the month it occurs.

    For each year in 2001–2008, iterates over monthly OSTIA SST fields and
    records the lowest SST value and the corresponding month at every grid
    cell.

    Returns
    -------
    minarr : numpy.ndarray
        3-D array of shape (n_years, 3600, 7200) with the minimum monthly SST
        (°C) for each year.  Cells that were never updated are NaN.
    monarr : numpy.ndarray of uint16
        3-D array of shape (n_years, 3600, 7200) with the calendar month
        (1–12) at which the minimum SST was reached.
    """
    years = range(2001,2009)

    minarr = np.zeros((len(years), 3600, 7200)) + 50
    monarr = np.zeros((len(years), 3600, 7200), np.uint16)

    pos = 0
    for year in tqdm(years):
        dtmvec = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="MS")
        for dtm in tqdm(dtmvec, leave=False):
            ds = ostia.open_dataset(dtm)
            mask = ds.sst.values[0,:,:] < minarr[pos]
            minarr[pos][mask] = ds.sst.values[0,:,:][mask]
            monarr[pos][mask] = dtm.month
        pos += 1
    minarr[minarr==minarr.max()] = np.nan
    return minarr, monarr


def clim_max_month():
    """Compute the climatological maximum monthly SST and the warmest month.

    Averages OSTIA SST over 2001–2008 for each calendar month, then takes the
    maximum across months at every grid cell to produce a climatological peak
    SST and the corresponding month index.

    Returns
    -------
    maxarr : numpy.ndarray
        2-D array of shape (3600, 7200) with the climatological maximum SST
        (°C) across all calendar months.
    monarr : numpy.ndarray
        2-D array of shape (3600, 7200) with the calendar month (1–12) of the
        climatological maximum.
    """
    years = range(2001,2009)

    sumarr = np.zeros((12, 3600, 7200))
    cntarr = np.zeros((12, 3600, 7200), np.uint16)
    for month in tqdm(range(12)):
        for year in tqdm(years, leave=False):
            ds = ostia.open_dataset(f"{year}-{month+1}-01")
            sumarr[month,:,:] += ds.sst.values[0,:,:]
            cntarr[month,:,:] += np.isfinite(ds.sst.values[0,:,:])
    maxarr = np.max(sumarr / cntarr, axis=0)
    monarr = np.argmax(sumarr / cntarr, axis=0) + 1
    return maxarr, monarr


def regrid(arr):
    """Regrid a global SST field from 0.05° to 2° resolution.

    Uses a bucket-average resampler to aggregate the full-resolution OSTIA
    grid (3600 × 7200) down to a 2°×2° plate carrée grid (90 × 180).

    Parameters
    ----------
    arr : numpy.ndarray
        2-D array of shape (3600, 7200) in the OSTIA 0.05° plate carrée
        projection.

    Returns
    -------
    numpy.ndarray
        2-D array of shape (90, 180) with bucket-averaged SST values.
        Cells with no valid source pixels are NaN.
    """
    grid = grids.setup_plate_carree_grid(-180,180,-90,90,7200,3600)
    grid2 = grids.setup_plate_carree_grid(-180,180,-90,90,180,90)
    resampler = bilinear.NumpyBilinearResampler(grid, grid2, radius_of_influence=50000, epsilon=1, reduce_data=False)

    lons,lats = grid.get_lonlats()
    lons = dask.array.array(lons)
    lats = dask.array.array(lats)
    resampler = BucketResampler(grid2, lons, lats)
    return np.array(resampler.get_average(arr, skipna=True))


def nearest_ostia(return_dist=False):
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

    ds = ostia.open_dataset(timetype="monthly", chunks="auto")
    lons,lats = np.meshgrid(ds.lon,ds.lat)
    latlon2 = np.deg2rad(np.array([lats.flatten(),lons.flatten()]).T)
    dist,ij = tree.query(latlon2)
    if return_dist:
        return dist,ij
    return ij


def ostia_to_2deg(arr, ij=None):
    """Resample an OSTIA field to the 2°-grid and assign values to regions.

    First regrids ``arr`` from 0.05° to 2° resolution using bucket averaging,
    then maps each 2°-grid cell to its corresponding mintmat region index.

    Parameters
    ----------
    arr : numpy.ndarray
        2-D array of shape (3600, 7200) in the OSTIA 0.05° grid.
    ij : numpy.ndarray, optional
        Precomputed flat indices from :func:`nearest`.  Computed on-the-fly if
        not provided.

    Returns
    -------
    numpy.ndarray
        1-D array of SST values (°C) with one entry per mintmat region
        (including a leading NaN for region 0).
    """
    ij = nearest_ostia() if ij is None else ij
    arr = regrid(arr)
    sst = arr.flatten()[ij]
    sst = np.hstack([np.array((np.nan,)), np.squeeze(sst)])
    return sst

def to_reg(arr, ij=None):
    """Average OSTIA full-resolution SST values into mintmat regions.

    Maps each OSTIA pixel to the nearest mintmat region centroid and computes
    the mean SST within each region.

    Parameters
    ----------
    arr : numpy.ndarray
        2-D array of shape (3600, 7200) with SST values (°C) on the full
        OSTIA 0.05° grid.
    ij : numpy.ndarray, optional
        Precomputed region indices from :func:`nearest_ostia`.  Computed
        on-the-fly if not provided.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by region number (starting at 0) with a single
        column ``sst`` containing the mean SST (°C) for each region.
        Region 0 is set to NaN.
    """
    ij = nearest_ostia() if ij is None else ij
    df = pd.DataFrame({"ij":np.squeeze(ij+1), "sst":arr.flatten()})
    svec = df.groupby("ij").mean().reset_index()
    svec = pd.concat([pd.DataFrame({"ij":0, "sst":[np.nan]}),svec], axis=0)
    svec.set_index("ij", inplace=True)
    return svec


def calc_reg_maxmon(year1=2001, year2=2009):
    """Compute regional maximum SST statistics and write them to the mintmat file.

    Calculates four per-region SST summary variables over 2001–2008 and
    appends them to ``mintmat_2001-2009.nc``:

    - ``abs_max_sst``: SST at the absolute warmest month across all years.
    - ``abs_mon_sst``: Calendar month of the absolute warmest SST.
    - ``clm_max_sst``: SST at the climatologically warmest month.
    - ``clm_mon_sst``: Calendar month of the climatological maximum.

    The output file is written with zlib compression (level 9).

    Returns
    -------
    None
    """
    ij = nearest_ostia()

    maxarr, mon = max_month(year1, year2)
    yr_idx = np.argmax(maxarr, axis=0)
    sp_idx = np.indices(yr_idx.shape)
    abs_max_sst = to_reg(maxarr[yr_idx, sp_idx[0], sp_idx[1]], ij)
    abs_mon_sst = to_reg(mon[yr_idx, sp_idx[0], sp_idx[1]], ij)

    climarr,mon = clim_max_month()
    clm_max_sst = to_reg(climarr, ij)
    clm_mon_sst = to_reg(mon, ij)

    def add_reg_sst(name, arr, desc):
        gl[name] = (("to_reg"), arr.sst)
        gl[name].attrs["Description"] = desc

    abs_mon_sst = abs_mon_sst.rename(columns={"sst":"month"})
    return pd.concat([abs_max_sst,abs_mon_sst], axis=1)
