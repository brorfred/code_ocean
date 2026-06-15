
import pathlib

import numpy as np
import xarray as xr

from abstemp.seagrid import gridtools

from .sst_ostia import (
    nearest_ostia,
    to_reg,
    max_month,
    min_month,
    clim_max_month
)
import abstemp
from . import longhurst

_DATADIR = pathlib.Path(__file__).parent.parent / "data"

lon = np.linspace(-179,179,180)
lat = np.linspace(-89,89,90)
lons, lats = np.meshgrid(lon, lat)

def add_maxmon():
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

    maxarr, mon = max_month()
    yr_idx = np.argmax(maxarr, axis=0)
    sp_idx = np.indices(yr_idx.shape)
    abs_max_sst = to_reg(maxarr[yr_idx, sp_idx[0], sp_idx[1]], ij)
    abs_max_mon = to_reg(mon[yr_idx, sp_idx[0], sp_idx[1]], ij)

    minarr, mon = min_month()
    yr_idx = np.argmin(minarr, axis=0)
    abs_min_sst = to_reg(minarr[yr_idx, sp_idx[0], sp_idx[1]], ij)
    abs_min_mon = to_reg(mon[yr_idx, sp_idx[0], sp_idx[1]], ij)

    climarr, mon = clim_max_month()
    clm_max_sst = to_reg(climarr, ij)
    clm_max_mon = to_reg(mon, ij)

    def add_reg_sst(name, arr, desc):
        gl[name] = (("to_reg"), arr.sst)
        gl[name].attrs["Description"] = desc

    arr_list = [abs_max_sst, abs_max_mon, abs_min_sst, abs_min_mon, clm_max_sst, clm_max_mon]
    name_list = ["abs_max_sst", "abs_max_mon", "abs_min_sst", "abs_min_mon", "clm_max_sst", "clm_max_mon"]
    desc_list = [
        "Monthly SST (°C) at the warmest month 2001-2009",
        "Warmest month 2001-2009",
        "Monthly SST (°C) at the coldest month 2001-2009",
        "Coldest month 2001-2009",
        "Monthly SST (°C) at the warmest climatological month 2001-2009",
        "Warmest climatological month 2001-2009",
    ]

    gl = abstemp.open_mintmat_ds()
    gl["nreg"] = (("to_reg"), np.arange(1,11117))
    for arr, name, desc in zip(arr_list, name_list, desc_list):
        add_reg_sst(name, arr, desc)
    encoding = {var: {"zlib": True, "complevel": 9} for var in gl.data_vars}
    gl.to_netcdf(_DATADIR / "mintmat_2001-2009.nc", encoding=encoding)


def add_longhurst():

    gl = abstemp.open_mintmat_ds()
    lh = longhurst.lh

    arr = lh.basins
    arr = arr.where(arr!=0)
    rg = longhurst.lh_to_reg(arr)
    gl["longhurst_basins"] = (("to_reg"), rg[arr.name])

    arr = lh.regions
    arr = arr.where(arr!=0)
    rg = longhurst.lh_to_reg(arr)
    gl["longhurst_regions"] = (("to_reg"), rg[arr.name])

    arr = lh.biomes
    arr = arr.where(arr!=0)
    rg = longhurst.lh_to_reg(arr)
    gl["longhurst_biomes"] = (("to_reg"), rg[arr.name])
    encoding = {var: {"zlib": True, "complevel": 9} for var in gl.data_vars}
    gl.to_netcdf(_DATADIR / "mintmat_2001-2009.nc", encoding=encoding)


def get_mode(resampler, arr, fill_value=np.nan):
    """Return the most common value in each bucket of a BucketResampler.

    Analogous to ``BucketResampler.get_average``, but returns the mode
    (most frequently occurring value) instead of the mean.  Uses
    ``get_fractions`` internally, so *arr* should contain categorical
    (integer) values.  Ties are broken by the last tied category encountered
    in sorted order (i.e. the largest tied value wins).  Empty buckets are
    set to *fill_value*.

    Parameters
    ----------
    resampler : pyresample.bucket.BucketResampler
        An already-initialised resampler.
    arr : numpy.ndarray or dask array
        Categorical source data to be binned.
    fill_value : float, optional
        Value used for empty buckets (no valid source pixels).
        Default: ``np.nan``.

    Returns
    -------
    numpy.ndarray
        2-D array with the same shape as ``resampler.target_area.shape``
        where each cell holds the modal value of the source pixels that
        fell into it, or *fill_value* for empty buckets.
    """
    out_shape = resampler.target_area.shape
    # Use fill_value=0 so empty buckets get fraction 0 and are never selected
    fractions = resampler.get_fractions(arr, categories=np.unique(arr), fill_value=0.0)

    mode = np.full(out_shape, fill_value)
    best_frac = np.zeros(out_shape)
    for cat, frac_arr in fractions.items():
        frac_arr = np.array(frac_arr)
        update = frac_arr > best_frac
        mode = np.where(update, cat, mode)
        best_frac = np.maximum(best_frac, frac_arr)
    return mode

def nearest(return_dist=False):
    """Find the nearest 2°-grid cell index for each mintmat region centroid.

    Builds a haversine BallTree on the 2°×2° plate carrée grid and queries it
    with the region centroids from ``data/mintmat_2001-2009.nc`` to map every
    region to the closest 2°-grid flat index.

    Parameters
    ----------
    return_dist : bool, optional
        If True, also return the haversine distances (in radians) to the
        nearest grid cell.  Default is False.

    Returns
    -------
    ij : numpy.ndarray
        1-D array of flat indices into the 2°-grid (shape 90 × 180) for each
        region centroid.
    dist : numpy.ndarray
        Haversine distances in radians (only returned when
        ``return_dist=True``).
    """
    lat = np.linspace(-89,89,90)
    lon = np.linspace(-179,179,180)
    lons,lats = np.meshgrid(lon,lat)
    ds = abstemp.open_mintmat_ds()
    from sklearn.neighbors import BallTree
    latlon = np.deg2rad(np.array([lats.flatten(),lons.flatten()]).T)
    tree = BallTree(latlon, metric="haversine")
    latlon2 = np.deg2rad(np.array([ds.reglat[1:],ds.reglon[1:]]).T)
    dist,ij = tree.query(latlon2)
    if return_dist:
        return dist,ij
    return ij


def nreg_arr(ij=None):
    """Build a 2-D array mapping each 2°-grid cell to its mintmat region index.

    Parameters
    ----------
    ij : numpy.ndarray, optional
        Precomputed flat indices from :func:`nearest`.  Computed on-the-fly if
        not provided.

    Returns
    -------
    numpy.ndarray
        2-D array of shape (90, 180) where each element holds the mintmat
        region index (1-based) of the corresponding 2°-grid cell.  Cells not
        assigned to any region retain the value 0.
    """
    ij = nearest() if ij is None else ij
    ijarr = np.zeros([90,180])
    ijarr.flat[ij.flatten()] = np.arange(1,11116).astype(int)

    return ijarr

def arr_to_regvec(arr):
    """Map a 2°×2° spatial array onto a per-region vector.

    Inverse of :func:`regvec_to_arr`.  For each mintmat region, picks the
    value of the 2°-grid cell whose centre is nearest to the region centroid.

    Parameters
    ----------
    arr : array-like, shape (90, 180)
        2-D array on the 2°×2° plate carrée grid.

    Returns
    -------
    numpy.ndarray, shape (11116,)
        1-D array indexed by mintmat region number (0-based).  Index 0 is
        ``NaN`` (no region); indices 1–11115 hold the value of the
        corresponding grid cell.
    """
    ij = nearest()
    values = np.asarray(arr).flatten()[ij.flatten()]
    return np.concatenate([[np.nan], values])


def regvec_to_arr(da):
    """Map a per-region DataArray back onto the 2°×2° spatial grid.

    Parameters
    ----------
    da : xarray.DataArray
        1-D DataArray indexed by mintmat region number (0-based), with one
        value per region.

    Returns
    -------
    numpy.ndarray
        2-D array of shape (90, 180) where each cell holds the value of the
        mintmat region it belongs to.
    """
    ijarr = nreg_arr()
    return da[ijarr.flatten().astype(int)].values.reshape(ijarr.shape)


def latlon_to_reg(lat_val, lon_val):
    """Return the mintmat region index for the 2°-grid cell nearest to (lat, lon).

    Parameters
    ----------
    lat_val : float
        Latitude in degrees (−90 to 90).
    lon_val : float
        Longitude in degrees (−180 to 180).

    Returns
    -------
    int
        Mintmat region index (1-based).  Returns 0 for cells with no assigned
        region.
    """
    i_lat = np.argmin(np.abs(lat - lat_val))
    i_lon = np.argmin(np.abs(lon - lon_val))
    return int(nreg_arr()[i_lat, i_lon])


def reg_area() -> np.ndarray:
    """Compute the area (km²) of each 2°×2° grid cell on the plate carrée grid.

    Returns
    -------
    numpy.ndarray
        2-D array of shape (90, 180) with grid-cell areas in km².
    """
    return gridtools.area(lats=lats, lons=lons)
