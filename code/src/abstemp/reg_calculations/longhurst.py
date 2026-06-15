import numpy as np
import xarray as xr
import pandas as pd
import dask
from tqdm import tqdm
from pyresample import bilinear
from pyresample.bucket import BucketResampler
from sklearn.neighbors import BallTree

import abstemp

def nearest_longhurst(mask=None, return_dist=False):
    """Find the nearest Longhurst pixel for each mintmat region centroid.

    Builds a haversine BallTree on the non-zero Longhurst grid pixels and
    queries it with each mintmat region centroid to obtain a Longhurst pixel
    assignment for every region.

    Parameters
    ----------
    mask : array_like of bool, optional
        2-D boolean mask selecting which Longhurst pixels to include in the
        tree.  Defaults to pixels where ``regions != 0``.
    return_dist : bool, optional
        If True, also return the haversine distances (in radians) to the
        nearest Longhurst pixel.  Default is False.

    Returns
    -------
    ij : numpy.ndarray
        1-D array of flat Longhurst-pixel indices, one per mintmat region
        (regions 1 … N, skipping index 0).
    dist : numpy.ndarray
        Haversine distances in radians (only returned when
        ``return_dist=True``).
    """
    lh = abstemp.open_longhurst()
    lons,lats = np.meshgrid(lh.lon,lh.lat)
    mask= lh.regions!=0 if mask is None else mask
    latlon = np.deg2rad(np.array([lats[mask], lons[mask]]).T)
    tree = BallTree(latlon, metric="haversine")

    gl = abstemp.open_mintmat_ds()
    latlon2 = np.deg2rad(np.array([gl.reglat[1:],gl.reglon[1:]]).T)
    dist,ij = tree.query(latlon2)
    if return_dist:
        return dist,ij
    return ij


def lh_to_reg(arr, ij=None):
    """Assign a Longhurst categorical value to each mintmat region.

    For each mintmat region centroid, looks up the value of the nearest
    non-zero Longhurst grid pixel using a precomputed nearest-neighbour index.

    Parameters
    ----------
    arr : xarray.DataArray or array_like
        2-D Longhurst field (e.g. ``basins``, ``regions``, or ``biomes``)
        defined on the Longhurst lat/lon grid.  If a DataArray, its
        ``.name`` attribute is used as the output column name.
    ij : numpy.ndarray, optional
        Precomputed pixel indices from :func:`nearest_longhurst`.  Computed
        on-the-fly if not provided.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by region number (0 … N) with one column (named
        after *arr*) containing the Longhurst category for each mintmat region.
        Index 0 is set to NaN (no region).
    """
    if isinstance(arr, xr.DataArray):
        name = arr.name
        arr = arr.data
    else:
        name = "data"

    lh = abstemp.open_longhurst()
    mask= lh.regions!=0
    ij = nearest_longhurst(mask=mask) if ij is None else ij
    vec = np.squeeze(arr[mask][ij].T)
    df = pd.DataFrame({"ij":np.arange(len(vec))+1, name:vec})
    df = pd.concat([pd.DataFrame({"ij":0, name:[np.nan]}),df], axis=0)
    df.set_index("ij", inplace=True)
    return df
    #agg(pd.Series.mode)
