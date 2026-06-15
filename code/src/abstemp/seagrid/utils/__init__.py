"""Utility functions for seagrid data sources.

Provides helpers for date ranges, directory management, netCDF validation,
and spatial subsetting.
"""

import os, pathlib
import warnings


import xarray as xr

#from .config import settings

import pandas as pd

from .lldist import ll2dist2vec

def dtm_vec(settings):
    """Generate a pandas date range from a settings dictionary.

    Parameters
    ----------
    settings : dict
        Dictionary containing the keys ``dtm_min``, ``dtm_max``, and
        ``dtm_freq`` which are forwarded to ``pandas.date_range``.

    Returns
    -------
    pandas.DatetimeIndex
        Date range spanning from ``dtm_min`` to ``dtm_max`` at the
        frequency given by ``dtm_freq``.
    """
    return pd.date_range(settings["dtm_min"], settings["dtm_max"],
                         freq=settings["dtm_freq"])

def datadir(settings):
    """Create and return the data directory path from a settings dictionary.

    The directory (and any missing parents) is created if it does not
    already exist.

    Parameters
    ----------
    settings : dict
        Dictionary containing the key ``datadir`` with the desired
        directory path.

    Returns
    -------
    str
        The data directory path stored in ``settings["datadir"]``.
    """
    pathlib.Path(settings["datadir"]).mkdir(parents=True, exist_ok=True)
    return settings["datadir"]

def is_nc_file(full_filename):
    """Check if a file exists and is a valid netCDF file.

    Attempts to open the file with ``xarray.open_dataset``.  If the file
    does not exist or cannot be opened, ``False`` is returned.

    Parameters
    ----------
    full_filename : str
        Absolute or relative path to the file to check.

    Returns
    -------
    bool
        ``True`` if the file exists and can be opened as a netCDF
        dataset, ``False`` otherwise.
    """
    if not os.path.isfile(full_filename):
        return False
    try:
        with xr.open_dataset(full_filename) as ds: pass
    except OSError as err:
        #os.unlink(full_filename)
        return False
    return True

def has_nc_data_var(full_filename, data_var):
    """Check if a netCDF file contains a specific data variable.

    The file is opened with ``xarray.open_dataset`` and the presence of
    *data_var* is verified.

    Parameters
    ----------
    full_filename : str
        Path to the netCDF file.
    data_var : str
        Name of the variable to look for in the dataset.

    Returns
    -------
    bool
        ``True`` if the file exists, is a valid netCDF file, and
        contains *data_var*; ``False`` otherwise.
    """
    if not os.path.exists(full_filename):
        return False
    try:
        warnings.simplefilter("ignore")
        ds = xr.open_dataset(full_filename)
    except (OSError, ValueError) as err:
        #os.unlink(full_filename)
        return False
    try:
        ds[data_var]
    except KeyError as err:
        return False
    return True



def vprint(text, verbose=False):
    """Conditionally print text based on a verbose flag.

    Prints *text* when either ``settings["verbose"]`` is truthy or the
    *verbose* parameter is ``True``.

    Parameters
    ----------
    text : str
        The text to print.
    verbose : bool, optional
        If ``True``, force printing regardless of the global settings.
        Default is ``False``.
    """
    if settings["verbose"] or verbose:
        print(text)

def slice_dataframe(ds, settings, latname="lat", lonname="lon"):
    """Subset an xarray Dataset to a lat/lon bounding box from settings.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to subset.
    settings : dict
        Dictionary with keys ``latmin``, ``latmax``, ``lonmin``, and
        ``lonmax`` defining the bounding box in degrees.
    latname : str, optional
        Name of the latitude coordinate in *ds*. Default is ``"lat"``.
    lonname : str, optional
        Name of the longitude coordinate in *ds*. Default is ``"lon"``.

    Returns
    -------
    xarray.Dataset
        The subset of *ds* whose coordinates fall within the bounding
        box defined by *settings*.
    """
    la = ds.coords[latname]
    lc = ds.coords[lonname]
    return ds.loc[{lonname:lc[(lc >= settings["lonmin"]) & 
                              (lc <= settings["lonmax"])], 
                   latname:la[(la >= settings["latmin"]) & 
                              (la <= settings["latmax"])]}]
