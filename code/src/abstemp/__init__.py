import os
import pathlib

import numpy as np
import xarray as xr
import pandas as pd

from abstemp.data import (
    open_warmest_2095,
)
from abstemp import reg_calculations

vector_figs = False


def _nc_dir() -> pathlib.Path:
    """Return the directory that holds large NetCDF files.

    Respects the ``ABSTEMP_DATA_DIR`` environment variable so the same code
    works both inside a Code Ocean capsule (``/data``) and in a local
    development environment (package data directory).
    """
    env = os.environ.get("ABSTEMP_DATA_DIR")
    return pathlib.Path(env) if env else pathlib.Path(__file__).parent / "data"

def open_mintmat_ds() -> xr.Dataset:
    """Open the mintmat connectivity/SST dataset bundled with the package.

    Returns
    -------
    xarray.Dataset
        Dataset containing the Dijkstra minimum-travel-time matrix
        (``dijkmintmat``), per-region coordinates (``reglat``, ``reglon``),
        and pre-computed SST statistics (``abs_max_sst``, ``abs_mon_sst``,
        ``clm_max_sst``, ``clm_mon_sst``, ``longhurst_basins``, etc.).
    """
    return xr.open_dataset(
        _nc_dir() / "mintmat_2001-2009.nc",
        decode_timedelta=False,
    )

def read_regdegvel() -> pd.DataFrame:
    """Read the pre-computed per-region degree-velocity parquet file.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by region number with columns including
        ``reglat``, ``reglon``, ``{year}_max_sst``, ``{year}_regdegvel``,
        ``{year}_regdegbas``, ``{year}_regdegfrm`` for years 1985, 2019,
        and 2095.
    """
    return pd.read_parquet(pathlib.Path(__file__).parent / "data/abstemp_reg_degvel.parquet")

def read_cmip6_hists(experiment: str = "ssp585") -> pd.DataFrame:
    """Read pre-computed CMIP6 SST histograms for a given SSP experiment.

    Parameters
    ----------
    experiment : str, optional
        SSP scenario identifier.  One of ``"ssp585"`` (default) or
        ``"ssp245"``.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by SST bin centre (°C) with one column per CMIP6
        model, containing area-weighted histogram counts (km²).
    """
    return pd.read_csv(pathlib.Path(__file__).parent / f"data/all_cmip6_hists_{experiment}.csv", index_col=0)

def read_ostia_hists() -> pd.DataFrame:
    """Read pre-computed OSTIA SST histograms for 1985–1990 and 2019–2023.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by SST bin centre (°C) with columns
        ``"1985-1990"`` and ``"2019-2023"``, containing area-weighted
        histogram counts (km²).
    """
    return pd.read_csv(pathlib.Path(__file__).parent / f"data/all_ostia_hists.csv", index_col=0)

def open_longhurst() -> xr.Dataset:
    """Open the Longhurst biogeochemical province dataset.

    Returns
    -------
    xarray.Dataset
        Dataset on a regular lat/lon grid with variables ``basins``,
        ``regions``, and ``biomes`` (integer province codes).
    """
    return xr.open_dataset(_nc_dir() / "Longhurst_Regions_2007.nc")
