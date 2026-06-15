
import os
import pathlib

import xarray as xr
import pandas as pd

from abstemp import warmest_month
from abstemp.tempvel import generate_regdegvel_df
from abstemp.seagrid import cmip6
from . import download


def _nc_dir() -> pathlib.Path:
    """Return the directory that holds large NetCDF files.

    Respects the ``ABSTEMP_DATA_DIR`` environment variable so the same code
    works both inside a Code Ocean capsule (``/data``) and in a local
    development environment (package data directory).
    """
    env = os.environ.get("ABSTEMP_DATA_DIR")
    return pathlib.Path(env) if env else pathlib.Path(__file__).parent

def max_min_month(ds):
    """Compute per-pixel maximum and minimum monthly SST statistics.

    Parameters
    ----------
    ds : xarray.Dataset
        Monthly SST dataset with a ``sst`` variable and ``time`` dimension.
        Must start in January and end in December.

    Returns
    -------
    xarray.Dataset
        Dataset with ``lat``/``lon`` dimensions containing:
        ``maxarr``, ``maxmon``, ``minarr``, ``minmon`` (absolute extremes),
        and ``clim_maxarr``, ``clim_maxmon``, ``clim_minarr``, ``clim_minmon``
        (climatological extremes).
    """
    if pd.Timestamp(ds.time[0].values).month != 1:
        raise ValueError("First month must be January")
    if pd.Timestamp(ds.time[-1].values).month != 12:
        raise ValueError("Last month must be December")

    maxarr = ds['sst'].fillna(-999).max(dim='time')
    maxmon = (ds['sst'].fillna(-999).argmax(dim='time') % 12) + 1
    minarr = ds['sst'].fillna(-999).min(dim='time')
    minmon = (ds['sst'].fillna(-999).argmin(dim='time') % 12) + 1

    clim = ds['sst'].groupby('time.month').mean('time')
    clim_maxmon = clim.idxmax('month')
    clim_minmon = clim.idxmin('month')
    clim_maxarr = clim.max('month')
    clim_minarr = clim.min('month')

    dds = xr.Dataset()
    dds["maxarr"] = maxarr
    dds["maxmon"] = maxmon
    dds["minarr"] = minarr
    dds["minmon"] = minmon
    dds["clim_maxarr"] = clim_maxarr
    dds["clim_maxmon"] = clim_maxmon
    dds["clim_minarr"] = clim_minarr
    dds["clim_minmon"] = clim_minmon

    return dds

def generate_ostia_maxmonsst_files() -> None:
    """Compute and save max-month SST statistics for OSTIA.

    Reads the raw monthly OSTIA SST files via :func:`open_ostia_1985` and
    :func:`open_ostia_2019`, applies :func:`max_min_month` to each, and writes
    the results to ``ostia_maxmonsst_1985-1990.nc`` and
    ``ostia_maxmonsst_2019-2023.nc`` in the package data directory.

    Requires the raw OSTIA monthly SST files to be present (see
    :func:`abstemp.data.download.ostia_sst_fields`).

    Returns
    -------
    None
    """
    datadir = pathlib.Path(__file__).parent
    ds = open_ostia_1985()
    dds = max_min_month(ds)
    dds.to_netcdf(datadir / "ostia_maxmonsst_1985-1990.nc")
    ds = open_ostia_2019()
    dds = max_min_month(ds)
    dds.to_netcdf(datadir / "ostia_maxmonsst_2019-2023.nc")

def generate_cmip6_maxmonsst_files() -> None:
    """Compute and save max-month SST statistics for all CMIP6 models.

    Iterates over all models in :data:`abstemp.seagrid.cmip6.cmip6_sst_models`
    and both SSP experiments, applies :func:`max_min_month` to each, and
    writes the results to
    ``maxmonsst_cmip6/{model}_{experiment}_maxmonsst.nc`` in the package data
    directory.  The output subdirectory is created if it does not exist.

    Requires the CMIP6 NetCDF files to be present in
    ``src/abstemp/data/cmip6/`` (see
    :func:`abstemp.seagrid.cmip6.retrieve_all_files`).

    Returns
    -------
    None
    """
    datadir = pathlib.Path(__file__).parent / "maxmonsst_cmip6"
    datadir.mkdir(exist_ok=True)
    for model in cmip6.cmip6_sst_models:
        for experiment in cmip6.experiments:
            ds = cmip6.open_dataset(model, experiment)
            try:
                dds = max_min_month(ds)
            except TypeError:
                print(f"{model}_{experiment} failed.")
                continue
            dds.to_netcdf(datadir / f"{model}_{experiment}_maxmonsst.nc")
            print(f"   Generated {model}_{experiment}_maxmonsst.nc")


def open_ostia_1985() -> xr.Dataset:
    """Open the raw monthly OSTIA reanalysis SST dataset for 1985–1990.

    Returns
    -------
    xarray.Dataset
        Monthly SST fields with ``time``, ``lat``, ``lon`` dimensions.
    """
    return xr.open_dataset(_nc_dir() / "ostia_sst_1985-1990.nc")

def open_ostia_2019() -> xr.Dataset:
    """Open the raw monthly OSTIA NRT SST dataset for 2019–2023.

    Returns
    -------
    xarray.Dataset
        Monthly SST fields with ``time``, ``lat``, ``lon`` dimensions.
    """
    return xr.open_dataset(_nc_dir() / "ostia_sst_2019-2023.nc")

def open_cmip6_2095(model="ecearth") -> xr.Dataset:
    """Open the raw monthly CMIP6 SST dataset for 2095–2100.

    Parameters
    ----------
    model : str, optional
        Model name prefix used to construct the filename
        ``{model}_sst_2095-2100.nc``.  Default is ``"ecearth"``.

    Returns
    -------
    xarray.Dataset
        Monthly SST fields with ``time``, ``lat``, ``lon`` dimensions.
    """
    return xr.open_dataset(_nc_dir() / f"{model}_sst_2095-2100.nc")


def open_warmest_1985() -> xr.Dataset:
    """Open the pre-computed OSTIA max-month SST statistics for 1985–1990.

    Returns
    -------
    xarray.Dataset
        Per-pixel SST statistics on a ``lat``/``lon`` grid; see
        :func:`max_min_month` for variable descriptions.
    """
    return xr.open_dataset(_nc_dir() / "ostia_maxmonsst_1985-1990.nc")

def open_warmest_2019() -> xr.Dataset:
    """Open the pre-computed OSTIA NRT max-month SST statistics for 2019–2023.

    Returns
    -------
    xarray.Dataset
        Per-pixel SST statistics on a ``lat``/``lon`` grid; see
        :func:`max_min_month` for variable descriptions.
    """
    return xr.open_dataset(_nc_dir() / "ostia_maxmonsst_2019-2023.nc")

def open_warmest_2095(model="cnrm_cm6_1_hr", experiment="ssp5_8_5") -> xr.Dataset:
    """Open pre-computed CMIP6 max-month SST statistics for 2095–2100.

    Parameters
    ----------
    model : str, optional
        Model key matching the filename stem in ``maxmonsst_cmip6/``.
        Default is ``"cnrm_cm6_1_hr"``.
    experiment : str, optional
        SSP experiment identifier.  Default is ``"ssp5_8_5"``.

    Returns
    -------
    xarray.Dataset
        Per-pixel SST statistics on a ``lat``/``lon`` grid, centred on the
        Greenwich meridian; see :func:`max_min_month` for variable
        descriptions.
    """
    datadir = _nc_dir() / "maxmonsst_cmip6"
    fn = f"{model}_{experiment}_maxmonsst.nc"
    return cmip6.center_on_gmt(xr.open_dataset(datadir / fn))

def open_longhurst() -> xr.Dataset:
    """Open the Longhurst province dataset bundled with the package.

    Returns
    -------
    xarray.Dataset
        Dataset on a regular lat/lon grid with variables ``basins``,
        ``regions``, and ``biomes``.
    """
    return xr.open_dataset(_nc_dir() / "Longhurst_Regions_2007.nc")


def all() -> None:
    """Regenerate all derived data files in the package data directory.

    Requires OSTIA access (via njord) for the 1985 and 2019 SST files and
    the CMIP6 monthly NetCDF files to be present (see
    :func:`abstemp.seagrid.cmip6.retrieve_all_files`).

    Returns
    -------
    None
    """
    _datadir = pathlib.Path(__file__).parent
    warmest_month.save_ostia_files()
    generate_ostia_maxmonsst_files()
    generate_cmip6_maxmonsst_files()
    df = generate_regdegvel_df()
    df.to_parquet(_datadir / "abstemp_reg_degvel.parquet")

def setup() -> None:
    """Download all files required for Use Case 1 (bundled-data workflow).

    Fetches the large SST NetCDF files that are not tracked in git
    (``ostia_maxmonsst_1985-1990.nc``, ``ostia_maxmonsst_2019-2023.nc``,
    and CNRM-CM6-1-HR max-month files in ``maxmonsst_cmip6/``) and the
    Longhurst province mask, plus the mintmat connectivity file from Zenodo.
    All other files (parquet, CSVs) are already present in the repository.

    Returns
    -------
    None
    """
    download.maxmonsst_fields()
    download.longhurst_regions()
    download.mintmat()
