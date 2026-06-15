from pathlib import Path

import numpy as np
import xarray as xr
import pandas as pd


from abstemp import warmest_month, read_ostia_hists, read_cmip6_hists
from abstemp.seagrid import cmip6
from abstemp.data import open_ostia_1985, open_ostia_2019

sstvec = np.arange(-4,40,0.25)

def process(ds: xr.Dataset, clim: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Compute an area-weighted histogram of the warmest monthly SST.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with an ``sst`` variable and spatial dimensions
        ``lat`` / ``lon``.
    clim : bool, optional
        If True, histogram all monthly SST values in ``ds.sst`` (shape
        ``(time, lat, lon)``) and normalise counts by the number of time
        steps.  If False (default), compute the climatological warmest-month
        SST via :func:`abstemp.warmest_month.warmest_monthly_sst` and
        histogram the resulting 2-D field.

    Returns
    -------
    x : numpy.ndarray
        Bin edges (°C), shape ``(n_bins + 1,)``.
    y : numpy.ndarray
        Area-weighted bin counts (km²), shape ``(n_bins,)``.
        Bins with zero area are set to NaN.
    """
    if clim:
        sst = ds.sst.values  # (time, lat, lon)
    else:
        sst = warmest_month.warmest_monthly_sst(ds)  # (lat, lon)
    mask = np.isfinite(sst)
    area = warmest_month.haversine_area(ds)
    area_3d = np.broadcast_to(area, sst.shape)
    y, x = np.histogram(sst[mask], sstvec, weights=area_3d[mask])
    y[y==0] = np.nan
    if clim:
        y = y / sst.shape[0]
    return x, y


def cesm2() -> tuple[np.ndarray, np.ndarray]:
    """Return warmest-month SST histogram for CESM2 SSP5-8.5 (2090–2100).

    Returns
    -------
    x : numpy.ndarray
        Bin edges (°C).
    y : numpy.ndarray
        Area-weighted bin counts (km²).
    """
    ds = cmip6.open_dataset(model="cesm2", experiment="ssp5_8_5")
    return process(ds)


def ostia85() -> tuple[np.ndarray, np.ndarray]:
    """Return warmest-month SST histogram from OSTIA reanalysis 1985–1990.

    Returns
    -------
    x : numpy.ndarray
        Bin edges (°C).
    y : numpy.ndarray
        Area-weighted bin counts (km²).
    """
    ds = open_ostia_1985()
    return process(ds)

def ostia19() -> tuple[np.ndarray, np.ndarray]:
    """Return warmest-month SST histogram from OSTIA NRT product 2019–2023.

    Returns
    -------
    x : numpy.ndarray
        Bin edges (°C).
    y : numpy.ndarray
        Area-weighted bin counts (km²).
    """
    ds = open_ostia_2019()
    return process(ds)

def all_cmip(experiment: str = "ssp5_8_5") -> pd.DataFrame:
    """Compute warmest-month SST histograms for all CMIP6 models.

    Parameters
    ----------
    experiment : str, optional
        SSP experiment identifier (e.g. ``"ssp5_8_5"`` or ``"ssp2_4_5"``).
        Default is ``"ssp5_8_5"``.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by SST bin centre (°C) with one column per
        CMIP6 model, containing area-weighted counts (km²).
    """
    mlist = list(cmip6.cmip6_sst_models)
    ydict = {}
    for model in mlist:
        x,y = process(cmip6.open_dataset(model, experiment=experiment))
        ydict[model] = y
    df = pd.DataFrame(ydict).set_index(sstvec[:-1])
    df.index.name = "SST"
    return df

def save_cmip_hists() -> None:
    """Compute and save CMIP6 SST histograms for both SSP scenarios.

    Writes two CSV files to the package data directory:
    ``all_cmip6_hists_ssp585.csv`` and ``all_cmip6_hists_ssp245.csv``.

    Returns
    -------
    None
    """
    datadir = Path(__file__).parent / "data"
    df = all_cmip(experiment="ssp585")
    df.to_csv(datadir / "all_cmip6_hists_ssp585.csv")
    df = all_cmip(experiment="ssp245")
    df.to_csv(datadir / "all_cmip6_hists_ssp245.csv")

def save_ostia_hists() -> pd.DataFrame:
    """Compute and save OSTIA SST histograms for 1985–1990 and 2019–2023.

    Writes ``all_ostia_hists.csv`` to the package data directory.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by SST bin centre (°C) with columns
        ``"1985-1990"`` and ``"2019-2023"``.
    """
    sst,area85 = process(open_ostia_1985())
    sst,area19 = process(open_ostia_2019())
    df = pd.DataFrame({"1985-1990":area85, "2019-2023":area19}).set_index(sst[:-1])
    df.index.name = "SST"
    datadir = Path(__file__).parent / "data"
    df.to_csv(datadir / "all_ostia_hists.csv")
    return df


def significant_round(x: float | np.ndarray, n_figs: int) -> float | np.ndarray:
    """Round *x* to *n_figs* significant figures.

    Parameters
    ----------
    x : float or array_like
        Value(s) to round.
    n_figs : int
        Number of significant figures to keep.

    Returns
    -------
    float or numpy.ndarray
        Rounded value(s).
    """
    power = 10 ** np.floor(np.log10(np.abs(x).clip(1e-200)))
    return np.round(x / power, n_figs - 1) * power


def stats_table():

    c24 = read_cmip6_hists(experiment="ssp245")
    c24 = c24.fillna(0).cumsum()/c24.fillna(0).sum(axis=0)
    c58 = read_cmip6_hists(experiment="ssp585")
    c58 = c58.fillna(0).cumsum()/c58.fillna(0).sum(axis=0)

    ost = read_ostia_hists()
    ost = ost.fillna(0).cumsum()/ost.fillna(0).sum(axis=0)
    df = ost.copy()
    df["2095-2100 ssp245"] = c24.mean(axis=1)
    df["2095-2100 ssp585"] = c58.mean(axis=1)
    tab = (1-df.loc[[24.75, 29.75, 30.75, 31.75, 33.75, 34.75]]).transpose()
    tab.rename(columns={24.75:">=25°C", 29.75:">=30°C", 30.75:">=31°C", 31.75:">=32°C", 33.75:">=34°C", 34.75:">=35°C"}, inplace=True)
    return significant_round(tab*100, 2).map(lambda x: f"{x:.2g}")
