"""Maps of maximum monthly sea-surface temperature across climate periods.

Functions here produce global maps of the warmest monthly SST, either for a
single time slice or in a three-panel comparison across 1985–1990,
2019–2024, and 2095–2100.  CMIP6 model maps can also be generated in bulk.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import xarray as xr
import pandas as pd

import cmap
import projmap

import abstemp
from abstemp import vector_figs
from abstemp import (
    data,
    warmest_month,
    global_sst_histograms,
    tempvel,
)
from abstemp.seagrid import cmip6
from .figpref import title_kws,lat,lon

def warmest_map(ds=None, mp=None, ax=None, colorbar=True, title="2019$-$2023"):
    """
    Plot a categorical map of the maximum monthly SST.

    SST values are binned into discrete categories at fixed thresholds
    (25, 30–37 °C) and rendered as a filled-contour map using a sequential
    colour scale.

    Parameters
    ----------
    ds : xarray.Dataset, optional
        Dataset containing an ``sst`` variable with a ``time`` dimension.
        If None, the 2019–2023 NRT OSTIA dataset is loaded automatically.
    mp : projmap.Map, optional
        Map projection object to draw on. If None, a new global ``"glob"``
        map is created after clearing the current figure.
    ax : int or matplotlib.axes.Axes, optional
        Axes index (for multi-panel maps) or Axes object passed to
        ``mp.contourf`` and ``mp.nice``. Default is None.
    colorbar : bool, optional
        Whether to draw a colourbar with labelled SST thresholds.
        Default is True.
    title : str, optional
        Annotation text placed in the upper-right of the map.
        Default is ``"2019$-$2023"``.

    Returns
    -------
    None
    """


    if ds is None:
        ds = data.open_warmest_2019()
    da = ds.maxarr if "maxarr" in ds else ds.sst.max(dim="time")
    cat = np.zeros(da.shape)
    levels = [25, 30, 31, 32, 33, 34, 35, 36, 37]
    for lev in levels:
        cat[da.values>=lev] += 1.0
    cat = np.where(np.isfinite(da.values), cat, -1)
    colors = ["#6baed6", "#addd8e", "#ffffb2", "#fed976",
              "#feb24c", "#fd8d3c", "#f03b20", "#bd0026",
              "darkviolet", "indigo"]
    if mp is None:
        plt.clf()
        mp = projmap.Map("glob")
    if "lons" in ds:
        lons = ds.lons.values
        lats = ds.lats.values
    else:
        lons = ds.lon.values
        lats = ds.lat.values
    mp.contourf(lons, lats, cat, np.arange(-1,10)+0.5, colors=colors, colorbar=colorbar, zorder=0, ax=ax)

    if colorbar:
        mp._cb.set_ticks(np.arange(0.5,9.5,1))
        mp._cb.set_ticklabels(["25", "30", "31", "32", "33", "34", "35", "36", "37"])
        mp._cb.ax.tick_params(labelsize=8)
        mp._cb.ax.set_xlabel("°C")
    mp.nice(borders=False, ax=ax)
    mp.text(150, 60, title, **title_kws)

def warmest_map_pcolor(ds=None, mp=None, ax=None, colorbar=True, clim=False, title="2019$-$2023"):
    """
    Plot a categorical map of the maximum monthly SST using pcolor.

    Equivalent to ``warmest_map`` but renders with ``pcolor`` rather than
    ``contourf``. SST values are binned into discrete categories at fixed
    thresholds (25, 30–37 °C) and a ``ListedColormap`` with ``BoundaryNorm``
    is used to map each integer category to the same sequential colour scale.

    Parameters
    ----------
    ds : xarray.Dataset, optional
        Dataset containing an ``sst`` variable with a ``time`` dimension.
        If None, the 2019–2023 NRT OSTIA dataset is loaded automatically.
    mp : projmap.Map, optional
        Map projection object to draw on. If None, a new global ``"glob"``
        map is created after clearing the current figure.
    ax : int or matplotlib.axes.Axes, optional
        Axes index (for multi-panel maps) or Axes object passed to
        ``mp.pcolor`` and ``mp.nice``. Default is None.
    colorbar : bool, optional
        Whether to draw a colourbar with labelled SST thresholds.
        Default is True.
    title : str, optional
        Annotation text placed in the upper-right of the map.
        Default is ``"2019$-$2023"``.

    Returns
    -------
    None

    See Also
    --------
    warmest_map : Equivalent function using ``contourf``.
    """
    if ds is None:
        ds = data.open_warmest_2019()
    if clim:
        da = ds.clim_maxarr if "clim_maxarr" in ds else ds.sst.max(dim="time")
    else:
        da = ds.maxarr if "maxarr" in ds else ds.sst.max(dim="time")
    cat = np.zeros(da.shape)
    levels = [25, 30, 31, 32, 33, 34, 35, 36, 37]
    for lev in levels:
        cat[da.values >= lev] += 1.0
    cat = np.ma.masked_where(~np.isfinite(da.values), cat)

    colors = ["#6baed6", "#addd8e", "#ffffb2", "#fed976",
              "#feb24c", "#fd8d3c", "#f03b20", "#bd0026",
              "darkviolet", "indigo"]
    cmap_cat = ListedColormap(colors)
    bounds = np.arange(-0.5, 10.5, 1)
    norm = BoundaryNorm(bounds, cmap_cat.N)

    if mp is None:
        plt.clf()
        mp = projmap.Map("glob")
    if "lons" in ds:
        lons = ds.lons.values
        lats = ds.lats.values
    else:
        lons = ds.lon.values
        lats = ds.lat.values

    pc = mp.pcolor(lons, lats, cat, cmap=cmap_cat, norm=norm, colorbar=colorbar, zorder=0, ax=ax)

    if colorbar:
        mp._cb.set_ticks(np.arange(0.5, 9.5, 1))
        mp._cb.set_ticklabels(["25", "30", "31", "32", "33", "34", "35", "36", "37"])
        mp._cb.ax.tick_params(labelsize=8)
        mp._cb.ax.set_xlabel("°C")
    mp.nice(borders=False, ax=ax)
    mp.text(150, 60, title, **title_kws)


def warmest_months():
    """
    Produce a three-panel map comparing maximum monthly SST across time periods.

    Draws panels for 1985–1990, 2019–2024, and 2095–2100 using
    ``warmest_map``, then saves the figure to
    ``figs/warmest_sst_map.pdf`` and ``figs/warmest_sst_map.png``
    at 600 dpi.

    Returns
    -------
    None
    """
    plt.clf()
    mp = projmap.Map()
    mp.subplots(3,1)
    for i, (year, cb) in enumerate([(1985, None), (2019, None), (2095, True)]):
        ds = getattr(data, f"open_warmest_{year}")()
        warmest_map_pcolor(ds, mp, i, colorbar=cb, title=f"{year}$-${year+5}")
    if vector_figs:
        plt.savefig("figs/warmest_sst_map.pdf", dpi=600, bbox_inches="tight")
    plt.savefig("figs/warmest_sst_map.png", dpi=600, bbox_inches="tight")


def all_model_maps() -> None:
    """Generate maximum monthly SST maps for every available CMIP6 model.

    Iterates over all models listed in :data:`abstemp.seagrid.cmip6.cmip6_sst_models`
    and both ``ssp5_8_5`` and ``ssp2_4_5`` experiments.  Each map is saved as
    a 600 dpi PNG under ``figs/model_maps/max/<experiment>/figs_<model>_<experiment>.png``.
    Output subdirectories are created automatically.

    Returns
    -------
    None
    """
    datadir = Path("figs/model_maps")
    datadir.mkdir(exist_ok=True)
    (datadir / "clm" / "ssp5_8_5").mkdir(exist_ok=True, parents=True)
    (datadir / "clm" / "ssp2_4_5").mkdir(exist_ok=True, parents=True)
    (datadir / "max" / "ssp5_8_5").mkdir(exist_ok=True, parents=True)
    (datadir / "max" / "ssp2_4_5").mkdir(exist_ok=True, parents=True)


    plt.clf()
    for model in cmip6.cmip6_sst_models:
        for experiment in ["ssp5_8_5", "ssp2_4_5"]:
            ds = data.open_warmest_2095(model=model, experiment=experiment)
            warmest_map_pcolor(ds, title=f"{model} ({experiment})")
            plt.savefig(datadir / "max" / experiment / f"figs_{model}_{experiment}.png", dpi=600, bbox_inches="tight")
            warmest_map_pcolor(ds, title=f"{model} ({experiment})", clim=True)
            plt.savefig(datadir / "clm" / experiment / f"figs_{model}_{experiment}.png", dpi=600, bbox_inches="tight")
            plt.close("all")
            print(model, experiment)
