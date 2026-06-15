"""Maps and histograms of regional temperature velocity.

Functions here visualise how fast SST is shifting in degrees per year
(temperature velocity) across the global ocean for three climate periods:
1985–1990, 2019–2024, and 2095–2100.
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd

import cmap
import projmap

import abstemp
from abstemp import tempvel
from abstemp.reg_calculations import regvec_to_arr
from .figpref import title_kws,lat,lon

#regvelcm = cmap.Colormap("colorcet:CET_L18").to_mpl()
_reds9 = cmap.Colormap("colorbrewer:YlOrBr_9").to_mpl()
_colors = np.vstack([[1, 1, 1, 1], _reds9(np.linspace(0, 1, 9))])
regvelcm = plt.matplotlib.colors.ListedColormap(_colors)



def tempvel_map(df: "pd.DataFrame", year: int = 1985) -> None:
    """Plot a single-panel temperature-velocity map for one climate period.

    Converts the regional-degree-velocity column to a 2-D array and renders
    it as a global pcolor map (units: °C yr⁻¹).  Longhurst-basin outliers
    (basins 2/3/4) are overlaid as cyan markers, and the highest-SST grid
    cells are shown as scatter points.

    Saves the figure to ``figs/regdegvel_map.pdf`` and
    ``figs/regdegvel_map.png`` at 600 dpi.

    Parameters
    ----------
    df : pandas.DataFrame
        Regional data frame as returned by :func:`abstemp.read_regdegvel`,
        containing ``{year}_regdegvel``, ``{year}_regdegbas``,
        ``{year}_max_sst``, ``longhurst_basins``, ``reglon``, and
        ``reglat`` columns.
    year : int, optional
        Start year of the climate period to plot. Default is 1985.

    Returns
    -------
    None
    """
    pltkw = dict(marker='x', c="cyan", linestyle='none', markeredgewidth=0.5, markersize=2)
    vel = df[f"{year}_regdegvel"]
    bas = df[f"{year}_regdegbas"]
    sst = df[f"{year}_max_sst"]
    darr = abstemp.regvec_to_arr(vel)

    plt.clf()
    mp = projmap.Map()
    mp.pcolor(lon, lat, darr/365, vmax=5, cmap=regvelcm, colorbar=True, rasterized=True)
    mask = (df.longhurst_basins==2) & (bas==3)
    mp.plot(df.reglon[mask], df.reglat[mask], **pltkw)
    mask = (df.longhurst_basins==4) & (bas==3)
    mp.plot(df.reglon[mask], df.reglat[mask], **pltkw)

    pltkw["marker"] = '+'
    mask = (df.longhurst_basins==2) & (bas==4)
    mp.plot(df.reglon[mask], df.reglat[mask], **pltkw)
    mask = (df.longhurst_basins==3) & (bas==4)
    mp.plot(df.reglon[mask], df.reglat[mask], **pltkw)

    maxsst = sst.max()
    hivec = sst>=maxsst-2
    mp.scatter(df.reglon[hivec], df.reglat[hivec], 10)
    mp._cb.ax.set_xlabel("Time to travel 1°C (years)")
    mp.nice()
    plt.savefig("figs/regdegvel_map.pdf", dpi=600, bbox_inches="tight")
    plt.savefig("figs/regdegvel_map.png", dpi=600, bbox_inches="tight")


def tempvel_maps(df=None) -> None:
    """Produce a three-panel temperature-velocity comparison map.

    Reads the regional-degree-velocity data frame and draws panels for
    1985–1990, 2019–2024, and 2095–2100 on a shared map projection.
    The colour bar (°C yr⁻¹) is shown only on the bottom panel.
    Longhurst-basin outliers are overlaid as green dots.

    Saves the figure to ``figs/regdegvel_maps.pdf`` and
    ``figs/regdegvel_maps.png`` at 600 dpi.

    Returns
    -------
    None
    """
    if df is None:
        df = abstemp.read_regdegvel()

    def map(year: int = 1985, mp: "projmap.Map | None" = None,
            ax: "int | None" = None, colorbar: "bool | None" = True,
            title: "str | None" = None) -> None:
        """Draw one temperature-velocity panel onto *mp*.

        Parameters
        ----------
        year : int, optional
            Start year of the climate period. Default is 1985.
        mp : projmap.Map or None, optional
            Shared map projection object. If None, a new one is created.
        ax : int or None, optional
            Sub-panel index for ``mp``. Default is None.
        colorbar : bool or None, optional
            Whether to draw the colour bar. Default is True.
        title : str or None, optional
            Annotation text placed in the upper-right of the panel.

        Returns
        -------
        None
        """
        if mp is None:
            mp = projmap.Map()
        pltkw = dict(marker='o', c="limegreen", linestyle='none', markeredgewidth=0.2, markersize=0.2, ax=ax)
        vel = df[f"{year}_regdegvel"]
        bas = df[f"{year}_regdegbas"]
        sst = df[f"{year}_max_sst"]
        darr = regvec_to_arr(vel)

        #plt.clf()
        #mp = projmap.Map()
        mp.pcolor(lon, lat, darr/365, vmax=5, cmap=regvelcm, colorbar=colorbar, rasterized=True, ax=ax)
        mask = (df.longhurst_basins==2) & (bas==3)
        mp.plot(df.reglon[mask], df.reglat[mask], **pltkw)
        mask = (df.longhurst_basins==4) & (bas==3)
        mp.plot(df.reglon[mask], df.reglat[mask], **pltkw)

        pltkw["marker"] = '+'
        mask = (df.longhurst_basins==2) & (bas==4)
        #mp.plot(df.reglon[mask], df.reglat[mask], **pltkw)
        mask = (df.longhurst_basins==3) & (bas==4)
        mp.plot(df.reglon[mask], df.reglat[mask], **pltkw)

        maxsst = sst.max()
        hivec = sst>=maxsst-1.7
        mp.scatter(df.reglon[hivec], df.reglat[hivec], 5, ax=ax)
        if colorbar:
            mp._cb.ax.set_xlabel("Time to travel 1°C (years)")
        mp.nice(ax=ax)
        mp.text(140, 52, title, **title_kws)


    plt.clf()
    mp = projmap.Map()
    mp.subplots(3, 1)
    map(1985, mp=mp, ax=0, colorbar=None, title="1985-1990")
    map(2019, mp=mp, ax=1, colorbar=None, title="2019-2024")
    map(2095, mp=mp, ax=2, colorbar=True, title="2095-2100")
    plt.savefig("figs/regdegvel_maps.pdf", dpi=600, bbox_inches="tight")
    plt.savefig("figs/regdegvel_maps.png", dpi=600, bbox_inches="tight")

def tempvel_histogram() -> None:
    """Generate and save the temperature-velocity histogram figure.

    Delegates to :func:`abstemp.tempvel.regdegvel_hist` and saves the result
    to ``figs/regdegvel_hist.pdf`` and ``figs/regdegvel_hist.png`` at 600 dpi.

    Returns
    -------
    None
    """
    tempvel.regdegvel_hist()
    plt.savefig("figs/regdegvel_hist.pdf", dpi=600, bbox_inches="tight")
    plt.savefig("figs/regdegvel_hist.png", dpi=600, bbox_inches="tight")



