"""Figures illustrating the methods used in the abstemp analysis.

Includes a checkerboard grid diagnostic and multi-panel maps of the OSTIA
maximum, minimum, and seasonal-range SST fields.
"""

import projmap
import numpy as np
import matplotlib.pyplot as plt

import abstemp
from abstemp import data, vector_figs
from abstemp.figure_scripts.figpref import lon, lat
from abstemp.reg_calculations import regvec_to_arr

def checkerboard() -> None:
    """Plot and save a 2°×2° checkerboard grid as a map diagnostic.

    Creates a binary checkerboard pattern on the module-level ``lon``/``lat``
    grid and renders it as a global pcolor map.  Saves the result to
    ``figs/checkerboard.pdf`` and ``figs/checkerboard.png`` at 600 dpi.

    Returns
    -------
    None
    """
    lons,lats = np.meshgrid(lon,lat)
    check = np.zeros(lons.shape)
    check[:-1:2,:-1:2] = 1
    check[1::2,1::2] = 1
    plt.clf()
    mp = projmap.Map()
    mp.pcolor(lons, lats, check, cmap="gray", vmin=-1)
    mp.nice()
    if vector_figs:
        plt.savefig("figs/checkerboard.pdf", dpi=600, bbox_inches="tight")
    plt.savefig("figs/checkerboard.png", dpi=600, bbox_inches="tight")



def sst_maps() -> None:
    """Plot and save three-panel OSTIA SST maps (max, min, and range).

    Reads the mintmat dataset via :func:`abstemp.open_mintmat_ds`, extracts
    absolute maximum and minimum SST fields, and produces a three-panel
    global map: maximum SST, minimum SST (shared ``Spectral_r`` scale,
    −5 to 35 °C), and the seasonal range (``Reds`` scale).  Colour bars for
    each scale are placed to the right of the panel array.

    Saves the figure to ``figs/ostia_sst.pdf`` and ``figs/ostia_sst.png``
    at 600 dpi.

    Returns
    -------
    None
    """
    ds = data.open_warmest_2019()

    maxsst = ds.maxmon.values
    minsst = ds.minmon.values
    lon = ds.lon.values
    lat = ds.lat.values

    cmap_name = "Spectral_r"
    plt.clf()
    mp = projmap.Map()
    mp.subplots(3, 1)
    mp.fig.subplots_adjust(right=0.85)

    mp.pcolor(lon, lat, maxsst, cmap=cmap_name, vmin=-5, vmax=35, ax=0)
    mp.pcolor(lon, lat, minsst, cmap=cmap_name, vmin=-5, vmax=35, ax=1)
    im_sst = mp._im
    mp.pcolor(lon, lat, maxsst-minsst, cmap="Reds", ax=2)
    im_range = mp._im

    for ax in range(3):
        mp.nice(ax=ax)

    mp.fig.canvas.draw()
    pos0 = mp.axes.flat[0].get_position()
    pos1 = mp.axes.flat[1].get_position()
    pos2 = mp.axes.flat[2].get_position()

    cb_x = pos0.x1 + 0.01
    cb_w = 0.03

    cax1 = mp.fig.add_axes([cb_x, pos1.y0, cb_w, pos0.y1 - pos1.y0])
    cb1 = mp.fig.colorbar(im_sst, cax=cax1, orientation="vertical")
    cb1.ax.set_ylabel("OSTIA SST (°C)")

    cax2 = mp.fig.add_axes([cb_x, pos2.y0, cb_w, pos2.height])
    cb2 = mp.fig.colorbar(im_range, cax=cax2, orientation="vertical")
    cb2.ax.set_ylabel("Max $-$ Min SST (°C)")
    if vector_figs:
        plt.savefig("figs/ostia_sst.pdf", dpi=600, bbox_inches="tight")
    plt.savefig("figs/ostia_sst.png", dpi=600, bbox_inches="tight")
