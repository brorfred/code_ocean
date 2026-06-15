
import matplotlib.pyplot as plt
import numpy as np

from abstemp.global_sst_histograms import read_cmip6_hists, read_ostia_hists

sstvec = np.arange(-4,40,0.25)

def main_histograms(include_c24: bool = False) -> None:
    """Plot and save the global SST histogram comparison figure.

    Draws area-weighted warmest-month SST histograms for three climate
    periods on a single panel: OSTIA 1985–1990, OSTIA 2019–2023, and the
    CMIP6 SSP5-8.5 ensemble mean for 2095–2100 (individual model curves
    shown in light red).  Optionally overlays the SSP2-4.5 ensemble mean.

    Saves the figure to ``figs/sst_hist_1985_2019_2095.pdf`` at 900 dpi.

    Parameters
    ----------
    include_c24 : bool, optional
        If True, also plot the SSP2-4.5 ensemble mean.  Default is False.

    Returns
    -------
    None
    """
    c58 = read_cmip6_hists(experiment="ssp585")
    c24 = read_cmip6_hists(experiment="ssp245")
    ost = read_ostia_hists()

    def ynorm(vec):
        return vec

    plt.clf()
    plt.stairs(ynorm(ost["1985-1990"]), sstvec, lw=1, color="tab:purple", label="OSTIA 1985-90", baseline=None)
    plt.stairs(ynorm(ost["2019-2023"]), sstvec, lw=1, color="tab:blue", label="OSTIA 2019-23", baseline=None)
    for model in c58.keys():
        plt.stairs(ynorm(c58[model]), sstvec, lw=0.5, alpha=0.1, color="tab:red")
    #for y in yarr: plt.stairs(y, x, color="0.5", alpha=0.5, lw=0.5)
    plt.stairs(ynorm(np.nanmean(ynorm(c58), axis=1)), sstvec, lw=1, color="tab:red", label="CMIP6 2095-2100")
    if include_c24:
        plt.stairs(ynorm(np.nanmean(ynorm(c24), axis=1)), sstvec, lw=1, color="tab:orange", label="CMIP6 2095-2100 ssp245")
    plt.legend()

    plt.xlim(-2,40)
    plt.ylim(-0.1e7, 1e7)
    plt.ylabel("Area (km$^2$)")
    plt.xlabel("Max monthly SST (°C)")

    plt.savefig("figs/sst_hist_1985_2019_2095.pdf", dpi=900)
#ds = warmest_month.open_dataset(f"1985-01-01", f"1990-12-31", #ds"sst", nrt=True)
#ds = warmest_month.open_dataset(f"2018-01-01", f"2023-12-31", "sst", nrt=True)
