
import numpy as np
from tqdm import tqdm
from abstemp.seagrid import ostia
import pandas as pd
import matplotlib.pyplot as plt

import abstemp
from abstemp import data as abstemp_data
from abstemp.seagrid import cmip6
from abstemp.reg_calculations import sst_cmip6, sst_ostia

def generate_regdegvel_df() -> pd.DataFrame:
    """Build the full per-region degree-velocity DataFrame for three time periods.

    Loads the mintmat dataset, extracts region metadata, and computes
    degree-velocity columns for 1985–1990, 2019–2023, and 2095–2100 by
    calling :func:`regvel_1985`, :func:`regvel_2019`, and
    :func:`regvel_2095`.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by region number with columns:
        ``nreg``, ``reglat``, ``reglon``, ``longhurst_basins``,
        ``longhurst_regions``, ``longhurst_biomes``, and for each year
        ``{year}_max_sst``, ``{year}_regdegvel``, ``{year}_regdegbas``,
        ``{year}_regdegfrm``.
    """
    gl = abstemp.open_mintmat_ds()
    df = gl[[
        "nreg",
        "reglat", "reglon",
        "longhurst_basins",
        "longhurst_regions",
        "longhurst_biomes"
    ]].to_dataframe()
    df =  regvel_1985(df)
    df =  regvel_2019(df)
    df =  regvel_2095(df)
    return df

def movedegree(degree=1, sst=None):
    """Calculate the time it takes to move x degree"""
    gl = abstemp.open_mintmat_ds()
    mintmat = gl.dijkmintmat
    mintmat = mintmat.where(mintmat>=3, np.nan)
    sst = gl.abs_max_sst if sst is None else sst
    regdegvel = sst * 100
    regdegbas = sst * 0
    regdegfrm = sst * 0
    regids = np.arange(len(regdegvel))
    basins = gl.longhurst_basins.values

    def find_nearest_warmer_pos(r):
        sst_r = float(sst.iloc[r] if hasattr(sst, "iloc") else sst[r])
        if not np.isfinite(sst_r):
            return np.nan, np.nan, np.nan
        try:
            mask = (sst >= sst_r + degree).values
            travel = mintmat.isel(from_reg=r).values
            idx = np.nanargmin(travel[mask])  # raises ValueError if empty or all-NaN
            pos = regids[mask][idx]
            return travel[pos], basins[pos], pos
        except ValueError:
            return np.nan, np.nan, np.nan

    for r in tqdm(regids):
        vel, bas, pos = find_nearest_warmer_pos(r)
        regdegvel[r] = vel
        regdegbas[r] = bas
        regdegfrm[r] = pos
    regdegvel.attrs["Description"] = f"Shortest time to move {degree}°C"
    regdegbas.attrs["Description"] = f"Basin particle came from"
    regdegfrm.attrs["Description"] = f"Region particle came from"
    return regdegvel,regdegbas,regdegfrm

def regvel_2095(df: "pd.DataFrame | None" = None) -> pd.DataFrame:
    """Add 2095–2100 degree-velocity columns to *df* using EC-Earth3-CC SSP5-8.5.

    Parameters
    ----------
    df : pandas.DataFrame or None, optional
        Per-region DataFrame to extend.  If None, :func:`abstemp.read_regdegvel`
        is used.

    Returns
    -------
    pandas.DataFrame
        Input DataFrame with four new columns:
        ``2095_max_sst``, ``2095_regdegvel``, ``2095_regdegbas``,
        ``2095_regdegfrm``.
    """
    df = abstemp.read_regdegvel() if df is None else df
    sstreg = sst_cmip6.to_reg(abstemp_data.open_warmest_2095())
    df["2095_max_sst"] = sstreg["sst"]
    regdegvel,regdegbas,regdegfrm = movedegree(1, df["2095_max_sst"])
    df["2095_regdegvel"] = regdegvel
    df["2095_regdegbas"] = regdegbas
    df["2095_regdegfrm"] = regdegfrm
    return df

def regvel_1985(df: "pd.DataFrame | None" = None) -> pd.DataFrame:
    """Add 1985–1990 degree-velocity columns to *df* using OSTIA reanalysis.

    Parameters
    ----------
    df : pandas.DataFrame or None, optional
        Per-region DataFrame to extend.  If None, :func:`abstemp.read_regdegvel`
        is used.

    Returns
    -------
    pandas.DataFrame
        Input DataFrame with four new columns:
        ``1985_max_sst``, ``1985_regdegvel``, ``1985_regdegbas``,
        ``1985_regdegfrm``.
    """
    df = abstemp.read_regdegvel() if df is None else df
    sstreg = sst_ostia.to_reg(abstemp_data.open_warmest_1985().maxarr.values)
    df["1985_max_sst"] = sstreg.sst
    regdegvel,regdegbas,regdegfrm = movedegree(1, df["1985_max_sst"])
    df["1985_regdegvel"] = regdegvel
    df["1985_regdegbas"] = regdegbas
    df["1985_regdegfrm"] = regdegfrm
    return df

def regvel_2019(df: "pd.DataFrame | None" = None) -> pd.DataFrame:
    """Add 2019–2023 degree-velocity columns to *df* using OSTIA NRT product.

    Parameters
    ----------
    df : pandas.DataFrame or None, optional
        Per-region DataFrame to extend.  If None, :func:`abstemp.read_regdegvel`
        is used.

    Returns
    -------
    pandas.DataFrame
        Input DataFrame with four new columns:
        ``2019_max_sst``, ``2019_regdegvel``, ``2019_regdegbas``,
        ``2019_regdegfrm``.
    """
    df = abstemp.read_regdegvel() if df is None else df
    sstreg = sst_ostia.to_reg(abstemp_data.open_warmest_2019().maxarr.values)
    df["2019_max_sst"] = sstreg.sst
    regdegvel,regdegbas,regdegfrm = movedegree(1, df["2019_max_sst"])
    df["2019_regdegvel"] = regdegvel
    df["2019_regdegbas"] = regdegbas
    df["2019_regdegfrm"] = regdegfrm
    return df

def regdegvel_hist(df: "pd.DataFrame | None" = None) -> None:
    """Plot area-weighted histograms of degree-velocity for 1–4 °C steps.

    Computes :func:`movedegree` for degree increments 1–4 using 2019
    SST values and draws a scatter+line figure with area-weighted counts
    on the y-axis and time-to-travel (years) on the x-axis.

    Parameters
    ----------
    df : pandas.DataFrame or None, optional
        Per-region DataFrame.  If None, :func:`abstemp.read_regdegvel`
        is used.

    Returns
    -------
    None
    """
    if df is None:
        df = abstemp.read_regdegvel()
    area = abstemp.reg_calculations.arr_to_regvec(abstemp.reg_calculations.reg_area())
    vel2019_1yr,_,_ = movedegree(1, df["2019_max_sst"])
    vel2019_2yr,_,_ = movedegree(2, df["2019_max_sst"])
    vel2019_3yr,_,_ = movedegree(3, df["2019_max_sst"])
    vel2019_4yr,_,_ = movedegree(4, df["2019_max_sst"])

    def plot_vel_hist(vel, color, label):
        counts, edges = np.histogram(vel/365,np.arange(0,10,0.25), weights=area)
        counts = np.array([counts[0]] + list(counts))/1000/1000
        plt.scatter(edges[1:], counts[1:], c=color, s=3)
        plt.plot(edges, counts, lw = 0.5, c=color, label=label)
        plt.ylabel("Area (km²)")
        plt.xlabel("Time to move from warmer to cooler (years)")

    plt.clf()
    plot_vel_hist(vel2019_1yr, "tab:blue", "1°C")
    plot_vel_hist(vel2019_2yr, "tab:orange", "2°C")
    plot_vel_hist(vel2019_3yr, "tab:green", "3°C")
    plot_vel_hist(vel2019_4yr, "tab:red", "4°C")
    plt.xlim(0,8)
    ylim = list(plt.gca().get_ylim())
    ylim[0] = 0
    plt.ylim(ylim)
    plt.legend()
