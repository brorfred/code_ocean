"""Phytoplankton growth-rate models and associated plotting helpers.

Provides parametric temperature–growth functions (Blackford, Norberg-Eppley)
together with their analytical first and second derivatives, Monte-Carlo noise
estimates, individual plotting routines, and a combined two-panel figure.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

NT = 200
T = np.linspace(0,50,NT)

def fN(x, Z=28, w=40.2):
    """
    Norberg unimodal growth model with fixed exponential scaling.

    Parameters
    ----------
    x : float or array_like
        Temperature in degrees Celsius.
    Z : float, optional
        Optimal temperature (centre of the unimodal curve) in degrees Celsius.
        Default is 28.
    w : float, optional
        Thermal niche width in degrees Celsius. Default is 40.2.

    Returns
    -------
    float or ndarray
        Growth rate at temperature *x*.
    """
    return (1 - ( (x-Z)/w)**2 )*0.59*np.exp(0.0633*x)

def norberg_eppley(T):
    """
    Norberg-Eppley growth model.

    Parameters
    ----------
    T : float or array_like
        Temperature in degrees Celsius.

    Returns
    -------
    float or ndarray
        Growth rate.

    References
    ----------
    Norberg (2004)
    Ribalt et al. (2025)
    """
    t0=20
    a=0.70
    b=0.21
    Tc=11.7
    w=40.2

    return a*np.exp(b*(T-t0)) * (1 - ((T - Tc) / (w / 2.0)) ** 2)

def norberg_eppley_dT(T, T0=20.0, a=0.70, b=0.21, Tc=11.7, w=40.2):
    """
    First derivative of the Norberg-Eppley growth model with respect to temperature.

    Parameters
    ----------
    T : float or array_like
        Temperature in degrees Celsius.
    T0 : float, optional
        Reference temperature for the exponential term in degrees Celsius.
        Default is 20.0.
    a : float, optional
        Scaling coefficient. Default is 0.70.
    b : float, optional
        Exponential growth rate coefficient. Default is 0.21.
    Tc : float, optional
        Optimal temperature (centre of the unimodal curve) in degrees Celsius.
        Default is 11.7.
    w : float, optional
        Thermal niche width in degrees Celsius. Default is 40.2.

    Returns
    -------
    float or ndarray
        First derivative of growth rate with respect to temperature at *T*.

    See Also
    --------
    norberg_eppley : The function whose derivative is computed here.
    norberg_eppley_dT2 : Second derivative of the Norberg-Eppley growth model.
    """
    exp_term = np.exp(b * (T - T0))
    quad_term = 1.0 - ((T - Tc) / (w / 2.0)) ** 2

    d_quad = -8.0 * (T - Tc) / (w ** 2)

    return a * exp_term * (
        b * quad_term + d_quad
    )


def norberg_eppley_dT2(T, T0=20.0, a=0.70, b=0.21, Tc=11.7, w=40.2):
    """
    Second derivative of the Norberg-Eppley growth model with respect to temperature.

    Useful for locating the inflection points of the growth curve.

    Parameters
    ----------
    T : float or array_like
        Temperature in degrees Celsius.
    T0 : float, optional
        Reference temperature for the exponential term in degrees Celsius.
        Default is 20.0.
    a : float, optional
        Scaling coefficient. Default is 0.70.
    b : float, optional
        Exponential growth rate coefficient. Default is 0.21.
    Tc : float, optional
        Optimal temperature (centre of the unimodal curve) in degrees Celsius.
        Default is 11.7.
    w : float, optional
        Thermal niche width in degrees Celsius. Default is 40.2.

    Returns
    -------
    float or ndarray
        Second derivative of growth rate with respect to temperature at *T*.

    See Also
    --------
    norberg_eppley : The function whose second derivative is computed here.
    norberg_eppley_dT : First derivative of the Norberg-Eppley growth model.
    """
    exp_term = np.exp(b * (T - T0))
    quad_term = 1.0 - ((T - Tc) / (w / 2.0)) ** 2

    d_quad = -8.0 * (T - Tc) / (w ** 2)
    d2_quad = -8.0 / (w ** 2)

    return a * exp_term * (
        b**2 * quad_term
        + 2.0 * b * d_quad
        + d2_quad
    )

def norberg_eppley_noise(T, T0=20.0, a=0.70, b=0.21, Tc=11.7, w=40.2):
    """
    Mean Norberg-Eppley growth rate under temperature variability.

    For each point in a 200-step temperature grid from 0 to 50 °C, draws
    500 000 normally distributed temperature samples (std = 5 °C) and returns
    the mean growth rate, capturing Jensen's inequality effects.

    Parameters
    ----------
    T : array_like
        Nominal temperature array (unused directly; grid is rebuilt internally).
    T0 : float, optional
        Reference temperature for the exponential term in degrees Celsius.
        Default is 20.0.
    a : float, optional
        Scaling coefficient. Default is 0.70.
    b : float, optional
        Exponential growth rate coefficient. Default is 0.21.
    Tc : float, optional
        Optimal temperature in degrees Celsius. Default is 11.7.
    w : float, optional
        Thermal niche width in degrees Celsius. Default is 40.2.

    Returns
    -------
    ndarray of shape (200,)
        Mean growth rate at each temperature grid point.

    See Also
    --------
    norberg_eppley : Deterministic version of this growth model.
    """
    NT = 200
    fENm = np.zeros(NT)
    fENstd = np.zeros(NT)
    Tstd = 5

    for y in range(0,NT):
        Trand = np.random.normal(loc=y*50/NT, scale=Tstd, size=500000)
        GREN =  norberg_eppley(Trand)
        #fBm[y] = np.mean( GRB[GRB > 0] )
        fENm[y] = np.mean( GREN )
    return fENm

def norberg_eppley_plot(ax=None):
    """
    Plot the Norberg-Eppley growth model with its derivatives and noise envelope.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None, a new 10×6 inch figure is created.
        Default is None.

    Returns
    -------
    None
    """
    NT = 200
    T = np.linspace(0,50,NT)

    if ax is None:
        plt.figure(figsize=(10,6))
        ax = plt.gca()

    ax.axhline(y=0, color='lightgrey', linestyle='--', linewidth=2)
    ax.plot(T,norberg_eppley(T), lw = 2, label = r'$\mu$')
    ax.plot(T,10*norberg_eppley_dT(T), lw = 1, label = 'First derivative')
    ax.plot(T,10*norberg_eppley_dT2(T), lw = 1, label = 'Second derivative')
    ax.plot(T,norberg_eppley_noise(T), lw = 2, c = 'red', label = r'$\mu$ with noise')
    ax.set_ylabel(r"Relative growth rate ($\mu$)")
    ax.set_xlabel("Temperature [°C]")
    ax.set_title("Epply Norberg")
    ax.legend(loc = 'lower left', fontsize=8, ncols=2)
    ax.set_xlim(0,35)
    ax.set_ylim(-2,3);
    #plt.savefig("EpplyNorberg.pdf")


def eppley(T, a, b, Tc, w):
    """
    Eppley-style unimodal growth model.

    Parameters
    ----------
    T : float or array_like
        Temperature in degrees Celsius.
    a : float
        Scaling coefficient (maximum growth rate at reference temperature).
    b : float
        Exponential growth rate coefficient.
    Tc : float
        Optimal temperature (centre of the unimodal curve) in degrees Celsius.
    w : float
        Thermal niche width in degrees Celsius.

    Returns
    -------
    float or ndarray
        Growth rate at temperature *T*.
    """
    return a*np.exp(b*(T-28)) * (1 - ((T - Tc) / (w / 2.0)) ** 2)

def blackford(x, Z=28, w=15, Q10=2):
    """
    Blackford growth model based on Q10 temperature scaling.

    Computes the difference between a slow anabolic term and a fast catabolic
    term, producing a unimodal growth-versus-temperature curve.

    Parameters
    ----------
    x : float or array_like
        Temperature in degrees Celsius.
    Z : float, optional
        Reference temperature for the catabolic (loss) term in degrees
        Celsius. Default is 28.
    w : float, optional
        Thermal niche width in degrees Celsius. Currently unused in the
        formula but retained for interface consistency. Default is 15.
    Q10 : float, optional
        Q10 temperature coefficient. Default is 2.

    Returns
    -------
    float or ndarray
        Net growth rate at temperature *x*.
    """
    return Q10**((x-10)/10) - Q10**((x-Z)/3)

def blackford_dT(x, Z=28, Q10=2):
    """
    First derivative of the Blackford growth model.

    Parameters
    ----------
    x : float or array_like
        Temperature in degrees Celsius.
    Z : float, optional
        Reference temperature for the catabolic term in degrees Celsius.
        Default is 28.
    Q10 : float, optional
        Q10 temperature coefficient. Default is 2.

    Returns
    -------
    float or ndarray
        First derivative of net growth rate with respect to temperature at *x*.

    See Also
    --------
    blackford : The function whose derivative is computed here.
    blackford_dT2 : Second derivative of the Blackford growth model.
    """
    lnQ10 = np.log(Q10)

    return (
        (lnQ10 / 10) * Q10**((x - 10) / 10)
        - (lnQ10 / 3) * Q10**((x - Z) / 3)
    )

def blackford_dT2(x, Z=28, Q10=2):
    """
    Second derivative of the Blackford growth model.

    Useful for locating the inflection points of the growth curve.

    Parameters
    ----------
    x : float or array_like
        Temperature in degrees Celsius.
    Z : float, optional
        Reference temperature for the catabolic term in degrees Celsius.
        Default is 28.
    Q10 : float, optional
        Q10 temperature coefficient. Default is 2.

    Returns
    -------
    float or ndarray
        Second derivative of net growth rate with respect to temperature at *x*.

    See Also
    --------
    blackford : The function whose second derivative is computed here.
    blackford_dT : First derivative of the Blackford growth model.
    """
    lnQ10 = np.log(Q10)

    return (
        (lnQ10**2 / 100) * Q10**((x - 10) / 10)
        - (lnQ10**2 / 9) * Q10**((x - Z) / 3)
    )

def blackford_noise(x, Z=28, Q10=2):
    """
    Mean Blackford growth rate under temperature variability.

    For each point in the module-level temperature grid, draws 200 000 normally
    distributed temperature samples (std = 5 °C) and returns the mean growth
    rate, capturing Jensen's inequality effects.

    Parameters
    ----------
    x : array_like
        Nominal temperature array (unused directly; the module-level grid is
        used internally).
    Z : float, optional
        Reference temperature for the catabolic term in degrees Celsius.
        Default is 28.
    Q10 : float, optional
        Q10 temperature coefficient. Default is 2.

    Returns
    -------
    ndarray of shape (NT,)
        Mean growth rate at each temperature grid point.

    See Also
    --------
    blackford : Deterministic version of this growth model.
    """
    fBm = np.zeros(NT)
    #fBstd = np.zeros(NT)
    #fNm = np.zeros(NT)
    #fNstd = np.zeros(NT)
    Tstd = 5

    for y in range(0,NT):
        Trand = np.random.normal(loc=y*50/NT, scale=Tstd, size=200000)
        GRB =  blackford(Trand)
        #fBm[y] = np.mean( GRB[GRB > 0] )
        fBm[y] = np.mean( GRB )
        #fBstd[y] = np.std( fB(Trand) )
        GRN =  fN(Trand)
        #fNm[y] = np.mean( GRN )
        #fNstd[y] = np.std( fN(Trand) )
    return fBm

def blackford_plot(ax=None):
    """
    Plot the Blackford growth model together with its first and second derivatives.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None, a new 10×6 inch figure is created.
        Default is None.

    Returns
    -------
    None
    """
    if ax is None:
        plt.figure(figsize=(10,6))
        ax = plt.gca()
    ax.axhline(y=0, color='lightgrey', linestyle='--', linewidth=2)
    ax.plot(T,blackford(T), lw = 2, label = r'$\mu$')
    ax.plot(T,10*blackford_dT(T), lw = 1, label = 'First derivative')
    ax.plot(T,10*blackford_dT2(T), lw = 1, label = 'Second derivative')
    ax.plot(T,blackford_noise(T), lw = 2, c = 'red', label = r'$\mu$ with noise')
    ax.set_ylabel(r"Relative growth rate ($\mu$)")
    ax.set_xlabel("Temperature [°C]")
    ax.set_title("Blackford et al. (2004)")
    ax.legend(loc = 'lower left', fontsize=8, ncols=2)
    ax.set_xlim(0,40)
    ax.set_ylim(-2,3);
    #plt.savefig("Blackford.pdf")

def plot():
    """
    Produce a combined two-panel figure of the Blackford and Norberg-Eppley growth models.

    Saves the figure to ``growth_models.pdf`` in the ``figs`` directory.

    Returns
    -------
    None
    """
    fig,axes = plt.subplots(2,1, sharex=True, sharey=True)
    blackford_plot(ax=axes[0])
    axes[0].set_xlabel("")
    norberg_eppley_plot(ax=axes[1])
    plt.savefig("figs/growth_models.pdf", bbox_inches='tight')


def growth_observations() -> None:
    """Plot median phytoplankton growth rates alongside the Blackford model.

    Reads per-temperature growth rates from
    ``src/abstemp/data/growth_rates.csv``, computes the per-integer-degree
    median, and draws a two-panel figure: scatter of observed medians (top)
    and the Blackford growth-model curve (bottom).

    Returns
    -------
    None
    """

    from pathlib import Path
    df = pd.read_csv(Path(__file__).parent.parent / "data" / "growth_rates.csv")
    df["temp"] = df["temperature"].astype(int)
    gr = df[["temp", "r"]].groupby("temp").median()


    fig,axes = plt.subplots(2,1, sharex=True)
    axes[0].scatter(gr.index, gr.r, 5)
    axes[0].set_ylabel("Growth rate")
    plot(ax=axes[1])
    axes[1].set_xlabel("Temperature (°C)")
