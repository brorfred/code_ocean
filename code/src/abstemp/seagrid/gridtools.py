"""Grid setup and metric calculations for pyresample grids.

Provides functions for creating ``pyresample.GridDefinition`` objects
and computing approximate grid-cell spacing and areas.
"""

import numpy as np
import pyresample as pr

from .utils import ll2dist2vec


def setup_grid_from_latlon_vecs(latvec, lonvec, **kw):
    """Create a pyresample grid from 1-D latitude and longitude vectors.

    Parameters
    ----------
    latvec : array_like
        1-D array of latitude values.
    lonvec : array_like
        1-D array of longitude values.
    **kw : dict
        Additional keyword arguments (unused, reserved).

    Returns
    -------
    pyresample.geometry.GridDefinition
        Grid with ``ivec``, ``jvec``, ``iarr``, and ``jarr`` attributes.
    """
    lons, lats = np.meshgrid(lonvec, latvec)
    grid = pr.geometry.GridDefinition(lons=lons, lats=lats)
    grid.ivec = np.arange(grid.shape[1])
    grid.jvec = np.arange(grid.shape[0])
    grid.iarr, grid.jarr = np.meshgrid(grid.ivec, grid.jvec)
    return grid


def setup_nasa_grid(resolution="4km", **kw):
    """Create a pyresample grid using the standard NASA equal-angle projection.

    Parameters
    ----------
    resolution : {'18km','9km', '4km', '2km', '1km'}, optional
        Grid resolution.
    **kw : dict
        Additional keyword arguments (unused, reserved).

    Returns
    -------
    pyresample.geometry.GridDefinition
        Grid with ``ivec``, ``jvec``, ``iarr``, and ``jarr`` attributes.

    Raises
    ------
    ValueError
        If ``resolution`` is not one of the supported values.
    """
    if resolution == "18km":
        i0t, imt, j0t, jmt = (0000, 2160, 0, 1080)
    elif resolution == "9km":
        i0t, imt, j0t, jmt = (0000, 4320, 0, 2160)
    elif resolution == "4km":
        i0t, imt, j0t, jmt = (0000, 8640, 0, 4320)
    elif resolution == "2km":
        i0t, imt, j0t, jmt = (0000, 17280, 0, 8640)
    elif resolution == "1km":
        i0t, imt, j0t, jmt = (0000, 34560, 0, 17280)
    else:
        raise ValueError("Wrong resolution")
    incr = 360.0 / imt
    jR = np.arange(j0t, jmt)
    iR = np.arange(i0t, imt)
    latvec = (90 - jR * incr - incr / 2)[::-1]
    lonvec = -180 + iR * incr + incr / 2
    lons, lats = np.meshgrid(lonvec, latvec)
    grid = pr.geometry.GridDefinition(lons=lons, lats=lats)
    grid.ivec = np.arange(grid.shape[1])
    grid.jvec = np.arange(grid.shape[0])
    grid.iarr, grid.jarr = np.meshgrid(grid.ivec, grid.jvec)
    return grid


def dx_approx(*, lats=None, lons=None, grid=None):
    """Calculate the approximate zonal (x-direction) grid spacing.

    Uses haversine distances between adjacent grid cells in the
    longitudinal direction.

    Parameters
    ----------
    lats : numpy.ndarray, optional
        2-D array of latitudes. Required if ``grid`` is not provided.
    lons : numpy.ndarray, optional
        2-D array of longitudes. Required if ``grid`` is not provided.
    grid : pyresample.geometry.GridDefinition, optional
        Grid definition. If provided, ``lats`` and ``lons`` are
        extracted from it.

    Returns
    -------
    numpy.ndarray
        Zonal spacing in metres, same shape as the input grid.

    Raises
    ------
    TypeError
        If neither ``grid`` nor both ``lats`` and ``lons`` are provided.
    """
    if isinstance(grid, pr.geometry.GridDefinition):
        lats = grid.lats
        lons = grid.lons
    elif (lats is None) or (lons is None):
        raise TypeError("Both lats and lons must be provided")
    dx = lons * np.nan
    for i in range(dx.shape[0]):
        latvec1 = (lats[i, 0:-2] + lats[i, 1:-1]) / 2
        latvec2 = (lats[i, 2:] + lats[i, 1:-1]) / 2
        lonvec1 = (lons[i, 0:-2] + lons[i, 1:-1]) / 2
        lonvec2 = (lons[i, 2:] + lons[i, 1:-1]) / 2
        dx[i, 1:-1] = ll2dist2vec(lonvec1, latvec1, lonvec2, latvec2)
    dx[:, 0] = 2 * dx[:, 1] - dx[:, 2]
    dx[:, -1] = 2 * dx[:, -2] - dx[:, -3]
    return dx


def dy_approx(*, lats=None, lons=None, grid=None):
    """Calculate the approximate meridional (y-direction) grid spacing.

    Uses haversine distances between adjacent grid cells in the
    latitudinal direction.

    Parameters
    ----------
    lats : numpy.ndarray, optional
        2-D array of latitudes. Required if ``grid`` is not provided.
    lons : numpy.ndarray, optional
        2-D array of longitudes. Required if ``grid`` is not provided.
    grid : pyresample.geometry.GridDefinition, optional
        Grid definition. If provided, ``lats`` and ``lons`` are
        extracted from it.

    Returns
    -------
    numpy.ndarray
        Meridional spacing in metres, same shape as the input grid.

    Raises
    ------
    TypeError
        If neither ``grid`` nor both ``lats`` and ``lons`` are provided.
    """
    if isinstance(grid, pr.geometry.GridDefinition):
        lats = grid.lats
        lons = grid.lons
    elif (lats is None) or (lons is None):
        raise TypeError("Both lats and lons must be provided")
    dy = lats * np.nan
    for j in range(dy.shape[1]):
        latvec1 = (lats[0:-2, j] + lats[1:-1, j]) / 2
        latvec2 = (lats[2:, j] + lats[1:-1, j]) / 2
        lonvec1 = (lons[0:-2, j] + lons[1:-1, j]) / 2
        lonvec2 = (lons[2:, j] + lons[1:-1, j]) / 2
        dy[1:-1, j] = ll2dist2vec(lonvec1, latvec1, lonvec2, latvec2)
    dy[0, :] = 2 * dy[1, :] - dy[2, :]
    dy[-1, :] = 2 * dy[-2, :] - dy[-3, :]
    return dy


def area(*, lats=None, lons=None, grid=None):
    """Calculate the approximate area of each grid cell.

    Computes ``dx * dy`` from :func:`dx_approx` and :func:`dy_approx`.

    Parameters
    ----------
    lats : numpy.ndarray, optional
        2-D array of latitudes.
    lons : numpy.ndarray, optional
        2-D array of longitudes.
    grid : pyresample.geometry.GridDefinition, optional
        Grid definition.

    Returns
    -------
    numpy.ndarray
        Cell areas in square metres, same shape as the input grid.
    """
    dx = dx_approx(lats=lats, lons=lons, grid=grid)
    dy = dy_approx(lats=lats, lons=lons, grid=grid)
    return dx * dy
