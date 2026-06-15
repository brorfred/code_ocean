"""Grid generation utilities for pyresample.

Provides functions to create standard NASA grids and grids from xarray
lat/lon coordinates.
"""

import numpy as np
import pyresample as pr


def setup_plate_carree_grid(lon1, lon2, lat1, lat2, width, height, **kw):
    area_def = pr.geometry.AreaDefinition(
        area_id="global_005deg",
        description="Global 0.05-degree lat/lon grid",
        proj_id="longlat_wgs84",
        projection={
            "proj": "longlat",
            "datum": "WGS84",
            "no_defs": True,
        },
        width=width,  # xsize: 360 / 0.05
        height=height,  # ysize: 180 / 0.05
        area_extent=(lon1, lat1, lon2, lat2),
    )
    return area_def


def setup_nasa_grid(resolution="4km", **kw):
    """Create a pyresample grid using the standard NASA equal-angle projection.

    Builds a global equal-angle latitude/longitude grid at the requested
    resolution and returns it as a ``pyresample.geometry.GridDefinition``
    with additional index attributes attached.

    Parameters
    ----------
    resolution : {'9km', '4km', '1km'}, optional
        Desired grid resolution. Default is ``'4km'``.
    **kw
        Additional keyword arguments (currently unused).

    Returns
    -------
    pyresample.geometry.GridDefinition
        Grid definition with the following extra attributes:

        - **ivec** -- 1-D array of column indices.
        - **jvec** -- 1-D array of row indices.
        - **iarr** -- 2-D meshgrid of column indices.
        - **jarr** -- 2-D meshgrid of row indices.

    Raises
    ------
    ValueError
        If *resolution* is not one of the accepted values.
    """
    if resolution == "9km":
        i0t, imt, j0t, jmt = (0000, 4320, 0, 2160)
    elif resolution == "4km":
        i0t, imt, j0t, jmt = (0000, 8640, 0, 4320)
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


def setup_latlon_grid(ds, **kw):
    """Create a pyresample grid from an xarray Dataset with lat/lon coordinates.

    Uses the ``lat`` and ``lon`` coordinates of *ds* to build a meshgrid
    and wraps it in a ``pyresample.geometry.GridDefinition``.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset containing ``lat`` and ``lon`` coordinate variables.
    **kw
        Additional keyword arguments (currently unused).

    Returns
    -------
    pyresample.geometry.GridDefinition
        Grid definition with the following extra attributes:

        - **ivec** -- 1-D array of column indices.
        - **jvec** -- 1-D array of row indices.
        - **iarr** -- 2-D meshgrid of column indices.
        - **jarr** -- 2-D meshgrid of row indices.
    """
    lons, lats = np.meshgrid(ds.lon, ds.lat)
    grid = pr.geometry.GridDefinition(lons=lons, lats=lats)
    grid.ivec = np.arange(grid.shape[1])
    grid.jvec = np.arange(grid.shape[0])
    grid.iarr, grid.jarr = np.meshgrid(grid.ivec, grid.jvec)
    return grid
