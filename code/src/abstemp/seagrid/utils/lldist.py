"""Haversine and spherical distance calculations on the Earth's surface.

All distances returned in metres unless noted.
"""

from numpy import *
import numpy as np


pi180 =  pi/180
earth_radius = 6371.0*1000
 

def lldist(lon,lat):
    """Calculate cumulative distance along a 1-D trajectory.

    Uses the haversine formula to compute the great-circle distance
    between consecutive coordinate pairs.

    Parameters
    ----------
    lon : array_like
        1-D array of longitudes in degrees.
    lat : array_like
        1-D array of latitudes in degrees.

    Returns
    -------
    numpy.ndarray
        Array of cumulative distances in metres, same shape as the
        input. The first element is 0.
    """
    lat = lat * pi180
    dlon = (lon[1:] - lon[:-1])*pi180 
    dlat = (lat[1:] - lat[:-1])

    a = (sin(dlat/2))**2 + cos(lat[:-1]) * cos(lat[1:]) * (sin(dlon/2))**2
    angles = lon * 0
    angles[1:] = 2 * arctan2( sqrt(a), sqrt(1-a) )
    return earth_radius * angles

def ll2dist(lon,lat):
    """Calculate haversine distance between two rows of a 2-D array.

    Computes the great-circle distance from row 0 to row 1 for each
    column of the input arrays.

    Parameters
    ----------
    lon : array_like
        2-D array of longitudes in degrees with at least two rows.
    lat : array_like
        2-D array of latitudes in degrees with at least two rows.

    Returns
    -------
    numpy.ndarray
        1-D array of distances in metres, one value per column.
    """
    lat = lat * pi180
    dlon = (lon[1,:] - lon[0,:])*pi180 
    dlat = (lat[1,:] - lat[0,:])

    a = (sin(dlat/2))**2 + cos(lat[0,:]) * cos(lat[1,:]) * (sin(dlon/2))**2
    angles = 2 * arctan2( sqrt(a), sqrt(1-a) )
    return earth_radius * angles

def ll2dist2vec(lonvec1, latvec1, lonvec2, latvec2):
    """Calculate pairwise haversine distances between corresponding elements of two coordinate vectors.

    Parameters
    ----------
    lonvec1 : array_like
        Longitudes of the first set of points, in degrees.
    latvec1 : array_like
        Latitudes of the first set of points, in degrees.
    lonvec2 : array_like
        Longitudes of the second set of points, in degrees.
    latvec2 : array_like
        Latitudes of the second set of points, in degrees.

    Returns
    -------
    numpy.ndarray
        Distances in metres between each pair of corresponding points.
    """
    latvec1 = np.deg2rad(latvec1)
    latvec2 = np.deg2rad(latvec2)
    lonvec1 = np.deg2rad(lonvec1)
    lonvec2 = np.deg2rad(lonvec2)

    dlon = (lonvec2 - lonvec1)
    dlat = (latvec2 - latvec1)

    a = (sin(dlat/2))**2 + cos(latvec1) * cos(latvec2) * (sin(dlon/2))**2
    angles = 2 * arctan2( sqrt(a), sqrt(1-a) )
    return earth_radius * angles


def spherical_dist(pos1, pos2, r=3958.75):
    """Calculate spherical distance using the law of cosines.

    Despite the default value of *r* (in miles), the function uses the
    module-level ``earth_radius`` constant (in metres) for the
    calculation, so distances are returned in metres.

    Parameters
    ----------
    pos1 : array_like
        Array with latitude and longitude in degrees along the last
        axis (i.e. ``pos1[..., 0]`` is latitude, ``pos1[..., 1]`` is
        longitude).
    pos2 : array_like
        Array with latitude and longitude in degrees, same layout as
        *pos1*.
    r : float, optional
        Earth radius parameter (default ``3958.75``). Note that the
        actual computation uses the ``earth_radius`` module constant
        (6 371 000 m).

    Returns
    -------
    numpy.ndarray
        Distances in metres.
    """
    pos1 = np.deg2rad(pos1)
    pos2 = np.deg2rad(pos2)
    cos_lat1 = np.cos(pos1[..., 0])
    cos_lat2 = np.cos(pos2[..., 0])
    cos_lat_d = np.cos(pos1[..., 0] - pos2[..., 0])
    cos_lon_d = np.cos(pos1[..., 1] - pos2[..., 1])
    return earth_radius * np.arccos(cos_lat_d - cos_lat1 *
                                    cos_lat2 * (1 - cos_lon_d))


def distance_matrix(lonvec, latvec):
    """Calculate a full distance matrix for all pairwise combinations.

    Parameters
    ----------
    lonvec : array_like
        1-D array of longitudes in degrees.
    latvec : array_like
        1-D array of latitudes in degrees.

    Returns
    -------
    numpy.ndarray
        2-D square array of distances in metres where element
        ``(i, j)`` is the distance between point *i* and point *j*.
    """
    latmat1,latmat2 = np.meshgrid(latvec, latvec)
    lonmat1,lonmat2 = np.meshgrid(lonvec, lonvec)
    pos1 = np.vstack((latmat1.flatten(), lonmat1.flatten())).T
    pos2 = np.vstack((latmat2.flatten(), lonmat2.flatten())).T
    
    dist =  spherical_dist(pos1,pos2)
    return dist.reshape(int(np.sqrt(len(dist))),-1)
