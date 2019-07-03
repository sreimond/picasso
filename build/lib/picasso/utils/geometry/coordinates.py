# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 13:34:18 2018

@author: sreimond
"""

from ..constants import defaults
import numpy as np

def geocentric2geodetic( x, y, z, 
                        a=defaults.EARTH_EQUATORIAL_RADIUS(),
                        b=defaults.EARTH_POLAR_RADIUS() ):
    """
    The `geocentric2geodetic` function converts geocentric  cartesian 
    coordinates (x,y,z) into geodetic coordinates (lon,lat,h). The reference 
    ellipsoid can be specified with the optional arguments a and b.
    
    Returns three values: lon (float), lat (float), height (float)
    
    Requires the `constants` package.
    """ 
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    e = np.sqrt((a*a-b*b)/(a*a))
    p = np.sqrt(x*x+y*y)
    lon = np.arctan2(y,x)
    lat1 = np.arctan2(z,p*(1.0-(e*e)))
    dlat = np.inf
    it = 0
    while np.any(dlat>1e-10) and (it<1e4):
        it += 1
        N = a/np.sqrt(1.0-e*e*np.sin(lat1)*np.sin(lat1))
        h = p/np.cos(lat1)-N
        lat2 = np.arctan2(z,p*(1.0-e*e*N/(N+h)))
        dlat = np.fabs(lat1-lat2)
        lat1 = lat2
    lon = lon  * 180.0 / np.pi
    lat = lat1 * 180.0 / np.pi
    h = p/np.cos(lat1)-N
    return lon, lat, h

    
def geodetic2geocentric( lon, lat, h, 
                        a=defaults.EARTH_EQUATORIAL_RADIUS(),
                        b=defaults.EARTH_POLAR_RADIUS() ):
    """
    The `geodetic2geocentric` function converts geodetic coordinates 
    (lon,lat,h) into geocentric cartesian coordinates (x,y,z). The reference 
    ellipsoid can be specified with the optional arguments a and b.
    
    Returns three values: x (float), y (float), z (float)
    
    Requires the `constants` package.
    """     
    lon = np.array(lon) * np.pi/180.0
    lat = np.array(lat) * np.pi/180.0
    h = np.array(h)
    e = np.sqrt((a*a-b*b)/(a*a))
    N = a/np.sqrt(1.0-e*e*np.sin(lat)*np.sin(lat))
    x = (N+h) * np.cos(lat) * np.cos(lon)
    y = (N+h) * np.cos(lat) * np.sin(lon)
    z = (b*b/(a*a)*N+h) * np.sin(lat)
    return x, y, z
    
    