# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 13:40:54 2018

@author: sreimond
"""

from . import GRS80

def EARTH_EQUATORIAL_RADIUS():
    """
    The `EARTH_EQUATORIAL_RADIUS` function returns the default Earth's 
    equatorial radius (unit is meters).
    """
    return GRS80.EQUATORIAL_RADIUS()

def EARTH_GM():
    """
    The `EARTH_GM` function returns the default geocentric gravitational 
    constant (unit is m3/s2).
    """
    return GRS80.GM()    
    
def EARTH_INVERSE_FLATTENING():
    """
    The `EARTH_INVERSE_FLATTENING` function returns the default Earth's inverse
    flattening factor 1/f.
    """
    return GRS80.INVERSE_FLATTENING()
    
def EARTH_J2():
    """
    The `EARTH_J2` function returns the default dynamic form factor J2.
    """
    return GRS80.J2()
    
def EARTH_POLAR_RADIUS():
    """
    The `EARTH_POLAR_RADIUS` function returns the default Earth's polar radius 
    (unit is meters).
    """
    return GRS80.POLAR_RADIUS()
    
def EARTH_MEAN_RADIUS():
    """
    The `EARTH_MEAN_RADIUS` function returns the default Earth's mean radius 
    (unit is meters).
    """
    return GRS80.MEAN_RADIUS()
    
def EARTH_OMEGA():
    """
    The `EARTH_OMEGA` function returns the default Earth's angular velocity 
    (unit is rad/s).
    """
    return GRS80.OMEGA()

    
    