# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 12:18:30 2018

@author: sreimond
"""

def EQUATORIAL_RADIUS():
    """
    The `EQUATORIAL_RADIUS` function returns the Earth's equatorial radius 
    (unit is meters) according to the GRS80 reference system.
    """
    return 6378137.0
    
def GM():
    """
    The `GM` function returns the geocentric gravitational constant (unit 
    is m3/s2) according to the GRS80 reference system.
    """
    return 3.986005e14
    
def INVERSE_FLATTENING():
    """
    The `INVERSE_FLATTENING` function returns the Earth's inverse flattening 
    factor 1/f according to the GRS80 reference system (unitless).
    """
    return 298.257222101
    
def J2():
    """
    The `J2` function returns the dynamic form factor J2 according to the 
    GRS80 reference system.
    """
    return 1.08263e-3
    
def POLAR_RADIUS():
    """
    The `POLAR_RADIUS` function returns the Earth's polar radius (unit is 
    meters) according to the GRS80 reference system.
    """
    return EQUATORIAL_RADIUS()*(1.0-1.0/INVERSE_FLATTENING())
    
def MEAN_RADIUS():
    """
    The `MEAN_RADIUS` function returns the Earth's mean radius (unit is meters)
    based on IUGG definition and the GRS80 reference system.
    """
    return (2.0*EQUATORIAL_RADIUS()+POLAR_RADIUS())/3.0
    
def OMEGA():
    """
    The `OMEGA` function returns the Earth's angular velocity (unit is rad/s) 
    according to the GRS80 reference system.
    """
    return 7.292115e-05
    
    
    
    
    
    