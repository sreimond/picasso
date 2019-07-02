# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 16:32:20 2018

@author: sreimond
"""

#import os, sys
#sys.path.append(os.path.join(os.path.dirname(__file__),'..','constants'))
from ..constants import defaults
import numpy as np

def compute_kepler_orbit( a, e, inc_0, omega, Omega_0, M_0, t,
                          delta_n=0, inc_dot=0, Omega_dot=0, C_uc=0, C_us=0,
                          C_ic=0, C_is=0, C_rc=0, C_rs=0, 
                          ref_sys='ECEF',
                          GM=defaults.EARTH_GM(),
                          omega_e=defaults.EARTH_OMEGA()):
    """
    The `compute_kepler_orbit` function computes the satellite positions with 
    respect to the Earth Centered Earth Fixed (ECEF) or Space Fixed (ECSF) 
    reference frame. 
    
    Required input:
        - a       ... semi-major axis [m]
        - e       ... eccentricity
        - i       ... inclination [°]
        - omega   ... argument of perigee [°]
        - Omega_0 ... right ascension of ascending node [°]
        - M       ... mean anomaly [°]
        - t       ... time, can either be scalar or vector
    
    Optional input can be correction terms including
        - delta_n    ... mean motion correction
        - inc_dot    ... rate of inclination angle
        - Omega_dot  ... drift of node's right ascension per second
        - C_uc, C_us ... correction coefficients (argument of perigee)
        - C_rc, C_rs ... correction coefficients (geometric distance)
        - C_ic, C_is ... correction coefficients (inclination)
    or the desired output coordinate system
        - ref_sys    ... can be either 'ECEF' (default) or 'ECSF'.
    
    Requires the `constants` package.
    """ 
    # convert to radians
    deg2rad = lambda x: x*np.pi/180.0
    inc_0 = deg2rad(inc_0)
    omega = deg2rad(omega)
    Omega_0 = deg2rad(Omega_0)
    M_0 = deg2rad(M_0)
    t = np.array(t,ndmin=1)
    # computed mean motion
    n_0 = np.sqrt(GM/a**3.0)
    # corrected mean motion
    n = n_0 + delta_n
    # mean anomaly
    M = M_0+n*t
    # eccentric anomaly
    epsilon = np.random.random((np.shape(M)[0],2))
    tol = 1e-5
    j=0
    while np.max(np.ptp(epsilon,axis=1))>tol:
        j += 1
        if j==1:
            E = M
        else:
            epsilon[:,0] = E
            E = M+e*np.sin(E)
            epsilon[:,1] = E
    # true anomaly
    cos_nu = np.cos(E)-e
    sin_nu = np.sqrt(1.0-e**2.0)*np.sin(E)
    nu = np.arctan2(sin_nu,cos_nu)
    # argument of latitude
    Phi = nu+omega
    # argument of latitude correction
    delta_u = C_uc*np.cos(2.0*Phi)+C_us*np.sin(2.0*Phi)
    # radius correction
    delta_r = C_rc*np.cos(2.0*Phi)+C_rs*np.sin(2.0*Phi)
    # inclination correction
    delta_inc = C_ic*np.cos(2.0*Phi)+C_is*np.sin(2.0*Phi)
    # corrected argument of latitude
    u = Phi+delta_u
    # corrected radius
    r = a*(1.0-e*np.cos(E))+delta_r
    # corrected inclination
    inc = inc_0+inc_dot*t+delta_inc
    # corrected longitude of ascending node
    Omega = Omega_0+(Omega_dot-omega_e)*t-omega_e*t[0]
    # position in the orbital plane
    X_op = r*np.cos(u)
    Y_op = r*np.sin(u)
    # compute satellite coordinates (ECEF)
    X_ECEF = X_op*np.cos(Omega)-Y_op*np.sin(Omega)*np.cos(inc)
    Y_ECEF = X_op*np.sin(Omega)+Y_op*np.cos(Omega)*np.cos(inc)
    Z_ECEF = Y_op*np.sin(inc)
    # compute satellite coordinates (ECSF)
    X_ECSF = X_op*np.cos(Omega_0)-Y_op*np.sin(Omega_0)*np.cos(inc)
    Y_ECSF = X_op*np.sin(Omega_0)+Y_op*np.cos(Omega_0)*np.cos(inc)
    Z_ECSF = Z_ECEF
    # output
    if ref_sys.lower()=='ecef':
        X = X_ECEF
        Y = Y_ECEF
        Z = Z_ECEF
    elif ref_sys.lower()=='ecsf':
        X = X_ECSF
        Y = Y_ECSF
        Z = Z_ECSF
    return X, Y, Z
    
    
    
    
    
    