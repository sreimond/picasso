# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:32:15 2018

@author: sreimond
"""

from ..constants import defaults
import numpy as np


def lambert_azimuthal_equal_area( lon, lat,
                                  central_point=None,
                                  R=defaults.EARTH_MEAN_RADIUS()):                                          
    lon = np.array(lon) * np.pi/180.0
    lat = np.array(lat) * np.pi/180.0
    if not central_point:
        central_point = (np.mean(lon),np.mean(lat))
    lon_0 = central_point[0]
    lat_0 = central_point[1]
    slat1 = np.sin(lat_0)
    clat1 = np.cos(lat_0)
    slat = np.sin(lat)
    clat = np.cos(lat)
    sdlon = np.sin(lon-lon_0)
    cdlon = np.cos(lon-lon_0)
    k = (2.0 / (1.0 + slat1 * slat + clat1 * clat * cdlon ))**0.5
    u = R * k * clat * sdlon
    v = R * k * ( clat1 * slat - slat1 * clat * cdlon )
    return u, v
    
def lambert_azimuthal_equal_area_inverse( x, y, central_point,
                                          R=defaults.EARTH_MEAN_RADIUS()):                                    
    x = np.array(x) 
    y = np.array(y)
    lon_0 = central_point[0] * np.pi/180.0
    lat_0 = central_point[1] * np.pi/180.0
    rho = (x**2.0+y**2.0)**0.5
    c = 2.0*(np.arcsin(rho/(2.0*R)))
    cc = np.cos(c)
    sc = np.sin(c)
    slat1 = np.sin( lat_0 )
    clat1 = np.cos( lat_0 )
    lat = np.arcsin(cc*slat1+y*sc*clat1/rho)
    lon = lon_0 + np.arctan2(x*sc,rho*clat1*cc-y*slat1*sc)
    lon *= 180.0/np.pi
    lat *= 180.0/np.pi
    return lon, lat
    
def utm( lon, lat,
         zone=None,
         a=defaults.EARTH_EQUATORIAL_RADIUS(),
         b=defaults.EARTH_POLAR_RADIUS(),         
         hemisphere='N'):            
    if not zone:
        zone=np.floor(((np.floor(np.mean(lon))+180)/6.0)+1.0) 
    lon = np.array(lon) * np.pi/180.0
    lat = np.array(lat) * np.pi/180.0
    # auxiliary quantities
    eprime = np.sqrt((a**2.0-b**2.0)/b**2.0)
    eprime_sq = eprime**2.0
    n = (a-b)/(a+b)
    alpha = ((a+b)*(1.0 + n**2.0/4.0 + n**4.0/64.0))/2.0
    beta = -3.0*n/2.0 + 9.0*n**3.0/16.0-3.0*n**5.0/32.0
    gamma = 15*n**2.0/16.0 - 15.0*n**4.0/32.0
    delta = -35.0*n**3.0/48.0 + 105.0*n**5.0/256.0
    epsilon = 315.0/512.0*n**4.0 
    # reference meridian
    l0 = (-183.0+(zone*6.0))*np.pi/180.0
    # base latitude
    Bphi = (alpha * (lat + beta * np.sin(2.0*lat) + gamma * np.sin(4.0*lat) +
            delta * np.sin(6.0*lat) + epsilon * np.sin(8.0*lat)))
    # auxiliary quantities
    eta = eprime_sq * np.cos(lat)**2.0
    N = (a**2.0)/(b*np.sqrt(1.0+eta))
    t = np.tan(lat)
    # difference to reference meridian
    l = lon - l0
    # projection in east-direction
    u = ((N*np.cos(lat)*l + (1.0/6.0)*N*np.cos(lat)**3.0 * (1.0 - t**2.0 + eta)*l**3.0 + 
         (1.0/120.0)*N*np.cos(lat)**5.0 * (5.0 - 18.0*t**2.0 + t**4.0 + 14.0*eta - 58.0*t**2.0*eta)*l**5.0 + 
         (1.0/5040.0)*N*np.cos(lat)**7.0 * (61.0 - 479.0*t**2.0 + 179.0*t**4.0 - t**6.0)*l**7.0)*0.9996)
    u = u + 500000;
    # projection in north-direction
    v = ((Bphi + (t/2.0) * N * np.cos(lat)**2.0 * l**2.0 + 
         (t/24.0) * N * np.cos(lat)**4.0 * (5.0 - t**2.0 + 9.0*eta + 4.0*eta**2.0)*l**4.0 + 
         (t/720.0) * N * np.cos(lat)**6.0 * (61.0 - 58.0*t**2.0 + t**4.0 + 270.0*eta - 330.0*eta*t**2.0)*l**6.0 + 
         (t/40320.0) * N * np.cos(lat)**8.0 * (1385.0 - 3111.0*t**2.0 + 543.0*t**4.0 - t**2.0)*l**8.0)*0.9996)
    # check if negative
    if hemisphere.lower()=='s':
        v = v + 10000000;
    return u, v, zone
    
def utm_inverse( x, y, zone,
                 a=defaults.EARTH_EQUATORIAL_RADIUS(),
                 b=defaults.EARTH_POLAR_RADIUS()):
    # shift and rescale                                 
    x = np.array(x) 
    y = np.array(y)
    x = (x-500000.0)/0.9996
    y = y/0.9996
    # reference meridian
    l0 = (-183.0+(zone*6.0))*np.pi/180.0
    # auxiliary quantities
    eprime = np.sqrt((a**2.0-b**2.0)/b**2.0)
    eprime_sq = eprime**2.0
    n = (a-b)/(a+b)
    alpha = (a+b)/2.0*(1.0+1.0/4.0*n**2.0+1.0/64.0*n**4.0)
    beta = 3.0/2.0*n-27.0/32.0*n**3.0+269.0/512.0*n**5.0
    gamma = 21.0/16.0*n**2.0-55.0/32.0*n**4.0
    delta = 151.0/96.0*n**3.0-417.0/128.0*n**5.0
    epsilon = 1097.0/512.0*n**4.0
    # footpoint latitude
    y_dash = y/alpha
    phi_f = (y_dash + beta * np.sin(2.0*y_dash) + gamma * np.sin(4.0*y_dash) + 
             delta * np.sin(6.0*y_dash) + epsilon * np.sin(8.0*y_dash))
    # auxiliary quantities
    eta_f = np.sqrt(eprime_sq*np.cos(phi_f)**2.0)
    N_f = a**2.0/(b*np.sqrt(1.0+eta_f**2.0))
    t_f = np.tan(phi_f)
    # conversion to latitude
    phi = (phi_f + t_f/(2.0*N_f**2.0) * (-1.0-eta_f**2.0) * x**2.0 + 
           t_f/(24.0*N_f**4.0) * (5.0+3.0*t_f**2.0+6.0*eta_f**2.0-6.0*t_f**2.0*eta_f**2.0-3.0*eta_f**4.0-9.0*t_f**2.0*eta_f**4.0) * x**4.0 + 
           t_f/(760.0*N_f**6.0) * (-61.0-90.0-t_f**2.0-107.0*eta_f**2.0+162.0*t_f**2.0*eta_f**2.0+45.0*t_f**4.0*eta_f**2.0) * x**6.0 + 
           t_f/(40320.0*N_f**8.0) * (1385.0+3633.0*t_f**2.0+4095.0*t_f**4.0+1575.0*t_f**6.0) * x**8.0)
    # conversion to longitude
    lam = (l0 + 1.0/(N_f*np.cos(phi_f)) * x + 
           1.0/(6.0*N_f**3.0*np.cos(phi_f)) * (-1.0-2.0*t_f**2.0-eta_f**2.0) * x**3.0 + 
           1.0/(120.0*N_f**5.0*np.cos(phi_f)) * (5.0+28.0*t_f**2.0+24.0*t_f**4.0+6.0*eta_f**2.0+8.0*t_f**2.0*eta_f**2.0) * x**5.0 + 
           1.0/(5040.0*N_f**7.0*np.cos(phi_f)) * (-61.0-662.0*t_f**2.0-1320.0*t_f**4.0-720.0*t_f**6.0) * x**7.0)
    # rad to deg
    lon = lam * 180.0/np.pi
    lat = phi * 180.0/np.pi
    return lon, lat
    

def azimuthal_equidistant( lon, lat,
                            central_point=None,
                            R=defaults.EARTH_MEAN_RADIUS()):
    lon = np.array(lon) * np.pi/180.0
    lat = np.array(lat) * np.pi/180.0
    if not central_point:
        central_point = (np.mean(lon),np.mean(lat))
    lon_0 = central_point[0]
    lat_0 = central_point[1]
    slat1 = np.sin(lat_0)
    clat1 = np.cos(lat_0)
    slat = np.sin(lat)
    clat = np.cos(lat)
    sdlon = np.sin(lon-lon_0)
    cdlon = np.cos(lon-lon_0)
    c = np.arccos( slat1 * slat + clat1 * clat * cdlon )
    k_prime = c/np.sin(c)
    u = R * k_prime * clat * sdlon
    v = R * k_prime * ( clat1 * slat - slat1 * clat * cdlon )
    return u, v
    
def azimuthal_equidistant_inverse( x, y, central_point,
                                          R=defaults.EARTH_MEAN_RADIUS()):                                    
    x = np.array(x) 
    y = np.array(y)
    lon_0 = central_point[0] * np.pi/180.0
    lat_0 = central_point[1] * np.pi/180.0
    rho = (x**2.0+y**2.0)**0.5
    c = rho/R
    cc = np.cos(c)
    sc = np.sin(c)
    slat1 = np.sin( lat_0 )
    clat1 = np.cos( lat_0 )
    lat = np.arcsin(cc*slat1+y*sc*clat1/rho)
    lat *= 180.0/np.pi
    lon = lon_0 + np.arctan2(x*sc,rho*clat1*cc-y*slat1*sc)    
    ix90 = lat==90
    ixm90 = lat==-90    
    lon[ix90] = lon_0 + np.arctan2(-x[ix90],y[ix90])
    lon[ixm90] = lon_0 + np.arctan2(x[ix90],y[ix90])
    lon *= 180.0/np.pi
    return lon, lat
    