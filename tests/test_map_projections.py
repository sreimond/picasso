# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
import numpy as np
from picasso.utils.geometry.points import Point2D
from picasso.utils.geometry import map_projections

def test_lambert_azimuthal_equal_area():
	lon = [15.61,15.22,16.00]
	lat = [22.16,22.90,21.00]
	R = 6371008.0
	u_true = [0,-39949.26325187112,40487.21380985146]
	v_true = [15567.30540523261,97902.15736458405,-113366.3751508453]
	u_test, v_test = map_projections.lambert_azimuthal_equal_area(lon,lat,R=R)
	assert u_test[0] == pytest.approx(u_true[0],1e-9)
	assert u_test[1] == pytest.approx(u_true[1],1e-9)
	assert u_test[2] == pytest.approx(u_true[2],1e-9)
	assert v_test[0] == pytest.approx(v_true[0],1e-9)
	assert v_test[1] == pytest.approx(v_true[1],1e-9)
	assert v_test[2] == pytest.approx(v_true[2],1e-9)
        
def test_lambert_azimuthal_equal_area_inverse():
	u = [1734.0,-35333.0,45123.0]
	v = [14882.0,99948.0,-121255.0]
	central_point = (15.0,22.0)
	R = 6371008.0
	lon_true = [15.016834827398265,14.655070183549958,15.434388348979571]
	lat_true = [22.133836044309135,22.898500257207296,20.908943005438836]
	lon_test, lat_test = map_projections.lambert_azimuthal_equal_area_inverse(u,v,R=R,central_point=central_point)
	assert lon_test[0] == pytest.approx(lon_true[0],1e-9)
	assert lon_test[1] == pytest.approx(lon_true[1],1e-9)
	assert lon_test[2] == pytest.approx(lon_true[2],1e-9)
	assert lat_test[0] == pytest.approx(lat_true[0],1e-9)
	assert lat_test[1] == pytest.approx(lat_true[1],1e-9)
	assert lat_test[2] == pytest.approx(lat_true[2],1e-9)
        
def test_utm():
	lon = [15.61,15.22,16.00]
	lat = [22.16,22.90,21.00]
	a = 6378137.0
	b = a*(1.0-1.0/298.257222101)
	u_true = [562894.6954475424,522562.5488405955,603932.761363141]
	v_true = [2450663.177064531,2532466.717880427,2322472.695370913]
	z_true = 33
	u_test, v_test, z_test = map_projections.utm(lon,lat,a=a,b=b)
	assert u_test[0] == pytest.approx(u_true[0],1e-9)
	assert u_test[1] == pytest.approx(u_true[1],1e-9)
	assert u_test[2] == pytest.approx(u_true[2],1e-9)
	assert v_test[0] == pytest.approx(v_true[0],1e-9)
	assert v_test[1] == pytest.approx(v_true[1],1e-9)
	assert v_test[2] == pytest.approx(v_true[2],1e-9)
	assert z_test == z_true
	
def test_utm_inverse():
	u = [562894.6954475424,522562.5488405955,603932.761363141]
	v = [2450663.177064531,2532466.717880427,2322472.695370913]
	zone = 33
	a = 6378137.0
	b = a*(1.0-1.0/298.257222101)
	lon_true = [15.61,15.22,16.00]
	lat_true = [22.16,22.90,21.00]
	lon_test, lat_test = map_projections.utm_inverse(u,v,zone=zone,a=a,b=b)
	assert lon_test[0] == pytest.approx(lon_true[0],1e-9)
	assert lon_test[1] == pytest.approx(lon_true[1],1e-9)
	assert lon_test[2] == pytest.approx(lon_true[2],1e-9)
	assert lat_test[0] == pytest.approx(lat_true[0],1e-9)
	assert lat_test[1] == pytest.approx(lat_true[1],1e-9)
	assert lat_test[2] == pytest.approx(lat_true[2],1e-9)
	
def test_azimuthal_equidistant():
	lon = [100.0]
	lat = [-20.0]
	R = 3.0
	u_true = [-5.8311398]
	v_true = [5.5444634]
	central_point = (-100.0*np.pi/180.0,40.0*np.pi/180.0)
	u_test, v_test = map_projections.azimuthal_equidistant(lon,lat,R=R,central_point=central_point)        
	assert u_test[0] == pytest.approx(u_true[0],1e-7)
	assert v_test[0] == pytest.approx(v_true[0],1e-7)
	
def test_azimuthal_equidistant_inverse():
	u = [-5.8311398]
	v = [5.5444634]
	central_point = (-100.0,40.0)
	R = 3.0
	lon_true = [-260.0]
	lat_true = [-20.0]
	lon_test, lat_test = map_projections.azimuthal_equidistant_inverse(u,v,R=R,central_point=central_point)
	assert lon_test[0] == pytest.approx(lon_true[0],1e-5)
	assert lat_test[0] == pytest.approx(lat_true[0],1e-5)
