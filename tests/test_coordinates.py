# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
from picasso.utils.geometry import coordinates


def test_geocentric2geodetic():
	"""
	Test if the function `geocentric2geodetic` correctly transforms 
	cartesian coordinates into geodetic coordinates.
	
	Reference: Hofmann-Wellenhof, Moritz 2006 (ISBN 978-3-211-33545-1)
			   MATLAB
			   You (2000)
	"""
	x = [4.2010e+06,-2259148.993]
	y = [1.7246e+05,3912960.837]
	z = [4.7801e+06,4488055.516]
	lon_true =  [2.3508,120.000]
	lat_true = [48.8562,45.000]
	h_true = [67.3700,1000.000]
	a_wgs84 = 6378137.0
	b_wgs84 = 6378137.0*(1.0-(1.0/298.257223563))
	lon_test,lat_test,h_test = coordinates.geocentric2geodetic( x, y, z,
																a=a_wgs84,
																b=b_wgs84 )
	assert lon_test[0] == pytest.approx(lon_true[0],1e-4)
	assert lat_test[0] == pytest.approx(lat_true[0],1e-4)
	assert h_test[0] == pytest.approx(h_true[0],1e-4)
	assert lon_test[1] == pytest.approx(lon_true[1],1e-3)
	assert lat_test[1] == pytest.approx(lat_true[1],1e-3)
	assert h_test[1] == pytest.approx(h_true[1],1e-3)
	
def test_geodetic2geocentric():
	"""
	Test if the function `geodetic2geocentric` correctly transforms 
	geodetic coordinates into geocentric coordinates.
	
	Reference: Hofmann-Wellenhof, Moritz 2006 (ISBN 978-3-211-33545-1)
			   You (2000)
	"""
	x_true = [-2265866.507,-2259148.993]
	y_true = [3924595.914,3912960.837]
	z_true = [4501490.544,4488055.516]
	lon = [120.0,120.0]
	lat = [45.0,45.0]
	h = [20000.0,1000.0]
	a_wgs84 = 6378137.0
	b_wgs84 = 6378137.0*(1.0-(1.0/298.257223563))
	x_test, y_test, z_test = coordinates.geodetic2geocentric( lon, lat, h,
															  a=a_wgs84,
															  b=b_wgs84 )
	assert x_test[0] == pytest.approx(x_true[0],1e-3)
	assert y_test[0] == pytest.approx(y_true[0],1e-3)
	assert z_test[0] == pytest.approx(z_true[0],1e-3)
	assert x_test[1] == pytest.approx(x_true[1],1e-3)
	assert y_test[1] == pytest.approx(y_true[1],1e-3)
	assert z_test[1] == pytest.approx(z_true[1],1e-3)
	
