# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
from picasso.utils.analysis import orbits
import numpy as np
    
def test_compute_kepler_orbit():
	"""
	Test if the function `compute_kepler_orbit` correctly computes the 
	orbital coordinates from Keplerian elements.
	"""
	a = 6775000.0
	e = 0.0015306
	inc = 89.0239
	M = 249.6038
	omega = 106.0444
	e = 0.0011067
	inc = 51.6276
	Omega = 176.0525        
	t = np.arange(0,25,5,dtype=float)
	x_true = [-6.7183914059392540e+06,
			  -6.7225791726657506e+06,
			  -6.7265679502936183e+06,
			  -6.7303576233946532e+06,
			  -6.7339480824924260e+06]
	y_true = [7.9234276181779266e+05,
			  7.7130436140848417e+05,
			  7.5024434625745448e+05,
			  7.2916324259044463e+05,
			  7.0806157724666782e+05]
	z_true = [-4.1418331512639782e+05,
			  -3.8421008459302568e+05,
			  -3.5422455643525091e+05,
			  -3.2422769021651061e+05,
			  -2.9422044589666737e+05]
	x_test, y_test, z_test = orbits.compute_kepler_orbit(a,
														 e,
														 inc,
														 omega,
														 Omega,
														 M,                                                             
														 t,
														 Omega_dot=2.5e-07,
														 GM=3.9860044180e+14,
														 omega_e=7.2921150e-05)
	assert x_test[0] == pytest.approx(x_true[0],1e-12)
	assert x_test[1] == pytest.approx(x_true[1],1e-12)
	assert x_test[2] == pytest.approx(x_true[2],1e-12)
	assert x_test[3] == pytest.approx(x_true[3],1e-12)
	assert x_test[4] == pytest.approx(x_true[4],1e-12)
	assert y_test[0] == pytest.approx(y_true[0],1e-12)
	assert y_test[1] == pytest.approx(y_true[1],1e-12)
	assert y_test[2] == pytest.approx(y_true[2],1e-12)
	assert y_test[3] == pytest.approx(y_true[3],1e-12)
	assert y_test[4] == pytest.approx(y_true[4],1e-12)
	assert z_test[0] == pytest.approx(z_true[0],1e-12)
	assert z_test[1] == pytest.approx(z_true[1],1e-12)
	assert z_test[2] == pytest.approx(z_true[2],1e-12)
	assert z_test[3] == pytest.approx(z_true[3],1e-12)
	assert z_test[4] == pytest.approx(z_true[4],1e-12)
        
