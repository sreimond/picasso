# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
from picasso.utils.analysis import special_functions
import numpy as np

def test_legendre1_degree2_order1_x02():
	"""
	Test if the fully normalized associated Legendre function of the first 
	kind of degree 2 and order 1 for the argument x = 0.2 is equal to the 
	analytical expression (within 12 digits accuracy).
	
	Reference: Hofmann-Wellenhof, Moritz 2006 (ISBN 978-3-211-33545-1)
	"""
	p = special_functions.legendre1(2,0.2)
	p21_test = p[0,4]
	p21_true = np.sqrt(15.0) * np.sin(0.2) * np.cos(0.2)
	assert p21_test == pytest.approx(p21_true,1e-12)
