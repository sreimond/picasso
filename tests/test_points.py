# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
from picasso.utils.geometry.points import Point2D

def test_init_distance():
	"""
	Tests the basic and distance methods of the Point2D class.
	"""
	p1 = Point2D(1,5)
	p2 = Point2D(-2,2)
	d_true = 4.242640687119286 
	d_test = p1.calculate_distance(p2)
	assert d_test == pytest.approx(d_true,1e-12)

def test_left_of():
	"""
	Test the `is_left_of` method.
	"""
	p1 = Point2D(-2,2)
	p2 = Point2D(2,-2)
	p = Point2D(0,2)
	assert p.is_left_of(p1,p2)
	p = Point2D(0,-2)
	assert not p.is_left_of(p1,p2)
	