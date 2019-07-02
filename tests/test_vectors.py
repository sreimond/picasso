# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
from picasso.utils.algebra import vectors
import numpy as np


def test_init_add_count():
	"""
	Test general functionality of the Vector class.
	"""
	v = vectors.Vector()
	v.add_element(3)
	v.add_element(5)
	v.add_element(-10)
	c = v.count_elements()
	assert c == 3

def test_cross():
	"""
	Test the cross product method.
	Reference: 
	https://math.oregonstate.edu/home/programs/undergrad/CalculusQuestStudyGuides/vcalc/crossprod/crossprod.html
	"""
	u = vectors.Vector()
	u.add_element(3)
	u.add_element(-2)
	u.add_element(-2)
	v = vectors.Vector()
	v.add_element(-1)
	v.add_element(0)
	v.add_element(5)
	c_true = [-10,-13,-2]
	c_test = u.cross(v).tolist()
	assert c_test == c_true
	
def test_dot():
	"""
	Test the dot product method.
	"""
	u = vectors.Vector()
	u.add_element(3)
	u.add_element(-2)
	u.add_element(-2)
	v = vectors.Vector()
	v.add_element(-1)
	v.add_element(0)
	v.add_element(5)
	c_true = -13
	c_test = u.dot(v)
	assert c_test == c_true

def test_norm():
	"""
	Test the norm method.
	"""
	u = vectors.Vector()
	u.add_element(3.0)
	u.add_element(7.9)
	u.add_element(5.2)
	u.add_element(-2.1)
	u.add_element(-4.4)
	n_true = 11.055315463612967
	n_test = u.norm()
	assert n_test == pytest.approx(n_true,1e-12)

def test_angle_between():
	"""
	Test the angles between method.
	Reference: http://onlinemschool.com/math/library/vector/angl/
	"""
	u = vectors.Vector()
	u.add_element(3)
	u.add_element(4)
	v = vectors.Vector()
	v.add_element(4)
	v.add_element(3)
	a_true = np.arccos(0.96)
	a_test = u.angle_between(v)
	assert a_test == pytest.approx(a_true,1e-12)
	  