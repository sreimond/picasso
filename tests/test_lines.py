# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
from picasso.utils.geometry.points import Point2D
from picasso.utils.geometry.lines import Line2D
from picasso.utils.algebra.vectors import Vector

def test_init_add_vector():
	"""
	Tests the basic methods of the Line2D class.
	"""
	p1 = Point2D(1,5)
	p2 = Point2D(-2,2)
	l = Line2D(point1=p1)
	l.add_endpoint(p2)
	v = l.as_vector()
	assert isinstance(v,Vector)
        
def test_length_epsilon():
	p1 = Point2D(-4,2)
	p2 = Point2D(6,2)
	l = Line2D()
	l.add_endpoint(p1)
	l.add_endpoint(p2)
	d_true = 10
	d_test = l.calculate_length()
	assert d_test == d_true
	e_true = 1e-4
	e_test = l.calculate_epsilon()
	assert e_test == e_true
        
def test_point_distance():
	p1 = Point2D(-4,2)
	p2 = Point2D(6,2)
	p3 = Point2D(0,5.75)
	l = Line2D()
	l.add_endpoint(p1)
	l.add_endpoint(p2)
	d_true = 3.75
	d_test = l.calculate_point_distance(p3)
	assert d_test == d_true

def test_point_location():
	p1 = Point2D(-4,2)
	p2 = Point2D(6,2)
	l = Line2D()
	l.add_endpoint(p1)
	l.add_endpoint(p2)        
	p = Point2D(0,5.75)
	assert l.determine_point_location(p) == 0
	p = Point2D(0,2)
	assert l.determine_point_location(p) == 1
	p = Point2D(-100,2)
	assert l.determine_point_location(p) == 2
    
def test_line_parallel():
	p1 = Point2D(-4,2)
	p2 = Point2D(6,2)
	l = Line2D()
	l.add_endpoint(p1)
	l.add_endpoint(p2) 
	l2 = Line2D(Point2D(5,10),Point2D(-10,10))
	assert l.determine_line_parallel(l2)
	l2 = Line2D(Point2D(0,0),Point2D(10,10))
	assert not l.determine_line_parallel(l2)
    
def test_line_intersection_point():
	p1 = Point2D(-4,2)
	p2 = Point2D(6,2)
	l = Line2D()
	l.add_endpoint(p1)
	l.add_endpoint(p2) 
	l2 = Line2D(Point2D(2.45,-5),Point2D(2.45,10))
	p_true = Point2D(2.45,2)
	p_test = l.calculate_line_intersection_point(l2)
	assert p_test.x == p_true.x
	assert p_test.y == p_true.y
        
def test_line_intersection():
	p1 = Point2D(-4,2)
	p2 = Point2D(6,2)
	l = Line2D()
	l.add_endpoint(p1)
	l.add_endpoint(p2) 
	l2 = Line2D(Point2D(-2,-2),Point2D(2,3))
	assert l.determine_line_intersection(l2) == 1
	l2 = Line2D(Point2D(7,-2),Point2D(2,-3))
	assert l.determine_line_intersection(l2) == 2
	l2 = Line2D(Point2D(5,10),Point2D(-10,10))
	assert l.determine_line_parallel(l2)

