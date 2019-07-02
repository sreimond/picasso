# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
import numpy as np
from picasso.utils.geometry.polygons import Polygon, MultiPolygon
from picasso.utils.geometry.points import Point2D
from picasso.utils.geometry.lines import Line2D


def test_triangle():
	P1 = Point2D(5,0)
	P2 = Point2D(0,5)
	P3 = Point2D(-5,0)
	points = (P1,P2,P3)
	a = 10.0
	b = np.sqrt(5.0**2.0+5.0**2.0)
	c = np.sqrt(5.0**2.0+5.0**2.0)
	triangle = Polygon(points)
	# test perimeter
	p_true = a+b+c
	p_test = triangle.calculate_perimeter()
	assert p_test == pytest.approx(p_true,1e-12)
	# test area
	s = (a+b+c)/2.0
	A_true = np.sqrt(s*(s-a)*(s-b)*(s-c))
	A_test = triangle.calculate_area()
	assert A_test == pytest.approx(A_true,1e-12)
	# test centroid
	x_true = (5.0+0-5.0)/3.0
	y_true = (0+5.0+0)/3.0
	c = triangle.calculate_centroid()
	x_test = c.x
	y_test = c.y
	assert x_test == pytest.approx(x_true,1e-12)
	assert y_test == pytest.approx(y_true,1e-12)
	# test orientation
	o_true = 'ccw'
	o_test = triangle.determine_orientation()
	assert o_test == o_true
	# test convexity
	c_true = True
	c_test = triangle.determine_convexity()
	assert c_test == c_true
	# test complexity
	c_true = False
	c_test = triangle.determine_complexity()
	assert c_test == c_true
	# tests for bounding box
	bbox = triangle.determine_bounding_box()
	p_true = 5.0*2.0+10.0*2.0
	p_test = bbox.calculate_perimeter()
	assert p_test == pytest.approx(p_true,1e-12)
	A_true = 5.0*10.0
	A_test = bbox.calculate_area()
	assert A_test == pytest.approx(A_true,1e-12)
	x_true = (5.0-5.0)/2.0
	y_true = (0+5.0)/2.0
	c = bbox.calculate_centroid()
	x_test = c.x
	y_test = c.y
	assert x_test == pytest.approx(x_true,1e-12)
	assert y_test == pytest.approx(y_true,1e-12)
	o_true = 'ccw'
	o_test = bbox.determine_orientation()
	assert o_test == o_true
	c_true = True
	c_test = bbox.determine_convexity()
	assert c_test == c_true
	c_true = False
	c_test = bbox.determine_complexity()
	assert c_test == c_true
	# test points in polygon
	point = Point2D(0,1)
	assert triangle.determine_point_location(point) == 1
	point = Point2D(-4,0.5)
	assert triangle.determine_point_location(point) == 1
	point = Point2D(4,1)
	assert triangle.determine_point_location(point) == 2
	point = Point2D(0,0)
	assert triangle.determine_point_location(point) == 2
	point = Point2D(-1,-1)
	assert triangle.determine_point_location(point) == 0

def test_concave_polygon():
	polygon = Polygon()
	polygon.add_vertex(Point2D(10,-10))
	polygon.add_vertex(Point2D(5,0))
	polygon.add_vertex(Point2D(10,10))
	polygon.add_vertex(Point2D(-10,10))
	polygon.add_vertex(Point2D(-5,0))
	polygon.add_vertex(Point2D(-10,-10))
	# test convexity
	c_true = False
	c_test = polygon.determine_convexity()
	assert c_test == c_true
	# test complexity
	c_true = False
	c_test = polygon.determine_complexity()
	assert c_test == c_true
	# test points in polygon
	point = Point2D(0,1)
	assert polygon.determine_point_location(point) == 1
	point = Point2D(-8,-10)
	assert polygon.determine_point_location(point) == 2
	point = Point2D(10,10)
	assert polygon.determine_point_location(point) == 2
	point = Point2D(10,0)
	assert polygon.determine_point_location(point) == 0
	
def test_complex_polygon():
	polygon = Polygon()
	polygon.add_vertex(Point2D(10,0))
	polygon.add_vertex(Point2D(0,10))
	polygon.add_vertex(Point2D(0,-10))
	polygon.add_vertex(Point2D(-10,10))
	polygon.add_vertex(Point2D(-10,0))
	# test convexity
	c_true = False
	c_test = polygon.determine_convexity()
	assert c_test == c_true
	# test complexity
	c_true = True
	c_test = polygon.determine_complexity()
	assert c_test == c_true
	# test points in polygon
	point = Point2D(5,2)
	assert polygon.determine_point_location(point) == 1
	point = Point2D(-1,-2)
	assert polygon.determine_point_location(point) == 1
	point = Point2D(-8,3)
	assert polygon.determine_point_location(point) == 1
	point = Point2D(-2,3)
	assert polygon.determine_point_location(point) == 0
	point = Point2D(-5,0)
	assert polygon.determine_point_location(point) == 2
	
def test_multipolygon():
	p1 = Polygon()
	p1.add_vertex(Point2D(-1,-1))
	p1.add_vertex(Point2D(-2,1))
	p1.add_vertex(Point2D(-1,5))
	p1.add_vertex(Point2D(-4,3))
	p1.add_vertex(Point2D(-5,-3))
	p1.add_vertex(Point2D(-3,0))
	p2 = Polygon()
	p2.add_vertex(Point2D(5,5))
	p2.add_vertex(Point2D(6,3))
	p2.add_vertex(Point2D(3,-3))
	p2.add_vertex(Point2D(2,2))        
	multi_polygon = MultiPolygon()
	multi_polygon.add_polygon(p1)
	multi_polygon.add_polygon(p2)        
	# test points in polygon
	point = Point2D(5,2)
	assert multi_polygon.determine_point_location(point) == 1
	point = Point2D(-4,-1)
	assert multi_polygon.determine_point_location(point) == 1
	point = Point2D(-2,0)
	assert multi_polygon.determine_point_location(point) == 1
	point = Point2D(-2,1)
	assert multi_polygon.determine_point_location(point) == 2
	point = Point2D(0,0)
	assert multi_polygon.determine_point_location(point) == 0
	point = Point2D(4,-3)
	assert multi_polygon.determine_point_location(point) == 0
	point = Point2D(-6,1)
	assert multi_polygon.determine_point_location(point) == 0
