# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
import os
from picasso.utils.files.polygon_files import read_polygon_file, export_polygon
from picasso.utils.files.polygon_files import convert_polygon_file_txt2xml
from picasso.utils.files.polygon_files import convert_polygon_file_xml2txt
from picasso.utils.geometry.polygons import Polygon, MultiPolygon
from picasso.utils.geometry.points import Point2D


def test_read_polygon_file(tmp_path):
	test_multi_polygon = ((0,1,2),
						  (0,4,5),
						  (0,7,8),
						  (1,-1,-2),
						  (1,-4,-5),
						  (1,-7,-8))
	f = open(os.path.join(tmp_path , 'test.txt'), 'w')
	for data in test_multi_polygon:
		f.write('%d\t%.20e\t%.20e\n' % (data[0],data[1],data[2]))
	f.close()
	multi_polygon = read_polygon_file(os.path.join(tmp_path, 'test.txt'))
	assert multi_polygon.polygons[0].vertices[0].x == 1
	assert multi_polygon.polygons[0].vertices[0].y == 2
	assert multi_polygon.polygons[0].vertices[1].x == 4
	assert multi_polygon.polygons[0].vertices[1].y == 5
	assert multi_polygon.polygons[0].vertices[2].x == 7
	assert multi_polygon.polygons[0].vertices[2].y == 8
	assert multi_polygon.polygons[1].vertices[0].x == -1
	assert multi_polygon.polygons[1].vertices[0].y == -2
	assert multi_polygon.polygons[1].vertices[1].x == -4
	assert multi_polygon.polygons[1].vertices[1].y == -5
	assert multi_polygon.polygons[1].vertices[2].x == -7
	assert multi_polygon.polygons[1].vertices[2].y == -8
        
def test_export_polygon(tmp_path):
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
	export_polygon(multi_polygon,os.path.join(tmp_path, 'test.txt'))
	multi_polygon_test = read_polygon_file(os.path.join(tmp_path, 'test.txt'))
	assert multi_polygon_test.polygons[0].vertices[0].x == -1
	assert multi_polygon_test.polygons[1].vertices[2].y == -3
        
def test_file_conversions(tmp_path):
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
	export_polygon(multi_polygon,os.path.join(tmp_path, 'test.txt'))
	convert_polygon_file_txt2xml(os.path.join(tmp_path, 'test.txt'))
	convert_polygon_file_xml2txt(os.path.join(tmp_path, 'test.xml'),
								 os.path.join(tmp_path, 'conv.txt'))
	multi_polygon_test = read_polygon_file(os.path.join(tmp_path, 'conv.txt'))
	assert multi_polygon_test.polygons[0].vertices[0].x == -1
	assert multi_polygon_test.polygons[1].vertices[2].y == -3
      
	  