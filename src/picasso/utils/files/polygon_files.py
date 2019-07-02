# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 11:23:44 2018

@author: sreimond
"""

from ..geometry.polygons import Polygon, MultiPolygon
from ..geometry.points import Point2D
import numpy as np
import os
import warnings

def read_polygon_file( input_filename ):
    """
    The function `read_polygon_file` inside the polygon_files module reads a 
    polygon in the txt file format and returns an instane of the MultiPolygon 
    class.
    
    The polygon file must consist of three columns (without header):
        1) Polygon ID (non-neagtive IDs only)
        2) x coordinate
        3) y coordinate
    """
    if not os.path.exists( input_filename ):
        warnings.warn('Cannot open TXT file!')                
        return 
    ix_polygon, x_polygon, y_polygon = np.genfromtxt( input_filename, unpack=True )
    ix_polygon = np.array(ix_polygon,dtype=int)
    multi_polygon = MultiPolygon()
    [multi_polygon.add_polygon(Polygon()) for _ in range(max(ix_polygon)+1)]
    for ix,x,y in zip(ix_polygon, x_polygon, y_polygon):
        multi_polygon.polygons[ix].add_vertex(Point2D(x,y))
    return multi_polygon

def export_polygon( polygon, output_filename ):
    """
    The function `export_polygon` inside the polygon_files module lets you 
    export an instance of the Polygon OR MultiPolygon class as a text file.
    """    
    if os.path.exists( output_filename ):
        warnings.warn('Existing output file will be overwritten!')
    f = open( output_filename, 'w' )
    if isinstance(polygon,Polygon):
        for vertex in polygon.vertices:
            f.write('%d\t%.20e\t%.20e\n' % (0,vertex.x,vertex.y))
    elif isinstance(polygon,MultiPolygon):
        for ix,polygon in enumerate(polygon.polygons):
            for vertex in polygon.vertices:
                f.write('%d\t%.20e\t%.20e\n' % (ix,vertex.x,vertex.y))
    f.close()

def convert_polygon_file_txt2xml( input_filename, output_filename=None ):
    """
    The function `convert_polygon_file_txt2xml` inside the polygon_files module 
    converts a polygon file from the txt to the GROOPS xml format.
    """    
    if not os.path.exists( input_filename ):
        warnings.warn('Cannot open TXT file!')                
        return
    if not output_filename:
        output_filename = input_filename[0:-4] + '.xml'
    if os.path.exists( output_filename ):
        warnings.warn('Existing output file will be overwritten!')
    polygons, x, y = np.genfromtxt( input_filename, unpack=True )
    polygon_count = int(polygons[-1] + 1)
    f = open(output_filename,'w')
    f.write('<?xml version="1.0" encoding="UTF-8" ?>\n')
    f.write('<groops>\n')
    f.write('\t<polygonList>\n')
    f.write('\t\t<polygonCount>%d</polygonCount>\n' % polygon_count)
    for i in range(polygon_count):
        ix = np.where(polygons==i)
        xi = x[ix]
        yi = y[ix]
        xi = np.append(xi,xi[-1])
        yi = np.append(yi,yi[-1])
        point_count = len(xi)
        f.write('\t\t<polygon>\n')        
        f.write('\t\t\t<pointCount>%d</pointCount>\n' % point_count)
        for j in range(point_count):
            f.write('\t\t\t<point>\n')
            f.write('\t\t\t\t<longitude>%.20e</longitude>\n' % xi[j])
            f.write('\t\t\t\t<latitude>%.20e</latitude>\n' % yi[j])
            f.write('\t\t\t</point>\n')
        f.write('\t\t</polygon>\n')    
    f.write('\t</polygonList>\n')
    f.write('</groops>\n')
    f.close()

def convert_polygon_file_xml2txt( input_filename, output_filename=None ):
    """
    The function `convert_polygon_file_xml2txt` inside the polygon_files module 
    converts a polygon file from the GROOPS xml format to the txt format.
    """    
    if not os.path.exists( input_filename ):
        warnings.warn('Cannot open XML file!')                
        return
    if not output_filename:
        output_filename = input_filename[0:-4] + '.txt'
    if os.path.exists( output_filename ):
        warnings.warn('Existing output file will be overwritten!')
    fw = open(output_filename,'w')
    fr = open(input_filename,'r')
    ix = -1
    for line in fr:
        if '<polygon>' in line:
            ix += 1
        if '<longitude>' in line:
            xi = float(line[line.index("<longitude>") + 
                       11:line.rindex("</longitude>")])
        if '<latitude>' in line:
            yi = float(line[line.index("<latitude>") + 
                       10:line.rindex("</latitude>")])
        if '</point>' in line:
            fw.write('%d\t%.20e\t%.20e\n' % (ix,xi,yi))
    fr.close()
    fw.close()




    