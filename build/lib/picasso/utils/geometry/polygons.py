# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 11:40:02 2018

@author: sreimond
"""

from .points import Point2D
from .lines import Line2D
import numpy as np
import os
import warnings

class Polygon( object ):
    """
    The `Polygon` class allows to represent a collection of 2D points in terms
    of a Polygon object. Methods include area calculation, point in polygon 
    check and convexity analysis.
    """
    def __init__(self,vertices=[]):
        self.vertices = []
        for vertex in vertices:
            if isinstance(vertex,tuple):
                point = Point2D(*vertex)
                self.add_vertex(point)
            elif isinstance(vertex,Point2D):
                self.add_vertex(vertex)
    
    def add_vertex(self,point):
        """
        Add a Point2D object to the list of vertices (if it differs from 
        the previous vertex).
        """     
        if len(self.vertices)==0:
            self.vertices.append((point))
            return
        if (self.vertices[-1].calculate_distance(point) > 0) and (self.vertices[0].calculate_distance(point) > 0):
            self.vertices.append(point)  
        self.edges = []
        self.update_edges()        
    
    def update_edges(self):
        """
        Update list of edges, i.e. collection of Line2D objects.
        """
        if len(self.vertices)<2:
            warnings.warn('This polygon consits of only one vertex!')
            return
        for ix,vertex in enumerate(self.vertices):
            edge = Line2D(self.vertices[ix-1],vertex)
            self.edges.append(edge)
    
    def calculate_perimeter(self):
        """ Calculates the total perimeter based on the list of vertices. """
        perimeter = 0
        for edge in self.edges:
            perimeter += edge.calculate_length()
        return perimeter
        
    def calculate_area(self):
        """ Calculates the total area based on the list of vertices. """
        area = 0
        for ix,vertex in enumerate(self.vertices):
            area += ( (vertex.y*self.vertices[ix-1].x)
                    - (vertex.x*self.vertices[ix-1].y) )
        return 0.5 * area  
    
    def calculate_centroid(self):
        """
        Calculates the centroid point based on the list of vertices and the
        polygon area.
        """
        area = self.calculate_area()
        xc = 0
        yc = 0        
        for ix,vertex in enumerate(self.vertices):
            xc += ((self.vertices[ix-1].x+vertex.x)*(self.vertices[ix-1].x*vertex.y-vertex.x*self.vertices[ix-1].y))
            yc += ((self.vertices[ix-1].y+vertex.y)*(self.vertices[ix-1].x*vertex.y-vertex.x*self.vertices[ix-1].y))
        xc /= (6.0*area)        
        yc /= (6.0*area)
        return Point2D(xc,yc)
        
    def determine_orientation(self):
        """ 
        Determines whether the polygon vertices are ordered clockwise or 
        counterclockwise or if they only form a line. 
        """
        result = 0
        for ix,vertex in enumerate(self.vertices):
            result += ( (vertex.x-self.vertices[ix-1].x)
                      * (vertex.y+self.vertices[ix-1].y) )
        if result > 0:
            return 'cw'
        elif result < 0:
            return 'ccw'
        else:
            return 'line'
    
    def determine_convexity(self):
        """ Determines whether the polygon is convex or not. """
        result = []
        for ix,vertex in enumerate(self.vertices):
            result.append(vertex.is_left_of(self.vertices[ix-2],self.vertices[ix-1]))
        return (np.all(result)) or (not np.any(result)) and (not self.determine_complexity())
    
    def determine_complexity(self):
        """ 
        Determines whether the polygon is complex or not (self intersecting). 
        """        
        for edge_i in self.edges:
            for edge_j in self.edges:
                if edge_i == edge_j:
                    continue
                if edge_i.determine_line_parallel(edge_j):
                    continue
                point = edge_i.calculate_line_intersection_point(edge_j)
                on_i = (edge_i.determine_point_location(point) == 1)
                on_j = (edge_j.determine_point_location(point) == 1)
                if on_i and on_j:
                    di1 = point.calculate_distance(edge_i.endpoints[0])
                    di2 = point.calculate_distance(edge_i.endpoints[1])
                    dj1 = point.calculate_distance(edge_j.endpoints[0])
                    dj2 = point.calculate_distance(edge_j.endpoints[1])
                    if ((di1 > edge_i.calculate_epsilon()) and
                        (di2 > edge_i.calculate_epsilon()) and
                        (dj1 > edge_j.calculate_epsilon()) and
                        (dj2 > edge_j.calculate_epsilon())):
                            return True
        return False
        
    def determine_bounding_box(self):
        """ Returns the bounding box of the polygon as a Polygon instance. """
        min_x = np.min( [vertex.x for vertex in self.vertices] )
        min_y = np.min( [vertex.y for vertex in self.vertices] )
        max_x = np.max( [vertex.x for vertex in self.vertices] )
        max_y = np.max( [vertex.y for vertex in self.vertices] )
        bounding_box = BoundingBox()
        bounding_box.add_vertex(Point2D(min_x,min_y))
        bounding_box.add_vertex(Point2D(max_x,min_y))
        bounding_box.add_vertex(Point2D(max_x,max_y))
        bounding_box.add_vertex(Point2D(min_x,max_y))
        return bounding_box
        
    def determine_point_location(self,point):
        """
        Returns 1 if the point lies inside the polygon, 2 if it is on the 
        boundary and 0 if it is outside the polygon.
        
        This algorithm is described in Galetzka and Glauner 2012.
        """
        # check if inside bounding box        
        if self.determine_bounding_box().determine_point_location(point) == 0:
            return 0
        # if yes, check if point lies on one of the edges
        for edge in self.edges:
            if edge.determine_point_location(point) == 1:
                return 2
        # if not, check if the polygon is convex so that a simple method can be 
        # applied
        if self.determine_convexity():
            left = []
            for edge in self.edges:
                left.append(point.is_left_of(edge.endpoints[0],edge.endpoints[1]))            
            return int((np.all(left)) or (not np.any(left)))        
        # if the polygon is not convex, use the algorithm cited above:        
        vertices = []
        for vertex in self.vertices:
            dx = vertex.x - point.x
            dy = vertex.y - point.y
            vertices.append(Point2D(dx,dy))
#        max_x = np.max( [vertex.x for vertex in vertices] ) * 1e6 + 1.0
        max_x = np.max( [vertex.x for vertex in vertices] ) + 1.0
        min_x = np.min( [vertex.x for vertex in vertices] ) - 1.0
#        max_x = 1e12
        x_axis_pos = Line2D(Point2D(0,0),Point2D(max_x,0))
        x_axis_tot = Line2D(Point2D(min_x,0),Point2D(max_x,0))
        indices_s = np.nonzero( [vertex.y for vertex in vertices] )[0]   
        if indices_s.size == 0:
            return 0
        count_intersections = 0
        for ix,index_si in enumerate(indices_s):
            index_s = indices_s[ix-1]
            point_s = Point2D(vertices[index_s].x,vertices[index_s].y)
            point_si = Point2D(vertices[index_si].x,vertices[index_si].y)
            edge = Line2D(point_s,point_si)
            count_skipped = (index_si - index_s) - 1
            if count_skipped==0 or count_skipped==-len(vertices):
                if edge.determine_line_intersection(x_axis_pos) == 1:
                    count_intersections += 1
            else:
                indices_skipped = []
                if index_s > index_si:
                    indices_skipped += range(0,index_si)
                    indices_skipped += range(index_s+1,len(vertices))
                else:
                    indices_skipped += range(index_s+1,index_si)
                x_skipped = [vertices[i].x for i in indices_skipped]
                count_skipped_x_pos = sum([x > 0 for x in x_skipped])                
                if count_skipped_x_pos > 0:
                    if edge.determine_line_intersection(x_axis_tot) == 1:
                        count_intersections += 1       
        return int((count_intersections % 2) != 0)
        
        
class BoundingBox( Polygon ):
    """
    The `BoundingBox` class is a subclass of the Polygon class and allows 
    faster computations.
    """
    def __init__(self):
        super(BoundingBox, self).__init__()
    
    def determine_point_location(self,point):
        """
        Returns 1 if the point lies inside the polygon, 2 if it is on the 
        boundary and 0 if it is outside the polygon.
        """
        if ( (point.x < self.vertices[0].x) or 
             (point.x > self.vertices[1].x) or
             (point.y < self.vertices[0].y) or
             (point.y > self.vertices[2].y) ):
            return 0
        elif ( (point.x == self.vertices[0].x) and
               (point.y >= self.vertices[0].y) and
               (point.y <= self.vertices[2].y) ):
            return 2
        elif ( (point.x == self.vertices[1].x) and
               (point.y >= self.vertices[0].y) and
               (point.y <= self.vertices[2].y) ):
            return 2
        elif ( (point.y == self.vertices[0].y) and
               (point.x >= self.vertices[0].x) and
               (point.x <= self.vertices[1].x) ):
            return 2
        elif ( (point.y == self.vertices[2].y) and
               (point.x >= self.vertices[0].x) and
               (point.x <= self.vertices[1].x) ):
            return 2
        else:
            return 1
            

class MultiPolygon( object ):
    """
    The `MultiPolygon` class allows to represent a collection of polygons in 
    terms of a MultiPolygon object.
    """
    def __init__(self,polygons=[]):
        self.polygons = []
        for polygon in polygons:
            if isinstance(polygon,Polygon):
                self.add_polygon(polygon)
    
    def add_polygon(self,polygon):
        """
        Add a Polygon object to the list of polygons.
        """     
        self.polygons.append(polygon)
    
    def calculate_area(self):
        """
        Sums up the areas of all polygons.
        """
        A = 0
        for polygon in self.polygons:
            A += polygon.calculate_area()
        return A
    
    def determine_bounding_box(self):
        """ 
        Returns the bounding box of the MultiPolygon as a 
        Polygon instance. 
        """
        min_x = np.min( [vertex.x for poly in self.polygons for vertex in poly.vertices] )
        min_y = np.min( [vertex.y for poly in self.polygons for vertex in poly.vertices] )
        max_x = np.max( [vertex.x for poly in self.polygons for vertex in poly.vertices] )
        max_y = np.max( [vertex.y for poly in self.polygons for vertex in poly.vertices] )
        bounding_box = BoundingBox()
        bounding_box.add_vertex(Point2D(min_x,min_y))
        bounding_box.add_vertex(Point2D(max_x,min_y))
        bounding_box.add_vertex(Point2D(max_x,max_y))
        bounding_box.add_vertex(Point2D(min_x,max_y))
        return bounding_box
    
    def determine_point_location(self,point):
        """
        Returns 1 if the point lies inside the MultiPolygon, 2 if it is on the 
        boundary and 0 if it is outside the MultiPolygon.
        """
        for polygon in self.polygons:
            result = polygon.determine_point_location(point)
            if result > 0:
                break
        return result
        
        

