# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 08:59:18 2018

@author: sreimond
"""

from .points import Point2D
from ..algebra.vectors import Vector
import numpy as np
import warnings

class Line2D( object ):
    """
    The `Line2D` class allows to represent two or more Point2D objects in terms
    of a Line (line segment) object. 
    """
    def __init__(self,point1=[],point2=[]):
        self.endpoints = []
        if isinstance(point1,tuple):
            point = Point2D(*point1)
            self.add_endpoint(point)
        elif isinstance(point1,Point2D):
            self.add_endpoint(point1)
        if isinstance(point2,tuple):
            point = Point2D(*point2)
            self.add_endpoint(point)
        elif isinstance(point2,Point2D):
            self.add_endpoint(point2)
    
    def add_endpoint(self,point):
        """
        Add a Point2D object to the list of endpoint (max. 2).
        """     
        if len(self.endpoints)==2:
            warnings.warn('Tried to add a third point to a Line2D object. '
                          'Ignored. Create a new Line2D object instead.')
            return
        if len(self.endpoints)==0:
            self.endpoints.append((point))
            return
        if self.endpoints[-1].calculate_distance(point) > 0:
            self.endpoints.append((point))   
        else:
            warnings.warn('The second endpoint must differ from the first!')    
            
    def as_vector(self):
        """ Returns the difference vector B-A as a Vector object. """
        if len(self.endpoints)<2:
            warnings.warn('Only one endpoint defined!')
            return 
        dx = self.endpoints[1].x-self.endpoints[0].x
        dy = self.endpoints[1].y-self.endpoints[0].y
        v = Vector()
        v.add_element(dx)
        v.add_element(dy)
        return v
        
    def calculate_length(self):
        return np.sqrt(  (self.endpoints[0].x-self.endpoints[1].x)**2.0
                       + (self.endpoints[0].y-self.endpoints[1].y)**2.0 )
                       
    def calculate_epsilon(self):
#        return 1e-12
        return 1e-5 * self.calculate_length()

    def calculate_point_distance(self,point):
        """
        Computes the perpendicular distance from a point to the line.
        """
        x1 = self.endpoints[0].x
        y1 = self.endpoints[0].y
        x2 = self.endpoints[1].x
        y2 = self.endpoints[1].y
        x0 = point.x
        y0 = point.y
        return np.fabs((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1)/self.calculate_length()

    def determine_point_location(self,point):
        """
        Returns 1 if the Point2D lies on the line and is inbetween the two 
        endpoints, 2 if the point is on the line but not between the two 
        endpoints and 0 if the point is not on the line at all.
        """
        if (point.calculate_distance(self.endpoints[0]) == 0 or
            point.calculate_distance(self.endpoints[1]) == 0):
                return 1
        dAB = self.as_vector()
        dX = Line2D(self.endpoints[0],point).as_vector()
        cp = dAB.cross(dX)
        is_on_line = np.fabs(cp) < self.calculate_epsilon()
        is_between_x = (((self.endpoints[0].x <= point.x) and
                         (point.x <= self.endpoints[1].x))
                        or
                        ((self.endpoints[1].x <= point.x) and
                        (point.x <= self.endpoints[0].x)))
        is_almost_x = ((abs(self.endpoints[0].x-point.x) < self.calculate_epsilon()) 
                       or
                       (abs(self.endpoints[1].x-point.x)< self.calculate_epsilon()))      
        is_between_y = (((self.endpoints[0].y <= point.y) and
                         (point.y <= self.endpoints[1].y))
                        or
                        ((self.endpoints[1].y <= point.y) and
                        (point.y <= self.endpoints[0].y)))
        is_almost_y = ((abs(self.endpoints[0].y-point.y) < self.calculate_epsilon()) 
                       or
                       (abs(self.endpoints[1].y-point.y)< self.calculate_epsilon()))
        if is_on_line and (is_between_x or is_almost_x) and (is_between_y or is_almost_y):
            return 1
        if is_on_line:
            return 2
        else:
            return 0
    
    def determine_line_parallel(self,line):
        """
        Checks if the two lines are parallel. Returns True if they are and #
        False they are not parallel.
        """        
        dAB = self.as_vector()
        dCD = line.as_vector()
        alpha = np.fabs(dAB.angle_between(dCD))
        if (alpha < self.calculate_epsilon() or 
            np.fabs(alpha-np.pi) < self.calculate_epsilon()):
            return True
        else:            
            return False
    
    def calculate_line_intersection_point(self,line):
        """
        Returns a Point2D object if the two line intersect. If the lines are
        parallel, this object is empty and a warning is produced.        
        """
        if self.determine_line_parallel(line):
            warnings.warn('Lines are parellel! No intersection point.')
            return Point2D(np.nan,np.nan)
        x1 = self.endpoints[0].x
        y1 = self.endpoints[0].y
        x2 = self.endpoints[1].x
        y2 = self.endpoints[1].y
        x3 = line.endpoints[0].x
        y3 = line.endpoints[0].y
        x4 = line.endpoints[1].x
        y4 = line.endpoints[1].y
        xi = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
        yi = ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
        return Point2D(xi,yi)
    
    def determine_line_intersection(self,line):
        """
        Returns 1 if the segments of the Line2D and the line intersect, 2 if 
        the two lines intersect anywhere, 0 if they are parallel.
        """
        if self.determine_line_parallel(line):
            return 0            
        point = self.calculate_line_intersection_point(line)
        on_self = self.determine_point_location(point) == 1
        on_line = line.determine_point_location(point) == 1
        if on_self and on_line:
            return 1
        else:
            return 2
        
        
        
    
            
    






