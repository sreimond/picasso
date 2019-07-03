# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 13:47:07 2018

@author: sreimond
"""

import numpy as np

class Point2D( object ):
    """
    The `Point2D` class allows to represent a 2D position of terms of a Point 
    object.
    """
    def __init__(self,x,y):
        self.x = float(x)
        self.y = float(y)
        
    def calculate_distance(self,point):
        return np.sqrt((self.x-point.x)**2.0+(self.y-point.y)**2.0)
        
    def is_left_of(self,point1,point2):
        return ((  ((point2.x-point1.x)*(self.y-point1.y)) 
                 - ((point2.y-point1.y)*(self.x-point1.x)) )
                 >= 0 )
                 
                 