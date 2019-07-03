# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 10:09:07 2018

@author: sreimond
"""

import numpy as np

class Vector( object ):
    """
    The `Vector` class allows to represent two 2D position of terms of a Vector 
    object.
    """
    def __init__(self):
        self.elements = []
    
    def add_element(self,element):
        self.elements.append(float(element))
    
    def count_elements(self):
        return len(self.elements)
    
    def cross(self,vector):
        return np.cross(self.elements,vector.elements)
        
    def dot(self,vector):
        return np.dot(self.elements,vector.elements)

    def norm(self):
        return np.sqrt(self.dot(self))

    def angle_between(self,vector):
        return np.arccos(self.dot(vector)/(self.norm()*vector.norm()))
        
        
      