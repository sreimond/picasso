# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 11:15:36 2018

@author: sreimond
"""

import numpy as np
import warnings
from .vectors import Vector

class Matrix( object ):
    """
    The `Matrix` class returns a Matrix object (all based on numpy).
    """
    def __init__(self,m,n):
        self.m = m
        self.n = n
        self.elements = np.zeros((m,n))
        
    def set_element(self,ixm,ixn,value):
        self.elements[ixm,ixn] = value
        
    def return_transpose(self):
        matrix = np.transpose(self.elements)
        return return_matrix_from_numpy_array(matrix)        
    
    def return_inverse(self):
        try:            
            matrix = np.linalg.inv(self.elements)
        except np.linalg.LinAlgError as err:
            if 'Singular matrix' in str(err):
                matrix = np.ones((self.m,self.n))*np.nan
            else:
                raise
        return return_matrix_from_numpy_array(matrix)
        
    def return_diagonal_vector(self):
        matrix = np.diag(self.elements)
        return return_matrix_from_numpy_array(matrix)
    
    def return_trace(self):
        return np.trace(self.elements)
        
    def return_normal_matrix(self):
        self_t = self.return_transpose()
        return self_t.multiply(self)   
            
    def determine_rank(self):
        return np.linalg.matrix_rank(self.elements)
    
    def is_square(self):
        return (self.m == self.n)

    def add(self,matrix):
        matrix = self.elements + matrix.elements
        return return_matrix_from_numpy_array(matrix)  
        
    def subtract(self,matrix):
        matrix = self.elements - matrix.elements
        return return_matrix_from_numpy_array(matrix)  
        
    def multiply(self,matrix):
        matrix = np.dot(self.elements, matrix.elements)
        return return_matrix_from_numpy_array(matrix)  
        
    def multiply_elementwise(self,matrix):
        matrix = self.elements * matrix.elements
        return return_matrix_from_numpy_array(matrix)  

class IdentityMatrix( Matrix ):
    """
    The `IdentityMatrix` class is a subclass of the Matrix class.
    """
    def __init__(self,m,n=None):
        n = m if n is None else n
        super(IdentityMatrix, self).__init__(m,n)
        self.initialize()
    
    def initialize(self):
        self.elements[:] = np.eye(N=self.m,M=self.n)
#        for ixm in list(range(self.m)):
#            for ixn in list(range(self.n)):
#                if ixm==ixn:
#                    self.set_element(ixm,ixn,1.0)

class RandomMatrix( Matrix ):
    """
    The `RandomMatrix` class is a subclass of the Matrix class.
    """
    def __init__(self,m,n):
        super(RandomMatrix, self).__init__(m,n)
        self.initialize()
    
    def initialize(self):
        for ixm in list(range(self.m)):
            for ixn in list(range(self.n)):
                self.set_element(ixm,ixn,np.random.random())                
        
class UpperTriangularMatrix( object ):
    """
    The `UpperTriangularMatrix` class returns an UpperTriangularMatrix 
    object (all based on numpy).
    """
    def __init__(self,n):
        self.n = n
        self.elements = np.zeros((n,n))
        self.indices = np.ones((n,n)) * np.nan
        self.initialize_indices()
    
    def initialize_indices(self):
        lin_ix = -1
        row_ix = -1        
        for ixm in list(range(self.n)):
            row_ix += 1
            for ixn in list(range(self.n)):
                if ixn < row_ix: 
                    continue
                lin_ix += 1
                self.indices[ixm,ixn] = lin_ix                
    
    def set_element(self,ix,value):
        indices = np.where(self.indices==ix)
        ixm = indices[0][0]
        ixn = indices[1][0]
        self.elements[ixm,ixn] = value
    
    def return_full_matrix(self):
        matrix = self.elements
        matrix += np.transpose(matrix) - np.eye(self.n) * np.diag(matrix)
        return return_matrix_from_numpy_array(matrix)      
    
class LowerTriangularMatrix( object ):
    """
    The `LowerTriangularMatrix` class returns an LowerTriangularMatrix 
    object (all based on numpy).
    """
    def __init__(self,n):
        self.n = n
        self.elements = np.zeros((n,n))
        self.indices = np.ones((n,n)) * np.nan
        self.initialize_indices()
    
    def initialize_indices(self):
        lin_ix = -1
        row_ix = -1        
        for ixm in list(range(self.n)):
            row_ix += 1
            for ixn in list(range(row_ix+1)):
                lin_ix += 1
                self.indices[ixm,ixn] = lin_ix                
    
    def set_element(self,ix,value):
        indices = np.where(self.indices==ix)
        ixm = indices[0][0]
        ixn = indices[1][0]
        self.elements[ixm,ixn] = value
    
    def return_full_matrix(self):
        matrix = self.elements
        matrix += np.transpose(matrix) - np.eye(self.n) * np.diag(matrix)
        return return_matrix_from_numpy_array(matrix)  
    
class MatrixEquation(object):
    """
    The `MatrixEquation` class returns a MatrixEquation object (all based on 
    numpy, scipy).Ax=b.
    """
    def __init__(self,A=None,b=None,x=None):
        self.A = A
        self.b = b
        self.x = x
    
    def define_coefficient_matrix(self,matrix):
        if not isinstance(matrix,Matrix):
            warnings.warn("Please use the Matrix class.")
            return
        self.A = matrix
        
    def define_right_hand_side(self,matrix):
        if not (isinstance(matrix,Matrix) or isinstance(matrix,Vector)):
            warnings.warn("Please use the Matrix or Vector class.")
            return
        self.b = matrix
    
    def solve(self):
        if not self.A.is_square():
            warnings.warn("Coefficient matrix must be square! Consider least squares instead.")
            return
        A = self.A.elements
        b = self.b.elements
        try:            
            x = np.linalg.solve(A,b)
        except np.linalg.LinAlgError as err:
            if 'Singular matrix' in str(err):
                x = np.ones(b.shape)*np.nan
            else:
                raise
        self.x = return_matrix_from_numpy_array(x)    
   
def return_matrix_from_numpy_array(numpy_array):
    if numpy_array.ndim<=1:
        numpy_array = np.transpose(np.array(numpy_array,ndmin=2))
    m = numpy_array.shape[0]
    n = numpy_array.shape[1]        
    matrix = Matrix(m,n)
    matrix.elements[:] = numpy_array[:]
#    for ixm in list(range(m)):
#        for ixn in list(range(n)):
#            matrix.set_element(ixm,ixn,numpy_array[ixm,ixn])
    return matrix   







        