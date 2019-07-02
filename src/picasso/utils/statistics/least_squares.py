# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 14:04:05 2018

@author: sreimond
"""

import numpy as np
from copy import deepcopy
import itertools
from ..algebra.matrices import Matrix
from ..algebra.matrices import IdentityMatrix
from ..algebra.matrices import MatrixEquation
from ..algebra.matrices import return_matrix_from_numpy_array 
from ..algebra.vectors import Vector

class OrdinaryLeastSquares( object ):
    """ 
    The `OrdinaryLeastSquares` class.
    """    
    def __init__(self):
        self.A = None
        self.b = None
        self.x = None
        
    def define_design_matrix(self,matrix):
        if not isinstance(matrix,Matrix):
            warnings.warn("Please use the Matrix class.")
            return
        self.A = matrix
        self.update_dimensions()
        
    def define_right_hand_side(self,matrix):
        if not (isinstance(matrix,Matrix)):
            warnings.warn("Please use the Matrix class.")
            return
        self.b = matrix
        self.bb = self.b.return_transpose().multiply(self.b).elements

    def update_dimensions(self):
        self.row_count = self.A.m
        self.col_count = self.A.n
        self.redundancy = self.row_count - self.col_count
        
    def determine_normal_equation(self):
        self.normal_equation = MatrixEquation()
        N = self.A.return_normal_matrix()
        self.normal_equation.define_coefficient_matrix(N)
        n = self.A.return_transpose().multiply(self.b)
        self.normal_equation.define_right_hand_side(n)        
    
    def solve(self):
        self.normal_equation.solve()
        self.x = self.normal_equation.x
        xNx = self.x.return_transpose().multiply(self.normal_equation.A).multiply(self.x).elements
        nx2 = 2.0 * self.normal_equation.b.return_transpose().multiply(self.x).elements        
        self.rss = xNx - nx2 + self.bb   # residual sum of squares
        self.variance_unit_weight = self.rss/self.redundancy
        self.parameter_covariance_matrix = self.normal_equation.A.return_inverse()
        self.parameter_covariance_matrix.elements *= self.variance_unit_weight
        self.sigma_x = self.parameter_covariance_matrix.return_diagonal_vector()
        self.sigma_x.elements **= 0.5
    
    def synthesis(self):
        ANA = self.A.multiply(self.A.return_normal_matrix().return_inverse()).multiply(self.A.return_transpose())
        self.adjusted_right_hand_side = self.A.multiply(self.x)
        self.adjusted_right_hand_side_covariance_matrix = ANA
        self.adjusted_right_hand_side_covariance_matrix.elements *= self.variance_unit_weight
        self.residuals = IdentityMatrix(self.row_count).subtract(ANA).multiply(self.b)
        self.residuals.elements *= -1.0
        self.residuals_covariance_matrix = IdentityMatrix(self.row_count).subtract(ANA)
        self.residuals_covariance_matrix.elements *= self.variance_unit_weight        

class GeneralizedLeastSquares( object ):
    """ 
    The `GeneralizedLeastSquares` class. Supporting covariance matrix of b.
    """    
    def __init__(self):
        self.A = None
        self.right_hand_side_covariance_matrix = None
        self.b = None
        self.x = None
        
    def define_design_matrix(self,matrix):
        if not isinstance(matrix,Matrix):
            warnings.warn("Please use the Matrix class.")
            return
        self.A = matrix
        self.update_dimensions()
        
    def define_right_hand_side_covariance_matrix(self,matrix):
        if not (isinstance(matrix,Matrix)):
            warnings.warn("Please use the Matrix class.")
            return
        self.right_hand_side_covariance_matrix = matrix
        self.P = matrix.return_inverse()
        
    def define_right_hand_side(self,matrix):
        if not (isinstance(matrix,Matrix)):
            warnings.warn("Please use the Matrix class.")
            return
        self.b = matrix
        self.bPb = self.b.return_transpose().multiply(self.P).multiply(self.b).elements

    def update_dimensions(self):
        self.row_count = self.A.m
        self.col_count = self.A.n
        self.redundancy = self.row_count - self.col_count
        
    def determine_normal_equation(self):
        self.normal_equation = MatrixEquation()
        N = self.A.return_transpose().multiply(self.P).multiply(self.A)
        self.normal_equation.define_coefficient_matrix(N)
        n = self.A.return_transpose().multiply(self.P).multiply(self.b)
        self.normal_equation.define_right_hand_side(n)        
    
    def solve(self):
        self.normal_equation.solve()
        self.x = self.normal_equation.x
        xNx = self.x.return_transpose().multiply(self.normal_equation.A).multiply(self.x).elements
        nx2 = 2.0 * self.normal_equation.b.return_transpose().multiply(self.x).elements        
        self.rss = xNx - nx2 + self.bPb   # residual sum of squares
        self.variance_unit_weight = self.rss/self.redundancy
        self.parameter_covariance_matrix = self.normal_equation.A.return_inverse()
        self.parameter_covariance_matrix.elements *= self.variance_unit_weight
        self.sigma_x = self.parameter_covariance_matrix.return_diagonal_vector()
        self.sigma_x.elements **= 0.5

    def synthesis(self):
        ANA = self.A.multiply(self.normal_equation.A.return_inverse()).multiply(self.A.return_transpose())
        self.adjusted_right_hand_side = ANA.multiply(self.P).multiply(self.b)
        self.adjusted_right_hand_side_covariance_matrix = ANA
        self.adjusted_right_hand_side_covariance_matrix.elements *= self.variance_unit_weight
        self.residuals = IdentityMatrix(self.row_count).subtract(ANA.multiply(self.P)).multiply(self.b)
        self.residuals.elements *= -1.0
        self.residuals_covariance_matrix = self.right_hand_side_covariance_matrix.subtract(ANA)
        self.residuals_covariance_matrix.elements *= self.variance_unit_weight         
        
class RegularizedOrdinaryLeastSquares( object ):
    """ 
    The `RegularizedLeastSquares` class. Tikhonov.
    """    
    def __init__(self):
        self.A = None
        self.b = None
        self.x = None
        self.R = None
        
    def define_design_matrix(self,matrix):
        if not isinstance(matrix,Matrix):
            warnings.warn("Please use the Matrix class.")
            return
        self.A = matrix
        self.update_dimensions()
        
    def define_right_hand_side(self,matrix):
        if not (isinstance(matrix,Matrix)):
            warnings.warn("Please use the Matrix class.")
            return
        self.b = matrix
        self.bb = self.b.return_transpose().multiply(self.b).elements

    def define_regularization_matrix(self,matrix):
        if not isinstance(matrix,Matrix):
            warnings.warn("Please use the Matrix class.")
            return
        self.R = matrix
        
    def update_dimensions(self):
        self.row_count = self.A.m
        self.col_count = self.A.n
        self.redundancy = self.row_count - self.col_count
        
    def determine_normal_equation(self):
        self.normal_equation = MatrixEquation()
        N = self.A.return_normal_matrix()
        self.normal_equation.define_coefficient_matrix(N)
        n = self.A.return_transpose().multiply(self.b)
        self.normal_equation.define_right_hand_side(n)        
    
    def determine_regularized_normal_equation(self):
        self.regularized_normal_equation = MatrixEquation()
        N = self.A.return_normal_matrix()
        N_regularized = N.add(self.R)
        self.regularized_normal_equation.define_coefficient_matrix(N_regularized)
        n = self.A.return_transpose().multiply(self.b)
        self.regularized_normal_equation.define_right_hand_side(n)        

    def determine_vce_regularization(self,sig1=1.0,sigu=1e2):
        iteration_count = 20
        for ix in list(range(iteration_count)):
            if ix==0:
                r1 = self.row_count - self.normal_equation.A.multiply(self.normal_equation.A.return_inverse()).return_trace()  
                ru = self.col_count - self.normal_equation.A.return_inverse().return_trace()  
            else:
                r1 = self.row_count - self.regularized_normal_equation.A.multiply(self.regularized_normal_equation.A.return_inverse()).return_trace()/sig1
                ru = self.col_count - self.regularized_normal_equation.A.return_inverse().return_trace()/sigu
                sig1 = self.rss/r1
                sigu = self.xss/ru
            if sigu<np.finfo(float).tiny:
                alpha = 0
            else:
                alpha = sig1/sigu
            R = alpha * np.eye(self.col_count)
            self.define_regularization_matrix( return_matrix_from_numpy_array(R) )
            self.regularized_normal_equation = MatrixEquation()
            N_regularized = self.normal_equation.A.add(self.R)
            self.regularized_normal_equation.define_coefficient_matrix(N_regularized)
            self.regularized_normal_equation.define_right_hand_side(self.normal_equation.b)            
            self.solve()            
            
    def solve(self):
        self.regularized_normal_equation.solve()
        self.x = self.regularized_normal_equation.x
        xNx = self.x.return_transpose().multiply(self.normal_equation.A).multiply(self.x).elements
        nx2 = 2.0 * self.normal_equation.b.return_transpose().multiply(self.x).elements   
#        xRx = self.x.return_transpose().multiply(self.R).multiply(self.x).elements
        xRx = self.x.return_transpose().multiply(self.x).elements
        self.rss = xNx - nx2 + self.bb   # residual sum of squares
        self.xss = xRx                   # regularization term: parameter sum of squares
        self.variance_unit_weight = self.rss/self.redundancy
        self.parameter_covariance_matrix = self.regularized_normal_equation.A.return_inverse().multiply(self.normal_equation.A).multiply(self.regularized_normal_equation.A.return_inverse()) # from quadrature error propagation (E.g.: https://de.scribd.com/document/288575917/Regularization)
        self.parameter_covariance_matrix.elements *= self.variance_unit_weight
        self.sigma_x = self.parameter_covariance_matrix.return_diagonal_vector()
        self.sigma_x.elements **= 0.5
        self.bias = self.regularized_normal_equation.A.return_inverse().multiply(self.R).multiply(self.x)   # Sciences of Geodesy - II: Innovations and Future Developments

#    def synthesis(self):
#        ANA = self.A.multiply(self..A.return_inverse()).multiply(self.A.return_transpose())
#        self.adjusted_right_hand_side = self.A.multiply(self.x)
#        self.adjusted_right_hand_side_covariance_matrix = ANA
#        self.adjusted_right_hand_side_covariance_matrix.elements *= self.variance_unit_weight
##        self.regularized_residuals = ANA.subtract(self.A).multiply(self.b)
#        self.residuals = IdentityMatrix(self.row_count).subtract(ANA).multiply(self.b)
#        self.residuals.elements *= -1.0
#        self.residuals_covariance_matrix = IdentityMatrix(self.row_count).subtract(ANA)
#        self.residuals_covariance_matrix.elements *= self.variance_unit_weight        

class RegularizedGeneralizedLeastSquares( object ):
    """ 
    The `RegularizedGeneralizedLeastSquares` class. Tikhonov.
    """    
    def __init__(self):
        self.A = None
        self.b = None
        self.x = None
        self.R = None
        self.right_hand_side_covariance_matrix = None
        
    def define_design_matrix(self,matrix):
        if not isinstance(matrix,Matrix):
            warnings.warn("Please use the Matrix class.")
            return
        self.A = matrix
        self.update_dimensions()

    def define_right_hand_side_covariance_matrix(self,matrix):
        if not (isinstance(matrix,Matrix)):
            warnings.warn("Please use the Matrix class.")
            return
        self.right_hand_side_covariance_matrix = matrix
        self.P = matrix.return_inverse()        
        
    def define_right_hand_side(self,matrix):
        if not (isinstance(matrix,Matrix)):
            warnings.warn("Please use the Matrix class.")
            return
        self.b = matrix
        self.bPb = self.b.return_transpose().multiply(self.P).multiply(self.b).elements
        
    def define_regularization_matrix(self,matrix):
        if not isinstance(matrix,Matrix):
            warnings.warn("Please use the Matrix class.")
            return
        self.R = matrix
        
    def update_dimensions(self):
        self.row_count = self.A.m
        self.col_count = self.A.n
        self.redundancy = self.row_count - self.col_count
        
    def determine_normal_equation(self):
        self.normal_equation = MatrixEquation()
        N = self.A.return_transpose().multiply(self.P).multiply(self.A)
        self.normal_equation.define_coefficient_matrix(N)
        n = self.A.return_transpose().multiply(self.P).multiply(self.b)
        self.normal_equation.define_right_hand_side(n)        
    
    def determine_regularized_normal_equation(self):
        self.regularized_normal_equation = MatrixEquation()
        N = self.A.return_normal_matrix()
        N_regularized = N.add(self.R)
        self.regularized_normal_equation.define_coefficient_matrix(N_regularized)
        n = self.A.return_transpose().multiply(self.b)
        self.regularized_normal_equation.define_right_hand_side(n)        
    
    def determine_vce_regularization(self,sig1=1.0,sigu=1e2):
        iteration_count = 100
        for ix in list(range(iteration_count)):
            if ix==0:
                r1 = self.row_count - self.normal_equation.A.multiply(self.normal_equation.A.return_inverse()).return_trace()  
                ru = self.col_count - self.normal_equation.A.return_inverse().return_trace()  
            else:
                r1 = self.row_count - self.regularized_normal_equation.A.multiply(self.regularized_normal_equation.A.return_inverse()).return_trace()/sig1
                ru = self.col_count - self.regularized_normal_equation.A.return_inverse().return_trace()/sigu
                sig1 = self.rss/r1
                sigu = self.xss/ru
            if sigu<np.finfo(float).tiny:
                alpha = 0
            else:
                alpha = sig1/sigu
            R = alpha * np.eye(self.col_count)
            self.define_regularization_matrix( return_matrix_from_numpy_array(R) )
            self.regularized_normal_equation = MatrixEquation()
            N_regularized = self.normal_equation.A.add(self.R)
            self.regularized_normal_equation.define_coefficient_matrix(N_regularized)
            self.regularized_normal_equation.define_right_hand_side(self.normal_equation.b)            
            self.solve()        
            
    def solve(self):
        self.regularized_normal_equation.solve()
        self.x = self.regularized_normal_equation.x
        xNx = self.x.return_transpose().multiply(self.normal_equation.A).multiply(self.x).elements
        nx2 = 2.0 * self.normal_equation.b.return_transpose().multiply(self.x).elements
#        xRx = self.x.return_transpose().multiply(self.R).multiply(self.x).elements
        xRx = self.x.return_transpose().multiply(self.x).elements
        self.rss = xNx - nx2 + self.bPb   # residual sum of squares
        self.xss = xRx              # regularization term: parameter sum of squares
        self.variance_unit_weight = self.rss/self.redundancy
        self.parameter_covariance_matrix = self.regularized_normal_equation.A.return_inverse().multiply(self.normal_equation.A).multiply(self.regularized_normal_equation.A.return_inverse())
        self.parameter_covariance_matrix.elements *= self.variance_unit_weight
        self.sigma_x = self.parameter_covariance_matrix.return_diagonal_vector()
        self.sigma_x.elements **= 0.5
        self.bias = self.regularized_normal_equation.A.return_inverse().multiply(self.R).multiply(self.x)

#    def synthesis(self):
#        ANA = self.A.multiply(self.normal_equation.A.return_inverse()).multiply(self.A.return_transpose())
#        self.adjusted_right_hand_side = ANA.multiply(self.P).multiply(self.b)
#        self.adjusted_right_hand_side_covariance_matrix = ANA
#        self.adjusted_right_hand_side_covariance_matrix.elements *= self.variance_unit_weight
##        self.regularized_residuals = ANA.subtract(self.A).multiply(self.b)
#        self.residuals = IdentityMatrix(self.row_count).subtract(ANA).multiply(self.b)
#        self.residuals.elements *= -1.0
#        self.residuals_covariance_matrix = IdentityMatrix(self.row_count).subtract(ANA)
#        self.residuals_covariance_matrix.elements *= self.variance_unit_weight    

