# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
from picasso.utils.algebra import matrices
import numpy as np

@pytest.fixture
def matrix1():
    m = matrices.Matrix(3,3)
    m.set_element(0,0,0)
    m.set_element(1,0,1)
    m.set_element(2,0,-3)
    m.set_element(0,1,-3)
    m.set_element(1,1,-4)
    m.set_element(2,1,4)
    m.set_element(0,2,-2)
    m.set_element(1,2,-2)
    m.set_element(2,2,1)
    return m
	
@pytest.fixture
def matrix2():
    m = matrices.Matrix(3,1)
    m.set_element(0,0,4)
    m.set_element(1,0,7)
    m.set_element(2,0,-3)
    return m
    

def test_matrix(matrix1):
    """
    Test general functionality of the Matrix class.
    """
    assert matrix1.elements[1,1] == -4
    assert matrix1.elements[0,2] == -2

def test_transpose(matrix1):
    """
    Test transpose functionality of the Matrix class.
    """
    mt = matrix1.return_transpose()
    assert mt.elements[1,1] == -4
    assert mt.elements[0,2] == -3
    
def test_inverse(matrix1):
    """
    Test inverse functionality of the Matrix class.
    """
    mt = matrix1.return_inverse()
    assert mt.elements[1,1] == pytest.approx(-6,1e-15)
    assert mt.elements[0,2] == pytest.approx(-2,1e-15)
	
def test_diagonal(matrix1):
    """
    Test diagonal functionality of the Matrix class.
    """
    mt = matrix1.return_diagonal_vector()
    assert mt.elements[0,0] == 0
    assert mt.elements[1,0] == -4
    assert mt.elements[2,0] == 1
    
def test_trace(matrix1):
    """
    Test trace functionality of the Matrix class.
    """
    mt = matrix1.return_trace()
    assert mt == -3

def test_normal(matrix1):
    """
    Test normal matrix functionality of the Matrix class.
    """
    mt = matrix1.return_normal_matrix()
    assert mt.elements[1,1] == 41
    assert mt.elements[0,2] == -5
	
def test_rank(matrix1):
    """
    Test rank functionality of the Matrix class.
    """
    mt = matrix1.determine_rank()
    assert mt == 3
	
def test_is_square(matrix1):
    """
    Test is_square functionality of the Matrix class.
    """
    mt = matrix1.is_square()
    assert mt
	
def test_add(matrix1):
    """
    Test add matrix functionality of the Matrix class.
    """
    mt = matrix1.add(matrix1)
    assert mt.elements[1,1] == -8
    assert mt.elements[0,2] == -4
	
def test_subtract(matrix1):
    """
    Test subtract matrix functionality of the Matrix class.
    """
    mt = matrix1.subtract(matrix1)
    assert mt.elements[1,1] == 0
    assert mt.elements[0,2] == 0
	
def test_multiply(matrix1):
    """
    Test multiply matrix functionality of the Matrix class.
    """
    mt = matrix1.multiply(matrix1)
    assert mt.elements[1,1] == 5
    assert mt.elements[0,2] == 4
	
def test_multiply_elementwise(matrix1):
    """
    Test multiply_elementwise matrix functionality of the Matrix class.
    """
    mt = matrix1.multiply_elementwise(matrix1)
    assert mt.elements[1,1] == 16
    assert mt.elements[0,2] == 4
	
	
def test_identity(matrix1):
    """
    Test identity matrix functionality of the Matrix class.
    """
    mt = matrices.IdentityMatrix(4)
    assert mt.elements[1,1] == 1
    assert mt.return_trace() == 4
	
def test_random(matrix1):
    """
    Test random matrix functionality of the Matrix class.
    """
    mt = matrices.RandomMatrix(4,4)
    assert mt.m == 4
    assert mt.elements[1,1] != 0
	
def test_upper_triangular(matrix1):
    """
    Test upper triangular matrix functionality of the Matrix class.
    """
    mt = matrices.UpperTriangularMatrix(4)
    assert mt.n == 4
	
def test_lower_triangular(matrix1):
    """
    Test lower triangular matrix functionality of the Matrix class.
    """
    mt = matrices.LowerTriangularMatrix(4)
    assert mt.n == 4

def test_matrix_equation(matrix1,matrix2):
    """
    Test matrix equation functionality of the Matrix class.
    """
    mt = matrices.MatrixEquation(A=matrix1,b=matrix2)
    mt.solve()
    assert mt.x.elements[0,0] == pytest.approx(-13,1e-15)
    assert mt.x.elements[1,0] == pytest.approx(-16,1e-15)
    assert mt.x.elements[2,0] == pytest.approx(22,1e-15)

		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
