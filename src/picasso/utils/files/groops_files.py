# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 09:33:41 2018

@author: sreimond
"""

from ..statistics import least_squares
from ..algebra import matrices
import numpy as np
import os, sys
import inspect
import warnings

def read_groops_normals( file_path ):
    """
    Imports normal equation and right hand side from GROOPS output path and
    returns LeastSquares class instance. 
    """
    # LS class
    LS = least_squares.RegularizedOrdinaryLeastSquares()    
    # full file names    
    f_norm = file_path + '-normals.txt'
    f_rghs = file_path + '-normals.rightHandSide.txt'
    f_info = file_path + '-normals.info.xml'
    f_x = file_path + '-x.txt'
    f_sigmax = file_path + '-sigmax.txt'
    # check for existence
    c1 = not os.path.exists( f_norm )
    c2 = not os.path.exists( f_rghs )
    c3 = not os.path.exists( f_info )
    if c1 or c2 or c3:
        warnings.warn('Groops files not found!')
        return LS
    # read info file
    f = open( f_info, 'r' )
    for line in f:
        if '<lPl>' in line:
            lPl = float(line[line.index("<lPl>") + 5:line.rindex("</lPl>")])
        if '<parameterCount>' in line:
            nx = int(line[line.index("<parameterCount>") + 
                     16:line.rindex("</parameterCount>")])
        if '<observationCount>' in line:
            nl = int(line[line.index("<observationCount>") + 
                     18:line.rindex("</observationCount>")])
    f.close()
    LS.bb = lPl
    LS.row_count = nl
    LS.col_count = nx
    LS.redundancy = nl-nx    
    # read normals file
    N = np.zeros( (nx,nx), dtype=np.float_ )
    f = open( f_norm, 'r' )
    next(f), next(f)
    ii, ij = -1, -1   
    for line in f:
        ii += 1
        ij += 1
        N[ii,ij:nx] = [float(ni) for ni in line.split()]
    f.close()  
    N += np.triu( N, 1).T        
    LS.normal_equation = matrices.MatrixEquation()
    LS.normal_equation.define_coefficient_matrix(matrices.return_matrix_from_numpy_array(N))
#    # read normals file (slow)
#    N_triu = matrices.UpperTriangularMatrix(LS.col_count)
#    f = open( f_norm, 'r' )
#    next(f), next(f)
#    ix = -1
#    for line in f:
#        values = [float(ni) for ni in line.split()]
#        for value in values:
#            ix += 1
#            N_triu.set_element(ix,value)
#    f.close()
#    LS.normal_equation = matrices.MatrixEquation()
#    LS.normal_equation.define_coefficient_matrix(N_triu.return_full_matrix())
    # read right hand side file
    n = np.genfromtxt( f_rghs, unpack=True, skip_header=2 )
    LS.normal_equation.define_right_hand_side(matrices.return_matrix_from_numpy_array(n))
    # x and sigmax
    x, sigmax = np.ones(LS.col_count)*np.nan, np.ones(LS.col_count)*np.nan
    if os.path.exists( f_x ):
        x = np.genfromtxt( f_x, unpack=True, skip_header=2 )
    if os.path.exists( f_sigmax ):
        sigmax = np.genfromtxt( f_sigmax, unpack=True, skip_header=2 )
    LS.normal_equation.x = matrices.return_matrix_from_numpy_array(x)
    LS.x = matrices.return_matrix_from_numpy_array(x)
    LS.sigma_x = matrices.return_matrix_from_numpy_array(sigmax)
    xNx = LS.x.return_transpose().multiply(LS.normal_equation.A).multiply(LS.x).elements
    nx2 = 2.0 * LS.normal_equation.b.return_transpose().multiply(LS.x).elements   
    xRx = 0
    LS.rss = xNx - nx2 + LS.bb 
    LS.xss = xRx
    LS.variance_unit_weight = LS.rss/LS.redundancy
    LS.parameter_covariance_matrix = LS.normal_equation.A.return_inverse()
    LS.parameter_covariance_matrix.elements *= LS.variance_unit_weight
    LS.sigma_x = LS.parameter_covariance_matrix.return_diagonal_vector()
    LS.sigma_x.elements **= 0.5
    LS.bias = None
    return LS



