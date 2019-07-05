# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:35:48 2018

@author: sreimond
"""

from picasso.utils.algebra import matrices
from picasso.utils.files import groops_files
import numpy as np
import pkg_resources
import os

def regularization(normals_path,vce=True):
    LS = groops_files.read_groops_normals(normals_path)
    if vce:
        LS = _vce(LS)
    else:
        LS = _no_regularization(LS)
    return LS

def regularization_matlab(normals_path,method='lcurve'):
    matlab_path = pkg_resources.resource_filename('picasso.data', 'MATLAB')
    LS = groops_files.read_groops_normals(normals_path)
    LS = _no_regularization(LS)
    sys_str = "matlab -nodisplay -r \""
    sys_str += "addpath('%s'); " % matlab_path
    sys_str += "addpath('%s/regularization_tools'); " % matlab_path
    sys_str += "normals2regtools('%s','%s'); " % (normals_path,method)
    sys_str += "exit;\""
    os.system(sys_str)
    x = np.genfromtxt(normals_path+'-x_regtools.txt')
    sigma_x = np.genfromtxt(normals_path+'-sigmax_regtools.txt')
    LS.x = matrices.return_matrix_from_numpy_array( np.array(x) )    
    LS.sigma_x = matrices.return_matrix_from_numpy_array( np.array(sigma_x) ) 
    return LS

def _vce(LS):
    LS.determine_vce_regularization(sig1=1.0,sigu=1e2)
    return LS

def _no_regularization(LS):
    R = matrices.IdentityMatrix(LS.col_count)
    R.elements *= 0
    LS.define_regularization_matrix(R)
    LS.regularized_normal_equation = matrices.MatrixEquation()
    N_regularized = LS.normal_equation.A.add(LS.R)
    LS.regularized_normal_equation.define_coefficient_matrix(N_regularized)
    LS.regularized_normal_equation.define_right_hand_side(LS.normal_equation.b)
    LS.solve()
    return LS
