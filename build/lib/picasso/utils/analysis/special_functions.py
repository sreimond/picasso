# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 15:27:59 2018

@author: sreimond
"""

import numpy as np

def legendre1( max_degree, x ):
    """
    The `legendre1` function returns the fully normalized associated Legendre 
    functions of the first kind up to degree max_degree for arguments x.
    
    Requires the functions: `_ravel_ix`, `_unravel_ix`.
    
    Reference: Fukushima 2012 (DOI: 10.1007/s00190-012-0561-8)
    """
    # initialize array
    x_count = np.array(x).size
    coeff_count = int((max_degree+1) * (max_degree+2) * 0.5)
    P = np.zeros((x_count,coeff_count),dtype=np.float_)
    # auxiliary indices 
    n_array = np.array( range(0,max_degree+1) )  
    l_array = np.array( range(0,coeff_count) )  
    nl,ml   = _unravel_ix( l_array )
    # auxiliary variables
    s = np.sin(np.array(x))
    t = np.cos(np.array(x))  
    with np.errstate(divide='ignore', invalid='ignore'):
        n  = np.array( nl, dtype=np.float_ )
        m  = np.array( ml, dtype=np.float_ )
        an = (2.0 * n + 1.0) * (2.0 * n - 1.0)
        ad = (n + m) * (n - m)
        bn = (2.0 * n + 1.0) * (n + m - 1.0) * (n - m - 1.0)
        bd = (2.0 * n - 3.0) * (n + m) * (n - m)
        a  = np.sqrt(an/ad)
        b  = np.sqrt(bn/bd)
    # seed values
    P[:,0] = 1.0
    P[:,1] = np.sqrt(3.0) * t
    P[:,2] = np.sqrt(3.0) * s  
    # sectorials
    imm = _ravel_ix(n_array,n_array)  
    i1  = imm[2:]
    i2  = imm[1:-1]  
    m   = n_array[2:]
    d   = ((2.0*m+1.0)/(2.0*m)) ** 0.5  
    P[:,i1] = np.outer( s, d )
    for i in range(i2.size):
        P[:,i1[i]] *= P[:,i2[i]]
    # semi-sectorials
    i1 = _ravel_ix(n_array[1:-1]+1,n_array[1:-1])
    P[:,i1] = np.outer( t, a[i1] ) * P[:,i2]
    # all others
    i1      = np.array([l for l in l_array if nl[l] >= (ml[l] + 2)], dtype=np.int_)
    n1, m1  = _unravel_ix( i1 )
    i2      = _ravel_ix( n1-1, m1 )
    i3      = _ravel_ix( n1-2, m1 )  
    for i in range(i1.size):
        P[:,i1[i]] = t * a[i1[i]] * P[:,i2[i]] - P[:,i3[i]] * b[i1[i]]
    return P
  
def _ravel_ix( n, m ):
    return np.array(0.5 * n * (n+1) + m,dtype=np.int_)

def _unravel_ix( l ):
    n = np.array( -0.5 + (0.25 + 2.0*l)**0.5, dtype=np.int_)
    m = np.array( l - n * 0.5 * (n+1), dtype=np.int_)
    return n, m