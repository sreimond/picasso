# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 11:32:25 2018

@author: sreimond
"""

import numpy as np
from .objective_functions import ObjectiveFunction, ObjectiveFunctionParameter
from .objective_functions import MultiObjectiveFunction

def six_hump_camel_back():
    f = ObjectiveFunction( _six_hump_camel_back )
    f.global_minimum = -1.031628453
    xst = lambda x:(x>=-3 and x<=3)
    yst = lambda y:(y>=-2 and y<=2)
    xm = lambda x: 3.0*np.sin(x*2.0*np.pi)
    ym = lambda y: 2.0*np.sin(y*2.0*np.pi)
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p1.true_value = 0.089842
    p1.bits = 4
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    p2.true_value = -0.712656
    p2.bits = 4
    f.add_parameter(p1)
    f.add_parameter(p2)
    return f

def _six_hump_camel_back( parameters ):
    x = [p.value for p in parameters]
    return (4.0 - 2.1*x[0]**2.0 + (x[0]**4.0)/3.0 ) * x[0]**2.0 + x[0]*x[1] + (-4.0+4.0*x[1]**2.0)*x[1]**2.0    

def constrained_rosenbrock1():
    f = ObjectiveFunction( _constrained_rosenbrock1 )
    f.global_minimum = 0
    xst = lambda x:(x>=-1.5 and x<=1.5)
    yst = lambda y:(y>=-0.5 and y<=2.5)
    xm = lambda x: (1.5-(-1.5))*0.5*np.sin(x*2.0*np.pi)+(1.5+(-1.5))*0.5
    ym = lambda y: (2.5-(-0.5))*0.5*np.sin(y*2.0*np.pi)+(2.5+(-0.5))*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p1.true_value = 1.0
    p1.bits = 4
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    p2.true_value = 1.0
    p2.bits = 4
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_subject(_subject_constrained_rosenbrock1a)
    f.add_subject(_subject_constrained_rosenbrock1b)
    return f
    
def _constrained_rosenbrock1( parameters ):
    x = [p.value for p in parameters]    
    return (1.0-x[0])**2.0 + 100.0*(x[1]-x[0]**2.0)**2.0

def _subject_constrained_rosenbrock1a( parameters ):
    x = [p.value for p in parameters]
    return ((x[0]-1.0)**3.0 - x[1] + 1) < 0
    
def _subject_constrained_rosenbrock1b( parameters ):
    x = [p.value for p in parameters]
    return (x[0]+x[1]-2.0) < 0
    
def constrained_rosenbrock2():
    f = ObjectiveFunction( _constrained_rosenbrock2 )
    f.global_minimum = 0
    xst = lambda x:(x>=-1.5 and x<=1.5)
    yst = lambda y:(y>=-1.5 and y<=1.5)
    xm = lambda x: (1.5-(-1.5))*0.5*np.sin(x*2.0*np.pi)+(1.5+(-1.5))*0.5
    ym = lambda y: (1.5-(-1.5))*0.5*np.sin(y*2.0*np.pi)+(1.5+(-1.5))*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p1.true_value = 1.0
    p1.bits = 4
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    p2.true_value = 1.0
    p2.bits = 4
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_subject(_subject_constrained_rosenbrock2a)
    return f
    
def _constrained_rosenbrock2( parameters ):
    x = [p.value for p in parameters]    
    return (1.0-x[0])**2.0 + 100.0*(x[1]-x[0]**2.0)**2.0

def _subject_constrained_rosenbrock2a( parameters ):
    x = [p.value for p in parameters]
    return (x[0]**2.0 + x[1]**2.0)<2 

def townsend_modified():
    f = ObjectiveFunction( _townsend_modified )
    f.global_minimum = -2.0239884
    xst = lambda x:(x>=-2.25 and x<=2.5)
    yst = lambda y:(y>=-2.5 and y<=1.75)
    xm = lambda x: (2.5-(-2.25))*0.5*np.sin(x*2.0*np.pi)+(2.5+(-2.25))*0.5
    ym = lambda y: (1.75-(-2.5))*0.5*np.sin(y*2.0*np.pi)+(1.75+(-2.5))*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p1.true_value = 2.0052938
    p1.bits = 8
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    p2.true_value = 1.1944509
    p2.bits = 8
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_subject(_subject_townsend_modified)
    return f
    
def _townsend_modified( parameters ):
    x = [p.value for p in parameters]    
    return -(np.cos((x[0]-0.1)*x[1])**2.0) - x[0]*np.sin(3.0*x[0]+x[1])

def _subject_townsend_modified( parameters ):
    x = [p.value for p in parameters]
    t = np.arctan2(x[0],x[1])
    l = x[0]**2.0 + x[1]**2.0
    r1 = (2.0*np.cos(t)-0.5*np.cos(2.0*t)-0.25*np.cos(3.0*t)-0.125*np.cos(4.0*t))**2.0
    r2 = 2.0 * np.sin(t)**2.0
    return l < (r1+r2)


def binh_korn():
    f1 = ObjectiveFunction( _binh_korna )
    f2 = ObjectiveFunction( _binh_kornb )
    xst = lambda x:(x>=0 and x<=5)
    yst = lambda y:(y>=0 and y<=3)
    xm = lambda x: (5.0-0)*0.5*np.sin(x*2.0*np.pi)+(5.0+0)*0.5
    ym = lambda y: (3.0-0)*0.5*np.sin(y*2.0*np.pi)+(3.0+0)*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    f1.add_parameter(p1)
    f1.add_parameter(p2)
    f2.add_parameter(p1)
    f2.add_parameter(p2)
    f = MultiObjectiveFunction()
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_objective_function(f1)
    f.add_objective_function(f2)
    f.add_subject(_subject_binh_korna)
    f.add_subject(_subject_binh_kornb)
    return f
    
def _binh_korna( parameters ):
    x = [p.value for p in parameters]
    return 4.0*x[0]**2.0 + 4.0*x[1]**2.0

def _binh_kornb( parameters ):
    x = [p.value for p in parameters]
    return (x[0]-5.0)**2.0 + (x[1]-5.0)**2.0

def _subject_binh_korna( parameters ):
    x = [p.value for p in parameters]
    return ((x[0]-5.0) + x[1]**2.0) <= 25
    
def _subject_binh_kornb( parameters ):
    x = [p.value for p in parameters]
    return ((x[0]-8.0)**2.0 + (x[1]+3.0)**2.0) >= 7.7
    

# function 1 - 7 (and exam_soo, exam_moo) from lecture 
# Stochastic Optimization Methods (P. Alotto)
def function1():
    f1 = ObjectiveFunction( _function1a )
    f2 = ObjectiveFunction( _function1b )    
    xst = lambda x:(x>=3 and x<=6)
    yst = lambda y:(y>=4 and y<=7)
    xm = lambda x: (6.0-3.0)*0.5*np.sin(x*2.0*np.pi)+(6.0+3.0)*0.5
    ym = lambda y: (7.0-4.0)*0.5*np.sin(y*2.0*np.pi)+(7.0+4.0)*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    f1.add_parameter(p1)
    f1.add_parameter(p2)
    f2.add_parameter(p1)
    f2.add_parameter(p2)
    f = MultiObjectiveFunction()
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_objective_function(f1)
    f.add_objective_function(f2)
    return f

def _function1a( parameters ):
    x = [p.value for p in parameters]
    return -np.pi*x[0]**2.0*x[1]
    
def _function1b( parameters ):
    x = [p.value for p in parameters]
    return np.pi*x[0]**2.0 + 2.0*np.pi*x[0]*x[1]

def function2():
    f1 = ObjectiveFunction( _function2a )
    f2 = ObjectiveFunction( _function2b )    
    xst = lambda x:(x>=0 and x<=1)
    yst = lambda y:(y>=0 and y<=1)
    xm = lambda x: abs(np.sin(x*2.0*np.pi))
    ym = lambda y: abs(np.sin(y*2.0*np.pi))
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    f1.add_parameter(p1)
    f1.add_parameter(p2)
    f2.add_parameter(p1)
    f2.add_parameter(p2)
    f = MultiObjectiveFunction()
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_objective_function(f1)
    f.add_objective_function(f2)
    f.add_subject(_subject2a)
    return f

def _function2a( parameters ):
    x = [p.value for p in parameters]
    return np.sqrt(x[0]**2.0+x[1]**2.0)
    
def _function2b( parameters ):
    x = [p.value for p in parameters]    
    return 1.0/np.sqrt(2.0) * (1.0-x[0]-x[1])
    
def _subject2a( parameters ):
    x = [p.value for p in parameters]
    return x[1] < (1.0 - x[0])
    
def function3():
    f1 = ObjectiveFunction( _function3a )
    f2 = ObjectiveFunction( _function3b )    
    xll = 0
    xul = 10.0
    yll = -10.0
    yul = 10.0
    xst = lambda x:(x>=xll and x<=xul)
    yst = lambda y:(y>=yll and y<=yul)
    xm = lambda x: (xul-xll)*0.5*np.sin(x*2.0*np.pi)+(xul+xll)*0.5
    ym = lambda y: (yul-yll)*0.5*np.sin(y*2.0*np.pi)+(yul+yll)*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    f1.add_parameter(p1)
    f1.add_parameter(p2)
    f2.add_parameter(p1)
    f2.add_parameter(p2)
    f = MultiObjectiveFunction()
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_objective_function(f1)
    f.add_objective_function(f2)
    return f

def _function3a( parameters ):
    x = [p.value for p in parameters]
    Ri = 1.0
    return x[0]/(Ri+x[0])
    
def _function3b( parameters ):
    x = [p.value for p in parameters]
    V = 1.0
    Ri = 1.0
    Xi = 1.0    
    return V**2-0 * x[0]/((Ri+x[0])**2.0+(Xi+x[1])**2.0)
        
def function4():
    f1 = ObjectiveFunction( _function4a )
    f2 = ObjectiveFunction( _function4b )
    p1 = ObjectiveFunctionParameter()
    p2 = ObjectiveFunctionParameter()
    f1.add_parameter(p1)
    f1.add_parameter(p2)
    f2.add_parameter(p1)
    f2.add_parameter(p2)
    f = MultiObjectiveFunction()
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_objective_function(f1)
    f.add_objective_function(f2)
    return f

def _function4a( parameters ):
    x = [p.value for p in parameters]
    return x[0]
    
def _function4b( parameters ):
    x = [p.value for p in parameters]
    n = len(x)
    g = 1.0 + 9.0/(n-1.0) * np.sum(x[1:])
    f1 =  _function4a( parameters )
    t = (1.0-np.sqrt(f1/g)-f1/g*np.sin(10.0*np.pi*f1))
    return g*t   
    
def function5():
    f1 = ObjectiveFunction( _function5a )
    f2 = ObjectiveFunction( _function5b )
    xst = lambda x:(x>=0 and x<=np.pi)
    yst = lambda y:(y>=0 and y<=np.pi)
    xm = lambda x: np.pi*abs(np.sin(x*2.0*np.pi))
    ym = lambda y: np.pi*abs(np.sin(y*2.0*np.pi))
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    f1.add_parameter(p1)
    f1.add_parameter(p2)
    f2.add_parameter(p1)
    f2.add_parameter(p2)
    f = MultiObjectiveFunction()
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_objective_function(f1)
    f.add_objective_function(f2)
    f.add_subject(_subject5a)
    f.add_subject(_subject5b)
    return f
    
def _function5a( parameters ):
    x = [p.value for p in parameters]
    return x[0]

def _function5b( parameters ):
    x = [p.value for p in parameters]
    return x[1]

def _subject5a( parameters ):
    x = [p.value for p in parameters]
    return (x[0]-0.5)**2.0 - 5.0*(x[1]-0.5)**2.0 <= 0
    
def _subject5b( parameters ):
    x = [p.value for p in parameters]
    c = -(x[0]**2.0+x[1]**2.0) + 1.0 + 0.1 * np.cos(16.0*np.arctan2(x[0],x[1]))
    return c <= 0

def function6():
    f1 = ObjectiveFunction( _function6a )
    f2 = ObjectiveFunction( _function6b )
    xst = lambda x:(x>=0 and x<=1)
    yst = lambda y:(y>=0 and y<=1)
    xm = lambda x: abs(np.sin(x*2.0*np.pi))
    ym = lambda y: abs(np.sin(y*2.0*np.pi))
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    f1.add_parameter(p1)
    f1.add_parameter(p2)
    f2.add_parameter(p1)
    f2.add_parameter(p2)
    f = MultiObjectiveFunction()
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_objective_function(f1)
    f.add_objective_function(f2)
    return f

def _function6a( parameters ):
    x = [p.value for p in parameters]
    return np.sin(0.5*np.pi*x[0])

def _function6b( parameters ):
    x = [p.value for p in parameters]
    num = (1.0-np.exp(-(x[1]-0.1)**2.0/(0.0001)))+(1.0-0.5*np.exp(-(x[1]-0.8)**2.0/(0.8)))
    denom = np.arctan(100.0*x[0])
    return num/denom

def function7():
    f1 = ObjectiveFunction( _function7a )
    f2 = ObjectiveFunction( _function7b )
    f3 = ObjectiveFunction( _function7c )
    xst = lambda x:(x>=-5 and x<=5)
    yst = lambda y:(y>=-5 and y<=5)
    xm = lambda x: 5.0*np.sin(x*2.0*np.pi)
    ym = lambda y: 5.0*np.sin(y*2.0*np.pi)
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    f1.add_parameter(p1)
    f1.add_parameter(p2)
    f2.add_parameter(p1)
    f2.add_parameter(p2)
    f3.add_parameter(p1)
    f3.add_parameter(p2)
    f = MultiObjectiveFunction()
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_objective_function(f1)
    f.add_objective_function(f2)
    f.add_objective_function(f3)
    return f
    
def _function7a( parameters ):
    x = [p.value for p in parameters]
    return 0.5 * (x[0]**2.0+x[1]**2.0) + np.sin(x[0]**2.0+x[1]**2.0)
    
def _function7b( parameters ):
    x = [p.value for p in parameters]
    t1 = ((3.0*x[0]-2.0*x[1]+4.0)**2.0)/8.0
    t2 = ((x[0]-x[1]+1.0)**2.0)/27.0
    return t1 + t2 + 15.0
    
def _function7c( parameters ):
    x = [p.value for p in parameters]
    t1 = 1.0/(x[0]**2.0+x[1]**2.0+1.0)
    t2 = 1.1 * np.exp(-x[0]**2.0-x[1]**2.0)
    return t1 - t2
    
def exam_soo():
    f = ObjectiveFunction( _exam_soo )
    xll = -4.5
    xul = 4.5
    xst = lambda x:(x>=xll and x<=xul)
    xm = lambda x: (xul-xll)*0.5*np.sin(x*2.0*np.pi)+(xul+xll)*0.5
    yll = -4.5
    yul = 4.5
    yst = lambda y:(y>=yll and y<=yul)
    ym = lambda y: (yul-yll)*0.5*np.sin(y*2.0*np.pi)+(yul+yll)*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p1.bits = 4
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    p2.bits = 4
    f.add_parameter(p1)
    f.add_parameter(p2)
    return f

def _exam_soo( parameters ):
    x = [p.value for p in parameters]
    return (x[0]**2.0+x[0])*np.cos(x[0]) 
    
def exam_moo():
    f1 = ObjectiveFunction( _exam_moo1 )
    f2 = ObjectiveFunction( _exam_moo2 )
    xll = -4.0
    xul = 4.0
    xst = lambda x:(x>=xll and x<=xul)
    xm = lambda x: (xul-xll)*0.5*np.sin(x*2.0*np.pi)+(xul+xll)*0.5
    yll = -4.0
    yul = 4.0
    yst = lambda y:(y>=yll and y<=yul)
    ym = lambda y: (yul-yll)*0.5*np.sin(y*2.0*np.pi)+(yul+yll)*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p1.bits = 4
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    p2.bits = 4
    f1.add_parameter(p1)
    f1.add_parameter(p2)
    f2.add_parameter(p1)
    f2.add_parameter(p2)
    f = MultiObjectiveFunction()
    f.add_parameter(p1)
    f.add_parameter(p2)
    f.add_objective_function(f1)
    f.add_objective_function(f2)
    return f

def _exam_moo1( parameters ):
    x = np.array([p.value for p in parameters])
    return (1.0 - np.exp(-np.sum((x-1.0/np.sqrt(2.0))**2.0)))
    
def _exam_moo2( parameters ):
    x = np.array([p.value for p in parameters])
    return (1.0 - np.exp(-np.sum((x+1.0/np.sqrt(2.0))**2.0)))

#F1-F16 from textbook Practival Genetic Algorithms (Haupt & Haupt)
def F6():
    f = ObjectiveFunction( _F6 )
    f.global_minimum = -100.22
    xll = -10.0
    xul = 10.0
    xst = lambda x:(x>=xll and x<=xul)
    xm = lambda x: (xul-xll)*0.5*np.sin(x*2.0*np.pi)+(xul+xll)*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p1.true_value = 9.6204
    p1.bits = 4
    f.add_parameter(p1)
    return f

def _F6( parameters ):
    x = [p.value for p in parameters]
    return (x[0]**2.0+x[0])*np.cos(x[0]) 
    
def F7():
    f = ObjectiveFunction( _F7 )
    f.global_minimum = -18.5547
    xll = 0
    xul = 10.0
    yll = 0
    yul = 10.0
    xst = lambda x:(x>=xll and x<=xul)
    yst = lambda y:(y>=yll and y<=yul)
    xm = lambda x: (xul-xll)*0.5*np.sin(x*2.0*np.pi)+(xul+xll)*0.5
    ym = lambda y: (yul-yll)*0.5*np.sin(y*2.0*np.pi)+(yul+yll)*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p1.true_value = 0.9039
    p1.bits = 4
    f.add_parameter(p1)
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    p2.true_value = 0.8668
    p2.bits = 4
    f.add_parameter(p2)
    return f

def _F7( parameters ):
    x = [p.value for p in parameters]
    return (x[0]*np.sin(4.0*x[0])+1.1*x[1]*np.sin(2.0*x[1]))
    
    
    
# This list is due www.sfu.ca/~ssurjano/optimization.html
def ackley( dimensions=2 ):
    f = ObjectiveFunction( _ackley )
    f.global_minimum = 0
    for _ in list(range(dimensions)):
        xll = -32.768
        xul = 32.768
        xst = lambda x:(x>=xll and x<=xul)
        xm = lambda x: (xul-xll)*0.5*np.sin(x*2.0*np.pi)+(xul+xll)*0.5
        p = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
        p.true_value = 0
        p.bits = 4
        f.add_parameter(p)
    return f
    
def _ackley( parameters ):
    x = np.array([p.value for p in parameters])
    a = 20.0
    b = 0.2
    c = 2.0*np.pi
    d = len(x)
    f1 = -a * np.exp( -b * np.sqrt( 1.0/d * np.sum(x**2.0) ))
    f2 = -np.exp(1.0/d * np.sum(np.cos(c*x))) + a + np.exp(1.0)
    return f1+f2
    
def bukin6():
    f = ObjectiveFunction( _bukin6 )
    f.global_minimum = 0
    xll = -15.0
    xul = -5.0
    yll = -3.0
    yul = 3.0
    xst = lambda x:(x>=xll and x<=xul)
    yst = lambda y:(y>=yll and y<=yul)
    xm = lambda x: (xul-xll)*0.5*np.sin(x*2.0*np.pi)+(xul+xll)*0.5
    ym = lambda y: (yul-yll)*0.5*np.sin(y*2.0*np.pi)+(yul+yll)*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p1.true_value = -10.0
    p1.bits = 4
    f.add_parameter(p1)
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    p2.true_value = 1.0
    p2.bits = 4
    f.add_parameter(p2)
    return f

def _bukin6( parameters ):
    x = [p.value for p in parameters]
    return 100.0 * np.sqrt( np.fabs(x[1] - 0.01*x[0]**2.0) + 0.01 * np.fabs(x[0] + 10.0))

def eggholder():
    f = ObjectiveFunction( _eggholder )
    f.global_minimum = -959.6407
    xll = -512.0
    xul = 512.0
    yll = -512.0
    yul = 512.0
    xst = lambda x:(x>=xll and x<=xul)
    yst = lambda y:(y>=yll and y<=yul)
    xm = lambda x: (xul-xll)*0.5*np.sin(x*2.0*np.pi)+(xul+xll)*0.5
    ym = lambda y: (yul-yll)*0.5*np.sin(y*2.0*np.pi)+(yul+yll)*0.5
    p1 = ObjectiveFunctionParameter(subject_to=xst,mapping=xm)
    p1.true_value = 512.0
    p1.bits = 4
    f.add_parameter(p1)
    p2 = ObjectiveFunctionParameter(subject_to=yst,mapping=ym)
    p2.true_value = 404.2319
    p2.bits = 4
    f.add_parameter(p2)
    return f

def _eggholder( parameters ):
    x = [p.value for p in parameters]
    return -(x[1]+47.0) * np.sin(np.sqrt(np.fabs(x[1]+0.5*x[0]+47.0))) - x[0] * np.sin(np.sqrt(np.fabs(x[0]-x[1]-47.0)))


  
    
    
    