# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 10:57:10 2018

@author: sreimond
"""

import matplotlib
matplotlib.use('agg')
import sys
sys.path.insert(0,'../optimization')
import numpy as np
from copy import deepcopy
import pylab
from mpl_toolkits.mplot3d import Axes3D
import warnings
#import SOBGA, SOCGA, MOCGA, downhill_simplex
from ..optimization import SOBGA
from ..optimization import SOCGA
from ..optimization import MOCGA
from ..optimization import downhill_simplex

def plot_convergence( optimization_result ):
    """
    The function `plot_convergence` function inside the optimization_plots 
    module visualizes the convergence of a GA result as obtained via the 
    optimization.SOBGA or optimization.SOCGA modules.
    """
    if (isinstance(optimization_result,downhill_simplex.Simplex)):
        return _plot_convergence_ds( optimization_result )
    if not (isinstance(optimization_result,SOBGA.Population) or isinstance(optimization_result,SOCGA.Population)):
        warnings.warn('Please use the SOBGA or SOCGA class.')
        return
    x = optimization_result.development['generation']
    y_best = optimization_result.development['best_fitness']
    y_mean = optimization_result.development['mean_fitness']    
    # settings
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    pylab.rcParams.update(params)
    fig,ax = pylab.subplots()
    ax.plot(x,y_best,'-k',label='max. fitness')
    ax.plot(x,y_mean,'--k',label='mean fitness')
    ax.set_xlabel('generation')
    ax.set_ylabel('fitness')
    ax.grid()
    ax.legend()
    return fig, ax

def _plot_convergence_ds( optimization_result ):
    """
    The function `_plot_convergence_ds` function inside the optimization_plots 
    module visualizes the convergence of a Downhill Simplex result as obtained 
    via the optimization.downhill_simplex module.
    """
    x = optimization_result.development['iteration']
    y_best = optimization_result.development['min_cost']
    y_mean = optimization_result.development['mean_cost']    
    # settings
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    pylab.rcParams.update(params)
    fig,ax = pylab.subplots()
    ax.plot(x,y_best,'-k',label='min. cost')
    ax.plot(x,y_mean,'--k',label='mean cost')
    ax.set_xlabel('iteration')
    ax.set_ylabel('cost')
    ax.grid()
    ax.legend()
    return fig, ax
    

def plot_objective_function( optimization_result ):
    """
    The function `plot_objective_function` function inside the optimization_plots 
    module visualizes 1D line or the 2D cost function and the final population 
    of a GA result as obtained via the optimization.SOBGA or optimization.SOCGA
    modules.
    """    
    if (isinstance(optimization_result,downhill_simplex.Simplex)):
        return _plot_objective_function_ds( optimization_result )
    if not (isinstance(optimization_result,SOBGA.Population) or isinstance(optimization_result,SOCGA.Population)):
        warnings.warn('Please use the SOBGA or SOCGA class.')
        return
    f = optimization_result.config['objective_function']
    if f.dimension > 2:
        warnings.warn('Only 1D or 2D objective functions are supported.')
        return    
    elif f.dimension == 1:
        return _plot_objective_function1D( optimization_result )
    elif f.dimension == 2:
        return _plot_objective_function2D( optimization_result )
    
def _plot_objective_function1D( optimization_result ):
    f = optimization_result.config['objective_function']
    t = np.linspace(-1.0,1.0,num=100)
    parameters = []
    for parameter in f.parameters:
        parameters.append( deepcopy(parameter) )
    x = np.zeros(len(t))
    for ix,ti in enumerate(t):
        parameters[0].set_mapping_parameter( ti )
        x[ix] = parameters[0].value
    ix = np.argsort(x)
    xs = [x[jx] for jx in ix]
    ts = [t[jx] for jx in ix]
    ys = np.ones(len(t)) * np.nan
    for ix,ti in enumerate(ts):
        parameters[0].set_mapping_parameter( ti )
        value = f.evaluate(parameters)
        if f.determine_feasibility(parameters) == 1:
            ys[ix] = value
    x_pop = []
    y_pop = []
    for individual in optimization_result.individuals:
        x_pop.append( individual.genes[0].value )
        value = f.evaluate( individual.genes )
        y_pop.append(value)
    # settings
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    pylab.rcParams.update(params)
    fig,ax = pylab.subplots()    
    ax.plot(xs,ys,'-k',label='f(x)')
    ax.plot(x_pop,y_pop,'xk',label='final population')
    ax.plot(x_pop[0],y_pop[0],'or',label='best individual')  
    x_true = parameters[0].true_value 
    y_true = f.global_minimum
    if (x_true is not None) and (y_true is not None):
        ax.plot(x_true,y_true,'*b',label='global minimum')  
    ax.set_xlabel('x')
    ax.set_ylabel('f(x)')
    ax.grid()
    ax.legend()
    return fig, ax
            
def _plot_objective_function2D( optimization_result ):
    f = optimization_result.config['objective_function']
    tx = np.linspace(-1.0,1.0,num=100)
    ty = np.linspace(-1.0,1.0,num=100)
    parameters = []
    for parameter in f.parameters:
        parameters.append( deepcopy(parameter) )
    x = np.zeros(len(tx))
    y = np.zeros(len(ty))
    for ix,_ in enumerate(tx):
        parameters[0].set_mapping_parameter( tx[ix] )
        x[ix] = parameters[0].value
        parameters[1].set_mapping_parameter( ty[ix] )
        y[ix] = parameters[0].value
    ix = np.argsort(x)
    txs = [tx[jx] for jx in ix]
    ix = np.argsort(y)
    tys = [ty[jx] for jx in ix]
    txg,tyg = np.meshgrid(txs,tys)
    x = np.zeros(txg.shape)
    y = np.zeros(txg.shape)
    c = np.ones(txg.shape)*np.nan
    for ix in list(range(txg.shape[0])):
        for iy in list(range(txg.shape[1])):
            parameters[0].set_mapping_parameter( txg[ix,iy] )
            parameters[1].set_mapping_parameter( tyg[ix,iy])
            x[ix,iy] = parameters[0].value
            y[ix,iy] = parameters[1].value
            value = f.evaluate(parameters)
            if f.determine_feasibility(parameters) == 1:
                c[ix,iy] = value   
    x_pop = []
    y_pop = []
    c_pop = []
    for individual in optimization_result.individuals:
        x_pop.append( individual.genes[0].value )
        y_pop.append( individual.genes[1].value )
        value = f.evaluate( individual.genes )
        c_pop.append(value)    
    # settings
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    pylab.rcParams.update(params)
    fig,ax = pylab.subplots()    
    ax.contourf(x,y,c,cmap=pylab.cm.rainbow)
#    ax.imshow(c)
#    ax.pcolor(x,y,c)
#    ax.scatter(x,y,c=c)
    ax.plot(x_pop,y_pop,'xk',label='final population')
    ax.plot(x_pop[0],y_pop[0],'or',label='best individual')
    x_true = parameters[0].true_value 
    y_true = parameters[1].true_value 
    if (x_true is not None) and (y_true is not None):
        ax.plot(x_true,y_true,'*b',label='global minimum')
    ax.set_xlabel('x')
    ax.set_ylabel('y')    
    ax.legend()
    return fig, ax


def _plot_objective_function_ds( optimization_result ):
    """
    The function `_plot_objective_function_ds` function inside the optimization_plots 
    module visualizes 1D line or the 2D cost function and the final population 
    of a DS result as obtained via the optimization.downhill_simplex module.
    """    
    f = optimization_result.objective_function
    if f.dimension > 2:
        warnings.warn('Only 1D or 2D objective functions are supported.')
        return    
    elif f.dimension == 1:
        return _plot_objective_function_ds_1D( optimization_result )
    elif f.dimension == 2:
        return _plot_objective_function_ds_2D( optimization_result )
    
def _plot_objective_function_ds_1D( optimization_result ):
    f = optimization_result.objective_function
    t = np.linspace(-1.0,1.0,num=100)
    parameters = []
    for parameter in f.parameters:
        parameters.append( deepcopy(parameter) )
    x = np.zeros(len(t))
    for ix,ti in enumerate(t):
        parameters[0].set_mapping_parameter( ti )
        x[ix] = parameters[0].value
    ix = np.argsort(x)
    xs = [x[jx] for jx in ix]
    ts = [t[jx] for jx in ix]
    ys = np.ones(len(t)) * np.nan
    for ix,ti in enumerate(ts):
        parameters[0].set_mapping_parameter( ti )
        value = f.evaluate(parameters)
        if f.determine_feasibility(parameters) == 1:
            ys[ix] = value
    x_pop = []
    y_pop = []
    for vertex in optimization_result.vertices:
        x_pop.append( vertex.parameters[0].value )
        value = f.evaluate( vertex.parameters )
        y_pop.append(value)
    # settings
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    pylab.rcParams.update(params)
    fig,ax = pylab.subplots()    
    ax.plot(xs,ys,'-k',label='f(x)')
    ax.plot(x_pop,y_pop,'xk',label='final simplex')
    ax.plot(x_pop[0],y_pop[0],'or',label='best vertex')  
    x_true = parameters[0].true_value 
    y_true = f.global_minimum
    if (x_true is not None) and (y_true is not None):
        ax.plot(x_true,y_true,'*b',label='global minimum')  
    ax.set_xlabel('x')
    ax.set_ylabel('f(x)')
    ax.grid()
    ax.legend()
    return fig, ax
            
def _plot_objective_function_ds_2D( optimization_result ):
    f = optimization_result.objective_function
    tx = np.linspace(-1.0,1.0,num=100)
    ty = np.linspace(-1.0,1.0,num=100)
    parameters = []
    for parameter in f.parameters:
        parameters.append( deepcopy(parameter) )
    x = np.zeros(len(tx))
    y = np.zeros(len(ty))
    for ix,_ in enumerate(tx):
        parameters[0].set_mapping_parameter( tx[ix] )
        x[ix] = parameters[0].value
        parameters[1].set_mapping_parameter( ty[ix] )
        y[ix] = parameters[0].value
    ix = np.argsort(x)
    txs = [tx[jx] for jx in ix]
    ix = np.argsort(y)
    tys = [ty[jx] for jx in ix]
    txg,tyg = np.meshgrid(txs,tys)
    x = np.zeros(txg.shape)
    y = np.zeros(txg.shape)
    c = np.ones(txg.shape)*np.nan
    for ix in list(range(txg.shape[0])):
        for iy in list(range(txg.shape[1])):
            parameters[0].set_mapping_parameter( txg[ix,iy] )
            parameters[1].set_mapping_parameter( tyg[ix,iy])
            x[ix,iy] = parameters[0].value
            y[ix,iy] = parameters[1].value
            value = f.evaluate(parameters)
            if f.determine_feasibility(parameters) == 1:
                c[ix,iy] = value   
    x_pop = []
    y_pop = []
    c_pop = []
    for vertex in optimization_result.vertices:
        x_pop.append( vertex.parameters[0].value )
        y_pop.append( vertex.parameters[1].value )
        value = f.evaluate( vertex.parameters )
        c_pop.append(value)    
    # settings
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    pylab.rcParams.update(params)
    fig,ax = pylab.subplots()    
    ax.contourf(x,y,c,cmap=pylab.cm.rainbow)
#    ax.imshow(c)
#    ax.pcolor(x,y,c)
#    ax.scatter(x,y,c=c)
    ax.plot(x_pop,y_pop,'xk',label='final simplex')
    ax.plot(x_pop[0],y_pop[0],'or',label='best vertex')
    x_true = parameters[0].true_value 
    y_true = parameters[1].true_value 
    if (x_true is not None) and (y_true is not None):
        ax.plot(x_true,y_true,'*b',label='global minimum')
    ax.set_xlabel('x')
    ax.set_ylabel('y')    
    ax.legend()
    return fig, ax
    
    
def plot_pareto_front( optimization_result ):
    """
    The function `plot_pareto_front` function inside the optimization_plots 
    module visualizes the final population of a GA result as obtained via the 
    optimization.MOCGA module and highlights the obtained Pareto front.
    """
    if not (isinstance(optimization_result,MOCGA.Population)):
        warnings.warn('Please use the MOCGA class.')
    f = optimization_result.config['objective_function']
    if f.objective_dimension == 2:
        return _plot_pareto_front2D( optimization_result )
    elif f.objective_dimension == 3:
        return _plot_pareto_front3D( optimization_result )
    else:
        warnings.warn('Only multi objective functions with two or three objectives are supported.')
        return
       
def _plot_pareto_front2D( optimization_result ):
    f = optimization_result.config['objective_function'] 
    rank = [ind.rank for ind in optimization_result.archive]
    x = [-ind.fitness[0] for ind in optimization_result.archive]
    y = [-ind.fitness[1] for ind in optimization_result.archive]
    ix = np.lexsort((x,y,rank))
    rs = [rank[jx] for jx in ix]
    xs = [x[jx] for jx in ix]
    ys = [y[jx] for jx in ix]
    xs0 = [x[jx] for jx in ix if rank[jx]==0]
    ys0 = [y[jx] for jx in ix if rank[jx]==0]    
#    xx = [ind.genes[0].value for ind in optimization_result.archive]
#    xxs0 = [('%.2f' % xx[jx]) for jx in ix if rank[jx]==0]
    # settings
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    pylab.rcParams.update(params)
    fig,ax = pylab.subplots()    
    ax.plot(xs0,ys0,'-r',label='Pareto front')
#    for ix in list(range(len(xxs0))):
#        ax.text(xs0[ix],ys0[ix],xxs0[ix])
#    ax.plot(xs0,ys0,'og',label='Pareto front')
    ax.plot(xs,ys,'xk',label='final archive')
    ax.set_xlabel('f1')
    ax.set_ylabel('f2')
    ax.grid()
    ax.legend()
    return fig, ax
    
def _plot_pareto_front3D( optimization_result ):
    f = optimization_result.config['objective_function'] 
    rank = [ind.rank for ind in optimization_result.archive]
    x = [-ind.fitness[0] for ind in optimization_result.archive]
    y = [-ind.fitness[1] for ind in optimization_result.archive]
    z = [-ind.fitness[2] for ind in optimization_result.archive]
    ix = np.lexsort((x,y,z,rank))
    rs = [rank[jx] for jx in ix]
    xs = [x[jx] for jx in ix]
    ys = [y[jx] for jx in ix]
    zs = [z[jx] for jx in ix]
    xs0 = [x[jx] for jx in ix if rank[jx]==0]
    ys0 = [y[jx] for jx in ix if rank[jx]==0]    
    zs0 = [z[jx] for jx in ix if rank[jx]==0]    
    # settings
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    pylab.rcParams.update(params)
    fig = pylab.figure()    
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs0,ys0,zs0,'-r',label='Pareto front')
    ax.plot(xs,ys,zs,'xk',label='final archive')
    ax.set_xlabel('f1')
    ax.set_ylabel('f2')
    ax.set_zlabel('f3')
    ax.grid()
    ax.legend()
    return fig, ax
    
    
    
def plot_pareto_set( optimization_result ):
    """
    The function `plot_pareto_set` function inside the optimization_plots 
    module visualizes the final population of a GA result as obtained via the 
    optimization.MOCGA module and highlights the obtained Pareto set.
    """
    if not (isinstance(optimization_result,MOCGA.Population)):
        warnings.warn('Please use the MOCGA class.')
    f = optimization_result.config['objective_function']
    if f.dimension == 1:
        return _plot_pareto_set1D( optimization_result )
    elif f.dimension == 2:
        return _plot_pareto_set2D( optimization_result )
    elif f.dimension == 3:
        return _plot_pareto_set3D( optimization_result )
    else:
        warnings.warn('Only multi objective functions with one, two or three dimensions are supported.')
        return
        
def _plot_pareto_set1D( optimization_result ):
    f = optimization_result.config['objective_function'] 
    rank = [ind.rank for ind in optimization_result.archive]
    x = [ind.genes[0].value for ind in optimization_result.archive]
    y = [ind.genes[0].value for ind in optimization_result.archive]
    ix = np.lexsort((x,y,rank))
    rs = [rank[jx] for jx in ix]
    xs = [x[jx] for jx in ix]
    ys = [y[jx] for jx in ix]
    xs0 = [x[jx] for jx in ix if rank[jx]==0]
    ys0 = [y[jx] for jx in ix if rank[jx]==0]    
    # settings
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    pylab.rcParams.update(params)
    fig,ax = pylab.subplots()    
    ax.plot(xs0,ys0,'or',label='Pareto set')
    ax.plot(xs,ys,'xk',label='final archive')
    ax.set_xlabel('x')
    ax.set_ylabel('x')
    ax.grid()
    ax.legend()
    return fig, ax

def _plot_pareto_set2D( optimization_result ):
    f = optimization_result.config['objective_function'] 
    rank = [ind.rank for ind in optimization_result.archive]
    x = [ind.genes[0].value for ind in optimization_result.archive]
    y = [ind.genes[1].value for ind in optimization_result.archive]
    ix = np.lexsort((x,y,rank))
    rs = [rank[jx] for jx in ix]
    xs = [x[jx] for jx in ix]
    ys = [y[jx] for jx in ix]
    xs0 = [x[jx] for jx in ix if rank[jx]==0]
    ys0 = [y[jx] for jx in ix if rank[jx]==0]    
    # settings
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    pylab.rcParams.update(params)
    fig,ax = pylab.subplots()    
    ax.plot(xs0,ys0,'or',label='Pareto set')
    ax.plot(xs,ys,'xk',label='final archive')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid()
    ax.legend()
    return fig, ax
        
def _plot_pareto_set3D( optimization_result ):
    f = optimization_result.config['objective_function'] 
    rank = [ind.rank for ind in optimization_result.archive]
    x = [ind.genes[0].value for ind in optimization_result.archive]
    y = [ind.genes[1].value for ind in optimization_result.archive]
    z = [ind.genes[2].value for ind in optimization_result.archive]
    ix = np.lexsort((x,y,z,rank))
    rs = [rank[jx] for jx in ix]
    xs = [x[jx] for jx in ix]
    ys = [y[jx] for jx in ix]
    zs = [z[jx] for jx in ix]
    xs0 = [x[jx] for jx in ix if rank[jx]==0]
    ys0 = [y[jx] for jx in ix if rank[jx]==0]    
    zs0 = [z[jx] for jx in ix if rank[jx]==0]    
    # settings
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    pylab.rcParams.update(params)
    fig = pylab.figure()    
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs0,ys0,zs0,'-r',label='Pareto set')
    ax.plot(xs,ys,zs,'xk',label='final archive')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.grid()
    ax.legend()
    return fig, ax
