# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 14:43:34 2018

@author: sreimond
"""

import numpy as np
from copy import deepcopy


class Simplex( object ):
    """ 
    The `Simplex` class holds vertices and optimizes the objective 
    function.
    """    
    def __init__( self, objective_function, initial_parameters=[] ):
        self.rho = 1.0
        self.chi = 2.0
        self.gamma = 0.5
        self.sigma = 0.5
        self.max_iterations = 10
        self.development = {'iteration':[],
                            'function_calls':0,
                            'min_cost':[],
                            'max_cost':[],
                            'mean_cost':[],
                            'std_cost':[]}
        self.objective_function = objective_function
        self.initial_parameters = initial_parameters
        self.vertices = []
        self.optimized_parameters = []
        self.iteration = 0
    
    def initialize_vertices(self):
        for ix in list(range(self.objective_function.dimension+1)):
            self.add_vertex()
            if ix < len(self.initial_parameters):
                self.vertices[ix].parameters = deepcopy( self.initial_parameters[ix] )
                self.vertices[ix].update_coordinates()
            
    def add_vertex(self):
        vertex = Vertex(self)
        vertex.initialize()
        self.vertices.append(vertex)
        
    def sort_vertices(self):
        cost = [vertex.cost for vertex in self.vertices]
        neg_feasibility = [-vertex.feasibility for vertex in self.vertices]     
        ix = np.lexsort((cost,neg_feasibility))
        return [self.vertices[jx] for jx in ix]
        
    def update_statistics(self):        
        cost = [vertex.cost for vertex in self.vertices]
        self.development['iteration'].append( self.iteration )
        self.development['min_cost'].append( np.min(cost) )
        self.development['max_cost'].append( np.max(cost) )
        self.development['mean_cost'].append( np.mean(cost) )
        self.development['std_cost'].append( np.std(cost) )    
    
    def test_vertices(self):
        for vertex in self.vertices:
            vertex.determine_cost()
            self.development['function_calls'] += 1
    
    def determine_centroid(self):
        self.centroid = Vertex(self)
        for ix,parameter in enumerate(self.centroid.parameters):
            x = []
            for vertex in self.vertices[0:-1]:
                x.append( vertex.coordinates[ix] )
            parameter.set_mapping_parameter( np.mean(x) ) 
        self.centroid.update_coordinates()
        
    def determine_reflection(self):
        x_mean = np.array(self.centroid.coordinates)
        x_n1 = np.array(self.vertices[-1].coordinates)
        coordinates = (1.0+self.rho)*x_mean - self.rho*x_n1
        self.reflection = Vertex(self)
        self.reflection.coordinates = [c for c in coordinates]
        self.reflection.update_parameters()
        self.reflection.determine_cost()
        self.development['function_calls'] += 1
        
    def determine_expansion(self):
        x_mean = np.array(self.centroid.coordinates)
        x_n1 = np.array(self.vertices[-1].coordinates)
        coordinates = (1.0+self.rho*self.chi)*x_mean - self.rho*self.chi*x_n1
        self.expansion = Vertex(self)
        self.expansion.coordinates = [c for c in coordinates]
        self.expansion.update_parameters()
        self.expansion.determine_cost()
        self.development['function_calls'] += 1
        
    def determine_contraction_outside(self):
        x_mean = np.array(self.centroid.coordinates)
        x_n1 = np.array(self.vertices[-1].coordinates)
        coordinates = (1.0+self.rho*self.gamma)*x_mean - self.rho*self.gamma*x_n1
        self.contraction_outside = Vertex(self)
        self.contraction_outside.coordinates = [c for c in coordinates]
        self.contraction_outside.update_parameters()
        self.contraction_outside.determine_cost()
        self.development['function_calls'] += 1
        
    def determine_contraction_inside(self):
        x_mean = np.array(self.centroid.coordinates)
        x_n1 = np.array(self.vertices[-1].coordinates)
        coordinates = (1.0-self.gamma)*x_mean - self.gamma*x_n1
        self.contraction_inside = Vertex(self)
        self.contraction_inside.coordinates = [c for c in coordinates]
        self.contraction_inside.update_parameters()
        self.contraction_inside.determine_cost()
        self.development['function_calls'] += 1
    
    def perform_shrink(self):
        ix_shrink = list(range(1,len(self.vertices)))
        x_1 = np.array(self.vertices[0].coordinates)
        for ix in ix_shrink:
            x_i = np.array(self.vertices[ix].coordinates)
            coordinates = x_1 + self.sigma*(x_i-x_1)
            self.vertices[ix] = Vertex(self)
            self.vertices[ix].coordinates = [c for c in coordinates]
            self.vertices[ix].update_parameters()
            self.vertices[ix].determine_cost()
            self.development['function_calls'] += 1
    
    def start_search(self):        
        self.initialize_vertices()
        self.test_vertices()
        self.vertices = self.sort_vertices()
        self.update_statistics()
        for _ in range(self.max_iterations-1):
            self.iteration += 1
            self.determine_centroid()
            self.determine_reflection()
            f1 = self.vertices[0].cost
            fr = self.reflection.cost
            fn = self.vertices[-2].cost
            fn1 = self.vertices[-1].cost
            if (f1<=fr) and (fr<fn):
                del self.vertices[-1]
                self.vertices.append(self.reflection)
            elif (fr<f1):
                self.determine_expansion()
                fe = self.expansion.cost
                del self.vertices[-1]
                if (fe<fr):                    
                    self.vertices.append(self.expansion)
                else:
                    self.vertices.append(self.reflection)
            elif (fn<=fr) and (fr<fn1):
                self.determine_contraction_outside()
                fc = self.contraction_outside.cost
                if (fc<=fr):
                    del self.vertices[-1]
                    self.vertices.append(self.contraction_outside)
                else:
                    self.perform_shrink()
            elif (fr>fn1):
                self.determine_contraction_inside()
                fcc = self.contraction_inside.cost
                if (fcc<fn1):
                    del self.vertices[-1]
                    self.vertices.append(self.contraction_inside)
                else:
                    self.perform_shrink()
            self.vertices = self.sort_vertices()
            self.update_statistics()
        self.optimized_parameters = self.vertices[0].parameters

            
class Vertex( object ):
    """ 
    The `Vertex` class holds parameters and values.
    """    
    def __init__( self, simplex ):
        self.simplex = simplex
        self.cost = None
        self.feasibility = None
        self.initialize()
        
    def initialize(self):
        self.parameters = []        
        for p in self.simplex.objective_function.parameters:
            parameter = deepcopy( p )
            parameter.set_random_value()
            self.parameters.append( parameter )
        self.update_coordinates()
        
    def update_coordinates(self):
        self.coordinates = np.array([parameter.mapping_parameter for parameter in self.parameters])
        
    def update_parameters(self):
        for ix,parameter in enumerate(self.parameters):
            value = self.coordinates[ix]
            parameter.set_mapping_parameter(value)
        
    def determine_cost(self):
        self.cost = self.simplex.objective_function.evaluate( self.parameters )
        self.feasibility = self.simplex.objective_function.determine_feasibility( self.parameters )
        
        
        
        
        
        