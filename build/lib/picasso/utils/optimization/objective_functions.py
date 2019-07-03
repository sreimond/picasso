# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 11:42:57 2018

@author: sreimond
"""

import numpy as np
import warnings

class ObjectiveFunction( object ):
    """ 
    The `ObjectiveFunction` class defines an objective function.
    """   
    def __init__(self,definition=None):
        self.definition = definition
        self.parameters = []
        self.global_minimum = None
        self.subject_to = []
    
    def add_parameter(self,parameter):
        if not isinstance(parameter,ObjectiveFunctionParameter):
            warnings.warn('Please use the ObjectiveFunctionParameter class.')
            return
        self.parameters.append(parameter)
        self.update_dimension()
    
    def add_subject(self,subject):
        if not callable(subject):
            warnings.warn('Subject must be a callable function returning a boolean.')
            return
        self.subject_to.append(subject)
        
    def update_dimension(self):
        self.dimension = len(self.parameters)
        
    def evaluate(self,parameters=[]):
        if not parameters:
            parameters = self.parameters
        return self.definition( parameters )
        
    def determine_feasibility(self,parameters=[]):
        if not parameters:
            parameters = self.parameters
        feasibility = []
        for parameter in parameters:
            feasibility.append( parameter.is_feasible() )
        for subject in self.subject_to:
            feasibility.append( subject(parameters) )
        return np.mean(feasibility)
        
             
class MultiObjectiveFunction( ObjectiveFunction ):
    """
    The `MultiObjectiveFunction` class defines various objective functions.
    """
    def __init__(self):
        super(MultiObjectiveFunction, self).__init__()
        self.objective_functions = []
    
    def add_objective_function(self,objective_function):
        if not isinstance(objective_function,ObjectiveFunction):
            warnings.warn('Please use the ObjectiveFunction class.')
            return
        self.objective_functions.append(objective_function)
        self.update_objective_dimension()
        
    def update_objective_dimension(self):
        self.objective_dimension = len(self.objective_functions)
        
    def evaluate(self,parameters=[]):
        if not parameters:
            parameters = self.parameters
        f = []
        for objective_function in self.objective_functions:
            f.append( objective_function.evaluate(parameters=parameters) )
        return f
        
    def determine_feasibility(self,parameters=[]):
        if not parameters:
            parameters = self.parameters
        feasibility = []
        for parameter in parameters:
            feasibility.append( parameter.is_feasible() )
        for subject in self.subject_to:
            feasibility.append( subject(parameters) )
        feasibility = [np.mean(feasibility)]
        for objective_function in self.objective_functions:
            feasibility.append( objective_function.determine_feasibility(parameters=parameters) )
        return np.mean(feasibility)
        
        
class ObjectiveFunctionParameter( object ):
    """ 
    The `ObjectiveFunctionParameter` class defines a parameter for 
    the ObjectiveFunction.
    """   
    def __init__(self,
                 subject_to=lambda x:(x>-np.inf and x<np.inf),
                 mapping=lambda x:x):
        self.subject_to = subject_to
        self.mapping = mapping
        self.value = None
        self.mapping_parameter = None # between 0 and 1 (for binary optimzers)
        self.true_value = None
        self.bits = 8
        self.code = None
        self.label = ''
        
    def set_random_value(self):
        self.mapping_parameter = np.random.random()
        self.value = self.mapping( self.mapping_parameter )    
        self.update_binary_representation()
    
    def set_mapping_parameter(self,value):
        self.mapping_parameter = value
        self.value = self.mapping( self.mapping_parameter )    
        self.update_binary_representation()
        
    def set_binary_code(self,code):
        self.code = code
        self.update_values()
    
    def update_values(self):
        # https://de.mathworks.com/matlabcentral/answers/25549-convert-floating-point-to-binary
        self.mapping_parameter = np.dot(self.code,[2.0**p for p in np.arange(-1,-(self.bits+1),-1)])
        self.value = self.mapping( self.mapping_parameter )
    
    def update_binary_representation(self):
        # https://de.mathworks.com/matlabcentral/answers/25549-convert-floating-point-to-binary
        self.code = [np.fix(np.fmod(self.mapping_parameter * 2.0**p, 2)) for p in np.arange(1,self.bits+1)]
    
    def is_feasible(self):
        return self.subject_to( self.value )
        
        
        
        
        
        
