# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 14:03:18 2018

@author: sreimond
"""

from ..optimization import MOCGA
from ..optimization import SOBGA
from ..optimization import SOCGA
from ..optimization import downhill_simplex
import numpy as np
import os
import inspect
import warnings

def export_optimization( optimization_result, output_filename ):
    """
    The function `export_optimization` inside the optimization_files module 
    lets you export an instance of the SOBGA, SOCGA OR MOCGA class as a text 
    file.
    """    
    if os.path.exists( output_filename ):
        warnings.warn('Existing output file will be overwritten!')
    if isinstance(optimization_result,SOBGA.Population):
        _export_SOGA( optimization_result, output_filename, 'SOBGA' )
    elif isinstance(optimization_result,SOCGA.Population ):
        _export_SOGA( optimization_result, output_filename, 'SOCGA' )
    elif isinstance(optimization_result,MOCGA.Population):
        _export_MOGA( optimization_result, output_filename, 'MOCGA' )
    elif isinstance(optimization_result,downhill_simplex.Simplex):
        _export_downhill_simplex( optimization_result, output_filename, 'downhill_simplex' )        
    else:
        warnings.warn('Please use the SOBGA, SOCGA or MOCGA class.')
        return

def _export_SOGA( optimization_result, output_filename, optimizer ):
    f = optimization_result.config['objective_function']
    x = optimization_result.optimized_parameters
    fw = open( output_filename, 'w' )
    fw.write('%s\n' %('-'*80))
    fw.write('RESULTS\n')
    fw.write('minimum             %+15.12e\n' %(optimization_result.development['best_fitness'][-1]))    
    blanks = ' '*10
    for ix,parameter in enumerate(x):
        if ix>9:
            blanks = ' '*9
        fw.write('parameter%d%s%+15.12e\n' %(ix+1,blanks,parameter.value))
    fw.write('feasibility         %+15.12e\n' %(f.determine_feasibility(x)))
    fw.write('optimizer           %s\n' %(optimizer))
    fw.write('%s\n' %('-'*80))
    fw.write('OBJECTIVE FUNCTION\n')
    fw.write('function name       %s\n' %(f.definition.__name__))
    fw.write('dimension           %d\n' %(f.dimension))    
    if f.global_minimum is not None:
        fw.write('global minimum      %+15.12e\n' %(f.global_minimum))
    else:
        fw.write('global minimum      %s\n' %('unknown'))
    for subject in f.subject_to:
        fw.write('subject to          %s\n' %(subject.__name__))        
    fw.write('%s\n' %('-'*80))
    fw.write('OBJECTIVE FUNCTION PARAMETERS\n')      
    for ix,parameter in enumerate(x):
        fw.write('parameter%d         \n' %(ix+1))
        fw.write('subject to          %s\n' %(inspect.getsourcelines(parameter.subject_to)[0][0]).strip())
        fw.write('mapping             %s\n' %(inspect.getsourcelines(parameter.mapping)[0][0]).strip())
        fw.write('value               %+15.12e\n' %(parameter.value))
        if parameter.true_value is not None:
            fw.write('true value          %+15.12e\n' %(parameter.true_value))
        else:
            fw.write('true value          %s\n' %('unknown'))
        fw.write('mapping parameter   %+15.12e\n' %(parameter.mapping_parameter))
        fw.write('bits                %d\n' %(parameter.bits))    
        fw.write('code                %s\n' %(''.join(str(int(x)) for x in parameter.code)))    
        fw.write('feasible            %r\n' %(parameter.is_feasible()))   
    fw.write('%s\n' %('-'*80))
    fw.write('GENETIC ALGORITHM SETTINGS\n')
    fw.write('pop_size            %d\n' %(optimization_result.config['pop_size']))
    fw.write('max_iterations      %d\n' %(optimization_result.config['max_iterations']))
    fw.write('nsel_method         %s\n' %(optimization_result.config['nsel_method']))
    fw.write('nsel_rate           %4.2f\n' %(optimization_result.config['nsel_rate']))
    fw.write('nsel_thresh         %4.2f\n' %(optimization_result.config['nsel_thresh']))
    fw.write('msel_method         %s\n' %(optimization_result.config['msel_method']))
    fw.write('msel_n_participants %d\n' %(optimization_result.config['msel_n_participants']))
    fw.write('mat_method          %s\n' %(optimization_result.config['mat_method']))
    if 'mat_blend' in optimization_result.config:
        fw.write('mat_blend           %r\n' %(optimization_result.config['mat_blend']))
    fw.write('mut_rate            %4.2f\n' %(optimization_result.config['mut_rate']))
    fw.write('mut_elitism         %r\n' %(optimization_result.config['mut_elitism']))
    fw.write('%s\n' %('-'*80))
    fw.write('POPULATION\n')
    fw.write('function calls      %d\n' %(optimization_result.development['function_calls']))
    fw.write('best fitness        %+15.12e\n' %(optimization_result.development['best_fitness'][-1]))
    fw.write('worst fitness       %+15.12e\n' %(optimization_result.development['worst_fitness'][-1]))
    fw.write('mean fitness        %+15.12e\n' %(optimization_result.development['mean_fitness'][-1]))
    fw.write('std fitness         %+15.12e\n' %(optimization_result.development['std_fitness'][-1]))
    fw.write('fittest indvidual\n')
    fw.write('born in generation  %d\n' %(optimization_result.individuals[0].development['generation'][0]))
    fw.write('worst indvidual\n')
    fw.write('born in generation  %d\n' %(optimization_result.individuals[-1].development['generation'][0]))
    fw.write('all individuals of final population\n')   
    fw.write('%12s' % 'individual')
    fw.write('%22s' % 'f')
    for ix,fun in enumerate(f.parameters):
        fw.write('%22s' % ('x%d' % (ix+1)))
    fw.write('%12s' % 'feasibility')
    fw.write('\n')    
    for ix,ind in enumerate(optimization_result.individuals):
        fw.write('%12d' % ix)
        fw.write('%+22.14e' % ind.development['fitness'][-1])
        x = ind.genes
        for xi in x:
            fw.write('%+22.14e' % xi.value)
        fw.write('%12.2f' % ind.development['feasibility'][-1])    
        fw.write('\n')
    fw.close()
    
def _export_MOGA( optimization_result, output_filename, optimizer ):
    f = optimization_result.config['objective_function']
    fw = open( output_filename, 'w' )    
    fw.write('%s\n' %('-'*80))
    fw.write('RESULTS\n')
    fw.write('Excerpt of Pareto front sorted by crowding distance\n')
    fw.write('%12s' % 'individual')
    for ix,fun in enumerate(f.objective_functions):
        fw.write('%12s' % ('f%d' % (ix+1)))
    for ix,fun in enumerate(f.parameters):
        fw.write('%12s' % ('x%d' % (ix+1)))
    fw.write('%12s' % 'feasibility')
    fw.write('\n')
    for ix in list(range(5)):
        if ix>=len(optimization_result.archive):
            break
        fw.write('%12d' % ix)
        x = optimization_result.archive[ix].genes
        minima = f.evaluate( x )
        for minimum in minima:
            fw.write('%+12.4e' % minimum)
        for xi in x:
            fw.write('%+12.4e' % xi.value)
        fw.write('%12.2f' % optimization_result.archive[ix].feasibility)
        fw.write('\n')
    fw.write('optimizer           %s\n' %(optimizer))
    fw.write('%s\n' %('-'*80))
    fw.write('OBJECTIVE FUNCTION\n')
#    fw.write('function name       %s\n' %(f))
    fw.write('dimension           %d\n' %(f.dimension))   
    fw.write('objectives          %d\n' %(f.objective_dimension))   
    for subject in f.subject_to:
        fw.write('subject to          %s\n' %(subject.__name__))              
    for ix,fi in enumerate(f.objective_functions):
        fw.write('objective%d         \n' %(ix+1))
        fw.write('function name       %s\n' %(fi.definition.__name__))
        fw.write('dimension           %d\n' %(fi.dimension))    
        if fi.global_minimum is not None:
            fw.write('global minimum      %+15.12e\n' %(fi.global_minimum))
        else:
            fw.write('global minimum      %s\n' %('unknown'))
        for subject in fi.subject_to:
            fw.write('subject to          %s\n' %(subject.__name__))                
    fw.write('%s\n' %('-'*80))
    fw.write('OBJECTIVE FUNCTION PARAMETERS\n')      
    for ix,parameter in enumerate(x):
        fw.write('parameter%d         \n' %(ix+1))
        fw.write('subject to          %s\n' %(inspect.getsourcelines(parameter.subject_to)[0][0]).strip())
        fw.write('mapping             %s\n' %(inspect.getsourcelines(parameter.mapping)[0][0]).strip())
        fw.write('bits                %d\n' %(parameter.bits))    
    fw.write('%s\n' %('-'*80))
    fw.write('GENETIC ALGORITHM SETTINGS\n')
    fw.write('pop_size            %d\n' %(optimization_result.config['pop_size']))
    fw.write('max_iterations      %d\n' %(optimization_result.config['max_iterations']))
    fw.write('archive_size        %d\n' %(optimization_result.config['archive_size']))    
    fw.write('nsel_rate           %4.2f\n' %(optimization_result.config['nsel_rate']))    
    fw.write('msel_method         %s\n' %(optimization_result.config['msel_method']))
    fw.write('msel_n_participants %d\n' %(optimization_result.config['msel_n_participants']))
    fw.write('mat_method          %s\n' %(optimization_result.config['mat_method']))
    if 'mat_blend' in optimization_result.config:
        fw.write('mat_blend           %r\n' %(optimization_result.config['mat_blend']))
    fw.write('mut_rate            %4.2f\n' %(optimization_result.config['mut_rate']))
    fw.write('mut_elitism         %r\n' %(optimization_result.config['mut_elitism']))
    fw.write('%s\n' %('-'*80))
    fw.write('ARCHIVE\n')
    fw.write('function calls      %d\n' %(optimization_result.function_calls))    
    fw.write('%12s' % 'individual')
    for ix,fun in enumerate(f.objective_functions):
        fw.write('%12s' % ('f%d' % (ix+1)))
    for ix,fun in enumerate(f.parameters):
        fw.write('%12s' % ('x%d' % (ix+1)))
    fw.write('%12s' % 'feasibility')
    fw.write('%12s' % 'rank')
    fw.write('\n')
    for ix in list(range(len(optimization_result.archive))):
        fw.write('%12d' % ix)
        x = optimization_result.archive[ix].genes
        minima = f.evaluate( x )
        for minimum in minima:
            fw.write('%+12.4e' % minimum)
        for xi in x:
            fw.write('%+12.4e' % xi.value)
        fw.write('%12.2f' % optimization_result.archive[ix].feasibility)
        fw.write('%12d' % optimization_result.archive[ix].rank)
        fw.write('\n')
    fw.close()
    
        
def _export_downhill_simplex( optimization_result, output_filename, optimizer ):
    f = optimization_result.objective_function
    x = optimization_result.optimized_parameters
    fw = open( output_filename, 'w' )
    fw.write('%s\n' %('-'*80))
    fw.write('RESULTS\n')
    fw.write('minimum             %+15.12e\n' %(f.evaluate(x)))    
    blanks = ' '*10
    for ix,parameter in enumerate(x):
        if ix>9:
            blanks = ' '*9
        fw.write('parameter%d%s%+15.12e\n' %(ix+1,blanks,parameter.value))
    fw.write('feasibility         %+15.12e\n' %(f.determine_feasibility(x)))
    fw.write('optimizer           %s\n' %(optimizer))
    fw.write('%s\n' %('-'*80))
    fw.write('OBJECTIVE FUNCTION\n')
    fw.write('function name       %s\n' %(f.definition.__name__))
    fw.write('dimension           %d\n' %(f.dimension))    
    if f.global_minimum is not None:
        fw.write('global minimum      %+15.12e\n' %(f.global_minimum))
    else:
        fw.write('global minimum      %s\n' %('unknown'))
    for subject in f.subject_to:
        fw.write('subject to          %s\n' %(subject.__name__))        
    fw.write('%s\n' %('-'*80))
    fw.write('OBJECTIVE FUNCTION PARAMETERS\n')      
    for ix,parameter in enumerate(x):
        fw.write('parameter%d         \n' %(ix+1))
        fw.write('subject to          %s\n' %(inspect.getsourcelines(parameter.subject_to)[0][0]).strip())
        fw.write('mapping             %s\n' %(inspect.getsourcelines(parameter.mapping)[0][0]).strip())
        fw.write('value               %+15.12e\n' %(parameter.value))
        if parameter.true_value is not None:
            fw.write('true value          %+15.12e\n' %(parameter.true_value))
        else:
            fw.write('true value          %s\n' %('unknown'))
        fw.write('mapping parameter   %+15.12e\n' %(parameter.mapping_parameter))
        fw.write('bits                %d\n' %(parameter.bits))    
        fw.write('code                %s\n' %(''.join(str(int(x)) for x in parameter.code)))    
        fw.write('feasible            %r\n' %(parameter.is_feasible()))   
    fw.write('%s\n' %('-'*80))    
    fw.close()
    
    