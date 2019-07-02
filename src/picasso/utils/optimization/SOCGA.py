# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 11:05:59 2018

This module is an implementation of the (S)ingle (O)bjective (C)ontinuous 
(G)enetic (A)lgorithm.

@author: sreimond
"""

import numpy as np
from copy import deepcopy
import itertools

class Population( object ):
    """ 
    The `Population` class holds Individuals and optimizes the objective 
    function as defined in the config dictionary.
    """    
    def __init__( self, config ):
        self.config = _default_config(config)
        self.individuals = []
        self.development = {'generation':[],
                            'function_calls':0,
                            'best_fitness':[],
                            'worst_fitness':[],
                            'mean_fitness':[],
                            'std_fitness':[]}
        self.optimized_parameters = []
        self.iteration = 0
    
    def initialize_individuals(self):
        for _ in range(self.config['pop_size']):
            self.add_individual()
    
    def add_individual(self):
        d_min = -np.inf
        for ix in range(10):
            individual = Individual(self)
            if ix==0:
                self.individuals.append(individual)
                continue
            d = []
            for individual2 in self.individuals:
                if individual == individual2:
                    continue
                d.append( individual.determine_distance(individual2) )            
            if np.min(d) > d_min:
                d_min = np.min(d)
                del self.individuals[-1]
                self.individuals.append(individual)        
        
    def sort_individuals(self):
        # fittest indivdual has max. fitness and is located at index 0
        neg_fitness = [-individual.development['fitness'][-1] for individual in self.individuals]
        neg_feasibility = [-individual.development['feasibility'][-1] for individual in self.individuals]        
        ix = np.lexsort((neg_fitness,neg_feasibility))
#        ix = np.argsort(neg_fitness)
        return [self.individuals[jx] for jx in ix]
         
    def update_statistics(self):        
        fitness = [individual.development['fitness'][-1] for individual in self.individuals]
        self.development['generation'].append( self.iteration )
        self.development['best_fitness'].append( np.max(fitness) )
        self.development['worst_fitness'].append( np.min(fitness) )
        self.development['mean_fitness'].append( np.mean(fitness) )
        self.development['std_fitness'].append( np.std(fitness) )
    
    def test_individuals(self):
        for individual in self.individuals:
            individual.determine_fitness()
            self.development['function_calls'] += 1
    
    def dismiss_individuals(self):
        method = self.config['nsel_method']
        pop_size = self.config['pop_size']
        rate = self.config['nsel_rate']
        thresh = self.config['nsel_thresh']
        if method == 'sorted_list':
            ix = int(pop_size - pop_size * rate )
            ix_dismiss = list(range(ix,pop_size))
        elif method == 'thresholding':
            ix_dismiss = []
            for ix,individual in self.individuals:
                if (-individual.development['fitness'][-1]) > thresh:
                    ix_dismiss.append(ix)
        ix_dismiss.sort(reverse=True)
        for ix in ix_dismiss:
            del self.individuals[ix]
        while len(self.individuals)<2:
            self.add_individual()
        if (pop_size - len(self.individuals)) % 2 != 0:
            del self.individuals[-1]
    
    def select_parents(self):
        pop_size = int(self.config['pop_size'])
        individual_count = len(self.individuals)
        missing_count = pop_size-individual_count
        parents_count = int(missing_count/2)
        participants_count = int(self.config['msel_n_participants'])
        method = self.config['msel_method']
        do_replace = False
        if ((individual_count < missing_count) or
            (individual_count < participants_count * 2)):
                do_replace = True
        if method == 'sorted_list':
            ix = range(missing_count)
            self.parents = (ix[::2],ix[1::2])
        elif method == 'random_picking':
            ix = np.random.choice(range(individual_count),
                                  replace=do_replace,
                                  size=missing_count)
            self.parents = (ix[0:parents_count],ix[parents_count:missing_count])
        elif method == 'tournament':
            parent1 = []
            parent2 = []
            for _ in range(parents_count):
                ix = np.random.choice(range(individual_count),
                                      replace=do_replace,
                                      size=participants_count*2)
                mums = ix[0:participants_count]
                dads = ix[participants_count:participants_count*2]
                mums_fitness = [self.individuals[mum].development['fitness'][-1] for mum in mums]
                dads_fitness = [self.individuals[dad].development['fitness'][-1] for dad in dads]
                parent1.append(mums[np.argmax(mums_fitness)])
                parent2.append(dads[np.argmax(dads_fitness)])
            self.parents = (parent1,parent2)
            
    def mate_individuals(self):
        ix_parents1, ix_parents2 = self.parents
        for ix,ix_parent1 in enumerate(ix_parents1):
            ix_parent2 = ix_parents2[ix]
            parent1 = self.individuals[ix_parent1]
            parent2 = self.individuals[ix_parent2]
            children = parent1.mate(parent2)
            self.individuals.append(children[0])
            self.individuals.append(children[1])
    
    def mutate_individuals(self):
        ix = 0
        if self.config['mut_elitism']:
            ix = 1
        for jx,individual in enumerate(self.individuals):
            if jx < ix:
                continue
            individual.mutate()
    
    def start_evolution(self):
        self.initialize_individuals()
        self.test_individuals()
        self.individuals = self.sort_individuals()
        self.update_statistics()
        for _ in range(self.config['max_iterations']-1):
            self.iteration += 1
            self.dismiss_individuals()
            self.select_parents()
            self.mate_individuals()
            self.mutate_individuals()
            self.test_individuals()
            self.individuals = self.sort_individuals()
            self.update_statistics()
        self.optimized_parameters = self.individuals[0].genes

class Individual( object ):
    """ 
    The `Individual` class defines Individuals and assigns Genes.
    """   
    def __init__(self,population):
        self.population = population
        self.development = {'generation':[],'fitness':[],'feasibility':[]}
        self.initialize()
    
    def initialize(self):
        self.genes = []        
        for parameter in self.population.config['objective_function'].parameters:
            gene = deepcopy( parameter )
            gene.set_random_value()
            self.genes.append( gene )
        self.update_dna()
    
    def determine_distance(self,individual):
        d = 0
        for ix,gene in enumerate(self.genes):
            d += (gene.value - individual.genes[ix].value)**2.0
        return np.sqrt(d)
    
    def update_dna(self):
        self.dna = np.array([gene.mapping_parameter for gene in self.genes])
    
    def update_genes(self):
        for ix,gene in enumerate(self.genes):
            value = self.dna[ix]
            gene.set_mapping_parameter(value)
            
    def determine_fitness(self):
        cost = self.population.config['objective_function'].evaluate( self.genes )
        fitness = -cost
        feasibility = self.population.config['objective_function'].determine_feasibility( self.genes )
        self.development['generation'].append(self.population.iteration)
        self.development['fitness'].append(fitness)
        self.development['feasibility'].append(feasibility)
    
    def mutate(self):
        for gene in self.genes:
            if np.random.random() > (1.0-self.population.config['mut_rate']):
                gene.set_random_value()
        self.update_dna()        
            
    def mate(self,individual):
        total_bits = len(self.dna)        
        parents_dna = (self.dna,individual.dna)
        children_dna = parents_dna
        mask = np.zeros(total_bits,dtype=bool)
        beta = 0
        if self.population.config['mat_method'] == 'single_point':
            ix = np.random.randint(0,total_bits)
            mask[0:ix] = True
        elif self.population.config['mat_method'] == 'two_point':
            ix = np.sort(np.random.choice(total_bits, 2, replace=False))
            mask[ix[0]:ix[1]] = True
        elif self.population.config['mat_method'] == 'uniform':
            mask[:] = np.random.randint(0,2,total_bits)            
        if self.population.config['mat_blend']:
            beta = np.random.random()   
        children_dna[0][mask] *= beta + (1.0-beta) * parents_dna[1][mask]
        children_dna[1][mask] *= beta + (1.0-beta) * parents_dna[0][mask]
        children = (Individual(self.population),Individual(self.population))
        for ix,child in enumerate(children):
            child.dna = children_dna[ix]
            child.update_genes()
        return children                
            
def _default_config(config={}):
    if not ('pop_size' in config):
        config['pop_size'] = 20
    if not ('max_iterations' in config):
        config['max_iterations'] = 30
    if not ('nsel_method' in config):
        config['nsel_method'] = 'sorted_list'
    if not ('nsel_rate' in config):
        config['nsel_rate'] = 0.7
    if not ('nsel_thresh' in config):
        config['nsel_thresh'] = 0
    if not ('msel_method' in config):
        config['msel_method'] = 'tournament'
    if not ('msel_n_participants' in config):
        config['msel_n_participants'] = 3
    if not ('mat_method' in config):
        config['mat_method'] = 'uniform'
    if not ('mat_blend' in config):
        config['mat_blend'] = True
    if not ('mut_rate' in config):
        config['mut_rate'] = 0.3
    if not ('mut_elitism' in config):
        config['mut_elitism'] = True
    return config       














