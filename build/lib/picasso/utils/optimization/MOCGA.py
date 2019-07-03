# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 11:56:46 2018

This module is an implementation of the (M)ingle (O)bjective (C)ontinuous 
(G)enetic (A)lgorithm. Using uon-dominated sorting, archive and crowding.

@author: sreimond
"""


import numpy as np
from copy import deepcopy
import itertools

class Population( object ):
    """ 
    The `Population` class holds Individuals and optimizes the multi-objective 
    function as defined in the config dictionary.
    """    
    def __init__( self, config ):
        self.config = _default_config(config)
        self.individuals = []
        self.archive = []
        self.iteration = 0
        self.function_calls = 0
    
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
                
    def test_individuals(self):
        for individual in self.individuals:
            individual.determine_fitness()
            self.function_calls += 1
                                
    def rank_individuals(self,group):
        for individual in group:
            individual.rank = None
        current_rank = -1
        none_count  = np.inf
        while none_count > 0:
            current_rank += 1
            none_count = sum([ind.rank is None for ind in group])
            for individual1 in group:
                if (individual1.rank is not None):
                    continue
                dominated_by = []
                for individual2 in group:
                    if (individual2.rank is not None) and (individual2.rank < current_rank):
                        continue
                    if individual1==individual2 and none_count > 1:
                        continue
                    dominations_count = 0
                    for ix,fitness1 in enumerate(individual1.fitness):
                        if fitness1 < individual2.fitness[ix]:
                            dominations_count += 1
                    dominated_by.append( dominations_count==len(individual1.fitness) )
                on_pareto_front = not np.any( dominated_by )
                if on_pareto_front:
                    individual1.rank = deepcopy( current_rank )
    
    def crowd_individuals(self,group):
        for individual in group:
            individual.crowding_distance = None
        for individual1 in group:
            d = 0
            for ix,fi in enumerate(individual1.fitness):
                neighbor_pos_ix = None
                neighbor_neg_ix = None
                neighbor_pos = np.inf
                neighbor_neg = np.inf
                min_fitness = np.inf
                max_fitness = -np.inf
                for jx,individual2 in enumerate(group):
                    fj = individual2.fitness[ix]
                    delta = fi-fj
                    if (np.sign(delta)==1) and (abs(delta)<neighbor_pos):
                        neighbor_pos = abs(delta)
                        neighbor_pos_ix = jx
                    if (np.sign(delta)==-1) and (abs(delta)<neighbor_neg):
                        neighbor_neg = abs(delta)
                        neighbor_neg_ix = jx
                    if fj < min_fitness:
                        min_fitness = fj
                    if fj > max_fitness:
                        max_fitness = fj
                d += (neighbor_pos+neighbor_neg)/(max_fitness-min_fitness+1.0)
            individual1.crowding_distance = d        
        
    def sort_individuals(self,group):
        rank = [individual.rank for individual in group]
        neg_crowding = [-individual.crowding_distance for individual in group]
        neg_feasibility = [-individual.feasibility for individual in group]
#        ix = np.lexsort((neg_crowding,neg_feasibility,rank))
        ix = np.lexsort((neg_crowding,rank,neg_feasibility))
        return [group[jx] for jx in ix]       
    
    def dismiss_individuals(self):
        pop_size = self.config['pop_size']
        rate = self.config['nsel_rate']
        ix = int(pop_size - pop_size * rate )
        ix_dismiss = list(range(ix,pop_size))
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
                mums_rank = [self.individuals[mum].rank for mum in mums]
                dads_rank = [self.individuals[dad].rank for dad in dads]
                parent1.append(mums[np.argmin(mums_rank)])
                parent2.append(dads[np.argmin(dads_rank)])
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
    
    def archive_individuals(self):
        for current_individual in self.individuals:
            individual = Individual(self)
            individual.dna = current_individual.dna
            individual.update_genes()
            individual.fitness = current_individual.fitness
            individual.feasibility = current_individual.feasibility
            self.archive.append(individual)
            
    def declutter_archive(self):
        while len(self.archive) > self.config['archive_size']:
            del self.archive[-1] 
            
    def start_evolution(self):
        self.initialize_individuals()
        self.test_individuals()
        self.rank_individuals(group=self.individuals)
        self.crowd_individuals(group=self.individuals)
        self.individuals = self.sort_individuals(group=self.individuals)
        for _ in range(self.config['max_iterations']-1):
            self.iteration += 1
            self.dismiss_individuals()
            self.select_parents()
            self.mate_individuals()
            self.mutate_individuals()
            self.test_individuals()
            self.rank_individuals(group=self.individuals)
            self.crowd_individuals(group=self.individuals)
            self.individuals = self.sort_individuals(group=self.individuals)
            self.archive_individuals()
            self.rank_individuals(group=self.archive)
            self.crowd_individuals(group=self.archive)
            self.archive = self.sort_individuals(group=self.archive)    
            self.declutter_archive()       


class Individual( object ):
    """ 
    The `Individual` class defines Individuals and assigns Genes.
    """   
    def __init__(self,population):
        self.population = population
        self.fitness = None
        self.feasibility = None
        self.rank = None
        self.crowding_distance = None
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
        self.fitness = [-c for c in cost]        
        self.feasibility = self.population.config['objective_function'].determine_feasibility( self.genes )
    
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
        config['pop_size'] = 40
    if not ('max_iterations' in config):
        config['max_iterations'] = 40
    if not ('archive_size' in config):
        config['archive_size'] = 500
    if not ('nsel_rate' in config):
        config['nsel_rate'] = 0.7
    if not ('msel_method' in config):
        config['msel_method'] = 'random_picking'
    if not ('msel_n_participants' in config):
        config['msel_n_participants'] = 3
    if not ('mat_method' in config):
        config['mat_method'] = 'single_point'
    if not ('mat_blend' in config):
        config['mat_blend'] = True
    if not ('mut_rate' in config):
        config['mut_rate'] = 0.3
    if not ('mut_elitism' in config):
        config['mut_elitism'] = True
    return config       

