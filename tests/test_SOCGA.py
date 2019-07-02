# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
import numpy as np
from picasso.utils.optimization import SOCGA
from picasso.utils.optimization import benchmarks


def test_methods():
	f = benchmarks.six_hump_camel_back()
	config = {'objective_function':f,
			  'pop_size':10,
			  'mut_elitism':True}
	camels = SOCGA.Population(config)
	camels.initialize_individuals()
	camels.test_individuals()
	camels.individuals = camels.sort_individuals()
	camels.update_statistics()
	feasibility = [camel.development['feasibility'][-1] for camel in camels.individuals]
	ix = np.argmax(feasibility)
	assert ix == 0
	fitness = [camel.development['fitness'][-1] for camel in camels.individuals]
	ix = np.argmax(fitness)
	assert ix == 0
	ix = np.argmin(fitness)
	assert ix == 9
	camels.iteration += 1
	camels.dismiss_individuals()
	assert not len(camels.individuals) == 10
	fitness2 = [camel.development['fitness'][-1] for camel in camels.individuals]
	assert fitness[0] == fitness2[0]
	camels.select_parents()
	camels.mate_individuals()
	camels.test_individuals()
	camels.individuals = camels.sort_individuals()
	camels.update_statistics()
	fitness3 = [camel.development['fitness'][-1] for camel in camels.individuals]
	assert len(camels.individuals) == 10
	camels.mutate_individuals()
	camels.test_individuals()
	camels.individuals = camels.sort_individuals()
	camels.update_statistics()
	fitness4 = [camel.development['fitness'][-1] for camel in camels.individuals]
	assert not np.mean(fitness3) == np.mean(fitness4)
	ix = np.argmax(fitness4)
	assert ix == 0
	ix = np.argmin(fitness4)
	assert ix == 9       

def test_six_hump_camel_back():
	f = benchmarks.six_hump_camel_back()
	config = {'objective_function':f}
	minima = []
	for ix in range(100):
		camels = SOCGA.Population(config)
		camels.start_evolution()
		minima.append( f.evaluate( camels.optimized_parameters ) )
	res = [abs(minimum-f.global_minimum) for minimum in minima]
	assert np.amax(res) < 1e-1
	

 