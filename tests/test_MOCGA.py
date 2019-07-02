# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
from picasso.utils.optimization import MOCGA
from picasso.utils.optimization import benchmarks

def test_methods():
	f = benchmarks.function1()
	config = {'objective_function':f,
			  'pop_size':10}
	moof = MOCGA.Population(config)
	moof.start_evolution()
	assert moof.archive[0].rank == 0
