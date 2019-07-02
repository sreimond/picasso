# -*- coding: utf-8 -*-
"""
@author: sreimond
"""

import pytest
from picasso.utils.data_visualization import colors

    
def test_rgb2hsv():
	rgb = [26,8,100]
	hsv_true = [251.74,0.92,0.3922]
	hsv_test = colors.rgb2hsv(rgb)
	assert hsv_test[0] == pytest.approx(hsv_true[0],1e-2)
	assert hsv_test[1] == pytest.approx(hsv_true[1],1e-2)
	assert hsv_test[2] == pytest.approx(hsv_true[2],1e-2)

def test_hsv2rgb():
	hsv = [251.74,0.92,0.3922]
	rgb_true = [26.0,8.0,100.0]
	rgb_test = colors.hsv2rgb(hsv)
	assert rgb_test[0] == pytest.approx(rgb_true[0],1e-1)
	assert rgb_test[1] == pytest.approx(rgb_true[1],1e-1)
	assert rgb_test[2] == pytest.approx(rgb_true[2],1e-1)

def test_rgb2yuv():
	rgb = [26,8,100]
	yuv_true = [23.9,37.5,1.9]
	yuv_test = colors.rgb2yuv(rgb)
	assert yuv_test[0] == pytest.approx(yuv_true[0],1e-1)
	assert yuv_test[1] == pytest.approx(yuv_true[1],1e-1)
	assert yuv_test[2] == pytest.approx(yuv_true[2],1e-1)
	
def test_yuv2rgb():
	yuv = [23.87,37.46,1.87]
	rgb_true = [26.0,8.0,100.0]
	rgb_test = colors.yuv2rgb(yuv)
	assert rgb_test[0] == pytest.approx(rgb_true[0],1e-1)
	assert rgb_test[1] == pytest.approx(rgb_true[1],1e-1)
	assert rgb_test[2] == pytest.approx(rgb_true[2],1e-1)

def test_rgb2hex():
	rgb = [26,8,100]
	hex_true = '#1a0864'
	hex_test = colors.rgb2hex(rgb)
	assert hex_test == hex_true
	
def test_hex2rgb():
	hex_code = '#1a0864'
	rgb_true = [26.0,8.0,100.0]
	rgb_test = colors.hex2rgb(hex_code)
	assert rgb_test[0] == pytest.approx(rgb_true[0],1e-1)
	assert rgb_test[1] == pytest.approx(rgb_true[1],1e-1)
	assert rgb_test[2] == pytest.approx(rgb_true[2],1e-1)
	
def test_hex2hsv():
	hex_code = '#1a0864'
	hsv_true = [251.74,0.92,0.3922]
	hsv_test = colors.hex2hsv(hex_code)
	assert hsv_test[0] == pytest.approx(hsv_true[0],1e-2)
	assert hsv_test[1] == pytest.approx(hsv_true[1],1e-2)
	assert hsv_test[2] == pytest.approx(hsv_true[2],1e-2)

def test_hsv2hex():
	hsv = [251.74,0.92,0.3922]
	hex_true = '#1a0864'
	hex_test = colors.hsv2hex(hsv)
	assert hex_test == hex_true
	
def test_hsv2yuv():
	hsv = [251.74,0.92,0.3922]
	yuv_true = [23.9,37.5,1.9]
	yuv_test = colors.hsv2yuv(hsv)
	assert yuv_test[0] == pytest.approx(yuv_true[0],1e-1)
	assert yuv_test[1] == pytest.approx(yuv_true[1],1e-1)
	assert yuv_test[2] == pytest.approx(yuv_true[2],1e-1)

def test_yuv2hsv():
	yuv = [23.87,37.46,1.87]
	hsv_true = [251.74,0.92,0.3922]
	hsv_test = colors.yuv2hsv(yuv)
	assert hsv_test[0] == pytest.approx(hsv_true[0],1e-2)
	assert hsv_test[1] == pytest.approx(hsv_true[1],1e-2)
	assert hsv_test[2] == pytest.approx(hsv_true[2],1e-2)

def test_hex2yuv():
	hex_code = '#1a0864'
	yuv_true = [23.9,37.5,1.9]
	yuv_test = colors.hex2yuv(hex_code)
	assert yuv_test[0] == pytest.approx(yuv_true[0],1e-1)
	assert yuv_test[1] == pytest.approx(yuv_true[1],1e-1)
	assert yuv_test[2] == pytest.approx(yuv_true[2],1e-1)

def test_yuv2hex():
	yuv = [23.87,37.463,1.869]
	hex_true = '#1a0864'
	hex_test = colors.yuv2hex(yuv)
	assert hex_test == hex_true

def test_name2rgb():
	color_name = 'goldenrod'
	rgb_true = [218.0,165.0,32.0]
	rgb_test = colors.name2rgb(color_name)
	assert rgb_test[0] == pytest.approx(rgb_true[0],1e-1)
	assert rgb_test[1] == pytest.approx(rgb_true[1],1e-1)
	assert rgb_test[2] == pytest.approx(rgb_true[2],1e-1)

def test_name2hsv():
	color_name = 'goldenrod'
	hsv_true = [42.9,0.8532,0.8549]
	hsv_test = colors.name2hsv(color_name)
	assert hsv_test[0] == pytest.approx(hsv_true[0],1e-2)
	assert hsv_test[1] == pytest.approx(hsv_true[1],1e-2)
	assert hsv_test[2] == pytest.approx(hsv_true[2],1e-2)

def test_name2yuv():
	color_name = 'goldenrod'
	yuv_true = [165.7,-65.8,45.9]
	yuv_test = colors.name2yuv(color_name)
	assert yuv_test[0] == pytest.approx(yuv_true[0],1e-1)
	assert yuv_test[1] == pytest.approx(yuv_true[1],1e-1)
	assert yuv_test[2] == pytest.approx(yuv_true[2],1e-1)

def test_rgb2name():
	rgb = [220.0,163.0,33.0]
	name_true = 'goldenrod'
	name_test = colors.rgb2name(rgb)
	assert name_true == name_test

def test_hsv2name():
	hsv = [43,0.83,0.84]
	name_true = 'goldenrod'
	name_test = colors.hsv2name(hsv)
	assert name_true == name_test
	
def test_yuv2name():
	yuv = [164,-65,45]
	name_true = 'goldenrod'
	name_test = colors.yuv2name(yuv)
	assert name_true == name_test
	
def test_color_gradient():
	colors_rgb = ([1,15,61],[87,1,66],[100,100,200])
	gradient_rgb = colors.color_gradient(colors_rgb,
										 color_count=10,
										 color_representation='rgb',
										 interp='rgb')
	colors_hsv = ()
	for rgb in colors_rgb:
		colors_hsv += (colors.rgb2hsv(rgb),)
	gradient_hsv = colors.color_gradient(colors_hsv,
										 color_count=10,
										 color_representation='hsv',
										 interp='hsv')
	colors_yuv = ()
	for rgb in colors_rgb:
		colors_yuv += (colors.rgb2yuv(rgb),)
	gradient_yuv = colors.color_gradient(colors_yuv,
										 color_count=10,
										 color_representation='yuv',
										 interp='yuv')      
	colors_hex = ()
	for rgb in colors_rgb:
		colors_hex += (colors.rgb2hex(rgb),)
	gradient_hex = colors.color_gradient(colors_hex,
										 color_count=10,
										 color_representation='hex',
										 interp='rgb')
	colors_name = ()
	for rgb in colors_rgb:
		colors_name += (colors.rgb2name(rgb),)
	gradient_name = colors.color_gradient(colors_name,
										  color_count=10,
										  color_representation='name',
										  interp='rgb')
	assert len(gradient_rgb) == 10
	assert len(gradient_hsv) == 10
	assert len(gradient_yuv) == 10
	assert len(gradient_hex) == 10
	assert len(gradient_name) == 10
	
def test_color_schemes_tol():
	cmap = colors.color_schemes_tol(3,'qualitative')
	assert cmap[1][1] == 204
	cmap = colors.color_schemes_tol(30,'sequential')
	assert len(cmap) == 30
	cmap = colors.color_schemes_tol(30,'diverging')
	assert len(cmap) == 30
	cmap = colors.color_schemes_tol(30,'rainbow')
	assert len(cmap) == 30
	cmap = colors.color_schemes_tol(3,'hues')
	assert cmap[1][1] == 119
        
