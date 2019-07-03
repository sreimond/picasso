# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 11:33:00 2018

@author: sreimond
"""

import numpy as np
import csv
from multiprocessing import Pool
from functools import partial
import pkg_resources
from picasso.utils.files import polygon_files as pof
from picasso.utils.geometry import points

def wgms_region_code_to_rgi_region_id(wgms_code):
    wgms_code = wgms_code.lower()
    regions = {}
    regions['acn'] = 3
    regions['acs'] = 4
    regions['ala'] = 1
    regions['ant'] = 19
    regions['asc'] = 13
    regions['ase'] = 15
    regions['asn'] = 10
    regions['asw'] = 14
    regions['cau'] = 12
    regions['ceu'] = 11
    regions['grl'] = 5
    regions['isl'] = 6
    regions['nzl'] = 18
    regions['rua'] = 9
    regions['san'] = 17
    regions['sca'] = 8
    regions['sjm'] = 7
    regions['trp'] = 16
    regions['wna'] = 2
    return regions[wgms_code]  

def rgi_region_id_to_area(rgi_id):
    # RGI Consortium, 2017
    # [km2]
    areas = dict([(1,86725.053),
                  (2,14524.224),
                  (3,105110.642),
                  (4,40888.228),
                  (5,89717.066),
                  (6,11059.7),
                  (7,33958.934),
                  (8,2949.103),
                  (9,51591.6),
                  (10,2410.051),
                  (11,2092.146),
                  (12,1306.992),
                  (13,49303.415),
                  (14,33568.298),
                  (15,14734.204),
                  (16,2341.036),
                  (17,29429.08),
                  (18,1161.801),
                  (19,132867.22)])
    return areas[rgi_id]

def rgi_region_id_to_rgi_data_name(rgi_id):    
    regions = dict([(1,'01_rgi60_Alaska'),
                    (2,'02_rgi60_WesternCanadaUS'),
                    (3,'03_rgi60_ArcticCanadaNorth'),
                    (4,'04_rgi60_ArcticCanadaSouth'),
                    (5,'05_rgi60_GreenlandPeriphery'),
                    (6,'06_rgi60_Iceland'),
                    (7,'07_rgi60_Svalbard'),
                    (8,'08_rgi60_Scandinavia'),
                    (9,'09_rgi60_RussianArctic'),
                    (10,'10_rgi60_NorthAsia'),
                    (11,'11_rgi60_CentralEurope'),
                    (12,'12_rgi60_CaucasusMiddleEast'),
                    (13,'13_rgi60_CentralAsia'),
                    (14,'14_rgi60_SouthAsiaWest'),
                    (15,'15_rgi60_SouthAsiaEast'),
                    (16,'16_rgi60_LowLatitudes'),
                    (17,'17_rgi60_SouthernAndes'),
                    (18,'18_rgi60_NewZealand'),
                    (19,'19_rgi60_AntarcticSubantarctic')])
    return regions[rgi_id]

def rgi_region_id_to_polygon(rgi_id):
    data_name = rgi_region_id_to_rgi_data_name(rgi_id)    
    polygon_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'polygons/%s.txt' % data_name)
    polygon = pof.read_polygon_file( polygon_file )
    return polygon

def rgi_region_id_to_glacier_attributes(rgi_id):
    data_name = rgi_region_id_to_rgi_data_name(rgi_id)
    attributes_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'RGI/%s.csv' % data_name)
    with open(attributes_file) as File:  
        reader = csv.reader(File,delimiter=',')
        next(reader)
        lon = []
        lat = []
        area = []
        for row in reader:
            lon.append(row[4])
            lat.append(row[5])
            area.append(row[8])
    return np.array(lon), np.array(lat), np.array(area)

def rgi_region_id_to_grid(rgi_id,ignore_zeros=True):
    data_name = rgi_region_id_to_rgi_data_name(rgi_id)
    grid_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'RGI/%s.grid' % data_name)
    grid = np.array(np.genfromtxt(grid_file,skip_header=2),ndmin=2)
    if ignore_zeros:
        ix = np.where(grid[:,3]>1e-6)
        grid = grid[ix[0],:]
    return grid[:,0], grid[:,1], grid[:,3]

def glacier_attributes_to_grid():
    lon = ()
    lat = ()
    area = ()
    for ix in list(range(1,20)):
        loni, lati, areai = rgi_region_id_to_glacier_attributes(ix)
        lon += (loni,)
        lat += (lati,)
        area += (areai,)
    lon = np.concatenate(lon).ravel()
    lat = np.concatenate(lat).ravel()
    area = np.concatenate(area).ravel()
    grid = np.zeros((np.size(lon),4),dtype=np.float_)
    grid[:,0] = lon
    grid[:,1] = lat
    grid[:,3] = area
    return grid
    
def glacier_attributes_in_polygon(polygon):
    grid = glacier_attributes_to_grid()
    if polygon is not None:
        ix = list(range(grid.shape[0]))
        p = Pool(None)
        pip = p.map(partial(_in_polygon_aux, polygon=polygon, grid=grid), ix, chunksize=25)
        pip = np.nonzero(pip)[0]
        grid = grid[pip,:]
    return grid[:,0], grid[:,1], grid[:,3]    

def rgi_grid_in_polygon(polygon,ignore_zeros=True):
    grid_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'RGI/00_rgi60_30-30grid_ascii.txt')
    grid = np.array(np.genfromtxt(grid_file,skip_header=4),ndmin=2)
    if ignore_zeros:
        ix = np.where(grid[:,3]>1e-6)
        grid = grid[ix[0],:]
    if polygon is not None:
        ix = list(range(grid.shape[0]))
        p = Pool(None)
        pip = p.map(partial(_in_polygon_aux, polygon=polygon, grid=grid), ix, chunksize=25)
        pip = np.nonzero(pip)[0]        
        grid = grid[pip,:]
    return grid[:,0], grid[:,1], grid[:,3]    
        
def _in_polygon_aux(ix,polygon,grid):
    point = points.Point2D(grid[ix,0],grid[ix,1])
    return polygon.determine_point_location(point)




