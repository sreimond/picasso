# -*- coding: utf-8 -*-
"""
Created on Wed May 23 13:33:08 2018

@author: sreimond
"""

import os
import uuid
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
from functools import partial
import glob
import pkg_resources
import shutil, tempfile
from multiprocessing import Pool
from PyQt5.QtWidgets import QMessageBox, QInputDialog, QLineEdit
from picasso.utils.constants import defaults
from picasso.utils.dates_and_time import date_functions
from picasso.utils.geometry import map_projections
from picasso.utils.files import polygon_files as pof
from picasso.utils.geometry import polygons, points
from picasso.data import data_functions_rgi as dfr
from picasso.gui import project_data_functions
from picasso.gui.programs import groops_interface
from picasso.gui.programs import compute_grace
import logging

# global variables
data_path = None
grid_file = None
buffer_grid_file = None
min_degree = None
max_degree = None
gauss = None
groops_bin = None

def compute_internal_validation(ui):
    # change button text
    logging.info('Internal validation computations started.')
    ui.pb_id003.setText('In progress...')
    ui.pb_id003.setStyleSheet("background-color: rgb(255, 250, 205); color: rgb(255,0,0);")
    # set data path
    global data_path    
    data_path = ui.le_id223.text()
    # global vars
    global min_degree
    global max_degree
    global gauss
    min_degree = ui.sb_id018.value()
    max_degree = ui.sb_id019.value()
    gauss = ui.dsb_id014.value()
    # groops bin
    global groops_bin
    groops_bin = ui.le_id218.text()
    if ui.cb_id079.isChecked():
        groops_bin = ui.le_id220.text()
    # polygon
    global grid_file
    global buffer_grid_file
    temp_dir = tempfile.mkdtemp()
    grid_file = os.path.join(temp_dir, 'tmp.grid')
    buffer_grid_file = os.path.join(temp_dir, 'tmp_buffer.grid')
    if ui.rb_id057.isChecked():
        polygon = compute_grace._region_settings(ui)
        lon,lat,area = dfr.rgi_grid_in_polygon(polygon,ignore_zeros=ui.cb_id075.isChecked())
        grid_file_global = pkg_resources.resource_filename('picasso.data', 'grids/geographical_30-30.grid')
        lon0,lat0,h0,area0 = np.genfromtxt(grid_file_global,skip_header=2,unpack=True)
        area_i = griddata((lon0,lat0), (area0), (lon, lat), method='nearest')
        groops_interface.make_grid_file(groops_bin,grid_file,lon,lat,0,area_i)
    if ui.rb_id100.isChecked():
        grid_file = ui.le_id236.text()
    # dirs
    project_data_functions.create_project_directories(ui)
    project_data_functions.clean_up_internal_validation_directories(ui)
#    # M1, M2
    _compute_M1_grids(ui)
    _compute_M2_grids(ui)
    _compute_M1_syn_grids(ui)
    _time_series(ui)
    shutil.rmtree(temp_dir)
    # change back button text
    ui.pb_id003.setText('Compute Internal Validation')
    ui.pb_id003.setStyleSheet("background-color: rgb(239, 240, 241); color: rgb(0,0,0);")
    logging.info('Internal validation computations finished.')
    
def _compute_M1_grids(ui):
    logging.info('Computation of M1 grids started.')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    grids_path = os.path.join(grace_path,'grids')
    internal_validation_path = os.path.join(project_dir,project_name,'validation','internal')
    # constant density
    density = 1025.0
    # loop over all dates
    mjd_start, mjd_end = compute_grace._date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")  
        grid_file_pm = os.path.join(grids_path,'%s-grid.txt' % (ymstr))
        pointmass_file = os.path.join(internal_validation_path,'grids','m1_%s-pointmass.txt' % (ymstr))
        temp_dir = tempfile.mkdtemp()
        tmp_grid = os.path.join(temp_dir, 'tmp.grid')
        try:
            if ui.rb_id058.isChecked():
                _compute_gfc2grid(ui,groops_bin,grid_file_pm,tmp_grid,mjd1,min_degree,max_degree,gauss)
            elif (ui.rb_id057.isChecked() or ui.rb_id100.isChecked()):
                _compute_gfc2grid(ui,groops_bin,grid_file,tmp_grid,mjd1,min_degree,max_degree,gauss)
        except:
            continue
        if not os.path.isfile(tmp_grid):
            continue
        lon,lat,h,area_i,tws = np.genfromtxt(tmp_grid,skip_header=2,unpack=True)
        area = (area_i*defaults.EARTH_EQUATORIAL_RADIUS()*defaults.EARTH_EQUATORIAL_RADIUS())*1e-6
        kg = tws * (area*1e6) * density
        groops_interface.make_grid_file(groops_bin,pointmass_file,lon,lat,h,area,*(kg,np.array(0)))
        shutil.rmtree(temp_dir)
        logging.info('Computing grid: %s',pointmass_file)
    logging.info('Computation of M1 grids finished.')

def _compute_M2_grids(ui):
    logging.info('Computation of M2 grids started.')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    grids_path = os.path.join(grace_path,'grids')
    internal_validation_path = os.path.join(project_dir,project_name,'validation','internal')
    # constant density
    density = 1025.0
    # which grid
    if (ui.rb_id057.isChecked() or ui.rb_id100.isChecked()):
        lon,lat,_,_ = np.genfromtxt(grid_file,skip_header=2,unpack=True)
        _compute_buffer_grid(ui,lon,lat)
    # loop over all dates
    mjd_start, mjd_end = compute_grace._date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")  
        grid_file_pm = os.path.join(grids_path,'%s-grid.txt' % (ymstr))
        pointmass_file = os.path.join(internal_validation_path,'grids','m2_%s-pointmass.txt' % (ymstr))
        temp_dir = tempfile.mkdtemp()
        tmp_grid = os.path.join(temp_dir, 'tmp.grid')
        _compute_gfc2grid(ui,groops_bin,buffer_grid_file,tmp_grid,mjd1,min_degree,max_degree,gauss)
        try:
            if ui.rb_id058.isChecked():
                lon,lat,_,_ = np.genfromtxt(grid_file_pm,skip_header=2,unpack=True)
                _compute_buffer_grid(ui,lon,lat)
            _compute_gfc2grid(ui,groops_bin,buffer_grid_file,tmp_grid,mjd1,min_degree,max_degree,gauss)
        except:
            continue
        if not os.path.isfile(tmp_grid):
            continue
        lon,lat,h,area_i,tws = np.genfromtxt(tmp_grid,skip_header=2,unpack=True)
        area = (area_i*defaults.EARTH_EQUATORIAL_RADIUS()*defaults.EARTH_EQUATORIAL_RADIUS())*1e-6
        kg = tws * (area*1e6) * density
        groops_interface.make_grid_file(groops_bin,pointmass_file,lon,lat,h,area,*(kg,np.array(0)))
        shutil.rmtree(temp_dir)
        logging.info('Computing grid: %s',pointmass_file)
    logging.info('Computation of M2 grids finished.')

def _compute_M1_syn_grids(ui):
    logging.info('Computation of M1_syn grids started.')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    internal_validation_path = os.path.join(project_dir,project_name,'validation','internal')
    # global grid
    grid_file_global = pkg_resources.resource_filename('picasso.data', 'grids/geographical_30-30.grid')
    grid0 = np.array(np.genfromtxt(grid_file_global,skip_header=2),ndmin=2)
    # indices
    if (ui.rb_id057.isChecked() or ui.rb_id100.isChecked()):
        lon,lat,_,_ = np.genfromtxt(grid_file,skip_header=2,unpack=True)
        ix = []
        for jx in list(range(np.size(lon))):            
            ix.append(_return_index_closest_point(np.array([lon[jx],lat[jx]]),grid0[:,0:2]))
    # constant density
    density = 1025.0
    # loop over all dates
    mjd_start, mjd_end = compute_grace._date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")  
        m1_file = os.path.join(internal_validation_path,'grids','m1_%s-pointmass.txt' % (ymstr))
        m2_file = os.path.join(internal_validation_path,'grids','m2_%s-pointmass.txt' % (ymstr))
        m1_star_file = os.path.join(internal_validation_path,'grids','m1_star_%s-pointmass.txt' % (ymstr))
        m2_syn_file = os.path.join(internal_validation_path,'grids','m2_syn_%s-pointmass.txt' % (ymstr))
        m1_syn_file = os.path.join(internal_validation_path,'grids','m1_syn_%s-pointmass.txt' % (ymstr))
        try:
            lon1,lat1,h1,area1,x1,sigma_x1 = np.genfromtxt(m1_file,skip_header=2,unpack=True)
            lon2,lat2,h2,area2,x2,sigma_x2 = np.genfromtxt(m2_file,skip_header=2,unpack=True)
        except:
            continue
        kg1 = np.sum(x1)
        kg2 = np.sum(x2)
        IAF = kg2/kg1
        m1_star = kg1*IAF
        groops_interface.make_grid_file(groops_bin,m1_star_file,lon1,lat1,h1,area1,*(x1*IAF,np.array(0)))
        if (ui.rb_id057.isChecked() or ui.rb_id100.isChecked()):
            m1_star_i = np.zeros(grid0.shape[0])
            m1_star_i[ix] = x1*IAF/((area1*1e6) * density)
        else:
            m1_star_i = griddata((lon1,lat1), x1*IAF/((area1*1e6) * density), (grid0[:,0], grid0[:,1]), method='cubic')
        temp_dir = tempfile.mkdtemp()
        tmp_grid1 = os.path.join(temp_dir, 'tmp1.grid')
        tmp_grid2 = os.path.join(temp_dir, 'tmp2.grid')
        groops_interface.make_grid_file(groops_bin,tmp_grid1,grid0[:,0],grid0[:,1],grid0[:,2],grid0[:,3],*(m1_star_i,np.array(0)))
        groops_interface.compute_grid_to_gfc_to_grid(groops_bin,tmp_grid1,tmp_grid2,buffer_grid_file,gauss=gauss)
        if not os.path.isfile(tmp_grid2):
            continue
        lon,lat,h,area_i,tws = np.genfromtxt(tmp_grid2,skip_header=2,unpack=True)
        area = (area_i*defaults.EARTH_EQUATORIAL_RADIUS()*defaults.EARTH_EQUATORIAL_RADIUS())*1e-6
        kg = tws * (area*1e6) * density
        kg2_syn = np.sum(kg)
        groops_interface.make_grid_file(groops_bin,m2_syn_file,lon,lat,h,area,*(kg,np.array(0)))
        signal_loss = 100.0 * (1.0 - kg2_syn/kg2)
        kg2_syn_star = kg2 * (1.0 + (1.0 - kg2_syn/kg2))
        FAF = kg2_syn_star/kg1
        kg1_syn = FAF * kg1
        groops_interface.make_grid_file(groops_bin,m1_syn_file,lon1,lat1,h1,area1,*(x1*FAF,np.array(0)))
        shutil.rmtree(temp_dir)
        logging.info('Computing grid: %s',m1_syn_file)
    logging.info('Computation of M1_syn grids finished.')

def _time_series(ui):
    logging.info('Creating time series files from point mass grids...')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    internal_validation_path = os.path.join(project_dir,project_name,'validation','internal')
    # loop over element types
    for element in ('m1','m2','m1_star','m2_syn','m1_syn'):
        # time series files
        time_series_path = os.path.join(internal_validation_path,'time_series')
        time_series_file = os.path.join(time_series_path,'%s.txt' % element)
        with open(time_series_file,'w') as f:
            f.write('%25s %25s %25s\n' % ('time (mjd)','mass (kg)','sigma (kg)'))
        # loop over all dates
        mjd_start, mjd_end = compute_grace._date_settings(ui)
        mjd = mjd_start
        while mjd<=mjd_end:
            mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
            mjd = (mjd2+1)
            ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")  
            pointmass_file = os.path.join(internal_validation_path,'grids','%s_%s-pointmass.txt' % (element,ymstr))
            try:
                lon,lat,h,area,x,sigma_x = np.genfromtxt(pointmass_file,skip_header=2,unpack=True)
                kg = np.sum(x)
                sigma_kg = np.sqrt(np.sum(sigma_x**2.0)/sigma_x.size)
                with open(time_series_file,'a+') as f:
                    f.write('%+25.16e %+25.16e %+25.16e\n' % (mjd1,kg,sigma_kg))
            except:
                continue
            logging.info('Creating time series file: %s_%s',element,ymstr)
    mjd,m1,_ = np.genfromtxt(os.path.join(time_series_path,'m1.txt'),skip_header=1,unpack=True)
    _,m2,_ = np.genfromtxt(os.path.join(time_series_path,'m2.txt'),skip_header=1,unpack=True)
    _,m2_syn,_ = np.genfromtxt(os.path.join(time_series_path,'m2_syn.txt'),skip_header=1,unpack=True)
    data = {}
    data['iaf'] = m2/m1
    data['signal_loss'] = 100.0 * (1.0 - m2_syn/m2)
    data['m2_syn_star'] = m2 * (1.0 + (1.0 - m2_syn/m2))
    data['faf'] = data['m2_syn_star']/m1
    for key in ('iaf','signal_loss','m2_syn_star','faf'):
        time_series_file = os.path.join(time_series_path,'%s.txt' % key)
        with open(time_series_file,'w') as f:
            f.write('%25s %25s %25s\n' % ('time (mjd)','value','dummy'))
        for jx in list(range(np.size(mjd))):
            with open(time_series_file,'a+') as f:
                f.write('%+25.16e %+25.16e %+25.16e\n' % (mjd[jx],data[key][jx],0))        
    logging.info('Time series files created.')
                
def _compute_gfc2grid(ui,groops_bin,input_grid,output_grid,mjd,min_degree,max_degree,gauss):
    if ui.rb_id059.isChecked():
        groops_interface.compute_gfc_grid(groops_bin,input_grid,output_grid,mjd,min_degree,max_degree,gauss,data_center='itsg')
    elif ui.rb_id060.isChecked():
        groops_interface.compute_gfc_grid(groops_bin,input_grid,output_grid,mjd,min_degree,max_degree,gauss,data_center='gfz')
    elif ui.rb_id061.isChecked():
        groops_interface.compute_gfc_grid(groops_bin,input_grid,output_grid,mjd,min_degree,max_degree,gauss,data_center='csr')
    elif ui.rb_id062.isChecked():
        groops_interface.compute_gfc_grid(groops_bin,input_grid,output_grid,mjd,min_degree,max_degree,gauss,data_center='jpl')
        
def _compute_buffer_grid(ui,lon0,lat0):
    buffer_size_km = ui.dsb_id015.value()
    grid_file_global = pkg_resources.resource_filename('picasso.data', 'grids/geographical_30-30.grid')
    lon,lat,h,area = np.genfromtxt(grid_file_global,skip_header=2,unpack=True)
    central_point = (0,0)
    u,v = map_projections.azimuthal_equidistant(lon,lat,central_point=central_point)
    ix = np.zeros(len(lon),dtype=bool)
    for jx in list(range(len(lon0))):
        point = points.Point2D(lon0[jx],lat0[jx])
        ui,vi = map_projections.azimuthal_equidistant(point.x,point.y,central_point=central_point)
        d_test = np.sqrt((u-ui)**2.0+(v-vi)**2.0)
        dx = d_test<=(buffer_size_km*1e3)
        ix[dx] = True    
    groops_interface.make_grid_file(groops_bin,buffer_grid_file,lon[ix],lat[ix],h[ix],area[ix])
    
def _return_index_closest_point(point, points):    
    deltas = points - point
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)
    
    
    
    
    
