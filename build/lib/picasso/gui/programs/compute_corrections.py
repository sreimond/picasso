# -*- coding: utf-8 -*-
"""
Created on Thu May 17 15:19:12 2018

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
from picasso.utils.files import polygon_files as pof
from picasso.data import data_functions_rgi as dfr
from picasso.gui import project_data_functions
from picasso.gui.programs import groops_interface
from picasso.gui.programs import compute_grace
import logging

# global variables
data_path = None
polygon = None
pip_ix = None
groops_bin = None


def compute_corrections(ui):
    # change button text
    logging.info('Corrections computations started.')
    ui.pb_id002.setText('In progress...')
    ui.pb_id002.setStyleSheet("background-color: rgb(255, 250, 205); color: rgb(255,0,0);")
    # set data path
    global data_path    
    data_path = ui.le_id223.text()
    # groops bin
    global groops_bin
    groops_bin = ui.le_id218.text()
    if ui.cb_id079.isChecked():
        groops_bin = ui.le_id220.text()
    # polygon
    global polygon
    global pip_ix
    if ui.rb_id050.isChecked():
        polygon = compute_grace._region_settings(ui)
        reference_grid = pkg_resources.resource_filename('picasso.src.picasso.data', 'grids/geographical_30-30.grid')
        grid = np.array(np.genfromtxt(reference_grid,skip_header=2),ndmin=2)
        pip_ix = list(range(grid.shape[0]))
        if polygon is not None:
            ix = pip_ix
            p = Pool(None)
            pip_ix = p.map(partial(compute_grace._in_polygon_aux, polygon=polygon, grid=grid), ix, chunksize=25)            
            pip_ix = np.nonzero(pip_ix)[0]
    # dirs
    project_data_functions.create_project_directories(ui)
    project_data_functions.clean_up_correction_directories(ui)
    # hydrology
    _compute_hydrology_grids(ui)
    if ui.cb_id015.isChecked():
        _compute_danubia_grid(ui)
    if ui.rb_id110.isChecked():
        _compute_goco_grids_spline(ui)
    # GIA        
    _compute_gia_grids(ui)
    # create time series files
    _time_series(ui)
    # change back button text
    ui.pb_id002.setText('Compute Corrections')
    ui.pb_id002.setStyleSheet("background-color: rgb(239, 240, 241); color: rgb(0,0,0);")
    logging.info('Corrections computations finished.')

def _compute_hydrology_grids(ui):
    logging.info('Computing hydrological grids...')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    grids_path = os.path.join(grace_path,'grids')
    corrections_path = os.path.join(project_dir,project_name,'corrections')
    # model names
    models = ('gldas','lsdm','wghm')
    if ui.rb_id109.isChecked():
        models += ('goco',)
    # loop over all dates
    mjd_start, mjd_end = compute_grace._date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")  
        grid_file = os.path.join(grids_path,'%s-grid.txt' % (ymstr))
        for name in models:
            hydro_grid = os.path.join(data_path,name.upper(),'%s_TWS_%s.txt' % (name,ymstr))
            if name is 'goco':
                hydro_grid = os.path.join(data_path,'GOCO','grids','GOCO-%s.txt' % (ymstr))
            pointmass_file = os.path.join(corrections_path,name,'grids','%s_%s-pointmass.txt' % (name,ymstr))
            # positions
            if ui.rb_id050.isChecked():
                try:
                    grid = np.array(np.genfromtxt(hydro_grid,skip_header=2),ndmin=2)
                    grid = grid[pip_ix,:]
                    lon,lat,h,area_i,tws_i = grid[:,0],grid[:,1],grid[:,2],grid[:,3],grid[:,4]       
                except:
                    continue
            elif ui.rb_id051.isChecked():                
                try:
                    lon0,lat0,h0,area0,tws0 = np.genfromtxt(hydro_grid,skip_header=2,unpack=True)
                    lon,lat,h,area = np.genfromtxt(grid_file,skip_header=2,unpack=True)
                    tws_i = griddata((lon0,lat0), (tws0), (lon, lat), method=ui.cob_id003.currentText())
                    area_i = griddata((lon0,lat0), (area0), (lon, lat), method=ui.cob_id003.currentText())
                except:
                    continue
            point_count = np.size(lon)
            # area
            if ui.rb_id052.isChecked():
                area = ui.dsb_id001.value()/point_count
            elif ui.rb_id053.isChecked():
                area = (area_i*defaults.EARTH_EQUATORIAL_RADIUS()*defaults.EARTH_EQUATORIAL_RADIUS())*1e-6
            elif ui.rb_id054.isChecked():
                area = ui.dsb_id012.value()/point_count
            # density
            density = 1025.0
            if ui.rb_id056.isChecked():
                density = ui.dsb_id013.value()
            # conversion to pointmass
            kg_i = tws_i * (area*1e6) * density
            # write file
            groops_interface.make_grid_file(groops_bin,pointmass_file,lon,lat,h,area,*(kg_i,np.array(0)))
            logging.info('Current file: %s',pointmass_file)
    logging.info('Computing hydrological grids finished.')
            
def _compute_danubia_grid(ui):
    logging.info('Computing PROMET/DANUBIA grids...')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    grids_path = os.path.join(grace_path,'grids')
    corrections_path = os.path.join(project_dir,project_name,'corrections')
    fh = Dataset(os.path.join(data_path,'DANUBIA','Data-Danubia_update','Danubia.nc'), mode='r')
    fh.set_auto_mask(False)
    lon0 = fh.variables['longitude'][:]
    lat0 = fh.variables['latitude'][:]
    LON0,LAT0 = np.meshgrid(lon0,lat0)
    mjd0 = fh.variables['time'][:]
    # loop over all dates
    mjd_start, mjd_end = compute_grace._date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")  
        grid_file = os.path.join(grids_path,'%s-grid.txt' % (ymstr))
        pointmass_file = os.path.join(corrections_path,'promet','grids','promet_%s-pointmass.txt' % (ymstr))
        ix = np.argwhere(mjd0==mjd)
        if not np.any(ix):
            continue
        ix = int(np.asscalar(ix))
        cr = fh.variables['channel_runoff'][ix,:]
        et = fh.variables['evapotranspiration'][ix,:]
        gs = fh.variables['glacier_snow'][ix,:]
        im = fh.variables['ice_melt'][ix,:]
        pc = fh.variables['precipitation'][ix,:]
        sm = fh.variables['soil_moisture'][ix,:]
        tr = fh.variables['total_runoff'][ix,:]
        ts = fh.variables['total_snow'][ix,:]
        tws0 = eval(ui.le_id225.text())*1e-3
        # mkae grid
        grid = np.zeros((np.size(tws0),5))
        grid[:,0] = LON0.flatten()
        grid[:,1] = LAT0.flatten()
        grid[:,4] = tws0.flatten()
        # positions
        if ui.rb_id050.isChecked():
            try:
                if (polygon is not None) and (pip is not None):
                    grid = grid[pip,:]
                elif (polygon is not None) and (pip is None):    
                    ix = list(range(grid.shape[0]))
                    p = Pool(None)
                    pip = p.map(partial(compute_grace._in_polygon_aux, polygon=polygon, grid=grid), ix, chunksize=25)
                    pip = np.nonzero(pip)[0]
                    grid = grid[pip,:]
                lon,lat,h,area_i,tws_i = grid[:,0],grid[:,1],grid[:,2],grid[:,3],grid[:,4]       
            except:
                continue
        elif ui.rb_id051.isChecked():                
            try:
                lon0,lat0,h0,area0,tws0 = grid[:,0],grid[:,1],grid[:,2],grid[:,3],grid[:,4]
                lon,lat,h,area = np.genfromtxt(grid_file,skip_header=2,unpack=True)
                tws_i = griddata((lon0,lat0), (tws0), (lon, lat), method=ui.cob_id003.currentText())
                area_i = griddata((lon0,lat0), (area0), (lon, lat), method=ui.cob_id003.currentText())
            except:
                continue
        point_count = np.size(lon)
        # area
        if ui.rb_id052.isChecked():
            area = ui.dsb_id001.value()/point_count
        elif ui.rb_id053.isChecked():
            area = 1.0 # 1 km^2 per pixel
        elif ui.rb_id054.isChecked():
            area = ui.dsb_id012.value()/point_count
        # density
        density = 1025.0
        if ui.rb_id056.isChecked():
            density = ui.dsb_id013.value()
        # conversion to pointmass
        kg_i = tws_i * (area*1e6) * density
        # write file
        groops_interface.make_grid_file(groops_bin,pointmass_file,lon,lat,h,area,*(kg_i,np.array(0)))
        logging.info('Current file: %s',pointmass_file)
    fh.close()
    logging.info('Computing PROMET/DANUBIA grids finished.')
    
def _compute_gia_grids(ui):
    logging.info('Computing GIA grids...')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    grids_path = os.path.join(grace_path,'grids')
    corrections_path = os.path.join(project_dir,project_name,'corrections')
    # loop over all dates
    mjd_start, mjd_end = compute_grace._date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")  
        grid_file = os.path.join(grids_path,'%s-grid.txt' % (ymstr))
        for name in ('a','ice5','ice6','ice6_grace','klemann'):
            try:
                gia_grid = glob.glob(os.path.join(data_path,'GIA','grids',name.upper(),'*-%s.txt' % (ymstr)))[0]
            except IndexError:
                continue
            pointmass_file = os.path.join(corrections_path,name,'grids','%s_%s-pointmass.txt' % (name,ymstr))
            # positions
            if ui.rb_id050.isChecked():
                try:
                    grid = np.array(np.genfromtxt(gia_grid,skip_header=2),ndmin=2)
                    grid = grid[pip_ix,:]
                    lon,lat,h,area_i,tws_i = grid[:,0],grid[:,1],grid[:,2],grid[:,3],grid[:,4]       
                except:
                    continue
            elif ui.rb_id051.isChecked():             
                try:
                    lon0,lat0,h0,area0,tws0 = np.genfromtxt(gia_grid,skip_header=2,unpack=True)
                    lon,lat,h,area = np.genfromtxt(grid_file,skip_header=2,unpack=True)
                    tws_i = griddata((lon0,lat0), (tws0), (lon, lat), method=ui.cob_id003.currentText())
                    area_i = griddata((lon0,lat0), (area0), (lon, lat), method=ui.cob_id003.currentText())
                except:
                    continue
            point_count = np.size(lon)
            # area
            if ui.rb_id052.isChecked():
                area = ui.dsb_id001.value()/point_count
            elif ui.rb_id053.isChecked():
                area = (area_i*defaults.EARTH_EQUATORIAL_RADIUS()*defaults.EARTH_EQUATORIAL_RADIUS())*1e-6
            elif ui.rb_id054.isChecked():
                area = ui.dsb_id012.value()/point_count
            # density
            density = 1025.0
            if ui.rb_id056.isChecked():
                density = ui.dsb_id013.value()
            # conversion to pointmass
            kg_i = tws_i * (area*1e6) * density
            # write file
            groops_interface.make_grid_file(groops_bin,pointmass_file,lon,lat,h,area,*(kg_i,np.array(0)))
            logging.info('Current file: %s',pointmass_file)
    logging.info('Computing GIA grids finished.')

def _compute_goco_grids_spline(ui):
    logging.info('Computing GOCO grids (splines)...')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    grids_path = os.path.join(grace_path,'grids')
    corrections_path = os.path.join(project_dir,project_name,'corrections','goco')
    # loop over all dates
    mjd_start, mjd_end = compute_grace._date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")  
        grid_file = os.path.join(grids_path,'%s-grid.txt' % (ymstr))
        goco_grid = os.path.join(data_path,'GOCO','grids','GOCO-%s.txt' % (ymstr))
        pointmass_file = os.path.join(corrections_path,'grids','goco_%s-pointmass.txt' % (ymstr))
        # positions              
        try:                    
            temp_dir = tempfile.mkdtemp()
            tmp_input_grid_file = os.path.join(temp_dir, 'grid_in.txt')
            tmp_output_grid_file = os.path.join(temp_dir, 'grid_out.txt')
            lon,lat,h,area = np.genfromtxt(grid_file,skip_header=2,unpack=True)
            groops_interface.make_grid_file(groops_bin,tmp_input_grid_file,lon,lat,h,area)
            groops_interface.compute_goco_grid(groops_bin,tmp_input_grid_file,tmp_output_grid_file,mjd1)
            lon,lat,h,area_i,tws_i = np.genfromtxt(tmp_output_grid_file,skip_header=2,unpack=True)
            shutil.rmtree(temp_dir)
            lon0,lat0,h0,area0,tws0 = np.genfromtxt(goco_grid,skip_header=2,unpack=True)
            area_i = griddata((lon0,lat0), (area0), (lon, lat), method=ui.cob_id003.currentText())
        except:
            continue                        
        point_count = np.size(lon)
        # area
        if ui.rb_id052.isChecked():
            area = ui.dsb_id001.value()/point_count
        elif ui.rb_id053.isChecked():
            area = (area_i*defaults.EARTH_EQUATORIAL_RADIUS()*defaults.EARTH_EQUATORIAL_RADIUS())*1e-6
        elif ui.rb_id054.isChecked():
            area = ui.dsb_id012.value()/point_count
        # density
        density = 1025.0
        if ui.rb_id056.isChecked():
            density = ui.dsb_id013.value()
        # conversion to pointmass
        kg_i = tws_i * (area*1e6) * density
        # write file
        groops_interface.make_grid_file(groops_bin,pointmass_file,lon,lat,h,area,*(kg_i,np.array(0)))
        logging.info('Current file: %s',pointmass_file)
    logging.info('Computing GOCO grids (splines) done.')
            
def _time_series(ui):
    logging.info('Creating time series files from point mass grids...')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    corrections_path = os.path.join(project_dir,project_name,'corrections')
    # loop over correction types
    for correction in ('a','gldas','ice5','ice6','ice6_grace','klemann','lsdm','promet','wghm','goco'):
        # time series files
        time_series_path = os.path.join(corrections_path,correction,'time_series')
        time_series_file = os.path.join(time_series_path,'%s.txt' % correction)
        with open(time_series_file,'w') as f:
            f.write('%25s %25s %25s\n' % ('time (mjd)','mass (kg)','sigma (kg)'))
        # loop over all dates
        mjd_start, mjd_end = compute_grace._date_settings(ui)
        mjd = mjd_start
        while mjd<=mjd_end:
            mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
            mjd = (mjd2+1)
            ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")  
            pointmass_file = os.path.join(corrections_path,correction,'grids','%s_%s-pointmass.txt' % (correction,ymstr))
            try:
                lon,lat,h,area,x,sigma_x = np.genfromtxt(pointmass_file,skip_header=2,unpack=True)
                kg = np.sum(x)
                sigma_kg = np.sqrt(np.sum(sigma_x**2.0)/sigma_x.size)
                with open(time_series_file,'a+') as f:
                    f.write('%+25.16e %+25.16e %+25.16e\n' % (mjd1,kg,sigma_kg))
            except:
                continue
            logging.info('Creating time series file: %s_%s',correction,ymstr)
    logging.info('Time series files created.')
            