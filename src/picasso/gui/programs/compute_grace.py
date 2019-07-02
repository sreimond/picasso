# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 10:38:45 2018

@author: sreimond
"""

import os
import uuid
import numpy as np
import glob
from copy import deepcopy
from multiprocessing import Pool
import subprocess
import time
import pkg_resources
from functools import partial
import shutil, tempfile
from PyQt5.QtWidgets import QMessageBox, QInputDialog, QLineEdit
from picasso.utils.dates_and_time import date_functions
from picasso.utils.files import groops_files as gof
from picasso.utils.files import polygon_files as pof
from picasso.utils.files import optimization_files as fof
from picasso.utils.optimization import SOCGA
from picasso.utils.optimization import MOCGA
from picasso.utils.optimization import downhill_simplex
from picasso.utils.optimization import objective_functions as of
from picasso.utils.geometry import points
from picasso.data import data_functions_rgi as dfr
from picasso.gui import project_data_functions
from picasso.gui.programs import groops_interface
from picasso.gui.programs import regularization
from picasso.gui.programs import optimization
import logging

# global variables
x2kg = 1.0/(4.0*np.pi) # conversion factor (groops density parameter --> kg)
pm_count_is_parameter = False
pm_positions_are_parameters = False
pm_grid_type_is_parameter = False
pm_grid_resolution_is_parameter = False
pm_depths_are_parameters = False
pm_magnitudes_are_parameters = False
regularization_is_parameter = False
current_normals_path = ''
current_mjd = None
current_polygon = None
current_polygon_file = None
current_pm_count = None
current_pm_lon = None
current_pm_lat = None
current_pm_area = None
current_pm_depths = None
current_pm_magnitudes = None
current_regularization = None
current_groops_bin = None
current_groops_bin_cluster = None
current_frontend_cores = None
current_pm_grid_type = None
current_pm_geo_ll = None
current_pm_geo_ul = None
current_pm_trivert_ll = None
current_pm_trivert_ul = None
current_pm_tricent_ll = None
current_pm_tricent_ul = None
current_pm_gauss_ll = None
current_pm_gauss_ul = None
current_pm_reuter_ll = None
current_pm_reuter_ul = None
current_pm_corput_ll = None
current_pm_corput_ul = None
current_pm_driscoll_ll = None
current_pm_driscoll_ul = None
picasso_on_leo = False
ga_solutions = ()
data_path = None

def compute_grace(ui):
    # change button text
    logging.info('GRACE computations started.')
    ui.pb_id001.setText('In progress...')
    ui.pb_id001.setStyleSheet("background-color: rgb(255, 250, 205); color: rgb(255,0,0);")
    # dirs
    project_data_functions.create_project_directories(ui)
    # set data path
    global data_path
    data_path = os.path.join(ui.le_id223.text(),'IFG','raw')
    # get objective function
    objective_function = _determine_objective_function(ui)    
    # export polygon to temporary file
    global current_polygon_file
    temp_dir = tempfile.mkdtemp()
    tmp_txt = os.path.join(temp_dir,'tmppoly.txt')
    current_polygon_file = os.path.join(temp_dir,'tmppoly.xml')
    if current_polygon is not None:
        pof.export_polygon( current_polygon, tmp_txt )
        pof.convert_polygon_file_txt2xml( tmp_txt, current_polygon_file )
    # start computations
    if objective_function is None:
        _no_optimization(ui)
    else:
        _optimization(ui,objective_function)
    shutil.rmtree( temp_dir )
    # change back button text
    ui.pb_id001.setText('Compute GRACE')
    ui.pb_id001.setStyleSheet("background-color: rgb(239, 240, 241); color: rgb(0,0,0);")
    logging.info('GRACE computations finished.')

def _show_dialog_grace_normals():
    msg = QMessageBox()
    msg.setText("Should GRACE normals be reused (if possible)?")
    msg.setWindowTitle("GRACE normals")
    msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
    ret = msg.exec_()
    reuse = None
    if ret == QMessageBox.Yes:
        reuse = True
        logging.info('Only existing normals files will be used.')
    elif ret == QMessageBox.No:
        reuse = False
        logging.info('Normals will be computed.')
    return reuse
            
def _determine_objective_function(ui):
    logging.info('Gathering information concerning the objective function...')
    # global variables
    global pm_count_is_parameter
    global pm_positions_are_parameters
    global pm_grid_type_is_parameter
    global pm_grid_resolution_is_parameter
    global pm_depths_are_parameters
    global pm_magnitudes_are_parameters
    global regularization_is_parameter
    global current_polygon
    global current_pm_count
    global current_pm_lon
    global current_pm_lat
    global current_pm_area
    global current_pm_depths
    global current_pm_magnitudes
    global current_regularization
    global current_pm_grid_type
    global current_pm_geo_ll
    global current_pm_geo_ul
    global current_pm_trivert_ll
    global current_pm_trivert_ul
    global current_pm_tricent_ll
    global current_pm_tricent_ul
    global current_pm_gauss_ll
    global current_pm_gauss_ul
    global current_pm_reuter_ll
    global current_pm_reuter_ul
    global current_pm_corput_ll
    global current_pm_corput_ul
    global current_pm_driscoll_ll
    global current_pm_driscoll_ul
    global current_groops_bin 
    global current_groops_bin_cluster
    global current_frontend_cores
    global picasso_on_leo
    # groops bin
    current_groops_bin = ui.le_id218.text()
    if ui.cb_id079.isChecked():
        picasso_on_leo = True
        current_groops_bin = ui.le_id220.text()
    current_groops_bin_cluster = ui.le_id222.text()
    current_frontend_cores = ui.sb_id016.value()
    # date range configuration
    mjd_start, mjd_end = _date_settings(ui)
    # region of interest
    current_polygon = _region_settings(ui)
    # number of point masses
    current_pm_count = _pm_count_settings(ui)
    pm_count_is_parameter = np.all(np.array(current_pm_count)==-100)
    # positions of point masses
    current_pm_lon, current_pm_lat, current_pm_area = _pm_position_settings(ui,current_polygon)
    pm_positions_are_parameters = np.all(np.array(current_pm_lon)==-100) and np.all(np.array(current_pm_lat)==-100) and np.all(np.array(current_pm_area)==-100)
    pm_grid_type_is_parameter = np.all(np.array(current_pm_lon)==-300) and np.all(np.array(current_pm_lat)==-300) and np.all(np.array(current_pm_area)==-300)
    if pm_grid_type_is_parameter or (np.all(np.array(current_pm_lon)==-200) and np.all(np.array(current_pm_lat)==-200) and np.all(np.array(current_pm_area)==-200)):
        pm_grid_resolution_is_parameter = True
    if pm_grid_resolution_is_parameter:
        pm_count_is_parameter = False
        logging.info('PM count, magnitudes and simultaenous grid refinement optimization not supported. PM counts and magnitudes settings ignored. Single layer depths will be used.')
        ui.rb_id032.setChecked(True)
        ui.rb_id018.setChecked(True)
    # depths of point masses
    current_pm_depths = _pm_depth_settings(ui)
    pm_depths_are_parameters = np.all(np.array(current_pm_depths)==-100)
    # magnitudes of point masses
    current_pm_magnitudes = _pm_magnitude_settings(ui)
    pm_magnitudes_are_parameters = np.all(np.array(current_pm_magnitudes)==-100)
    # regularization
    current_regularization = _regularization_settings(ui)
    regularization_is_parameter = np.all(np.array(current_regularization)==-100)
    # if no optimization takes place, return here
    if not (pm_count_is_parameter or pm_positions_are_parameters or pm_depths_are_parameters or pm_magnitudes_are_parameters or regularization_is_parameter or pm_grid_type_is_parameter or pm_grid_resolution_is_parameter):
        logging.info('No optimization will be applied.')
        return None
    # objective function definition
    objective_function = _objective_function_settings(ui)
    # initial positions, possible point mass counts
    ga_lon, ga_lat, ga_area = _ga_initial_position_settings(ui,current_polygon)
    ga_lon, ga_lat, ga_area = np.array(ga_lon,dtype=float), np.array(ga_lat,dtype=float), np.array(ga_area,dtype=float)
    ga_positions_are_random = np.all(np.array(ga_lon)==-100) and np.all(np.array(ga_lat)==-100) and np.all(np.array(ga_area)==-100)
    if ga_positions_are_random and pm_positions_are_parameters:
        if pm_count_is_parameter:
            ga_count = ui.sb_id008.value()
        else:
            ga_count = ui.sb_id005.value()
    else:
        ga_count = np.size(ga_lon)
        if pm_count_is_parameter and (ui.sb_id008.value()>ga_count):
            ui.sb_id008.setValue(ga_count)
        if ui.rb_id010.isChecked() and (ui.sb_id005.value()<=ga_count):
            ga_count = ui.sb_id005.value()
        elif ui.rb_id010.isChecked() and (ui.sb_id005.value()>ga_count):
            ui.sb_id005.setValue(ga_count) 
    current_pm_count = ga_count    
    # add objective function parameters: pm count
    if pm_count_is_parameter:
        ll, ul = ui.sb_id007.value(), ui.sb_id008.value()
        objective_function.add_parameter(_define_objective_parameter(ui,ll,ul,'pm_count'))
        if isinstance(objective_function,of.MultiObjectiveFunction):
            for f in objective_function.objective_functions:
                f.add_parameter(_define_objective_parameter(ui,ll,ul,'pm_count'))
    # grid types and resolutions
    if pm_grid_resolution_is_parameter:
        ll, ul = 0, 1.0
        current_pm_grid_type = str(ui.cob_id039.currentText())
        current_pm_geo_ll, current_pm_geo_ul = ui.dsb_id037.value(), ui.dsb_id038.value()
        current_pm_trivert_ll, current_pm_trivert_ul = ui.sb_id021.value(), ui.sb_id022.value()
        current_pm_tricent_ll, current_pm_tricent_ul = ui.sb_id023.value(), ui.sb_id024.value()
        current_pm_gauss_ll, current_pm_gauss_ul = ui.sb_id025.value(), ui.sb_id026.value()
        current_pm_reuter_ll, current_pm_reuter_ul = ui.sb_id027.value(), ui.sb_id028.value()
        current_pm_corput_ll, current_pm_corput_ul = ui.sb_id029.value(), ui.sb_id030.value()
        current_pm_driscoll_ll, current_pm_driscoll_ul = ui.sb_id031.value(), ui.sb_id032.value()
        objective_function.add_parameter(_define_objective_parameter(ui,ll,ul,'pm_grid_resolution'))
        if isinstance(objective_function,of.MultiObjectiveFunction):
            for f in objective_function.objective_functions:
                f.add_parameter(_define_objective_parameter(ui,ll,ul,'pm_grid_resolution'))            
    if pm_grid_type_is_parameter:
        ll = 0
        ul = 7.0 
        objective_function.add_parameter(_define_objective_parameter(ui,ll,ul,'pm_grid_type'))
        if isinstance(objective_function,of.MultiObjectiveFunction):
            for f in objective_function.objective_functions:
                f.add_parameter(_define_objective_parameter(ui,ll,ul,'pm_grid_type')) 
    # adjust pm position search if random initial positions are used
    if ga_positions_are_random:
        ui.rb_id030.setChecked(True)
    # add objective function parameters: pm positions
    if pm_positions_are_parameters:
        if ui.rb_id030.isChecked():
            bbox = current_polygon.determine_bounding_box()
            x_min = np.ones(ga_count) * bbox.vertices[0].x
            x_max = np.ones(ga_count) * bbox.vertices[1].x
            y_min = np.ones(ga_count) * bbox.vertices[0].y
            y_max = np.ones(ga_count) * bbox.vertices[2].y
        else:
            delta_max_lon = ui.dsb_id002.value()
            delta_max_lat = ui.dsb_id002.value()
            x_min = ga_lon-delta_max_lon
            x_max = ga_lon+delta_max_lon
            y_min = ga_lat-delta_max_lat
            y_max = ga_lat+delta_max_lat
        for ix in list(range(ga_count)):
            x_mini = x_min[ix]
            x_maxi = x_max[ix]
            y_mini = y_min[ix]
            y_maxi = y_max[ix]
            objective_function.add_parameter(_define_objective_parameter(ui,x_mini,x_maxi,'pm_lon'))
            objective_function.add_parameter(_define_objective_parameter(ui,y_mini,y_maxi,'pm_lat'))
            if isinstance(objective_function,of.MultiObjectiveFunction):
                for f in objective_function.objective_functions:
                    f.add_parameter(_define_objective_parameter(ui,x_mini,x_maxi,'pm_lon'))
                    f.add_parameter(_define_objective_parameter(ui,y_mini,y_maxi,'pm_lat'))
            if ui.rb_id030.isChecked():
                objective_function.add_subject( optimization._subject_pointmass_positions )
                if isinstance(objective_function,of.MultiObjectiveFunction):
                    for f in objective_function.objective_functions:
                        f.add_subject( optimization._subject_pointmass_positions )
    # add objective function parameters: pm depths
    if pm_depths_are_parameters:
        if ui.rb_id033.isChecked():
            depths_count = ga_count
        else:
            depths_count = 1
        for _ in list(range(depths_count)):
            ll, ul = ui.dsb_id003.value(), ui.dsb_id004.value()
            objective_function.add_parameter(_define_objective_parameter(ui,ll,ul,'pm_depth'))
            if isinstance(objective_function,of.MultiObjectiveFunction):
                for f in objective_function.objective_functions:
                    f.add_parameter(_define_objective_parameter(ui,ll,ul,'pm_depth'))
    # add objective function parameters: pm magnitudes
    if pm_magnitudes_are_parameters:
        for _ in list(range(ga_count)):
            ll, ul = ui.dsb_id005.value(), ui.dsb_id006.value()
            ll *= 1e12 * 4.0 * np.pi
            ul *= 1e12 * 4.0 * np.pi
            objective_function.add_parameter(_define_objective_parameter(ui,ll,ul,'pm_magnitude'))
            if isinstance(objective_function,of.MultiObjectiveFunction):
                for f in objective_function.objective_functions:
                    f.add_parameter(_define_objective_parameter(ui,ll,ul,'pm_magnitude'))        
    # add objective function parameters: regularization
    if regularization_is_parameter:
        ll, ul = ui.dsb_id007.value(), ui.dsb_id008.value()
        objective_function.add_parameter(_define_objective_parameter(ui,ll,ul,'regularization'))
        if isinstance(objective_function,of.MultiObjectiveFunction):
            for f in objective_function.objective_functions:
                f.add_parameter(_define_objective_parameter(ui,ll,ul,'regularization'))
    # return instance
    logging.info('Objective function initialized.')
    return objective_function
            

def _date_settings(ui):
    date_begin = ui.de_id001.date().toString("yyyy-MM-dd")
    date_end = ui.de_id002.date().toString("yyyy-MM-dd")
    if ui.rb_id003.isChecked():
        date_begin = "2002-04-01"
        date_end = "2017-12-01"
    mjd_start = date_functions.ymstring2mjd(date_begin)
    mjd_end = date_functions.ymstring2mjd(date_end)
    logging.info('Date range (MJD): %f - %f',mjd_start,mjd_end)
    return mjd_start, mjd_end

def _region_settings(ui):
    if ui.rb_id005.isChecked():
        wgms_id = str(ui.cob_id002.currentText())
        rgi_id = dfr.wgms_region_code_to_rgi_region_id(wgms_id)
        polygon = dfr.rgi_region_id_to_polygon(rgi_id)
        logging.info('Polygon instance based on RGI definitions created.')
    elif ui.rb_id006.isChecked():
        polygon_file = ui.le_id004.text()
        if polygon_file.endswith('.xml'):
            pof.convert_polygon_file_xml2txt( polygon_file )
            polygon_file = polygon_file.replace('.xml','.txt')
        polygon = pof.read_polygon_file( polygon_file )
        logging.info('Polygon instance created based on file: %s',polygon_file)
    elif ui.rb_id007.isChecked():
        polygon = None
        logging.info('No border (polygon) defined.')
    return polygon
    
def _pm_count_settings(ui):
    # -100 = GA
    # -200 = based on positions
    if ui.rb_id008.isChecked():
        pm_count = -100
        logging.info('Point mass count will be optimized using GA.')
    elif ui.rb_id009.isChecked():
        pm_count = -200
        logging.info('Point mass count will determined based on positions.')
    elif ui.rb_id010.isChecked():
        pm_count = ui.sb_id005.value()
        logging.info('Point mass count will be fixed to: %d',pm_count)
    return pm_count        
    
def _pm_position_settings(ui,polygon):
    # -100 = GA
    if ui.rb_id011.isChecked():
        if ui.rb_id112.isChecked(): 
            pm_lon, pm_lat, pm_area = -200, -200, -200
            logging.info('Point mass grid resolution will be optimized using GA.')
        elif ui.rb_id114.isChecked(): 
            pm_lon, pm_lat, pm_area = -300, -300, -300
            logging.info('Point mass grid type and resolution will be optimized using GA.')
        else:
            pm_lon, pm_lat, pm_area = -100, -100, -100
            logging.info('Point mass positions will be optimized using GA.')
    elif ui.rb_id012.isChecked():
        ignore_zeros = ui.cb_id072.isChecked()
        if ui.rb_id005.isChecked():
            wgms_id = str(ui.cob_id002.currentText())
            rgi_id = dfr.wgms_region_code_to_rgi_region_id(wgms_id)
            pm_lon, pm_lat, pm_area = dfr.rgi_region_id_to_grid(rgi_id,ignore_zeros=ignore_zeros)                    
        else:
            pm_lon, pm_lat, pm_area = dfr.rgi_grid_in_polygon(polygon,ignore_zeros=ignore_zeros)
        log_str = 'Point mass positions are fixed based on the RGI grid definition inside the polygon.'
        if ignore_zeros:
            log_str = log_str.replace('polygon.','polygon (zero-area points excluded).')
        logging.info(log_str)
        ui.rb_id027.setChecked(True)
        ui.cb_id073.setChecked(ignore_zeros)
    elif ui.rb_id013.isChecked():
        if ui.rb_id005.isChecked():
            wgms_id = str(ui.cob_id002.currentText())
            rgi_id = dfr.wgms_region_code_to_rgi_region_id(wgms_id)
            pm_lon, pm_lat, pm_area = dfr.rgi_region_id_to_glacier_attributes(rgi_id)
        else:
            pm_lon, pm_lat, pm_area = dfr.glacier_attributes_in_polygon(polygon)
        logging.info('Point mass positions are fixed based on the RGI glacier coordinates inside the polygon.')
        ui.rb_id028.setChecked(True)
    elif ui.rb_id014.isChecked():
        grid_file = ui.le_id005.text()
        grid = np.array(np.genfromtxt(grid_file,skip_header=2),ndmin=2)
        if polygon is not None:
            ix = list(range(grid.shape[0]))
            p = Pool(None)
            pip = p.map(partial(_in_polygon_aux, polygon=polygon, grid=grid), ix, chunksize=25)
            pip = np.nonzero(pip)[0]
            grid = grid[pip,:]
        pm_lon, pm_lat, pm_area = grid[:,0], grid[:,1], grid[:,3]
        logging.info('Point mass positions are fixed inside the polygon based on the grid file: %s',grid_file)
        ui.rb_id029.setChecked(True)
        ui.le_id006.setText(grid_file)
    return pm_lon, pm_lat, pm_area
        
def _in_polygon_aux(ix,polygon,grid):
    point = points.Point2D(grid[ix,0],grid[ix,1])
    return polygon.determine_point_location(point)

def _pm_depth_settings(ui):
    # -100 = GA
    if ui.rb_id015.isChecked():
        pm_depth = -100
        logging.info('Point mass depths will be optimized using GA.')
    else:
        pm_depth = ui.sb_id006.value()
        logging.info('Point mass depths will be fixed to: %f',pm_depth)
    return pm_depth

def _pm_magnitude_settings(ui):
    # -100 = GA
    # -200 = LSA
    if ui.rb_id017.isChecked():
        pm_magnitude = -100
        logging.info('Point mass magnitudes will be optimized using GA.')
    elif ui.rb_id018.isChecked():
        pm_magnitude = -200
        logging.info('Point mass magnitudes will be determined via LSA.')
    return pm_magnitude

def _regularization_settings(ui):
    # -100 = GA
    # -200 = VCE
    # -300 = L-Curve
    # -400 = GCV
    # -500 = No regularization
    if ui.rb_id019.isChecked():
        regularization = -100
        logging.info('Regularization will be optimized using GA.')
    elif ui.rb_id020.isChecked():
        regularization = -200
        logging.info('Regularization will be done via VCE.')
    elif ui.rb_id021.isChecked():
        regularization = -300
        logging.info('Regularization will be done via L-Curve.')
    elif ui.rb_id022.isChecked():
        regularization = -400
        logging.info('Regularization will be done via GCV.')
    elif ui.rb_id023.isChecked():
        regularization = -500
        logging.info('No regularization will be applied.')
    return regularization

def _objective_function_settings(ui):
    # number of objectives
    objectives_count = 0
    for ix in list(range(5,12)):
        key = 'cb_id%03d' % ix
        item = getattr(ui,key)
        objectives_count += int(item.isChecked())
    objectives_count += int(ui.cb_id084.isChecked())
    # determinte objective functions
    if ui.rb_id024.isChecked() and objectives_count>=2:
        ui.rb_id035.setChecked(True) # no binary implementation
        ui.rb_id049.setChecked(True) # no hybrid implementation
        f = _multi_objective(ui)
    else:
        f = _single_objective(ui)
    return f    
            
def _single_objective(ui):
    # accumulate objectives
    logging.info('Single objective optimization will be used.')
    ofs = ()
    if ui.cb_id005.isChecked():
        ofs += ( optimization._minimize_alpha, )
        logging.info('Objective (alpha) added.')
    if ui.cb_id006.isChecked():
        ofs += ( optimization._minimize_alphaxTx, )
        logging.info('Objective (alpha xTx) added.')
    if ui.cb_id007.isChecked():
        ofs += ( optimization._minimize_condN, )
        logging.info('Objective (cond(N)) added.')
    if ui.cb_id008.isChecked():
        ofs += ( optimization._minimize_ePe, )
        logging.info('Objective (ePe) added.')
    if ui.cb_id009.isChecked():
        ofs += ( optimization._minimize_ePealphaxTx, )
        logging.info('Objective (ePe alpha xTx) added.')
    if ui.cb_id010.isChecked():
        ofs += ( optimization._minimize_sigmaxTsigmax, )
        logging.info('Objective (sigmaxTsigmax) added.')
    if ui.cb_id011.isChecked():
        ofs += ( optimization._minimize_xTx, )
        logging.info('Objective (xTx) added.')
    if ui.cb_id084.isChecked():
        ofs += ( optimization._minimize_negatives2nratio, )
        logging.info('Objective (negative signal-to-noise ratio) added.')
    sum_of_ofs = lambda parameters: sum( _of(parameters) for _of in ofs )
    # return an instance of the objective_function class
    f = of.ObjectiveFunction( sum_of_ofs )
    return f

def _multi_objective(ui):
    # add objectives to MulitObjective object
    logging.info('Multi objective optimization (Pareto) will be used.')
    f = of.MultiObjectiveFunction()
    if ui.cb_id005.isChecked():
        f.add_objective_function( of.ObjectiveFunction( optimization._minimize_alpha ) )
        logging.info('Objective (alpha) added.')
    if ui.cb_id006.isChecked():
        f.add_objective_function( of.ObjectiveFunction( optimization._minimize_alphaxTx ) )
        logging.info('Objective (alpha xTx) added.')
    if ui.cb_id007.isChecked():
        f.add_objective_function( of.ObjectiveFunction( optimization._minimize_condN ) )     
        logging.info('Objective (cond(N)) added.')   
    if ui.cb_id008.isChecked():
        f.add_objective_function( of.ObjectiveFunction( optimization._minimize_ePe ) )  
        logging.info('Objective (ePe) added.')              
    if ui.cb_id009.isChecked():
        f.add_objective_function( of.ObjectiveFunction( optimization._minimize_ePealphaxTx ) )
        logging.info('Objective (ePe alpha xTx) added.')
    if ui.cb_id010.isChecked():
        f.add_objective_function( of.ObjectiveFunction( optimization._minimize_sigmaxTsigmax ) )   
        logging.info('Objective (sigmaxTsigmax) added.')     
    if ui.cb_id011.isChecked():
        f.add_objective_function( of.ObjectiveFunction( optimization._minimize_xTx ) )
        logging.info('Objective (xTx) added.')
    if ui.cb_id084.isChecked():
        f.add_objective_function( of.ObjectiveFunction( optimization._minimize_negatives2nratio ) )
        logging.info('Objective (negative signal-to-noise ratio) added.')
    return f

def _ga_initial_position_settings(ui,polygon):
    # -100 = random
    if ui.rb_id026.isChecked():
        ga_lon, ga_lat, ga_area = -100, -100, -100
        logging.info('Initial point mass positions (GA) will be distributed randomly.')
    elif ui.rb_id027.isChecked():
        ignore_zeros = ui.cb_id073.isChecked()
        if ui.rb_id005.isChecked():
            wgms_id = str(ui.cob_id002.currentText())
            rgi_id = dfr.wgms_region_code_to_rgi_region_id(wgms_id)
            ga_lon, ga_lat, ga_area = dfr.rgi_region_id_to_grid(rgi_id,ignore_zeros=ignore_zeros)        
        else:
            ga_lon, ga_lat, ga_area = dfr.rgi_grid_in_polygon(polygon,ignore_zeros=ignore_zeros)
        log_str = 'Initial point mass positions (GA) are fixed based on the RGI grid definition inside the polygon.'
        if ignore_zeros:
            log_str = log_str.replace('polygon.','polygon (zero-area points excluded).')
        logging.info(log_str)
    elif ui.rb_id028.isChecked():
        if ui.rb_id005.isChecked():
            wgms_id = str(ui.cob_id002.currentText())
            rgi_id = dfr.wgms_region_code_to_rgi_region_id(wgms_id)
            ga_lon, ga_lat, ga_area = dfr.rgi_region_id_to_glacier_attributes(rgi_id)
        else:
            ga_lon, ga_lat, ga_area = dfr.glacier_attributes_in_polygon(polygon)
        logging.info('Initial point mass positions (GA) are fixed based on the RGI glacier coordinates inside the polygon.')
    elif ui.rb_id029.isChecked():
        grid_file = ui.le_id006.text()
        grid = np.array(np.genfromtxt(grid_file,skip_header=2),ndmin=2)
        if polygon is not None:
            ix = list(range(grid.shape[0]))
            p = Pool(None)
            pip = p.map(partial(_in_polygon_aux, polygon=polygon, grid=grid), ix, chunksize=25)
            pip = np.nonzero(pip)[0]
            grid = grid[pip,:]
        ga_lon, ga_lat, ga_area = grid[:,0], grid[:,1], grid[:,3]
        logging.info('Initial point mass positions (GA) are fixed inside the polygon based on the grid file: %s',grid_file)
    return ga_lon, ga_lat, ga_area
    
def _define_objective_parameter(ui,ll,ul,label):
    st = lambda value:(value>=ll and value<=ul)
    m  = lambda value: (ul-ll)*0.5*np.sin(value*2.0*np.pi)+(ul+ll)*0.5
    p = of.ObjectiveFunctionParameter(subject_to=st,mapping=m)
    if ui.rb_id035.isChecked():
        p.bits = 0
    p.label = label
    logging.info('Objective parameter added (label, lower bound, upper bound): %s, %f, %f',label,ll,ul)
    return p

def _no_optimization(ui):
    reuse_normals = _show_dialog_grace_normals()
    if reuse_normals is False:
        _make_grids(ui)
        _build_pointmass_normals(ui)
    elif reuse_normals is None:
        return
    _regularization(ui)
    _pointmass_grid(ui)
    _time_series(ui)

def _make_grids(ui):
    logging.info('Creating monthly single solution grids...')
    # output path
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    # groops bin
    groops_bin = ui.le_id218.text()
    if ui.cb_id079.isChecked():
        groops_bin = ui.le_id220.text()
    # grid settings
    polygon = _region_settings(ui)
    pm_lon, pm_lat, pm_area = _pm_position_settings(ui,polygon)
    pm_depths = _pm_depth_settings(ui)
    pm_h = (-pm_depths*1e3)
    # loop over all dates
    mjd_start, mjd_end = _date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")
        grid_file = os.path.join(grace_path,'grids','%s-grid.txt' % ymstr)
        groops_interface.make_grid_file(groops_bin,grid_file,pm_lon,pm_lat,pm_h,pm_area)   
        logging.info('%s created!',grid_file)
    logging.info('Monthly single solution grids created.')

def _build_pointmass_normals(ui,goce_only=False):
    if ui.rb_id002.isChecked():
        _build_pointmass_normals_leo(ui)
        return
    elif ui.rb_id001.isChecked() and ui.cb_id079.isChecked():
        ui.rb_id107.setChecked(True)
        _build_pointmass_normals_leo(ui)
        return
    logging.info('Building monthly single solution normals...')
    # compute goce?
    compute_goce = ui.cb_id074.isChecked()
    if (not compute_goce) and goce_only:
        logging.info('Neither GRACE nor GOCE is selected, nothing is done.')
        return
    # output path
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    normals_path = os.path.join(grace_path,'normals')
    # groops bin
    groops_bin = ui.le_id218.text()
    if ui.cb_id079.isChecked():
        groops_bin = ui.le_id220.text()
    # loop over all dates
    mjd_start, mjd_end = _date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")
        grid_file = os.path.join(grace_path,'grids','%s-grid.txt' % ymstr)  
        if not os.path.exists(grid_file):
            continue
        groops_interface.build_pointmass_normals(groops_bin,mjd1,mjd2,grid_file,normals_path,leo=False,compute_goce=compute_goce,goce_only=goce_only)
        groops_interface.combine_eliminate_solve_pointmass_normals(groops_bin,mjd1,mjd2,grid_file,normals_path,leo=False,compute_goce=compute_goce,goce_only=goce_only)
        logging.info('%s normals built!',ymstr)
    logging.info('Monthly single solution normals built.')

def _build_pointmass_normals_leo(ui,goce_only=False):
    logging.info('Building monthly single solution normals on LEO...')
    # compute goce?
    compute_goce = ui.cb_id074.isChecked()
    if (not compute_goce) and goce_only:
        logging.info('Neither GRACE nor GOCE is selected, nothing is done.')
        return
    # on nodes?
    on_nodes = False
    if ui.rb_id108.isChecked():
        on_nodes = True
        logging.info('Computations will be distributed on LEO nodes.')
    # output path
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    project_dir_leo = str(ui.le_id219.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    normals_path = os.path.join(grace_path,'normals')
    # get temporary file and save 
    uuid_str = str(uuid.uuid4())
    tmp_dir = os.path.join(grace_path,'normals',uuid_str)
    tmp_dir_leo = os.path.join(project_dir_leo,uuid_str)
    os.system("mkdir -p %s" % tmp_dir)
    logging.info('Temporary working directory is: %s',tmp_dir)
    logging.info('Temporary working directory on LEO is: %s',tmp_dir_leo)
    # groops bin
    groops_bin = ui.le_id218.text()
    groops_bin_cluster = ui.le_id222.text()
    if ui.cb_id079.isChecked():
        groops_bin = ui.le_id220.text()
    # xml file
    xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/build_pointmass_normals_grace_only.xml')
    if compute_goce and goce_only:
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/build_pointmass_normals_goce_only.xml')
    elif compute_goce:
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/build_pointmass_normals_grace_goce.xml')
    logging.info('XML file is: %s',xml_file)
    # copy xml file
    xml_file_tmp = os.path.join(tmp_dir,'tmp.xml')
    os.system("cp %s %s" % (xml_file,xml_file_tmp))
    # leo settings
    leo_queue = ui.cob_id001.currentText()
    tmp = ui.te_id001.time().toString("HH:mm:ss")
    h = int(tmp[0:2]) + 24*ui.sb_id003.value()
    leo_walltime = '%02d%s' % (h,tmp[2:])
    leo_scratch = ui.le_id226.text()
    leo_nodes = ui.sb_id001.value()
    leo_cores = ui.sb_id002.value()
    if not on_nodes:
        leo_nodes = 1
        leo_cores = ui.sb_id016.value()
    leo_notifications = ''
    if ui.cb_id002.isChecked():
        leo_notifications += 'a'
    if ui.cb_id003.isChecked():
        leo_notifications += 'b'
    if ui.cb_id004.isChecked():
        leo_notifications += 'e'
    leo_mail = ui.le_id003.text()
    logging.info('LEO queue: %s',leo_queue)
    logging.info('LEO walltime: %s',leo_walltime)
    logging.info('LEO scratch: %s',leo_scratch)
    logging.info('LEO nodes: %d',leo_nodes)
    logging.info('LEO cores: %d',leo_cores)
    logging.info('LEO notifications: %s',leo_notifications)
    # loop over all dates
    pbs_tmp = ()
    mjd_start, mjd_end = _date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")
        grid_file = os.path.join(grace_path,'grids','%s-grid.txt' % ymstr)   
        grid_file_tmp = os.path.join(tmp_dir,'%s-grid.txt' % ymstr)
        os.system("cp %s %s" % (grid_file,grid_file_tmp))
        sys_str = groops_interface.build_pointmass_normals(groops_bin,mjd1,mjd2,grid_file_tmp,tmp_dir,leo=True,compute_goce=compute_goce,goce_only=goce_only)
        pbs_tmp += (sys_str,)
    commands_count = len(pbs_tmp)
    total_cores = leo_nodes*leo_cores
    files_count = ui.sb_id004.value()    
    commands_per_file = int(commands_count/files_count)
    simulataneous_commands = ui.sb_id020.value()
    logging.info('LEO commands count: %d',commands_count)
    logging.info('LEO files count: %d',files_count)
    logging.info('LEO commands per file: %d',commands_per_file)
    logging.info('LEO simulataneous commands (if on frontend): %d',simulataneous_commands)
    groops_command = 'mpiexec -n %d %s ' % (total_cores,groops_bin_cluster)
    # pbs
    pbs_line = ()
    pbs_line += ('#!/bin/bash',)
    if on_nodes:
        pbs_line += ('#',
                     '#PBS -N %s' % project_name,
                     '#PBS -q %s' % leo_queue,
                     '#PBS -l walltime=%s' % leo_walltime,
                     '#PBS -d %s' % leo_scratch,
                     '#PBS -l nodes="%d:ppn=%d"' % (leo_nodes,leo_cores),    
                     '#PBS -m %s' % leo_notifications,
                     '#PBS -M %s' % leo_mail,
                     '')
    pbs_line += ('module load mpich','')
    pbs_files = ()
    for ix in list(range(files_count)):
        tmp = os.path.join(tmp_dir,'pbs%d.pbs' % ix)
        pbs_files += (tmp,)
        fw = open(tmp,'w')
        for line in pbs_line:
            fw.write('%s\n' % line)
        ixi = ix*commands_per_file
        ixj = (ix+1)*commands_per_file
        if ix==files_count-1:
            ixj = commands_count
        simul_count = 1
        for jx in list(range(ixi,ixj)):            
            simul_count += 1
            tmp = pbs_tmp[jx]
            tmp = tmp.replace(groops_bin,groops_command)
            tmp = tmp.replace(xml_file,xml_file_tmp)
            tmp = tmp.replace(data_path,os.path.join(ui.le_id224.text(),'IFG','raw'))
            tmp = tmp.replace(os.path.join(grace_path,'grids'),tmp_dir)
            tmp = tmp.replace(tmp_dir,tmp_dir_leo)            
            if simul_count<=simulataneous_commands-1 and not on_nodes:
                tmp += ' &'
            else:
                simul_count = 0
            fw.write('%s\n' % tmp)
        fw.close()     
    if picasso_on_leo:
        user = ''
        pw = ''
    else:
        user,ok = QInputDialog.getText(ui,"Enter user name for LEO","User",QLineEdit.Normal)
        pw,ok = QInputDialog.getText(ui,"Enter password for LEO","Password",QLineEdit.Password)        
    if picasso_on_leo:
        os.system("cp -r %s %s" % (tmp_dir,tmp_dir_leo)) 
    else:
        os.system("sshpass -p'%s' scp -r %s %s@leo1.iwf.oeaw.ac.at:%s" % (pw,tmp_dir,user,tmp_dir_leo)) 
    if on_nodes:
        command = '/usr/bin/qsub'
    else:
        command = '/usr/bin/bash'
    for pbs_file in pbs_files:
        if picasso_on_leo:
            os.system("%s %s" % (command,pbs_file.replace(tmp_dir,tmp_dir_leo)))
        else:
            os.system("sshpass -p'%s' ssh -X -l %s@oeaw.ads leo1.iwf.oeaw.ac.at -o StrictHostKeyChecking=no '%s %s'" % (pw,user,command,pbs_file.replace(tmp_dir,tmp_dir_leo)))
    if on_nodes:
        command = '/usr/bin/qstat'
    else:
        command = '/usr/bin/pgrep groops'
    logging.info('LEO computations started...')
    qstat = 'busy'
    while qstat:
        try:
            if picasso_on_leo:
                qstat = subprocess.check_output(command)
            else:
                qstat = subprocess.check_output(["sshpass", '-p%s' % pw, "ssh", "-X", "-l", "%s@oeaw.ads" % user, "leo1.iwf.oeaw.ac.at", "-o", "StrictHostKeyChecking=no", "%s" % command])
        except subprocess.CalledProcessError as grepexc:     
            print("error code", grepexc.returncode, grepexc.output)
            qstat = False
        time.sleep(60*5)        
    logging.info('LEO computations finished.')
    if picasso_on_leo:
        os.system("cp -r %s/* %s" % (tmp_dir_leo,normals_path))
        os.system("rm -rf %s" % (tmp_dir_leo))
    else:    
        os.system("sshpass -p'%s' scp -r %s@leo1.iwf.oeaw.ac.at:%s/* %s" % (pw,user,tmp_dir_leo,normals_path))    
        os.system("sshpass -p'%s' ssh -X -l %s@oeaw.ads leo1.iwf.oeaw.ac.at -o StrictHostKeyChecking=no 'rm -rf %s'" % (pw,user,tmp_dir_leo))
    # remove dir
    os.system("rm -rf %s" % (tmp_dir))    
    logging.info('LEO files (normals) copied back to local directory.')
    # loop over all dates
    mjd_start, mjd_end = _date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")
        grid_file = os.path.join(grace_path,'grids','%s-grid.txt' % ymstr)        
        groops_interface.combine_eliminate_solve_pointmass_normals(groops_bin,mjd1,mjd2,grid_file,normals_path,leo=False,compute_goce=compute_goce,goce_only=goce_only)    
        logging.info('Solving normals: %s.',ymstr)

def _regularization(ui,goce_only=False):
    logging.info('Regularization of normal equations...')
    # output path
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    normals_path = os.path.join(grace_path,'normals')
    # groops bin
    groops_bin = ui.le_id218.text()
    if ui.cb_id079.isChecked():
        groops_bin = ui.le_id220.text()
    # sats
    sats = ()
    if not goce_only:
        sats += ('grace',)
    sats += ('goce','grace_goce')
    # vce
    vce = True
    if not goce_only:
        vce = ui.rb_id020.isChecked()
    # matlab regularization
    matlab_regularization = None
    if current_regularization==-300:
        matlab_regularization = 'lcurve'
    elif current_regularization==-400:
        matlab_regularization = 'gcv'
    # loop over all dates
    mjd_start, mjd_end = _date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")
        for sat in sats:
            file_path = os.path.join(normals_path,'%s_%s' % (sat,ymstr))
            out_x = os.path.join(normals_path,'%s_%s-x.txt' % (sat,ymstr))
            out_sigmax = os.path.join(normals_path,'%s_%s-sigmax.txt' % (sat,ymstr))   
            try:                
                if current_regularization==-300 or current_regularization==-400:
                    LS = regularization.regularization_matlab(file_path,method=matlab_regularization)
                else:
                    LS = regularization.regularization(file_path,vce=vce)
            except:
                continue
            groops_interface.make_matrix_file(groops_bin,out_x,LS.x.elements)
            groops_interface.make_matrix_file(groops_bin,out_sigmax,LS.sigma_x.elements)
            logging.info('Regularizing normals: %s_%s',sat,ymstr)
    logging.info('Regularization of normal equations done.')

def _pointmass_grid(ui,goce_only=False):
    logging.info('Creating point mass grid files from solutions...')
    # output path
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    normals_path = os.path.join(grace_path,'normals')
    grids_path = os.path.join(grace_path,'grids')
    # groops bin
    groops_bin = ui.le_id218.text()
    if ui.cb_id079.isChecked():
        groops_bin = ui.le_id220.text()
    # sats
    sats = ()
    if not goce_only:
        sats += ('grace',)
    sats += ('grace_goce','goce')
    # loop over all dates
    mjd_start, mjd_end = _date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")
        for sat in sats:
            x_file = os.path.join(normals_path,'%s_%s-x.txt' % (sat,ymstr))
            sigmax_file = os.path.join(normals_path,'%s_%s-sigmax.txt' % (sat,ymstr))
            pointmass_file = os.path.join(grids_path,'%s_%s-pointmass.txt' % (sat,ymstr))
            grid_file = os.path.join(grids_path,'%s-grid.txt' % (ymstr))
            try:
                lon,lat,h,area = np.genfromtxt(grid_file,skip_header=2,unpack=True)
                x = np.genfromtxt(x_file,skip_header=2)
                x *= x2kg
                sigma_x = np.genfromtxt(sigmax_file,skip_header=2)
                sigma_x *= x2kg
                groops_interface.make_grid_file(groops_bin,pointmass_file,lon,lat,h,area,*(x,sigma_x))
            except:
                continue    
            logging.info('Creating grid: %s_%s-pointmass',sat,ymstr)
    logging.info('Point mass grid files created.')
    
def _time_series(ui):
    logging.info('Creating time series files from point mass grids...')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','single_solution')
    grids_path = os.path.join(grace_path,'grids')
    # time series files
    time_series_path = os.path.join(grace_path,'time_series')
    grace_file = os.path.join(time_series_path,'grace.txt')
    goce_file = os.path.join(time_series_path,'goce.txt')
    grace_goce_file = os.path.join(time_series_path,'grace_goce.txt')  
    for sat in ('grace','goce','grace_goce'):
        current_file = os.path.join(time_series_path,'%s.txt' % sat)
        with open(current_file,'w') as f:
            f.write('%25s %25s %25s\n' % ('time (mjd)','mass (kg)','sigma (kg)'))
    # loop over all dates
    mjd_start, mjd_end = _date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")  
        for sat in ('grace','goce','grace_goce'):  
            pointmass_file = os.path.join(grids_path,'%s_%s-pointmass.txt' % (sat,ymstr))
            try:
                lon,lat,h,area,x,sigma_x = np.genfromtxt(pointmass_file,skip_header=2,unpack=True)
                kg = np.sum(x)
               # sigma_kg = np.sqrt(np.sum(sigma_x**2.0)/sigma_x.size)
                sigma_kg = np.sqrt(np.sum(sigma_x**2.0))
                current_file = os.path.join(time_series_path,'%s.txt' % sat)
                with open(current_file,'a+') as f:
                    f.write('%+25.16e %+25.16e %+25.16e\n' % (mjd1,kg,sigma_kg))
            except:
                continue
            logging.info('Creating time series file: %s_%s',sat,ymstr)
    logging.info('Time series files created.')
                
def _time_series_ga(ui):
    logging.info('Creating time series files from point mass grids (GA)...')
    # number of individuals
    pop_size = ui.sb_id009.value()
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','genetic_algorithm')
    grids_path = os.path.join(grace_path,'grids')
    # time series files
    time_series_path = os.path.join(grace_path,'time_series')
    grace_file = os.path.join(time_series_path,'grace.txt')
    goce_file = os.path.join(time_series_path,'goce.txt')
    grace_goce_file = os.path.join(time_series_path,'grace_goce.txt')  
    for sat in ('grace','goce','grace_goce'):
        for individual_ix in list(range(pop_size)):
            current_file = os.path.join(time_series_path,'%s_individual%d.txt' % (sat,individual_ix))
            with open(current_file,'w') as f:
                f.write('%25s %25s %25s\n' % ('time (mjd)','mass (kg)','sigma (kg)'))
    # loop over all dates
    mjd_start, mjd_end = _date_settings(ui)
    mjd = mjd_start
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(mjd1,'%Y-%m')[0].decode("utf-8")  
        for sat in ('grace','goce','grace_goce'):  
            for individual_ix in list(range(pop_size)):
                pointmass_file = os.path.join(grids_path,'%s_%s-individual%d-pointmass.txt' % (sat,ymstr,individual_ix))
                try:
                    lon,lat,h,area,x,sigma_x = np.genfromtxt(pointmass_file,skip_header=2,unpack=True)
                    kg = np.sum(x)
                    #sigma_kg = np.sqrt(np.sum(sigma_x**2.0)/sigma_x.size)
                    sigma_kg = np.sqrt(np.sum(sigma_x**2.0))
                    current_file = os.path.join(time_series_path,'%s_individual%d.txt' % (sat,individual_ix))
                    with open(current_file,'a+') as f:
                        f.write('%+25.16e %+25.16e %+25.16e\n' % (mjd1,kg,sigma_kg))
                except:
                    continue
                logging.info('Creating time series file: %s_%s-individual%d',sat,ymstr,individual_ix)
    logging.info('Time series files (GA) created.')
        
def _optimization(ui,objective_function): 
    logging.info('Optimization started...')
    # global vars
    global current_normals_path
    global current_mjd
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace','genetic_algorithm')
    # compute normals if fixed ga, continue afterwards
    if (not pm_count_is_parameter) and (not pm_depths_are_parameters) and (not pm_positions_are_parameters) and (not pm_grid_type_is_parameter) and (not pm_grid_resolution_is_parameter):
        _no_optimization(ui)
    # set up GA configuration 
    config = {'pop_size':ui.sb_id009.value(),
              'max_iterations':ui.sb_id010.value(),
              'nsel_rate':ui.dsb_id009.value(),
              'nsel_thresh':ui.dsb_id010.value(),
              'msel_n_participants':ui.sb_id012.value(),
              'mut_rate':ui.dsb_id011.value(),
              'archive_size':ui.sb_id017.value(),
             }
    if ui.rb_id036.isChecked():
        config['nsel_method'] = 'sorted_list'
    if ui.rb_id037.isChecked():
        config['nsel_method'] = 'thresholding'
    if ui.rb_id038.isChecked():
        config['msel_method'] = 'tournament'
    if ui.rb_id039.isChecked():
        config['msel_method'] = 'sorted_list'
    if ui.rb_id040.isChecked():
        config['msel_method'] = 'random_picking'
    if ui.rb_id041.isChecked():
        config['mat_blend'] = True
    if ui.rb_id043.isChecked():
        config['mat_method'] = 'single_point'
    if ui.rb_id044.isChecked():
        config['mat_method'] = 'two_point'
    if ui.rb_id045.isChecked():
        config['mat_method'] = 'uniform'
    if ui.rb_id046.isChecked():
        config['mut_elitism'] = True
    logging.info('GA configuration:')
    for key,item in config.items():
        logging.info('%s:%s',key,item)
    # loop over all dates
    mjd_start, mjd_end = _date_settings(ui)
    mjd = mjd_start
    year,month,day = date_functions.mjd2ymd(mjd)
    while mjd<=mjd_end:
        current_mjd, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        ymstr = date_functions.mjd2ymstring(current_mjd,'%Y-%m')[0].decode("utf-8")  
        current_normals_path = os.path.join(grace_path,'normals','grace_%s-' % ymstr)
        logging.info('Optimizing month: %s',ymstr)
        config['objective_function'] = objective_function
        if isinstance(objective_function,of.MultiObjectiveFunction):
            _ga_multi_objective(ui,config)
        else:
            _ga_single_objective(ui,config)
    #    try:
    #        config['objective_function'] = objective_function
    #        if isinstance(objective_function,of.MultiObjectiveFunction):
    #            _ga_multi_objective(ui,config)
    #        else:
    #            _ga_single_objective(ui,config)
    #    except:
    #        continue
    if pm_count_is_parameter or pm_depths_are_parameters or pm_positions_are_parameters or pm_grid_type_is_parameter or pm_grid_resolution_is_parameter:
        goce_only = (pm_magnitudes_are_parameters or regularization_is_parameter)
        _build_pointmass_normals(ui,goce_only=goce_only)
        _regularization(ui,goce_only=goce_only)
        _pointmass_grid(ui,goce_only=goce_only)
    _time_series_ga(ui)
  #  if not isinstance(objective_function,of.MultiObjectiveFunction):
  #      _time_series(ui)
    _time_series(ui)
    logging.info('Optimization finished.')

def _ga_single_objective(ui,config):
    logging.info('Single Objective Optimization (SOO) started.')
    # global vars
    global ga_solutions
    soof = None
    # repeat n times until best mean population is found
    iterations = ui.sb_id011.value()
    best_population = -np.inf
    for _ in range(iterations):
        if ui.rb_id034.isChecked():   
            tmp = SOBGA.Population(config) # binary GA
        elif ui.rb_id035.isChecked(): 
            tmp = SOCGA.Population(config) # continuous GA
        tmp.start_evolution() # run GA
        fitness = tmp.development['mean_fitness'][-1]    
        if (fitness > best_population) or np.isinf(best_population): # check if current population is better than all others before
            best_population = fitness
            soof = deepcopy(tmp)
            # copy/rename/delete files 
            for ls_dict in ga_solutions:                
                if ls_dict['ga_id'] is not None:
                    normals_path = current_normals_path + ls_dict['ga_id']
                    grid_file = current_normals_path.replace('grace_','').replace('normals','grids')
                    grid_file += ls_dict['ga_id'] + '-grid.txt'
                    pointmass_file = normals_path.replace('normals','grids') + '-pointmass.txt'
                    lon,lat,h,area = np.genfromtxt(grid_file,skip_header=2,unpack=True)                    
                    x = ls_dict['LS'].x.elements.flatten()
                    x *= x2kg
                    sigma_x = ls_dict['LS'].sigma_x.elements.flatten()
                    sigma_x *= x2kg
                    ret = groops_interface.make_grid_file(current_groops_bin,pointmass_file,lon,lat,h,area,*(x,sigma_x))                    
                    if len(ls_dict['parameters'])==len(soof.optimized_parameters):
                        if np.all([ls_dict['parameters'][jx].value==soof.optimized_parameters[jx].value for jx in range(len(ls_dict['parameters']))]):
                            for file_name in glob.iglob('%s*' % normals_path):
                                shutil.copy2(file_name,file_name.replace('genetic_algorithm','single_solution'))
                            for file_name in glob.iglob('%s*' % normals_path.replace('genetic_algorithm','single_solution')):
                                os.system("mv %s %s" % (file_name,file_name.replace(ls_dict['ga_id']+'-','')))
                                if file_name.endswith('.log'):
                                    os.system("rm %s" % file_name)
                            for file_name in (grid_file,pointmass_file):
                                file_name2 = file_name.replace('genetic_algorithm','single_solution')
                                shutil.copy2(file_name,file_name2)
                                os.system("mv %s %s" % (file_name2,file_name2.replace(ls_dict['ga_id']+'-','')))
                    for individual_ix,individual in enumerate(soof.individuals):
                        if len(ls_dict['parameters'])==len(individual.genes):
                            if np.all([ls_dict['parameters'][jx].value==individual.genes[jx].value for jx in range(len(ls_dict['parameters']))]):
                                for file_name in glob.iglob('%s*' % normals_path):     
                                    os.system("mv %s %s" % (file_name,file_name.replace(ls_dict['ga_id'],'individual%d' % individual_ix)))
                                for file_name in (grid_file,pointmass_file):                                    
                                    os.system("mv %s %s" % (file_name,file_name.replace(ls_dict['ga_id'],'individual%d' % individual_ix)))
                else:
                    normals_path = current_normals_path.replace('genetic_algorithm','single_solution')
                    grid_file = normals_path.replace('grace_','').replace('normals','grids')
                    grid_file += 'grid.txt'
                    lon,lat,h,area = np.genfromtxt(grid_file,skip_header=2,unpack=True)
                    x = ls_dict['LS'].x.elements.flatten()
                    x *= x2kg
                    sigma_x = ls_dict['LS'].sigma_x.elements.flatten()
                    sigma_x *= x2kg
                    if len(ls_dict['parameters'])==len(soof.optimized_parameters):
                        if np.all([ls_dict['parameters'][jx].value==soof.optimized_parameters[jx].value for jx in range(len(ls_dict['parameters']))]):
                            pointmass_file = normals_path.replace('normals','grids') + 'pointmass.txt'
                            groops_interface.make_grid_file(current_groops_bin,pointmass_file,lon,lat,h,area,*(x,sigma_x))                    
                    for individual_ix,individual in enumerate(soof.individuals):
                        if len(ls_dict['parameters'])==len(individual.genes):
                            if np.all([ls_dict['parameters'][jx].value==individual.genes[jx].value for jx in range(len(ls_dict['parameters']))]):
                                pointmass_file = current_normals_path.replace('normals','grids') + 'individual%d-pointmass.txt' % individual_ix
                                groops_interface.make_grid_file(current_groops_bin,pointmass_file,lon,lat,h,area,*(x,sigma_x))
        for ls_dict in ga_solutions:
            if ls_dict['ga_id'] is not None:
                normals_path = current_normals_path + ls_dict['ga_id']
                grid_file = current_normals_path.replace('grace_','').replace('normals','grids')
                grid_file += ls_dict['ga_id'] + '-grid.txt'
                pointmass_file = normals_path.replace('normals','grids') + '-pointmass.txt'
                for file_name in glob.iglob('%s*' % normals_path):
                    os.system("rm %s" % file_name)
                for file_name in (grid_file,pointmass_file):
                    os.system("rm %s" % file_name)
        ga_solutions = ()    
    if not soof:
        soof = tmp       
    try:         
        fof.export_optimization( soof, current_normals_path.replace('normals','log')+'log.txt' )
    except:
        pass
    logging.info('Single Objective Optimization (SOO) finished.')

def _ga_multi_objective(ui,config):
    logging.info('Multi Objective Optimization (MOO) started.')
    # global vars
    global ga_solutions
    moof = None
    # repeat n times until best mean population is found
    iterations = ui.sb_id011.value()
    best_archive = -np.inf
    for _ in range(iterations):
        tmp = MOCGA.Population(config)
        tmp.start_evolution()
        rank0 = sum([1 for individual in tmp.archive if individual.rank==0])
        if rank0 > best_archive or np.isinf(best_archive):
            best_archive = rank0
            moof = deepcopy(tmp)
            # copy/rename/delete files 
            for ls_dict in ga_solutions:
                if ls_dict['ga_id'] is not None:
                    normals_path = current_normals_path + ls_dict['ga_id']
                    grid_file = current_normals_path.replace('grace_','').replace('normals','grids')
                    grid_file += ls_dict['ga_id'] + '-grid.txt'
                    pointmass_file = normals_path.replace('normals','grids') + '-pointmass.txt'
                    lon,lat,h,area = np.genfromtxt(grid_file,skip_header=2,unpack=True)
                    x = ls_dict['LS'].x.elements.flatten()
                    x *= x2kg
                    sigma_x = ls_dict['LS'].sigma_x.elements.flatten()
                    sigma_x *= x2kg
                    ret = groops_interface.make_grid_file(current_groops_bin,pointmass_file,lon,lat,h,area,*(x,sigma_x))
                    for individual_ix,individual in enumerate(moof.archive):
                        if len(ls_dict['parameters'])==len(individual.genes):
                            if np.all([ls_dict['parameters'][jx].value==individual.genes[jx].value for jx in range(len(ls_dict['parameters']))]):
                                for file_name in glob.iglob('%s*' % normals_path):
                                    os.system("mv %s %s" % (file_name,file_name.replace(ls_dict['ga_id'],'individual%d' % individual_ix)))           
                                for file_name in (grid_file,pointmass_file):
                                    os.system("mv %s %s" % (file_name,file_name.replace(ls_dict['ga_id'],'individual%d' % individual_ix)))
                else:
                    normals_path = current_normals_path.replace('genetic_algorithm','single_solution')
                    grid_file = normals_path.replace('grace_','').replace('normals','grids')
                    grid_file += 'grid.txt'
                    lon,lat,h,area = np.genfromtxt(grid_file,skip_header=2,unpack=True)
                    x = ls_dict['LS'].x.elements.flatten()
                    x *= x2kg
                    sigma_x = ls_dict['LS'].sigma_x.elements.flatten()
                    sigma_x *= x2kg       
                    for individual_ix,individual in enumerate(moof.archive):
                        if len(ls_dict['parameters'])==len(individual.genes):
                            if np.all([ls_dict['parameters'][jx].value==individual.genes[jx].value for jx in range(len(ls_dict['parameters']))]):
                                pointmass_file = current_normals_path.replace('normals','grids') + 'individual%d-pointmass.txt' % individual_ix
                                groops_interface.make_grid_file(current_groops_bin,pointmass_file,lon,lat,h,area,*(x,sigma_x))
#                median_parameters = []
#                for med_ix in range(len(ls_dict['parameters'])):
#                    med_tmp = []
#                    for individual in moof.archive: 
#                        if len(ls_dict['parameters'])==len(individual.genes):
#                            med_tmp.append(individual.genes[med_ix].value)
#                    median_parameters.append(np.median(med_tmp))
#                for individual_ix,individual in enumerate(moof.archive):
#                    if np.all([ls_dict['parameters'][jx].value==median_parameters[jx] for jx in range(len(ls_dict['parameters']))]):
#                        pointmass_file = current_normals_path.replace('normals','grids').replace('genetic_algorithm','single_solution') + '-pointmass.txt'
#                        groops_interface.make_grid_file(current_groops_bin,pointmass_file,lon,lat,h,area,*(x,sigma_x))
        grids = ()
        for individual_ix,individual in enumerate(moof.archive):
            pointmass_file = current_normals_path.replace('normals','grids') + 'individual%d-pointmass.txt' % individual_ix
            gridi = np.genfromtxt(pointmass_file,skip_header=2,unpack=False)
            grids += (gridi,)
        grid = np.nanmean(grids,axis=0)
        pointmass_file = current_normals_path.replace('normals','grids').replace('genetic_algorithm','single_solution') + 'pointmass.txt'
        groops_interface.make_grid_file(current_groops_bin,pointmass_file,grid[:,0],grid[:,1],grid[:,2],grid[:,3],*(grid[:,4],grid[:,5]))
        for ls_dict in ga_solutions:
            if ls_dict['ga_id'] is not None:
                normals_path = current_normals_path + ls_dict['ga_id']
                grid_file = current_normals_path.replace('grace_','').replace('normals','grids')
                grid_file += ls_dict['ga_id'] + '-grid.txt'
                pointmass_file = normals_path.replace('normals','grids') + '-pointmass.txt'
                for file_name in glob.iglob('%s*' % normals_path):
                    os.system("rm %s" % file_name)
                for file_name in (grid_file,pointmass_file):
                    os.system("rm %s" % file_name)
        ga_solutions = ()    
    if not moof:
        moof = tmp     
    if (not pm_count_is_parameter) and (not pm_depths_are_parameters) and (not pm_positions_are_parameters):
        try:
            fof.export_optimization( moof, current_normals_path.replace('normals','log')+'log.txt' )
        except:
            pass
    logging.info('Multi Objective Optimization (MOO) finished.')
        
def _update_region_area(ui):
    wgms_id = ui.cob_id002.currentText()
    rgi_id = dfr.wgms_region_code_to_rgi_region_id(wgms_id)
    area = dfr.rgi_region_id_to_area(rgi_id)
    ui.dsb_id001.setValue(area)
    
    
    
    
    
