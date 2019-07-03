# -*- coding: utf-8 -*-
"""
Created on Wed May  9 14:59:20 2018

@author: sreimond
"""

import os
import numpy as np
import uuid
import time
import subprocess
import shutil, tempfile
from picasso.utils.dates_and_time import date_functions
from picasso.gui.programs import groops_interface
from picasso.utils.algebra import matrices
from picasso.utils.files import groops_files
from picasso.utils.files import polygon_files as pof
from picasso.utils.geometry import points
from picasso.gui.programs import compute_grace as cg
from picasso.gui.programs import regularization as rg
LS = None

# define the objective functions
def _update_LS( parameters ):
    # global vars
    global LS
    mjd1, mjd2 = date_functions.mjd2mjdmrange(cg.current_mjd)
    # return if already updated
    LS_check = False
    for ls_dict_i in cg.ga_solutions:
        if len(ls_dict_i['parameters'])==len(parameters):
            LS_check = np.all([parameters[ix].value==ls_dict_i['parameters'][ix].value for ix in range(len(parameters))])
    if LS_check:
        LS = ls_dict_i['LS']
        return
    # ls dictionary
    ls_dict = {}
    ls_dict['parameters'] = parameters
    # temporary directory
    ls_dict['ga_id'] = str(uuid.uuid4())
    normals_path = cg.current_normals_path + ls_dict['ga_id']
    grid_file = cg.current_normals_path.replace('grace_','').replace('normals','grids')
    grid_file += ls_dict['ga_id'] + '-grid.txt'
    # possible optimization candidates
    pm_count = None
    pm_lon = []
    pm_lat = []
    pm_depth = []
    pm_magnitude = []
    pm_grid_resolution = None
    pm_grid_type = None
    regularization = 0
    # loop over parameters
    for p in parameters:
        if p.label == 'pm_count':
            pm_count = np.round(p.value).astype(int)
        elif p.label == 'pm_lon':
            pm_lon.append(p.value)
        elif p.label == 'pm_lat':
            pm_lat.append(p.value)
        elif p.label == 'pm_depth':
            pm_depth.append(p.value)
        elif p.label == 'pm_magnitude':
            pm_magnitude.append(p.value)
        elif p.label == 'regularization':
            regularization = (10.0**p.value)
        elif p.label == 'pm_grid_resolution':
            pm_grid_resolution = p.value
        elif p.label == 'pm_grid_type':
            pm_grid_type = np.floor(p.value).astype(int)
            pm_grid_type = 6 if pm_grid_type>6 else pm_grid_type
    # reassign values
    if not pm_grid_type:
        pm_grid_type = cg.current_pm_grid_type
    if pm_grid_resolution is not None:
        resolution = _map_grid_resolution( pm_grid_type, pm_grid_resolution )
        temp_dir = tempfile.mkdtemp()
        tmp_grid = os.path.join(temp_dir,'tmp.grid')
        groops_interface.make_specific_grid_in_polygon_file(cg.current_groops_bin,tmp_grid,cg.current_polygon_file,pm_grid_type,resolution)
        loni,lati,_,_ = np.genfromtxt(tmp_grid,skip_header=2,unpack=True)
        pm_lon, pm_lat = loni.tolist(), lati.tolist()
        pm_count = np.size(pm_lon)
        shutil.rmtree(temp_dir)
    if not pm_count:
        pm_count = cg.current_pm_count
    if not pm_lon:
        pm_lon = cg.current_pm_lon
    if not pm_lat:
        pm_lat = cg.current_pm_lat
    if not pm_depth:
        pm_depth = cg.current_pm_depths
    # check if (re-)computation of normals is necessary
    if (not cg.pm_count_is_parameter) and (not cg.pm_positions_are_parameters) and (not cg.pm_depths_are_parameters) and (not cg.pm_grid_type_is_parameter) and (not cg.pm_grid_resolution_is_parameter) :
        normals_path = normals_path.replace('genetic_algorithm','single_solution').replace('-'+ls_dict['ga_id'],'')
        ls_dict['ga_id'] = None
        LS = groops_files.read_groops_normals(normals_path)
    else:
        groops_interface.make_grid_file(cg.current_groops_bin,grid_file,pm_lon,pm_lat,(-np.array(pm_depth)*1e3),cg.current_pm_area,**{'point_count':pm_count})
        if not cg.picasso_on_leo:
            groops_interface.build_pointmass_normals_ga(cg.current_groops_bin,mjd1,mjd2,grid_file,normals_path,love_enabled=cg.pm_magnitudes_are_parameters)
        else:
            ret = groops_interface.build_pointmass_normals_ga(cg.current_groops_bin,mjd1,mjd2,grid_file,normals_path,leo=True,love_enabled=cg.pm_magnitudes_are_parameters)
            #groops_command = 'module load mpich ; mpiexec -n %d %s ' % (cg.current_frontend_cores,ret)
            #groops_command = 'module load mpich ; mpiexec -n %d %s ' % (40,ret.replace('bin/groops','bin/pgroops'))
            groops_command = 'mpiexec -n %d %s ' % (40,ret.replace('bin/groops','bin/pgroops'))
            #print(groops_command)
            temp_dir = tempfile.mkdtemp()
            pbs_file = os.path.join(temp_dir, 'tmp.pbs')
            #pbs_file = '/home/sreimond/data/scenario/projects/tmp.pbs'
            f = open(pbs_file,"w") 
            f.write('#!/bin/bash\n')
            f.write('#PBS -N picassoga\n')
            f.write('#PBS -q normal\n')
            f.write('#PBS -l walltime=23:59:00\n')
            f.write('#PBS -d /scratch/sreimond/myrun\n')
            f.write('#PBS -l nodes="1:ppn=40"\n')
            f.write('module load mpich\n')
            f.write(groops_command)
            f.close()            
            os.system('/usr/bin/qsub %s' % pbs_file)
            qstat = 'busy'
            while qstat:
        #        qstat = subprocess.check_output("/usr/bin/qstat")
        #        print(qstat)
                try:
                    #qstat = subprocess.check_output(["/usr/bin/pgrep","groops"])
                    qstat = subprocess.check_output("/usr/bin/qstat")
                except subprocess.CalledProcessError as grepexc:     
                    print("error code", grepexc.returncode, grepexc.output)
                    qstat = False
                time.sleep(60) 
            shutil.rmtree(temp_dir)
        groops_interface.eliminate_solve_pointmass_normals_ga(cg.current_groops_bin,mjd1,mjd2,grid_file,normals_path)
        LS = groops_files.read_groops_normals(normals_path)
    # solve
    R = matrices.IdentityMatrix(LS.col_count)
    R.elements *= regularization
    LS.define_regularization_matrix(R)
    LS.regularized_normal_equation = matrices.MatrixEquation()
    N_regularized = LS.normal_equation.A.add(LS.R)
    LS.regularized_normal_equation.define_coefficient_matrix(N_regularized)
    LS.regularized_normal_equation.define_right_hand_side(LS.normal_equation.b)
    LS.solve()
    if cg.current_regularization==-200:
        LS.determine_vce_regularization(sig1=1.0,sigu=1e2)
    elif cg.current_regularization==-300:
        LS = rg.regularization_matlab(normals_path,method='lcurve')
    elif cg.current_regularization==-400:
        LS = rg.regularization_matlab(normals_path,method='gcv')
    # update parameter vector if optimized
    if pm_magnitude:
        LS.x = matrices.return_matrix_from_numpy_array( np.array(pm_magnitude) )
        xNx = LS.x.return_transpose().multiply(LS.normal_equation.A).multiply(LS.x).elements
        nx2 = 2.0 * LS.normal_equation.b.return_transpose().multiply(LS.x).elements   
        LS.rss = xNx - nx2 + LS.bb
        xRx = LS.x.return_transpose().multiply(LS.x).elements
        LS.xss = xRx
    ls_dict['LS'] = LS
    cg.ga_solutions += (ls_dict,)
    
def _minimize_alpha( parameters ):
    try:
        _update_LS( parameters )
        alpha = LS.R.elements[0,0]
        cost = np.asscalar(alpha)
    except:
        cost = np.inf
    return cost

def _minimize_alphaxTx( parameters ):
    try:
        _update_LS( parameters )
        alpha = LS.R.elements[0,0]
        cost = alpha * np.asscalar(LS.xss)
    except:
        cost = np.inf
    return cost
    
def _minimize_condN( parameters ):
    try:
        _update_LS( parameters )
        c = np.linalg.cond(LS.regularized_normal_equation.A.elements)
        cost = np.asscalar(c)        
    except:
        cost = np.inf
    return cost

def _minimize_ePe( parameters ):
    try:
        _update_LS( parameters )
        cost = np.asscalar(LS.rss)
    except:
        cost = np.inf
    return cost
    
def _minimize_ePealphaxTx( parameters ):
    try:
        _update_LS( parameters )
        alpha = LS.R.elements[0,0]
        cost = np.asscalar(LS.rss) + alpha * np.asscalar(LS.xss)
    except:
        cost = np.inf
    return cost
       
def _minimize_sigmaxTsigmax( parameters ):
    try:
        _update_LS( parameters )
        rms = np.sqrt((LS.sigma_x.return_transpose().multiply(LS.sigma_x).elements)/LS.col_count)
        cost = np.asscalar(rms)
    except:
        cost = np.inf
    return cost

def _minimize_negatives2nratio( parameters ):
    try:
        _update_LS( parameters )
        s = np.sqrt((LS.x.return_transpose().multiply(LS.x).elements)/LS.col_count)
        n = np.sqrt((LS.sigma_x.return_transpose().multiply(LS.sigma_x).elements)/LS.col_count)
        #s = LS.x.elements
        #n = LS.sigma_x.elements
        s2n = (-1.0) * (s/n)
        #s2n = (-1.0) * np.sum(s)/np.sqrt(np.sum(n**2.0))
        #s2n = (-1.0) * np.sqrt(np.sum((s/n)**2.0)/LS.col_count)
        cost = np.asscalar(s2n)
    except:
        cost = np.inf
    return cost
    
def _minimize_xTx( parameters ):
    try:
        _update_LS( parameters )
        cost = np.asscalar(LS.xss)
    except:
        cost = np.inf
    return cost
    
def _subject_pointmass_positions( parameters ):
    pm_lon = []
    pm_lat = []
    # loop over parameters
    for p in parameters:
        if p.label == 'pm_lon':
            pm_lon.append(p.value)
        elif p.label == 'pm_lat':
            pm_lat.append(p.value)
    pip = []
    for lon,lat in zip(pm_lon,pm_lat):
        point = points.Point2D(lon,lat)
        pip.append(cg.current_polygon.determine_point_location(point))
    point_count = len(pip)
    return (np.count_nonzero( pip ) == point_count)

def _map_grid_resolution( grid_type, normalized_value ):
    do_round = True
    if grid_type==0 or grid_type=='geographical':
        do_round = False
        ll, ul = cg.current_pm_geo_ll, cg.current_pm_geo_ul
    elif grid_type==1 or grid_type=='triangleVertex':
        ll, ul = cg.current_pm_trivert_ll, cg.current_pm_trivert_ul
    elif grid_type==2 or grid_type=='triangleCenter':
        ll, ul = cg.current_pm_tricent_ll, cg.current_pm_tricent_ul
    elif grid_type==3 or grid_type=='gauss':
        ll, ul = cg.current_pm_gauss_ll, cg.current_pm_gauss_ul
    elif grid_type==4 or grid_type=='reuter':
        ll, ul = cg.current_pm_reuter_ll, cg.current_pm_reuter_ul
    elif grid_type==5 or grid_type=='corput':
        ll, ul = cg.current_pm_corput_ll, cg.current_pm_corput_ul
    elif grid_type==6 or grid_type=='driscoll':
        ll, ul = cg.current_pm_driscoll_ll, cg.current_pm_driscoll_ul
    resolution = (ul-ll)*0.5*np.sin(normalized_value*2.0*np.pi)+(ul+ll)*0.5
    if do_round:
        resolution = np.round(resolution).astype(int)
    return resolution
