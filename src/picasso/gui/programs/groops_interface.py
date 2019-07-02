# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 12:37:04 2018

@author: sreimond
"""

import os
import numpy as np
import pkg_resources
import shutil, tempfile
from picasso.utils.dates_and_time import date_functions
from picasso.gui.programs import compute_grace as cg
from picasso.gui.programs import compute_corrections as cc
from picasso.gui.programs import compute_internal_validation as civ

def make_grid_file(groops_bin,grid_file,lon,lat,h,area,*args,**kwargs):
    # determine dimensions
    point_count = np.size(lon)
    args_count = len(args)
    if args_count>5: # max. 5 args allowed
        args_count=5
    # determine xml file location
    xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/ascii2grid_args%d.xml' % args_count)   
    # create and fill output array
    output_ascii = np.ones((point_count,4+args_count)) * np.nan     
    output_ascii[:,0] = lon
    output_ascii[:,1] = lat
    output_ascii[:,2] = h
    output_ascii[:,3] = area    
    for ix in list(range(args_count)):
        output_ascii[:,4+ix] = args[ix]
    # replace nans and infs
    ix_nan = np.isnan(output_ascii)
    output_ascii[ix_nan] = 1e30    
    ix_inf = np.isinf(output_ascii)
    output_ascii[ix_inf] = 1e30
    # limit point count?
    if 'point_count' in kwargs:
        new_point_count = kwargs.get("point_count")
        if new_point_count>point_count or new_point_count<1:
            new_point_count = point_count
        output_ascii = output_ascii[0:new_point_count,:]
    # get temporary file and save 
    temp_dir = tempfile.mkdtemp()
    ascii_file = os.path.join(temp_dir, 'ascii.txt')
    np.savetxt(ascii_file,output_ascii,'%+.16e')
    # call groops
    sys_str = ""
    sys_str += groops_bin
    sys_str += " -g outputfileGriddedData=%s" % grid_file
    sys_str += " -g inputfile=%s" % ascii_file
    sys_str += " %s 2>/dev/null" %(xml_file)
    ret = os.system(sys_str)
    shutil.rmtree(temp_dir)
    return ret

def build_pointmass_normals(groops_bin,mjd_start,mjd_end,grid_file,output_path,leo=False,compute_goce=False,goce_only=False):
    year,month,day = date_functions.mjd2ymd(mjd_start)
    gridi = np.genfromtxt(grid_file,skip_header=2)
    gridi = np.array(gridi,ndmin=2)
    point_count = gridi.shape[0]    
    xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/build_pointmass_normals_grace_only.xml')
    if compute_goce and goce_only:
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/build_pointmass_normals_goce_only.xml')
    elif compute_goce:
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/build_pointmass_normals_grace_goce.xml')
    # system string template
    sys_str = ""
    sys_str += "xxx_grops_bin"
    sys_str += " -g timeStart=xxx_mjd_time_start"
    sys_str += " -g timeEnd=xxx_mjd_time_end"
    sys_str += " -g numberOfPoints=xxx_number_of_points"
    sys_str += " -g inputfileGrid=xxx_grid_file"
    sys_str += " -g outputfileNormalequationGraceDat=xxx_output_path_grace_xxxx-xx-normals.dat"
    sys_str += " -g outputfileNormalequationGoceDat=xxx_output_path_goce_xxxx-xx-normals.dat"
    sys_str += " -g groopsPath=xxx_groops_path"
    sys_str += " xxx_xml_file 2>/dev/null"
    # replace placeholder 
    sys_str_i = sys_str
    sys_str_i = sys_str_i.replace('xxx_grops_bin',groops_bin)
    sys_str_i = sys_str_i.replace('xxx_mjd_time_start','%08.2f' % mjd_start)
    sys_str_i = sys_str_i.replace('xxx_mjd_time_end','%08.2f' % mjd_end)
    sys_str_i = sys_str_i.replace('xxxx-xx',('%4d-%02d' % (year,month)))
    sys_str_i = sys_str_i.replace('xxx_number_of_points','%d' % point_count)
    sys_str_i = sys_str_i.replace('xxx_grid_file',grid_file)
    sys_str_i = sys_str_i.replace('xxx_output_path_',output_path+'/')
    sys_str_i = sys_str_i.replace('xxx_xml_file',xml_file)
    sys_str_i = sys_str_i.replace('xxx_groops_path',cg.data_path)
    # execute command if local, otherwise (if leo) return command string
    if not leo:
        ret = os.system(sys_str_i)
    else:
        ret = sys_str_i
    return ret

def build_pointmass_normals_ga(groops_bin,mjd_start,mjd_end,grid_file,output_path,leo=False,love_enabled=False):
    year,month,day = date_functions.mjd2ymd(mjd_start)
    gridi = np.genfromtxt(grid_file,skip_header=2)
    gridi = np.array(gridi,ndmin=2)
    point_count = gridi.shape[0]
    xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/build_pointmass_normals_grace_only_ga.xml')
    if love_enabled:
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/build_pointmass_normals_grace_only_love_enabled_ga.xml')
    # system string template
    sys_str = ""
    sys_str += "xxx_grops_bin"
    sys_str += " -g timeStart=xxx_mjd_time_start"
    sys_str += " -g timeEnd=xxx_mjd_time_end"
    sys_str += " -g numberOfPoints=xxx_number_of_points"
    sys_str += " -g inputfileGrid=xxx_grid_file"
    sys_str += " -g outputfileNormalequationGraceDat=xxx_output_path_-normals.dat"
    sys_str += " -g groopsPath=xxx_groops_path"
    sys_str += " xxx_xml_file 2>/dev/null"
    # replace placeholder 
    sys_str_i = sys_str
    sys_str_i = sys_str_i.replace('xxx_grops_bin',groops_bin)
    sys_str_i = sys_str_i.replace('xxx_mjd_time_start','%08.2f' % mjd_start)
    sys_str_i = sys_str_i.replace('xxx_mjd_time_end','%08.2f' % mjd_end)
    sys_str_i = sys_str_i.replace('xxxx-xx',('%4d-%02d' % (year,month)))
    sys_str_i = sys_str_i.replace('xxx_number_of_points','%d' % point_count)
    sys_str_i = sys_str_i.replace('xxx_grid_file',grid_file)
    sys_str_i = sys_str_i.replace('xxx_output_path_',output_path)
    sys_str_i = sys_str_i.replace('xxx_xml_file',xml_file)
    sys_str_i = sys_str_i.replace('xxx_groops_path',cg.data_path)
    # execute command if local, otherwise (if leo) return command string
    if not leo:
        ret = os.system(sys_str_i)
    else:
        ret = sys_str_i
    return ret    
    
def combine_eliminate_solve_pointmass_normals(groops_bin,mjd_start,mjd_end,grid_file,output_path,leo=False,compute_goce=False,goce_only=False):
    year,month,day = date_functions.mjd2ymd(mjd_start)
    gridi = np.genfromtxt(grid_file,skip_header=2)
    gridi = np.array(gridi,ndmin=2)
    point_count = gridi.shape[0]
    xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/combine_eliminate_solve_pointmass_normals_grace_only.xml')
    if compute_goce and goce_only:
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/combine_eliminate_solve_pointmass_normals_goce_only.xml')
    elif compute_goce:
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/combine_eliminate_solve_pointmass_normals.xml')
    # system string template
    sys_str = ""
    sys_str += "xxx_grops_bin"
    sys_str += " -g timeStart=xxx_mjd_time_start"
    sys_str += " -g timeEnd=xxx_mjd_time_end"
    sys_str += " -g numberOfPoints=xxx_number_of_points"
    sys_str += " -g inputfileGrid=xxx_grid_file"
    sys_str += " -g outputfileNormalequationGraceGoceDat=xxx_output_path_grace_goce_xxxx-xx-normals.dat"
    sys_str += " -g outputfileNormalequationGraceGoceTxt=xxx_output_path_grace_goce_xxxx-xx-normals.txt"
    sys_str += " -g outputfileNormalequationRegularizedGraceGoceTxt=xxx_output_path_grace_goce_xxxx-xx-normalsRegularized.txt"
    sys_str += " -g outputfileSolutionGraceGoce=xxx_output_path_grace_goce_xxxx-xx-x.txt"
    sys_str += " -g outputfileSigmaxGraceGoce=xxx_output_path_grace_goce_xxxx-xx-sigmax.txt"
    sys_str += " -g outputfileNormalequationGraceDat=xxx_output_path_grace_xxxx-xx-normals.dat"
    sys_str += " -g outputfileNormalequationGraceTxt=xxx_output_path_grace_xxxx-xx-normals.txt"
    sys_str += " -g outputfileNormalequationRegularizedGraceTxt=xxx_output_path_grace_xxxx-xx-normalsRegularized.txt"
    sys_str += " -g outputfileSolutionGrace=xxx_output_path_grace_xxxx-xx-x.txt"
    sys_str += " -g outputfileSigmaxGrace=xxx_output_path_grace_xxxx-xx-sigmax.txt"
    sys_str += " -g outputfileNormalequationGoceDat=xxx_output_path_goce_xxxx-xx-normals.dat"
    sys_str += " -g outputfileNormalequationGoceTxt=xxx_output_path_goce_xxxx-xx-normals.txt"
    sys_str += " -g outputfileNormalequationRegularizedGoceTxt=xxx_output_path_goce_xxxx-xx-normalsRegularized.txt"
    sys_str += " -g outputfileSolutionGoce=xxx_output_path_goce_xxxx-xx-x.txt"
    sys_str += " -g outputfileSigmaxGoce=xxx_output_path_goce_xxxx-xx-sigmax.txt"
    sys_str += " -g groopsPath=xxx_groops_path"
    sys_str += " xxx_xml_file 2>/dev/null"
    # replace placeholder 
    sys_str_i = sys_str
    sys_str_i = sys_str_i.replace('xxx_grops_bin',groops_bin)
    sys_str_i = sys_str_i.replace('xxx_mjd_time_start','%08.2f' % mjd_start)
    sys_str_i = sys_str_i.replace('xxx_mjd_time_end','%08.2f' % mjd_end)
    sys_str_i = sys_str_i.replace('xxxx-xx',('%4d-%02d' % (year,month)))
    sys_str_i = sys_str_i.replace('xxx_number_of_points','%d' % point_count)
    sys_str_i = sys_str_i.replace('xxx_grid_file',grid_file)
    sys_str_i = sys_str_i.replace('xxx_output_path_',output_path+'/')
    sys_str_i = sys_str_i.replace('xxx_xml_file',xml_file)
    sys_str_i = sys_str_i.replace('xxx_groops_path',cg.data_path)
    # execute command if local, otherwise (if leo) return command string
    if not leo:
        ret = os.system(sys_str_i)
    else:
        ret = sys_str_i
    return ret

def eliminate_solve_pointmass_normals_ga(groops_bin,mjd_start,mjd_end,grid_file,output_path,leo=False):
    year,month,day = date_functions.mjd2ymd(mjd_start)
    gridi = np.genfromtxt(grid_file,skip_header=2)
    gridi = np.array(gridi,ndmin=2)
    point_count = gridi.shape[0]
    xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/eliminate_solve_pointmass_normals_grace_only_ga.xml')
    # system string template
    sys_str = ""
    sys_str += "xxx_grops_bin"
    sys_str += " -g timeStart=xxx_mjd_time_start"
    sys_str += " -g timeEnd=xxx_mjd_time_end"
    sys_str += " -g numberOfPoints=xxx_number_of_points"
    sys_str += " -g inputfileGrid=xxx_grid_file"
    sys_str += " -g outputfileNormalequationGraceDat=xxx_output_path_-normals.dat"
    sys_str += " -g outputfileNormalequationGraceTxt=xxx_output_path_-normals.txt"
    sys_str += " -g outputfileNormalequationRegularizedGraceTxt=xxx_output_path_-normalsRegularized.txt"
    sys_str += " -g outputfileSolutionGrace=xxx_output_path_-x.txt"
    sys_str += " -g outputfileSigmaxGrace=xxx_output_path_-sigmax.txt"
    sys_str += " -g groopsPath=xxx_groops_path"
    sys_str += " xxx_xml_file 2>/dev/null"
    # replace placeholder 
    sys_str_i = sys_str
    sys_str_i = sys_str_i.replace('xxx_grops_bin',groops_bin)
    sys_str_i = sys_str_i.replace('xxx_mjd_time_start','%08.2f' % mjd_start)
    sys_str_i = sys_str_i.replace('xxx_mjd_time_end','%08.2f' % mjd_end)
    sys_str_i = sys_str_i.replace('xxxx-xx',('%4d-%02d' % (year,month)))
    sys_str_i = sys_str_i.replace('xxx_number_of_points','%d' % point_count)
    sys_str_i = sys_str_i.replace('xxx_grid_file',grid_file)
    sys_str_i = sys_str_i.replace('xxx_output_path_',output_path)
    sys_str_i = sys_str_i.replace('xxx_xml_file',xml_file)
    sys_str_i = sys_str_i.replace('xxx_groops_path',cg.data_path)
    # execute command if local, otherwise (if leo) return command string
    if not leo:
        ret = os.system(sys_str_i)
    else:
        ret = sys_str_i
    return ret    

def make_matrix_file(groops_bin,matrix_file,array):
    # determine dimensions
    array = np.array(array,ndmin=2)
    row_count = array.shape[0]
    col_count = array.shape[1]
    ele_count = row_count * col_count
    # determine xml file location
    xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/ascii2matrix.xml')   
    # create and fill output array    
    output_ascii = np.zeros((ele_count,3))
    ix = -1
    for m in list(range(row_count)):
        for n in list(range(col_count)):
            ix += 1
            output_ascii[ix,0] = m
            output_ascii[ix,1] = n
            output_ascii[ix,2] = array[m,n]
    # get temporary file and save 
    temp_dir = tempfile.mkdtemp()
    ascii_file = os.path.join(temp_dir, 'ascii.txt')
    np.savetxt(ascii_file,output_ascii,'%d %d %+.16e')
    # call groops
    sys_str = ""
    sys_str += groops_bin    
    sys_str += " -g numberColumns=%d" % col_count
    sys_str += " -g numberRows=%d" % row_count
    sys_str += " -g inputfileAscii=%s" % ascii_file
    sys_str += " -g outputfileMatrix=%s" % matrix_file
    sys_str += " %s 2>/dev/null" %(xml_file)
    ret = os.system(sys_str)
    shutil.rmtree(temp_dir)
    return ret

def compute_goco_grid(groops_bin,input_grid,output_grid,mjd):
    # determine xml file location
    xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/goco2grid.xml')   
    # call groops
    sys_str = ""
    sys_str += groops_bin
    # system string template
    sys_str = ""
    sys_str += "xxx_grops_bin"
    sys_str += " -g timeStart=xxx_mjd_time_start"
    sys_str += " -g inputfileGrid=xxx_input_grid_file"
    sys_str += " -g outputfileGriddedData=xxx_output_grid_file"
    sys_str += " -g groopsPath=xxx_groops_path"
    sys_str += " xxx_xml_file 2>/dev/null"    
    # replace placeholder 
    sys_str_i = sys_str
    sys_str_i = sys_str_i.replace('xxx_grops_bin',groops_bin)
    sys_str_i = sys_str_i.replace('xxx_mjd_time_start','%08.2f' % mjd)
    sys_str_i = sys_str_i.replace('xxx_input_grid_file',input_grid)
    sys_str_i = sys_str_i.replace('xxx_output_grid_file',output_grid)
    sys_str_i = sys_str_i.replace('xxx_xml_file',xml_file)
    sys_str_i = sys_str_i.replace('xxx_groops_path',os.path.join(cc.data_path,'IFG','raw'))
    # execute command if local
    ret = os.system(sys_str_i)
    return ret
    
def make_grid_in_polygon_file(groops_bin,output_grid,polygon_file):
    # determine xml file location
    xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/grid_in_polygon.xml')  
    # call groops
    sys_str = ""
    sys_str += groops_bin
    sys_str += " -g outputfileGrid=%s" % output_grid
    sys_str += " -g inputfilePolygon=%s" % polygon_file
    sys_str += " %s 2>/dev/null" %(xml_file)
    ret = os.system(sys_str)
    return ret

def make_specific_grid_in_polygon_file(groops_bin,output_grid,polygon_file,grid_type,grid_resolution):
    # call groops
    sys_str = ""
    sys_str += groops_bin
    sys_str += " -g outputfileGrid=%s" % output_grid
    sys_str += " -g inputfilePolygon=%s" % polygon_file
    # determine xml file location
    if grid_type==0 or grid_type=='geographical':
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/grid_in_polygon_geographical.xml')
        sys_str += " -g delta=%.10f" % grid_resolution
    elif grid_type==1 or grid_type=='triangleVertex':
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/grid_in_polygon_triangleVertex.xml')
        sys_str += " -g level=%d" % grid_resolution
    elif grid_type==2 or grid_type=='triangleCenter':
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/grid_in_polygon_triangleCenter.xml')
        sys_str += " -g level=%d" % grid_resolution
    elif grid_type==3 or grid_type=='gauss':
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/grid_in_polygon_gauss.xml')
        sys_str += " -g parallelsCount=%d" % grid_resolution
    elif grid_type==4 or grid_type=='reuter':
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/grid_in_polygon_reuter.xml')
        sys_str += " -g gamma=%d" % grid_resolution
    elif grid_type==5 or grid_type=='corput':
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/grid_in_polygon_corput.xml')
        sys_str += " -g globalPointsCount=%d" % grid_resolution
    elif grid_type==6 or grid_type=='driscoll':
        xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/grid_in_polygon_driscoll.xml') 
        sys_str += " -g dimension=%d" % grid_resolution
    # call groops
    sys_str += " %s 2>/dev/null" %(xml_file)
    ret = os.system(sys_str)
    return ret
        
def compute_gfc_grid(groops_bin,input_grid,output_grid,mjd,min_degree,max_degree,gauss,data_center='itsg'):
    # determine xml file location
    xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/%s2grid.xml' % data_center)  
    # system string template
    sys_str = ""
    sys_str += "xxx_grops_bin"
    sys_str += " -g mjd=xxx_mjd"
    sys_str += " -g inputfileGrid=xxx_grid_file"
    sys_str += " -g outputfileGriddedData=xxx_output_grid"
    sys_str += " -g groopsPath=xxx_groops_path"
    sys_str += " -g minDegree=xxx_min_degree"
    sys_str += " -g maxDegree=xxx_max_degree"
    sys_str += " -g gaussRadius=xxx_gauss"
    sys_str += " xxx_xml_file 2>/dev/null"       
    # replace placeholder 
    sys_str_i = sys_str
    sys_str_i = sys_str_i.replace('xxx_grops_bin',groops_bin)
    sys_str_i = sys_str_i.replace('xxx_mjd','%08.2f' % mjd)
    sys_str_i = sys_str_i.replace('xxx_grid_file',input_grid)
    sys_str_i = sys_str_i.replace('xxx_output_grid',output_grid)
    sys_str_i = sys_str_i.replace('xxx_min_degree','%d' % min_degree)
    sys_str_i = sys_str_i.replace('xxx_max_degree','%d' % max_degree)
    sys_str_i = sys_str_i.replace('xxx_gauss','%d' % gauss)
    sys_str_i = sys_str_i.replace('xxx_xml_file',xml_file)
    sys_str_i = sys_str_i.replace('xxx_groops_path',civ.data_path)
    # execute command if local
    ret = os.system(sys_str_i)
    return ret
    

def compute_grid_to_gfc_to_grid(groops_bin,input_gridded_data,output_grid,input_grid,gauss,min_degree=0,max_degree=60):
    # determine xml file location
    xml_file = pkg_resources.resource_filename('picasso.src.picasso.data', 'GROOPS/grid2gfc2grid.xml')  
    # system string template
    sys_str = ""
    sys_str += "xxx_grops_bin"
    sys_str += " -g inputfileGriddedData=xxx_gridded_data"
    sys_str += " -g outputfileGriddedData=xxx_output_grid"
    sys_str += " -g outputfilePotentialCoefficients=xxx_gfc_file"
    sys_str += " -g inputfileGrid=xxx_grid_file"
    sys_str += " -g groopsPath=xxx_groops_path"
    sys_str += " -g minDegree=xxx_min_degree"
    sys_str += " -g maxDegree=xxx_max_degree"
    sys_str += " -g gaussRadius=xxx_gauss"
    sys_str += " xxx_xml_file 2>/dev/null"      
    # get temporary file and save 
    temp_dir = tempfile.mkdtemp()
    gfc_file = os.path.join(temp_dir, 'tmp.gfc') 
    # replace placeholder 
    sys_str_i = sys_str
    sys_str_i = sys_str_i.replace('xxx_grops_bin',groops_bin)
    sys_str_i = sys_str_i.replace('xxx_gridded_data',input_gridded_data)
    sys_str_i = sys_str_i.replace('xxx_output_grid',output_grid)
    sys_str_i = sys_str_i.replace('xxx_grid_file',input_grid)
    sys_str_i = sys_str_i.replace('xxx_gfc_file',gfc_file)
    sys_str_i = sys_str_i.replace('xxx_min_degree','%d' % min_degree)
    sys_str_i = sys_str_i.replace('xxx_max_degree','%d' % max_degree)
    sys_str_i = sys_str_i.replace('xxx_gauss','%d' % gauss)
    sys_str_i = sys_str_i.replace('xxx_xml_file',xml_file)
    sys_str_i = sys_str_i.replace('xxx_groops_path',civ.data_path)
    # execute command if local
    ret = os.system(sys_str_i)
    shutil.rmtree(temp_dir)
    return ret    
    
    
    
    
    
    
    
    
    
