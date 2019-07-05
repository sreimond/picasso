# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 13:53:28 2018

@author: sreimond
"""

import os
import numpy as np
import csv
from picasso.utils.dates_and_time import date_functions
from PyQt5.QtWidgets import *
from PyQt5.QtCore import QDate, QTime
import logging

def initialize_project_settings():
    project_data = {}
    # gui data
    for ix in list(range(1,1000)):
        key = 'le_id%03d' % ix
        project_data[key] = ''
        key = 'rb_id%03d' % ix
        project_data[key] = False
        key = 'cb_id%03d' % ix
        project_data[key] = False
        key = 'cob_id%03d' % ix
        project_data[key] = 0
        key = 'sb_id%03d' % ix
        project_data[key] = 0
        key = 'dsb_id%03d' % ix
        project_data[key] = 0
        key = 'de_id%03d' % ix
        project_data[key] = ''
        key = 'te_id%03d' % ix
        project_data[key] = ''      
    logging.info('Internal projet data dictionary initialized.')
    return project_data
    
def update_project_settings(ui,project_data):
    # gui data
    for ix in list(range(1,1000)):
        # line edits
        key = 'le_id%03d' % ix
        if hasattr(ui,key):
            le = getattr(ui,key)
            project_data[key] = le.text()
        # radio buttons
        key = 'rb_id%03d' % ix
        if hasattr(ui,key):
            rb = getattr(ui,key)
            project_data[key] = rb.isChecked()
        # check box
        key = 'cb_id%03d' % ix
        if hasattr(ui,key):
            cb = getattr(ui,key)
            project_data[key] = cb.isChecked()
        # combo box
        key = 'cob_id%03d' % ix   
        if hasattr(ui,key):
            cob = getattr(ui,key)
            project_data[key] = cob.currentIndex()
        # spin box
        key = 'sb_id%03d' % ix
        if hasattr(ui,key):
            sb = getattr(ui,key)
            project_data[key] = sb.value()
        # double spin box
        key = 'dsb_id%03d' % ix
        if hasattr(ui,key):
            dsb = getattr(ui,key)
            project_data[key] = dsb.value()
        # date edit
        key = 'de_id%03d' % ix
        if hasattr(ui,key):
            de = getattr(ui,key)
            project_data[key] = de.date().toString("yyyy-MM-dd")
        # time edit   
        key = 'te_id%03d' % ix
        if hasattr(ui,key):
            te = getattr(ui,key)
            project_data[key] = te.time().toString("HH:mm:ss")
    logging.info('Internal projet data dictionary updated.')
    return project_data

def write_project_settings_to_file(project_data):
    try:
        file_path = project_data['le_id221']
        f = open(file_path,'w')
        w = csv.writer(f)
        for key,value in project_data.items():
            w.writerow([key,value])
        f.close()
        logging.info('Project settings exported to: %s',file_path)
    except:
        logging.warning('Could not export project settings to: %s',file_path)

def read_project_settings(ui,file_path):
    try:
        f = open(file_path,'r')
        r = csv.reader(f)
        for row in r:
            key = row[0]
            value = row[1]
            if not hasattr(ui,key):
                continue
            item = getattr(ui,key)
            if key.startswith('le_'):
                item.setText(value)
            elif key.startswith('rb_'):
                item.setChecked(value.lower()=="true") 
            elif key.startswith('cb_'):
                item.setChecked(value.lower()=="true") 
            elif key.startswith('cob_'):
                ix = int(value)
                item.setCurrentIndex(ix)
            elif key.startswith('sb_'):
                val = int(value)
                item.setValue(val)
            elif key.startswith('dsb_'):
                val = float(value)
                item.setValue(val)
            elif key.startswith('de'):
                date = value.split('-')
                item.setDate(QDate(int(date[0]),int(date[1]),int(date[2])))
            elif key.startswith('te'):
                time = value.split(':')
                item.setTime(QTime(int(time[0]),int(time[1]),int(time[2])))
        f.close()
        logging.info('GUI updated according to file: %s',file_path)
    except:
        logging.warning('Could not update GUI according to file: %s',file_path)

def save_project_settings(ui):
    file_path = ui.le_id221.text()
    if not file_path:
        browse_save_file(ui)
        file_path = ui.le_id221.text()
    if not file_path:
        logging.warning('Could not save settings.')
        return
    project_data = initialize_project_settings()
    project_data = update_project_settings(ui,project_data)
    write_project_settings_to_file(project_data)       
    
def load_project_settings(ui):
    file_path = browse_load_file(ui)
    if file_path == '.picpj':
        logging.info('No settings file selected.')
        return
    read_project_settings(ui,file_path)

def browse_save_file(ui,filter_string="PICASSO projects (*.picpj)"):    
    show_dir = ui.le_id002.text()
    if not show_dir:
        show_dir = os.path.expanduser("~")
    file_name, _  = QFileDialog.getSaveFileName(ui,"Save file",show_dir,filter_string)
    if not file_name.lower().endswith('.picpj'):
        file_name += '.picpj'
    if file_name == '.picpj':
        logging.info('No settings file selected.')
        return
    ui.le_id221.setText(file_name)
    
def browse_load_file(ui,filter_string="PICASSO projects (*.picpj)"):    
    show_dir = ui.le_id002.text()
    show_dir = '/home/sreimond/work/projects/paper_new'
    if not show_dir:
        show_dir = os.path.expanduser("~")
    file_name, _  = QFileDialog.getOpenFileName(ui,"Select file",show_dir,filter_string)
    if (filter_string=="PICASSO projects (*.picpj)") and not file_name.lower().endswith('.picpj'):
        file_name += '.picpj'
    return file_name

def browse_any_file(ui,filter_string="All files (*)"):
    filenames = None
    dlg = QFileDialog(ui)
    dlg.setFileMode(QFileDialog.AnyFile)
    dlg.setNameFilters([filter_string])
    if dlg.exec_():
        filenames = dlg.selectedFiles()
    if filenames:
        fileName  = filenames[0]     
        return str(fileName)
    else:
        logging.info('No file selected.')

def export_log_file(ui):
    tmp_log = logging.getLogger().handlers[0].baseFilename
    file_name = browse_any_file(ui,"LOG file (*.log)")
    if not file_name:
        logging.warning('Could not export LOG file.')
        return
    try:
        os.system("cp %s %s" % (tmp_log,file_name))
        logging.info('LOG file exported to: %s', file_name)
        logging.getLogger().handlers[0].close()
        logging.getLogger().handlers[0].baseFilename = file_name
        logging.info('LOG location changed to: %s', file_name)
    except:
        logging.warning('Could not export LOG file.')
    
def browse_folder(ui):
    show_dir = ui.le_id002.text()
    if not show_dir:
        show_dir = os.path.expanduser("~")
    directory_name = QFileDialog.getExistingDirectory(ui,"Select a directory",show_dir)
    return directory_name
    
def pick_color(le):
    color = QColorDialog.getColor()
    if color.isValid():
        le.setText(color.name())
        logging.info('Color in %s set to %s.',le.objectName(),color.name())
    
def create_project_directories(ui):
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    project_full_path = os.path.join(project_dir,project_name)
    _try_create_directory(project_full_path + '/grace')
    for sub_dir1 in ('single_solution','genetic_algorithm'):
        for sub_dir2 in ('normals','grids','time_series','regression'):
            _try_create_directory(os.path.join(project_full_path,'grace',sub_dir1,sub_dir2))
    _try_create_directory(os.path.join(project_full_path,'grace','genetic_algorithm','log'))
    _try_create_directory(project_full_path + '/corrections')
    for sub_dir1 in ('gldas','lsdm','wghm','promet','a','ice5','ice6','ice6_grace','klemann','goco'):
        for sub_dir2 in ('grids','time_series','regression'):
            _try_create_directory(os.path.join(project_full_path,'corrections',sub_dir1,sub_dir2))
    _try_create_directory(project_full_path + '/validation')
    for sub_dir1 in ('internal','external'):
        for sub_dir2 in ('grids','time_series','regression'):
            _try_create_directory(os.path.join(project_full_path,'validation',sub_dir1,sub_dir2))
    _try_create_directory(os.path.join(project_full_path,'validation','external','mass_balances'))

def _try_create_directory(dir_name):
    try:
        os.system("mkdir -p %s" % dir_name)
        logging.info('Directory created or exists already: %s',dir_name)
    except:
        logging.info('Could not create directory: %s',dir_name)

def clean_up_grace_directories(ui):
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    grace_path = os.path.join(project_dir,project_name,'grace')
    for sub_dir1 in ('single_solution','genetic_algorithm'):
        for sub_dir2 in ('normals','grids','time_series','regression'):
            _try_clean_up_directory(os.path.join(grace_path,sub_dir1,sub_dir2))
    _try_clean_up_directory(os.path.join(grace_path,'genetic_algorithm','log'))

def clean_up_correction_directories(ui):
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    corrections_path = os.path.join(project_dir,project_name,'corrections')
    for sub_dir1 in ('gldas','lsdm','wghm','promet','a','ice5','ice6','ice6_grace','klemann','goco'):
        for sub_dir2 in ('grids','time_series','regression'):
            _try_clean_up_directory(os.path.join(corrections_path,sub_dir1,sub_dir2))
            
def clean_up_internal_validation_directories(ui):
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    validation_path = os.path.join(project_dir,project_name,'validation','internal')
    for sub_dir in ('grids','time_series','regression'):
        _try_clean_up_directory(os.path.join(validation_path,sub_dir))
        
def clean_up_external_validation_directories(ui):
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    validation_path = os.path.join(project_dir,project_name,'validation','external')
    for sub_dir in ('grids','time_series','regression','mass_balances'):
        _try_clean_up_directory(os.path.join(validation_path,sub_dir))

def _try_clean_up_directory(dir_name):
    try:
        os.system("rm -rf %s/*" % dir_name)
        logging.info('All files deleted in directory: %s',dir_name)
    except:
        logging.info('Could not delete files in directory: %s',dir_name)
