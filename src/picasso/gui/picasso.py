# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 12:20:17 2018

@author: sreimond
"""

import sys, os
import threading
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox
from PyQt5 import QtCore
import functools as ft
from picasso.gui import picasso_design
from picasso.gui import project_data_functions
from picasso.gui.programs import compute_grace
from picasso.gui.programs import compute_corrections
from picasso.gui.programs import compute_internal_validation
from picasso.gui.programs import compute_external_validation
from picasso.gui.programs import compute_regression
from picasso.gui.programs import make_plot
import shutil, tempfile
import logging

class Picasso(QMainWindow, picasso_design.Ui_MainWindow):
    def __init__(self, parent=None):
        super(Picasso, self).__init__(parent)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setupUi(self)
        self.a_id001.triggered.connect(lambda: project_data_functions.save_project_settings(self))
        self.a_id002.triggered.connect(lambda: project_data_functions.load_project_settings(self))
        self.a_id004.triggered.connect(lambda: project_data_functions.export_log_file(self))
        self.pb_id008.clicked.connect(self.button_browse_project_path)
        self.pb_id042.clicked.connect(self.button_save_location)
        self.pb_id043.clicked.connect(self.button_browse_data_path)
        self.pb_id041.clicked.connect(self.button_browse_groops_bin)
        self.pb_id048.clicked.connect(self.button_browse_polygon)
        self.pb_id009.clicked.connect(self.button_browse_grid)
        self.pb_id010.clicked.connect(self.button_browse_grid_initial)
        self.pb_id052.clicked.connect(self.button_browse_grid_validation_internal)
        self.pb_id045.clicked.connect(lambda: project_data_functions.clean_up_grace_directories(self))      
        self.pb_id044.clicked.connect(lambda: project_data_functions.clean_up_correction_directories(self))    
        self.pb_id046.clicked.connect(lambda: project_data_functions.clean_up_internal_validation_directories(self))      
        self.pb_id047.clicked.connect(lambda: project_data_functions.clean_up_external_validation_directories(self))  
        self.connect_color_buttons()    
        self.pb_id001.clicked.connect(self.button_grace)
        self.pb_id002.clicked.connect(self.button_corrections)
        self.pb_id003.clicked.connect(self.button_internal_validation)
        self.pb_id004.clicked.connect(self.button_external_validation)
        self.pb_id005.clicked.connect(self.button_regression)
        self.pb_id006.clicked.connect(lambda: make_plot.make_plot(self))
        self.buttonGroup_35.buttonClicked.connect(lambda: compute_regression._update_gui(self))
        self.cob_id002.currentIndexChanged.connect(lambda: compute_grace._update_region_area(self))    
        logging.info('GUI started')
        logging.info('Temporary log file: %s',logging.getLogger().handlers[0].baseFilename)

    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Quit', 'Are you sure you want to quit?', QMessageBox.No | QMessageBox.Yes)
        if reply == QMessageBox.Yes:
            reply_save = QMessageBox.question(self, 'Save settings', 'Do you want to save the current project settings before closing?', QMessageBox.No | QMessageBox.Yes)
            if reply_save == QMessageBox.Yes:
                project_data_functions.save_project_settings(self)
            logging.info('GUI closed')
            event.accept()
        else:
            event.ignore()
                
    def button_save_location(self):
        filter_projects = "PICASSO projects (*.picpj)"
        file_name = project_data_functions.browse_any_file( self, filter_projects ) 
        if not file_name:
            logging.warning('Could not set project file location.')
            return
        if not file_name.lower().endswith('.picpj'):
            file_name += '.picpj'
        self.le_id221.setText(file_name)
        logging.info('Project file location set: %s', file_name)
    
    def button_browse_project_path(self):
        directory_name = project_data_functions.browse_folder(self)
        if not directory_name:
            logging.warning('Could not set project working directory.')
            return
        self.le_id002.setText(directory_name)   
        logging.info('Project working directory set: %s', directory_name)
        
    def button_browse_data_path(self):
        directory_name = project_data_functions.browse_folder(self)
        if not directory_name:
            logging.warning('Could not set data directory.')
            return
        self.le_id223.setText(directory_name)  
        logging.info('Data directory set: %s', directory_name)
        
    def button_browse_groops_bin(self):
        filter_projects = "GROOPS executable (*)"
        file_name = project_data_functions.browse_any_file( self, filter_projects ) 
        if not file_name:
            logging.warning('Could not set GROOPS bin location.')
            return
        self.le_id218.setText(file_name)
        logging.info('GROOPS bin location set: %s', file_name)
        
    def button_browse_polygon(self):
        filter_polygons = "GROOPS polygon file (*.xml)"
        file_name = project_data_functions.browse_load_file( self, filter_polygons ) 
        if not file_name:
            logging.warning('Could not set polygon file location.')
            return
        self.le_id004.setText(file_name)
        logging.info('Polygon file location set: %s', file_name)
        
    def button_browse_grid(self):
        filter_grids = "GROOPS grid file (*.txt *.grid)"
        file_name = project_data_functions.browse_load_file( self, filter_grids ) 
        if not file_name:
            logging.warning('Could not set grid file location.')
            return
        self.le_id005.setText(file_name)
        logging.info('Grid file location set: %s', file_name)
        
    def button_browse_grid_initial(self):
        filter_grids = "GROOPS grid file (*.txt *.grid)"
        file_name = project_data_functions.browse_load_file( self, filter_grids ) 
        if not file_name:
            logging.warning('Could not set initial grid file location.')
            return
        self.le_id006.setText(file_name)
        logging.info('Initial grid file location set: %s', file_name)
        
    def button_browse_grid_validation_internal(self):
        filter_grids = "GROOPS grid file (*.txt *.grid)"
        file_name = project_data_functions.browse_load_file( self, filter_grids ) 
        if not file_name:
            logging.warning('Could not set internal validation grid file location.')
            return
        self.le_id236.setText(file_name)
        logging.info('Internal validation grid file location set: %s', file_name)
    
    def connect_color_buttons(self):
        pb_ids = list(range(11,41)) + list(range(49,52))
        le_ids = list(range(127,206,3)) + list(range(209,216,3)) + list(range(227,234,3))
        for ix,pb_id in enumerate(pb_ids):
            key = 'pb_id%03d' % pb_id
            pb = getattr(self,key)
            key = 'le_id%03d' % le_ids[ix]
            le = getattr(self,key)
            pb.clicked.connect(ft.partial(project_data_functions.pick_color,le))
    
    def button_grace(self):
        logging.info('GRACE button pressed.')
        if self.cb_id083.isChecked():
            self.thread_grace = threading.Thread()
            self.thread_grace = threading.Thread(target=compute_grace.compute_grace,args=(self,))
            self.thread_grace.start()        
        else:
            compute_grace.compute_grace(self)            
        
    def button_corrections(self):
        logging.info('Corrections button pressed.')
        if self.cb_id083.isChecked():
            self.thread_corrections = threading.Thread()
            self.thread_corrections = threading.Thread(target=compute_corrections.compute_corrections,args=(self,))
            self.thread_corrections.start()   
        else:
            compute_corrections.compute_corrections(self)
        
    def button_internal_validation(self):
        logging.info('Internal validation button pressed.')
        if self.cb_id083.isChecked():
            self.thread_internal_validation = threading.Thread()
            self.thread_internal_validation = threading.Thread(target=compute_internal_validation.compute_internal_validation,args=(self,))
            self.thread_internal_validation.start()  
        else:
            compute_internal_validation.compute_internal_validation(self)
        
    def button_external_validation(self):
        logging.info('External validation button pressed.')
        if self.cb_id083.isChecked():
            self.thread_external_validation = threading.Thread()
            self.thread_external_validation = threading.Thread(target=compute_external_validation.compute_external_validation,args=(self,))
            self.thread_external_validation.start()  
        else:
            compute_external_validation.compute_external_validation(self)
        
    def button_regression(self):
        logging.info('Regression button pressed.')
        compute_regression.compute_regression(self)
        compute_regression._update_gui(self)
            
            
def main():    
    temp_dir = tempfile.mkdtemp()
    log_file = os.path.join(temp_dir, 'log.txt')
    logging.basicConfig(format='PICASSO (%(asctime)s, %(threadName)s): %(message)s', 
                        datefmt='%Y-%m-%d %H:%M:%S', 
                        level=logging.DEBUG,
                        filename=log_file,
                        filemode='a')
    logging.getLogger().addHandler(logging.StreamHandler())
    app = QApplication(sys.argv)
    form = Picasso()
    form.show()
    app.exec_()
    shutil.rmtree(temp_dir)
    return
  
    
    

if __name__ == '__main__':
    main()
    
