# -*- coding: utf-8 -*-
"""
Created on Fri May 25 10:46:47 2018

@author: sreimond
"""

import os
import uuid
import numpy as np
from scipy import signal
from scipy import stats
import glob
import csv
from PyQt5.QtWidgets import QMessageBox, QInputDialog, QLineEdit
from picasso.utils.constants import defaults
from picasso.utils.dates_and_time import date_functions
from picasso.utils.files import polygon_files as pof
from picasso.utils.external import lowess
from picasso.utils.algebra import matrices
from picasso.utils.statistics import least_squares
from picasso.data import data_functions_rgi as dfr
from picasso.gui import project_data_functions
import logging

regression_data = {}

def compute_regression(ui):
    # change button text
    logging.info('Regression computations started.')
    ui.pb_id005.setText('In progress...')
    ui.pb_id005.setStyleSheet("background-color: rgb(255, 250, 205); color: rgb(255,0,0);")
    # start analysis
    _load_time_series(ui)
    _identify_outliers_remove_mean(ui)
    _smoothing(ui)
    _trend_modeling(ui)
    _evaluate_trend(ui)
    _variance_analysis(ui)
    _assess_uncertainties(ui)
    _write_files(ui)
    # change back button text
    ui.pb_id005.setText('Compute Regression')
    ui.pb_id005.setStyleSheet("background-color: rgb(239, 240, 241); color: rgb(0,0,0);")
    logging.info('Regression computations finished.')
    

def _load_time_series(ui):
    logging.info('Load time series files...')
    # globals
    global regression_data
    regression_data = {}
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    # load corrections
    keys = ('a','gldas','goco','ice5','ice6','ice6_grace','klemann','lsdm','promet','wghm')
    for key in keys:
        file_path = os.path.join(project_dir,project_name,'corrections',key,'time_series','%s.txt' % key)        
        try:
            mjd,kg,kg_sigma = np.genfromtxt(file_path,skip_header=1,unpack=True)
            log_txt = 'success'
        except:
            mjd,kg,kg_sigma = np.nan, np.nan, np.nan
            log_txt = 'falied'
        logging.info('Loading file %s: %s',file_path,log_txt)
        regression_data[key] = {'mjd':np.array(mjd,ndmin=1),
                                'kg':np.array(kg,ndmin=1),
                                'kg_sigma':np.array(kg_sigma,ndmin=1)}
    # load grace
    keys = ('goce','grace','grace_goce')
    for key in keys:
        file_path = os.path.join(project_dir,project_name,'grace','single_solution','time_series','%s.txt' % key)
        try:
            mjd,kg,kg_sigma = np.genfromtxt(file_path,skip_header=1,unpack=True)
            log_txt = 'success'
        except:
            mjd,kg,kg_sigma = np.nan, np.nan, np.nan
            log_txt = 'falied'
        logging.info('Loading file %s: %s',file_path,log_txt)
        regression_data[key] = {'mjd':np.array(mjd,ndmin=1),
                                'kg':np.array(kg,ndmin=1),
                                'kg_sigma':np.array(kg_sigma,ndmin=1)}
    # combine grace and grace_goce
    mjd,kg,kg_sigma = regression_data['grace']['mjd'], regression_data['grace']['kg'], regression_data['grace']['kg_sigma']
    for mjdi in mjd:
        ix = np.where(regression_data['grace_goce']['mjd']==mjdi)
        kgj = regression_data['grace_goce']['kg'][ix]
        kg_sigmaj = regression_data['grace_goce']['kg_sigma'][ix]
        if kgj and (kgj is not None) and (not np.isnan(kgj)):
            kg[ix] = kgj
            kg_sigma[ix] = kg_sigmaj
#    kg *= 5
    regression_data['grace_goce'] = {'mjd':np.array(mjd,ndmin=1),
                                     'kg':np.array(kg,ndmin=1),
                                     'kg_sigma':np.array(kg_sigma,ndmin=1)}
    # load grace: GA population
    pop_size = ui.sb_id009.value()
    keys = ('goce','grace','grace_goce')
    for key in keys:
        for ix in list(range(pop_size)):
            key_i = '%s_individual%d' % (key,ix)
            file_path = os.path.join(project_dir,project_name,'grace','genetic_algorithm','time_series','%s.txt' % key_i)
            try:
                mjd,kg,kg_sigma = np.genfromtxt(file_path,skip_header=1,unpack=True)
                log_txt = 'success'
            except:
                mjd,kg,kg_sigma = np.nan, np.nan, np.nan
                log_txt = 'falied'
            logging.info('Loading file %s: %s',file_path,log_txt)
            regression_data[key_i] = {'mjd':np.array(mjd,ndmin=1),
                                      'kg':np.array(kg,ndmin=1),
                                      'kg_sigma':np.array(kg_sigma,ndmin=1)}
    # combine grace and grace_goce
    for ix in list(range(pop_size)):
        key_i = '%s_individual%d' % ('grace',ix)
        mjd,kg,kg_sigma = regression_data[key_i]['mjd'], regression_data[key_i]['kg'], regression_data[key_i]['kg_sigma']
        key_i = '%s_individual%d' % ('grace_goce',ix)
        for mjdi in mjd:
            jx = np.where(regression_data[key_i]['mjd']==mjdi)
            kgj = regression_data[key_i]['kg'][jx]
            kg_sigmaj = regression_data[key_i]['kg_sigma'][jx]
            if kgj and (kgj is not None) and (not np.isnan(kgj)):
                kg[jx] = kgj
                kg_sigma[jx] = kg_sigmaj
        regression_data[key_i] = {'mjd':np.array(mjd,ndmin=1),
                                  'kg':np.array(kg,ndmin=1),
                                  'kg_sigma':np.array(kg_sigma,ndmin=1)}
    # load internal validation
    keys = ('faf','iaf','m1_star','m1_syn','m1','m2_syn_star','m2','signal_loss')
    for key in keys:
        file_path = os.path.join(project_dir,project_name,'validation','internal','time_series','%s.txt' % key) 
        try:
            mjd,kg,kg_sigma = np.genfromtxt(file_path,skip_header=1,unpack=True)
            log_txt = 'success'
        except:
            mjd,kg,kg_sigma = np.nan, np.nan, np.nan
            log_txt = 'falied'
        logging.info('Loading file %s: %s',file_path,log_txt)
        regression_data[key] = {'mjd':np.array(mjd,ndmin=1),
                                'kg':np.array(kg,ndmin=1),
                                'kg_sigma':np.array(kg_sigma,ndmin=1)}
    # load external validation
    keys = ('wgms','wgms_calibrated')
    for key in keys:
        file_path = os.path.join(project_dir,project_name,'validation','external','time_series','%s.txt' % key)      
        try:
            mjd,kg,kg_sigma = np.genfromtxt(file_path,skip_header=1,unpack=True)
            log_txt = 'success'
        except:
            mjd,kg,kg_sigma = np.nan, np.nan, np.nan
            log_txt = 'falied'
        logging.info('Loading file %s: %s',file_path,log_txt)
        regression_data[key] = {'mjd':np.array(mjd,ndmin=1),
                                'kg':np.array(kg,ndmin=1),
                                'kg_sigma':np.array(kg_sigma,ndmin=1)}
    logging.info('Time series files loaded.')

def _identify_outliers_remove_mean(ui):
    logging.info('Identification of outliers and removal of mean started...')
    # globals
    global regression_data
    # time data
    date_begin = ui.de_id003.date().toString("yyyy-MM-dd")
    date_end = ui.de_id004.date().toString("yyyy-MM-dd")
    mjd_start = date_functions.ymstring2mjd(date_begin)
    tmp = date_functions.ymstring2mjd(date_end)
    _, mjd_end = date_functions.mjd2mjdmrange(tmp)
    if ui.cb_id021.isChecked():
        mjd_start = -np.inf
        mjd_end = np.inf
    # analyze data
    for key,value in regression_data.items():
        logging.info('Current time series: %s', key)
        if np.all(np.isnan(regression_data[key]['mjd'])):
            logging.info('Keeping all data.')
            regression_data[key]['mjd_use'] = regression_data[key]['mjd']
            regression_data[key]['kg_use'] = regression_data[key]['kg']
            regression_data[key]['kg_sigma_use'] = regression_data[key]['kg_sigma']
            regression_data[key]['outliers'] = []
            continue
        ix1 = regression_data[key]['mjd'] >= mjd_start
        ix2 = regression_data[key]['mjd'] <= mjd_end        
        tix = np.logical_and(ix1,ix2)
        mjd = regression_data[key]['mjd'][tix]
        kg = regression_data[key]['kg'][tix]
        kg_sigma = regression_data[key]['kg_sigma'][tix]
        if ui.rb_id074.isChecked():
            oix = np.zeros(np.count_nonzero(tix),dtype=bool)
            regression_data[key]['outliers'] = oix            
        elif not np.any(tix): 
            logging.info('No data in date range.')
            oix = np.ones(np.count_nonzero(tix),dtype=bool)
            regression_data[key]['outliers'] = oix 
        else:
            ydt = signal.detrend(kg)
            if ui.rb_id075.isChecked() or ui.rb_id076.isChecked():
                m = np.nanmean(ydt)
            if ui.rb_id077.isChecked() or ui.rb_id078.isChecked():
                m = np.nanmedian(ydt)        
            if ui.rb_id075.isChecked() or ui.rb_id077.isChecked():
                s = 3.*np.nanstd(ydt)
            if ui.rb_id076.isChecked() or ui.rb_id078.isChecked():
                s = 2.*np.nanstd(ydt)
            ix1 = ydt > (m+s)
            ix2 = ydt < (m-s)
            oix = np.logical_or(ix1,ix2)
            regression_data[key]['outliers'] = oix
        kg_outliers = np.copy(kg)
        mjd[oix] = np.nan
        kg[oix] = np.nan
        kg_sigma[oix] = np.nan
        kg_outliers[~oix] = np.nan
        regression_data[key]['kg_outliers'] = kg_outliers
        regression_data[key]['mjd_use'] = mjd
        regression_data[key]['kg_use'] = kg
        regression_data[key]['kg_sigma_use'] = kg_sigma
        logging.info('Number of outliers: %d',np.count_nonzero(oix))
        if ui.cb_id077.isChecked():
            regression_data[key]['kg_use'] -= np.nanmean(kg)
            logging.info('Tempoaral mean removed')
    logging.info('Identification of outliers and removal of mean finished.')

def _smoothing(ui):
    logging.info('Smoothing started...')
    # globals
    global regression_data
    # smooth data
    for key,value in regression_data.items():   
        logging.info('Current time series: %s', key)     
        mjd = regression_data[key]['mjd_use']
        kg  = regression_data[key]['kg_use']
        if np.all(np.isnan(mjd)):
            regression_data[key]['kg_smooth'] = regression_data[key]['kg']
            continue
        if ui.rb_id079.isChecked():
            logging.info('No smoothing applied')    
            regression_data[key]['kg_smooth'] = kg
        elif ui.rb_id082.isChecked():
            regression_data[key]['kg_smooth'] = lowess.lowess(mjd,kg,f=ui.dsb_id018.value())
            logging.info('LOESS smoothing applied, factor: %.3f',ui.dsb_id018.value())
        elif ui.rb_id080.isChecked():
            w = ui.sb_id013.value()
            tmp = np.copy(kg) * np.nan
            for ix in list(range(kg.size)):
                ixi = ix-w
                ixj = ix+w
                if ixi<0 or ixj>(kg.size-1):
                    continue
                tmp[ix] = np.nanmean(kg[ixi:ixj])
            regression_data[key]['kg_smooth'] = tmp
            logging.info('Moving average filter applied, window size: %d',w)
        elif ui.rb_id081.isChecked():
            w = ui.sb_id014.value()
            tmp = np.copy(kg) * np.nan
            for ix in list(range(kg.size)):
                ixi = ix-w
                ixj = ix+w
                if ixi<0 or ixj>(kg.size-1):
                    continue
                tmp[ix] = np.nanmedian(kg[ixi:ixj])
            regression_data[key]['kg_smooth'] = tmp
            logging.info('Moving median filter applied, window size: %d',w)
    logging.info('Smoothing finished.')
    
def _trend_modeling(ui):
    logging.info('Trend modeling started...')
    # globals
    global regression_data    
    mjd_ref = 54466.0
    T = 365.25
    # perform regression
    for key,value in regression_data.items():      
        logging.info('Current time series: %s', key)     
        mjd = regression_data[key]['mjd_use']
        if np.all(np.isnan(mjd)):
            regression_data[key]['parameter_names'] = ()
            regression_data[key]['LS'] = None
            continue
        kg  = regression_data[key]['kg_smooth']
        kg_sigma = regression_data[key]['kg_sigma_use']
        regression_data[key]['parameter_names'] = ()
        regression_data[key]['parameter_names'] += ('beta_0',)
        # check data
        if (np.any(kg_sigma==0)) or (np.any(np.isnan(kg_sigma))) or (kg_sigma is None):
            kg_sigma = np.ones(mjd.size)
        ix1 = ~np.isnan(mjd)
        ix2 = ~np.isnan(kg)
        ix = np.logical_and(ix1,ix2)
        mjd = mjd[ix]
        kg = kg[ix]
        kg_sigma = kg_sigma[ix]           
        # design matrix
        A = np.ones(mjd.size)
        t  = (mjd-mjd_ref)/T
        if ui.rb_id083.isChecked():
            degree = 1
        elif ui.rb_id084.isChecked():
            degree = 2
        elif ui.rb_id085.isChecked():
            degree = 3
        elif ui.rb_id086.isChecked():
            degree = ui.sb_id015.value()
        logging.info('Polynomial degree: %d',degree)
        for p in np.arange(1,degree+1,dtype=np.float_):
            A = np.column_stack((A,t**p))
            regression_data[key]['parameter_names'] += ('beta_%d' % p,)
            logging.info('Parameter beta_%d added', p)
        # sinusiodals
        if ui.cb_id022.isChecked():
            f = 1.0
            A = np.column_stack((A,np.cos(2.0*np.pi*f*t)))
            A = np.column_stack((A,np.sin(2.0*np.pi*f*t)))
            regression_data[key]['parameter_names'] += ('C_0','S_0',)
            logging.info('Parameters C_0 and S_0 added.')
        if ui.cb_id023.isChecked():
            f = 0.5
            A = np.column_stack((A,np.cos(2.0*np.pi*f*t)))
            A = np.column_stack((A,np.sin(2.0*np.pi*f*t)))
            regression_data[key]['parameter_names'] += ('C_1','S_1',)
            logging.info('Parameters C_1 and S_1 added.')
        if ui.cb_id024.isChecked():
            f = 161.0/T
            A = np.column_stack((A,np.cos(2.0*np.pi*f*t)))
            A = np.column_stack((A,np.sin(2.0*np.pi*f*t)))
            regression_data[key]['parameter_names'] += ('C_2','S_2',)
            logging.info('Parameters C_2 and S_2 added.')
        if ui.cb_id025.isChecked():
            f = 3.73
            A = np.column_stack((A,np.cos(2.0*np.pi*f*t)))
            A = np.column_stack((A,np.sin(2.0*np.pi*f*t)))
            regression_data[key]['parameter_names'] += ('C_3','S_3',)
            logging.info('Parameters C_3 and S_3 added.')
        if ui.cb_id026.isChecked():
            periods = np.array(ui.le_id011.text().split(),dtype=np.float_)
            for jx,period in enumerate(periods):
                f = period/T
                A = np.column_stack((A,np.cos(2.0*np.pi*f*t)))
                A = np.column_stack((A,np.sin(2.0*np.pi*f*t)))
                regression_data[key]['parameter_names'] += ('C_%d' % (jx+4),'S_%d' % (jx+4),)
                logging.info('Parameters C_%d and S_%d added.',jx+4,jx+4)
        logging.info('Building and solving LS system.')
        # covariance matrix
        Q = np.diag(kg_sigma**2.0)
        if ui.rb_id087.isChecked():
            Q = np.eye(mjd.size)
        # LS object
        A = matrices.return_matrix_from_numpy_array(A)
        Q = matrices.return_matrix_from_numpy_array(Q)
        b = matrices.return_matrix_from_numpy_array(kg)
        LS = least_squares.GeneralizedLeastSquares()
        LS.define_design_matrix(A)
        LS.define_right_hand_side_covariance_matrix(Q)
        LS.define_right_hand_side(b)
        LS.determine_normal_equation()
        LS.solve()
        LS.synthesis()
        regression_data[key]['LS'] = LS
        logging.info('LS built and solved.')
    logging.info('Trend modeling finished...')
        
def _evaluate_trend(ui):
    logging.info('Trend evaluation started...')
    # globals
    global regression_data    
    mjd_ref = 54466.0
    T = 365.25
    mjd_eval = []
    date_begin = ui.de_id003.date().toString("yyyy-MM-dd")
    date_end = ui.de_id004.date().toString("yyyy-MM-dd")
    mjd_start = date_functions.ymstring2mjd(date_begin)
    tmp = date_functions.ymstring2mjd(date_end)
    _, mjd_end = date_functions.mjd2mjdmrange(tmp)
    mjd = mjd_start    
    while mjd<=mjd_end:
        mjd1, mjd2 = date_functions.mjd2mjdmrange(mjd)
        mjd = (mjd2+1)
        mjd_eval.append(mjd1)
    mjd = np.array(mjd_eval)
    # regression
    for key,value in regression_data.items():      
        logging.info('Current time series: %s', key)     
        regression_data[key]['mjd_regression'] = mjd
        if regression_data[key]['LS'] is None:
            regression_data[key]['kg_regression'] = None
            continue
        LS = regression_data[key]['LS']
       # regression_data[key]['parameter_names'] = ()
       # regression_data[key]['parameter_names'] += ('beta_0',)
        # design matrix
        A = np.ones(mjd.size)
        t  = (mjd-mjd_ref)/T
        if ui.rb_id083.isChecked():
            degree = 1
        elif ui.rb_id084.isChecked():
            degree = 2
        elif ui.rb_id085.isChecked():
            degree = 3
        elif ui.rb_id086.isChecked():
            degree = ui.sb_id015.value()
        logging.info('Polynomial degree: %d',degree)
        for p in np.arange(1,degree+1,dtype=np.float_):
            A = np.column_stack((A,t**p))
        # sinusiodals
        if ui.cb_id022.isChecked():
            f = 1.0
            A = np.column_stack((A,np.cos(2.0*np.pi*f*t)))
            A = np.column_stack((A,np.sin(2.0*np.pi*f*t)))
            logging.info('Parameters C_0 and S_0 added.')
        if ui.cb_id023.isChecked():
            f = 0.5
            A = np.column_stack((A,np.cos(2.0*np.pi*f*t)))
            A = np.column_stack((A,np.sin(2.0*np.pi*f*t)))
            logging.info('Parameters C_1 and S_1 added.')
        if ui.cb_id024.isChecked():
            f = 161.0/T
            A = np.column_stack((A,np.cos(2.0*np.pi*f*t)))
            A = np.column_stack((A,np.sin(2.0*np.pi*f*t)))
            logging.info('Parameters C_2 and S_2 added.')
        if ui.cb_id025.isChecked():
            f = 3.73
            A = np.column_stack((A,np.cos(2.0*np.pi*f*t)))
            A = np.column_stack((A,np.sin(2.0*np.pi*f*t)))
            logging.info('Parameters C_3 and S_3 added.')
        if ui.cb_id026.isChecked():
            periods = np.array(ui.le_id011.text().split(),dtype=np.float_)
            for jx,period in enumerate(periods):
                f = period/T
                A = np.column_stack((A,np.cos(2.0*np.pi*f*t)))
                A = np.column_stack((A,np.sin(2.0*np.pi*f*t)))
                logging.info('Parameters C_%d and S_%d added.',jx+4,jx+4)
        logging.info('Evaluationg LS system.')
        # LS evaluation
        A = matrices.return_matrix_from_numpy_array(A)
        regression_data[key]['kg_regression'] = A.multiply(LS.x).elements
        logging.info('LS system evaluated.')
    logging.info('Trend evaluation finished...')
    
def _variance_analysis(ui):
    logging.info('Variance analysis started...')
    # globals
    global regression_data 
    # confidence level
    alpha = 0.05
    # perform analysis
    for key,value in regression_data.items():
        logging.info('Current time series: %s', key)     
        LS = regression_data[key]['LS']
        if LS is None:
            continue
        # basic anova table
        regression_data[key]['SS_R'] = LS.x.return_transpose().multiply(LS.A.return_transpose()).multiply(LS.b).elements - LS.row_count*np.mean(LS.b.elements)**2.0
        regression_data[key]['SS_T'] = LS.b.return_transpose().multiply(LS.b).elements - LS.row_count*np.mean(LS.b.elements)**2.0
        regression_data[key]['SS_E'] = LS.b.return_transpose().multiply(LS.b).elements - LS.x.return_transpose().multiply(LS.A.return_transpose()).multiply(LS.b).elements
        regression_data[key]['df_R'] = LS.col_count - 1
        regression_data[key]['df_E'] = LS.row_count - LS.col_count
        regression_data[key]['df_T'] = LS.row_count - 1
        regression_data[key]['MS_R'] = regression_data[key]['SS_R']/regression_data[key]['df_R']
        regression_data[key]['MS_E'] = regression_data[key]['SS_E']/regression_data[key]['df_E']
        regression_data[key]['F_0'] = regression_data[key]['MS_R']/regression_data[key]['MS_E']
        regression_data[key]['pval'] = 1.0 - stats.f.cdf(regression_data[key]['F_0'],regression_data[key]['df_R'],regression_data[key]['df_E'])
        regression_data[key]['fv'] =  stats.f.isf(1.0-alpha,regression_data[key]['df_R'],regression_data[key]['df_E'])
        regression_data[key]['R_sq'] = regression_data[key]['SS_R']/regression_data[key]['SS_T']
        regression_data[key]['R_sq_adj'] = 1.0-(regression_data[key]['SS_E']/regression_data[key]['df_E'])/(regression_data[key]['SS_T']/regression_data[key]['df_T'])
        regression_data[key]['RMSE'] = np.sqrt(LS.variance_unit_weight)
        regression_data[key]['mean_of_response'] = np.nanmean(LS.adjusted_right_hand_side.elements)
        if regression_data[key]['F_0'] > regression_data[key]['fv']:
            regression_data[key]['h0f'] = 'rejected'
        else:
            regression_data[key]['h0f'] = 'cannot be rejected'            
        # parameter estimates
        regression_data[key]['tv'] =  stats.t.isf(1.0-alpha/2.0,regression_data[key]['df_E'])
        regression_data[key]['t_ratios'] = LS.x.elements/LS.sigma_x.elements
        regression_data[key]['t_pvals'] = 2.0*(1.0-stats.t.cdf(np.fabs(regression_data[key]['t_ratios']),regression_data[key]['df_E']))
        regression_data[key]['t_txt'] = ()
        for t_ratio in regression_data[key]['t_ratios']:
            if np.fabs(t_ratio)>regression_data[key]['tv']:
                regression_data[key]['t_txt'] += ('rejected',)
            else:
                regression_data[key]['t_txt'] += ('not rejected',)
        # confidence intervals
        regression_data[key]['cil'] = LS.x.elements - LS.sigma_x.elements * regression_data[key]['tv']
        regression_data[key]['ciu'] = LS.x.elements + LS.sigma_x.elements * regression_data[key]['tv']
    logging.info('Variance analysis finished.')

def _assess_uncertainties(ui):
    logging.info('Assessing uncertainties started...')
    # globals
    global regression_data 
    regression_data['final_results'] = {'grace':np.nan,
                                        'grace_uncertainty':np.nan,
                                        'hydrology':np.nan,
                                        'hydrology_uncertainty':np.nan,
                                        'gia':np.nan,
                                        'gia_uncertainty':np.nan,
                                        'internal_validation':np.nan,
                                        'internal_validation_uncertainty':np.nan,
                                        'external_validation':np.nan,
                                        'external_validation_uncertainty':np.nan,
                                        'total':np.nan,
                                        'total_uncertainty':np.nan,}
    # GRACE
    logging.info('GRACE uncertainty:')
    if ui.rb_id105.isChecked():
        key = 'grace_goce'
    else:
        key = 'grace'
    error_grace = ()
    if ui.cb_id078.isChecked():        
        f = ui.dsb_id019.value()        
        try:
            error_grace += (regression_data[key]['LS'].sigma_x.elements[1]*f,)
        except:
            error_grace += (np.nan,)
        logging.info('added: %.2f times std. error of regression = %.2f',f,error_grace[-1])
    if ui.cb_id027.isChecked():
        f = ui.dsb_id020.value()
        pop_size = ui.sb_id009.value()
        x1 = []
        for ix in list(range(pop_size)):
            key_i = '%s_individual%d' % (key,ix)
            try:
                x1.append(regression_data[key_i]['LS'].x.elements[1])
            except:
                x1.append(np.nan)
        error_grace += (np.nanstd(x1)*f,)
        logging.info('added: %.2f times std. of population differences = %.2f',f,error_grace[-1])
    try:
        regression_data['final_results']['grace'] = regression_data[key]['LS'].x.elements[1]
        regression_data['final_results']['grace'] += regression_data['goco']['LS'].x.elements[1]
    except:
        pass
    try:
        regression_data['final_results']['grace_uncertainty'] = np.sqrt(np.nansum(np.asarray(error_grace)**2.0))
    except:
        pass
    logging.info('final results: %f +- %f',regression_data['final_results']['grace'],regression_data['final_results']['grace_uncertainty'])
    # HYDROLOGY
    logging.info('Hydrology uncertainty:')
    keys = ()
    if ui.cb_id012.isChecked():
        keys += ('gldas',)
    if ui.cb_id013.isChecked():
        keys += ('lsdm',)
    if ui.cb_id014.isChecked():
        keys += ('wghm',)
    if ui.cb_id015.isChecked():
        keys += ('promet',)
    error_hydro = ()
    if ui.cb_id028.isChecked():
        f = ui.dsb_id021.value()
        x1 = []
        for key in keys:
            try:
                x1.append(regression_data[key]['LS'].x.elements[1])
            except:
                x1.append(np.nan)
        error_hydro += (np.nanstd(x1)*f,)
        logging.info('added: %.2f times std. of model differences = %.2f',f,error_hydro[-1])
    if ui.cb_id031.isChecked():
        f = ui.dsb_id024.value()
        sigma_x = 0
        for key in keys:
            try:
                sigma_x += (regression_data[key]['LS'].sigma_x.elements[1]**2.0)
            except:
                sigma_x += 0
        error_hydro += (np.sqrt(sigma_x)*f,)
        logging.info('added: %.2f times sum of std. errors of regression analysis = %.2f',f,error_hydro[-1])
    if ui.cb_id029.isChecked():
        f = ui.dsb_id022.value()/100.0
        x1 = []
        for key in keys:
            try:
                x1.append(regression_data[key]['LS'].x.elements[1])
            except:
                x1.append(np.nan)
        error_hydro += (np.nanmean(x1)*f,)
        logging.info('added: %.2f percent of the average signal = %.2f',f*100.0,error_hydro[-1])
    if ui.cb_id030.isChecked():
        error_hydro += (ui.dsb_id023.value()*1e12,)
        logging.info('added: fixed value = %.2f',error_hydro[-1])
    x1 = []
    for key in keys:
        try:
            x1.append(regression_data[key]['LS'].x.elements[1])
        except:
            x1.append(np.nan)
    try:
        regression_data['final_results']['hydrology'] = np.nanmean(x1)
    except:
        pass
    try:
        regression_data['final_results']['hydrology_uncertainty'] = np.sqrt(np.nansum(np.asarray(error_hydro)**2.0))
    except:
        pass
    logging.info('final results: %f +- %f',regression_data['final_results']['hydrology'],regression_data['final_results']['hydrology_uncertainty'])
    # GIA    
    logging.info('GIA uncertainty:')
    keys = ()
    if ui.cb_id016.isChecked():
        keys += ('a',)
    if ui.cb_id017.isChecked():
        keys += ('ice5',)
    if ui.cb_id018.isChecked():
        keys += ('ice6',)
    if ui.cb_id019.isChecked():
        keys += ('ice6_grace',)
    if ui.cb_id020.isChecked():
        keys += ('klemann',)    
    error_gia = ()
    if ui.checkBox_31.isChecked():
        f = ui.dsb_id025.value()
        x1 = []
        for key in keys:
            try:
                x1.append(regression_data[key]['LS'].x.elements[1])
            except:
                x1.append(np.nan)
        error_gia += (np.nanstd(x1)*f,)
        logging.info('added: %.2f times std. of model differences = %.2f',f,error_gia[-1])
    if ui.cb_id032.isChecked():
        f = ui.dsb_id026.value()/100.0
        x1 = []
        for key in keys:
            try:
                x1.append(regression_data[key]['LS'].x.elements[1])
            except:
                x1.append(np.nan)
        error_gia += (np.nanmean(x1)*f,)
        logging.info('added: %.2f percent of the average signal = %.2f',f*100.0,error_gia[-1])
    if ui.cb_id033.isChecked():
        error_gia += (ui.dsb_id027.value()*1e12,)    
        logging.info('added: fixed value = %.2f',error_gia[-1])
    x1 = []
    for key in keys:
        try:
            x1.append(regression_data[key]['LS'].x.elements[1])
        except:
            x1.append(np.nan)
    try:
        regression_data['final_results']['gia'] = np.nanmean(x1)
    except:
        pass
    try:
        regression_data['final_results']['gia_uncertainty'] = np.sqrt(np.nansum(np.asarray(error_gia)**2.0))
    except:
        pass
    logging.info('final results: %f +- %f',regression_data['final_results']['gia'],regression_data['final_results']['gia_uncertainty'])
    # internal validation
    logging.info('Internal validation uncertainty:')
    error_internal_validation = ()
    if ui.cb_id034.isChecked():
        key = 'm1_syn'
        f = ui.dsb_id028.value()
        try:
            error_internal_validation += (regression_data[key]['LS'].sigma_x.elements[1],)
        except:
            error_internal_validation += (0,)
        logging.info('added: %.2f times std. error of regression = %.2f',f,error_internal_validation[-1])
    if ui.cb_id035.isChecked():
        f = ui.dsb_id029.value()
        try:
            x1_m1 = regression_data['m1']['LS'].x.elements[1]
            x1_m1_syn = regression_data['m1_syn']['LS'].x.elements[1]
            error_internal_validation += (np.fabs(x1_m1-x1_m1_syn)*f,)
        except:
            error_internal_validation += (0,)
        logging.info('added: %.2f times difference of M1 and M1_syn = %.2f',f,error_internal_validation[-1])
    try:
        regression_data['final_results']['internal_validation'] = regression_data['m1_syn']['LS'].x.elements[1]
    except:
        pass
    try:
        regression_data['final_results']['internal_validation_uncertainty'] = np.sqrt(np.nansum(np.asarray(error_internal_validation)**2.0))
    except:
        pass
    logging.info('final results: %f +- %f',regression_data['final_results']['internal_validation'],regression_data['final_results']['internal_validation_uncertainty'])
    # external validation
    logging.info('External validation uncertainty:')
    error_external_validation = ()
    if ui.cb_id036.isChecked():
        key = 'wgms'
        f = ui.dsb_id030.value()
        try:
            error_external_validation += (regression_data[key]['LS'].sigma_x.elements[1],)
        except:
            error_external_validation += (0,)
        logging.info('added: %.2f times std. error of regression = %.2f',f,error_external_validation[-1])
    if ui.cb_id037.isChecked():
        f = ui.dsb_id031.value()
        try:
            x1_m1 = regression_data['wgms']['LS'].x.elements[1]
            x1_m1_syn = regression_data['wgms_calibrated']['LS'].x.elements[1]
            error_external_validation += (np.fabs(x1_m1-x1_m1_syn)*f,)
        except:
            error_external_validation += (0,)
        logging.info('added: %.2f times difference of calibrated and non-calibrated balances = %.2f',f,error_external_validation[-1])
    try:
        regression_data['final_results']['external_validation'] = regression_data['wgms']['LS'].x.elements[1]
    except:
        pass
    try:
        regression_data['final_results']['external_validation_uncertainty'] = np.sqrt(np.nansum(np.asarray(error_external_validation)**2.0))
    except:
        pass
    logging.info('final results: %f +- %f',regression_data['final_results']['external_validation'],regression_data['final_results']['external_validation_uncertainty'])
    # total uncertainty
    logging.info('Total uncertainty:')
    error_tuple = ()
    if ui.cb_id038.isChecked():
        error_tuple += (error_grace,)
        logging.info('added: GRACE errors')
    if ui.cb_id039.isChecked():
        error_tuple += (error_hydro,)
        logging.info('added: Hydrology errors')
    if ui.cb_id040.isChecked():
        error_tuple += (error_gia,)
        logging.info('added: GIA errors')
    if ui.cb_id041.isChecked():
        error_tuple += (error_internal_validation,)  
        logging.info('added: Internal validation errors')
    try:
        tmp = [0,regression_data['final_results']['grace'],-regression_data['final_results']['hydrology'],-regression_data['final_results']['gia']]
        regression_data['final_results']['total'] = np.nansum(tmp)
        regression_data['final_results']['total_uncertainty'] = np.sqrt(np.nansum(np.asarray(error_tuple)**2.0))
    except:
        pass
    logging.info('final results: %f +- %f',regression_data['final_results']['total'],regression_data['final_results']['total_uncertainty'])
    logging.info('Assessing uncertainties finished.')

def _write_files(ui):
    logging.info('Write results to files...')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    # final results
    file_name = os.path.join(project_dir,project_name,'final_results.txt')
    f = open(file_name,'w')
    w = csv.writer(f)
    for key_i,value_i in regression_data['final_results'].items():            
        w.writerow([key_i,value_i])
    f.close()
    # corrections
    file_path = os.path.join(project_dir,project_name,'corrections')
    keys = ('a','gldas','goco','ice5','ice6','ice6_grace','klemann','lsdm','promet','wghm')
    for key in keys:
        file_name =  os.path.join(file_path,key,'regression','%s.txt' % key)
        f = open(file_name,'w')
        w = csv.writer(f)
        for key_i,value_i in regression_data[key].items():            
            w.writerow([key_i,value_i])
            if key_i is 'LS' and value_i is not None:
                w.writerow(['x',value_i.x.elements])
                w.writerow(['sigma_x',value_i.sigma_x.elements])
                w.writerow(['row_count',value_i.row_count])
                w.writerow(['adjusted_right_hand_side',value_i.adjusted_right_hand_side.elements])
        f.close()
    # grace
    file_path = os.path.join(project_dir,project_name,'grace','single_solution')
    keys = ('goce','grace','grace_goce')
    for key in keys:
        file_name =  os.path.join(file_path,'regression','%s.txt' % key)
        f = open(file_name,'w')
        w = csv.writer(f)
        for key_i,value_i in regression_data[key].items():
            w.writerow([key_i,value_i])
            if key_i is 'LS' and value_i is not None:
                w.writerow(['x',value_i.x.elements])
                w.writerow(['sigma_x',value_i.sigma_x.elements])
                w.writerow(['row_count',value_i.row_count])
                w.writerow(['adjusted_right_hand_side',value_i.adjusted_right_hand_side.elements])
        f.close()
    # GA population
    pop_size = ui.sb_id009.value()
    keys = ('goce','grace','grace_goce')
    for key in keys:
        for ix in list(range(pop_size)):
            key_i = '%s_individual%d' % (key,ix)
            file_path = os.path.join(project_dir,project_name,'grace','genetic_algorithm','regression','%s.txt' % key_i)
            f = open(file_path,'w')
            w = csv.writer(f)
            for key_j,value_j in regression_data[key_i].items():
                w.writerow([key_j,value_j])
                if key_j is 'LS' and value_j is not None:
                    w.writerow(['x',value_j.x.elements])
                    w.writerow(['sigma_x',value_j.sigma_x.elements])
                    w.writerow(['row_count',value_j.row_count])
                    w.writerow(['adjusted_right_hand_side',value_j.adjusted_right_hand_side.elements])
            f.close()        
    # internal validation
    keys = ('faf','iaf','m1_star','m1_syn','m1','m2_syn_star','m2','signal_loss')
    for key in keys:
        file_path = os.path.join(project_dir,project_name,'validation','internal','regression','%s.txt' % key) 
        f = open(file_path,'w')
        w = csv.writer(f)
        for key_i,value_i in regression_data[key].items():
            w.writerow([key_i,value_i])
            if key_i is 'LS' and value_i is not None:
                w.writerow(['x',value_i.x.elements])
                w.writerow(['sigma_x',value_i.sigma_x.elements])
                w.writerow(['row_count',value_i.row_count])
                w.writerow(['adjusted_right_hand_side',value_i.adjusted_right_hand_side.elements])
        f.close()        
    # external validation
    keys = ('wgms','wgms_calibrated')
    for key in keys:
        file_path = os.path.join(project_dir,project_name,'validation','external','regression','%s.txt' % key)      
        f = open(file_path,'w')
        w = csv.writer(f)
        for key_i,value_i in regression_data[key].items():
            w.writerow([key_i,value_i])
            if key_i is 'LS' and value_i is not None:
                w.writerow(['x',value_i.x.elements])
                w.writerow(['sigma_x',value_i.sigma_x.elements])
                w.writerow(['row_count',value_i.row_count])
                w.writerow(['adjusted_right_hand_side',value_i.adjusted_right_hand_side.elements])
        f.close()        
    logging.info('Write results to files finished.')
    
def _update_gui(ui):
    logging.info('Update GUI with regression results...')
    # paths
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    # reset to None    
    for ix in list(range(12,124)):
        # line edits
        key = 'le_id%03d' % ix
        le = getattr(ui,key)
        le.setText('None')
    #
    if ui.rb_id089.isChecked():
        key = 'grace_goce'
        file_path = os.path.join(project_dir,project_name,'grace','single_solution','regression','%s.txt' % key)      
    if ui.rb_id104.isChecked():
        key = 'grace'
        file_path = os.path.join(project_dir,project_name,'grace','single_solution','regression','%s.txt' % key)      
    if ui.rb_id111.isChecked():
        key = 'goce'
        file_path = os.path.join(project_dir,project_name,'grace','single_solution','regression','%s.txt' % key)      
    if ui.rb_id090.isChecked():
        key = 'gldas'
        file_path = os.path.join(project_dir,project_name,'corrections',key,'regression','%s.txt' % key)      
    if ui.rb_id091.isChecked():
        key = 'lsdm'
        file_path = os.path.join(project_dir,project_name,'corrections',key,'regression','%s.txt' % key)      
    if ui.rb_id092.isChecked():
        key = 'wghm'
        file_path = os.path.join(project_dir,project_name,'corrections',key,'regression','%s.txt' % key)      
    if ui.rb_id093.isChecked():
        key = 'promet'
        file_path = os.path.join(project_dir,project_name,'corrections',key,'regression','%s.txt' % key)      
    if ui.rb_id094.isChecked():
        key = 'm1_syn'
        file_path = os.path.join(project_dir,project_name,'validation','internal','regression','%s.txt' % key)      
    if ui.rb_id095.isChecked():
        key = 'wgms'
        file_path = os.path.join(project_dir,project_name,'validation','external','regression','%s.txt' % key)  
    # load data    
    results = {}
    f = open(file_path,'r')
    r = csv.reader(f)    
    for row in r:
        key = row[0]
        value = row[1]
        results[key] = value
    f.close()
    # summary of fit
    if 'R_sq' in results.keys():
        tmp = [element for element in results['R_sq'].replace('[','').replace(']','').split()]
        ui.le_id012.setText(tmp[0])
    if 'R_sq_adj' in results.keys():
        tmp = [element for element in results['R_sq_adj'].replace('[','').replace(']','').split()]
        ui.le_id013.setText(tmp[0])
    if 'RMSE' in results.keys():
        tmp = [element for element in results['RMSE'].replace('[','').replace(']','').split()]
        ui.le_id014.setText(tmp[0])
    if 'mean_of_response' in results.keys():
        tmp = [element for element in results['mean_of_response'].replace('[','').replace(']','').split()]
        ui.le_id015.setText(tmp[0])
    if 'row_count' in results.keys():
        tmp = [element for element in results['row_count'].replace('[','').replace(']','').split()]
        ui.le_id016.setText(tmp[0])
    # basic anova table
    if 'SS_R' in results.keys():
        tmp = [element for element in results['SS_R'].replace('[','').replace(']','').split()]
        ui.le_id017.setText(tmp[0])
    if 'df_R' in results.keys():
        tmp = [element for element in results['df_R'].replace('[','').replace(']','').split()]
        ui.le_id018.setText(tmp[0])
    if 'MS_R' in results.keys():
        tmp = [element for element in results['MS_R'].replace('[','').replace(']','').split()]
        ui.le_id019.setText(tmp[0])
    if 'F_0' in results.keys():
        tmp = [element for element in results['F_0'].replace('[','').replace(']','').split()]
        ui.le_id020.setText(tmp[0])
    if 'SS_E' in results.keys():
        tmp = [element for element in results['SS_E'].replace('[','').replace(']','').split()]
        ui.le_id021.setText(tmp[0])
    if 'df_E' in results.keys():
        tmp = [element for element in results['df_E'].replace('[','').replace(']','').split()]
        ui.le_id022.setText(tmp[0])
    if 'MS_E' in results.keys():
        tmp = [element for element in results['MS_E'].replace('[','').replace(']','').split()]
        ui.le_id023.setText(tmp[0])
    if 'SS_T' in results.keys():
        tmp = [element for element in results['SS_T'].replace('[','').replace(']','').split()]
        ui.le_id024.setText(tmp[0])
    if 'df_T' in results.keys():
        tmp = [element for element in results['df_T'].replace('[','').replace(']','').split()]
        ui.le_id025.setText(tmp[0])
    if 'pval' in results.keys():
        tmp = [element for element in results['pval'].replace('[','').replace(']','').split()]
        ui.le_id026.setText(tmp[0])
    if 'h0f' in results.keys():
        ui.le_id027.setText(results['h0f'])
    # parameter estimates
    if 'x' in results.keys():
        x = [element for element in results['x'].replace('[','').replace(']','').split()]
    if 'sigma_x' in results.keys():
        sigma_x = [element for element in results['sigma_x'].replace('[','').replace(']','').split()]
    if 't_ratios' in results.keys():
        t_ratios = [element for element in results['t_ratios'].replace('[','').replace(']','').split()]
    if 't_pvals' in results.keys():
        t_pvals = [element for element in results['t_pvals'].replace('[','').replace(']','').split()]
    if 't_txt' in results.keys():
        t_txt = eval(results['t_txt'])
    le_ojects = {'beta_0':(ui.le_id028,ui.le_id029,ui.le_id030,ui.le_id031,ui.le_id032),
                 'beta_1':(ui.le_id033,ui.le_id034,ui.le_id035,ui.le_id036,ui.le_id037),
                 'beta_2':(ui.le_id038,ui.le_id039,ui.le_id040,ui.le_id041,ui.le_id042),
                 'beta_3':(ui.le_id043,ui.le_id044,ui.le_id045,ui.le_id046,ui.le_id047),
                 'C_0':(ui.le_id048,ui.le_id049,ui.le_id050,ui.le_id051,ui.le_id052),
                 'S_0':(ui.le_id053,ui.le_id054,ui.le_id055,ui.le_id056,ui.le_id057),
                 'C_1':(ui.le_id058,ui.le_id059,ui.le_id060,ui.le_id061,ui.le_id062),
                 'S_1':(ui.le_id063,ui.le_id064,ui.le_id065,ui.le_id066,ui.le_id067),
                 'C_2':(ui.le_id068,ui.le_id069,ui.le_id070,ui.le_id071,ui.le_id072),
                 'S_2':(ui.le_id073,ui.le_id074,ui.le_id075,ui.le_id076,ui.le_id077),
                 'C_3':(ui.le_id078,ui.le_id079,ui.le_id080,ui.le_id081,ui.le_id082),
                 'S_3':(ui.le_id083,ui.le_id084,ui.le_id085,ui.le_id086,ui.le_id087)}    
    for ix,parameter in enumerate(eval(results['parameter_names'])):
        if not (parameter in le_ojects.keys()):
            continue
        if x:
            le_ojects[parameter][0].setText(x[ix])
        if sigma_x:
            le_ojects[parameter][1].setText(sigma_x[ix])
        if t_ratios:
            le_ojects[parameter][2].setText(t_ratios[ix])
        if t_pvals:
            le_ojects[parameter][3].setText(t_pvals[ix])
        if t_txt:
            le_ojects[parameter][4].setText(t_txt[ix])
    # confidence intervals
    if 'cil' in results.keys():
        cil = [element for element in results['cil'].replace('[','').replace(']','').split()]
    if 'ciu' in results.keys():
        ciu = [element for element in results['ciu'].replace('[','').replace(']','').split()]
    le_ojects = {'beta_0':(ui.le_id088,ui.le_id089),
                 'beta_1':(ui.le_id090,ui.le_id091),
                 'beta_2':(ui.le_id092,ui.le_id093),
                 'beta_3':(ui.le_id094,ui.le_id095),
                 'C_0':(ui.le_id096,ui.le_id097),
                 'S_0':(ui.le_id098,ui.le_id099),
                 'C_1':(ui.le_id100,ui.le_id101),
                 'S_1':(ui.le_id102,ui.le_id103),
                 'C_2':(ui.le_id104,ui.le_id105),
                 'S_2':(ui.le_id106,ui.le_id107),
                 'C_3':(ui.le_id108,ui.le_id109),
                 'S_3':(ui.le_id110,ui.le_id111)}
    for ix,parameter in enumerate(eval(results['parameter_names'])):
        if not (parameter in le_ojects.keys()):
            continue
        if cil:
            le_ojects[parameter][0].setText(cil[ix])
        if ciu:
            le_ojects[parameter][1].setText(ciu[ix])
    # final results
    file_path = os.path.join(project_dir,project_name,'final_results.txt')  
    # load data    
    results = {}
    f = open(file_path,'r')
    r = csv.reader(f)    
    for row in r:
        key = row[0]
        value = row[1]
        tmp = [element for element in value.replace('[','').replace(']','').split()]
        results[key] = '%.2f' % (float(tmp[0])*1e-12)
    f.close()
    ui.le_id112.setText(results['grace'])
    ui.le_id113.setText(results['grace_uncertainty'])
    ui.le_id114.setText(results['hydrology'])
    ui.le_id115.setText(results['hydrology_uncertainty'])
    ui.le_id116.setText(results['gia'])
    ui.le_id117.setText(results['gia_uncertainty'])
    ui.le_id118.setText(results['total'])
    ui.le_id119.setText(results['total_uncertainty'])
    ui.le_id120.setText(results['internal_validation'])
    ui.le_id121.setText(results['internal_validation_uncertainty'])
    ui.le_id122.setText(results['external_validation'])
    ui.le_id123.setText(results['external_validation_uncertainty'])
    logging.info('Update GUI with regression results done.')

    
