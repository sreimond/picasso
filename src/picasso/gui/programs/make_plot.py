# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 09:49:50 2018

@author: sreimond
"""

import os
import uuid
import numpy as np
import glob
import csv
#import matplotlib.pyplot as plt
import pylab as plt
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

project_name = None
project_dir = None

def make_plot(ui,skip_try=False):
    if not skip_try:
        logging.info('New plot...')
        try:
            make_plot(ui,skip_try=True)
        except:
            logging.info('Could not create plots. Check settings and result files.')
        return
    # paths
    global project_name
    global project_dir
    project_name = str(ui.le_id001.text())
    project_dir = str(ui.le_id002.text())
    # upadte rc params
    _update_rc_params(ui)
    # initiate figure
    fig,ax = plt.subplots()  
    # plot data
    _plot_grace_goce(ui,ax)
    _plot_grace(ui,ax)
    _plot_goce(ui,ax)
    _plot_hydrology(ui,ax)
    _plot_gia(ui,ax)
    _plot_internal_validation(ui,ax)
    _plot_external_validation(ui,ax)
    # ax
    title = ui.le_id124.text()
    if title:
        ax.set_title(title)
    xlabel = ui.le_id125.text()
    if xlabel:
        ax.set_xlabel(xlabel)
    ylabel = ui.le_id126.text()
    if ylabel:
        ax.set_ylabel(ylabel)
    if ui.rb_id097.isChecked():
        xl = date_functions.mjd2datetimenumber(ui.dsb_id032.value())
        xr = date_functions.mjd2datetimenumber(ui.dsb_id033.value())
        ax.set_xlim(left=xl,right=xr)
    else:
        date_begin = ui.de_id003.date().toString("yyyy-MM-dd")
        date_end = ui.de_id004.date().toString("yyyy-MM-dd")
        mjd_start = date_functions.ymstring2mjd(date_begin)
        tmp = date_functions.ymstring2mjd(date_end)
        _, mjd_end = date_functions.mjd2mjdmrange(tmp)
        xl = date_functions.mjd2datetimenumber(mjd_start)
        xr = date_functions.mjd2datetimenumber(mjd_end)
        ax.set_xlim(left=xl,right=xr)
    if ui.rb_id099.isChecked():
        ax.set_ylim(bottom=ui.dsb_id034.value(),top=ui.dsb_id035.value())
    if ui.rb_id102.isChecked():
        ax.set_yscale('linear')
    elif ui.rb_id103.isChecked():
        ax.set_yscale('log')
    ax.legend()
    fig.autofmt_xdate() 
    # export plot data
    file_path = os.path.join(project_dir,project_name,'plot_data.txt')
    f = open(file_path,'w')
    w = csv.writer(f)
    for line in ax.lines:
        x = line.get_xdata()
        x = date_functions.datetimenumber2mjd(x)
        x = date_functions.mjd2ymstring(x,'%Y-%m-%d')
        y = line.get_ydata()
        l = line.get_label()
        w.writerow(['%s_x' % l,x])
        w.writerow(['%s_y' % l,y])
    f.close()
    # show figure
    fig.show()

def _update_rc_params(ui):
    params = {
        'axes.labelsize': 8,
        'font.size': 8,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'lines.linewidth':1,
        'grid.color':'k',
        'grid.linestyle':':',
        'grid.linewidth':0.5,
        'legend.facecolor':'w',
        'legend.edgecolor':'w',
        'legend.numpoints':1,
        'legend.framealpha':1,
        'text.usetex': False,
        'figure.figsize': [6,6]
    }
    try:
        paramse = eval(le_id208.text())
        for key,value in paramse.items():
            params[key] = value
    except:
        pass
    plt.rcParams.update(params)

def _load_data(file_path):
    results = {}
    f = open(file_path,'r')
    r = csv.reader(f)    
    for row in r:
        key = row[0]
        value = row[1]
        results[key] = value
    f.close() 
    # assign arrays
    try:
        mjd = [float(element) for element in results['mjd'].replace('[','').replace(']','').split()]
    except:
        mjd = np.nan
    try:
        kg = [float(element) for element in results['kg'].replace('[','').replace(']','').split()]
    except:
        kg = np.nan
    try:
        mjd_use = [float(element) for element in results['mjd_use'].replace('[','').replace(']','').split()]
    except:
        mjd_use = np.nan
    try:
        mjd_regression = [float(element) for element in results['mjd_regression'].replace('[','').replace(']','').split()]
    except:
        mjd_regression = np.nan
    try:
        kg_use = [float(element) for element in results['kg_use'].replace('[','').replace(']','').split()]
    except:
        kg_use = np.nan
    try:
        kg_smooth = [float(element) for element in results['kg_smooth'].replace('[','').replace(']','').split()]
    except:
        kg_smooth = np.nan
    try:
        kg_regression = [float(element) for element in results['kg_regression'].replace('[','').replace(']','').split()]
    except:
        kg_regression = np.nan
    try:
        kg_adjusted = [float(element) for element in results['adjusted_right_hand_side'].replace('[','').replace(']','').split()]
    except:
        kg_adjusted = np.nan
    try:
        outliers = [float(element) for element in results['kg_outliers'].replace('[','').replace(']','').split()]
    except:
        outliers = np.nan     
    mjd = np.array(np.squeeze(mjd),ndmin=1)
    kg = np.array(np.squeeze(kg),ndmin=1)*1e-12
    mjd_use = np.array(np.squeeze(mjd_use),ndmin=1)
    mjd_regression = np.array(np.squeeze(mjd_regression),ndmin=1)
    kg_use = np.array(np.squeeze(kg_use),ndmin=1)*1e-12
    kg_smooth = np.array(np.squeeze(kg_smooth),ndmin=1)*1e-12
    kg_adjusted = np.array(np.squeeze(kg_adjusted),ndmin=1)*1e-12
    kg_regression = np.array(np.squeeze(kg_regression),ndmin=1)*1e-12
    outliers = np.array(np.squeeze(outliers),ndmin=1)*1e-12
    return mjd, kg, mjd_use, mjd_regression, kg_use, kg_smooth, kg_adjusted, kg_regression, outliers

def _plot_grace_goce(ui,ax):
    if not (ui.cb_id042.isChecked() or ui.cb_id043.isChecked() or ui.cb_id044.isChecked()):
        return     
    # load grace+goco data   
    key = 'grace_goce'
    file_path = os.path.join(project_dir,project_name,'grace','single_solution','regression','%s.txt' % key)  
    mjd, kg, mjd_use, mjd_regression, kg_use, kg_smooth, kg_adjusted, kg_regression, outliers = _load_data(file_path)
    key = 'grace'
    file_path = os.path.join(project_dir,project_name,'grace','single_solution','regression','%s.txt' % key)  
    mjd_grace, kg_grace, mjd_use_grace, mjd_regression_grace, kg_use_grace, kg_smooth_grace, kg_adjusted_grace, kg_regression_grace, outliers_grace = _load_data(file_path)
    key = 'goco'
    file_path = os.path.join(project_dir,project_name,'corrections',key,'regression','%s.txt' % key)
    mjd_goco, kg_goco, mjd_use_goco, mjd_regression_goco, kg_use_goco, kg_smooth_goco, kg_adjusted_goco, kg_regression_goco, outliers_goco = _load_data(file_path)
    # original data
    if ui.cb_id042.isChecked():
        x = np.unique(np.concatenate((mjd,mjd_grace,mjd_goco)))
        y_gg, y_g, y_go, y = np.ones(x.size)*np.nan, np.ones(x.size)*np.nan, np.ones(x.size)*np.nan, np.ones(x.size)*np.nan
        for ix,xi in enumerate(x):
            if np.isnan(xi):
                continue
            jx = mjd==xi
            if kg[jx]:
                y_gg[ix] = kg[jx]            
            jx = mjd_grace==xi
            if kg_grace[jx]:
                y_g[ix] = kg_grace[jx]
            jx = mjd_goco==xi
            if kg_goco[jx]:
                y_go[ix] = kg_goco[jx]
        for ix,xi in enumerate(x):
            if np.isnan(y_gg[ix]):
                y[ix] = y_go[ix] + y_g[ix]
            else:
                y[ix] = y_go[ix] + y_gg[ix]
        color = ui.le_id127.text()
        if not color:
            color = tuple(np.random.random(3))
        linestyle = ui.cob_id006.currentText()
        if not linestyle:
            linestyle = '-'
        label = ui.le_id128.text()
        if not label:
            label = 'GRACE+GOCO (+GOCE)'
        x = [date_functions.mjd2datetimenumber(xj) for xj in x]
        if ui.le_id129.text():
            eval(ui.le_id129.text())
        else:
            ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)
    # smoothed data
    if ui.cb_id043.isChecked():
        x = np.unique(np.concatenate((mjd_use,mjd_use_grace,mjd_use_goco)))
        y_gg, y_g, y_go, y = np.ones(x.size)*np.nan, np.ones(x.size)*np.nan, np.ones(x.size)*np.nan, np.ones(x.size)*np.nan
        for ix,xi in enumerate(x):
            if np.isnan(xi):
                continue
            jx = mjd_use==xi
            if kg_smooth[jx]:
                y_gg[ix] = kg_smooth[jx]            
            jx = mjd_use_grace==xi
            if kg_smooth_grace[jx]:
                y_g[ix] = kg_smooth_grace[jx]
            jx = mjd_use_goco==xi
            if kg_smooth_goco[jx]:
                y_go[ix] = kg_smooth_goco[jx]
        for ix,xi in enumerate(x):
            if np.isnan(y_gg[ix]):
                y[ix] = y_go[ix] + y_g[ix]
            else:
                y[ix] = y_go[ix] + y_gg[ix]
        color = ui.le_id130.text()
        if not color:
            color = tuple(np.random.random(3))
        linestyle = ui.cob_id007.currentText()
        if not linestyle:
            linestyle = '-'
        label = ui.le_id131.text()
        if not label:
            label = 'GRACE+GOCO (+GOCE) smoothed'
        x = [date_functions.mjd2datetimenumber(xj) for xj in x]
        if ui.le_id132.text():
            eval(ui.le_id132.text())
        else:
            ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)
    # regression model
    if ui.cb_id044.isChecked():
        x = np.unique(np.concatenate((mjd_regression,mjd_regression_grace,mjd_regression_goco)))
        y_gg, y_g, y_go, y = np.ones(x.size)*np.nan, np.ones(x.size)*np.nan, np.ones(x.size)*np.nan, np.ones(x.size)*np.nan
        for ix,xi in enumerate(x):
            if np.isnan(xi):
                continue
            jx = mjd_regression==xi
            if kg_regression[jx]:
                y_gg[ix] = kg_regression[jx]   
            jx = mjd_regression_grace==xi
            if kg_regression_grace[jx]:
                y_g[ix] = kg_regression_grace[jx]
            jx = mjd_regression_goco==xi
            if kg_regression_goco[jx]:
                y_go[ix] = kg_regression_goco[jx]
        for ix,xi in enumerate(x):
            if np.isnan(y_gg[ix]):
                y[ix] = y_go[ix] + y_g[ix]
            else:
                y[ix] = y_go[ix] + y_gg[ix]
        color = ui.le_id133.text()
        if not color:
            color = tuple(np.random.random(3))
        linestyle = ui.cob_id008.currentText()
        if not linestyle:
            linestyle = '-'
        label = ui.le_id134.text()
        if not label:
            label = 'GRACE+GOCO (+GOCE) regression'
        x = [date_functions.mjd2datetimenumber(xj) for xj in x]
        if ui.le_id135.text():
            eval(ui.le_id135.text())
        else:
            ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)

def _plot_grace(ui,ax):
    if not (ui.cb_id046.isChecked() or ui.cb_id047.isChecked() or ui.cb_id048.isChecked() or ui.cb_id045.isChecked()):
        return     
    # load grace+goco data 
    key = 'grace'
    file_path = os.path.join(project_dir,project_name,'grace','single_solution','regression','%s.txt' % key)      
    mjd, kg, mjd_use, mjd_regression, kg_use, kg_smooth, kg_adjusted, kg_regression, outliers = _load_data(file_path)
    mjd_ga, kg_ga, mjd_use_ga, mjd_regression_ga, kg_use_ga, kg_smooth_ga, kg_adjusted_ga, kg_regression_ga, outliers_ga = (), (), (), (), (), (), (), (), ()
    pop_size = ui.sb_id009.value()
    for ix in list(range(pop_size)):
        key_i = '%s_individual%d' % (key,ix)
        file_path = os.path.join(project_dir,project_name,'grace','genetic_algorithm','regression','%s.txt' % key_i)  
        mjd_i, kg_i, mjd_use_i, mjd_regression_i, kg_use_i, kg_smooth_i, kg_adjusted_i, kg_regression_i, outliers_i = _load_data(file_path)
        mjd_ga += (mjd_i,)
        kg_ga += (kg_i,)
        mjd_use_ga += (mjd_use_i,)
        mjd_regression_ga += (mjd_regression_i,)
        kg_use_ga += (kg_use_i,)
        kg_smooth_ga += (kg_smooth_i,)
        kg_adjusted_ga += (kg_adjusted_i,)
        kg_regression_ga += (kg_regression_i,)
        outliers_ga += (outliers_i,)
    # original data
    if ui.cb_id046.isChecked():
        x = mjd
        y = kg
        xo = mjd_use
        yo = outliers
        ix = np.isnan(xo)
        xo = xo[~ix]
        yo = yo[~ix]
        color = ui.le_id209.text()
        if not color:
            color = tuple(np.random.random(3))
        linestyle = ui.cob_id010.currentText()
        if not linestyle:
            linestyle = '-'
        label = ui.le_id210.text()
        if not label:
            label = 'GRACE'
        x = [date_functions.mjd2datetimenumber(xj) for xj in x]
        if ui.le_id211.text():
            eval(ui.le_id211.text())
        else:
            ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)
            if not np.all(np.isnan(outliers)):
                xo = [date_functions.mjd2datetimenumber(xj) for xj in xo]
                ax.plot_date(xo,yo,color=color,linestyle='',marker='o',markerfacecolor='w',label=label+' (outliers)')
    # smoothed data
    if ui.cb_id047.isChecked():
        x = mjd_use
        y = kg_smooth
        color = ui.le_id212.text()
        if not color:
            color = tuple(np.random.random(3))
        linestyle = ui.cob_id011.currentText()
        if not linestyle:
            linestyle = '-'
        label = ui.le_id213.text()
        if not label:
            label = 'GRACE smoothed'
        x = [date_functions.mjd2datetimenumber(xj) for xj in x]
        if ui.le_id214.text():
            eval(ui.le_id214.text())
        else:
            ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)
    # regression model
    if ui.cb_id048.isChecked():
        x = mjd_regression
        y = kg_regression
        color = ui.le_id215.text()
        if not color:
            color = tuple(np.random.random(3))
        linestyle = ui.cob_id012.currentText()
        if not linestyle:
            linestyle = '-'
        label = ui.le_id216.text()
        if not label:
            label = 'GRACE regression'
        x = [date_functions.mjd2datetimenumber(xj) for xj in x]
        if ui.le_id217.text():
            eval(ui.le_id217.text())
        else:
            ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)
    # GA population
    if ui.cb_id045.isChecked():
        color = ui.le_id136.text()
        if not color:
            color = tuple(np.random.random(3))
        linestyle = ui.cob_id009.currentText()
        if not linestyle:
            linestyle = '-'
        label = ui.le_id137.text()
        if not label:
            label = 'GRACE GA population'
        for ix,(x,y) in enumerate(zip(mjd_ga,kg_ga)):
            x = [date_functions.mjd2datetimenumber(xj) for xj in x]
            if ix>0:
                label = '_nolegend_'
            if ui.le_id138.text():
                eval(ui.le_id138.text())
            else:
                ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)
        
def _plot_goce(ui,ax):
    if not (ui.cb_id080.isChecked() or ui.cb_id081.isChecked() or ui.cb_id082.isChecked()):
        return     
    # load goce data
    key = 'goce'
    file_path = os.path.join(project_dir,project_name,'grace','single_solution','regression','%s.txt' % key)  
    mjd, kg, mjd_use, mjd_regression, kg_use, kg_smooth, kg_adjusted, kg_regression, outliers = _load_data(file_path)
    # original data
    if ui.cb_id080.isChecked():
        x = mjd
        y = kg
        xo = mjd_use
        yo = outliers
        ix = np.isnan(xo)
        xo = xo[~ix]
        yo = yo[~ix]
        color = ui.le_id227.text()
        if not color:
            color = tuple(np.random.random(3))
        linestyle = ui.cob_id010.currentText()
        if not linestyle:
            linestyle = '-'
        label = ui.le_id228.text()
        if not label:
            label = 'GOCE'
        x = [date_functions.mjd2datetimenumber(xj) for xj in x]
        if ui.le_id229.text():
            eval(ui.le_id229.text())
        else:
            ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)
            if not np.all(np.isnan(outliers)):
                xo = [date_functions.mjd2datetimenumber(xj) for xj in xo]
                ax.plot_date(xo,yo,color=color,linestyle='',marker='o',markerfacecolor='w',label=label+' (outliers)')
    # smoothed data
    if ui.cb_id081.isChecked():
        x = mjd_use
        y = kg_smooth
        color = ui.le_id230.text()
        if not color:
            color = tuple(np.random.random(3))
        linestyle = ui.cob_id011.currentText()
        if not linestyle:
            linestyle = '-'
        label = ui.le_id231.text()
        if not label:
            label = 'GOCE smoothed'
        x = [date_functions.mjd2datetimenumber(xj) for xj in x]
        if ui.le_id232.text():
            eval(ui.le_id232.text())
        else:
            ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label) 
    # regression model
    if ui.cb_id082.isChecked():
        x = mjd_regression
        y = kg_regression
        color = ui.le_id233.text()
        if not color:
            color = tuple(np.random.random(3))
        linestyle = ui.cob_id012.currentText()
        if not linestyle:
            linestyle = '-'
        label = ui.le_id234.text()
        if not label:
            label = 'GOCE regression'
        x = [date_functions.mjd2datetimenumber(xj) for xj in x]
        if ui.le_id235.text():
            eval(ui.le_id235.text())
        else:
            ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)   
            
def _plot_hydrology(ui,ax):
    keys = ('gldas','lsdm','wghm','promet')
    cb = {'gldas':(ui.cb_id049,ui.cb_id050,ui.cb_id051),
          'lsdm':(ui.cb_id052,ui.cb_id053,ui.cb_id054),
          'wghm':(ui.cb_id055,ui.cb_id056,ui.cb_id057),
          'promet':(ui.cb_id058,ui.cb_id059,ui.cb_id060)}
    col_le = {'gldas':(ui.le_id139,ui.le_id142,ui.le_id145),
              'lsdm':(ui.le_id148,ui.le_id151,ui.le_id154),
              'wghm':(ui.le_id157,ui.le_id160,ui.le_id163),
              'promet':(ui.le_id166,ui.le_id169,ui.le_id172)}
    ls_le = {'gldas':(ui.cob_id013,ui.cob_id014,ui.cob_id015),
              'lsdm':(ui.cob_id016,ui.cob_id017,ui.cob_id018),
              'wghm':(ui.cob_id019,ui.cob_id020,ui.cob_id021),
              'promet':(ui.cob_id022,ui.cob_id023,ui.cob_id024)}
    label_le = {'gldas':(ui.le_id140,ui.le_id143,ui.le_id146),
              'lsdm':(ui.le_id149,ui.le_id152,ui.le_id155),
              'wghm':(ui.le_id158,ui.le_id161,ui.le_id164),
              'promet':(ui.le_id167,ui.le_id170,ui.le_id173)}
    cmd_le = {'gldas':(ui.le_id141,ui.le_id144,ui.le_id147),
              'lsdm':(ui.le_id150,ui.le_id153,ui.le_id156),
              'wghm':(ui.le_id159,ui.le_id162,ui.le_id165),
              'promet':(ui.le_id168,ui.le_id171,ui.le_id174)}
    for key in keys:
        # load correction data
        file_path = os.path.join(project_dir,project_name,'corrections',key,'regression','%s.txt' % key)  
        mjd, kg, mjd_use, mjd_regression, kg_use, kg_smooth, kg_adjusted, kg_regression, outliers = _load_data(file_path)
        # original data
        if cb[key][0].isChecked():
            x = mjd
            y = kg
            xo = mjd_use
            yo = outliers
            ix = np.isnan(xo)
            xo = xo[~ix]
            yo = yo[~ix]
            color = col_le[key][0].text()
            if not color:
                color = tuple(np.random.random(3))
            linestyle = ls_le[key][0].currentText()
            if not linestyle:
                linestyle = '-'
            label = label_le[key][0].text()
            if not label:
                label = key.upper()
            x = [date_functions.mjd2datetimenumber(xj) for xj in x]
            if cmd_le[key][0].text():
                eval(cmd_le[key][0].text())
            else:
                ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)
                if not np.all(np.isnan(outliers)):
                    xo = [date_functions.mjd2datetimenumber(xj) for xj in xo]
                    ax.plot_date(xo,yo,color=color,linestyle='',marker='o',markerfacecolor='w',label=label+' (outliers)')
        # smoothed data
        if cb[key][1].isChecked():
            x = mjd_use
            y = kg_smooth
            color = col_le[key][1].text()
            if not color:
                color = tuple(np.random.random(3))
            linestyle = ls_le[key][1].currentText()
            if not linestyle:
                linestyle = '-'
            label = label_le[key][1].text()
            if not label:
                label = key.upper() + ' smoothed'
            x = [date_functions.mjd2datetimenumber(xj) for xj in x]
            if cmd_le[key][1].text():
                eval(cmd_le[key][1].text())
            else:
                ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label) 
        # regression model
        if cb[key][2].isChecked():
            x = mjd_regression
            y = kg_regression
            color = col_le[key][2].text()
            if not color:
                color = tuple(np.random.random(3))
            linestyle = ls_le[key][2].currentText()
            if not linestyle:
                linestyle = '-'
            label = label_le[key][2].text()
            if not label:
                label = key.upper() + ' regression'
            x = [date_functions.mjd2datetimenumber(xj) for xj in x]
            if cmd_le[key][2].text():
                eval(cmd_le[key][2].text())
            else:
                ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)  
                
def _plot_gia(ui,ax):
    keys = ('a','ice5','ice6','ice6_grace','klemann')
    cb = {'a':ui.cb_id061,
          'ice5':ui.cb_id062,
          'ice6':ui.cb_id063,
          'ice6_grace':ui.cb_id064,
          'klemann':ui.cb_id065}
    col_le = {'a':ui.le_id175,
          'ice5':ui.le_id178,
          'ice6':ui.le_id181,
          'ice6_grace':ui.le_id184,
          'klemann':ui.le_id187}
    ls_le = {'a':ui.cob_id025,
          'ice5':ui.cob_id026,
          'ice6':ui.cob_id027,
          'ice6_grace':ui.cob_id028,
          'klemann':ui.cob_id029}
    label_le = {'a':ui.le_id176,
          'ice5':ui.le_id179,
          'ice6':ui.le_id182,
          'ice6_grace':ui.le_id185,
          'klemann':ui.le_id188}
    cmd_le = {'a':ui.le_id177,
          'ice5':ui.le_id180,
          'ice6':ui.le_id183,
          'ice6_grace':ui.le_id186,
          'klemann':ui.le_id189}
    for key in keys:
        # load correction data
        file_path = os.path.join(project_dir,project_name,'corrections',key,'regression','%s.txt' % key)  
        mjd, kg, mjd_use, mjd_regression, kg_use, kg_smooth, kg_adjusted, kg_regression, outliers = _load_data(file_path)
        # original data
        if cb[key].isChecked():
            x = mjd
            y = kg
            color = col_le[key].text()
            if not color:
                color = tuple(np.random.random(3))
            linestyle = ls_le[key].currentText()
            if not linestyle:
                linestyle = '-'
            label = label_le[key].text()
            if not label:
                label = key.upper()
            x = [date_functions.mjd2datetimenumber(xj) for xj in x]
            if cmd_le[key].text():
                eval(cmd_le[key].text())
            else:
                ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)
       
def _plot_internal_validation(ui,ax):
    keys = ('m1','m1_syn')
    cb = {'m1':(ui.cb_id066,None,None),
          'm1_syn':(ui.cb_id067,ui.cb_id068,ui.cb_id069)}
    col_le = {'m1':(ui.le_id190,None,None),
              'm1_syn':(ui.le_id193,ui.le_id196,ui.le_id199)}
    ls_le = {'m1':(ui.cob_id030,None,None),
             'm1_syn':(ui.cob_id031,ui.cob_id032,ui.cob_id033)}
    label_le = {'m1':(ui.le_id191,None,None),
                'm1_syn':(ui.le_id194,ui.le_id197,ui.le_id200)}
    cmd_le = {'m1':(ui.le_id192,None,None),
              'm1_syn':(ui.le_id195,ui.le_id198,ui.le_id201)}
    for key in keys:
        # load validation data
        file_path = os.path.join(project_dir,project_name,'validation','internal','regression','%s.txt' % key)          
        mjd, kg, mjd_use, mjd_regression, kg_use, kg_smooth, kg_adjusted, kg_regression, outliers = _load_data(file_path)
        # original data
        if cb[key][0].isChecked():
            x = mjd
            y = kg
            xo = mjd_use
            yo = outliers
            ix = np.isnan(xo)
            xo = xo[~ix]
            yo = yo[~ix]
            color = col_le[key][0].text()
            if not color:
                color = tuple(np.random.random(3))
            linestyle = ls_le[key][0].currentText()
            if not linestyle:
                linestyle = '-'
            label = label_le[key][0].text()
            if not label:
                label = key.upper()
            x = [date_functions.mjd2datetimenumber(xj) for xj in x]
            if cmd_le[key][0].text():
                eval(cmd_le[key][0].text())
            else:
                ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)
                if not np.all(np.isnan(outliers)):
                    xo = [date_functions.mjd2datetimenumber(xj) for xj in xo]
                    ax.plot_date(xo,yo,color=color,linestyle='',marker='o',markerfacecolor='w',label=label+' (outliers)')
        # smoothed data
        if not cb[key][1]:
            continue
        if cb[key][1].isChecked():
            x = mjd_use
            y = kg_smooth
            color = col_le[key][1].text()
            if not color:
                color = tuple(np.random.random(3))
            linestyle = ls_le[key][1].currentText()
            if not linestyle:
                linestyle = '-'
            label = label_le[key][1].text()
            if not label:
                label = key.upper() + ' smoothed'
            x = [date_functions.mjd2datetimenumber(xj) for xj in x]
            if cmd_le[key][1].text():
                eval(cmd_le[key][1].text())
            else:
                ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label) 
        # regression model
        if cb[key][2].isChecked():
            x = mjd_regression
            y = kg_regression
            color = col_le[key][2].text()
            if not color:
                color = tuple(np.random.random(3))
            linestyle = ls_le[key][2].currentText()
            if not linestyle:
                linestyle = '-'
            label = label_le[key][2].text()
            if not label:
                label = key.upper() + ' regression'
            x = [date_functions.mjd2datetimenumber(xj) for xj in x]
            if cmd_le[key][2].text():
                eval(cmd_le[key][2].text())
            else:
                ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)  

def _plot_external_validation(ui,ax):
    keys = ('wgms','wgms_calibrated')
    cb = {'wgms':ui.cb_id070,
          'wgms_calibrated':ui.cb_id071}
    col_le = {'wgms':ui.le_id202,
          'wgms_calibrated':ui.le_id205}
    ls_le = {'wgms':ui.cob_id034,
          'wgms_calibrated':ui.cob_id035}
    label_le = {'wgms':ui.le_id203,
          'wgms_calibrated':ui.le_id206}
    cmd_le = {'wgms':ui.le_id204,
          'wgms_calibrated':ui.le_id207}
    for key in keys:
        # load external validation data
        file_path = os.path.join(project_dir,project_name,'validation','external','regression','%s.txt' % key)  
        mjd, kg, mjd_use, mjd_regression, kg_use, kg_smooth, kg_adjusted, kg_regression, outliers = _load_data(file_path)
        # original data
        if cb[key].isChecked():
            x = mjd
            y = kg
            xo = mjd_use
            yo = outliers
            ix = np.isnan(xo)
            xo = xo[~ix]
            yo = yo[~ix]
            color = col_le[key].text()
            if not color:
                color = tuple(np.random.random(3))
            linestyle = ls_le[key].currentText()
            if not linestyle:
                linestyle = '-'
            label = label_le[key].text()
            if not label:
                label = key.upper()
            x = [date_functions.mjd2datetimenumber(xj) for xj in x]
            if cmd_le[key].text():
                eval(cmd_le[key].text())
            else:
                ax.plot_date(x,y,color=color,linestyle=linestyle,marker=None,label=label)
                if not np.all(np.isnan(outliers)):
                    xo = [date_functions.mjd2datetimenumber(xj) for xj in xo]
                    ax.plot_date(xo,yo,color=color,linestyle='',marker='o',markerfacecolor='w',label=label+' (outliers)')


