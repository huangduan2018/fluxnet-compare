#!/usr/bin/env python

import sys, os, rc
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from cPickle import load as pickle_load
from cPickle import dump as pickle_dump

# This script reads messy FluxNet obs files, cleans and regroups the information

#===============================================================================
def main():
#===============================================================================

#-------------------------------------------------------------------------------
# ================================= USER INPUT =================================

# read the settings from the rc file
    rcdict     = rc.read('settings.rc')

# ==============================================================================
#-------------------------------------------------------------------------------
# extract the needed information from the rc file
    sites      = [s.strip(' ') for s in rcdict['sites'].split(',')]
    resolution = rcdict['resolution']  # can be hourly or daily

    # directory paths
    fluxnetdir = rcdict['obsdir']
    obsdir     = os.path.join(fluxnetdir, 'regrouped_data')
 
#-------------------------------------------------------------------------------
    if resolution == 'daily':
        filelist   = ['BE-Lon_FLUXNET2015_FULLSET_DD_2004-2014.csv',
                      'FR-Gri_FLUXNET2015_FULLSET_DD_2004-2013.csv',
                      'DE-Kli_FLUXNET2015_FULLSET_DD_2004-2014.csv',
                      'IT-BCi_mais_2004-2009_daily.csv']
    elif resolution == '3-hourly':
        filelist   = ['BE-Lon_FLUXNET2015_FULLSET_HH_2004-2014.csv',
                      'FR-Gri_FLUXNET2015_FULLSET_HH_2004-2013.csv',
                      'DE-Kli_FLUXNET2015_FULLSET_HH_2004-2014.csv',
                      'IT-BCi_mais_2004-2009_daily.csv']

#-------------------------------------------------------------------------------
# Extract timeseries for the different sites

    # read files for the diferent sites
    f = open_csv(obsdir,filelist,convert_to_float=True)

    series = dict()
    filepath = os.path.join(fluxnetdir,'%s_timeseries_OBS.pickle'%resolution)

    for fnam,site in zip(filelist, sites):

        print site
       
        # TA_F_DAY: average daytime Ta_day from meas and ERA (*C)
        # SW_IN_F: SWin from meas and ERA (W.m-2) 
        # VPD_F: VPD consolidated from VPD_F_MDS and VPD_F_ERA (hPa)
        # TS_F_MDS_1 to 4: Tsoil of 4 soil layers (*C)
        # SWC_F_MDS_1 to 4: soil water content (%) of 4 layers (1=shallow) 
        # NT = night-time partitioning method (gC m-2 s-1)
        # VUT: variable ref u* between years
        FLUX_variables = ['TA_F_DAY', 'SW_IN_F', 'VPD_F', 'TS_F_MDS_1', 
                          'TS_F_MDS_2', 'TS_F_MDS_3', 'SWC_F_MDS_1', 'SWC_F_MDS_2',
                          'SWC_F_MDS_3', 'GPP_NT_VUT_REF', 'RECO_NT_VUT_REF',
                          'NEE_VUT_REF', 'crop', 'LAI', 'AGB', 'C_height']
        FLUX_varnames  = ['Ta_day', 'SWin', 'VPD', 'Ts_1', 'Ts_2', 'Ts_3', 'SWC_1',
                          'SWC_2', 'SWC_3', 'GPP', 'TER', 'NEE', 'crop_no', 'LAI',
                          'AGB', 'CHT']
        IT_variables = ['SWC_avg', 'GPP', 'Reco', 'NEE', 'crop', 'GLAI', 'AGB', 
                        'C_height']
        IT_varnames  = ['SWC', 'GPP', 'TER', 'NEE', 'crop_no', 'LAI', 'AGB', 
                        'CHT']

        # timestamps for all daily timeseries
        startyear = str(f[fnam]['TIMESTAMP'][0])[0:4]
        endyear   = str(f[fnam]['TIMESTAMP'][-1])[0:4]
        startdate = '%s-01-01 00:00:00'%startyear
        enddate   = '%s-12-31 23:30:00'%endyear
        if site=='DE-Kli': enddate = '%s-12-31 23:00:00'%endyear

        series[site] = dict()

        if resolution == '3-hourly':
            tm = pd.date_range(startdate, enddate, freq='30min')
            if (site!='IT-BCi'):
                for var,varname in zip(FLUX_variables[:12], FLUX_varnames[:12]):
                    # if the fluxes are half-hourly, I convert them to 3-hourly
                    if varname == 'Ta_day':
                        series[site]['Ta'] = pd.Series(f[fnam]['TA_F'], index=tm)
                    elif ((varname == 'SWC_2' or varname == 'SWC_3') and 
                    site == 'FR-Gri'):
                        series[site][varname] = pd.Series([-9999.]*len(tm), index=tm)
                    else: 
                        series[site][varname] = pd.Series(f[fnam][var], index=tm)
                    print varname

        elif resolution == 'daily':
            tm = pd.date_range(startdate, enddate, freq='1d')
            if (site!='IT-BCi'):
                for var,varname in zip(FLUX_variables, FLUX_varnames):
                    series[site][varname] = pd.Series(f[fnam][var], index=tm)
                    print varname
            else:
                tm_irreg = [pd.to_datetime('%s-%s-%s'%(str(t)[0:4],str(t)[4:6],
                                    str(t)[6:8])) for t in f[fnam]['TIMESTAMP']]
                # since the time records has gaps in the IT-BCi data, we use a 
                # special function to fill the gaps with -9999. values and
                # convert it to pandas timeseries
                for var,varname in zip(IT_variables, IT_varnames):
                    #if varname == 'VPD':
                    #    ta     = f[fnam]['T_avg']
                    #    dayvar = f[fnam]['Rh_avg'] / 100. * 6.11 * np.exp(ta /\
                    #             (238.3 + ta) * 17.2694)
                    dayvar = f[fnam][var]
                    series[site][varname] = convert2pandas(tm_irreg, dayvar, tm)
                    print varname
        else:
            print "Wrong CO2 fluxes temporal resolution: must be either "+\
                  "'daily' or '3-hourly'"
            sys.exit()



    # we store the pandas series in one pickle file
    pickle_dump(series, open(filepath,'wb'))

#-------------------------------------------------------------------------------
# plot timeseries

    # Let's plot the available micromet variables that are important for WOFOST
    #plot_fluxnet_micromet(obsdir,sites,[2005,2005],'-')

    # Let's plot GPP, TER, NEE
    #plot_fluxnet_daily_c_fluxes(obsdir,sites,[2004,2014],'-')

    #plot_fluxnet_LAI_CHT_AGB(obsdir,sites,[2004,2014],'o')
#-------------------------------------------------------------------------------

    return series

#===============================================================================
def convert_halfhourly_to_threehourly(halfhourly_var):
#===============================================================================

    return None

#===============================================================================
def convert2pandas(irreg_tm, VAR, tm):
#===============================================================================
    """
    Function to make a regular pandas timeseries out of an irregular one

    Function arguments:
    ------------------
    irreg_tm: the irregular time axis (should be in the format of a datetime
              list, and there are several ways to do that)
    VAR:      the y axis that correspond to the irregular time axis(format is 
              a list of floats)
    tm :      the pandas timeframe on which to convert the irregular timeseries

    example of use:
    --------------
    irreg_tm = [pd.to_datetime('2004-01-01'), pd.to_datetime('2004-01-05'), 
                pd.to_datetime('2004-01-31'), pd.to_datetime('2004-02-03')]
    LAI      = [0.,1.,2.,0.5]
    tm       = pd.date_range('2004-01-01 00:00:00', '2004-12-31 23:59:59', freq='1d')

    pandas_LAI_series = convert2pandas(irreg_tm, LAI, tm)
    print pandas_LAI_series

    returns:
    -------
    2004-01-01       0
    2004-01-02   -9999
    2004-01-03   -9999
    2004-01-04   -9999
    2004-01-05       1
                  ...
    2004-01-31       2
                  ...
    2004-02-03     0.5
                  ...
    2004-12-28   -9999
    2004-12-29   -9999
    2004-12-30   -9999
    2004-12-31   -9999
    Freq: D, dtype: float64

    """

    # initialize converted VAR time series
    reg_var = []
    # we go through the pandas timeframe of the carbon fluxes
    for t in tm:
        if t not in irreg_tm:
            reg_var += [-9999.]
        # if the date matches a date reported in the irregular data
        else:
            ind_var = np.argmin(np.abs([i-t for i in irreg_tm]))
            reg_var += [VAR[ind_var]]
    # we store as a pandas timeseries object
    pdseries = pd.Series(reg_var, index=tm)

    return pdseries

#===============================================================================
def plot_fluxnet_LAI_CHT_AGB(directory, sites_list, timeframe, linestyle):
#===============================================================================

    """
    This function plots the LAI, CHT and AGB auxilliary data from FluxNet sites
 
    Function arguments:
    ------------------
    directory:       path to the directory where the daily variables files
                     ('FLX_sitename_FLUXNET2015_FULLSET_DD_2004-2014_1-1')
                     are stored
    sites_list: list of fluxnet sites names
    time frame:      list of two years (e.g. [2005, 2007]) gives the time frame 
                     to plot
    linestyle:       line style of the plot to make (can be '-', '--', '+', etc)

    """

    plt.close('all')

    print "\nPlotting timeseries of LAI, CHT and AGB\n"

    for site in sites_list:
        print site+'\n'

        filepath = os.path.join(directory,'timeseries_%s.pickle'%site)
        series = pickle_load(open(filepath,'rb'))
       
        # settings of the plot:
        figs, axes = plt.subplots(nrows=4, ncols=1, figsize=(8,10))
        figs.subplots_adjust(0.1,0.07,0.98,0.95,0.,0.)
        variables = ['crop_no','LAI','CHT','AGB']
        axlabels  = ['crop ID',r'LAI (m$^{-2}$ m$^{-2}$)','CHT (m)',r'AGB (kg$_{DM}$ m$^{-2}$)']
        ylims     = [(0.,14.),(0.,7.5),(0.,3.),(0.,3.)]
        start = str(int(timeframe[0]))
        end   = str(int(timeframe[1]))
        lab = ''
        ls = linestyle
        fsz = 14 # fonsize of x and y axis ticks

        for ax, var, axlabel, ylim in zip(axes,variables,axlabels,ylims):
            if var=='crop_no': 
                ls='-'
            else: 
                ls=linestyle
            series[site][var][start:end].plot(ax=ax, lw=2, style=ls,
            label=lab, fontsize=fsz)
            ax.axhline(y=0., c='k')
            ax.set_ylim(ylim)
            ax.set_ylabel(axlabel)
            if var != 'AGB': ax.get_xaxis().set_visible(False)
        figs.suptitle(site, fontsize=14)
    plt.show()

    return None

#===============================================================================
def plot_fluxnet_daily_c_fluxes(directory, sites_list, timeframe, linestyle):
#===============================================================================

    """
    This function plots the average daytime Tair, average SWin, average VPD, 
    average Tsoil and average volumetric soil moisture content, all taken from
    the formatted 2015 FluxNet DD (daily data) files.
 
    NB: Daily averages have been computed from the half-hourly data.

    Function arguments:
    ------------------
    directory:       path to the directory where the daily variables files
                     ('FLX_sitename_FLUXNET2015_FULLSET_DD_2004-2014_1-1')
                     are stored
    sites_list: list of fluxnet sites names
    time frame:      list of two years (e.g. [2005, 2007]) gives the time frame 
                     to plot
    linestyle:       line style of the plot to make (can be '-', '--', '+', etc)

    """

    plt.close('all')

    print "\nPlotting carbon fluxes timeseries\n"

    for site in sites_list:
        print site+'\n'

        filepath = os.path.join(directory,'timeseries_%s.pickle'%site)
        series = pickle_load(open(filepath,'rb'))
       
        # settings of the plot:
        figs, axes = plt.subplots(nrows=4, ncols=1, figsize=(8,10))
        figs.subplots_adjust(0.1,0.07,0.98,0.95,0.,0.)
        variables = ['crop_no','GPP','TER','NEE']
        axlabels  = ['crop ID',r'GPP (g m$^{-2}$ d$^{-1}$)',r'TER (g m$^{-2}$ d$^{-1}$)',r'NEE (g m$^{-2}$ d$^{-1}$)']
        ylims     = [(0.,14.),(-30.,5.),(-5.,20.),(-30.,30.)]
        start = str(int(timeframe[0]))
        end   = str(int(timeframe[1]))
        lab = ''
        ls = linestyle
        fsz = 14 # fonsize of x and y axis ticks

        for ax, var, axlabel, ylim in zip(axes,variables,axlabels,ylims):
            if (var=='TER' or var=='GPP'): lab='-'
            series[site][var][start:end].plot(ax=ax, lw=2, style=ls,
            label=lab, fontsize=fsz)
            if var=='GPP':
                series[site]['GPP_day'][start:end].plot(ax=ax, lw=2, style=ls, 
                label='day partition', fontsize=fsz)
                series[site]['GPP_night'][start:end].plot(ax=ax, lw=2, style=ls, 
                label='night partition', fontsize=fsz)
            if var=='TER':
                series[site]['TER_day'][start:end].plot(ax=ax, lw=2, style=ls, 
                label='day partition', fontsize=fsz)
                series[site]['TER_night'][start:end].plot(ax=ax, lw=2, style=ls, 
                label='night partition', fontsize=fsz)
            ax.axhline(y=0., c='k')
            ax.set_ylim(ylim)
            if (var=='GPP_day'): ax.legend(loc='lower left',prop={'size':12})
            if (var=='TER_day'): ax.legend(loc='upper left',prop={'size':12})
            ax.set_ylabel(axlabel)
            if var != 'NEE': ax.get_xaxis().set_visible(False)
        figs.suptitle(site, fontsize=14)
    plt.show()

    return None

#===============================================================================
def plot_fluxnet_micromet(directory, sites_list, timeframe, linestyle):
#===============================================================================

    """
    This function plots the average daytime Tair, average SWin, average VPD, 
    average Tsoil and average volumetric soil moisture content, all taken from
    the formatted 2015 FluxNet DD (daily data) files.
 
    NB: Daily averages have been computed from the half-hourly data.

    Function arguments:
    ------------------
    directory:       path to the directory where the daily variables files
                     ('FLX_sitename_FLUXNET2015_FULLSET_DD_2004-2014_1-1')
                     are stored
    sites_list: list of fluxnet sites names
    time frame:      list of two years (e.g. [2005, 2007]) gives the time frame 
                     to plot
    linestyle:       line style of the plot to make (can be '-', '--', '+', etc)

    """

    plt.close('all')

    for site in sites_list:
        print site+'\n'

        filepath = os.path.join(directory,'timeseries_%s.pickle'%site)
        series = pickle_load(open(filepath,'rb'))
       
        # settings of the plot:
        figs, axes = plt.subplots(nrows=5, ncols=1, figsize=(8,9))
        figs.subplots_adjust(0.1,0.07,0.98,0.95,0.,0.)
        variables = ['Ta_day','SWin','VPD','Ts_1','SWC_1']
        axlabels  = [r'T$_{air}$ day (deg C)', r'SW$_{in}$ (W m$^{-2}$)',
                     'VPD (hPa)',r'T$_{soil}$ (deg C)','SWC (%)']
        ylims     = [(-20.,40.),(0.,400.),(0.,17.),(-5.,35.),(0.,60.)]
        start = str(int(timeframe[0]))
        end   = str(int(timeframe[1]))
        ls = linestyle
        lab = 'bullshit'
        fsz = 14 # fonsize of x and y axis ticks

        for ax, var, axlabel, ylim in zip(axes,variables,axlabels,ylims):
            if var.endswith('_1'): lab='layer 1'
            series[site][var][start:end].plot(ax=ax, lw=2, style=ls,
            label=lab, fontsize=fsz)
            if var=='Ts_1':
                series[site]['Ts_2'][start:end].plot(ax=ax, lw=2, style=ls, 
                label='layer 2', fontsize=fsz)
                series[site]['Ts_3'][start:end].plot(ax=ax, lw=2, style=ls, 
                label='layer 3', fontsize=fsz)
            if var=='SWC_1':
                series[site]['SWC_2'][start:end].plot(ax=ax, lw=2, style=ls, 
                label='layer 2', fontsize=fsz)
                series[site]['SWC_3'][start:end].plot(ax=ax, lw=2, style=ls, 
                label='layer 3', fontsize=fsz)
            if (var=='Ts_1'): ax.legend(loc='upper left',prop={'size':12})
            if (var=='SWC_1'): ax.legend(loc='lower left',prop={'size':12})
            ax.axhline(y=0., c='k')
            ax.set_ylim(ylim)
            ax.set_ylabel(axlabel)
            if var=='Ta_day': ax.set_yticks([-10.,0.,10.,20.,30.,40.])
            if var=='SWin': ax.set_yticks([0.,100.,200.,300.,400.])
            if var=='VPD': ax.set_yticks([0.,5.,10.,15.])
            if var=='Ts_1': ax.set_yticks([0.,10.,20.,30.])
            if var != 'SWC_1': ax.get_xaxis().set_visible(False)
        figs.suptitle(site, fontsize=14)
    plt.show()

    return None

#===============================================================================
# Function to open normal csv files
def open_csv(inpath,filelist,convert_to_float=False):
#===============================================================================

    from csv import reader as csv_reader

    Dict = {}

    for i,namefile in enumerate(filelist):
         
        #print "\nOpening %s......"%(namefile)

        # open file, read all lines
        inputpath = os.path.join(inpath,namefile)
        f=open(inputpath,'rU') 
        reader=csv_reader(f, delimiter=',', skipinitialspace=True)
        lines=[]
        for row in reader:
            lines.append(row)
        f.close()

        # storing headers in list headerow
        headerow=lines[0]

        # deleting rows that are not data (first and last rows of the file)
        del lines[0]

        # transforming data from string to float type
        converted_data=[]
        for line in lines:
            if '' in line:
                newline = []
                for it in line:
                    if it=='': newline += ['-9999.']
                    if it!='': newline += [it]
                line = newline
            converted_data.append(map(float,line))
        data = np.array(converted_data)

        # creating one dictionnary and storing the float data in it
        dictnamelist= {}
        for j,varname in enumerate(headerow):
            dictnamelist[varname]=data[:,j]
        Dict[namefile] = dictnamelist
    
        #print "Dictionary created!"

    return Dict
#===============================================================================
if __name__=='__main__':
    main()
#===============================================================================
