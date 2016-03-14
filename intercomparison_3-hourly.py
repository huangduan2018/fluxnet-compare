#!/usr/bin/env python

import sys, os, rc, math
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from cPickle import load as pickle_load
from scipy.stats.stats import linregress as linreg

#===============================================================================
# This script does some standard analysis on FluxNet sites
def main():
#===============================================================================
    global wofostdir, sibcasadir, obsdir
#-------------------------------------------------------------------------------
# ================================= USER INPUT =================================

# read the settings from the rc file (mostly directory paths)
    rcdict    = rc.read('settings.rc')
    sites      = [s.strip(' ') for s in rcdict['sites'].split(',')]
    years      = [int(s.strip(' ')) for s in rcdict['years'].split(',')]
    TER_method = rcdict['TER_method']
    R10        = rcdict['R10']
    resolution = rcdict['resolution']

# specific plotting options:
    #TER_method = 'grow-only' # this is to select the corresponding WOFOST output file
    #R10        = '0.08' # this is to select the corresponding WOFOST output file


#===============================================================================
#-------------------------------------------------------------------------------
# extract the needed information

    # input data directory paths
    rootdir     = rcdict['rootdir']
    sibcasadir  = os.path.join(rootdir,'intercomparison_study/SiBCASA_runs')
    wofostdir   = rcdict['outputdir'] 
    obsdir      = rcdict['obsdir']
    figdir      = os.path.join(rootdir,'intercomparison_study/figures')

#-------------------------------------------------------------------------------
# Start a directory to store OBS, SIMW (wofost), SIMB (SiBCASA)

    # recover the FluxNet observed data from pickle files
    res_timeseries = dict()
    res_timeseries['OBS'] = dict()
    res_timeseries['SIMB'] = dict()
    res_timeseries['SIMW'] = dict()

    filename = os.path.join(obsdir, '%s_timeseries_OBS.pickle'%resolution)
    try:
        res_timeseries['OBS'] = pickle_load(open(filename,'rb'))
    except IOError:
        print 'could not find the observations output file %s'%filename
        res_timeseries['OBS'] = None

    # recover the SiBCASA runs
    filename = os.path.join(sibcasadir, '%s_timeseries_SiBCASA.pickle'%resolution)
    try:
        res_timeseries['SIMB'] = pickle_load(open(filename,'rb'))
    except IOError:
        print 'could not find the SiBCASA output file %s'%filename
        res_timeseries['SIMB'] = None

    # recover the WOFOST runs
    filename = os.path.join(wofostdir, '%s_timeseries_'%resolution +\
               '%s_R10=%s_WOFOST_crop_rotation.pickle'%(TER_method,R10))
    try:
        res_timeseries['SIMC'] = pickle_load(open(filename,'rb'))
    except IOError:
        print 'could not find the WOFOST output file %s'%filename
        res_timeseries['SIMC'] = None

#-------------------------------------------------------------------------------
# plot the observed and simulated timeseries with pandas library
# with pandas we plot all years one after another, and can zoom in on one 
# particular year

    plt.close('all')

    # create figure sub-folder if it doesn't already exists
    figsubdir = os.path.join(figdir,'R10=%s/TER_%s/'%(R10,TER_method)+\
                '3-hourly_fluxes_perf')
    if not os.path.exists(figsubdir):
        print 'creating new directory %s'%figsubdir
        os.makedirs(figsubdir)

#-------------------------------------------------------------------------------
# we plot the 3-hourly fluxes of simulations versus observations for the years 
# and sites that perform well on the daily scale

    # we plot years 2005, 2009, 2013 for site BE-Lon, which showed an extremely
    # good result on the SIM vs OBS daily fluxes comparison

    years = [2005,2009,2013]

    variables = ['GPP','TER','NEE']
    axlabels  = [r'GPP (g m$^{-2}$ d$^{-1}$)',
                  r'TER (g m$^{-2}$ d$^{-1}$)',r'NEE (g m$^{-2}$ d$^{-1}$)']
    ylims     = [(-60.,5.),(0.,20.),(-50.,15.)]
    one_to_one = np.arange(-100,100,10)

    for site in ['BE-Lon']:
        if site != 'IT-BCi':
            for year in years:
                timeframe = [year,year]
                start = str(int(timeframe[0]))+'-05-01'
                end   = str(int(timeframe[1]))+'-07-01'
                print site
                for var, axlabel, lim in zip(variables,axlabels,ylims):
                    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5))
                    fig.subplots_adjust(0.15,0.15,0.85,0.85,0.,0.)
                    # select every 6 half-hourly flux in the observations,
                    # to get to 3-hourly frequency
                    OBS = res_timeseries['OBS'][site][var][::6][start:end].dropna()
                    # convert the 3-hourly simulated fluxes from gC.m-2.s-1
                    # to micromol CO2 .m-2.s-1
                    SIM = res_timeseries['SIMC'][site][var][start:end].dropna()
                    SIM = SIM * 1000000. / 12.01 #micro mol CO2 per m2 per sec
                    # the observed GPP needs a minus sign for convention
                    if var=='GPP': OBS=-OBS
                    # use the min and max to frame the figure
                    print var, min(min(OBS), min(SIM)), max(max(OBS), max(SIM))
                    varmin = math.floor(min(min(OBS), min(SIM)))
                    varmax = math.ceil(max(max(OBS), max(SIM)))
                    ax.scatter(OBS, SIM, marker='o')

                    # fit a linear regression line in OBS/SIM scatter plot
                    # and plot line and R2
                    mask = ~np.isnan(SIM)
                    z = np.polyfit(OBS[mask], SIM[mask], 1)
                    p = np.poly1d(z)
       	            ax.plot(one_to_one,p(one_to_one),'r-')
                    slope, intercept, r_value, p_value, std_err = \
                                                    linreg(OBS[mask], SIM[mask])
                    ax.annotate(r'r$^2$ = %.2f'%r_value**2, xy=(0.95, 0.15), 
                       xytext=(0.15, 0.9), xycoords='axes fraction',
                       ha='center', va='center', color='r')

                    ax.plot(one_to_one,one_to_one, c='k', lw=1)
                    ax.set_xlabel('obs')
                    ax.set_ylabel('sim')
                    ax.set_xlim(varmin,varmax)
                    ax.set_ylim(varmin,varmax)
                    fig.suptitle(r'%s 3-hourly %s fluxes ($\mu$'%(site, var)+
                     r'mol m$^{-2}$ s$^{-1}$)'+'\nfrom %s to %s\n'%(start,end))
                    fig.savefig(os.path.join(figsubdir,
                                       '%s_%s_%s.png'%(site,year,var)), dpi=300)
                #plt.show()

#===============================================================================
if __name__=='__main__':
    main()
#===============================================================================
