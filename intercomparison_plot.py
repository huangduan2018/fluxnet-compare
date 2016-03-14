#!/usr/bin/env python

import sys, os, rc, math
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from cPickle import load as pickle_load

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

    filename = os.path.join(obsdir, 'timeseries_OBS.pickle')
    try:
        res_timeseries['OBS'] = pickle_load(open(filename,'rb'))
    except IOError:
        print 'could not find the observations output file %s'%filename
        res_timeseries['OBS'] = None

    # recover the SiBCASA runs
    filename = os.path.join(sibcasadir, 'timeseries_SiBCASA.pickle')
    try:
        res_timeseries['SIMB'] = pickle_load(open(filename,'rb'))
    except IOError:
        print 'could not find the SiBCASA output file %s'%filename
        res_timeseries['SIMB'] = None

    # recover the WOFOST runs
    filename = os.path.join(wofostdir, 'timeseries_%s_R10=%s_'%(TER_method,R10) +\
               'WOFOST_crop_rotation.pickle')
    try:
        print 'opening the WOFOST output file %s'%filename
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
    figsubdir = os.path.join(figdir,'R10=%s/TER_%s/'%(R10,TER_method))
    if not os.path.exists(figsubdir):
        print 'creating new directory %s'%figsubdir
        os.makedirs(figsubdir)

#-------------------------------------------------------------------------------

    years = np.arange(2004,2015,1)

    for site in sites:
        for year in years:
            timeframe = [year,year]
            print site
            figs, axes = plt.subplots(nrows=4, ncols=1, figsize=(8,10))
            figs.subplots_adjust(0.1,0.07,0.98,0.95,0.,0.)
            variables = ['crop_no','GPP','TER','NEE']
            axlabels  = ['crop ID',r'GPP (g m$^{-2}$ d$^{-1}$)',
                          r'TER (g m$^{-2}$ d$^{-1}$)',r'NEE (g m$^{-2}$ d$^{-1}$)']
            ylims     = [(0.,14.),(-18.,2.),(-1.,12.),(-10.,10.)]
            start = str(int(timeframe[0]))
            end   = str(int(timeframe[1]))
            print '[%s:%s]'%(start,end)
            fsz = 14 # fonsize of x and y axis ticks
         
            for ax, var, axlabel, ylim in zip(axes,variables,axlabels,ylims):
                if (var=='crop_no'): 
                    try:
                        OBS = res_timeseries['OBS'][site][var][start:end].dropna()
                        OBS[~(OBS==-9999.)].plot(ax=ax, lw=2, 
                        #OBS.plot(ax=ax, lw=2, 
                        style='-', label='obs', fontsize=fsz)
                        crop_no = OBS[0]
                        minobs = OBS[~(OBS==-9999.)].min()
                        maxobs = OBS[~(OBS==-9999.)].max()
                    except TypeError:
                        minobs = 0.
                        maxobs = 0.
                    minwof = 1.
                    maxwof = 1.
                elif (var=='TER'): 
                    # observations
                    try:
                        OBS = res_timeseries['OBS'][site][var][start:end].dropna()
                        OBS[~(OBS==-9999.)].plot(ax=ax, lw=2, 
                        #OBS.plot(ax=ax, lw=2, c='b',
                        style='+', label='obs', fontsize=fsz)
                        minobs = OBS[~(OBS==-9999.)].min()
                        maxobs = OBS[~(OBS==-9999.)].max()
                    except TypeError:
                        minobs = 0.
                        maxobs = 0.
                    # SIBCASA sims
                    try:
                        #res_timeseries['SIMB'][site]['Raut'][start:end].plot(ax=ax, 
                        #lw=2, c='g', style=':', label='SiBCASA Raut', fontsize=fsz)
                        res_timeseries['SIMB'][site][var][start:end].plot(ax=ax, 
                        lw=2, c='g', style='--', label='SiBCASA TER', fontsize=fsz)
                    except TypeError:
                        pass
                    # WOFOST sims
                    try:
                        #WOF = res_timeseries['SIMC'][site]['Raut'][start:end].dropna()
                        #WOF.plot(ax=ax, lw=2, c='r',
                        #style='_', label='WOFOST Raut', fontsize=fsz)
                        WOF = res_timeseries['SIMC'][site][var][start:end].dropna()
                        WOF.plot(ax=ax, lw=2, c='r',
                        style='x', label='WOFOST TER', fontsize=fsz)
                        minwof = WOF.min()
                        maxwof = WOF.max()
                    except TypeError:
                        minwof = 0.
                        maxwof = 0.
                        WOF = 0.
                else:
                    # observations
                    try:
                        OBS = res_timeseries['OBS'][site][var][start:end].dropna()
                        OBS[~(OBS==-9999.)].plot(ax=ax, lw=2, 
                        #OBS.plot(ax=ax, lw=2, c='b',
                        style='+', label='obs', fontsize=fsz)
                        minobs = OBS[~(OBS==-9999.)].min()
                        maxobs = OBS[~(OBS==-9999.)].max()
                    except TypeError:
                        minobs = 0.
                        maxobs = 0.
                    # SIBCASA sims
                    try:
                        res_timeseries['SIMB'][site][var][start:end].plot(ax=ax, lw=2, c='g',
                        style='--', label='SiBCASA', fontsize=fsz)
                    except TypeError:
                        pass
                    # WOFOST simsA
                    try:
                        WOF = res_timeseries['SIMC'][site][var][start:end].dropna()
                        #WOF[~(WOF==-9999.)].plot(ax=ax, lw=2, 
                        WOF.plot(ax=ax, lw=2, c='r',
                        style='x', label='WOFOST', fontsize=fsz)
                        minwof = WOF.min()
                        maxwof = WOF.max()
                    except TypeError:
                        minwof = 0.
                        maxwof = 0.
                        WOF = 0.
                ax.axhline(y=0., c='k')
                minvar = math.floor(min(minobs,minwof))-1.
                maxvar = math.ceil(max(maxobs,maxwof))+1.
                ax.set_ylim(minvar,maxvar)
                #ax.set_ylim(ylim)
                if (var=='GPP'): ax.legend(loc='lower left',prop={'size':12})
                #if (var=='TER'): ax.legend(loc='upper left',prop={'size':10})
                ax.set_ylabel(axlabel)
                if var != 'NEE': ax.get_xaxis().set_visible(False)
            figs.suptitle(site, fontsize=14)
            figs.savefig(os.path.join(figsubdir,'crop%i_%s_%i.png'%(crop_no,site,timeframe[0])))
            plt.close('all')
    #plt.show()

#-------------------------------------------------------------------------------

    timeframe = [2004,2014]
    for site in sites:

        print site
        figs, axes = plt.subplots(nrows=4, ncols=1, figsize=(15,10))
        figs.subplots_adjust(0.1,0.07,0.98,0.95,0.,0.)
        variables = ['crop_no','GPP','TER','NEE']
        axlabels  = ['crop ID',r'GPP (g m$^{-2}$ d$^{-1}$)',
                      r'TER (g m$^{-2}$ d$^{-1}$)',r'NEE (g m$^{-2}$ d$^{-1}$)']
        ylims     = [(0.,14.),(-30.,2.),(-2.,20.),(-20.,10.)]
        start = str(int(timeframe[0]))
        end   = str(int(timeframe[1]))
        print '[%s:%s]'%(start,end)
        fsz = 14 # fonsize of x and y axis ticks
        
        for ax, var, axlabel, ylim in zip(axes,variables,axlabels,ylims):
            if (var=='crop_no'): 
                try:
                    OBS = res_timeseries['OBS'][site][var][start:end].dropna()
                    OBS[~(OBS==-9999.)].plot(ax=ax, lw=2, 
                    #OBS.plot(ax=ax, lw=2, 
                    style='-', label='obs', fontsize=fsz)
                    crop_no = OBS[0]
                    minobs = OBS[~(OBS==-9999.)].min()
                    maxobs = OBS[~(OBS==-9999.)].max()
                except TypeError:
                    minobs = 0.
                    maxobs = 0.
                minwof = 1.
                maxwof = 1.
            else:
                # observations
                try:
                    OBS = res_timeseries['OBS'][site][var][start:end].dropna()
                    OBS[~(OBS==-9999.)].plot(ax=ax, lw=2, 
                    #OBS.plot(ax=ax, lw=2, c='b',
                    style='+', label='obs', fontsize=fsz)
                    minobs = OBS[~(OBS==-9999.)].min()
                    maxobs = OBS[~(OBS==-9999.)].max()
                except TypeError:
                    minobs = 0.
                    maxobs = 0.
                # SIBCASA sims
                try:
                    res_timeseries['SIMB'][site][var][start:end].plot(ax=ax, lw=2, c='g',
                    style='--', label='SiBCASA', fontsize=fsz)
                except TypeError:
                    pass
                # WOFOST simsA
                try:
                    WOF = res_timeseries['SIMC'][site][var][start:end].dropna()
                    #WOF[~(WOF==-9999.)].plot(ax=ax, lw=2, 
                    WOF.plot(ax=ax, lw=2, c='r',
                    style='x', label='WOFOST', fontsize=fsz)
                    minwof = WOF.min()
                    maxwof = WOF.max()
                except TypeError:
                    minwof = 0.
                    maxwof = 0.
                    WOF = 0.
            ax.axhline(y=0., c='k')
            minvar = math.floor(min(minobs,minwof))-1.
            maxvar = math.ceil(max(maxobs,maxwof))+1.
            #ax.set_ylim(minvar,maxvar)
            ax.set_ylim(ylim)
            if (var=='GPP'): ax.legend(loc='lower left',prop={'size':12})
            ax.set_ylabel(axlabel)
            if var != 'NEE': ax.get_xaxis().set_visible(False)
        figs.suptitle(site, fontsize=14)
        figs.savefig(os.path.join(figsubdir,'timeseries_crop%i_%s_%i-%i.png'%(crop_no,site,timeframe[0],timeframe[1])))

    plt.close('all')

#===============================================================================
if __name__=='__main__':
    main()
#===============================================================================
