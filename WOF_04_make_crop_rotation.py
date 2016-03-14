#!/usr/bin/env python

import sys, os, rc
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from cPickle import load as pickle_load
from cPickle import dump as pickle_dump

#===============================================================================
# This script does some standard analysis on FluxNet sites
def main():
#===============================================================================
    global outputdir, obsdir
#-------------------------------------------------------------------------------
# ================================= USER INPUT =================================

# read the settings from the rc file
    rcdict    = rc.read('settings.rc')

#===============================================================================
#-------------------------------------------------------------------------------
# extract the needed information for that script
    sites      = [s.strip(' ') for s in rcdict['sites'].split(',')]
    years      = [s.strip(' ') for s in rcdict['years'].split(',')]
    TER_method = rcdict['TER_method']
    R10        = rcdict['R10']
    resolution  = rcdict['resolution']  # can be hourly or daily
    if resolution=='daily': res='1d'
    elif resolution=='3-hourly': res='3H'

    # directory paths
    outputdir  = rcdict['outputdir']
    obsdir     = rcdict['obsdir']
    forwardir  = os.path.join(outputdir, 'forward_runs')

#-------------------------------------------------------------------------------
# load the WOFOST runs of all crops

    # we store the two pandas series in one pickle file
    filepath = os.path.join(forwardir,'%s_timeseries_'%resolution+\
                            '%s_WOFOST.pickle'%TER_method)
    series   = pickle_load(open(filepath,'rb'))

    filepath = os.path.join(obsdir,'daily_timeseries_OBS.pickle')
    obs      = pickle_load(open(filepath,'rb'))

    final_series = dict()

    for s,site in enumerate(sites):
        print site
        print obs[site].keys()
        final_series[site] = dict()

        # read the crop rotation from FluxNet file
        rotation = obs[site]['crop_no']

        # slice each year's required time series, append to final series
        for varname in ['GPP','TER','Raut','Rhet','NEE']:
            print 'variable %s'%varname
            var = []
            for year in years:
                
                # get the crop number for that year
                if site != 'IT-BCi':
                    try:
                        crop_no = rotation[year:year][0]
                    except IndexError: # index error occurs when the year is
                                       # not in the rotation time series
                        startdate = '%s-01-01 00:00:00'%year
                        enddate   = '%s-12-31 23:59:59'%year
                        dtimes    = pd.date_range(startdate, enddate, freq=res)
                        na_vals   = np.array(len(dtimes)*[np.nan])
                        var      += [pd.Series(na_vals, index=dtimes)]
                        print '   ',site, year, 'unknown crop cover: skip.'
                        continue
                elif site == 'IT-BCi':
                    if int(year) not in np.arange(2004,2010,1): 
                        startdate = '%s-01-01 00:00:00'%year
                        enddate   = '%s-12-31 23:59:59'%year
                        dtimes    = pd.date_range(startdate, enddate, freq=res)
                        na_vals   = np.array(len(dtimes)*[np.nan])
                        var      += [pd.Series(na_vals, index=dtimes)]
                        print '   ',site, year, 'unknown crop cover: skip.'
                        continue
                    else:
                        crop_no = 2

                # try slicing and concatenating that year's timeseries from file
                try:
                    # if the GPP = 0 (failed growing season), we set TER and 
                    # NEE to zero as well
                    if np.mean(series[site]['c%i'%crop_no]['GPP'][year:year]) == 0.:
                        startdate = '%s-01-01 00:00:00'%year
                        enddate   = '%s-12-31 23:59:59'%year
                        dtimes    = pd.date_range(startdate, enddate, freq=res)
                        zeros     = np.array(len(dtimes)*[0.])
                        var      += [pd.Series(zeros, index=dtimes)]
                    else:
                        var += [series[site]['c%i'%crop_no][varname][year:year]]
                    print '   ',site, year, '%2i'%crop_no, 'slicing'
                except KeyError: # key error occurs when we haven't ran a crop
                                 # or a year with WOFOST
                    startdate = '%s-01-01 00:00:00'%year
                    enddate   = '%s-12-31 23:59:59'%year
                    dtimes    = pd.date_range(startdate, enddate, freq=res)
                    na_vals   = np.array(len(dtimes)*[np.nan])
                    var      += [pd.Series(na_vals, index=dtimes)]
                    print '   ',site, year, '%2i'%crop_no, 'skip.'
                
            final_series[site][varname] = pd.concat(var)
        #final_series[site]['GPP'].plot()
        #plt.show()

    # store the final WOFOST timeseries
    filepath = os.path.join(outputdir,'%s_timeseries_'%resolution+\
               '%s_R10=%s_WOFOST_crop_rotation.pickle'%(TER_method,R10))
    pickle_dump(final_series, open(filepath,'wb'))
    print 'successfully dumped %s'%filepath

#===============================================================================
if __name__=='__main__':
    main()
#===============================================================================

