#!/usr/bin/env python

import sys, os
import numpy as np
import netCDF4 as cdf
import pandas as pd
from matplotlib import pyplot as plt
from cPickle import dump as pickle_dump

#===============================================================================
# This script analyzes if we have arrived at steady-state for the intial C pools
# of SiBCASA
def main():
#===============================================================================

    cwdir     = os.getcwd()
    SIBrundir = os.path.join(cwdir, 'forward_runs')

    sites = ['BE-Lon','DE-Kli','FR-Gri','IT-BCi','NL-Dij','NL-Lan']

    # time axis of all time series
    tm = pd.date_range('2000-01-01 00:00:00', '2010-12-31 23:59:59', freq='1d')
    series = dict()

    # strange thing is: SiBCASA ignores leap years and does not simulate 29 feb
    # we have to delete 3 dates on 3 leap years between 2000 and 2010
    new_tm = tm[0:59].union(tm[60:1520].union(tm[1521:2981].union(tm[2982:])))
    print new_tm, len(new_tm)

    for site in sites:

        # open all the years and store in one list
        namefile = '%s_2000-2010/'%(site) +'hsib_*.qp2.nc'
        pathfile = os.path.join(SIBrundir, namefile)

        # open all 11 years * 12 files
        f = cdf.MFDataset(pathfile)
        # get daily GPP and NEE (in micromoles/m2/s) and convert
        # the fluxes to gC/m2/d:
        fac = 0.000001*12. # conversion from micromoles to gC
        dt  = 3600. * 24.  # nb of seconds in a day
        Sib_gpp  = np.array(-f.variables['gpp'][:])*fac*dt
        Sib_ter  = np.array(f.variables['resp_tot'][:])*fac*dt
        Sib_rhet = np.array(f.variables['resp_het'][:])*fac*dt
        Sib_raut = np.array(f.variables['resp_auto'][:])*fac*dt
        Sib_nee  = np.array(f.variables['NEE_2'][:])*fac*dt
        # from moles/m2 to gC/m2
        Sib_csoil= np.array(f.variables['carb_soil'][:])*fac*1000000.*dt
        # close file
        f.close()

        series[site] = dict()
        series[site]['GPP'] = pd.Series([l[0] for l in Sib_gpp],  index=new_tm)
        series[site]['TER'] = pd.Series([l[0] for l in Sib_ter],  index=new_tm)
        series[site]['Rhet'] = pd.Series([l[0] for l in Sib_rhet], index=new_tm)
        series[site]['Raut'] = pd.Series([l[0] for l in Sib_raut], index=new_tm)
        series[site]['NEE'] = pd.Series([l[0] for l in Sib_nee],  index=new_tm)

        fig, ax = plt.subplots(nrows=1, ncols=1)
        fig.suptitle(site, fontsize=14)
        series[site]['GPP'].plot(label='GPP')
        series[site]['TER'].plot(label='TER')
        series[site]['Rhet'].plot(label='Rhet')
        series[site]['Raut'].plot(label='Raut')
        series[site]['NEE'].plot(label='NEE')
        ax.legend()

    # store the formatted pandas timeseries in a pickle file
    filepath = os.path.join(SIBrundir,'timeseries_SiBCASA.pickle')
    pickle_dump(series, open(filepath,'wb'))

    # preview the timeseries per site
    plt.show()

#===============================================================================
if __name__=='__main__':
    main()
#===============================================================================
