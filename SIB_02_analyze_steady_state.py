#!/usr/bin/env python

import sys, os
import numpy as np
import netCDF4 as cdf
import pandas as pd
from matplotlib import pyplot as plt

#===============================================================================
# This script analyzes if we have arrived at steady-state for the intial C pools
# of SiBCASA
def main():
#===============================================================================

    cwdir       = os.getcwd()
    stdstatedir = os.path.join(cwdir, 'search_steady_cropland')

    sites = ['BE-Lon','DE-Kli','FR-Gri','IT-BCi','NL-Dij','NL-Lan']

    # time axis of all time series
    plt.close('all')
    tm = pd.date_range('1999-01-01 00:00:00', '1999-12-31 23:59:59', freq='1d')

    for site in sites:

        fig, ax = plt.subplots(nrows=1, ncols=1)
        fig.suptitle(site, fontsize=14)

        for it in range(0,10):
            namefile = '%s_run%i/'%(site,it) +'hsib_1999*.qp2.nc'
            pathfile = os.path.join(stdstatedir, namefile)
            print pathfile
            # open all 12 files
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

            #carb_pool_soil = pd.Series([l[0] for l in Sib_csoil], index=tm)
            series_gpp = pd.Series([l[0] for l in Sib_gpp], index=tm)
            series_ter = pd.Series([l[0] for l in Sib_ter], index=tm)
            series_nee = pd.Series([l[0] for l in Sib_nee], index=tm)
            if it==0:
                series_gpp.plot()
                series_ter.plot()
                series_nee.plot(label=it)
            else:
                series_nee.plot(label=it)

            # calculate the sum (or average?) gpp and nee over the
            # year: we reach steady state when the nee/gpp < 1%
            sum1 = sum(Sib_gpp)[0]
            sum2 = sum(Sib_nee)[0]
            ave1 = np.mean(Sib_gpp)
            ave2 = np.mean(Sib_nee)
            print '\n',site,'run %02d'%it
            print "from DOY 150 to 200, average GPP, TER, NEE per "+\
                   "day (gC/m2/d) = %.2f ; %.2f ; %.2f"%(
                   np.mean(Sib_gpp[150:200]), 
                   np.mean(Sib_ter[150:200]), 
                   np.mean(Sib_nee[150:200]))
            print "Average daily NEE, GPP, NEE/GPP over the year"+\
                  "(gC/m2/d):", ave2, ave1, ave2/ave1*100.
            print "Cumulative NEE, GPP, NEE/GPP (gC/m2/y):", \
                  sum2, sum1, sum2/sum1*100.
            # if we meet our criteria, we have reached steady-state
            #if (sum2/sum1<0.01): 
            if (ave2/ave1<0.01): 
                print "\nWe have reached steady-state at run %02d"%it
                #break

        ax.legend()
        plt.show()


#===============================================================================
def build_years_dict():
#===============================================================================
    years = dict()
    years['Winter wheat'] = dict()
    years['Grain maize']  = dict()
    # list of years for winter wheat sites:
    years['Winter wheat']['BE-Lon'] = [2005,2007] #winter wheat rotation
    years['Winter wheat']['DE-Kli'] = [2006]
    years['Winter wheat']['FR-Aur'] = [2006] #not sure between 2005 and 2006
    years['Winter wheat']['FR-Gri'] = [2006]
    years['Winter wheat']['FR-Lam'] = [2007]
    # list of years for grain maize sites:
    years['Grain maize']['DE-Kli']  = [2007]
    years['Grain maize']['FR-Gri']  = [2005]
    years['Grain maize']['FR-Lam']  = [2006]
    years['Grain maize']['IT-BCi']  = [2004]
    years['Grain maize']['NL-Lan']  = [2005]
    years['Grain maize']['NL-Dij']  = [2007]

    return years


#===============================================================================
if __name__=='__main__':
    main()
#===============================================================================
