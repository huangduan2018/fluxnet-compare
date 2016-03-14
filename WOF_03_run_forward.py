#!/usr/bin/env python

import sys, os, rc
import numpy as np
import pandas as pd
import netCDF4 as cdf
from cPickle import load as pickle_load
from cPickle import dump as pickle_dump

#===============================================================================
# This script does some standard analysis on FluxNet sites
def main():
#===============================================================================
    global inputdir, codedir, outputdir, CGMSdir, ECMWFdir, optimidir, forwardir,\
           EUROSTATdir, mmC, mmCO2, mmCH2O
#-------------------------------------------------------------------------------
# fixed molar masses for unit conversion of carbon fluxes
    mmC    = 12.01
    mmCO2  = 44.01
    mmCH2O = 30.03 

# ================================= USER INPUT =================================

# read the settings from the rc file
    rcdict     = rc.read('settings.rc')

#===============================================================================
#-------------------------------------------------------------------------------
# extract the needed information from the rc file
    sites      = [s.strip(' ') for s in rcdict['sites'].split(',')]
    #site_lons  = [float(s.strip(' ')) for s in rcdict['site_lons'].split(',')]
    #site_lats  = [float(s.strip(' ')) for s in rcdict['site_lats'].split(',')]
    #gridcells  = [float(s.strip(' ')) for s in rcdict['gridcells'].split(',')]
    #NUTS_reg   = [s.strip(' ') for s in rcdict['NUTS_reg'].split(',')]
    crops      = [s.strip(' ') for s in rcdict['crops'].split(',')]
    crop_nos   = [int(s.strip(' ')) for s in rcdict['crop_nos'].split(',')]
    years      = [int(s.strip(' ')) for s in rcdict['years'].split(',')]

    # forward runs settings
    force_forwardsim = str_to_bool(rcdict['force_forwardsim'])
    selec_method  = rcdict['selec_method']
    ncells        = int(rcdict['ncells'])
    nsoils        = int(rcdict['nsoils'])
    weather       = rcdict['weather']

    # carbon cycle settings
    TER_method  = rcdict['TER_method'] # if grow-only: NEE = GPP + Rgrow + Rsoil
    Eact0       = float(rcdict['Eact0'])
    R10         = float(rcdict['R10'])
    resolution  = rcdict['resolution']  # can be hourly or daily

    # directory paths
    outputdir  = rcdict['outputdir']
    inputdir   = rcdict['inputdir']
    codedir    = rcdict['codedir']
    CGMSdir     = os.path.join(inputdir, 'CGMS')
    ECMWFdir    = os.path.join(inputdir, 'ECMWF')
    EUROSTATdir = os.path.join(inputdir, 'EUROSTATobs')

#-------------------------------------------------------------------------------
    # get the sites longitude and latitudes
    from WOF_00_retrieve_input_data import open_csv
    sitdict = open_csv(inputdir, 'sites_info2.csv', convert_to_float=False)
    site_lons = [float(l) for l in sitdict['site_lons']]
    site_lats = [float(l) for l in sitdict['site_lats']]
    gridcells = [int(g) for g in sitdict['gridcells']]
    NUTS_reg  = sitdict['NUTS_reg']
#-------------------------------------------------------------------------------
# run WOFOST at the location / year / crops specified by user

    print '\nYLDGAPF(-),  grid_no,  year,  stu_no, stu_area(ha), '\
     +'TSO(kgDM.ha-1), TLV(kgDM.ha-1), TST(kgDM.ha-1), '\
     +'TRT(kgDM.ha-1), maxLAI(m2.m-2), rootdepth(cm), TAGP(kgDM.ha-1)'

    # we format the time series using the pandas python library, for easy plotting
    startdate = '%i-01-01 00:00:00'%years[0]
    enddate   = '%i-12-31 23:59:59'%years[-1]
    if resolution == 'daily':
        dtimes = pd.date_range(startdate, enddate, freq='1d')
    elif resolution == '3-hourly':
        dtimes = pd.date_range(startdate, enddate, freq='3H')
    else:
        print "Wrong CO2 fluxes temporal resolution: must be either 'daily' or '3-hourly'"
        sys.exit()

    series = dict()
    for s,site in enumerate(sites):
        lon = site_lons[s]
        lat = site_lats[s]
        grid_no = gridcells[s]
        NUTS_no = NUTS_reg[s]
        series[site] = dict()

        for c,crop_name in enumerate(crops):
            cpno = crop_nos[c]
            series[site]['c%i'%cpno] = dict()
            list_of_gpp  = np.array([])
            list_of_raut = np.array([])
            list_of_rhet = np.array([])
            list_of_ter  = np.array([])
            list_of_nee  = np.array([])

            for year in years:
                # create output folder if it doesn't already exists
                optimidir = os.path.join(outputdir,'fgap/%i/c%i/'%(year,cpno))

                # create output folder if it doesn't already exists
                forwardir = os.path.join(outputdir,'forward_runs/%i/c%i/'%(year,
                                                                        cpno))
                if not os.path.exists(forwardir):
                    os.makedirs(forwardir)

                print '\n', site, NUTS_no, year, crop_name

                # RETRIEVE OPTIMUM FGAP:
                # either the NUTS2 optimum if it exists
                ygf_path  = os.path.join(optimidir,'fgap_%s_optimized.pickle'%NUTS_no)
                # or the gapfilled version
                if not os.path.exists(ygf_path):
                    ygf_file  = [f for f in os.listdir(optimidir) if (NUTS_no in f) 
                                and ('_gapfilled' in f)][0]
                    ygf_path = os.path.join(optimidir, ygf_file)
                fgap_info = pickle_load(open(ygf_path,'rb'))
                yldgapf   = fgap_info[2]

                # FORWARD SIMULATIONS:
                perform_yield_sim(cpno, grid_no, int(year), yldgapf, 
                                  selec_method, nsoils, force_forwardsim)
                # POST-PROCESSING OF GPP, RAUTO, RHET, NEE:
                SimData = compute_timeseries_fluxes(cpno, grid_no, lon, lat, 
                                                    year, R10, Eact0, selec_method, 
                                                    nsoils, TER_method=TER_method,
                                                    scale=resolution)
                list_of_gpp  = np.concatenate([list_of_gpp,  SimData[1]], axis=0)
                list_of_raut = np.concatenate([list_of_raut, SimData[2]], axis=0)
                list_of_rhet = np.concatenate([list_of_rhet, SimData[3]], axis=0)
                list_of_ter  = np.concatenate([list_of_ter,  SimData[4]], axis=0)
                list_of_nee  = np.concatenate([list_of_nee,  SimData[5]], axis=0)

            print dtimes, list_of_gpp
            
            series[site]['c%i'%cpno]['GPP']  = pd.Series(list_of_gpp,  index=dtimes)
            series[site]['c%i'%cpno]['Raut'] = pd.Series(list_of_raut, index=dtimes)
            series[site]['c%i'%cpno]['Rhet'] = pd.Series(list_of_rhet, index=dtimes)
            series[site]['c%i'%cpno]['TER']  = pd.Series(list_of_ter,  index=dtimes)
            series[site]['c%i'%cpno]['NEE']  = pd.Series(list_of_nee,  index=dtimes)

    # we store the two pandas series in one pickle file
    filepath = os.path.join(outputdir,'forward_runs/'+\
               '%s_timeseries_%s_WOFOST.pickle'%(resolution,TER_method))
    pickle_dump(series, open(filepath,'wb'))



#===============================================================================
# Function to do forward simulations of crop yield for a given YLDGAPF and for a
# selection of grid cells x soil types within a NUTS region
def perform_yield_sim(crop_no, grid_no, year, fgap, selec_method, nsoils, force_forwardsim):
#===============================================================================
# Temporarily add code directory to python path, to be able to import pcse
# modules
    sys.path.insert(0, codedir) 
    sys.path.insert(0, os.path.join(codedir,'carbon_cycle')) 
#-------------------------------------------------------------------------------
    import glob
    from pcse.fileinput.cabo_weather import CABOWeatherDataProvider
    from maries_toolbox import select_soils
    from pcse.models import Wofost71_WLP_FD
    from pcse.exceptions import WeatherDataProviderError
#-------------------------------------------------------------------------------
    # fixed settings for these point simulations:
    weather = 'ECMWF'   
#-------------------------------------------------------------------------------
    # skipping already performed forward runs if required by user
    outlist = glob.glob(os.path.join(forwardir,'wofost_g%i*'%grid_no))
    if (len(outlist)==nsoils and force_forwardsim==False):
        print "        We have already done that forward run! Skipping."
        return None
#-------------------------------------------------------------------------------
    # Retrieve the weather data of one grid cell
    if (weather == 'CGMS'):
        filename = os.path.join(CGMSdir,'weatherobject_g%d.pickle'%grid_no)
        weatherdata = WeatherDataProvider()
        weatherdata._load(filename)
    if (weather == 'ECMWF'):
        weatherdata = CABOWeatherDataProvider('%i'%(grid_no), 
                                                         fpath=ECMWFdir)
    #print weatherdata(datetime.date(datetime(2006,4,1)))

    # Retrieve the soil data of one grid cell 
    filename = os.path.join(CGMSdir,'soildata_objects','soilobject_g%d.pickle'%grid_no)
    soil_iterator = pickle_load(open(filename,'rb'))

    # Retrieve calendar data of one year for one grid cell
    filename = os.path.join(CGMSdir,'timerdata_objects/%i/c%i/'%(year,crop_no),
                      'timerobject_g%d_c%d_y%d.pickle'%(grid_no, crop_no, year))
    timerdata = pickle_load(open(filename,'rb'))
                    
    # Retrieve crop data of one year for one grid cell
    filename = os.path.join(CGMSdir,'cropdata_objects/%i/c%i/'%(year,crop_no),
                         'cropobject_g%d_c%d_y%d.pickle'%(grid_no,crop_no,year))
    cropdata = pickle_load(open(filename,'rb'))

    # retrieve the fgap data of one year and one grid cell
    cropdata['YLDGAPF'] = fgap

    # Select soil types to loop over for the forward runs
    selected_soil_types = select_soils(crop_no, [grid_no], CGMSdir, 
                                       method=selec_method, n=nsoils)

    for smu, stu_no, stu_area, soildata in selected_soil_types[grid_no]:
        
        resfile = os.path.join(forwardir,"wofost_g%i_s%i.txt"%(grid_no,stu_no))

        # Retrieve the site data of one year, one grid cell, one soil type
        if str(grid_no).startswith('1'): 
            dum = str(grid_no)[0:2]
        else:
            dum = str(grid_no)[0:1]
        filename = os.path.join(CGMSdir,'sitedata_objects/%i/c%i/grid_%s'%(year,
                      crop_no,dum),'siteobject_g%d_c%d_y%d_s%d.pickle'%(grid_no,
                                                           crop_no,year,stu_no))
        sitedata = pickle_load(open(filename,'rb'))

        # run WOFOST
        wofost_object = Wofost71_WLP_FD(sitedata, timerdata, soildata, cropdata, 
                                                                    weatherdata)
        try:
            wofost_object.run_till_terminate() #will stop the run when DVS=2
        except WeatherDataProviderError:
            print 'Error with the weather data'
            return None

        # get time series of the output and take the selected variables
        wofost_object.store_to_file(resfile)

        # get major summary output variables for each run
        # total dry weight of - dead and alive - storage organs (kg/ha)
        TSO       = wofost_object.get_variable('TWSO')
        # total dry weight of - dead and alive - leaves (kg/ha) 
        TLV       = wofost_object.get_variable('TWLV')
        # total dry weight of - dead and alive - stems (kg/ha) 
        TST       = wofost_object.get_variable('TWST')
        # total dry weight of - dead and alive - roots (kg/ha) 
        TRT       = wofost_object.get_variable('TWRT')
        # maximum LAI
        MLAI      = wofost_object.get_variable('LAIMAX')
        # rooting depth (cm)
        RD        = wofost_object.get_variable('RD')
        # Total above ground dry matter (kg/ha)
        TAGP      = wofost_object.get_variable('TAGP')

        #output_string = '%10.3f, %8i, %5i, %7i, %15.2f, %12.5f, %14.2f, '
                        #%(yldgapf, grid_no, year, stu_no, arable_area/10000.,stu_area/10000.,TSO) 
        output_string = '%10.3f, %8i, %5i, %7i, %12.5f, %14.2f, '%(fgap,
                         grid_no, year, stu_no, stu_area/10000., TSO) + \
                        '%14.2f, %14.2f, %14.2f, %14.2f, %13.2f, %15.2f'%(TLV,
                         TST, TRT, MLAI, RD, TAGP)
        print output_string

    return None


#===============================================================================
def compute_timeseries_fluxes(crop_no, grid_no, lon, lat, year, R10, Eact0, 
    selec_method, nsoils, TER_method='grow-only',scale='daily'):
# possible methods: 'grow-only': NEE = GPP + Rgrow + Rsoil
#                   'rauto':     NEE = GPP + Rgrow + Rmaint + Rsoil
#===============================================================================
        
    import math
    import pandas as pd
    import datetime as dt
    from maries_toolbox import open_pcse_csv_output, select_soils

    print '- grid cell no %i: lon = %.2f , lat = %.2f'%(grid_no,lon,
                                                                lat)

    prod_figure = False

    # we retrieve the tsurf and rad variables from ECMWF
    filename_rad = 'rad_ecmwf_%i_lon%.2f_lat%.2f.pickle'%(year,lon,lat)
    path_rad     = os.path.join(ECMWFdir,filename_rad)
    rad = pickle_load(open(path_rad, 'rb'))

    filename_ts = 'ts_ecmwf_%i_lon%.2f_lat%.2f.pickle'%(year,lon,lat)
    path_ts     = os.path.join(ECMWFdir,filename_ts)
    ts = pickle_load(open(path_ts, 'rb'))

    
    # we initialize the timeseries for the grid cell
    # time list for timeseries
    time_cell_persec_timeseries = rad[0]
    time_cell_perday_timeseries = rad[0][::8]/(3600.*24.)
    # length of all carbon fluxes timeseries
    len_persec = len(rad[0])
    len_perday = len(rad[0][::8])
    # GPP timeseries
    gpp_cell_persec_timeseries  = np.array([0.]*len_persec)
    gpp_cell_perday_timeseries  = np.array([0.]*len_perday)
    # autotrophic respiration timeseries
    raut_cell_persec_timeseries = np.array([0.]*len_persec)
    raut_cell_perday_timeseries = np.array([0.]*len_perday)
    # heterotrophic respiration timeseries
    rhet_cell_persec_timeseries = np.array([0.]*len_persec)
    rhet_cell_perday_timeseries = np.array([0.]*len_perday)

    # we initialize some variables
    sum_stu_areas = 0. # sum of soil types areas
    delta = 3600. * 3. # number of seconds in delta (here 3 hours)

    if (prod_figure == True):
        from matplotlib import pyplot as plt
        plt.close('all')
        fig1, axes = plt.subplots(nrows=3, ncols=1, figsize=(14,10))
        fig1.subplots_adjust(0.1,0.1,0.98,0.9,0.2,0.2)

    # Select soil types to loop over
    soilist = select_soils(crop_no, [grid_no],
                            CGMSdir, method=selec_method, n=nsoils)

    #---------------------------------------------------------------
    # loop over soil types
    #---------------------------------------------------------------
    for smu, stu_no, stu_area, soildata in soilist[grid_no]:

        # We open the WOFOST results file
        filename    = 'wofost_g%i_s%i.txt'%(grid_no, stu_no) 
        results_set = open_pcse_csv_output(forwardir, [filename])
        wofost_data = results_set[0]

        # We apply the short wave radiation diurnal cycle on the GPP 
        # and R_auto

        # we create empty time series for this specific stu
        gpp_cycle_timeseries   = np.array([])
        raut_cycle_timeseries  = np.array([])
        gpp_perday_timeseries  = np.array([])
        raut_perday_timeseries = np.array([])
     
        # we compile the sum of the stu areas to do a weighted average of
        # GPP and Rauto later on
        sum_stu_areas += stu_area 
     
        #-----------------------------------------------------------
        # loop over days of the year
        #-----------------------------------------------------------
        for DOY, timeinsec in enumerate(time_cell_persec_timeseries[::8]):
            # conversion of current time in seconds into date
            time = dt.date(year,1,1) + dt.timedelta(DOY)
            #print 'date:', time
     
            # we test to see if we are within the growing season
            test_sow = (time - wofost_data[filename]['day'][0]).total_seconds()
            test_rip = (time - wofost_data[filename]['day'][-1]).total_seconds() 
            #print 'tests:', test_sow, test_rip
     
            # if the day of the time series is before sowing date: plant 
            # fluxes are set to zero
            if test_sow < 0.: 
                gpp_day  = 0.
                raut_day = 0.
            # or if the day of the time series is after the harvest date: 
            # plant fluxes are set to zero
            elif test_rip > 0.: 
                gpp_day  = 0.
                raut_day = 0.
            # else we get the daily total GPP and Raut in kgCH2O/ha/day
            # from wofost, and we weigh it with the stu area to later on 
            # calculate the weighted average GPP and Raut in the grid cell
            else: 
                # index of the sowing date in the time_cell_timeseries:
                if (test_sow == 0.): DOY_sowing  = DOY
                if (test_rip == 0.): DOY_harvest = DOY
                #print 'DOY sowing:', DOY_sowing
                # translation of cell to stu timeseries index
                index_day_w  = DOY - DOY_sowing
                #print 'index of day in wofost record:', index_day_w
     
                # unit conversion: from kgCH2O/ha/day to gC/m2/day
                gpp_day    = - wofost_data[filename]['GASS'][index_day_w] * \
                                                        (mmC / mmCH2O) * 0.1
                maint_resp = wofost_data[filename]['MRES'][index_day_w] * \
                                                        (mmC / mmCH2O) * 0.1
                try: # if there are any available assimilates for growth
                    growth_fac = (wofost_data[filename]['DMI'][index_day_w]) / \
                             (wofost_data[filename]['GASS'][index_day_w] - 
                              wofost_data[filename]['MRES'][index_day_w])
                    growth_resp = (1.-growth_fac)*(-gpp_day-maint_resp) 
                except ZeroDivisionError: # otherwise there is no crop growth
                    growth_resp = 0.
                if TER_method == 'rauto':
                    raut_day   = growth_resp + maint_resp
                elif TER_method == 'grow-only':
                    raut_day   = growth_resp
     
            # we select the radiation diurnal cycle for that date
            # NB: the last index is ignored in the selection, so we DO have
            # 8 time steps selected only (it's a 3-hourly dataset)
            rad_cycle      = rad[1][DOY*8:DOY*8+8] 
     
            # we apply the radiation cycle on the GPP and Rauto
            # and we transform the daily integral into a rate
            weights        = rad_cycle / sum(rad_cycle)
            # the sum of the 8 rates is equal to total/delta:
            sum_gpp_rates  = gpp_day   / delta
            sum_raut_rates = raut_day  / delta
            # the day's 8 values of actual gpp and raut rates per second:
            gpp_cycle      = weights * sum_gpp_rates
            raut_cycle     = weights * sum_raut_rates
            # NB: we check if the applied diurnal cycle is correct
            assert (sum(weights)-1. < 0.000001), "wrong radiation kernel"
            assert (len(gpp_cycle)*int(delta) == 86400), "wrong delta in diurnal cycle"
            assert ((sum(gpp_cycle)*delta-gpp_day) < 0.00001), "wrong diurnal cycle "+\
                "applied on GPP: residual=%.2f "%(sum(gpp_cycle)*delta-gpp_day) +\
                "on DOY %i"%DOY
            assert ((sum(raut_cycle)*delta-raut_day) < 0.00001), "wrong diurnal cycle "+\
                "applied on Rauto: residual=%.2f "%(sum(raut_cycle)*delta-raut_day) +\
                "on DOY %i"%DOY
     
            # if the applied diurnal cycle is ok, we append that day's cycle
            # to the yearly record of the stu
            gpp_cycle_timeseries  = np.concatenate((gpp_cycle_timeseries, 
                                                   gpp_cycle), axis=0)
            raut_cycle_timeseries = np.concatenate((raut_cycle_timeseries,
                                                   raut_cycle), axis=0)
            # we also store the carbon fluxes per day, for comparison with fluxnet
            gpp_perday_timeseries = np.concatenate((gpp_perday_timeseries,
                                                   [gpp_day]), axis=0) 
            raut_perday_timeseries = np.concatenate((raut_perday_timeseries,
                                                   [raut_day]), axis=0)

        #-----------------------------------------------------------
        # end of day nb loop
        #-----------------------------------------------------------

        # plot the soil type timeseries if requested by the user
        if (prod_figure == True):
            for ax, var, name, lims in zip(axes.flatten(), 
            [gpp_perday_timeseries, raut_perday_timeseries, 
            gpp_perday_timeseries + raut_perday_timeseries],
            ['GPP', 'Rauto', 'NPP'], [[-20.,0.],[0.,10.],[-15.,0.]]):
                ax.plot(time_cell_perday_timeseries, var, 
                                              label='stu %i'%stu_no)
                #ax.set_xlim([40.,170.])
                #ax.set_ylim(lims)
                ax.set_ylabel(name + r' (g$_{C}$ m$^{-2}$ d$^{-1}$)', 
                                                        fontsize=14)

        # We compile time series of carbon fluxes in units per day and per second
        # a- sum the PER SECOND timeseries
        gpp_cell_persec_timeseries  = gpp_cell_persec_timeseries + \
                                      gpp_cycle_timeseries*stu_area
        raut_cell_persec_timeseries = raut_cell_persec_timeseries + \
                                      raut_cycle_timeseries*stu_area

        # b- sum the PER DAY timeseries
        gpp_cell_perday_timeseries  = gpp_cell_perday_timeseries + \
                                      gpp_perday_timeseries*stu_area
        raut_cell_perday_timeseries = raut_cell_perday_timeseries + \
                                      raut_perday_timeseries*stu_area
    #---------------------------------------------------------------
    # end of soil type loop
    #---------------------------------------------------------------

    # finish ploting the soil type timeseries if requested by the user
    if (prod_figure == True):
        plt.xlabel('time (DOY)', fontsize=14)
        plt.legend(loc='upper left', ncol=2, fontsize=10)
        fig1.suptitle('Daily carbon fluxes of %s for all '%crop+\
                     'soil types of grid cell %i in %i'%(grid_no,
                                              year), fontsize=18)
        figname = 'GPP_allsoils_%i_c%i_g%i.png'%(year,crop_no,\
                                                            grid_no)
        #plt.show()
        fig1.savefig(os.path.join(analysisdir,figname))

    # compute the weighted average of GPP, Rauto over the grid cell
    # a- PER SECOND
    gpp_cell_persec_timeseries  = gpp_cell_persec_timeseries  / sum_stu_areas
    raut_cell_persec_timeseries = raut_cell_persec_timeseries / sum_stu_areas
   
    # b- PER DAY
    gpp_cell_perday_timeseries  = gpp_cell_perday_timeseries  / sum_stu_areas
    raut_cell_perday_timeseries = raut_cell_perday_timeseries / sum_stu_areas

    # compute the heterotrophic respiration with the A-gs equation
    # NB: we assume here Rhet only dependant on tsurf, not soil moisture
    #fw = Cw * wsmax / (wg + wsmin)
    tsurf_inter = Eact0 / (283.15 * 8.314) * (1 - 283.15 / ts[1])
    # a- PER SEC:
    rhet_cell_persec_timeseries = R10 * np.array([ math.exp(t) for t in tsurf_inter ]) 
    # b- PER DAY:
    for i in range(len(rhet_cell_perday_timeseries)):
        rhet_cell_perday_timeseries[i] = rhet_cell_persec_timeseries[i*8] * delta +\
                                       rhet_cell_persec_timeseries[i*8+1] * delta +\
                                       rhet_cell_persec_timeseries[i*8+2] * delta +\
                                       rhet_cell_persec_timeseries[i*8+3] * delta +\
                                       rhet_cell_persec_timeseries[i*8+4] * delta +\
                                       rhet_cell_persec_timeseries[i*8+5] * delta +\
                                       rhet_cell_persec_timeseries[i*8+6] * delta +\
                                       rhet_cell_persec_timeseries[i*8+7] * delta 
   
    # conversion from mgCO2 to gC
    conversion_fac = (mmC / mmCO2) * 0.001
    rhet_cell_persec_timeseries = rhet_cell_persec_timeseries * conversion_fac
    rhet_cell_perday_timeseries = rhet_cell_perday_timeseries * conversion_fac

    # calculate TER:
    ter_cell_persec_timeseries = raut_cell_persec_timeseries + \
                                 rhet_cell_persec_timeseries
    ter_cell_perday_timeseries = raut_cell_perday_timeseries + \
                                 rhet_cell_perday_timeseries

    # calculate NEE:
    nee_cell_persec_timeseries = gpp_cell_persec_timeseries  + \
                                 raut_cell_persec_timeseries + \
                                 rhet_cell_persec_timeseries
    nee_cell_perday_timeseries = gpp_cell_perday_timeseries  + \
                                 raut_cell_perday_timeseries + \
                                 rhet_cell_perday_timeseries

    # here we choose to return the carbon fluxes PER DAY
    if scale=='daily':
        return time_cell_perday_timeseries, gpp_cell_perday_timeseries, \
               raut_cell_perday_timeseries, rhet_cell_perday_timeseries, \
               ter_cell_perday_timeseries, nee_cell_perday_timeseries
    elif scale=='3-hourly':
        return time_cell_persec_timeseries, gpp_cell_persec_timeseries, \
               raut_cell_persec_timeseries, rhet_cell_persec_timeseries, \
               ter_cell_persec_timeseries, nee_cell_persec_timeseries


#===============================================================================
# function that will retrieve the surface temperature from the ECMWF data
# (ERA-interim). It will return two arrays: one of the time in seconds since
# 1st of Jan, and one with the tsurf variable in K.
def retrieve_ecmwf_tsurf(year, lon, lat):
#===============================================================================

    import netCDF4 as cdf

    tsurf = np.array([])
    time  = np.array([])

    for month in range (1,13):
        #print year, month
        for day in range(1,32):

            # open file if it exists
            namefile = 't_%i%02d%02d_00p03.nc'%(year,month,day)
            if (os.path.exists(os.path.join(ecmwfdir_tsurf,'%i/%02d'%(year,month),
                                                             namefile))==False):
                #print 'cannot find %s'%namefile
                continue
            pathfile = os.path.join(ecmwfdir_tsurf,'%i/%02d'%(year,month),namefile)
            f = cdf.Dataset(pathfile)

            # retrieve closest latitude and longitude index of desired location
            lats = f.variables['lat'][:]
            lons = f.variables['lon'][:]
            latindx = np.argmin( np.abs(lats - lat) )
            lonindx = np.argmin( np.abs(lons - lon) )

            # retrieve the temperature at the highest pressure level, at that 
            # lon,lat location 
            #print f.variables['ssrd'] # to get the dimensions of the variable
            tsurf = np.append(tsurf, f.variables['T'][0:8, 0, latindx, lonindx])
            # retrieve the nb of seconds on day 1 of that year
            if (month ==1 and day ==1):
                convtime = f.variables['time'][0]
            # NB: the file has 8 time steps (3-hourly)
            time = np.append(time, f.variables['time'][:] - convtime)  
            f.close()

    if (len(time) < 2920):
        print '!!!WARNING!!!'
        print 'there are less than 365 days of data that we could retrieve'
        print 'check the folder %s for year %i'%(ecmwfdir_tsurf, year)
 
    return time, tsurf

#===============================================================================
# function that will retrieve the incoming surface shortwave radiation from the
# ECMWF data (ERA-interim). It will return two arrays: one of the time in
# seconds since 1st of Jan, and one with the ssrd variable in W.m-2.
def retrieve_ecmwf_ssrd(year, lon, lat):
#===============================================================================

    import netCDF4 as cdf

    ssrd = np.array([])
    time = np.array([])

    for month in range (1,13):
        #print year, month
        for day in range(1,32):

            # open file if it exists
            namefile = 'ssrd_%i%02d%02d_00p03.nc'%(year,month,day)
            if (os.path.exists(os.path.join(ecmwfdir_ssrd,'%i/%02d'%(year,month),
                                                             namefile))==False):
                #print 'cannot find %s'%namefile
                continue
            pathfile = os.path.join(ecmwfdir_ssrd,'%i/%02d'%(year,month),namefile)
            f = cdf.Dataset(pathfile)

            # retrieve closest latitude and longitude index of desired location
            lats = f.variables['lat'][:]
            lons = f.variables['lon'][:]
            latindx = np.argmin( np.abs(lats - lat) )
            lonindx = np.argmin( np.abs(lons - lon) )

            # retrieve the shortwave downward surface radiation at that location 
            #print f.variables['ssrd'] # to get the dimensions of the variable
            ssrd = np.append(ssrd, f.variables['ssrd'][0:8, latindx, lonindx])
            # retrieve the nb of seconds on day 1 of that year
            if (month ==1 and day ==1):
                convtime = f.variables['time'][0]
            # NB: the file has 8 time steps (3-hourly)
            time = np.append(time, f.variables['time'][:] - convtime)  
            f.close()

    if (len(time) < 2920):
        print '!!!WARNING!!!'
        print 'there are less than 365 days of data that we could retrieve'
        print 'check the folder %s for year %i'%(ecmwfdir_ssrd, year)
 
    return time, ssrd

#===============================================================================
def str_to_bool(s):
#===============================================================================
    if s.strip(' ') == 'True':
         return True
    elif s.strip(' ') == 'False':
         return False
    else:
         raise ValueError

#===============================================================================
if __name__=='__main__':
    main()
#===============================================================================
