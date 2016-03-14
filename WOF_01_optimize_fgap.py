#!/usr/bin/env python

import sys, os, rc
import numpy as np
from cPickle import load as pickle_load
from cPickle import dump as pickle_dump

#===============================================================================
# This script does some standard analysis on FluxNet sites
def main():
#===============================================================================
    global inputdir, codedir, outputdir, CGMSdir, ECMWFdir, optimidir,\
           EUROSTATdir, custom_yns
#-------------------------------------------------------------------------------
# ================================= USER INPUT =================================

# read the settings from the rc file
    rcdict     = rc.read('settings.rc')

#===============================================================================
#-------------------------------------------------------------------------------
# extract the needed information from the rc file
    sites      = [s.strip(' ') for s in rcdict['sites'].split(',')]
    crops      = [s.strip(' ') for s in rcdict['crops'].split(',')]
    crop_nos   = [int(s.strip(' ')) for s in rcdict['crop_nos'].split(',')]
    years      = [int(s.strip(' ')) for s in rcdict['years'].split(',')]

    # optimization settings
    force_optimization = str_to_bool(rcdict['force_optimization'])
    selec_method  = rcdict['selec_method']
    ncells        = int(rcdict['ncells'])
    nsoils        = int(rcdict['nsoils'])
    weather       = rcdict['weather']

    # directory paths
    outputdir  = rcdict['outputdir']
    inputdir   = rcdict['inputdir']
    codedir    = rcdict['codedir']
    CGMSdir     = os.path.join(inputdir, 'CGMS')
    ECMWFdir    = os.path.join(inputdir, 'ECMWF')
    EUROSTATdir = os.path.join(inputdir, 'EUROSTATobs')

#-------------------------------------------------------------------------------
    # get the list of NUTS 2 region names associated to the list of FluxNet sites
    from WOF_00_retrieve_input_data import open_csv
    sitdict = open_csv(inputdir, 'sites_info2.csv', convert_to_float=False)
    NUTS_reg  = sitdict['NUTS_reg']
    #custom_yieldnsow =
#-------------------------------------------------------------------------------
    # get local yield and sowing date information
    import xlrd
    from xlrd.xldate import xldate_as_datetime

    xl_workbook=xlrd.open_workbook(os.path.join(inputdir,'site_yields.xlsx'))
    sheet_names = xl_workbook.sheet_names()
    xl_sheet = xl_workbook.sheet_by_name(sheet_names[0])
    xl_sites = xl_sheet.col(0)
    xl_years = xl_sheet.col(1)
    xl_crops = xl_sheet.col(2)
    xl_yield = xl_sheet.col(3)
    xl_sowda = xl_sheet.col(9)
    datemode = xl_workbook.datemode
    custom_yns = []
    for si,ye,cr,so,yi in zip(xl_sites[1:38], xl_years[1:38], xl_crops[1:38], 
                                                 xl_sowda[1:38], xl_yield[1:38]): 
        sit = str(si.value) 
        yea = int(ye.value) 
        cro = int(cr.value)
        if int(so.value) != -9999: sow = xldate_as_datetime(so.value, datemode)
        else: sow = np.nan
        if int(yi.value) != -9999.: yie = yi.value
        else: yie = np.nan
        custom_yns += [(sit, yea, cro, sow, yie)]
    
    for row in custom_yns: print row
#-------------------------------------------------------------------------------
# optimize fgap at the location / year / crops specified by user

    for s,site in enumerate(sites):

        for c,crop_name in enumerate(crops):
            crop_no = crop_nos[c]

            for year in years:

                # create output folder if it doesn't already exists
                optimidir = os.path.join(outputdir,'fgap/%i/c%i/'%(year,crop_no))
                if not os.path.exists(optimidir):
                    print 'creating new directory %s'%optimidir
                    os.makedirs(optimidir)

                # we try to optimize fgap for the NUTS 2, 1, 0 regions 
                for NUTS_level in range(3):
                    NUTS_no =  NUTS_reg[s][0:4-NUTS_level]

                    print '\n', site, NUTS_no, year, crop_name
                    
                    # OPTIMIZATION OF FGAP:
                    yldgapf = optimize_fgap(site, crop_no, crop_name, year, NUTS_no, 
                                            selec_method, ncells, nsoils, 
                                            weather, force_optimization)


#===============================================================================
# Function to optimize fgap for one NUTS region
def optimize_fgap(site, crop_no, crop_name, year, NUTS_no, selec_method, ncells, 
                                           nsoils, weather, force_optimization):
#===============================================================================
# Temporarily add code directory to python path, to be able to import pcse
# modules
    sys.path.insert(0, codedir) 
    sys.path.insert(0, os.path.join(codedir,'carbon_cycle')) 
#-------------------------------------------------------------------------------
    import glob
    from maries_toolbox import define_opti_years,\
                               select_cells, select_soils
#-------------------------------------------------------------------------------
    # if the optimization has already been performed and we don't want
    # to redo it, we skip that region
    filepath = os.path.join(optimidir,'fgap_%s_optimized.pickle'%NUTS_no)
    if (os.path.exists(filepath) and force_optimization==False):
        optimum = pickle_load(open(filepath,'rb'))
        print "We have already calculated the optimum fgap for that "+\
              "year and crop: fgap=%.2f"%optimum[2]
        return optimum[2]

#-------------------------------------------------------------------------------
    # we select the grid cell of the FluxNet site
    gridlist = pickle_load(open(os.path.join(CGMSdir,
                                 'gridlist_objects/shortgridlist.pickle'),'rb'))
    selected_grid_cells = gridlist[NUTS_no]
#-------------------------------------------------------------------------------
    # where possible, we retrieve local information about yield and sowing date
    local_sowda = None
    local_yield = None
    for row in custom_yns:
        if row[0]==site and row[1]==year and row[2]==crop_no:
            local_sowda = row[3]
            local_yield = row[4]
            print 'We recovered local info from site %s:'%site
            print 'sowing date of %s:'%crop_name, local_sowda, 'grain yield: %.3f'%local_yield
            break
    if local_sowda==None and local_yield==None:
        print 'No local information on sowing date and yield.'
#-------------------------------------------------------------------------------
# we retrieve the EUROSTAT pre-processed yield observations:
    if local_sowda==None and local_yield==None:
    try:
        filename1   = os.path.join(EUROSTATdir, 'preprocessed_yields.pickle')
        yields_dict = pickle_load(open(filename1,'rb'))
    except IOError:
        print '\nYou have not preprocessed the EUROSTAT observations'
        print 'Run the script 03_preprocess_obs.py first!\n'
        sys.exit() 
    # NB: we do NOT detrend the yields anymore, since fgap is not supposed to be
    # representative of multi-annual gap
    obs_yields = yields_dict[crop_name][NUTS_no]
    return None
#-------------------------------------------------------------------------------
	# if there were no reported yield on the year X, we skip that region
    if (year not in obs_yields[1]):
        print 'No reported yield, we have to gap-fill later'
        filename = os.path.join(optimidir,'fgap_%s_tobegapfilled.pickle'%NUTS_no)
        outlist = [NUTS_no, 2, 1., selected_grid_cells]
        pickle_dump(outlist, open(filename,'wb'))
        return 1.
#-------------------------------------------------------------------------------
    # NB: in the optimization routine, we use the observed cultivation
    # fraction of the crop to calculate the soil cultivated areas, and
    # thus to compute the regional yields (= weighted average of yields
    # using soil cultivated areas)

    # if the observed cultivated fraction is zero, we skip that region
    selected_soil_types = select_soils(crop_no,[g for g,a in selected_grid_cells],
                                       CGMSdir, method=selec_method, n=nsoils)
    print 'we selected grid cell %i, top %i soil types, for optimization'%(
                                              selected_grid_cells[0][0], nsoils)

#-------------------------------------------------------------------------------
    # we set the optimization code (gives us info on how we optimize)
    opti_code = 1 # 1= observations are available for optimization
                  # 2= no obs available 

    #print obs_yields[1], obs_yields[0]
    # in all other cases, we optimize the yield gap factor
    optimum = optimize_regional_yldgapf_dyn(NUTS_no, obs_yields,
                                                            crop_no,
                                                            selected_grid_cells,
                                                            selected_soil_types,
                                                            weather,
                                                            CGMSdir,
                                                            [year],
                                                            obs_type='yield',
                                                            plot_rmse=False)

    # pickle the information per NUTS region
    outlist = [NUTS_no, opti_code, optimum, selected_grid_cells]
    filename = os.path.join(optimidir,'fgap_%s_optimized.pickle'%NUTS_no)
    pickle_dump(outlist, open(filename,'wb'))

    return optimum

#===============================================================================
# Function to optimize the regional yield gap factor using the difference
# between the regional simulated and the observed harvest or yield (ie. 1 gap to
# optimize per NUTS region). This function iterates dynamically to find the
# optimum YLDGAPF.
def optimize_regional_yldgapf_dyn(NUTS_no_, detrend, crop_no_, 
    selected_grid_cells_, selected_soil_types_, weather, inputdir, opti_years_, 
    obs_type='yield', plot_rmse=False):
#===============================================================================

    import math
    from operator import itemgetter as operator_itemgetter
    from matplotlib import pyplot as plt
    from pcse.models import Wofost71_WLP_FD
    from pcse.base_classes import WeatherDataProvider
    from pcse.fileinput.cabo_weather import CABOWeatherDataProvider

    # aggregated yield method:
    
    # 2- we construct a 2D array with same dimensions as TSO_regional,
    # containing the observed yields
    row = [] # this list will become the row of the 2D array
    for y,year in enumerate(opti_years_):
        index_year = np.argmin(np.absolute(detrend[1]-year))
        row = row + [detrend[0][index_year]]
    OBS = np.tile(row, (5,1)) # repeats the list as a row 3 times, to get a 
                              # 2D array

    # 3- we calculate all the individual yields from the selected grid cells x
    # soils combinations

    # NB: we explore the range of yldgapf between 0.1 and 1.
    f0  = 0.
    f2  = 0.5
    f4  = 1.
    f_step  = 0.25 
    # Until the precision of the yield gap factor is good enough (i.e. < 0.02)
    # we loop over it. We do 12 iterations in total with this method.
    iter_no = 0
    RMSE_stored = list()
    while (f_step >= 0.02):

        iter_no = iter_no + 1
        # sub-method: looping over the yield gap factors

        # we build a range of 3 yield gap factors to explore one low bound, one
        # high bound, one in the middle
        f_step = (f4 - f0)/4.
        f1 = f0 + f_step
        f3 = f2 + f_step
        f_range = [f0, f1, f2, f3, f4]

        RES = [] # list in which we will store the yields of the combinations

        counter=0
        for grid, arable_land in selected_grid_cells_:
 
            frac_arable = arable_land / 625000000.

            # Retrieve the weather data of one grid cell (all years are in one
            # file) 
            if (weather == 'CGMS'):
                filename = os.path.join(inputdir,'weather_objects/',
                           'weatherobject_g%d.pickle'%grid)
                weatherdata = WeatherDataProvider()
                weatherdata._load(filename)
            if (weather == 'ECMWF'):
                weatherdata = CABOWeatherDataProvider('%i'%grid,fpath=ECMWFdir)
                        
            # Retrieve the soil data of one grid cell (all possible soil types) 
            filename = os.path.join(inputdir,'soildata_objects/',
                       'soilobject_g%d.pickle'%grid)
            soil_iterator = pickle_load(open(filename,'rb'))

            for smu, stu_no, weight, soildata in selected_soil_types_[grid]:

                # TSO will store all the yields of one grid cell x soil 
                # combination, for all years and all 3 yldgapf values
                TSO = np.zeros((len(f_range), len(opti_years_)))

                counter +=1
        
                for y, year in enumerate(opti_years_): 

                    # Retrieve yearly data 
                    filename = os.path.join(inputdir,
                               'timerdata_objects/%i/c%i/'%(year,crop_no_),
                               'timerobject_g%d_c%d_y%d.pickle'\
                                                           %(grid,crop_no_,year))
                    timerdata = pickle_load(open(filename,'rb'))
                    filename = os.path.join(inputdir,
                               'cropdata_objects/%i/c%i/'%(year,crop_no_),
                               'cropobject_g%d_c%d_y%d.pickle'\
                                                           %(grid,crop_no_,year))
                    cropdata = pickle_load(open(filename,'rb'))
                    if str(grid).startswith('1'):
                        dum = str(grid)[0:2]
                    else:
                        dum = str(grid)[0]
                    filename = os.path.join(inputdir,
                               'sitedata_objects/%i/c%i/grid_%s/'
                                                          %(year,crop_no_,dum),
                               'siteobject_g%d_c%d_y%d_s%d.pickle'\
                                                   %(grid,crop_no_,year,stu_no))
                    sitedata = pickle_load(open(filename,'rb'))

                    for f,factor in enumerate(f_range):
            
                        cropdata['YLDGAPF']=factor
                       
                        # run WOFOST
                        wofost_object = Wofost71_WLP_FD(sitedata, timerdata,
                                                soildata, cropdata, weatherdata)
                        wofost_object.run_till_terminate()
        
                        # get the yield (in kgDM.ha-1) 
                        TSO[f,y] = wofost_object.get_variable('TWSO')

                    #print grid, stu_no, year, counter, [y[0] for y in TSO], OBS[0]
                RES = RES + [(grid, stu_no, weight*frac_arable, TSO)]

        # 4- we aggregate the yield or harvest into the regional one with array
        # operations

        sum_weighted_vals = np.zeros((len(f_range), len(opti_years_)))
                                    # empty 2D array with same dimension as TSO
        sum_weights       = 0.
        for grid, stu_no, weight, TSO in RES:
            # adding weighted 2D-arrays in the empty array sum_weighted_yields
            # NB: variable 'weight' is actually the cultivated area in m2
            sum_weighted_vals   = sum_weighted_vals + (weight/10000.)*TSO 
            # computing the total sum of the cultivated area in ha 
            sum_weights         = sum_weights       + (weight/10000.) 

        if (obs_type == 'harvest'):
            TSO_regional = sum_weighted_vals / 1000000. # sum of the individual 
                                                        # harvests in 1000 tDM
        elif (obs_type == 'yield'):
            TSO_regional = sum_weighted_vals / sum_weights # weighted average of 
                                                        # all yields in kgDM/ha

        # 5- we compute the (sim-obs) differences.
        DIFF = TSO_regional - OBS
        if (TSO_regional[-1][0] <= 0.):
            print 'WARNING: no simulated crop growth. We set the optimum fgap to 1.'
            return 1.
        if (TSO_regional[-1] <= OBS[-1]):
            print 'WARNING: obs yield > sim yield. We set optimum to 1.'
            return 1.
        
        # 6- we calculate the RMSE (root mean squared error) of the 3 yldgapf
        # The RMSE of each yldgapf is based on N obs-sim differences for the N
        # years looped over

        RMSE = np.zeros(len(f_range))
        for f,factor in enumerate(f_range):
            list_of_DIFF = []
            for y, year in enumerate(opti_years_):
                list_of_DIFF = list_of_DIFF + [DIFF[f,y]]
            RMSE[f] = np.sqrt(np.mean( [ math.pow(j,2) for j in
                                                           list_of_DIFF ] ))
        #print RMSE, f_range
        # We store the value of the RMSE for plotting purposes
        RMSE_stored = RMSE_stored + [(f_range[1], RMSE[1]), (f_range[3], RMSE[3])]
        if (iter_no == 1):
            RMSE_stored = RMSE_stored + [(f_range[0], RMSE[0]), 
                                         (f_range[2], RMSE[2]),
                                         (f_range[4], RMSE[4])]

        # 7- We update the yldgapf range to explore for the next iteration. 
        # For this we do a linear interpolation of RMSE between the 3 yldgapf
        # explored here, and the next range to explore is the one having the
        # smallest interpolated RMSE

        index_new_center = RMSE.argmin()
        # if the optimum is close to 1:
        if index_new_center == len(f_range)-1:
            f0 = f_range[index_new_center-2]
            f2 = f_range[index_new_center-1]
            f4 = f_range[index_new_center]
        # if the optimum is close to 0:
        elif index_new_center == 0:
            f0 = f_range[index_new_center]
            f2 = f_range[index_new_center+1]
            f4 = f_range[index_new_center+2]
        else:
            f0 = f_range[index_new_center-1]
            f2 = f_range[index_new_center]
            f4 = f_range[index_new_center+1]

	# when we are finished iterating on the yield gap factor range, we plot the
    # RMSE as a function of the yield gap factor
    if (plot_rmse == True):
        RMSE_stored  = sorted(RMSE_stored, key=operator_itemgetter(0))
        x,y = zip(*RMSE_stored)
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,5))
        fig.subplots_adjust(0.15,0.16,0.95,0.96,0.4,0.)
        ax.plot(x, y, c='k', marker='o')
        ax.set_xlabel('yldgapf (-)')
        ax.set_ylabel('RMSE')
        fig.savefig('%s_opti_fgap.png'%NUTS_no_)
        #pickle_dump(RMSE_stored,open('%s_RMSE.pickle'%NUTS_no_,'wb'))

    # 8- when we are finished iterating on the yield gap factor range, we return
    # the optimum value. We look for the yldgapf with the lowest RMSE
    index_optimum   = RMSE.argmin()
    optimum_yldgapf = f_range[index_optimum] 

    print 'optimum found: %.2f +/- %.2f'%(optimum_yldgapf, f_step)

    # 10- we return the optimized YLDGAPF
    return optimum_yldgapf


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
