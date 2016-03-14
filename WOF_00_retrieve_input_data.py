#!/usr/bin/env python

import sys, os, rc
import numpy as np
from cPickle import load as pickle_load
from cPickle import dump as pickle_dump


#===============================================================================
# This script does some standard analysis on FluxNet sites
def main():
#===============================================================================
    global inputdir, codedir, outputdir, CGMSdir, obsdir\
#-------------------------------------------------------------------------------
    import cx_Oracle
    import sqlalchemy as sa
    from datetime import datetime
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

    obsdir     = rcdict['obsdir']
    inputdir   = rcdict['inputdir']
    CGMSdir     = os.path.join(inputdir, 'CGMS')
    codedir    = rcdict['codedir']
#-------------------------------------------------------------------------------
# get the closest CGMS grid cell id number for each FluxNet site

    # get the sites longitude and latitudes
    sitdict = open_csv(os.path.join(obsdir,'regrouped_data'), 'sites_info.txt',
                       convert_to_float=False)
    site_lons = sitdict['site_lons']
    site_lats = sitdict['site_lats']

    # we read the CGMS grid cells coordinates from file
    CGMS_cells = open_csv(CGMSdir, 'CGMS_grid_list.csv', convert_to_float=True)
    all_grids  = CGMS_cells['GRID_NO']
    all_lons   = CGMS_cells['LONGITUDE']
    all_lats   = CGMS_cells['LATITUDE']

    flux_gri = dict()
    for i,site in enumerate(sitdict['sites']):
        lon = float(site_lons[i])
        lat = float(site_lats[i])
        # compute the distance to site for all CGMS grid cells
        dist_list = list()
        for j,grid_no in enumerate(all_grids):
            distance = ((all_lons[j]-lon)**2. + (all_lats[j]-lat)**2.)**(1./2.)
            dist_list += [distance] 
        # select the closest grid cell
        indx = np.argmin(np.array(dist_list))
        flux_gri[site] = all_grids[indx]

        print 'FluxNet site %s with lon=%5.2f, lat=%5.2f: closest grid cell is %i'%(site, lon, lat, all_grids[indx])

#-------------------------------------------------------------------------------
# create new file with grid cell number in it

    filename = os.path.join(inputdir,'sites_info2.csv')
    newres = open(filename,'wb')
    oldres = open(os.path.join(obsdir,'regrouped_data/sites_info.txt'),'rU') 
    reader = oldres.readlines()
    oldres.close()
    for l,line in enumerate(reader):
        site = line.split(',')[0].strip(' ')
        if l==0: line = line.strip('\n')+', gridcells\n'
        else: line = line.strip('\n') + ',%10i'%int(flux_gri[site]) + '\n'
        newres.write(line)
    newres.close()
    print '\nWe successfully created the input file with grid cell IDs:\n%s'%filename
    

#-------------------------------------------------------------------------------
# retrieve the necessary input data for all sites

    # settings of the connection
    user = "cgms12eu_select"
    password = "OnlySelect"
    tns = "EURDAS.WORLD"
    dsn = "oracle+cx_oracle://{user}:{pw}@{tns}".format(user=user,pw=password,tns=tns)
    engine = sa.create_engine(dsn)
    print engine

    # test the connection:
    try:
        connection = cx_Oracle.connect("cgms12eu_select/OnlySelect@eurdas.world")
    except cx_Oracle.DatabaseError:
        print '\nBEWARE!! The Oracle database is not responding. Probably, you are'
        print 'not using a computer wired within the Wageningen University network.'
        print '--> Get connected with ethernet cable before trying again!'
        sys.exit()

    for c,crop in enumerate(crops):
        crop_no = crop_nos[c]

        print '\nRetrieving input data for %s (CGMS id=%i)'%(crop,crop_no)
        # We add a timestamp at start of the retrieval
        start_timestamp = datetime.utcnow()
        
		# We retrieve the list of suitable soil types for the selected crop
		# species
        filename = os.path.join(CGMSdir, 'soildata_objects/',
                   'suitablesoilsobject_c%d.pickle'%(crop_no))
        if os.path.exists(filename):
            suitable_stu = pickle_load(open(filename,'rb'))
        else:
            from pcse.db.cgms11 import STU_Suitability
            suitable_stu = STU_Suitability(engine, crop_no)
            suitable_stu_list = []
            for item in suitable_stu:
                suitable_stu_list = suitable_stu_list + [item]
            suitable_stu = suitable_stu_list
            pickle_dump(suitable_stu,open(filename,'wb'))       
            print 'retrieving suitable soils for %s'%crop

        # WE LOOP OVER ALL YEARS:
        for y, year in enumerate(years): 
            print '\n######################## Year %i ##############'%year+\
            '##########\n'
        
            # if we do a serial iteration, we loop over the grid cells that 
            # contain arable land
            for grid in flux_gri.values():
                retrieve_CGMS_input(grid, year, crop_no, suitable_stu, engine)
        
        # We add a timestamp at end of the retrieval, to time the process
        end_timestamp = datetime.utcnow()
        print '\nDuration of the retrieval:', end_timestamp-start_timestamp

#===============================================================================
# function that will retrieve CGMS input data from Oracle database
def retrieve_CGMS_input(grid, year, crop_no, suitable_stu, engine, retrieve_weather=False):
#===============================================================================
# Temporarily add code directory to python path, to be able to import pcse
# modules
    sys.path.insert(0, codedir) 
#-------------------------------------------------------------------------------
    from pcse.exceptions import PCSEError 
    from pcse.db.cgms11 import TimerDataProvider, SoilDataIterator, \
                               CropDataProvider, STU_Suitability, \
                               SiteDataProvider, WeatherObsGridDataProvider
# if the retrieval does not raise an error, the crop was cultivated that year
    print '    - grid cell no %i'%grid
    try:
        # We retrieve the crop calendar (timerdata)
        filename = os.path.join(CGMSdir,
                   'timerdata_objects/%i/c%i/'%(year,crop_no),
                   'timerobject_g%d_c%d_y%d.pickle'%(grid, crop_no, year))
        if os.path.exists(filename):
            pass
        else:
            timerdata = TimerDataProvider(engine, grid, crop_no, year)
            pickle_dump(timerdata,open(filename,'wb'))    

        # If required by the user, we retrieve the weather data
        if retrieve_weather == True: 
            filename = os.path.join(CGMSdir, 'weather_objects/',
                       'weatherobject_g%d.pickle'%(grid))
            if os.path.exists(filename):
                pass
            else:
                weatherdata = WeatherObsGridDataProvider(engine, grid)
                weatherdata._dump(filename)

        # We retrieve the soil data (soil_iterator)
        filename = os.path.join(CGMSdir, 'soildata_objects/',
                   'soilobject_g%d.pickle'%(grid))
        if os.path.exists(filename):
            soil_iterator = pickle_load(open(filename,'rb'))
        else:
            soil_iterator = SoilDataIterator(engine, grid)
            pickle_dump(soil_iterator,open(filename,'wb'))       

        # We retrieve the crop variety info (crop_data)
        filename = os.path.join(CGMSdir,
                   'cropdata_objects/%i/c%i/'%(year,crop_no),
                   'cropobject_g%d_c%d_y%d.pickle'%(grid,crop_no,year))
        if os.path.exists(filename):
            pass
        else:
            cropdata = CropDataProvider(engine, grid, crop_no, year)
            pickle_dump(cropdata,open(filename,'wb'))     

        # WE LOOP OVER ALL SOIL TYPES LOCATED IN THE GRID CELL:
        for smu_no, area_smu, stu_no, percentage, soildata in soil_iterator:

            # NB: we remove all unsuitable soils from the iteration
            if (stu_no not in suitable_stu):
                pass
            else:
                print '        soil type no %i'%stu_no

                # We retrieve the site data (site management)
                if (str(grid)).startswith('1'):
                    dum = str(grid)[0:2]
                else:
                    dum = str(grid)[0]
                filename = os.path.join(CGMSdir,
                           'sitedata_objects/%i/c%i/grid_%s/'%(year,crop_no,dum),
                           'siteobject_g%d_c%d_y%d_s%d.pickle'%(grid, crop_no,
                                                                  year, stu_no))
                if os.path.exists(filename):
                    pass
                else:
                    sitedata = SiteDataProvider(engine,grid,crop_no,year,stu_no)
                    pickle_dump(sitedata,open(filename,'wb'))     

    # if an error is raised, the crop was not grown that year
    except PCSEError:
        print '        the crop was not grown that year in that grid cell'
    except Exception as e:
        print '        Unexpected error', e#sys.exc_info()[0]

    return None

#===============================================================================
# Function to open normal csv files
def open_csv(inpath, namefile, convert_to_float=False):
#===============================================================================

    from csv import reader as csv_reader

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
        if convert_to_float==True:
            converted_data.append(map(float,line))
        else:
            converted_data.append(line)
    data = np.array(converted_data)

    # creating one dictionnary and storing the float data in it
    dictnamelist= {}
    for j,varname in enumerate(headerow):
        dictnamelist[varname]=data[:,j]
    
    return dictnamelist

#===============================================================================
if __name__=='__main__':
    main()
#===============================================================================
