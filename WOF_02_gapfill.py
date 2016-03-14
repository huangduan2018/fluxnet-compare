#!/usr/bin/env python

import sys, os, rc
import numpy as np
from cPickle import load as pickle_load
from cPickle import dump as pickle_dump

#===============================================================================
# This script does some standard analysis on FluxNet sites
def main():
#===============================================================================
    global inputdir, outputdir, optimidir
#-------------------------------------------------------------------------------
# ================================= USER INPUT =================================

# read the settings from the rc file
    rcdict     = rc.read('settings.rc')

#===============================================================================
# extract the needed information from the rc file
    sites      = [s.strip(' ') for s in rcdict['sites'].split(',')]
    #NUTS_reg   = [s.strip(' ') for s in rcdict['NUTS_reg'].split(',')]
    crops      = [s.strip(' ') for s in rcdict['crops'].split(',')]
    crop_nos   = [int(s.strip(' ')) for s in rcdict['crop_nos'].split(',')]
    years      = [int(s.strip(' ')) for s in rcdict['years'].split(',')]

    # directory paths
    outputdir  = rcdict['outputdir']
    inputdir   = rcdict['inputdir']

#-------------------------------------------------------------------------------
    # get the list of NUTS 2 region names associated to the list of FluxNet sites
    from WOF_00_retrieve_input_data import open_csv
    sitdict = open_csv(inputdir, 'sites_info2.csv', convert_to_float=False)
    NUTS_reg  = sitdict['NUTS_reg']
#-------------------------------------------------------------------------------
# list the old gapfilled files to remove, and remove them all

    for s,site in enumerate(sites):

        for c,crop_name in enumerate(crops):
            crop_no = crop_nos[c]

            for year in years:
                optimidir = os.path.join(outputdir,'fgap/%i/c%i/'%(year,crop_no))

                files2remove = [f for f in os.listdir(optimidir) if '_gapfilled' in f]
                for f in files2remove:
                    os.remove(os.path.join(optimidir,f))

#-------------------------------------------------------------------------------
# gap fill

    for s,site in enumerate(sites):
        NUTS_no = NUTS_reg[s]

        for c,crop_name in enumerate(crops):
            crop_no = crop_nos[c]

            for year in years:
                # create output folder if it doesn't already exists
                optimidir = os.path.join(outputdir,'fgap/%i/c%i/'%(year,crop_no))

                # detect if there is this year needs to be gapfilled
                f2gapfill = [f for f in os.listdir(optimidir) if ('_tobegapfilled' 
                             in f) and (NUTS_no in f)]
                if len(f2gapfill)==0:
                    continue

                print '\nWe gap fill:', site, NUTS_no, year, crop_name

                # GAP-FILLING YLDGAPF for NUTS2 level:
                prevyear = os.path.join(optimidir.replace('%04d'%year, 
                         '%04d'% (year-1 )),'fgap_%s_optimized.pickle'%NUTS_no)
                nextyear = os.path.join(optimidir.replace('%04d'%year, 
                         '%04d'% (year+1 )),'fgap_%s_optimized.pickle'%NUTS_no)
                availfiles = []
                availyears = []
                for yr in range(1995,2020):
                    searchyear = os.path.join(optimidir.replace('%04d'%year, 
                               '%04d'% yr ),'fgap_%s_optimized.pickle'%NUTS_no)
                    if os.path.exists(searchyear) :
                        availfiles.append(searchyear)
                        availyears.append(yr)
                print "%d years found for gap filling:" %len(availfiles), availyears

                # Use average from y-1 and y+1
                if prevyear in availfiles and nextyear in availfiles:
                    optimi_info = pickle_load(open(prevyear,'rb'))
                    ygf_prev    = optimi_info[2]
                    optimi_info = pickle_load(open(nextyear,'rb'))
                    ygf_next    = optimi_info[2]
                    ygf         = (ygf_prev+ygf_next)/2.0  # simply average
                    opt_code    = 'gapfilled02'
                    shortlist_cells = optimi_info[3]

                # Use previous year value
                elif prevyear in availfiles: 
                    optimi_info = pickle_load(open(prevyear,'rb'))
                    ygf         = optimi_info[2]
                    opt_code    = 'gapfilled03a'
                    shortlist_cells = optimi_info[3]
                    print shortlist_cells

                # Use next year value
                elif nextyear in availfiles:
                    optimi_info = pickle_load(open(nextyear,'rb'))
                    ygf         = optimi_info[2]
                    opt_code    = 'gapfilled03b'
                    shortlist_cells = optimi_info[3]

                # Use climatological average from other years if nyear > 2
                elif len(availfiles) > 2:  
                    ygf = 0.0
                    for filename in availfiles:
                        optimi_info = pickle_load(open(filename,'rb'))
                        ygf += optimi_info[2]
                    ygf = ygf/len(availfiles)
                    opt_code = 'gapfilled04'
                    shortlist_cells = optimi_info[3]
                # Use upper NUTS level optimum (NUTS1, or NUTS0 at worst)
                else:
                    try:
                        nuts1file = os.path.join(optimidir, 
                                     'fgap_%s_optimized.pickle'%NUTS_no[0:3])
                        data = pickle_load(open(nuts1file,'rb'))
                        ygf  = data[2]
                        opt_code = 'gapfilled05a'
                        shortlist_cells = data[3]
                    except IOError:
                        try:
                            nuts0file = os.path.join(optimidir, 
                                         'fgap_%s_optimized.pickle'%NUTS_no[0:2])
                            data = pickle_load(open(nuts0file,'rb'))
                            ygf  = data[2]
                            opt_code = 'gapfilled05b'
                            shortlist_cells = data[3]
                # Use default value if all previous methods fail
                        except IOError:
                            ygf  = 0.8
                            opt_code = 'gapfilled06'
                            shortlist_cells = []

                print "Using ygf of %5.2f and code of %s"%(ygf, opt_code)
                print "created file fgap_%s_%s.pickle"%(NUTS_no, opt_code)+\
                      " in folder %s"%optimidir
                currentyear = os.path.join(optimidir,'fgap_%s_%s.pickle'%(
                                                          NUTS_no, opt_code))
                pickle_dump([NUTS_no,opt_code,ygf,shortlist_cells],
                                                      open(currentyear,'wb'))



#===============================================================================
if __name__=='__main__':
    main()
#===============================================================================
