#!/bin/sh

# source the Modules package, and load it
source /usr/local/Modules/3.2.8/init/sh
module load intel/11.1.084/compiler
module load intel/11.1.084/netcdf/3.6.3

site[1]='BE-Lon'
site[2]='DE-Kli'
site[3]='FR-Gri'
site[4]='IT-BCi'
site[5]='NL-Dij'
site[6]='NL-Lan'

lon[1]=4.74
lon[2]=13.52
lon[3]=1.95
lon[4]=15.01 # we had to change it from 14.95 to 15.01 to pick a land - not ocean - grid cell
lon[5]=5.65
lon[6]=4.90

lat[1]=50.55
lat[2]=50.89
lat[3]=48.84
lat[4]=40.52
lat[5]=51.99
lat[6]=51.95

for i in 2 3 4 5 6
do
    ### make the default run (nb 0) ###
    echo ''
    echo ${site[$i]}

    # erase old initial run
    outdir=output_runs/search_steady/${site[$i]}_run0
    rm -r $outdir
    mkdir $outdir
    # start from the default files
    cp namel_sibdrv.default namel_sibdrv
    cp /Storage/SiBCASA/michiel/Siberia/default/sib_requib_00.nc input_files/sib_requib_0.nc

    # overwrite some settings in the namelist
    # output dir path, restart file path
    sed -i '' "9s#output_runs#${outdir}#1" namel_sibdrv
    sed -i '' "10s#michiel/Siberia/default/sib_requib_00.nc#marie/input_files/sib_requib_0.nc#1" namel_sibdrv
    # lon and lat of site location
    sed -i '' "25s/15.01/${lon[$i]}/1" namel_sibdrv
    sed -i '' "26s/15.01/${lon[$i]}/1" namel_sibdrv
    sed -i '' "27s/40.52/${lat[$i]}/1" namel_sibdrv
    sed -i '' "28s/40.52/${lat[$i]}/1" namel_sibdrv
    # start and end year: the decade before 2000-2014
    sed -i '' "39s/2004/1999/1" namel_sibdrv
    sed -i '' "41s/2004/1999/1" namel_sibdrv

    # run SiBCASA
    /Volumes/Storage/SiBCASA/marie/SiBD3

    # prepare to restart the following runs from here
    cp $outdir/sib_requib.nc input_files/sib_requib_1.nc

    ### make the following tryout runs (nb 1-10) ###
    for runnb in 1 2 3 4 5 6 7 8 9 10
    do
        echo '... restart nb: '$runnb

        # create an output directory for the run if it doesn't exist yet
        outdir=output_runs/search_steady/${site[$i]}_run${runnb}
        mkdir -p $outdir

        # overwrite some settings in the namelist: output dir path, restart file path
        oldnb=$((runnb-1))
        sed -i '' "9s#output_runs/search_steady/${site[$i]}_run${oldnb}#${outdir}#1" namel_sibdrv
        sed -i '' "10s#sib_requib_${oldnb}.nc#sib_requib_${runnb}.nc#1" namel_sibdrv
 
        # run SiBCASA
        /Volumes/Storage/SiBCASA/marie/SiBD3
 
        # prepare to restart the following runs from here
        newnb=$((runnb+1))
        cp $outdir/sib_requib.nc input_files/sib_requib_${newnb}.nc

    done
done

echo ''
echo 'Success! All repetitions of spin-up runs are done.'
echo 'Next: analyze steady-state.'
