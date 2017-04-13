#!/bin/bash - 
#===============================================================================
#
#          FILE: function.sh
# 
#         USAGE: ./function.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tang (Tang), tangchao90908@sina.com
#  ORGANIZATION: KLA
#       CREATED: 02/02/17 16:41:10 RET
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

prefix="rsds_AFR-44"

function getdata4titan()
{
    rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/*ensmean* ./
    rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/*timstd.nc ./
    rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/*timmean.nc ./
}
#=================================================== cal

function printmeanchanges()
{
    echo printmeanchanges

    #           region_1 region_2 ... region_i ...
    #  model_1
    #  model_2
    #  ...
    #  model_m
    #  ...

    for m in $(cat original.cal)
    do
        for i in 1 2 3 4 5 6 7 8 9
        do
            meanchanges=$(ncks -v rsds $m.fldmean.daychanges.timmean.$i.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            mean_ref=$(ncks -v rsds $m.1976-2005.nc.$i.nc.fldmean.nc.timmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            echo $meanchanges $mean_ref | awk '{printf "%5.3f\t",$1/$4}'
        done
        echo ""
    done

        for i in 1 2 3 4 5 6 7 8 9
        do
            ensmeanchanges=$(ncks -v rsds rsds_AFR-44.fldmean.daychanges.timmean.ensmean.$i.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            ensmean_ref=$(ncks -v rsds rsds_AFR.1976-2005.fldmean.timmean.ensmean.$i.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            echo $ensmeanchanges $ensmean_ref | awk '{printf "%5.3f\t",$1/$4}'
        done
        echo ""
}

function printmeanref()
{

    echo printmeanref
    #           region_1 region_2 ... region_i ...
    #  model_1
    #  model_2
    #  ...
    #  model_m
    #  ...

    echo ------------ ref mean ---------

    for m in $(cat model)
    do
        for i in 1 2 3 4 5 6 7 8 9
        do
            mean_ref=$(ncks -v rsds ${prefix}_$m.hist.day.1970-1999.$i.fldmean.timmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            #mean_future=$(ncks -v rsds ${prefix}_$m.rcp85.day.2070-2099.$i.timmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')

            echo $mean_ref | awk '{printf "%5.3f\t",$1}'
        done
        echo ""
    done

        for i in 1 2 3 4 5 6 7 8 9
        do
            ensmean_ref=$(ncks -v rsds rsds_AFR-44.hist.day.1970-1999.$i.fldmean.timmean.ensmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            echo $ensmean_ref | awk '{printf "%5.3f\t",$1}'
        done
        echo ""
}
function printmeanfuture()
{

    echo printmeanfuture
    #           region_1 region_2 ... region_i ...
    #  model_1
    #  model_2
    #  ...
    #  model_m
    #  ...

    echo ------------ future mean ---------

    for m in $(cat model)
    do
        for i in 1 2 3 4 5 6 7 8 9
        do
            mean_ref=$(ncks -v rsds ${prefix}_$m.rcp85.day.2070-2099.$i.fldmean.timmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            #mean_future=$(ncks -v rsds ${prefix}_$m.rcp85.day.2070-2099.$i.timmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')

            echo $mean_ref | awk '{printf "%5.3f\t",$1}'
        done
        echo ""
    done

        for i in 1 2 3 4 5 6 7 8 9
        do
            ensmean_ref=$(ncks -v rsds rsds_AFR-44.rcp85.day.2070-2099.$i.fldmean.timmean.ensmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            echo $ensmean_ref | awk '{printf "%5.3f\t",$1}'
        done
        echo ""
}
function printstdref()
{

    echo printstdref
    #           region_1 region_2 ... region_i ...
    #  model_1
    #  model_2
    #  ...
    #  model_m
    #  ...

    echo ------------ ref mean ---------
for t in daymean monmean yearmean
do

    echo ------------ $t -----------
    for m in $(cat model)
    do
        for i in 1 2 3 4 5 6 7 8 9
        do
            mean_ref=$(ncks -v rsds ${prefix}_$m.hist.day.1970-1999.$i.fldmean.detrend.nc.$t.mask*.timstd.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            #mean_future=$(ncks -v rsds ${prefix}_$m.rcp85.day.2070-2099.$i.timmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')

            echo $mean_ref | awk '{printf "%5.3f,",$1}'
        done
        echo "],\\"
    done

        for i in 1 2 3 4 5 6 7 8 9
        do
            ensmean_ref=$(ncks -v rsds rsds_AFR-44.hist.day.1970-1999.$i.fldmean.detrend.nc.$t.mask*.timstd.ensmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            echo $ensmean_ref | awk '{printf "%5.3f,",$1}'
        done
        echo ""
done
}
function printstdfuture()
{
    echo printstdfuture
    #           region_1 region_2 ... region_i ...
    #  model_1
    #  model_2
    #  ...
    #  model_m
    #  ...

    echo ------------ ref mean ---------
for t in daymean monmean yearmean
do

    echo ------------ $t -----------
    for m in $(cat model)
    do
        for i in 1 2 3 4 5 6 7 8 9
        do
            mean_ref=$(ncks -v rsds ${prefix}_$m.rcp85.day.2070-2099.$i.fldmean.detrend.nc.$t.mask*.timstd.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            #mean_future=$(ncks -v rsds ${prefix}_$m.rcp85.day.2070-2099.$i.timmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')

            echo $mean_ref | awk '{printf "%5.3f,",$1}'
        done
        echo "],\\"
    done

        for i in 1 2 3 4 5 6 7 8 9
        do
            ensmean_ref=$(ncks -v rsds rsds_AFR-44.rcp85.day.2070-2099.$i.fldmean.detrend.nc.$t.mask*.timstd.ensmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            echo $ensmean_ref | awk '{printf "%5.3f,",$1}'
        done
        echo ""
done
}
printmeanref
printmeanfuture
printstdref
printstdfuture
