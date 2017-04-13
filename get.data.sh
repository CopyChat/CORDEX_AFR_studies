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

    # for annual_series
    for v in rsds clt
    do
        rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CMIP5/monthly/$v/rcp8.5/*2099*fldmean.yearmean.nc ./data
        rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CMIP5/monthly/$v/rcp8.5/*anomaly.yearmean.nc ./data/

        rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/$v/*anomaly.yearmean.nc ./data
        rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/$v/*2099*fld*yearmean.nc ./data
    done

    # for time variability
    for v in rsds clt
    do
        rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/$v/*.?.*mask*std.nc ./data/
    done

    #rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/rsds/*SA.timmean.nc ./data/

    #rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/rsds/*99.?.timmean.nc ./data/

    #rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/rsds/*fldmean.timmean.nc ./data/

    #rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/rsds/*SA.ymonmean.nc ./data/
    #rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/rsds/*ensstd.nc ./data/
    #rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/rsds/*ensmean.nc ./data/
    #rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/rsds/*timstd.nc ./data/
    #rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/rsds/*fldmean.yearmean.nc ./data/

#    rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/tas/*SA.timmean.nc ./data/
    #rsync -arSzPH ctang@ccur.univ-reunion.fr:
    
    #rsync -arSzPH ctang@ccur.univ-reunion.fr:/worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/CORDEX-time-variability/sfcWind/*SA.timmean.nc ./data/
}

getdata4titan

