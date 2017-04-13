#!/bin/bash - 
#===============================================================================
#
#          FILE: data.mage.CORDEX-time-variability.sh
# 
#         USAGE: ./data.mage.CORDEX-time-variability.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tang (Tang), tangchao90908@sina.com
#  ORGANIZATION: KLA
#       CREATED: 02/02/17 16:40:47 RET
#      REVISION:  ---
#===============================================================================

set -o nounset               # Treat unset variables as an error
source /Users/ctang/Code/CORDEX_AFR_studies/cordex.function.sh

VAR='rsds'
YEAR1='1976 2005'
YEAR2='2041 2070'
#=================================================== 
copy data
#copydaily $VAR

#selyear $(eval echo $YEAR1 $YEAR2)

#selregion $VAR

#fldmean
#detrending 
#maskannual

#=================================================== calculate
