#!/usr/bin/env python
"""
========
Ctang, A bar plot of time variability changes projection 
        from CORDEX AFR-44, in Southern Africa
        Data was restored on titan
========
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import pandas as pd

import textwrap
import datetime
import ctang

from mpl_toolkits.basemap import Basemap 

from matplotlib.ticker import AutoMinorLocator
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange

import math
from scipy import stats
import subprocess

plt.close('all')

#=================================================== pre-defined
VAR='ssrd'
CORDEX_VAR='rsds'
#=================================================== 
DIR='/Users/ctang/Code/CORDEX_AFR_studies/data/cor.rsds_clt.fldcor/'
OBS='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/'

# input file
DJF='ERA_In.ssrd.clt.seasonal.fldcor.DJF.nc'
JJA='ERA_In.ssrd.clt.seasonal.fldcor.JJA.nc'
COR='ERA_In.ssrd.clt.seasonal.fldcor.nc'

#=================================================== read
# read time series
print(OBS+DJF)
TIME_OBS=ctang.get_netcdf_time(OBS+DJF)

# read correlation
COR_OBS_DJF=ctang.read_time_netcdf(VAR,OBS+DJF)
COR_OBS_JJA=ctang.read_time_netcdf(VAR,OBS+JJA)

CORDEX=(\
        'CCCma-CanESM2_SMHI-RCA4_v1',\
        'CNRM-CERFACS-CNRM-CM5_SMHI-RCA4_v1',\
        'CSIRO-QCCCE-CSIRO-Mk3-6-0_SMHI-RCA4_v1',\
        'ICHEC-EC-EARTH_SMHI-RCA4_v1',\
        'IPSL-IPSL-CM5A-MR_SMHI-RCA4_v1',\
        'MIROC-MIROC5_SMHI-RCA4_v1',\
        'MOHC-HadGEM2-ES_SMHI-RCA4_v1',\
        'MPI-M-MPI-ESM-LR_SMHI-RCA4_v1',\
        'NCC-NorESM1-M_SMHI-RCA4_v1',\
        'NOAA-GFDL-GFDL-ESM2M_SMHI-RCA4_v1',\
        \
        'CNRM-CERFACS-CNRM-CM5_CLMcom-CCLM4-8-17_v1',\
        'ICHEC-EC-EARTH_CLMcom-CCLM4-8-17_v1',\
        'MOHC-HadGEM2-ES_CLMcom-CCLM4-8-17_v1',\
        'MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17_v1',\
        'ICHEC-EC-EARTH_DMI-HIRHAM5_v2',\
        'NCC-NorESM1-M_DMI-HIRHAM5_v1',\
        'ICHEC-EC-EARTH_KNMI-RACMO22T_v1',\
        'MOHC-HadGEM2-ES_KNMI-RACMO22T_v2',\
        'ICHEC-EC-EARTH_MPI-CSC-REMO2009_v1',\
        'IPSL-IPSL-CM5A-LR_GERICS-REMO2009_v1',\
        'MPI-M-MPI-ESM-LR_MPI-CSC-REMO2009_v1',\
        )

N_model=len(CORDEX)

TIME_CORDEX=ctang.get_netcdf_time(\
        DIR+'AFR-44_'+CORDEX[2]+'.rsds.clt.seasonal.fldcor.DJF.nc')

COR_CORDEX_DJF=np.array([ctang.read_time_netcdf(CORDEX_VAR,\
        DIR+'AFR-44_'+CORDEX[i]+'.rsds.clt.seasonal.fldcor.DJF.nc')\
        for i in range(N_model)])

COR_CORDEX_JJA=np.array([ctang.read_time_netcdf(CORDEX_VAR,\
        DIR+'AFR-44_'+CORDEX[i]+'.rsds.clt.seasonal.fldcor.JJA.nc')\
        for i in range(N_model)])

print COR_CORDEX_JJA.shape

#=================================================== 
# define subplots
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,7),\
        facecolor='w', edgecolor='k') # figsize=(w,h)

# set limits
# ax.set_xlim( dates[0], dates[-1] )

# set grid
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.xaxis.grid(color='gray', linestyle='dashed')

# set title
ax.set_xlabel('time', fontsize=12)
ax.set_ylabel('Correlation of CLT vs SSR', fontsize=12)

# big title
Title='fldcor of seasonal series CLT vs SSR over SA between 1971-2099'
plt.suptitle(Title)

# The hour locator takes the hour or sequence of hours you want to
# tick, not the base multiple

ax.xaxis.set_major_formatter( DateFormatter('%Y') )
ax.fmt_xdata = DateFormatter('%Y')
fig.autofmt_xdate()

#=================================================== plot


# cordex in JJA
for m in range(N_model):
    ax.plot(TIME_CORDEX,COR_CORDEX_JJA[m],label='CORDEX JJA' if m == 0 else "",linestyle='--',color='r',alpha=0.5)

# cordex in DJF
    ax.plot(TIME_CORDEX,COR_CORDEX_DJF[m],label='CORDEX DJF' if m == 0 else "",linestyle='--',color='b',alpha=0.5)
#=================================================== multimodel mean
ax.plot(TIME_CORDEX,np.mean(COR_CORDEX_DJF,axis=0),\
        label='CORDEX mean DJF' if m == 0 else "",\
        linestyle='-',linewidth=6,color='b')
ax.plot(TIME_CORDEX,np.mean(COR_CORDEX_JJA,axis=0),\
        label='CORDEX mean DJF' if m == 0 else "",
        linestyle='-',linewidth=6,color='r')

# plot OBS
ax.plot(TIME_OBS,COR_OBS_DJF,label='ERA_In DJF',color='b',linewidth=2)
ax.plot(TIME_OBS,COR_OBS_JJA,label='ERA_In JJA',color='r',linewidth=2)
ax.legend()

#=================================================== output

# save image
Out_Image='seasonal.series.fldcor'

plt.savefig(Out_Image+'.eps',format='eps')
plt.savefig(Out_Image+'.png')

plt.show()
