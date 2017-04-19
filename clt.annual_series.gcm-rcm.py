#!/usr/bin/env python

########################################
#Globale Karte fuer tests
# from Rabea Amther
########################################
# http://gfesuite.noaa.gov/developer/netCDFPythonInterface.html

import math
import pdb
import numpy as np
import pandas as pd
import pylab as pl
import datetime
import Scientific.IO.NetCDF as IO
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.lines as lines
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.dates import YearLocator,MonthLocator,DateFormatter,drange
from matplotlib.ticker import AutoMinorLocator
import textwrap

# to load my functions
import sys 
sys.path.append('/Users/ctang/Code/Python/')
import ctang

pl.close('all')

########################## for CORDEX charactors
DIR='/Users/ctang/Code/CORDEX_AFR_studies/data/regional_annual_series/'
N_region = 7
N_model = 21
N_model_GCM = 8
VAR = 'clt'
YEAR1=1970
YEAR2=1999
YEAR3=2006
YEAR4=2099

MODEL=(\
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
	'MPI-M-MPI-ESM-LR_MPI-CSC-REMO2009_v1')

GCM=(\
    'clt_Amon_CNRM-CM5_historical-rcp85_r1i1p1',\
    'clt_Amon_CanESM2_historical-rcp85_r1i1p1',\
    'clt_Amon_GFDL-ESM2M_historical-rcp85_r1i1p1',\
    'clt_Amon_HadGEM2-ES_historical-rcp85_r1i1p1',\
    'clt_Amon_IPSL-CM5A-LR_historical-rcp85_r1i1p1',\
    'clt_Amon_MIROC5_historical-rcp85_r1i1p1',\
    'clt_Amon_MPI-ESM-LR_historical-rcp85_r1i1p1',\
    'clt_Amon_NorESM1-M_historical-rcp85_r1i1p1')
#---------------------------------------------------  data
# 11 * 9 table: generated by time.variability.function.sh

    #               region_1 region_2 ... region_N_region
    #  model_1
    #  model_2
    #  ...
    #  model_N_model
    #  ens
#--------------------------------------------------- 
# define the time series:
TIME=ctang.get_netcdf_time(\
        DIR+VAR+'_AFR-44_'+MODEL[0]+'.hist_rcp85.day.1970-2099.1'+\
        '.fldmean.yearmean.nc')
TIME=np.array(range(YEAR1,YEAR4+1,1))

#each map is mean_ref[0,0,:,:].shape
Anomaly = np.zeros((N_model, N_region, YEAR4-YEAR1+1))
Anomaly_GCM = np.zeros((N_model_GCM, N_region, YEAR4-YEAR1+1))
Changes = np.zeros((N_model, N_region, YEAR4-YEAR1+1))
Changes_GCM = np.zeros((N_model_GCM, N_region, YEAR4-YEAR1+1))
yearmean_future = np.zeros((N_model, N_region, YEAR4-YEAR1+1))
yearmean_future_GCM = np.zeros((N_model_GCM, N_region, YEAR4-YEAR1+1))

# in 3D: (N_model,N_region,130years)
for r in range(N_region):
    yearmean_future[:,r,:]=np.array([[ctang.read_time_netcdf(VAR,\
        DIR+VAR+'_AFR-44_'+MODEL[i]+'.hist_rcp85.day.1970-2099.'+str(r+1)+\
        '.fldmean.yearmean.nc')] for i in range(N_model)])[:,0,:]
    yearmean_future_GCM[:,r,:]=np.array([[ctang.read_time_netcdf(VAR,\
        DIR+GCM[i]+'.1970-2099.'+str(r+1)+'.fldmean.yearmean.nc')] \
        for i in range(N_model_GCM)])[:,0,:]
    Anomaly[:,r,:]=np.array([[ctang.read_time_netcdf(VAR,\
        DIR+VAR+'_AFR-44_'+MODEL[i]+'.hist_rcp85.day.1970-2099.'+str(r+1)+\
        '.fldmean.monmean.anomaly.yearmean.nc')] for i in range(N_model)])[:,0,:]
    Anomaly_GCM[:,r,:]=np.array([[ctang.read_time_netcdf(VAR,\
        DIR+GCM[i]+'.1970-2099.'+str(r+1)+'.fldmean.monmean.anomaly.yearmean.nc')]\
        for i in range(N_model_GCM)])[:,0,:]


# timmean_ref as ref of anomalies
timmean_ref=np.mean(yearmean_future[:,:,0:YEAR2-YEAR1],axis=2)
timmean_ref_GCM=np.mean(yearmean_future_GCM[:,:,0:YEAR2-YEAR1],axis=2)

# calculate changes:
for r in range(N_region):
    for m in range(N_model):
        Changes[m,r,:] = np.array([(t - timmean_ref[m,r])*100/timmean_ref[m,r] \
                for t in yearmean_future[m,r,:]])
    for m in range(N_model_GCM):
        Changes_GCM[m,r,:] = np.array( [ (t - timmean_ref_GCM[m,r])*100/timmean_ref_GCM[m,r]\
                for t in yearmean_future_GCM[m,r,:]])

# nomarlize by the timmean_ref
for r in range(N_region):
    for m in range(N_model):
        Anomaly[m,r,:]=np.array(\
            [ (t)*100/timmean_ref[m,r] for t in Anomaly[m,r,:]])
    for m in range(N_model_GCM):
        Anomaly_GCM[m,r,:]=np.array(\
            [ (t)*100/timmean_ref_GCM[m,r] for t in Anomaly_GCM[m,r,:]])

#=================================================== test
# fig1, ax = plt.subplots()
# for r in range(N_region):
    # for m in range(N_model):
        # plt.plot(range(2099-1969),Anomaly[m,4,:])
        # plt.plot(range(2099-1969),np.mean(Anomaly[:,4,:],axis=0),color='b')
        # # plt.plot(range(2099-1969),range(2099-1969))

# fig2, ax = plt.subplots()
# for r in range(N_region):
    # for m in range(N_model_GCM):
        # plt.plot(range(2099-1969),Anomaly_GCM[m,4,:])
        # plt.plot(range(2099-1969),np.mean(Anomaly_GCM[:,4,:],axis=0),color='b')
        # # plt.plot(range(2099-1969),range(2099-1969))
# # plt.show()

#=================================================== test



# Ensmean of timstd in ref estimation of noise of climate change signal
timstd_ref=np.std(Anomaly[:,:,0:YEAR2-YEAR1],axis=2)
timstd_ref_GCM=np.std(Anomaly_GCM[:,:,0:YEAR2-YEAR1],axis=2)
Ensmean_timstd_ref=np.mean(timstd_ref,axis=0)
Ensmean_timstd_ref_GCM=np.mean(timstd_ref_GCM,axis=0)

# running mean of each model in each region:
Anomaly2 = np.zeros((N_model, N_region, 101))
Anomaly_GCM2 = np.zeros((N_model_GCM, N_region, 101))

for r in range(N_region):
    for m in range(N_model):
        Anomaly2[m,r,:]=np.array(ctang.running_mean3(Anomaly[m,r,:],30))
    for m in range(N_model_GCM):
        Anomaly_GCM2[m,r,:]=np.array(ctang.running_mean3(Anomaly_GCM[m,r,:],30))

print Anomaly_GCM2.shape
#=================================================== test
# fig1, ax = plt.subplots()
# for r in range(N_region):
    # for m in range(N_model):
        # plt.plot(range(101),Anomaly2[m,4,:])
        # plt.plot(range(101),np.mean(Anomaly2[:,4,:],axis=0),color='b')
        # plt.plot(range(2099-1969),range(2099-1969))

# fig2, ax = plt.subplots()
# for r in range(N_region):
    # for m in range(N_model_GCM):
        # plt.plot(range(101),Anomaly_GCM2[m,4,:])
        # plt.plot(range(101),np.mean(Anomaly_GCM2[:,4,:],axis=0),color='b')
        # plt.plot(range(2099-1969),range(2099-1969))

#=================================================== test
# get ensmean & ensstd of anomaly
ensmean_anomaly=np.mean(Anomaly2,axis=0)
ensmean_anomaly_GCM=np.mean(Anomaly_GCM2,axis=0)
ensstd_anomaly=np.std(Anomaly2,axis=0)
ensstd_anomaly_GCM=np.std(Anomaly_GCM2,axis=0)

#=================================================== plot
Title='annual time series (30-year running mean) of clt along the 21st century under RCP8.5'


fig, axes = plt.subplots(nrows=3, ncols=3,figsize=(16, 9),facecolor='w', edgecolor='k')
fig.subplots_adjust(left=0.1,right=0.95, top=0.9,bottom=0.11,hspace=0.5,wspace=0.45)
axes = axes.flatten() # reshape plot (3*3) to 9*1

# set ticks and tick labels
for k in range(N_region):
    axes[k].set_xlim((1985, 2085))

    axes[k].set_xlabel('Year',fontsize=10)
    axes[k].tick_params(direction='out',length=6,width=2,labelsize=10)

    axes[k].set_title('Region '+str(k+1),fontsize=10)
    axes[k].set_ylabel('rsds anomalies in %', fontsize=10)
    
    axes[k].set_axisbelow(True)
    axes[k].yaxis.grid(color='gray', linestyle='dashed')
    axes[k].xaxis.grid(color='gray', linestyle='dashed')
    minorLocator = AutoMinorLocator()
    axes[k].xaxis.set_minor_locator(minorLocator)
    axes[k].yaxis.set_minor_locator(minorLocator)

    # Hide the right and top spines
    axes[k].spines['right'].set_visible(False)
    axes[k].spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    axes[k].yaxis.set_ticks_position('left')
    axes[k].xaxis.set_ticks_position('bottom')
    axes[k].axhline(y=0, xmin=0, xmax=N_region+1,color='black',linewidth=1)

#=================================================== plot every model, RCMs and GCMs
    
    # for m in range(5):
        # axes[k].plot(TIME[15:116], Anomaly2[m,k,:], '-',linewidth=0.5)
        #legend = axes[5].legend(loc='upper center', shadow=True)

    #for m in range(N_model_GCM):
        #axes[k].plot(TIME[15:116], ctang.running_mean(Anomaly_GCM[m,k,:],5),\
                #'-',linewidth=0.5)
        #legend = axes[5].legend(loc='upper center', shadow=True)

    # for m in range(N_model):
        # axes[k].plot(TIME, Changes[m,k,:], '-',linewidth=0.5)

    # for m in range(N_model_GCM):
        # axes[k].plot(TIME, Anomaly_GCM[1,k,:], '-',linewidth=2.5,color='green')

    # for m in range(N_model_GCM):
        # axes[k].plot(TIME, Changes_GCM[1,k,:], '-',linewidth=0.5,color='red')

#=================================================== plot ensmean of timstd ref

    axes[k].errorbar(1990,0,yerr=Ensmean_timstd_ref[k], color='red', ls='-',linewidth=2)
    axes[k].plot(1990,Ensmean_timstd_ref[k],'_',color='red',linewidth=8)
    axes[k].plot(1990,-1*Ensmean_timstd_ref[k],'_',color='red',linewidth=8)

#=================================================== plot runmean ensmean anomalies

    # axes[k].plot(TIME[15:116],ensmean_anomaly[k],\
    axes[k].plot(TIME[15:116],ensmean_anomaly[k],\
        '-', label='21 RCMs mean', color='blue', linewidth=2,zorder=2)

    # axes[k].plot(TIME[15:116],ensmean_anomaly_GCM[k],\
    axes[k].plot(TIME[15:116],ensmean_anomaly_GCM[k],\
        '-', label='8 GCMs mean', color='black', linewidth=2,zorder=2)

#=================================================== plot ensstd of anomalies

    # rcm
    axes[k].plot(TIME[15:116],np.subtract(ensmean_anomaly[k], ensstd_anomaly[k]),\
        '-',  color='blue', linewidth=0.1,zorder=1)
    axes[k].plot(TIME[15:116],np.add(ensmean_anomaly[k], ensstd_anomaly[k]),\
        '-',  color='blue', linewidth=0.1,zorder=1)
    axes[k].fill_between(TIME[15:116],\
        np.subtract(ensmean_anomaly[k], ensstd_anomaly[k]),\
        np.add(ensmean_anomaly[k], ensstd_anomaly[k]),\
        color='blue',alpha=0.3,zorder=3)

    # gcm
    axes[k].plot(TIME[15:116],np.subtract(ensmean_anomaly_GCM[k],ensstd_anomaly_GCM[k]),\
        '-', color='black', linewidth=0.1,zorder=1)
    axes[k].plot(TIME[15:116],np.add(ensmean_anomaly_GCM[k], ensstd_anomaly_GCM[k]),\
        '-', color='black', linewidth=0.1,zorder=9)
    axes[k].fill_between(TIME[15:116], \
            np.subtract(ensmean_anomaly_GCM[k],ensstd_anomaly_GCM[k]),\
            np.add(ensmean_anomaly_GCM[k],ensstd_anomaly_GCM[k]),\
            color='black',alpha=0.3,zorder=5)
    axes[k].legend(loc='lower left', shadow=False)

#=================================================== add legend
#legend = ax.legend(loc='upper center', shadow=True)

for kk in [7,8]:
    ctang.empty_plot(axes[kk])

plt.suptitle(Title)

#===================================================  end of  plot

plt.savefig('clt.annual_series.gcm-rcm.eps', format='eps')
plt.savefig('clt.annual_series.gcm-rcm.png')
plt.show()
quit()

