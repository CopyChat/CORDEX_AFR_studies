#!/usr/bin/env python
"""
========
Ctang, A bar plot of time variability changes projection 
        from CORDEX AFR-44, in Southern Africa
        Data was restored on titan
========
"""
import math
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.ticker import AutoMinorLocator


# to load my functions
import sys 
sys.path.append('/Users/ctang/Code/Python/')
import ctang

#=================================================== Definitions

#---------------------------------------------------  data
# 11 * 9 table: generated by time.variability.function.sh

    #               region_1 region_2 ... region_N_region
    #  model_1
    #  model_2
    #  ...
    #  model_N_model
    #  ens
#--------------------------------------------------- 

#=================================================== Definitions
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/time_variability/'
N_region = 7
N_model = 21
N_model_GCM = 10
VAR ='rsds'

N_plot = 3 # mean, monthly variability, annual variability

GCM_Model=(\
        'CNRM-CM5',\
        'CSIRO-Mk3-6-0',\
        'CanESM2',\
        'GFDL-ESM2M',\
        'HadGEM2-ES',\
        'IPSL-CM5A-LR',\
        'IPSL-CM5A-MR',\
        'MIROC5',\
        'MPI-ESM-LR',\
        'NorESM1-M')

RCM_Model=(\
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

#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

histfix='.hist.day.1970-1999.'
rcp85fix = '.rcp85.day.2070-2099.'

histfix_GCM ='_historical-rcp85_r1i1p1.1970-1999.'
rcp85fix_GCM = '_historical-rcp85_r1i1p1.2070-2099.'

filefix = (\
        '.fldmean.timmean.nc',\
        # '.fldmean.detrend.nc.daymean.maskannual.timstd.nc',\
        '.fldmean.detrend.nc.monmean.maskannual.timstd.nc',\
        '.fldmean.detrend.nc.yearmean.masknoooon.timstd.nc')

filefix_GCM = (\
        '.fldmean.timmean.nc',\
        '.fldmean.detrend.nc.monmean.maskannual.timstd.nc',\
        '.fldmean.detrend.nc.yearmean.masknoooon.timstd.nc')

Ref = np.zeros(( N_plot,N_model, N_region))
t_value = np.zeros(( N_model, N_region))
Future = np.zeros((N_plot, N_model, N_region))

Ref_GCM = np.zeros(( N_plot,N_model_GCM, N_region))
t_value_GCM = np.zeros(( N_model_GCM, N_region))
Future_GCM = np.zeros((N_plot, N_model_GCM, N_region))

# Read timmean for in each region, used for significance calculation
# the significance of mean changes over each region, this map of significance
# will be used for the time variability
# the loaded data will be in 4D map: 4Var, N_model, N_region, lat,lon

# print Data+VAR+'_AFR-44_'+RCM_Model[0]+histfix+str(1+1)+filefix[0]

# print Data+VAR+'_Amon_'+GCM_Model[0]+rcp85fix_GCM+str(3+1)+\
                # filefix_GCM[1]
#===================================================  reading functions
# Read for CORDEX 
def reading(Ref,Future):
    print("reading RCMs data ...")
    for v in range(N_plot):
        for i in range(N_model):
            Ref[v,i]=np.array([ctang.read_lonlatmap_netcdf(VAR,\
                Data+VAR+'_AFR-44_'+RCM_Model[i]+histfix+str(k+1)+filefix[v])[0,0]\
                for k in range(N_region)])
            Future[v,i]=np.array([ctang.read_lonlatmap_netcdf(VAR,\
                Data+VAR+'_AFR-44_'+RCM_Model[i]+rcp85fix+str(k+1)+filefix[v])[0,0]\
                for k in range(N_region)])
    return Ref,Future

#=================================================== 
# Read for CMPI5:
def reading_GCM(Ref_GCM,Future_GCM):
    print("reading GCMs data ...")
    for v in range(N_plot):
        for i in range(N_model_GCM):
            Ref_GCM[v,i]=np.array([ctang.read_lonlatmap_netcdf(VAR,\
                Data+VAR+'_Amon_'+GCM_Model[i]+histfix_GCM+str(k+1)+\
                filefix_GCM[v])[0,0]\
                for k in range(N_region)])
            Future_GCM[v,i]=np.array([ctang.read_lonlatmap_netcdf(VAR,\
                Data+VAR+'_Amon_'+GCM_Model[i]+rcp85fix_GCM+str(k+1)+\
                filefix_GCM[v])[0,0]\
                for k in range(N_region)])
    return Ref_GCM,Future_GCM

#=================================================== 
# get significance in earch region
def significanceMap(N_model,N_region,t_value):
    for m in range(N_model):
        for i in range(N_region):
            Ref_map=np.array([ctang.read_lonlatmap_netcdf(VAR,\
                Data+VAR+'_AFR-44_'+RCM_Model[m]+histfix+str(i+1)+'.timmean.nc')])
            Future_map=np.array([ctang.read_lonlatmap_netcdf(VAR,\
                Data+VAR+'_AFR-44_'+RCM_Model[m]+rcp85fix+str(i+1)+'.timmean.nc')])
            Change_map = Future_map - Ref_map
            if np.abs(stats.ttest_1samp(Change_map[0].flatten(),0)[0]) > \
                    ctang.get_T_value(len(Change_map[0].flatten())):
                t_value[m,i] = 1
    return t_value
#=================================================== 
# get significance in earch region
def significanceMap_GCM(N_model_GCM,N_region,t_value):
    for m in range(N_model_GCM):
        for i in range(N_region):
            Ref_map_GCM=np.array([ctang.read_lonlatmap_netcdf(VAR,\
                Data+VAR+'_Amon_'+GCM_Model[m]+histfix_GCM+str(i+1)+'.timmean.nc')])
            Future_map_GCM=np.array([ctang.read_lonlatmap_netcdf(VAR,\
                Data+VAR+'_Amon_'+GCM_Model[m]+rcp85fix_GCM+str(i+1)+'.timmean.nc')])
            Change_map_GCM = Future_map_GCM - Ref_map_GCM
            if np.abs(stats.ttest_1samp(Change_map_GCM[0].flatten(),0)[0]) > \
                    ctang.get_T_value(len(Change_map_GCM[0].flatten())):
                t_value[m,i] = 1
    return t_value
#=================================================== 

#--------------------------------------------------- t-test
t_value = significanceMap(N_model,N_region,t_value)
t_value_GCM = significanceMap_GCM(N_model_GCM,N_region,t_value_GCM)

print t_value_GCM.shape
print t_value.shape


Ref,Future = reading(Ref,Future)
Ref_GCM,Future_GCM = reading_GCM(Ref_GCM,Future_GCM)

# GCMs have missing values:
Ref_GCM[Ref_GCM > 999] = np.nan
Future_GCM[Future_GCM > 999] = np.nan

print Ref.shape
print Ref_GCM.shape

#ctang.Save2mat('Ref',Ref) # save to txt file to save time
#ctang.Save2mat('Future',Future) # save to txt file to save time

# Ref = ctang.Loadmat('Ref')
# Future = ctang.Loadmat('Future')

#=================================================== cal
#Ref[0][ Ref[0] == 0 ] = np.nan
Changes = np.array([(Future[i] - Ref[i])*100/Ref[0] for i in range(N_plot)])
Emean_Changes = np.mean( Changes, axis = 1)

Changes_GCM = np.array([(Future_GCM[i] - Ref_GCM[i])*100/Ref_GCM[0] \
        for i in range(N_plot)])
Emean_Changes_GCM = np.nanmean( Changes_GCM, axis = 1)
#=================================================== plot
ind = np.arange(1,N_region+1)   # num of regions
width = 0.35                    # the width of the bars


fig, axes = plt.subplots(nrows=2, ncols=2,figsize=(10, 6),facecolor='w', edgecolor='k')
fig.subplots_adjust(bottom=0.15,hspace=0.8,wspace=0.4)
ax0, ax1, ax2, ax3 = axes.flatten()
axes = axes.flatten() # reshape plot (2*2) to 4*1


print ind
# plot gcm:
rects0 = ax0.bar(ind-width, Emean_Changes_GCM[0], width, color='orange',align='edge',zorder=1)
rects1 = ax1.bar(ind-width, Emean_Changes_GCM[1], width, color='orange',align='edge',zorder=1)
rects2 = ax2.bar(ind-width, Emean_Changes_GCM[2], width, color='orange',align='edge',zorder=1)
# rects3 = ax3.bar(ind-width, Emean_Changes_GCM[3], width, color='blue',align='edge',zorder=1)


# plot rcm:
rects0 = ax0.bar(ind, Emean_Changes[0], width, color='green',align='edge',zorder=1)
rects1 = ax1.bar(ind, Emean_Changes[1], width, color='green',align='edge',zorder=1)
rects2 = ax2.bar(ind, Emean_Changes[2], width, color='green',align='edge',zorder=1)
# rects3 = ax3.bar(ind, Emean_Changes[3], width, color='blue',align='edge',zorder=1)
#-------------------- t-test in p<0.05,nof=20,2.086
def plot_individual_model(modelchange,t_value):
    """
    input array should be in 11 * 9 
    
    """
    for v in range(N_plot):
        for i in range(N_region):
            for j in range(N_model):
                if abs(t_value[j,i]) > 0:
                    axes[v].scatter((ind[i]+width/2),(modelchange[v,j,i]),s=10,\
                        facecolors='b',edgecolors='b',zorder=2)
                else:
                    axes[v].scatter((ind[i]+width/2),(modelchange[v,j,i]),s=10,\
                        facecolors='none',edgecolors='r',zorder=2)

#--------------------------------------------------- end of function
#-------------------- t-test in p<0.05,nof=20,2.086
def plot_individual_model_GCM(modelchange,t_value_GCM):
    """
    input array should be in 11 * 9 
    
    """
    for v in range(N_plot):
        for i in range(N_region):
            for j in range(N_model_GCM):
                if abs(t_value_GCM[j,i]) > 0:
                    axes[v].scatter((ind[i]-width/2),(modelchange[v,j,i]),s=10,\
                        facecolors='b',edgecolors='b',zorder=2)
                else:
                    axes[v].scatter((ind[i]-width/2),(modelchange[v,j,i]),s=10,\
                        facecolors='none',edgecolors='r',zorder=2)

#--------------------------------------------------- end of function
print t_value_GCM.shape
print t_value.shape
plot_individual_model(Changes,t_value)
plot_individual_model_GCM(Changes_GCM,t_value_GCM)


#=================================================== 
# add some text for labels, title and axes ticks

Title=['mean changes', 'monthly variability', 'annual variability']

#ax0.legend((rects1[0], ('Men')))

# set ticks and tick labels
for k in range(N_plot):
    axes[k].set_xlim((width, N_region+1))
    axes[k].set_xticks(ind )
    axes[k].set_xticklabels(('1', '2', '3', '4', '5', '6', '7' ))
    axes[k].set_xlabel('Region',fontsize=14)
    axes[k].spines['left'].set_linewidth(2)
    axes[k].spines['bottom'].set_linewidth(2)
    axes[k].tick_params(direction='out',length=6,width=2,labelsize=12)

    #axes[0].set_ylim((-10, 10))
    #axes[1].set_ylim((-15, 15))
    #axes[2].set_ylim((-15, 15))
    #axes[3].set_ylim((-35, 35))

    # plt.yscale('log')
    # axes[k].set_yscale('log')

    #axes[k].set_yticks([-1, 0, 1])
    axes[k].set_title(Title[k],fontsize=14)
    axes[k].set_ylabel('Changes   %', fontsize=14)
    #axes[k].set_yscale('symlog')
    
    axes[k].set_axisbelow(True)
    axes[k].yaxis.grid(color='gray', linestyle='dashed')
    axes[k].xaxis.grid(color='gray', linestyle='dashed')
    #axes[k].yaxis.set_minor_formatter(NullFormatter())
    minorLocator = AutoMinorLocator()
    axes[k].yaxis.set_minor_locator(minorLocator)

    # Hide the right and top spines
    axes[k].spines['right'].set_visible(False)
    axes[k].spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    axes[k].yaxis.set_ticks_position('left')
    axes[k].xaxis.set_ticks_position('bottom')

    axes[k].axhline(y=0, xmin=0, xmax=N_region+1,color='black',linewidth=2)


ctang.empty_plot(axes[3])
#=================================================== end of plot
plt.suptitle('Projected changes in the mean and the time variability of SSR',fontsize=16)

#plt.savefig('time.variability.eps', format='eps')
plt.savefig('time.variability.rsds.png')
plt.show()
