#!/usr/bin/env python
"""
========
Ctang, A map of mean max and min of ensembles
        from CORDEX AFR-44, in Southern Africa
        Data was restored on titan
========
"""
import pdb
import math
import subprocess
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.basemap import Basemap , addcyclic
import textwrap

# to load my functions
import sys 
sys.path.append('/Users/ctang/Code/Python/')
import ctang


#=================================================== I/O info:
# input data:
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/GCM-RCM_map/'


# output:
output='rsds.meanchange.rcm'
DataMat=output

#=================================================== Definitions
N_region = 9
N_model = 21
VAR = 'rsds' #,'tas','sfcWind') #,'PVpot')
#---------------------------------------------------  data
# 21 * 4 table: 21 models vs 4 vars

    #               region_1 region_2 ... region_N_region
    #  model_1
    #  model_2
    #  ...
    #  model_N_model
    #  ens
#--------------------------------------------------- 
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

#================================================== reading data
# 21 maps
filefix=(\
        '.hist.day.1970-1999.SA.timmean.nc' ,\
        '.rcp85.day.2070-2099.SA.timmean.nc')

# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf(\
        Data+VAR+'_AFR-44_'+RCM_Model[1]+filefix[0])

# Read Ensmean of timmean for CORDEX 
## reading 21 models : (21,90,137)
mean_ref=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[0])\
    for i in range(N_model)])
mean_future=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[1]) \
    for i in range(N_model)])
#=================================================== Cal

mean_change=(mean_future-mean_ref)*100/mean_ref

t_value=mean_change

#--------------------------------------------------- 
# for t test in each point
#--------------------------------------------------- 
#t_value in N_model,lat,lon
def significant_map(t_value):
    count=0
    for lat in range(t_value.shape[1]):
        for lon in range(t_value.shape[2]):
            grid = t_value[:,lat,lon]
            grid = grid[grid < 999]  # remove missing values
            # print grid

            print(lat,lon,len(grid))
            if len(grid) < 1:
                t_value[:,lat,lon]=np.NaN
                print lat,lon
            else:
                Sig = stats.ttest_1samp(grid,0)[0]
                if np.abs(Sig) > ctang.get_T_value(len(grid)):
                    count+=1
                else:
                    t_value[:,lat,lon]=np.NaN
    return t_value
#--------------------------------------------------- 
#--------------------------------------------------- 

def robust(t_value):
#   The following criteria are used to
#   consider them as either uncertain or negligible:
#       negligible: NaN >= N_model/2
#       uncertain:  at least 2 significant individual differ in sign
#       robust:     NaN < N_model/2, sum(abs(x_i)) = abs(sum(x_i))

#   input: 21 maps with t values
#   return: 21 maps with t values
    count=1
    for lat in range(t_value.shape[1]):
        for lon in range(t_value.shape[2]):
            if np.isnan(t_value[:,lat,lon]).sum() >= N_model/2:     # negligible
                #t_value[v,:,lat,lon]=[ 99999 for t in range(N_model)] 
                t_value[:,lat,lon]=np.nan
            else:
                if (N_model-np.isnan(t_value[:,lat,lon]).sum()) > np.abs(np.nansum(np.sign(t_value[:,lat,lon]))):
                    #t_value[v,:,lat,lon]=[-99999 for t in range(N_model)] # uncertain
                    t_value[:,lat,lon]=np.nan
                else:
                    t_value[:,lat,lon]=t_value[:,lat,lon]
                    #print "get robust point",count
                    count += 1
    return t_value
#========================================= end of function

# t_value=robust(significant_map(t_value))
# t_value=significant_map(t_value)

# save to txt file to save time
# ctang.Save2mat(output,t_value) 

# reading from mat file
t_value=ctang.Loadmat(output+'.mat')
print np.isnan(t_value).sum()

# plotting array:
Ensmean_change_ttest=np.nanmean(t_value,axis=0)
Ensmean_change=np.nanmean(mean_change,axis=0)

#=================================================== plot
degree_sign= u'\N{DEGREE SIGN}'
Title='rsds change %'
Unit=(('(W/m2)','(%)'))
TTT=(('Mean RSDS changes RCP8.5 ', 'PVpot changes RCP8.5 '))
#--------------------------------------------------- 
def PlotMap(array2D,lons,lats,axx,vmin,vmax):
    cmap = plt.cm.jet
    cmap = plt.cm.seismic
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,9)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    map=Basemap(projection='cyl',llcrnrlat=lats[:,0].min(),urcrnrlat=lats[:,0].max(),\
        llcrnrlon=lons[0,:].min(),urcrnrlon=lons[0,:].max(),resolution='l')
    ctang.setMap(map)

    x,y=map(lats,lons)

    map.pcolormesh(y,x,array2D,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    #axx.axis('off')
    axx.xaxis.set_visible(False)
    axx.yaxis.set_visible(False)

    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=0.8) 
    cb.ax.tick_params(labelsize=10) 
#--------------------------------------------------- 

#=================================================== main
fig, axes = plt.subplots(nrows=1, ncols=2,\
        figsize=(12, 6),facecolor='w', edgecolor='k') # figsize=(w,h)
fig.subplots_adjust(hspace=0.35,top=0.99,wspace=0.3,bottom=0.25)

plt.sca(axes[0]) # active shis subplot for GCM
ax=axes[0]
PlotMap(Ensmean_change,lons,lats,ax,-10,10)
plt.sca(axes[1]) # active shis subplot for GCM
ax=axes[1]
PlotMap(Ensmean_change_ttest,lons,lats,ax,-10,10)

plt.suptitle('21 RCM experiments mean changes of SSR (%)',fontsize=14)

print output+'.eps'
plt.savefig(output+'.eps',format='eps')
plt.show()

quit()
