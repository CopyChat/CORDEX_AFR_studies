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
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/'

# output:
output='rsds.meanchange.gcm'
DataMat=output

#=================================================== Definitions
N_region = 9
N_model = 11
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

GCM_Model=(\
    'CNRM-CM5',\
    'CSIRO-Mk3-6-0',\
    'CanESM2',\
    'EC-EARTH',\
    'GFDL-ESM2M',\
    'HadGEM2-ES',\
    'IPSL-CM5A-LR',\
    'IPSL-CM5A-MR',\
    'MIROC5',\
    'MPI-ESM-LR',\
    'NorESM1-M',\
    )

#================================================== reading data
# 21 maps
filefix=(\
        '.hist.1970-1999.SA.timmean.remap.nc' ,\
        '.rcp85.2070-2099.SA.timmean.remap.nc')

# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf_1D(\
        Data+VAR+'_Amon_'+GCM_Model[1]+filefix[0])

# Read Ensmean of timmean for CORDEX 
## reading 21 models : (21,90,137)
mean_ref=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+filefix[0])\
    for i in range(N_model)])
mean_future=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+filefix[1]) \
    for i in range(N_model)])
#=================================================== Cal

mean_change=(mean_future-mean_ref)*100/mean_ref

t_value=mean_change

#--------------------------------------------------- 
# for t test in each point
T=2.228 # dof=10 two tail 0.05
#--------------------------------------------------- 
#t_value in N_model,lat,lon
def significant_map(t_value):
    count=0
    for m in range(N_model):
        print 'calculating t test in dim: ',str(m)
        for lat in range(t_value.shape[1]):
            for lon in range(t_value.shape[2]):
                if abs(stats.ttest_1samp( t_value[:,lat,lon],0)[0]) > T:
                    count+=1
                else:
                    t_value[m,lat,lon]=np.NaN
    print("get good point in %: ", \
        str(count*100/N_model/t_value.shape[1]/t_value[2]))
    return t_value
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

plt.suptitle('11 CMPI5 models mean changes of SSR (%)',fontsize=14)

print output+'.eps'
plt.savefig(output+'.eps',format='eps')
plt.show()

quit()
