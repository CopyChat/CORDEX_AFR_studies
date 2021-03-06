#!/usr/bin/env python
"""
========
Ctang, A map of mean max and min of ensembles
        from CORDEX AFR-44, in Southern Africa
        Data was restored on titan
========
"""
import math
import subprocess
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.basemap import Basemap , addcyclic
import textwrap

# to load my functions
import sys 
sys.path.append('/Users/ctang/Code/Python/')
import ctang

#=================================================== Definitions
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/'
N_region = 9
N_model = 10
T=2.086 # nof=20 P<0.05
VAR = ('rsds','tas','sfcWind') #,'PVpot')
N_var=len(VAR)+1

#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

    #               region_1 region_2 ... region_N_region
    #  model_1
    #  model_2
    #  ...
    #  model_N_model
    #  ens
#--------------------------------------------------- 

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

#each map is mean_ref[0,0,:,:].shape

#=================================================== test
print np.nanmax(np.abs([np.nan, 1, 2,-3]))
#=================================================== test

# Read lon,lat
lons,lats=ctang.read_lonlat_netcdf(\
        Data+VAR[0]+'_AFR-44_'+MODEL[0]+'.hist.day.1970-1999.SA.timmean.nc')

print lons.shape

mean_ref=np.zeros((len(VAR),N_model,lons.shape[0],lons.shape[1]))
mean_future=np.zeros((len(VAR),N_model,lons.shape[0],lons.shape[1]))
mean_change=np.zeros((len(VAR)+1,N_model,lons.shape[0],lons.shape[1]))

for v in range(len(VAR)):
    mean_ref[v]=np.array([[ctang.read_lonlatmap_netcdf(VAR[v],\
            Data+VAR[v]+'_AFR-44_'+MODEL[i]+'.hist.day.1970-1999.SA.timmean.nc')] \
            for i in range(N_model)])[:,0,:,:]

    mean_future[v]=np.array([[ctang.read_lonlatmap_netcdf(VAR[v],\
            Data+VAR[v]+'_AFR-44_'+MODEL[i]+'.rcp85.day.2070-2099.SA.timmean.nc')] \
            for i in range(N_model)])[:,0,:,:]

print mean_ref.shape

# Calculate change  4D: len(VAR):N_model: lat : lon
mean_change[0]=np.array((mean_future[0,:,:,:]-mean_ref[0,:,:,:]))                 # RADS
mean_change[1]=np.array((ctang.PVpot(mean_future[0,:,:,:],mean_future[1,:,:,:])-\
        ctang.PVpot(mean_ref[0,:,:,:],mean_ref[1,:,:,:]))*100/\
        ctang.PVpot(mean_ref[0,:,:,:],mean_ref[1,:,:,:]))
mean_change[2]=np.array((mean_future[1,:,:,:]-mean_ref[1,:,:,:]))  # TAS
mean_change[3]=np.array((mean_future[2,:,:,:]-mean_ref[2,:,:,:])*100/mean_ref[2,:,:,:])  # VWS
print mean_change.shape
#=================================================== end of read



#--------------------------------------------------- 
# for t test in each point
t_value=mean_change
#--------------------------------------------------- 
def significant_map(mean_change,t_value):
    for v in range(N_var):
        for m in range(N_model):
            for lat in range(mean_change.shape[2]):
                for lon in range(mean_change.shape[3]):
                    if np.abs(stats.ttest_1samp(\
                            mean_change[v,:,lat,lon],mean_change[v,m,lat,lon])[0]) < T:
                        t_value[v,m,lat,lon]=np.NaN
    return t_value
#--------------------------------------------------- 
t_value=significant_map(mean_change,t_value)
#t_value=np.reshape(t_value,(N_model*(len(VAR)+1),lons.shape[0],lons.shape[1]))
print t_value.shape

Max=[np.nanmax(np.abs(t_value[t,:,:,:])) for t in range(N_var)]

#=================================================== plot
Title='Changes in RSDS(W/m2), PVpot(%), TAS($^\circ$C) and VWS(%) - RCP8.5'

# Define a plot of multi subplots
#fig, axes = plt.subplots(nrows=N_model, ncols=len(VAR),\
        #figsize=(len(VAR),N_model),facecolor='w', edgecolor='k')
#fig.subplots_adjust(bottom=0.15,hspace=0.8,wspace=0.4)
##axes = axes.flatten() # reshape plots to 1D if needed
#print axes.shape

#=================================================== 
# define the colormap
cmap = plt.cm.jet
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
#cmaplist[0] = (.5,.5,.5,1.0)
#cmaplist[-1] = (.8,1.0,1.0,1.0)
# create the new map
#cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
# define the bins and normalize
bounds = np.linspace(-21,21,11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

#=================================================== 
fig, axes = plt.subplots(nrows=N_model, ncols=N_var,\
        figsize=(11, 20),facecolor='w', edgecolor='k')
fig.subplots_adjust(bottom=0.15,hspace=0.5,wspace=0.4)

#fig=plt.figure(figsize=(24,16))
gs = gridspec.GridSpec(N_model, N_var)
gs.update(wspace=0.1, hspace=0.1, left=0.1, right=0.4, bottom=0.1, top=0.9) 

for v in range(N_var):
    for m in range(N_model):
        plt.sca(axes[m,v]) # active shis subplot
        map=Basemap(projection='cyl',llcrnrlat=lats.min(),urcrnrlat=lats.max(),\
            llcrnrlon=lons.min(),urcrnrlon=lons.max(),resolution='h')
        map.drawcoastlines(linewidth=1)
        map.drawparallels(np.arange(-90.,91.,10.),labels=[1,0,0,0],linewidth=0.5)
        map.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],linewidth=0.5)
        map.drawmapboundary()
        map.drawcountries()
        x,y=map(lats,lons)
        map.pcolormesh(y,x,t_value[v,m,:,:],cmap=cmap,\
                vmin=-Max[v],vmax=Max[v])
        #axes[m,v].set_xlabel('lon',fontsize=0.4)
                
        #axes[m,v].tick_params(
                #axis='both',          # changes apply to the x-axis
                #which='both',      # both major and minor ticks are affected
                #bottom='off',      # ticks along the bottom edge are off
                #top='off',         # ticks along the top edge are off
                #labelbottom='off') # labels along the bottom edge are off
        axes[m,v].set_xticks([])
        axes[m,v].set_yticks([])
        axes[m,v].xaxis.set_visible(False)
        plt.annotate(str(m+1), xy=(1, 0), fontsize=9,\
            horizontalalignment='left', verticalalignment='top')
        plt.colorbar(orientation='horizontal',shrink=0.9) # draw colorbar


#plt.title("\n".join(textwrap.wrap(Title,55)))

plt.tight_layout()
plt.savefig('mpa.add_value.eps',format='eps')
plt.savefig('map.add_value.pdf')
plt.savefig('map.add_value.png')
plt.show()

quit()
