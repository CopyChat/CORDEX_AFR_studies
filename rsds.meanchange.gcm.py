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
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/seasonal.meanchanges.rcp45.rcp85.mid.end.rsds/'

# output:
Output='rsds.meanchange.rcp45.rcp85.min.end.gcm'

#=================================================== Definitions
N_model = 11
VAR = 'rsds' #,'tas','sfcWind') #,'PVpot')
#---------------------------------------------------  data

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

RCP=('hist','rcp45','rcp85')

YEAR=('1970-1999','2035-2065','2070-2099')

Season='JJA'
Season='DJF'

#================================================== reading data
# 21 maps

# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf_1D(\
        Data+VAR+'_Amon_'+GCM_Model[1]+'.'+RCP[0]+'.'+YEAR[0]+'.SA.'+str(Season)+'.timmean.remap.nc')

# rsds_Amon_NorESM1-M.rcp45.2070-2099.SA.JJA.timmean.remap.nc
# print Data+VAR+'_Amon_'+GCM_Model[1]+'.'+RCP[0]+'.'+YEAR[0]+'.SA.'+str(Season)+'.timmean.remap.nc'


#===================================================  reading GCMs
# Read Ensmean of timmean for CORDEX 
GCM_hist_JJA=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+'.hist.1970-1999.SA.JJA.timmean.remap.nc')\
    for i in range(N_model)])

GCM_rcp45_JJA_mid=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+'.rcp45.2036-2065.SA.JJA.timmean.remap.nc')\
    for i in range(N_model)])

GCM_rcp45_JJA_end=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+'.rcp45.2070-2099.SA.JJA.timmean.remap.nc')\
    for i in range(N_model)])

GCM_rcp85_JJA_mid=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+'.rcp85.2036-2065.SA.JJA.timmean.remap.nc')\
    for i in range(N_model)])

GCM_rcp85_JJA_end=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+'.rcp85.2070-2099.SA.JJA.timmean.remap.nc')\
    for i in range(N_model)])

#=================================================== missing value
GCM_hist_JJA[GCM_hist_JJA < 0]=np.nan
GCM_hist_JJA[GCM_hist_JJA < 0]=np.nan

GCM_rcp45_JJA_mid[GCM_rcp45_JJA_mid < 0]=np.nan
GCM_rcp45_JJA_end[GCM_rcp45_JJA_end < 0]=np.nan

GCM_rcp85_JJA_mid[GCM_rcp85_JJA_mid < 0]=np.nan
GCM_rcp85_JJA_end[GCM_rcp85_JJA_end < 0]=np.nan

#=================================================== changes GCMs

GCM_rcp45_JJA_mid_Change=GCM_rcp45_JJA_mid-GCM_hist_JJA
GCM_rcp45_JJA_end_Change=GCM_rcp45_JJA_end-GCM_hist_JJA

GCM_rcp85_JJA_mid_Change=GCM_rcp85_JJA_mid-GCM_hist_JJA 
GCM_rcp85_JJA_end_Change=GCM_rcp85_JJA_end-GCM_hist_JJA

#=================================================== Cal significant
## input: changes map
## output: significant map, NB

ZEROS=np.zeros((\
        GCM_rcp45_JJA_end_Change.shape[1],\
        GCM_rcp45_JJA_end_Change.shape[2]))


GCM_rcp45_JJA_mid_NB=ZEROS
GCM_rcp45_JJA_end_NB=ZEROS
GCM_rcp85_JJA_mid_NB=ZEROS
GCM_rcp85_JJA_end_NB=ZEROS

# GCM_rcp45_JJA_mid_Sig=ZEROS
# GCM_rcp45_JJA_end_Sig=ZEROS
# GCM_rcp85_JJA_mid_Sig=ZEROS
# GCM_rcp85_JJA_end_Sig=ZEROS
#--------------------------------------------------- 
def Significant_map(multi_model_map,NB):
    Sig=np.zeros((\
            multi_model_map.shape[1],\
            multi_model_map.shape[2]))

    print Sig.shape,Sig
    for lat in range(multi_model_map.shape[1]):
        for lon in range(multi_model_map.shape[2]):

            # lat=34
            # lon=48
            print Sig.shape,Sig,lat,lon
            grid = multi_model_map[:,lat,lon]
            print grid
            print GCM_rcp45_JJA_end_Change[:,lat,lon]
            grid_nonnan=[value for value in grid if not math.isnan(value)]
            print len(grid_nonnan)
            print grid_nonnan
            print(lat,lon,len(grid_nonnan))

            if len(grid_nonnan) < 1:
                Sig[lat,lon]=0
                NB[lat,lon]=0
                print "----------------- bad point"
            else:
                # print grid
                Significance = stats.ttest_1samp(grid_nonnan,0)[0]
                T=ctang.get_T_value(len(grid))
                

                if np.abs(Significance) > T:

                    print "Significance = "+str(Significance)
                    Sig[lat,lon]=1
                else:
                    Sig[lat,lon]=0

                print Sig.shape,Sig,lat,lon,Sig[lat,lon],Significance,T

                # how many models have value significant
                cc=0
                for jk in range(len(grid_nonnan)):
                    print grid_nonnan[jk],np.mean(grid_nonnan)
                    if (grid_nonnan[jk]*np.mean(grid_nonnan) > 0):
                        if (np.abs(grid_nonnan[jk])>T):
                            print lat,lon,grid_nonnan[jk],np.mean(grid_nonnan),T
                            cc+=1
                NB[lat,lon]=cc
            print lat,lon,grid_nonnan[jk],np.mean(grid_nonnan),T,cc,NB[lat,lon]
    return Sig
#--------------------------------------------------- 
# GCM_rcp45_JJA_mid_Sig=\
        # Significant_map(GCM_rcp45_JJA_mid_Change,\
        # GCM_rcp45_JJA_mid_NB )

GCM_rcp45_JJA_end_Sig=\
        Significant_map(GCM_rcp45_JJA_end_Change,\
        GCM_rcp45_JJA_end_NB )

# GCM_rcp85_JJA_mid_Sig=\
        # Significant_map(GCM_rcp85_JJA_mid_Change,\
        # GCM_rcp85_JJA_mid_NB )

# GCM_rcp85_JJA_end_Sig=\
        # Significant_map(GCM_rcp85_JJA_end_Change,\
        # GCM_rcp85_JJA_end_NB )
        
# GCM_rcp45_JJA_end_Sig,GCM_rcp45_JJA_end_NB = \
        # Significant_map(GCM_rcp45_JJA_end_Change)
# GCM_rcp85_JJA_mid_Sig,GCM_rcp85_JJA_mid_NB = \
        # Significant_map(GCM_rcp85_JJA_mid_Change)
# GCM_rcp85_JJA_end_Sig,GCM_rcp85_JJA_end_NB = \
        # Significant_map(GCM_rcp85_JJA_end_Change)

#=================================================== check
print GCM_rcp45_JJA_end_NB
print GCM_rcp45_JJA_end_Sig
print GCM_rcp45_JJA_end_NB[34,44]
print GCM_rcp45_JJA_end_Sig[34,44]
print len(GCM_rcp45_JJA_end_Change[:,34,44])

print ctang.get_T_value(len(GCM_rcp45_JJA_end_Change[:,34,44]))
print GCM_rcp45_JJA_end_Change[:,34,44]

#=================================================== 

quit()

# plotting array:
Ensmean_change=np.nanmean(mean_change,axis=0)
Ensmean_change_ttest=np.nanmean(t_value,axis=0)

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

plt.suptitle('11 CMPI5 mean changes of SSR (%)',fontsize=14)

print output+'.eps'
plt.savefig(output+'.eps',format='eps')
plt.show()

quit()
