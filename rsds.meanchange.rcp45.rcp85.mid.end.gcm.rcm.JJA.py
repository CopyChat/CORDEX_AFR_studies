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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import textwrap

# to load my functions
import sys 
sys.path.append('/Users/ctang/Code/Python/')
import ctang


#=================================================== I/O info:
# input data:
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/seasonal.meanchanges.rcp45.rcp85.mid.end.rsds/'

# output:
Output='rsds.meanchange.rcp45.rcp85.mid.end'

#=================================================== Definitions
VAR = 'rsds' #,'tas','sfcWind') #,'PVpot')
#---------------------------------------------------  data
RCM_Model=(\
        'CCCma-CanESM2_SMHI-RCA4_v1',\
        'CNRM-CERFACS-CNRM-CM5_CLMcom-CCLM4-8-17_v1',\
        'CNRM-CERFACS-CNRM-CM5_SMHI-RCA4_v1',\
        'CSIRO-QCCCE-CSIRO-Mk3-6-0_SMHI-RCA4_v1',\
        'ICHEC-EC-EARTH_CLMcom-CCLM4-8-17_v1',\
        'ICHEC-EC-EARTH_DMI-HIRHAM5_v2',\
        'ICHEC-EC-EARTH_KNMI-RACMO22T_v1',\
        'ICHEC-EC-EARTH_MPI-CSC-REMO2009_v1',\
        'ICHEC-EC-EARTH_SMHI-RCA4_v1',\
        'IPSL-IPSL-CM5A-MR_SMHI-RCA4_v1',\
        'MIROC-MIROC5_SMHI-RCA4_v1',\
        'MOHC-HadGEM2-ES_CLMcom-CCLM4-8-17_v1',\
        'MOHC-HadGEM2-ES_KNMI-RACMO22T_v2',\
        'MOHC-HadGEM2-ES_SMHI-RCA4_v1',\
        'MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17_v1',\
        'MPI-M-MPI-ESM-LR_MPI-CSC-REMO2009_v1',\
        'MPI-M-MPI-ESM-LR_SMHI-RCA4_v1',\
        'NCC-NorESM1-M_DMI-HIRHAM5_v1',\
        'NCC-NorESM1-M_SMHI-RCA4_v1',\
        'NOAA-GFDL-GFDL-ESM2M_SMHI-RCA4_v1',\
        )


GCM_Model=(\
    'CNRM-CM5',\
    'CSIRO-Mk3-6-0',\
    'CanESM2',\
    'EC-EARTH',\
    # 'GFDL-ESM2M',\
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
lons_GCM,lats_GCM=ctang.read_lonlat_netcdf_1D(\
        Data+VAR+'_Amon_'+GCM_Model[1]+'.'+RCP[0]+'.'+YEAR[0]+'.SA.'+str(Season)+'.timmean.remap.nc')


lons_RCM,lats_RCM=ctang.read_lonlat_netcdf(\
        Data+VAR+'_AFR-44_'+RCM_Model[1]+'.rcp85.day.2036-2065.SA.'+str(Season)+'.timmean.nc')
# rsds_AFR-44_MIROC-MIROC5_SMHI-RCA4_v1.rcp85.day.2036-2065.SA.JJA.timmean.nc

#===================================================  reading RCMs
# Read Ensmean of timmean for CORDEX 
RCM_hist_JJA=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_AFR-44_'+RCM_Model[i]+\
    '.hist.day.1970-1999.SA.'+str(Season)+'.timmean.nc')\
    for i in range(len(RCM_Model))])

RCM_rcp45_JJA_mid=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_AFR-44_'+RCM_Model[i]+\
    '.rcp45.day.2036-2065.SA.'+str(Season)+'.timmean.nc')\
    for i in range(len(RCM_Model))])

RCM_rcp45_JJA_end=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_AFR-44_'+RCM_Model[i]+\
    '.rcp45.day.2070-2099.SA.'+str(Season)+'.timmean.nc')\
    for i in range(len(RCM_Model))])

RCM_rcp85_JJA_mid=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_AFR-44_'+RCM_Model[i]+\
    '.rcp85.day.2036-2065.SA.'+str(Season)+'.timmean.nc')\
    for i in range(len(RCM_Model))])

RCM_rcp85_JJA_end=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_AFR-44_'+RCM_Model[i]+\
    '.rcp85.day.2070-2099.SA.'+str(Season)+'.timmean.nc')\
    for i in range(len(RCM_Model))])
#===================================================  reading GCMs
# Read Ensmean of timmean for CORDEX 
GCM_hist_JJA=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+'.hist.1970-1999.SA.JJA.timmean.remap.nc')\
    for i in range(len(GCM_Model))])

GCM_rcp45_JJA_mid=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+'.rcp45.2036-2065.SA.JJA.timmean.remap.nc')\
    for i in range(len(GCM_Model))])

GCM_rcp45_JJA_end=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+'.rcp45.2070-2099.SA.JJA.timmean.remap.nc')\
    for i in range(len(GCM_Model))])

GCM_rcp85_JJA_mid=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+'.rcp85.2036-2065.SA.JJA.timmean.remap.nc')\
    for i in range(len(GCM_Model))])

GCM_rcp85_JJA_end=np.array([ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_Amon_'+GCM_Model[i]+'.rcp85.2070-2099.SA.JJA.timmean.remap.nc')\
    for i in range(len(GCM_Model))])

#=================================================== missing value
RCM_hist_JJA[RCM_hist_JJA < 0]=np.nan

RCM_rcp45_JJA_mid[RCM_rcp45_JJA_mid < 0]=np.nan
RCM_rcp45_JJA_end[RCM_rcp45_JJA_end < 0]=np.nan

RCM_rcp85_JJA_mid[RCM_rcp85_JJA_mid < 0]=np.nan
RCM_rcp85_JJA_end[RCM_rcp85_JJA_end < 0]=np.nan

GCM_hist_JJA[GCM_hist_JJA < 0]=np.nan

GCM_rcp45_JJA_mid[GCM_rcp45_JJA_mid < 0]=np.nan
GCM_rcp45_JJA_end[GCM_rcp45_JJA_end < 0]=np.nan

GCM_rcp85_JJA_mid[GCM_rcp85_JJA_mid < 0]=np.nan
GCM_rcp85_JJA_end[GCM_rcp85_JJA_end < 0]=np.nan

#=================================================== changes GCMs
RCM_rcp45_JJA_mid_Change=RCM_rcp45_JJA_mid-RCM_hist_JJA
RCM_rcp45_JJA_end_Change=RCM_rcp45_JJA_end-RCM_hist_JJA

RCM_rcp85_JJA_mid_Change=RCM_rcp85_JJA_mid-RCM_hist_JJA 
RCM_rcp85_JJA_end_Change=RCM_rcp85_JJA_end-RCM_hist_JJA

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
ZEROS2=np.zeros((\
        RCM_rcp45_JJA_end_Change.shape[1],\
        RCM_rcp45_JJA_end_Change.shape[2]))

RCM_rcp45_JJA_mid_NB=ZEROS2
RCM_rcp45_JJA_end_NB=ZEROS2
RCM_rcp85_JJA_mid_NB=ZEROS2
RCM_rcp85_JJA_end_NB=ZEROS2

GCM_rcp45_JJA_mid_NB=ZEROS
GCM_rcp45_JJA_end_NB=ZEROS
GCM_rcp85_JJA_mid_NB=ZEROS
GCM_rcp85_JJA_end_NB=ZEROS

#=================================================== 
def SigNB_map(multi_model_map):
    NB=np.zeros((\
            multi_model_map.shape[1],\
            multi_model_map.shape[2]))
    for lat in range(multi_model_map.shape[1]):
        for lon in range(multi_model_map.shape[2]):
            grid = multi_model_map[:,lat,lon]
            grid_nonnan=[value for value in grid if not math.isnan(value)]

            if len(grid_nonnan) < 1:
                Sig[lat,lon]=0
                NB[lat,lon]=0
                print "----------------- bad point"
            else:
                # print grid
                T=ctang.get_T_value(len(grid))
                cc=0
                for jk in range(len(grid_nonnan)):
                    if (grid_nonnan[jk]*np.mean(grid_nonnan) > 0):
                        if (np.abs(grid_nonnan[jk])>T):
                            cc+=1
                NB[lat,lon]=cc
            # print lat,lon,grid_nonnan[jk],np.mean(grid_nonnan),T,cc,NB[lat,lon]
    return NB
#--------------------------------------------------- 
#--------------------------------------------------- 
def Significant_map(multi_model_map):
    Sig=np.zeros((\
            multi_model_map.shape[1],\
            multi_model_map.shape[2]))

    for lat in range(multi_model_map.shape[1]):
        for lon in range(multi_model_map.shape[2]):

            grid = multi_model_map[:,lat,lon]
            grid_nonnan=[value for value in grid if not math.isnan(value)]

            if len(grid_nonnan) < 1:
                Sig[lat,lon]=0
                NB[lat,lon]=0
                print "----------------- bad point"
            else:
                Significance = stats.ttest_1samp(grid_nonnan,0)[0]
                T=ctang.get_T_value(len(grid))
                

                if np.abs(Significance) > T:

                    # print "Significance = "+str(Significance)
                    Sig[lat,lon]=1
                else:
                    Sig[lat,lon]=0

                # print Sig.shape,Sig,lat,lon,Sig[lat,lon],Significance,T

    return Sig
#--------------------------------------------------- 
RCM_rcp45_JJA_mid_Sig=\
        Significant_map(RCM_rcp45_JJA_mid_Change)

RCM_rcp45_JJA_end_Sig=\
        Significant_map(RCM_rcp45_JJA_end_Change)

RCM_rcp85_JJA_mid_Sig=\
        Significant_map(RCM_rcp85_JJA_mid_Change)

RCM_rcp85_JJA_end_Sig=\
        Significant_map(RCM_rcp85_JJA_end_Change)
RCM_rcp45_JJA_mid_Sig=\
        Significant_map(RCM_rcp45_JJA_mid_Change)
        
RCM_rcp45_JJA_mid_NB=SigNB_map(RCM_rcp45_JJA_mid_Change)
RCM_rcp45_JJA_end_NB=SigNB_map(RCM_rcp45_JJA_end_Change)

RCM_rcp85_JJA_mid_NB=SigNB_map(RCM_rcp85_JJA_mid_Change)
RCM_rcp85_JJA_end_NB=SigNB_map(RCM_rcp85_JJA_end_Change)
#=================================================== 

GCM_rcp45_JJA_mid_Sig=\
        Significant_map(GCM_rcp45_JJA_mid_Change)

GCM_rcp45_JJA_end_Sig=\
        Significant_map(GCM_rcp45_JJA_end_Change)

GCM_rcp85_JJA_mid_Sig=\
        Significant_map(GCM_rcp85_JJA_mid_Change)

GCM_rcp85_JJA_end_Sig=\
        Significant_map(GCM_rcp85_JJA_end_Change)
GCM_rcp45_JJA_mid_Sig=\
        Significant_map(GCM_rcp45_JJA_mid_Change)
        
GCM_rcp45_JJA_mid_NB=SigNB_map(GCM_rcp45_JJA_mid_Change)
GCM_rcp45_JJA_end_NB=SigNB_map(GCM_rcp45_JJA_end_Change)

GCM_rcp85_JJA_mid_NB=SigNB_map(GCM_rcp85_JJA_mid_Change)
GCM_rcp85_JJA_end_NB=SigNB_map(GCM_rcp85_JJA_end_Change)

#====================================== check
# print GCM_rcp45_JJA_end_NB[34,44]
# print GCM_rcp45_JJA_end_Sig[34,44]
# print len(GCM_rcp45_JJA_end_Change[:,34,44])

# print ctang.get_T_value(len(GCM_rcp45_JJA_end_Change[:,34,44]))
# print GCM_rcp45_JJA_end_Change[:,34,44]

#=================================================== 

#===================================== multi model mean:
RCM_rcp45_JJA_mid_Mean=np.nanmean(RCM_rcp45_JJA_mid_Change,\
        axis=0)
RCM_rcp45_JJA_end_Mean=np.nanmean(RCM_rcp45_JJA_end_Change,\
        axis=0)

RCM_rcp85_JJA_mid_Mean=np.nanmean(RCM_rcp85_JJA_mid_Change,\
        axis=0)
RCM_rcp85_JJA_end_Mean=np.nanmean(RCM_rcp85_JJA_end_Change,\
        axis=0)

GCM_rcp45_JJA_mid_Mean=np.nanmean(GCM_rcp45_JJA_mid_Change,\
        axis=0)
GCM_rcp45_JJA_end_Mean=np.nanmean(GCM_rcp45_JJA_end_Change,\
        axis=0)

GCM_rcp85_JJA_mid_Mean=np.nanmean(GCM_rcp85_JJA_mid_Change,\
        axis=0)
GCM_rcp85_JJA_end_Mean=np.nanmean(GCM_rcp85_JJA_end_Change,\
        axis=0)
print GCM_rcp45_JJA_mid_Mean.shape


#=========================================== mask the non-sig points

RCM_rcp45_JJA_mid_Mean[RCM_rcp45_JJA_mid_Sig == 0] = np.nan
RCM_rcp45_JJA_end_Mean[RCM_rcp45_JJA_end_Sig == 0] = np.nan
RCM_rcp85_JJA_mid_Mean[RCM_rcp85_JJA_mid_Sig == 0] = np.nan
RCM_rcp85_JJA_end_Mean[RCM_rcp85_JJA_end_Sig == 0] = np.nan

GCM_rcp45_JJA_mid_Mean[GCM_rcp45_JJA_mid_Sig == 0] = np.nan

GCM_rcp45_JJA_end_Mean[GCM_rcp45_JJA_end_Sig == 0] = np.nan
GCM_rcp85_JJA_mid_Mean[GCM_rcp85_JJA_mid_Sig == 0] = np.nan
GCM_rcp85_JJA_end_Mean[GCM_rcp85_JJA_end_Sig == 0] = np.nan


GCM_rcp45_JJA_mid_NB=GCM_rcp45_JJA_mid_NB*10
GCM_rcp45_JJA_end_NB=GCM_rcp45_JJA_end_NB*10

GCM_rcp85_JJA_mid_NB=GCM_rcp85_JJA_mid_NB*10
GCM_rcp85_JJA_end_NB=GCM_rcp85_JJA_end_NB*10

RCM_rcp45_JJA_mid_NB=RCM_rcp45_JJA_mid_NB*5
RCM_rcp45_JJA_end_NB=RCM_rcp45_JJA_end_NB*5

RCM_rcp85_JJA_mid_NB=RCM_rcp85_JJA_mid_NB*5
RCM_rcp85_JJA_end_NB=RCM_rcp85_JJA_end_NB*5

print GCM_rcp45_JJA_mid_NB
print RCM_rcp45_JJA_mid_NB
print GCM_rcp45_JJA_mid_NB.shape
print np.sum(GCM_rcp45_JJA_mid_NB==GCM_rcp85_JJA_mid_NB)

#=================================================== plot
degree_sign= u'\N{DEGREE SIGN}'
Title='rsds change %'
Unit=( '($\mathregular{W/m^{2}}$)','($\mathregular{W/m^{2}}$)')
TTT=(('Mean RSDS changes RCP4.5 ',\
    'Mean RSDS changes RCP8.5 '))


#=================================================== functions
def PlotSigMap(array2D,lons,lats,axx,vmin,vmax,title,number,cbar):
    cmap = plt.cm.jet
    cmap = plt.cm.seismic
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,13)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    map=Basemap(projection='cyl',llcrnrlat=-39,urcrnrlat=-1,\
        llcrnrlon=1,urcrnrlon=59,resolution='l')
    ctang.setMap_nostick(map)

    x,y=map(lats,lons)

    map.pcolormesh(y,x,array2D,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)

    plt.title(title,fontsize=9)

    # colorbar:
    divider = make_axes_locatable(axx)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    if (cbar==0):
        print ""
        ctang.empty_plot(cax)
    else:
        cb=plt.colorbar(cmap=plt.cm.jet,cax=cax) 
        # cb=plt.colorbar(cmap=plt.cm.jet,cax=cbaxes,orientation='vertical',shrink=0.6) 
        cb.ax.tick_params(labelsize=10) 
    ax.text(0.9, 0.9,str(Season), ha='center', va='center', transform=axx.transAxes)
    ax.text(0.1, 0.9,str(number), ha='center', va='center', transform=axx.transAxes)
#--------------------------------------------------- 
def PlotMap(array2D,lons,lats,axx,vmin,vmax,number,cbar):
    cmap = plt.cm.rainbow
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,11)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    map=Basemap(projection='cyl',llcrnrlat=lats[:,0].min(),urcrnrlat=lats[:,0].max(),\
        llcrnrlon=lons[0,:].min(),urcrnrlon=lons[0,:].max(),resolution='l')
    ctang.setMap_nostick(map)

    x,y=map(lats,lons)

    map.pcolormesh(y,x,array2D,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)


    # colorbar:
    divider = make_axes_locatable(axx)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    if (cbar==0):
        ctang.empty_plot(cax)
    else:
        cb=plt.colorbar(cmap=plt.cm.jet,cax=cax) 
        # cb=plt.colorbar(cmap=plt.cm.jet,cax=cbaxes,orientation='vertical',shrink=0.6) 
        cb.ax.tick_params(labelsize=10) 
    ax.text(0.9, 0.9,str(Season), ha='center', va='center', transform=axx.transAxes)
    # ax.text(0.1, 0.9,str(number), ha='center', va='center', transform=axx.transAxes)


title1='GCM trend to 2050 RCP4.5 '+str(Unit[0])
title2='GCM trend to 2050 RCP8.5 '+str(Unit[0])

title3='RCM trend to 2050 RCP4.5 '+str(Unit[0])
title4='RCM trend to 2050 RCP8.5 '+str(Unit[0])

title11='GCM trend to 2100 RCP4.5 '+str(Unit[0])
title22='GCM trend to 2100 RCP8.5 '+str(Unit[0])

title33='RCM trend to 2100 RCP4.5 '+str(Unit[0])
title44='RCM trend to 2100 RCP8.5 '+str(Unit[0])


title=''
#=================================================== main
fig, axes = plt.subplots(nrows=4, ncols=4,\
        figsize=(12, 8),facecolor='w', edgecolor='k') # figsize=(w,h)
fig.subplots_adjust(right=0.95,hspace=0.25,top=0.9,wspace=0.05,bottom=0.05)

plt.sca(axes[0,0]) # active shis subplot for GCM
ax=axes[0,0]
PlotSigMap(GCM_rcp45_JJA_mid_Mean,lons_GCM,lats_GCM,ax,-12,12,title1,'(a)',cbar=0)
plt.sca(axes[0,1]) # active shis subplot for GCM
ax=axes[0,1]
PlotSigMap(GCM_rcp85_JJA_mid_Mean,lons_GCM,lats_GCM,ax,-12,12,title2,'(b)',cbar=0)

plt.sca(axes[0,2]) # active shis subplot for GCM
ax=axes[0,2]
PlotSigMap(RCM_rcp45_JJA_mid_Mean,lons_RCM,lats_RCM,ax,-12,12,title3,'(c)',cbar=0)
plt.sca(axes[0,3]) # active shis subplot for GCM
ax=axes[0,3]
PlotSigMap(RCM_rcp85_JJA_mid_Mean,lons_RCM,lats_RCM,ax,-12,12,title4,'(d)',cbar=1)
#--------------------------------------------------- 

plt.sca(axes[1,0]) # active shis subplot for GCM
ax=axes[1,0]
PlotMap(GCM_rcp45_JJA_mid_NB,lons_GCM,lats_GCM,ax,-100,100,'(a)',cbar=0)
plt.sca(axes[1,1]) # active shis subplot for GCM
ax=axes[1,1]
PlotMap(GCM_rcp85_JJA_mid_NB,lons_GCM,lats_GCM,ax,-100,100,'(b)',cbar=0)

plt.sca(axes[1,2]) # active shis subplot for GCM
ax=axes[1,2]
PlotMap(RCM_rcp45_JJA_mid_NB,lons_RCM,lats_RCM,ax,0,100,'(c)',cbar=0)
plt.sca(axes[1,3]) # active shis subplot for GCM
ax=axes[1,3]
PlotMap(RCM_rcp85_JJA_mid_NB,lons_RCM,lats_RCM,ax,0,100,'(d)',cbar=1)
#--------------------------------------------------- 

plt.sca(axes[2,0]) # active shis subplot for GCM
ax=axes[2,0]
PlotSigMap(GCM_rcp45_JJA_end_Mean,lons_GCM,lats_GCM,ax,-12,12,title11,'(e)',cbar=0)
plt.sca(axes[2,1]) # active shis subplot for GCM
ax=axes[2,1]
PlotSigMap(GCM_rcp85_JJA_end_Mean,lons_GCM,lats_GCM,ax,-12,12,title22,'(f)',cbar=0)

plt.sca(axes[2,2]) # active shis subplot for GCM
ax=axes[2,2]
PlotSigMap(RCM_rcp45_JJA_end_Mean,lons_RCM,lats_RCM,ax,-12,12,title33,'(g)',cbar=0)
plt.sca(axes[2,3]) # active shis subplot for GCM
ax=axes[2,3]
PlotSigMap(RCM_rcp85_JJA_end_Mean,lons_RCM,lats_RCM,ax,-12,12,title44,'(h)',cbar=1)
#--------------------------------------------------- 

plt.sca(axes[3,0]) # active shis subplot for GCM
ax=axes[3,0]
PlotMap(GCM_rcp45_JJA_end_NB,lons_GCM,lats_GCM,ax,0,100,'(e)',cbar=0)
plt.sca(axes[3,1]) # active shis subplot for GCM
ax=axes[3,1]
PlotMap(GCM_rcp85_JJA_end_NB,lons_GCM,lats_GCM,ax,0,100,'(f)',cbar=0)

plt.sca(axes[3,2]) # active shis subplot for GCM
ax=axes[3,2]
PlotMap(RCM_rcp45_JJA_end_NB,lons_RCM,lats_RCM,ax,0,100,'(g)',cbar=0)
plt.sca(axes[3,3]) # active shis subplot for GCM
ax=axes[3,3]
PlotMap(RCM_rcp85_JJA_end_NB,lons_RCM,lats_RCM,ax,0,100,'(h)',cbar=1)
#--------------------------------------------------- 


plt.suptitle('model simulateed mean changes of SSR (W/m2)',fontsize=14)

print Output+'.eps'
plt.savefig(Output+'.eps',format='eps')
plt.show()

quit()
