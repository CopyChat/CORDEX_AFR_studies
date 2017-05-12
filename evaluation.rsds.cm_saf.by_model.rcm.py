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
import Taylor
import ctang

#=================================================== Definitions
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/validation_CM_SAF/'
OBS_Dir='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/'
N_model = 21
VAR ='rsds' # ,'tas','sfcWind') #,'PVpot')
OBSvar = 'SIS'
N_column = 2
N_row = 21
N_plot = N_column*N_row
#=================================================== test
##
#=================================================== end of test
# use CanESM2 instead of all not avail model:
RCM_name=(\
	'RCA4_v1',\
	'RCA4_v1',\
	'RCA4_v1',\
	'RCA4_v1',\
	'RCA4_v1',\
	'RCA4_v1',\
	'RCA4_v1',\
	'RCA4_v1',\
	'RCA4_v1',\
	'RCA4_v1',\
        \
	'CCLM4-8-17_v1',\
	'CCLM4-8-17_v1',\
	'CCLM4-8-17_v1',\
	'CCLM4-8-17_v1',\
	'HIRHAM5_v2',\
	'HIRHAM5_v2',\
	'RACMO22T_v1',\
	'RACMO22T_v1',\
	'REMO2009_v1',\
	'REMO2009_v1',\
	'REMO2009_v1')
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

GCM_name=(\
	'CanESM2',\
	'CNRM-CM5',\
	'CSIRO-Mk3-6-0',\
	'EC-EARTH',\
	'IPSL-CM5A-MR',\
	'MIROC5',\
	'HadGEM2-ES',\
	'MPI-ESM-LR',\
	'NorESM1-M',\
	'GFDL-ESM2M',\
        \
	'CNRM-CM5',\
	'EC-EARTH',\
	'HadGEM2-ES',\
	'MPI-ESM-LR',\
	'EC-EARTH',\
	'NorESM1-M',\
	'EC-EARTH',\
	'HadGEM2-ES',\
	'EC-EARTH',\
	'IPSL-CM5A-LR',\
	'M-MPI-ESM-LR')

#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

OBSfile=(\
        'SISmm.CDR.mon.mean.198301-200512.SA.timmean.nc',\
        'SISmm.CDR.mon.mean.198301-200512.SA.monmean.detrended.maskannual.timstd.nc',\
        'SISmm.CDR.mon.mean.198301-200512.SA.yearmean.detrended.masknoooon.timstd.nc')

OBS_remap=(\
        'SISmm.CDR.mon.mean.198301-200512.SA.timmean.remap.rcm.nc',\
        'SISmm.CDR.mon.mean.198301-200512.SA.monmean.detrended.maskannual.timstd.remap.rcm.nc',\
        'SISmm.CDR.mon.mean.198301-200512.SA.yearmean.detrended.masknoooon.timstd.remap.rcm.nc')

filefix=(\
        '.hist_rcp85.day.1983-2005.SA.timmean.nc',\
        '.hist_rcp85.day.1983-2005.SA.monmean.detrended.maskannual.timstd.nc',\
        '.hist_rcp85.day.1983-2005.SA.yearmean.detrended.masknoooon.timstd.nc')

# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf(\
        Data+VAR+'_AFR-44_'+RCM_Model[1]+filefix[0])


# Read lon,lat for OBS Plot the remap OBS, because time variability cannot be normalised by RCM in low resolution

lonsOBS,latsOBS=ctang.read_lonlat_netcdf(OBS_Dir+OBS_remap[0])



# Read Ensmean of timmean for CORDEX & OBS
timmean_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[0])\
        for i in range(N_model)])
timmean_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBSfile[0]))
timmean_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBS_remap[0]))

print timmean_OBS_remap
#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

OBSfile=(\
        'SISmm.CDR.mon.mean.198301-200512.SA.timmean.nc',\
        'SISmm.CDR.mon.mean.198301-200512.SA.monmean.detrended.maskannual.timstd.nc',\
        'SISmm.CDR.mon.mean.198301-200512.SA.yearmean.detrended.masknoooon.timstd.nc')

OBS_remap=(\
        'SISmm.CDR.mon.mean.198301-200512.SA.timmean.remap.rcm.nc',\
        'SISmm.CDR.mon.mean.198301-200512.SA.monmean.detrended.maskannual.timstd.remap.rcm.nc',\
        'SISmm.CDR.mon.mean.198301-200512.SA.yearmean.detrended.masknoooon.timstd.remap.rcm.nc')

filefix=(\
        '.hist_rcp85.day.1983-2005.SA.timmean.nc',\
        '.hist_rcp85.day.1983-2005.SA.monmean.detrended.maskannual.timstd.nc',\
        '.hist_rcp85.day.1983-2005.SA.yearmean.detrended.masknoooon.timstd.nc')

# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf(\
        Data+VAR+'_AFR-44_'+RCM_Model[1]+filefix[0])

# Read lon,lat for OBS Plot the remap OBS, because time variability cannot be normalised by RCM in low resolution

lonsOBS,latsOBS=ctang.read_lonlat_netcdf(OBS_Dir+OBS_remap[0])


# Read Ensmean of timmean for CMIP5 & OBS
timmean_CMIP5=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[0])\
        for i in range(N_model)])
timmean_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBSfile[0]))
timmean_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBS_remap[0]))

print timmean_OBS_remap


Bias=np.array([timmean_CMIP5[i]-timmean_OBS_remap for i in range(N_model)])
print Bias.shape

print timmean_OBS_remap.shape
print timmean_CMIP5.shape

#=================================================== plot
Title='Evaluation of the simulated SSR in the historical period 1983-2005'
#=================================================== 
Unit=( '(W/m2)','(W/m2)','(W/m2)')
#=================================================== 
TITLE2=('','bias vs CM_SAF')
def PlotMap(array2D,lons,lats,m,k,axx,vmin,vmax,cmap):
    # cmap = plt.cm.jet
    # cmap = plt.cm.bwr
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,11)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    map=Basemap(projection='cyl',llcrnrlat=lats[:,0].min(),urcrnrlat=lats[:,0].max(),\
            llcrnrlon=lons[0,:].min(),urcrnrlon=lons[0,:].max(),resolution='l')
    ctang.setMap(map)

    x,y=map(lats,lons)

    map.pcolormesh(y,x,array2D,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    #axx.axis('off')
    axx.xaxis.set_visible(False)
    axx.yaxis.set_visible(False)

    plt.title(str(m+1)+' '+GCM_name[m]+" "+RCM_name[m]+" "+ TITLE2[k]+" "+Unit[k],fontsize= 8)

    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=0.7) 
    cb.ax.tick_params(['{:.0f}'.format(x) for x in bounds ],labelsize=6) 
    #cbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in np.arange(cbar_min, cbar_max+cbar_step, cbar_step)], fontsize=16, weight='bold')
#=================================================== 
#=================================================== ploting
fig, axes = plt.subplots(nrows=N_row, ncols=N_column,\
        # sharex=True, sharey=True,\
        figsize=(6, 35),facecolor='w', edgecolor='k') # figsize=(w,h)
fig.subplots_adjust(hspace=0.3,top=0.96,wspace=0)
#=================================================== 
LIMIT=np.array([ [120,300],[-75,75]])

for m in range(N_row):
    if m == 0:
        for k in range(N_column):
            print 'm='+str(m),'k='+str(k)
            plt.sca(axes[m,k]) # active shis subplot for RCM
            axx=axes[m,k]
            if k == 0:
                cmap = plt.cm.jet
                PlotMap(timmean_CMIP5[m],lons,lats,m,k,axx,LIMIT[k][0],LIMIT[k][1],cmap)
            if k == 1:
                cmap = plt.cm.bwr
                PlotMap(Bias[m],lons,lats,m,k,axx,LIMIT[k][0],LIMIT[k][1],cmap)
    else:
        if RCM_Model[m] == 'CanESM2':
            ctang.NotAvailable(axes[m,0])
            ctang.NotAvailable(axes[m,1])
        else:
            for k in range(N_column):
                print 'm='+str(m),'k='+str(k)
                plt.sca(axes[m,k]) # active shis subplot for RCM
                axx=axes[m,k]
                if k == 0:
                    cmap = plt.cm.jet
                    PlotMap(timmean_CMIP5[m],lons,lats,m,k,axx,LIMIT[k][0],LIMIT[k][1],cmap)
                if k == 1:
                    cmap = plt.cm.bwr
                    PlotMap(Bias[m],lons,lats,m,k,axx,LIMIT[k][0],LIMIT[k][1],cmap)

#TaylorPlot(samples,refstd,fig,rect,ax4):


#=================================================== test

plt.suptitle(Title)

#plt.savefig('evaluation.eps',format='eps')
plt.savefig('evaluation.rsds.cm_saf.by_model.rcm.png')
plt.show()

quit()
