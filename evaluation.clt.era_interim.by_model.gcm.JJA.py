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
import pdb
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
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/validation_ERA_Interim/'
OBS_Dir='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/'
N_model = 21
VAR ='clt' # ,'tas','sfcWind') #,'PVpot')
OBS='ERA_Interim'
OBSvar = 'tcc'
N_column = 2
N_row = 21
N_plot = N_column*N_row
Season='JJA'
#=================================================== test
##
#=================================================== end of test
# use CanESM2 instead of all not avail model:
GCM_Model=(\
	'CanESM2',\
	'CNRM-CM5',\

	# 'CSIRO-Mk3-6-0',\
	'CanESM2',\

        'EC-EARTH',\
        'IPSL-CM5A-MR',\
	'MIROC5',\
	'HadGEM2-ES',\
	'MPI-ESM-LR',\
	'NorESM1-M',\
	'GFDL-ESM2M',\
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
	'MPI-ESM-LR')

GCM_name=GCM_Model
# real model:

# GCM_Model=(\
        # 'CNRM-CM5',\
        # 'CanESM2',\
        # 'GFDL-ESM2M',\
        # 'HadGEM2-ES',\
        # 'IPSL-CM5A-LR',\
        # 'MIROC5',\
        # 'MPI-ESM-LR',\
        # 'NorESM1-M')
#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

OBSfile=(\
        'ERA_In.clt.mon.mean.1979-2005.SA.'+str(Season)+'.timmean.nc',\
        'SISmm.CDR.mon.mean.197901-200512.SA.monmean.detrended.maskannual.timstd.remap.gcm.nc',\
        'SISmm.CDR.mon.mean.197901-200512.SA.yearmean.detrended.masknoooon.timstd.remap.gcm.nc')

OBS_remap=(\
        'ERA_In.clt.mon.mean.1979-2005.SA.'+str(Season)+'.timmean.remap.gcm.nc',\
        'SISmm.CDR.mon.mean.197901-200512.SA.monmean.detrended.maskannual.timstd.remap.gcm.nc',\
        'SISmm.CDR.mon.mean.197901-200512.SA.yearmean.detrended.masknoooon.timstd.remap.gcm.nc')

filefix=(\
        '_historical-rcp85_r1i1p1_1951-2099.nc.1979-2005.SA.'+str(Season)+'.timmean.nc',\
        '_historical-rcp85_r1i1p1.1970-2099.nc.1979-2005.SA.monmean.detrended.maskannual.timstd.remap.nc',\
        '_historical-rcp85_r1i1p1.1970-2099.nc.1979-2005.SA.yearmean.detrended.masknoooon.timstd.remap.nc')

filefix_remap=(\
        '_historical-rcp85_r1i1p1_1951-2099.nc.1979-2005.SA.'+str(Season)+'.timmean.remap.nc',\
        '_historical-rcp85_r1i1p1.1970-2099.nc.1979-2005.SA.monmean.detrended.maskannual.timstd.remap.nc',\
        '_historical-rcp85_r1i1p1.1970-2099.nc.1979-2005.SA.yearmean.detrended.masknoooon.timstd.remap.nc')

# Read lon,lat for model
# lons,lats=ctang.read_lonlat_netcdf_1D(\
        # Data+VAR+'_Amon_'+GCM_Model[1]+filefix[0])
lons=np.array([ctang.read_lon_netcdf_1D(\
        Data+VAR+'_Amon_'+GCM_Model[i]+filefix[0])\
        for i in range(N_model)])
lats=np.array([ctang.read_lat_netcdf_1D(\
        Data+VAR+'_Amon_'+GCM_Model[i]+filefix[0])\
        for i in range(N_model)])

# Read lon,lat for OBS Plot the remap OBS, because time variability cannot be normalised by GCM in low resolution

lonsOBS,latsOBS=ctang.read_lonlat_netcdf_1D(OBS_Dir+OBS_remap[0])

# Read Ensmean of timmean for CMIP5 & OBS
timmean_CMIP5=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_Amon_'+GCM_Model[i]+filefix[0])\
        for i in range(N_model)])

timmean_CMIP5_remap=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_Amon_'+GCM_Model[i]+filefix_remap[0])\
        for i in range(N_model)])

timmean_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBSfile[0])*100)
timmean_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBS_remap[0])*100)

print timmean_OBS_remap

Bias=np.array([timmean_CMIP5_remap[i]-timmean_OBS_remap for i in range(N_model)])
print Bias.shape

print timmean_OBS_remap.shape
print timmean_CMIP5.shape

#=================================================== plot
Title='Evaluation of the simulated '+str(VAR)+' in the historical period 1979-2005'
#=================================================== 
Unit=( '(%)','(%)','(%)')
#=================================================== 
TITLE2=('','- '+str(OBS))
def PlotMap(array2D,lons,lats,m,k,axx,vmin,vmax):
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,11)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    print lats.shape
    map=Basemap(projection='cyl',llcrnrlat=lats[:,0].min(),urcrnrlat=lats[:,0].max(),\
            llcrnrlon=lons[0,:].min(),urcrnrlon=lons[0,:].max(),resolution='l')
    ctang.setMap(map)

    x,y=map(lats,lons)

    map.pcolormesh(y,x,array2D,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    #axx.axis('off')
    axx.xaxis.set_visible(False)
    axx.yaxis.set_visible(False)

    plt.title(str(m+1)+' '+GCM_name[m]+" "+TITLE2[k]+" "+Unit[k],fontsize= 8)

    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=0.5) 
    cb.ax.tick_params(['{:.0f}'.format(x) for x in bounds ],labelsize=6) 
    axx.text(0.9, 0.9,str(Season), ha='center', va='center', transform=axx.transAxes)
    #cbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in np.arange(cbar_min, cbar_max+cbar_step, cbar_step)], fontsize=16, weight='bold')
#=================================================== 
#=================================================== ploting
fig, axes = plt.subplots(nrows=N_row, ncols=N_column,\
        # sharex=True, sharey=True,\
        figsize=(8, 40),facecolor='w', edgecolor='k') # figsize=(w,h)
fig.subplots_adjust(hspace=0.3,top=0.96,wspace=0)
#=================================================== 
LIMIT=np.array([ [0,100],[-50,50]])

for m in range(N_row):
    if m == 0:
        for k in range(N_column):
            print 'm='+str(m),'k='+str(k)
            plt.sca(axes[m,k]) # active shis subplot for GCM
            axx=axes[m,k]
            if k == 0:
                PlotMap(timmean_CMIP5[m],lons[m],lats[m],m,k,axx,LIMIT[k][0],LIMIT[k][1])
            if k == 1:
                PlotMap(Bias[m],lonsOBS,latsOBS,m,k,axx,LIMIT[k][0],LIMIT[k][1])
    else:
        if GCM_Model[m] == 'CanESM2':
            ctang.NotAvailable(axes[m,0])
            ctang.NotAvailable(axes[m,1])
        else:
            for k in range(N_column):
                print 'm='+str(m),'k='+str(k)
                plt.sca(axes[m,k]) # active shis subplot for GCM
                axx=axes[m,k]
                if k == 0:
                    PlotMap(timmean_CMIP5[m],lons[m],lats[m],m,k,axx,LIMIT[k][0],LIMIT[k][1])
                if k == 1:
                    PlotMap(Bias[m],lonsOBS,latsOBS,m,k,axx,LIMIT[k][0],LIMIT[k][1])

#TaylorPlot(samples,refstd,fig,rect,ax4):


#=================================================== test

plt.suptitle(Title)

OutputImage='evaluation.'+str(VAR)+'.'+str(OBS)+'.by_model.gcm'+str(Season)
#plt.savefig(OutputImage+'.eps',format='eps')
plt.savefig(OutputImage+'.png')
plt.show()

quit()
