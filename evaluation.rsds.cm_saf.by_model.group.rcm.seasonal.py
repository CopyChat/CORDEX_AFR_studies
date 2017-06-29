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
VAR ='rsds' # ,'tas','sfcWind') #,'PVpot')
OBS='CM SAF'
OBSvar = 'SIS'
Season='JJA'
#=================================================== test
##
#=================================================== end of test
# use CanESM2 instead of all not avail model:
RCM_name=(\
	'RCA4',\
	'CCLM4',\
	'HIRHAM5',\
	'RACMO22T',\
	'REMO',\
	'CORDEX-ENSEMBLE_5RCMmean',\
        )
RCM_Model=RCM_name
N_model = len(RCM_Model)

N_column = 2
N_row = N_model
N_plot = N_column*N_row
#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

OBSfile=(\
        'SISmm.CDR.mon.mean.198301-200512.SA.timmean.remap.rcm.nc',\
        'SISmm.CDR.mon.mean.197901-200512.SA.monmean.detrended.maskannual.timstd.remap.gcm.nc',\
        'SISmm.CDR.mon.mean.197901-200512.SA.yearmean.detrended.masknoooon.timstd.remap.gcm.nc')

OBS_remap=(\
        'SISmm.CDR.mon.mean.198301-200512.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
        'SISmm.CDR.mon.mean.197901-200512.SA.monmean.detrended.maskannual.timstd.remap.gcm.nc',\
        'SISmm.CDR.mon.mean.197901-200512.SA.yearmean.detrended.masknoooon.timstd.remap.gcm.nc')

filefix=(\
        '.hist_rcp85.day.1983-2005.SA.'+str(Season)+'.timmean.nc',\
        '_historical-rcp85_r1i1p1.1970-2099.nc.1979-2005.SA.monmean.detrended.maskannual.timstd.remap.nc',\
        '_historical-rcp85_r1i1p1.1970-2099.nc.1979-2005.SA.yearmean.detrended.masknoooon.timstd.remap.nc')



# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf(\
        Data+VAR+'_AFR-44_'+RCM_Model[1]+filefix[0])
# lons=np.array([ctang.read_lon_netcdf_1D(\
        # Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[0])\
        # for i in range(N_model)])
# lats=np.array([ctang.read_lat_netcdf_1D(\
        # Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[0])\
        # for i in range(N_model)])

# Read lon,lat for OBS Plot the remap OBS, because time variability cannot be normalised by RCM in low resolution
lonsOBS,latsOBS=ctang.read_lonlat_netcdf(OBS_Dir+OBS_remap[0])

# Read Ensmean of timmean for CORDEX & OBS
timmean_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[0])\
        for i in range(N_model)])
timmean_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBSfile[0]))
timmean_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBS_remap[0]))

print timmean_OBS_remap
Bias=np.array([timmean_CORDEX[i]-timmean_OBS_remap for i in range(N_model)])
print Bias.shape

print timmean_OBS_remap.shape
print timmean_CORDEX.shape

#=================================================== plot
Title='Evaluation of the simulated '+str(VAR)+' in the historical period 1983-2005'
#=================================================== 
Unit=( '(W/m2)','(W/m2)','(W/m2)')
#=================================================== 
TITLE2=('',' - '+str(OBS))
def PlotMap(array2D,lons,lats,m,k,axx,vmin,vmax,cmap):
    # cmap = plt.cm.jet
    # cmap = plt.cm.seismic
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,21)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    map=Basemap(projection='cyl',llcrnrlat=lats[:,0].min(),urcrnrlat=lats[:,0].max(),\
            llcrnrlon=lons[0,:].min(),urcrnrlon=lons[0,:].max(),resolution='l')
    ctang.setMap(map)

    x,y=map(lats,lons)

    map.pcolormesh(y,x,array2D,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    #axx.axis('off')
    axx.xaxis.set_visible(False)
    axx.yaxis.set_visible(False)

    plt.title(RCM_name[m]+" mean "+ TITLE2[k]+" "+Unit[k],fontsize= 8)

    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=0.8) 
    cb.ax.tick_params(['{:.0f}'.format(x) for x in bounds ],labelsize=6) 
    #cbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in np.arange(cbar_min, cbar_max+cbar_step, cbar_step)], fontsize=16, weight='bold')
    
    axx.text(0.9, 0.9,str(Season), ha='center', va='center', transform=axx.transAxes)

#=================================================== 
#=================================================== ploting
fig, axes = plt.subplots(nrows=N_row, ncols=N_column,\
        # sharex=True, sharey=True,\
        figsize=(9, 21),facecolor='w', edgecolor='k') # figsize=(w,h)
fig.subplots_adjust(hspace=0.1,top=0.96,wspace=0.1)
#=================================================== 
LIMIT=np.array([ [60,360],[-100,100]])

for m in range(N_row):
    if m == 0:
        for k in range(N_column):
            print 'm='+str(m),'k='+str(k)
            plt.sca(axes[m,k]) # active shis subplot for RCM
            axx=axes[m,k]
            if k == 0:
                cmap = plt.cm.jet
                PlotMap(timmean_CORDEX[m],lons,lats,m,k,axx,LIMIT[k][0],LIMIT[k][1],cmap)
            if k == 1:
                cmap = plt.cm.seismic
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
                    PlotMap(timmean_CORDEX[m],lons,lats,m,k,axx,LIMIT[k][0],LIMIT[k][1],cmap)
                if k == 1:
                    cmap = plt.cm.seismic
                    PlotMap(Bias[m],lons,lats,m,k,axx,LIMIT[k][0],LIMIT[k][1],cmap)

#TaylorPlot(samples,refstd,fig,rect,ax4):


#=================================================== test

plt.suptitle(Title)

OutputImage='evaluation.'+str(VAR)+'.'+str(OBS)+'.by_model.rcm.'+str(Season)
#plt.savefig(OutputImage+'.eps',format='eps')
plt.savefig(OutputImage+'.png')
plt.show()

quit()
