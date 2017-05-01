#!/usr/bin/env python
"""
========
Ctang, A map of mean max and min of ensembles
        from CORDEX AFR-44, in Southern Africa
        Data was restored on titan
========
"""
import math
import pdb
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
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/validation_CRU/'
OBS_Dir='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/'
N_model = 21
VAR ='clt'
OBSvar = 'cld'
N_column = 2
N_row = 21
N_plot = N_column*N_row
#=================================================== test
##
#=================================================== end of test
# use CanESM2 instead of all not avail model:
GCM_Model=(\
	'CanESM2',\
	'CNRM-CM5',\

	# 'CSIRO-Mk3-6-0',\
	'CanESM2',\

	# 'EC-EARTH',\
	'CanESM2',\
	# 'IPSL-CM5A-MR',\
	'CanESM2',\

	'MIROC5',\
	'HadGEM2-ES',\
	'MPI-ESM-LR',\
	'NorESM1-M',\
	'GFDL-ESM2M',\
	'CNRM-CM5',\

	# 'EC-EARTH',\
	'CanESM2',\

	'HadGEM2-ES',\
	'MPI-ESM-LR',\

	# 'EC-EARTH',\
	'CanESM2',\

	'NorESM1-M',\

	# 'EC-EARTH',\
	'CanESM2',\

	'HadGEM2-ES',\

	# 'EC-EARTH',\
	'CanESM2',\

	'IPSL-CM5A-LR',\
	'MPI-ESM-LR')

GCM_name=GCM_Model
# real model:
# GCM_Model=(\
        # 'CNRM-CM5',\
        # 'CSIRO-Mk3-6-0',\
        # 'CanESM2',\
        # 'GFDL-ESM2M',\
        # 'HadGEM2-ES',\
        # 'IPSL-CM5A-LR',\
        # 'IPSL-CM5A-MR',\
        # 'MIROC5',\
        # 'MPI-ESM-LR',\
        # 'NorESM1-M')
#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

OBSfile='cru_ts4.00.cld.1970-1999.SA.timmean.nc'

OBS_remap='cru_ts4.00.cld.1970-1999.SA.timmean.remap.gcm.nc'

filefix='_historical-rcp85_r1i1p1.1970-1999.SA.timmean.remap.nc'

print Data+VAR+'_Amon_'+GCM_Model[1]+filefix
# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf_1D(\
        Data+VAR+'_Amon_'+GCM_Model[1]+filefix)

# Read lon,lat for OBS Plot the remap OBS, because time variability cannot be normalised by GCM in low resolution

lonsOBS,latsOBS=ctang.read_lonlat_netcdf_1D(OBS_Dir+OBS_remap)

# Read Ensmean of timmean for CMIP5 & OBS
timmean_CMIP5=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_Amon_'+GCM_Model[i]+filefix)\
        for i in range(N_model)])
timmean_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBSfile))
timmean_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBS_remap))

print timmean_OBS_remap

Bias=np.array([timmean_CMIP5[i]-timmean_OBS_remap for i in range(N_model)])
print Bias.shape

print timmean_OBS_remap.shape
print timmean_CMIP5.shape

#=================================================== plot
Title='Evaluation of the simulated CLT in the historical period 1970 1999'
#=================================================== 
Unit=( '(%)','(%)','(%)')
#=================================================== 
TITLE2=('','bias vs CRU')
def PlotMap(array2D,lons,lats,m,k,axx,vmin,vmax):
    cmap = plt.cm.jet
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

    plt.title(str(m+1)+' '+GCM_name[m]+" "+TITLE2[k]+" "+Unit[k],fontsize= 8)

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
LIMIT=np.array([ [0,80],[-25,25]])

for m in range(N_row):
    if m == 0:
        for k in range(N_column):
            print 'm='+str(m),'k='+str(k)
            plt.sca(axes[m,k]) # active shis subplot for GCM
            axx=axes[m,k]
            if k == 0:
                PlotMap(timmean_CMIP5[m],lons,lats,m,k,axx,LIMIT[k][0],LIMIT[k][1])
            if k == 1:
                PlotMap(Bias[m],lons,lats,m,k,axx,LIMIT[k][0],LIMIT[k][1])
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
                    PlotMap(timmean_CMIP5[m],lons,lats,m,k,axx,LIMIT[k][0],LIMIT[k][1])
                if k == 1:
                    PlotMap(Bias[m],lons,lats,m,k,axx,LIMIT[k][0],LIMIT[k][1])

#TaylorPlot(samples,refstd,fig,rect,ax4):


#=================================================== test

plt.suptitle(Title)

#plt.savefig('evaluation.eps',format='eps')
plt.savefig('evaluation.rsds.cru.by_model.gcm.png')
plt.show()

quit()
