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

# commit
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
N_column = 4
N_row = 3
N_plot = N_column*N_row
#=================================================== test
##
#=================================================== end of test
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

# Read Ensmean of monthly std for CORDEX & OBS
monstd_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[1]) for i in range(N_model)])
monstd_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar,OBS_Dir+OBSfile[1]))
monstd_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(OBSvar,OBS_Dir+OBS_remap[1]))

# Read Ensmean of annual std for CORDEX & OBS
annualstd_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[2]) for i in range(N_model)])
annualstd_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBSfile[2]))
annualstd_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBS_remap[2]))

#=================================================== cal
# remove the missing points
timmean_OBS[timmean_OBS > 100055] = np.nan
monstd_OBS[monstd_OBS > 100000] = np.nan
annualstd_OBS[annualstd_OBS > 100000] = np.nan

timmean_OBS_remap[timmean_OBS_remap > 100055] = np.nan
monstd_OBS_remap[monstd_OBS_remap > 100000] = np.nan
annualstd_OBS_remap[annualstd_OBS_remap > 100000] = np.nan

timmean_CORDEX[timmean_CORDEX > 100055] = np.nan
monstd_CORDEX[monstd_CORDEX > 100000] = np.nan
annualstd_CORDEX[annualstd_CORDEX > 100000] = np.nan

# Ensmean of CORDEX 
Ensmean_timmean_CORDEX=np.mean(timmean_CORDEX,axis=0)
Ensmean_monstd_CORDEX=np.mean(monstd_CORDEX,axis=0)*100/Ensmean_timmean_CORDEX
Ensmean_annualstd_CORDEX=np.mean(annualstd_CORDEX,axis=0)*100/Ensmean_timmean_CORDEX

# for OBS
#monstd_OBS = monstd_OBS*100/Ensmean_timmean_CORDEX_remap
#annualstd_OBS = annualstd_OBS*100/Ensmean_timmean_CORDEX_remap

monstd_OBS_remap = monstd_OBS_remap*100/Ensmean_timmean_CORDEX
annualstd_OBS_remap = annualstd_OBS_remap*100/Ensmean_timmean_CORDEX


# Ensmean of Bias
Ensmean_timmean_Bias=(Ensmean_timmean_CORDEX-timmean_OBS_remap)*100/Ensmean_timmean_CORDEX
Ensmean_monstd_Bias=(Ensmean_monstd_CORDEX-monstd_OBS_remap)
Ensmean_annualstd_Bias=(Ensmean_annualstd_CORDEX-annualstd_OBS_remap)


Climatology=np.array([ Ensmean_timmean_CORDEX, Ensmean_monstd_CORDEX, Ensmean_annualstd_CORDEX])
OBSData=np.array([ timmean_OBS_remap, monstd_OBS_remap, annualstd_OBS_remap])
BiasData=np.array([ Ensmean_timmean_Bias, Ensmean_monstd_Bias, Ensmean_annualstd_Bias])

print OBSData[0]
print OBSData[1].shape
print OBSData[2].shape

#=================================================== input
## Reference dataset
RefStd = np.array([np.nanstd(var) for var in (\
        timmean_OBS_remap.flatten(),\
        (monstd_OBS_remap).flatten(),\
        (annualstd_OBS_remap).flatten())])

# to get std and corr: 1st, make nan = mask
Samples0 = np.array([[np.nanstd(m,ddof=1)/RefStd[0], \
        np.ma.corrcoef(timmean_OBS_remap.flatten(),m,allow_masked='Ture')[0,1]]\
        for m in [np.ma.array(timmean_CORDEX[i,:,:].flatten(), \
            mask=np.isnan(timmean_CORDEX[i,:,:].flatten())) \
            for i in range(N_model)]])

Samples1 = np.array([[np.nanstd(m,ddof=1)/RefStd[1], \
        np.ma.corrcoef(monstd_OBS_remap.flatten(),m,allow_masked='Ture')[0,1]] for m in \
        [np.ma.array(monstd_CORDEX[i,:,:].flatten(), \
            mask=np.isnan(monstd_CORDEX[i,:,:].flatten())) \
            for i in range(N_model)]])

Samples2 = np.array([[np.nanstd(m,ddof=1)/RefStd[2], \
        np.ma.corrcoef(annualstd_OBS_remap.flatten(),m,allow_masked='Ture')[0,1]] for m in \
        [np.ma.array(annualstd_CORDEX[i,:,:].flatten(), \
            mask=np.isnan(annualstd_CORDEX[i,:,:].flatten())) \
            for i in range(N_model)]])
SAMPLE=np.array([Samples0, Samples1, Samples2])

print RefStd
print SAMPLE[1]
#=================================================== end of cal
#=================================================== plot
Title='Evaluation of the simulated SSR in the historical period 1983-2005'
TTT1=('Climatology','OBS(CM_SAF','Bias','Taylor diagram')
TTT0=('Mean', 'Monthly variability', 'Annual variability')
#=================================================== 
Unit=(\
    ('(W/m2)','(W/m2)','(%)','(%)'),\
    ('(%)','(%)','(%)','(%)'),\
    ('(%)','(%)','(%)','(%)'))
#=================================================== 
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

    plt.title('SSR '+TTT0[m]+' '+TTT1[k]+' '+Unit[m][k],fontsize= 8)
    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=1) 
    cb.ax.tick_params(['{:.0f}'.format(x) for x in bounds ],labelsize=7) 
    #cbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in np.arange(cbar_min, cbar_max+cbar_step, cbar_step)], fontsize=16, weight='bold')
#=================================================== 
#=================================================== ploting
fig, axes = plt.subplots(nrows=N_row, ncols=N_column,\
        figsize=(30, 16),facecolor='w', edgecolor='k') # figsize=(w,h)
#fig.subplots_adjust(hspace=0.35,top=0.96,wspace=0.3)
#=================================================== 
LIMIT=np.array([\
        [[150,300],[150,300],[-45,45]],\
        [[4,12],[4,12],[-5,5]],\
        [[0,5],[0,5],[-2,2]]])

for m in range(N_row):
    for k in range(N_column):
        print 'm='+str(m),'k='+str(k)
        plt.sca(axes[m,k]) # active shis subplot for GCM
        axx=axes[m,k]
        if k == 0:
            if m == 0:
                PlotMap(Climatology[m],lons,lats,m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
            else:
                PlotMap(Climatology[m],lons,lats,m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
        if k == 1:
            if m == 0:
                PlotMap(OBSData[m],lonsOBS,latsOBS,m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
            else:
                PlotMap(OBSData[m],lonsOBS,latsOBS,m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
        if k == 2:
            PlotMap(BiasData[m],lonsOBS,latsOBS,m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
        if k == 3:
            Taylor.TaylorPlot(SAMPLE[m],1,fig,(N_row,N_column,m*N_column+k+1),axx)
#TaylorPlot(samples,refstd,fig,rect,ax4):




#=================================================== test
## Now adding the colorbar
#cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 
#cb = plt.colorbar(ax1, cax = cbaxes)  
#=================================================== end test
#ax1 = plt.subplot(11,2,22)
#ax1.axis('off')
#plt.colorbar(cax,cmap=plt.cm.bwr,orientation='horizontal',shrink=0.9) 

plt.suptitle(Title)

#plt.savefig('evaluation.eps',format='eps')
plt.savefig('evaluation.rsds.cm_saf.rcm.png')
plt.show()

quit()
