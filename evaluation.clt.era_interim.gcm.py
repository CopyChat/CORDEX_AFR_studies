#!/usr/bin/env python
"""
========
Ctang, A map of mean max and min of ensembles
        from CMIP5 AFR-44, in Southern Africa
        Data was restored on titan
========
"""
import math
import pdb
import pandas as pd
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
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/validation_ERA_Interim/'
OBS_Dir='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/'
N_model = 8
VAR ='clt' 
OBSvar = 'tcc'
N_column = 4
N_row = 3
N_plot = N_column*N_row
#=================================================== test
##
#=================================================== end of test
GCM_Model=(\
        'CNRM-CM5',\
        'CanESM2',\
        'GFDL-ESM2M',\
        'HadGEM2-ES',\
        'IPSL-CM5A-LR',\
        'MIROC5',\
        'MPI-ESM-LR',\
        'NorESM1-M')

#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

OBSfile=(\
        'ERA_In.clt.mon.mean.1979-2005.SA.timmean.nc',\
        'ERA_In.clt.mon.mean.1979-2005.SA.monmean.detrend.maskannual.timstd.nc',\
        'ERA_In.clt.mon.mean.1979-2005.SA.yearmean.detrend.masknoooon.timstd.nc')

OBS_remap=(\
        'ERA_In.clt.mon.mean.1979-2005.SA.timmean.remap.gcm.nc',\
        'ERA_In.clt.mon.mean.1979-2005.SA.monmean.detrend.maskannual.timstd.remap.gcm.nc',\
        'ERA_In.clt.mon.mean.1979-2005.SA.yearmean.detrend.masknoooon.timstd.remap.gcm.nc')

filefix=(\
        '_historical-rcp85_r1i1p1_1951-2099.nc.1979-2005.SA.timmean.remap.nc',\
        '_historical-rcp85_r1i1p1_1951-2099.nc.1979-2005.SA.yearmean.detrended.masknoooon.timstd.remap.nc',\
        '_historical-rcp85_r1i1p1_1951-2099.nc.1979-2005.SA.yearmean.detrended.masknoooon.timstd.remap.nc')
print Data+VAR+'_Amon_'+GCM_Model[1]+filefix[0]
# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf_1D(\
        Data+VAR+'_Amon_'+GCM_Model[1]+filefix[0])

# Read lon,lat for OBS Plot the remap OBS, because time variability cannot be normalised by GCM in low resolution

lonsOBS,latsOBS=ctang.read_lonlat_netcdf_1D(OBS_Dir+OBSfile[0])

# Read Ensmean of timmean for CMIP5 & OBS
timmean_CMIP5=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_Amon_'+GCM_Model[i]+filefix[0])\
        for i in range(N_model)])
timmean_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBSfile[0])*100)
timmean_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBS_remap[0])*100)

print timmean_OBS_remap
print Data+VAR+'_Amon_'+GCM_Model[1]+filefix[1]

# Read Ensmean of monthly std for CMIP5 & OBS
monstd_CMIP5=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_Amon_'+GCM_Model[i]+filefix[1]) for i in range(N_model)])

monstd_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar,OBS_Dir+OBSfile[1]))
monstd_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(OBSvar,OBS_Dir+OBS_remap[1]))

monstd_OBS=monstd_OBS*100
monstd_OBS_remap=monstd_OBS_remap*100


# Read Ensmean of annual std for CMIP5 & OBS
annualstd_CMIP5=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_Amon_'+GCM_Model[i]+filefix[2]) for i in range(N_model)])
annualstd_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBSfile[2]))
annualstd_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBS_remap[2]))

annualstd_OBS=annualstd_OBS*100
annualstd_OBS_remap=annualstd_OBS_remap*100
#=================================================== cal
# remove the missing points
timmean_OBS[timmean_OBS < -100] = np.nan
monstd_OBS[monstd_OBS < -100] = np.nan
annualstd_OBS[annualstd_OBS < -100] = np.nan

timmean_OBS_remap[timmean_OBS_remap < -100] = np.nan
monstd_OBS_remap[monstd_OBS_remap < -100] = np.nan
annualstd_OBS_remap[annualstd_OBS_remap < -100] = np.nan

timmean_CMIP5[timmean_CMIP5 > 1000] = np.nan
monstd_CMIP5[monstd_CMIP5 > 1000] = np.nan
annualstd_CMIP5[annualstd_CMIP5 > 1000] = np.nan

# Ensmean of CMIP5 
Ensmean_timmean_CMIP5=np.mean(timmean_CMIP5,axis=0)
Ensmean_monstd_CMIP5=np.mean(monstd_CMIP5,axis=0)
Ensmean_annualstd_CMIP5=np.mean(annualstd_CMIP5,axis=0)

# for OBS
#monstd_OBS = monstd_OBS
#annualstd_OBS = annualstd_OBS

# monstd_OBS_remap = monstd_OBS_remap
# annualstd_OBS_remap = annualstd_OBS_remap


# Ensmean of Bias
Ensmean_timmean_Bias=(Ensmean_timmean_CMIP5-timmean_OBS_remap)
Ensmean_monstd_Bias=(Ensmean_monstd_CMIP5-monstd_OBS_remap)
Ensmean_annualstd_Bias=(Ensmean_annualstd_CMIP5-annualstd_OBS_remap)


Climatology=np.array([ Ensmean_timmean_CMIP5, Ensmean_monstd_CMIP5, Ensmean_annualstd_CMIP5])
OBSData=np.array([ timmean_OBS, monstd_OBS, annualstd_OBS])
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

#=================================================== 
# output: 21 model x ( normalized std, corr)
def cal_corr_std(obs2D,N_model,modelarray):
    a=obs2D.flatten()
    for m in range(N_model):
        a=np.vstack((a, modelarray[m].flatten()))
    b=a.T

    df = pd.DataFrame(b)

    std1 = df.std()
    std0 = df.std()[0]
    print std0
    corr1 = df.corr()[0]

    return np.vstack((np.array(std1[1:]/std0),\
            np.array(corr1[1:]))).T
#=================================================== 

Samples0 = cal_corr_std(timmean_OBS_remap,N_model,timmean_CMIP5)
Samples1 = cal_corr_std(monstd_OBS_remap,N_model,monstd_CMIP5)
Samples2 = cal_corr_std(annualstd_OBS_remap,N_model,annualstd_CMIP5)

SAMPLE=np.array([Samples0, Samples1, Samples2])


print RefStd
print SAMPLE[0]
print SAMPLE[1]
print SAMPLE[2]
#=================================================== end of cal
#=================================================== plot
Title='Evaluation of the simulated CLT vs ERA_Interim in the historical period 1979-2005'
TTT1=('Climatology','OBS(ERA_Interim','Bias','Taylor diagram')
TTT0=('Mean', 'Monthly variability', 'Annual variability')
#=================================================== 
Unit=(\
    ('(abs %)','(abs %)','(abs %)','(%)'),\
    ('(abs %)','(abs %)','(abs %)','(abs %)'),\
    ('(abs %)','(abs %)','(abs %)','(abs %)'),\
    ('(%)','(%)','(%)','(abs %)'),\
    ('(W/m2)','(W/m2)','(W/m2)','(%)'))
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

    plt.title('CLT '+TTT0[m]+' '+TTT1[k]+' '+Unit[m][k],fontsize= 8)
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
        [[10,80],[10,80],[-25,25]],\
        [[0,15],[0,15],[-10,10]],\
        [[0,5],[0,5],[-5,5]]])

for m in range(N_row):
    for k in range(N_column):
        print 'm='+str(m),'k='+str(k)
        plt.sca(axes[m,k]) # active shis subplot for GCM
        axx=axes[m,k]
        if k == 0:
            if m == 0:
                print lons.shape
                PlotMap(Climatology[m],lons,lats,m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
            else:
                PlotMap(Climatology[m],lons,lats,m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
        if k == 1:
            if m == 0:
                PlotMap(OBSData[m],lonsOBS,latsOBS,m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
            else:
                PlotMap(OBSData[m],lonsOBS,latsOBS,m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
        if k == 2:
            PlotMap(BiasData[m],lons,lats,m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
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
plt.savefig('evaluation.clt.era_in.gcm.png')
plt.show()

quit()
