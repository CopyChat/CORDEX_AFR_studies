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
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.basemap import Basemap , addcyclic
import textwrap


# to load my functions
import sys 
sys.path.append('/Users/ctang/Code/My_Python_Code/')
import Taylor
import ctang

#=================================================== Definitions
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/'
OBS_Dir='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/'
N_model = 21
VAR ='rsds' # ,'tas','sfcWind') #,'PVpot')
OBSvar = 'ssrd'
N_column = 3
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
# reading GEBA

OBSfile='/Users/ctang/climate/GLOBALDATA/OBSDATA/GEBAdata/GEBA_5year_Southern_Africa'
GEBA = pd.DataFrame.from_csv(OBSfile,index_col=[0])
print GEBA.yearmean
print GEBA.columns.values
print GEBA.lat

# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf(\
        Data+VAR+'_AFR-44_'+RCM_Model[1]+'.hist.day.1970-1999.SA.timmean.nc')

# get the index of model map, which covers the station
idx=np.zeros((2,GEBA.lon.count()))
for i in range(GEBA.lon.count()):
    idx[0,i] = (np.abs(lons[0,:]-GEBA.lon[i])).argmin()
    idx[1,i] = (np.abs(lats[:,0]-GEBA.lat[i])).argmin()


# 21 * 4 table: 21 models vs 4 vars
filefix=(\
        '.hist.day.1970-1999.SA.timmean.nc',\
        '.hist.day.1970-1999.SA.detrended.daymean.maskannual.timstd.nc',\
        '.hist.day.1970-1999.SA.monmean.detrended.maskannual.timstd.nc',\
        '.hist.day.1970-1999.SA.yearmean.detrended.masknoooon.timstd.nc')

# Read lon,lat for OBS
lonsOBS = np.array(GEBA.lon)
latsOBS = np.array(GEBA.lat)



# read RCMs model data
#=================================================== 
# Read Ensmean of timmean for CORDEX & OBS
timmean_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[0])\
        for i in range(N_model)])

# Read Ensmean of daily std for CORDEX & OBS
#dailystd_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        #Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[1]) for i in range(N_model)])

# Read Ensmean of monthly std for CORDEX & OBS
monstd_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[2]) for i in range(N_model)])

# Read Ensmean of annual std for CORDEX & OBS
annualstd_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[3]) for i in range(N_model)])

#=================================================== cal
# remove the missing points
#timmean_CORDEX_remap[timmean_CORDEX_remap > 100055] = np.nan
#dailystd_CORDEX_remap[dailystd_CORDEX_remap > 100000] = np.nan
#monstd_CORDEX_remap[monstd_CORDEX_remap > 100000] = np.nan
#annualstd_CORDEX_remap[annualstd_CORDEX_remap > 100000] = np.nan

# Ensmean of CORDEX 
Ensmean_timmean_CORDEX=np.mean(timmean_CORDEX,axis=0)
#Ensmean_dailystd_CORDEX=np.mean(dailystd_CORDEX,axis=0)*100/Ensmean_timmean_CORDEX
Ensmean_monstd_CORDEX=np.mean(monstd_CORDEX,axis=0)*100/Ensmean_timmean_CORDEX
Ensmean_annualstd_CORDEX=np.mean(annualstd_CORDEX,axis=0)*100/Ensmean_timmean_CORDEX

# CORDEX over the satation
timmean_CORDEX_station=np.zeros((N_model,GEBA.lon.count()))
monstd_CORDEX_station=np.zeros((N_model,GEBA.lon.count()))
annualstd_CORDEX_station=np.zeros((N_model,GEBA.lon.count()))
for n in range(N_model):
    for i in range(GEBA.lon.count()):
        timmean_CORDEX_station[n,i]=timmean_CORDEX[n,idx[1,i],idx[0,i]]
        monstd_CORDEX_station[n,i]=monstd_CORDEX[n,idx[1,i],idx[0,i]]
        annualstd_CORDEX_station[n,i]=annualstd_CORDEX[n,idx[1,i],idx[0,i]]

# for OBS
#dailystd_OBS = dailystd_OBS*100/Ensmean_timmean_CORDEX_remap
timmean_OBS= np.array(GEBA.yearmean)
monstd_OBS = np.array(GEBA.mon_std)*100/np.mean(timmean_CORDEX_station,axis=0)
annualstd_OBS = np.array(GEBA.year_std)*100/np.mean(timmean_CORDEX_station,axis=0)

# Ensmean of Bias
Ensmean_timmean_Bias=np.zeros((GEBA.lon.count()))
Ensmean_monstd_Bias=np.zeros((GEBA.lon.count()))
Ensmean_annualstd_Bias=np.zeros((GEBA.lon.count()))

for i in range(GEBA.lon.count()):
    Ensmean_timmean_Bias[i]=Ensmean_timmean_CORDEX[idx[1,i],idx[0,i]] - timmean_OBS[i]
    Ensmean_monstd_Bias[i]=Ensmean_monstd_CORDEX[idx[1,i],idx[0,i]] - monstd_OBS[i]
    Ensmean_annualstd_Bias[i]=Ensmean_annualstd_CORDEX[idx[1,i],idx[0,i]] - annualstd_OBS[i]


Climatology=np.array([ Ensmean_timmean_CORDEX,Ensmean_monstd_CORDEX, Ensmean_annualstd_CORDEX])
OBSData=np.array([ timmean_OBS, monstd_OBS, annualstd_OBS])
BiasData=np.array([ Ensmean_timmean_Bias, Ensmean_monstd_Bias, Ensmean_annualstd_Bias])

print BiasData[2]
print OBSData[1]
print Climatology[2]

print BiasData.shape

print OBSData.shape

#=================================================== input
## Reference dataset

RefStd = np.array([np.nanstd(var) for var in (\
        timmean_OBS.flatten(),\
        (monstd_OBS).flatten(),\
        (annualstd_OBS).flatten())])
print RefStd
#RefStd = np.array([ RefStd[0],RefStd[1]*100, RefStd[2]*100, RefStd[3]*100])


# to get std and corr: 1st, make nan = mask
Samples0 = np.array([[np.nanstd(m,ddof=1)/RefStd[0], \
        np.ma.corrcoef(timmean_OBS,m,allow_masked='Ture')[0,1]]\
        for m in [np.ma.array(timmean_CORDEX_station[i,:], \
            mask=np.isnan(timmean_CORDEX_station[i,:])) \
            for i in range(N_model)]])

# time variability are needed to be normalised by Ensmean_timmean_CORDEX_remap
Samples1 = np.array([[np.nanstd(m,ddof=1)/RefStd[1], \
        np.ma.corrcoef(monstd_OBS,m,allow_masked='Ture')[0,1]] for m in \
        [np.ma.array(monstd_CORDEX_station[i,:], \
            mask=np.isnan(monstd_CORDEX_station[i,:])) \
            for i in range(N_model)]])

Samples2 = np.array([[np.nanstd(m,ddof=1)/RefStd[2], \
        np.ma.corrcoef(annualstd_OBS.flatten(),m,allow_masked='Ture')[0,1]] for m in \
        [np.ma.array(annualstd_CORDEX_station[i,:].flatten(), \
            mask=np.isnan(annualstd_CORDEX_station[i,:].flatten())) \
            for i in range(N_model)]])
SAMPLE=np.array([Samples0, Samples1, Samples2 ])

print timmean_CORDEX_station.shape
print RefStd.shape
print RefStd
print SAMPLE[1].shape
#=================================================== end of cal
#=================================================== plot
Title='Evaluation of the simulated RSDS in the historical period'
TTT1=('Climatology','Bias','Taylor diagram')
TTT0=('Mean', 'Monthly variability', 'Annual variability')
#=================================================== 
Unit=(\
    ('(W/m2)','(%)','(%)'),\
    ('(%)','(%)','(%)'),\
    ('(%)','(%)','(%)'))
#=================================================== 
cmap = plt.cm.jet
def PlotScatter(x,y,c,m,k,axx,vmin,vmax):
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,11)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    map=Basemap(projection='cyl',\
        llcrnrlat=lats[:,0].min(),urcrnrlat=lats[:,0].max(),\
        llcrnrlon=lons[0,:].min(),urcrnrlon=lons[0,:].max(),resolution='h')
    ctang.setMap(map)

    sc=plt.scatter(x, y, c=c, cmap=cmap,norm=norm,edgecolors='black',zorder=2,vmin=vmin, vmax=vmax,\
            s=35)
    
    axx.xaxis.set_visible(False)
    axx.yaxis.set_visible(False)

    plt.title('RSDS '+TTT0[m]+' '+TTT1[k]+' '+Unit[m][k],fontsize= 8)
    #plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=1) 
    cb=plt.colorbar(sc,orientation='horizontal',shrink=1)
    cb.ax.tick_params(labelsize=7) 

def PlotMap(array2D,lons,lats,m,k,axx,vmin,vmax):
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,10)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    map=Basemap(projection='cyl',llcrnrlat=lats[:,0].min(),urcrnrlat=lats[:,0].max(),\
            llcrnrlon=lons[0,:].min(),urcrnrlon=lons[0,:].max(),resolution='l')
    ctang.setMap(map)

    x,y=map(lats,lons)
    map.pcolormesh(y,x,array2D,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    axx.xaxis.set_visible(False)
    axx.yaxis.set_visible(False)

    plt.title('RSDS '+TTT0[m]+' '+TTT1[k]+' '+Unit[m][k],fontsize= 8)
    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=1) 
    cb.ax.tick_params(labelsize=7) 
#=================================================== 
#=================================================== ploting
fig, axes = plt.subplots(nrows=N_row, ncols=N_column,\
        figsize=(30, 16),facecolor='w', edgecolor='k') # figsize=(w,h)
fig.subplots_adjust(left=0.12,bottom=0.3,right=0.9,hspace=0.15,top=0.94,wspace=0.2)
#=================================================== 
LIMIT=np.array([\
        [[150,300],[-30,30]],\
        [[4,12],[-15,5]],\
        [[0,5],[-2,2]]])

for m in range(N_row):
    for k in range(N_column):
        print 'm='+str(m),'k='+str(k)
        plt.sca(axes[m,k]) # active shis subplot for GCM
        axx=axes[m,k]
        if k == 0:
            axx=axes[m,k]
            # mean
            PlotMap(Climatology[m],lons,lats,m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
            axx.scatter(lonsOBS, latsOBS, c=OBSData[m], \
                edgecolors='b',zorder=2,\
                vmin=LIMIT[m,k][0], vmax=LIMIT[m,k][1], s=35, cmap=cmap)
        if k == 1:
            # plot bias
            PlotScatter(lonsOBS,latsOBS,BiasData[m],m,k,axx,LIMIT[m,k][0],LIMIT[m,k][1])
        if k == 2:
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

#plt.tight_layout()
#plt.savefig('evaluation.eps',format='eps')
#plt.savefig('evaluation.pdf')
plt.savefig('evaluation.rsds.station.png')
plt.show()

quit()
