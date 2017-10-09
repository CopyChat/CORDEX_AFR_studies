#!/usr/bin/env python
"""
========
Ctang, A map of mean max and min of ensembles
        from CORDEX AFR-44, in Southern Africa
        Data was restored on titan
========
"""
import math
import os
import os.path as path
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
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/validation_CLARA-A2/'
OBS_Dir='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/validation_CLARA-A2/'
VAR ='clt' # ,'tas','sfcWind') #,'PVpot')
OBS='CLARA-A2'
OBSvar = 'cfc'
N_column = 2
Season='JJA'
Season='DJF'
#=================================================== test
##
#=================================================== end of test
# use CanESM2 instead of all not avail model:

# real model:
Model_name=(\
        'ERAINT',\
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

GCM_Model=(\
        # 'CMIP5-ENSEMBLE',\
        'ERA_In.clt.mon.mean.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # GCM:
        'clt_Amon_CNRM-CM5_historical-rcp85_r1i1p1_1951-2099.nc.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'clt_Amon_CSIRO-Mk3-6-0_historical-rcp85_r1i1p1_1951-2099.nc.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'clt_Amon_CanESM2_historical-rcp85_r1i1p1_1951-2099.nc.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'clt_Amon_EC-EARTH_historical-rcp85_r1i1p1_1951-2099.nc.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'clt_Amon_GFDL-ESM2M_historical-rcp85_r1i1p1_1951-2099.nc.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'clt_Amon_HadGEM2-ES_historical-rcp85_r1i1p1_1951-2099.nc.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'clt_Amon_IPSL-CM5A-LR_historical-rcp85_r1i1p1_1951-2099.nc.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'clt_Amon_IPSL-CM5A-MR_historical-rcp85_r1i1p1_1951-2099.nc.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'clt_Amon_MIROC5_historical-rcp85_r1i1p1_1951-2099.nc.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'clt_Amon_MPI-ESM-LR_historical-rcp85_r1i1p1_1951-2099.nc.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'clt_Amon_NorESM1-M_historical-rcp85_r1i1p1_1951-2099.nc.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        )

GCM_name=Model_name
N_model = len(GCM_Model)
N_row = N_model
N_plot = N_column*N_row
LABLE='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
LABLE2='abcdefghijklmnopqrstuvwxyz'
#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

OBS_remap=(\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.rea.nc',\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.CNRM-CM5.nc',\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.CSIRO-Mk3-6-0.nc',\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.CanESM2.nc',\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.EC-EARTH.nc',\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.GFDL-ESM2M.nc',\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.HadGEM2-ES.nc',\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.IPSL-CM5A-LR.nc',\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.IPSL-CM5A-MR.nc',\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.MIROC5.nc',\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.MPI-ESM-LR.nc',\
    'CFCmm.CM_SAF.CLARA-A2.1990-2005.GL.SA.'+str(Season)+'.timmean.remap.NorESM1-M.nc',\
    )


# Read lon,lat for model
# lons,lats=ctang.read_lonlat_netcdf_1D(\
        # Data+VAR+'_Amon_'+GCM_Model[1]+filefix[0])
print  Data+GCM_Model[1]
lons=np.array([ctang.read_lon_netcdf_1D(\
        Data+GCM_Model[i]) for i in range(N_model)])
lats=np.array([ctang.read_lat_netcdf_1D(\
        Data+GCM_Model[i]) for i in range(N_model)])

# Read timmean for CMIP5 & OBS
timmean_CMIP5=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+GCM_Model[i]) for i in range(N_model)])

timmean_CMIP5[0]=timmean_CMIP5[0]*100

print OBS_Dir+OBS_remap[0]
# SISmm.SARAH-2.1990-2005.SA.DJF.timmean.remap.rcm.nc


# for f in range(N_model):
    # PATH=OBS_Dir+OBS_remap[f]
    # if os.path.isfile(PATH):
        # print f
    # else:
        # print "find it!",f
        # print PATH

print OBS_Dir+OBS_remap[1]
        
timmean_OBS_remap=np.array([ctang.read_lonlatmap_netcdf(OBSvar,\
        OBS_Dir+OBS_remap[i]) for i in range(N_model)])


# remove missing values in SARAH-2: -999.0
for i in range(len(GCM_Model)):
    jj=timmean_OBS_remap[i]
    jj[jj == -999] = np.NAN
    timmean_OBS_remap[i]=jj

print timmean_OBS_remap[0]
print timmean_OBS_remap[1]

print timmean_OBS_remap.shape

Bias=np.array([timmean_CMIP5[i]-timmean_OBS_remap[i] for i in range(N_model)])
print Bias.shape
print Bias[1].shape


#=================================================== 
# statistic values:

COR=np.zeros((N_model))
RMSE=np.zeros((N_model))
MeanBias=np.zeros((N_model))

for n in range(N_model):
    OBS_remap_1d=timmean_OBS_remap[n].ravel()
    OBS_nonnan=np.array([value for value in OBS_remap_1d if not math.isnan(value)])
    CMIP5_nonnan=np.array([timmean_CMIP5[n].ravel()[i] for i in range(len(OBS_remap_1d)) if not math.isnan(OBS_remap_1d[i])])
    print OBS_remap_1d.shape,OBS_nonnan.shape,CMIP5_nonnan.shape

    # MeanBias:
    MeanBias[n]=np.mean(CMIP5_nonnan-OBS_nonnan) 

    # spatial cor:
    COR[n]= np.ma.corrcoef(OBS_nonnan,CMIP5_nonnan)[0,1]

    # RMSE
    RMSE[n]=np.sqrt(np.mean((CMIP5_nonnan-OBS_nonnan)**2))

print COR,RMSE,MeanBias

#=================================================== end of statistics




#=================================================== plot
Title='Evaluation of the simulated CLT in the historical period 1990-2005 in '\
        +str(Season)
#=================================================== 
# Unit=( '($\mathregular{W/m^{2}}$)','($\mathregular{W/m^{2}}$)')
Unit=( '(%)','(%)')
#=================================================== 
TITLE2=('',' - '+str(OBS))
def PlotMap(array2D,lons,lats,m,k,axx,vmin,vmax,cmap):
    # cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,21)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    print lats.shape
    map=Basemap(projection='cyl',llcrnrlat=-39,urcrnrlat=-0.1,\
            # llcrnrlon=lons[0,:].min(),urcrnrlon=lons[0,:].max(),resolution='l')
            llcrnrlon=1,urcrnrlon=59,resolution='l')
    ctang.setMap(map)

    x,y=map(lats,lons)

    map.pcolormesh(y,x,array2D,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    #axx.axis('off')
    axx.xaxis.set_visible(False)
    axx.yaxis.set_visible(False)

    plt.title(GCM_name[m]+TITLE2[k]+" "+Unit[k],fontsize= 8)

    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=0.8) 
    cb.ax.tick_params(['{:.0f}'.format(x) for x in bounds ],labelsize=6) 
    axx.text(0.9, 0.9,str(Season), ha='center', va='center', transform=axx.transAxes)
    axx.text(0.95, 0.25,'r='+str("%.2f" % COR[m]), ha='right',va='center', transform=axx.transAxes)
    axx.text(0.95, 0.15,'MBE='+str("%.2f" % MeanBias[m]), ha='right',va='center', transform=axx.transAxes)
    axx.text(0.95, 0.05,'RMSE='+str("%.2f" % RMSE[m]), ha='right',va='center', transform=axx.transAxes)
    # if k==0:
        # axx.text(0.1, 0.9,'('+str(LABLE[m])+')', ha='center', va='center', transform=axx.transAxes)
    # else:
        # axx.text(0.1, 0.9,'('+str(LABLE2[m])+')', ha='center', va='center', transform=axx.transAxes)

    #cbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in np.arange(cbar_min, cbar_max+cbar_step, cbar_step)], fontsize=16, weight='bold')
#=================================================== 
#=================================================== ploting
fig, axes = plt.subplots(nrows=N_row, ncols=N_column,\
        # sharex=True, sharey=True,\
        figsize=(9, 40),facecolor='w', edgecolor='k') # figsize=(w,h)
fig.subplots_adjust(hspace=0.2,top=0.96,wspace=0)
#=================================================== 
LIMIT=np.array([ [0,100],[-50,50]])

for m in range(N_row):
    if m < 1:
        for k in range(N_column):
            print 'm='+str(m),'k='+str(k)
            plt.sca(axes[m,k]) # active shis subplot for GCM
            axx=axes[m,k]
            if k == 0:
                cmap = plt.cm.jet
                PlotMap(timmean_CMIP5[m],lons[m],lats[m],m,k,axx,LIMIT[k][0],LIMIT[k][1],cmap)
            if k == 1:
                cmap = plt.cm.seismic
                PlotMap(Bias[m],lons[m],lats[m],m,k,axx,LIMIT[k][0],LIMIT[k][1],cmap)
    else:
        if GCM_Model[m] == 'CanESM2222':
            ctang.NotAvailable(axes[m,0])
            ctang.NotAvailable(axes[m,1])
        else:
            for k in range(N_column):
                print 'm='+str(m),'k='+str(k)
                plt.sca(axes[m,k]) # active shis subplot for gcm
                axx=axes[m,k]
                if k == 0:
                    cmap = plt.cm.jet
                    PlotMap(timmean_CMIP5[m],lons[m],lats[m],m,k,axx,LIMIT[k][0],LIMIT[k][1],cmap)
                if k == 1:
                    cmap = plt.cm.seismic
                    PlotMap(Bias[m],lons[m],lats[m],m,k,axx,LIMIT[k][0],LIMIT[k][1],cmap)

#TaylorPlot(samples,refstd,fig,rect,ax4):


#=================================================== test

plt.suptitle(Title)

OutputImage='evaluation.'+str(VAR)+'.'+str(OBS)+'.gcm.'+str(Season)
plt.savefig(OutputImage+'.eps',format='eps')
plt.savefig(OutputImage+'.png')
plt.show()

quit()
