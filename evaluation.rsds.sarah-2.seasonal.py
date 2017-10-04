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
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/validation_SARAH_2/'
OBS_Dir='/Users/ctang/Code/CORDEX_AFR_studies/data/validation_SARAH_2/'
VAR ='rsds' # ,'tas','sfcWind') #,'PVpot')
OBS='SARAH_2'
OBSvar = 'SIS'
N_column = 2
Season='DJF'
Season='JJA'
#=================================================== test
##
#=================================================== end of test
# use CanESM2 instead of all not avail model:

# real model:
Model_name=(\
        'CanESM2',\
        'CNRM-CM5',\
        'CSIRO-Mk3-6-0',\
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
        # GCM:
        'rsds_Amon_CNRM-CM5_historical-rcp85_r1i1p1.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_Amon_CSIRO-Mk3-6-0_historical-rcp85_r1i1p1.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_Amon_CanESM2_historical-rcp85_r1i1p1.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_Amon_EC-EARTH_historical-rcp85_r1i1p1.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_Amon_GFDL-ESM2M_historical-rcp85_r1i1p1.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_Amon_HadGEM2-ES_historical-rcp85_r1i1p1.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_Amon_IPSL-CM5A-LR_historical-rcp85_r1i1p1.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_Amon_IPSL-CM5A-MR_historical-rcp85_r1i1p1.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_Amon_MIROC5_historical-rcp85_r1i1p1.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_Amon_MPI-ESM-LR_historical-rcp85_r1i1p1.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_Amon_NorESM1-M_historical-rcp85_r1i1p1.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # evaluation:
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_198901-200812.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_DMI-HIRHAM5_v2_mon_198901-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22T_v1_mon_197901-201212.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_MPI-CSC-REMO2009_v1_mon_198902-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_SMHI-RCA4_v1_mon_198001-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_UQAM-CRCM5_v1_mon_197901-201212.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # RCM:
        # 'rsds_AFR-44_CCCma-CanESM2_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_CNRM-CERFACS-CNRM-CM5_CLMcom-CCLM4-8-17_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_CNRM-CERFACS-CNRM-CM5_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_CSIRO-QCCCE-CSIRO-Mk3-6-0_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_198901-200812.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_DMI-HIRHAM5_v2_mon_198901-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22T_v1_mon_197901-201212.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_MPI-CSC-REMO2009_v1_mon_198902-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_SMHI-RCA4_v1_mon_198001-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_UQAM-CRCM5_v1_mon_197901-201212.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ICHEC-EC-EARTH_CLMcom-CCLM4-8-17_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ICHEC-EC-EARTH_DMI-HIRHAM5_v2.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ICHEC-EC-EARTH_KNMI-RACMO22T_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ICHEC-EC-EARTH_MPI-CSC-REMO2009_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_ICHEC-EC-EARTH_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_IPSL-IPSL-CM5A-LR_GERICS-REMO2009_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_IPSL-IPSL-CM5A-MR_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_MIROC-MIROC5_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_MOHC-HadGEM2-ES_CLMcom-CCLM4-8-17_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_MOHC-HadGEM2-ES_KNMI-RACMO22T_v2.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_MOHC-HadGEM2-ES_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_MPI-M-MPI-ESM-LR_MPI-CSC-REMO2009_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_MPI-M-MPI-ESM-LR_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_NCC-NorESM1-M_DMI-HIRHAM5_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_NCC-NorESM1-M_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # 'rsds_AFR-44_NOAA-GFDL-GFDL-ESM2M_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        )

GCM_name=GCM_Model
N_model = len(GCM_Model)
N_row = N_model
N_plot = N_column*N_row
LABLE='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
LABLE2='abcdefghijklmnopqrstuvwxyz'
#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

OBS_remap=(\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.CMIP5-ENSEMBLE.nc',\
'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.CNRM-CM5.nc',\
'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.CSIRO-Mk3-6-0.nc',\
'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.CanESM2.nc',\
'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.EC-EARTH.nc',\
'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.GFDL-ESM2M.nc',\
'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.HadGEM2-ES.nc',\
'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.IPSL-CM5A-LR.nc',\
'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.IPSL-CM5A-MR.nc',\
'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.MIROC5.nc',\
'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.MPI-ESM-LR.nc',\
'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.NorESM1-M.nc',\
# #
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# #
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
# 'SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc',\
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

#=================================================== plot
Title='Evaluation of the simulated SSR in the historical period 1983-2005 in '\
        +str(Season)
#=================================================== 
Unit=( '(W/m2)','(W/m2)','(W/m2)')
#=================================================== 
TITLE2=('','- '+str(OBS))
def PlotMap(array2D,lons,lats,m,k,axx,vmin,vmax,cmap):
    # cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,21)
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

    plt.title(GCM_name[m]+" "+TITLE2[k]+" "+Unit[k],fontsize= 8)

    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=0.8) 
    cb.ax.tick_params(['{:.0f}'.format(x) for x in bounds ],labelsize=6) 
    axx.text(0.9, 0.9,str(Season), ha='center', va='center', transform=axx.transAxes)
    if k==0:
        axx.text(0.1, 0.9,'('+str(LABLE[m])+')', ha='center', va='center', transform=axx.transAxes)
    else:
        axx.text(0.1, 0.9,'('+str(LABLE2[m])+')', ha='center', va='center', transform=axx.transAxes)

    #cbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in np.arange(cbar_min, cbar_max+cbar_step, cbar_step)], fontsize=16, weight='bold')
#=================================================== 
#=================================================== ploting
fig, axes = plt.subplots(nrows=N_row, ncols=N_column,\
        # sharex=True, sharey=True,\
        figsize=(9, 40),facecolor='w', edgecolor='k') # figsize=(w,h)
fig.subplots_adjust(hspace=0.2,top=0.96,wspace=0)
#=================================================== 
LIMIT=np.array([ [60,360],[-100,100]])

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

OutputImage='evaluation.'+str(VAR)+'.'+str(OBS)+'.'+str(Season)
#plt.savefig(OutputImage+'.eps',format='eps')
plt.savefig(OutputImage+'.png')
plt.show()

quit()
