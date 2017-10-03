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
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/validation_SARAH_2/'
OBS_Dir='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/validation_SARAH-2/'
VAR ='rsds' # ,'tas','sfcWind') #,'PVpot')
OBS='SARAH-2'
OBSvar = 'SIS'
N_column = 2
N_row = 22
N_plot = N_column*N_row
Season='DJF'
Season='JJA'
#=================================================== test
##
#=================================================== end of test
# use CanESM2 instead of all not avail model:
RCM_name=(\
	'RCA4',\
	'RCA4',\
	'RCA4',\
	'RCA4',\
	'RCA4',\
	'RCA4',\
	'RCA4',\
	'RCA4',\
	'RCA4',\
	'RCA4',\
        \
	'CCLM4',\
	'CCLM4',\
	'CCLM4',\
	'CCLM4',\
	'HIRHAM5',\
	'HIRHAM5',\
	'RACMO22T',\
	'RACMO22T',\
	'REMO2009',\
	'REMO2009',\
	'REMO2009',\
	'21RCMs',\
        )
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
	'MPI-M-MPI-ESM-LR_MPI-CSC-REMO2009_v1',\
        
	'CORDEX-ENSEMBLE_21RCMs',\
        )
RCM_name=(\
        # evaluation:
        'ERAINT-CCLM4',\
        'ERAINT-HIRHAM5',\
        'ERAINT-RACMO22T',\
        'ERAINT-REMO2009',\
        'ERAINT-RCA4',\
        # RCM:
        'CanESM2-RCA4',\
        'CNRM-CM5-CCLM4',\
        'CNRM-CM5-RCA4',\
        'CSIRO-Mk3-6-0-RCA4',\
        'EC-EARTH-CCLM4',\
        'EC-EARTH-HIRHAM5',\
        'EC-EARTH-RACMO22T',\
        'EC-EARTH-REMO2009',\
        'EC-EARTH-RCA4',\
        'IPSL-CM5A-LR-REMO2009',\
        'IPSL-CM5A-MR-RCA4',\
        'MIROC5-RCA4',\
        'HadGEM2-ES-CCLM4',\
        'HadGEM2-ES-RACMO22T',\
        'HadGEM2-ES-RCA4',\
        'MPI-ESM-LR-CCLM4',\
        'MPI-ESM-LR-REMO2009',\
        'MPI-ESM-LR-RCA4',\
        'NorESM1-M-HIRHAM5',\
        'NorESM1-M-RCA4',\
        'GFDL-ESM2M-RCA4',\
        )

RCM_Model=(\
        # 'ERA_In.ssrd.mon.mean.1990-2005.SA.JJA.timmean.nc',\
        'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_198901-200812.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_DMI-HIRHAM5_v2_mon_198901-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22T_v1_mon_197901-201212.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_MPI-CSC-REMO2009_v1_mon_198902-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_SMHI-RCA4_v1_mon_198001-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        # RCM:
        'rsds_AFR-44_CCCma-CanESM2_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_CNRM-CERFACS-CNRM-CM5_CLMcom-CCLM4-8-17_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_CNRM-CERFACS-CNRM-CM5_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_CSIRO-QCCCE-CSIRO-Mk3-6-0_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_CLMcom-CCLM4-8-17_v1_mon_198901-200812.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_DMI-HIRHAM5_v2_mon_198901-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22T_v1_mon_197901-201212.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_MPI-CSC-REMO2009_v1_mon_198902-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_SMHI-RCA4_v1_mon_198001-201012.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_UQAM-CRCM5_v1_mon_197901-201212.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ICHEC-EC-EARTH_CLMcom-CCLM4-8-17_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ICHEC-EC-EARTH_DMI-HIRHAM5_v2.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ICHEC-EC-EARTH_KNMI-RACMO22T_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ICHEC-EC-EARTH_MPI-CSC-REMO2009_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_ICHEC-EC-EARTH_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_IPSL-IPSL-CM5A-LR_GERICS-REMO2009_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_IPSL-IPSL-CM5A-MR_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_MIROC-MIROC5_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_MOHC-HadGEM2-ES_CLMcom-CCLM4-8-17_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_MOHC-HadGEM2-ES_KNMI-RACMO22T_v2.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_MOHC-HadGEM2-ES_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_MPI-M-MPI-ESM-LR_MPI-CSC-REMO2009_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_MPI-M-MPI-ESM-LR_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_NCC-NorESM1-M_DMI-HIRHAM5_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_NCC-NorESM1-M_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        'rsds_AFR-44_NOAA-GFDL-GFDL-ESM2M_SMHI-RCA4_v1.hist_rcp85.day.1990-2005.SA.'+str(Season)+'.timmean.nc',\
        )


OBS_remap='SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc'


#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

OBSfile='SISmm.SARAH-2.1990-2005.SA.'+str(Season)+'.timmean.remap.rcm.nc'



print Data+RCM_Model[0]

N_model = len(RCM_Model)
#=================================================== 

# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf(Data+RCM_Model[0])

# Read lon,lat for OBS Plot the remap OBS, because time variability cannot be normalised by RCM in low resolution
lonsOBS,latsOBS=ctang.read_lonlat_netcdf(OBS_Dir+OBS_remap)

# Read Ensmean of timmean for CORDEX & OBS
timmean_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+RCM_Model[i]) for i in range(N_model)])
timmean_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBSfile))
timmean_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(OBSvar, OBS_Dir+OBS_remap))

# remove missing values in SARAH-2: -999.0
timmean_OBS_remap[timmean_OBS_remap == -999] = np.NAN

print timmean_OBS_remap
Bias=np.array([timmean_CORDEX[i]-timmean_OBS_remap for i in range(N_model)])
print Bias.shape


# statistic values:

OBS_remap_1d=timmean_OBS_remap.ravel()

MeanBias=np.array([np.nanmean(timmean_CORDEX[uuu].ravel()-OBS_remap_1d) for uuu in range(N_model)])


COR=np.array([np.ma.corrcoef(timmean_CORDEX[uuu].ravel(),OBS_remap_1d)[0,1] for uuu in range(N_model)])

print COR
print timmean_CORDEX[0]
aa=timmean_CORDEX[0].ravel()
print aa.shape
print OBS_remap_1d.shape
OBS_nonnan=OBS_remap_1d[OBS_remap_1d != np.nan]
aa_nonnan=aa[OBS_remap_1d != np.nan]

print OBS_nonnan
print aa_nonnan

print np.correlate(aa_nonnan,OBS_nonnan)

print np.ma.corrcoef(timmean_CORDEX[0].ravel(),OBS_remap_1d)

# RMSE
RMSE=np.array([np.sqrt(np.mean((timmean_CORDEX[jj].ravel()-OBS_remap_1d)**2)) for jj in range(N_model)])

print RMSE
# error=InputValue-OBS
# square=error**2
# mean=np.mean(square,dtype=np.float64)
# RMSE = "%.2f" %  np.sqrt(mean)


print timmean_OBS_remap.shape
print timmean_CORDEX.shape

#=================================================== plot
Title='Evaluation of the simulated '+str(VAR)+' in the historical period 1990-2005'
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

    plt.title(str(m+1)+' '+RCM_name[m]+" "+ TITLE2[k]+" "+Unit[k],fontsize= 8)

    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=0.8) 
    cb.ax.tick_params(['{:.0f}'.format(x) for x in bounds ],labelsize=6) 
    #cbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in np.arange(cbar_min, cbar_max+cbar_step, cbar_step)], fontsize=16, weight='bold')
    
    axx.text(0.9, 0.9,str(Season), ha='center', va='center', transform=axx.transAxes)
    axx.text(0.9, 0.1,'MB='+str("%.2f" % MeanBias[m]), ha='right',va='center', transform=axx.transAxes)

#=================================================== 
#=================================================== ploting
fig, axes = plt.subplots(nrows=N_row, ncols=N_column,\
        # sharex=True, sharey=True,\
        figsize=(9, 80),facecolor='w', edgecolor='k') # figsize=(w,h)
fig.subplots_adjust(hspace=0.2,top=0.96,wspace=0)
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
