#!/usr/bin/env python
"""
========
Ctang, A map of mean max and min of ensembles
        from CORDEX AFR-44, in Southern Africa
        Data was restored on titan
TODO:   1. calculate multiyear monthly mean for past 30 years: to defined seasons
        2. in future: to check if there's RSDS shift in monthly scale
        3. to plot the multi-model results of 1) & 2), map of std/mean in % is needed
        4. this script will generate 4 figs, each including 12 maps
========
"""
import math
import pdb
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.basemap import Basemap 

# to load my functions
import sys 
sys.path.append('/Users/ctang/Code/Python/')
import ctang

#=================================================== Definitions
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/'
N_model = 21
N_fig=2
VAR ='rsds' 
degree_sign= u'\N{DEGREE SIGN}'
#=================================================== reading data
GCM_Model=(\
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
GCM_resolution=(\
	str('2.8*2.8')+degree_sign,\
	str('1.4*1.4')+degree_sign,\
	str('1.865*1.875')+degree_sign,\
	str('1.12*1.12')+degree_sign,\
	str('1.267*2.5')+degree_sign,\
	str('1.4*1.4')+degree_sign,\
	str('1.25*1.875')+degree_sign,\
	str('1.865*1.875')+degree_sign,\
	str('1.895*2.5')+degree_sign,\
	str('2.02*2.5')+degree_sign,\
	str('1.4*1.4')+degree_sign,\
	str('1.12*1.12')+degree_sign,\
	str('1.25*1.875')+degree_sign,\
	str('1.865*1.875')+degree_sign,\
	str('1.12*1.12')+degree_sign,\
	str('1.895*2.5')+degree_sign,\
	str('1.12*1.12')+degree_sign,\
	str('1.25*1.875')+degree_sign,\
	str('1.12*1.12')+degree_sign,\
	str('1.895*3.75')+degree_sign,\
	str('1.865*1.875')+degree_sign)

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
RCM_Name=(\
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
	'RACMO22T_v2',\
	'MPI-CSC-REMO2009_v1',\
	'GERICS-REMO2009_v1',\
	'MPI-CSC-REMO2009_v1')

Monthstep=('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
#=================================================== test
#=================================================== end test

# Read lon,lat
lons,lats=ctang.read_lonlat_netcdf(\
        Data+VAR+'_AFR-44_'+RCM_Model[0]+'.hist.day.1970-1999.SA.ymonmean.nc')

# read RSDS in 3D: 12month,lat,lon
mean_ref=np.array(ctang.read_3D_netcdf(VAR,\
        Data+'rsds.hist.1970-1999.ymonmean.SA.ensmean.nc'))
std_ref=np.array(ctang.read_3D_netcdf(VAR,\
        Data+'rsds.hist.1970-1999.ymonmean.SA.ensstd.nc'))
mean_future=np.array(ctang.read_3D_netcdf(VAR,\
        Data+'rsds.rcp85.2070-2099.ymonmean.SA.ensmean.nc'))
std_future=np.array(ctang.read_3D_netcdf(VAR,\
        Data+'rsds.rcp85.2070-2099.ymonmean.SA.ensstd.nc'))

Changes = mean_future - mean_ref

Ref=np.array([mean_ref,Changes,\
        mean_future,std_future*100/mean_ref])   # (4,12,90,137)


# read changes_GCM

print mean_ref.shape

#=================================================== end of Reading

#=================================================== plot setting
LIMIT=[[150,320],[-20,20],\
        [150,320],[0,30]]   # 2 for mean 2 for std in %

CbarLabel=(\
        'rsds (W/m2)',\
        'rsds (W/m2)',\
        'rsds (W/m2)',\
        'std (%)',\
        'rsds (W/m2)',\
        'std (%)')
Title=(\
        'RCMs Simulated RSDS (W/m2) in eath month',\
        'RCMs Simulated RSDS mean changes (W/m2) in eath month',\

        'RCMs Simulated RSDS mean(W/m2) in eath month in 2070-2099 - RCP8.5',\
        'RCMs Simulated RSDS std(%) in eath month in 2070-2099 - RCP8.5',\

        'GCM Simulated RSDS mean(W/m2) in eath month in 1970-1999',\
        'GCM Simulated RSDS std(%) in eath month in 2070-2099 - RCP8.5')

#=================================================== ploting
def PlotMap(array2D,month,axx,vmin,vmax):

    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,9)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    map=Basemap(projection='cyl',llcrnrlat=lats[:,0].min(),urcrnrlat=lats[:,0].max(),\
            llcrnrlon=lons[0,:].min(),urcrnrlon=lons[0,:].max(),resolution='l')
    ctang.setMap(map)

    x,y=map(lats,lons)
    im=map.pcolormesh(y,x,array2D,cmap=cmap,vmin=vmin,vmax=vmax)
    #axx.axis('off')

    axx.xaxis.set_visible(False)
    axx.yaxis.set_visible(False)
    plt.title(str(Monthstep[month]),fontsize= 12)
    return im
#--------------------------------------------------- 
for f in [1]:  # plot the corresponding var in Ref
    plt.figure(f)
    fig,axes = plt.subplots(nrows=4, ncols=3,\
        figsize=(12, 12),facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace=0.3,bottom=0.3)
    #fig.subplots_adjust(0,0,1,1,0,0)
    #gs = gridspec.GridSpec(11,4)
    #gs.update(wspace=0.1, hspace=0.1, left=0.1, right=0.1, bottom=0.1, top=0.9) 
#=================================================== 
    for s in range(4): # in seasons
        for m in range(3): # 3 month per season
            print 'figure num = '+str(f),'season='+str(s),'month='+str(m)
            plt.sca(axes[s,m]) # active shis subplot for GCM
            axx=axes[s,m]
            # to calculate the index of 12 month 0-11
            # h=s+1;v=m+1; #month=(h-1)*3 + v # 1-12
            month=s*3 + m
            im=PlotMap(Ref[f,month,:,:],month,axx,vmin=LIMIT[f][0],vmax=LIMIT[f][1])
    cbaxes = fig.add_axes([0.2, 0.2, 0.6, 0.02]) 
    #                    [left, bottom, width, height]
    cb = plt.colorbar(im, cax = cbaxes,orientation='horizontal')  
    cb.ax.set_xlabel(str(CbarLabel[f]))
    plt.suptitle(Title[f])

    # plt.savefig('map.12month.rsds.'+str(f)+'.eps',format='eps')
    plt.savefig('map.12month.rsds.'+str(f)+'.png')

#=================================================== 


#=================================================== test
## Now adding the colorbar
#cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 
#=================================================== end test
#ax1 = plt.subplot(11,2,22)
#ax1.axis('off')
#plt.colorbar(cax,cmap=plt.cm.bwr,orientation='horizontal',shrink=0.9) 


plt.show()
quit()
