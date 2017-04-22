#!/usr/bin/env python
"""
========
Ctang, A map of mean max and min of ensembles
        from CORDEX AFR-44, in Southern Africa
        Data was restored on titan
TODO:   1. calculate multiyear monthly mean for past 30 years: to defined seasons
        2. in future: to check if there's CLT shift in monthly scale
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
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/map.12mon.SA/'
N_model = 21
VAR ='clt' 
output='map.12month.meanchanges.rcm.clt'
#=================================================== reading data
RCM_Model=(\
    'clt_AFR-44_CCCma-CanESM2_SMHI-RCA4_v1',\
    'clt_AFR-44_CNRM-CERFACS-CNRM-CM5_CLMcom-CCLM4-8-17_v1',\
    'clt_AFR-44_CNRM-CERFACS-CNRM-CM5_SMHI-RCA4_v1',\
    'clt_AFR-44_CSIRO-QCCCE-CSIRO-Mk3-6-0_SMHI-RCA4_v1',\
    'clt_AFR-44_ICHEC-EC-EARTH_CLMcom-CCLM4-8-17_v1',\
    'clt_AFR-44_ICHEC-EC-EARTH_DMI-HIRHAM5_v2',\
    'clt_AFR-44_ICHEC-EC-EARTH_KNMI-RACMO22T_v1',\
    'clt_AFR-44_ICHEC-EC-EARTH_MPI-CSC-REMO2009_v1',\
    'clt_AFR-44_ICHEC-EC-EARTH_SMHI-RCA4_v1',\
    'clt_AFR-44_IPSL-IPSL-CM5A-LR_GERICS-REMO2009_v1',\
    'clt_AFR-44_IPSL-IPSL-CM5A-MR_SMHI-RCA4_v1',\
    'clt_AFR-44_MIROC-MIROC5_SMHI-RCA4_v1',\
    'clt_AFR-44_MOHC-HadGEM2-ES_CLMcom-CCLM4-8-17_v1',\
    'clt_AFR-44_MOHC-HadGEM2-ES_KNMI-RACMO22T_v2',\
    'clt_AFR-44_MOHC-HadGEM2-ES_SMHI-RCA4_v1',\
    'clt_AFR-44_MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17_v1',\
    'clt_AFR-44_MPI-M-MPI-ESM-LR_MPI-CSC-REMO2009_v1',\
    'clt_AFR-44_MPI-M-MPI-ESM-LR_SMHI-RCA4_v1',\
    'clt_AFR-44_NCC-NorESM1-M_DMI-HIRHAM5_v1',\
    'clt_AFR-44_NCC-NorESM1-M_SMHI-RCA4_v1',\
    'clt_AFR-44_NOAA-GFDL-GFDL-ESM2M_SMHI-RCA4_v1')

Monthstep=('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

# Read lon,lat
lons,lats=ctang.read_lonlat_netcdf(Data+RCM_Model[1]+'.changes.ymon.mean.SA.nc')

# read CLT in 3D: 12month,lat,lon
Changes = np.array(\
        [ctang.read_3D_netcdf(VAR,Data+RCM_Model[i]+'.changes.ymon.mean.SA.nc')\
        for i in range(N_model)])
print Changes.shape


MeanChanges = np.mean(Changes, axis=0)
print MeanChanges.shape

#=================================================== end of Reading
#--------------------------------------------------- 
# for t test in each point
#--------------------------------------------------- 
#t_value in N_model,lat,lon
def significant_map(t_value):
    count=0
    for mon in range(12):
        for lat in range(t_value.shape[2]):
            for lon in range(t_value.shape[3]):
                grid = t_value[:,mon,lat,lon]
                grid = grid[grid < 999]  # remove missing values
                # print grid

                print(mon,lat,lon,len(grid))
                if len(grid) < 1:
                    t_value[:,mon,lat,lon]=np.NaN
                    print m,mon,lat,lon
                else:
                    Sig = stats.ttest_1samp(grid,0)[0]
                    if np.abs(Sig) > ctang.get_T_value(len(grid)):
                        count+=1
                    else:
                        t_value[:,mon,lat,lon]=np.NaN
    return t_value
#--------------------------------------------------- 

# MeanChanges = np.nanmean(significant_map(Changes),axis=0)

# save to txt file to save time
# ctang.Save2mat(output,MeanChanges) 

# reading from mat file
MeanChanges=ctang.Loadmat(output+'.mat')

print MeanChanges.shape


#=================================================== plot setting
LIMIT=[[-20,20]]

CbarLabel='clt (%)'

Title=str(N_model)+' RCMs Simulated mean CLT changes (absolute %) in eath month'

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
fig,axes = plt.subplots(nrows=4, ncols=3,\
    figsize=(12, 12),facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.3,bottom=0.3)
#fig.subplots_adjust(0,0,1,1,0,0)
#gs = gridspec.GridSpec(11,4)
#gs.update(wspace=0.1, hspace=0.1, left=0.1, right=0.1, bottom=0.1, top=0.9) 
#=================================================== 
for s in range(4): # in seasons
    for m in range(3): # 3 month per season
        print 'figure num = ','season='+str(s),'month='+str(m)
        plt.sca(axes[s,m]) # active shis subplot for GCM
        axx=axes[s,m]
        # to calculate the index of 12 month 0-11
        # h=s+1;v=m+1; #month=(h-1)*3 + v # 1-12
        month=s*3 + m

        im=PlotMap(MeanChanges[month,:,:],month,axx,vmin=-10,vmax=10)
cbaxes = fig.add_axes([0.2, 0.2, 0.6, 0.02]) 
#                    [left, bottom, width, height]
cb = plt.colorbar(im, cax = cbaxes,orientation='horizontal')  
cb.ax.set_xlabel(str(CbarLabel))
plt.suptitle(Title)

# plt.savefig('map.12month.clt.'+str(f)+'.eps',format='eps')
plt.savefig('map.12month.meanchanges.rcm.clt'+'.png')

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
