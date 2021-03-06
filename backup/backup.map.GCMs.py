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
sys.path.append('/Users/ctang/Code/My_Python_Code/')
import ctang

degree_sign= u'\N{DEGREE SIGN}'
#=================================================== Definitions
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/'
N_region = 9
N_model = 21
T=2.086 # nof=20 P<0.05
VAR ='rsds' # ,'tas','sfcWind') #,'PVpot')
N_var=len(VAR)+1

#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars

    #               region_1 region_2 ... region_N_region
    #  model_1
    #  model_2
    #  ...
    #  model_N_model
    #  ens
#--------------------------------------------------- 
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

#each map is mean_ref[0,0,:,:].shape

#=================================================== test
#quit()
#=================================================== test

# Read lon,lat
lons,lats=ctang.read_lonlat_netcdf(\
        Data+VAR+'_AFR-44_'+RCM_Model[0]+'.hist.day.1970-1999.SA.timmean.nc')

mean_ref=np.array([[ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_AFR-44_'+RCM_Model[i]+'.hist.day.1970-1999.SA.timmean.nc')] \
    for i in range(N_model)])[:,0]
mean_future=np.array([[ctang.read_lonlatmap_netcdf(VAR,\
    Data+VAR+'_AFR-44_'+RCM_Model[i]+'.rcp85.day.2070-2099.SA.timmean.nc')] \
    for i in range(N_model)])[:,0]

mean_ref_GCM=np.array([[ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_Amon_'+GCM_Model[i]+'.hist.1970-1999.SA.timmean.nc')] \
        for i in range(N_model)])[:,0]
mean_future_GCM=np.array([[ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_Amon_'+GCM_Model[i]+'.rcp85.2070-2099.SA.timmean.nc')] \
        for i in range(N_model)])[:,0]

lonsGCM=np.array([[ctang.read_lonlat_netcdf_1D(\
        Data+VAR+'_Amon_'+GCM_Model[i]+'.rcp85.2070-2099.SA.timmean.nc')[0]] \
        for i in range(N_model)])[:,0]
latsGCM=np.array([[ctang.read_lonlat_netcdf_1D(\
        Data+VAR+'_Amon_'+GCM_Model[i]+'.rcp85.2070-2099.SA.timmean.nc')[1]] \
        for i in range(N_model)])[:,0]

#print np.array(ctang.read_lonlat_netcdf_1D(\
        #Data+VAR+'_Amon_'+GCM_Model[0]+'.rcp85.2070-2099.SA.timmean.nc')).shape

#=================================================== end of read

# Calculate change  4D: len(VAR):N_model: lat : lon
mean_change=np.array((mean_future-mean_ref)*100/mean_ref)                 # RADS
mean_change_GCM=np.array((mean_future_GCM-mean_ref_GCM)*100/mean_ref_GCM) # RADS

print mean_change.shape
print mean_change_GCM.shape
#=================================================== plot
Title='Porjection Changes in RSDS(%) - RCP8.5'
#=================================================== 
# define the colormap
cmap = plt.cm.bwr
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
#cmaplist[0] = (.5,.5,.5,1.0)
#cmaplist[-1] = (.8,1.0,1.0,1.0)
# create the new map
#cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
# define the bins and normalize
bounds = np.linspace(-21,21,11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

#=================================================== 
fig, axes = plt.subplots(nrows=11, ncols=4,\
        figsize=(11, 20),facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.35,top=0.96,wspace=0.2)

#fig=plt.figure(figsize=(24,16))
gs = gridspec.GridSpec(11,4)
gs.update(wspace=0.1, hspace=0.1, left=0.1, right=0.4, bottom=0.1, top=0.9) 

#=================================================== 
def Plot_GCM(m):
    map=Basemap(projection='cyl',llcrnrlat=latsGCM[m].min(),urcrnrlat=latsGCM[m].max(),\
        llcrnrlon=lonsGCM[m].min(),urcrnrlon=lonsGCM[m].max(),resolution='h')
    ctang.setMap(map)
    x,y=map(latsGCM[m],lonsGCM[m])
    print latsGCM[m].shape, lonsGCM[m].shape
    print x.shape, y.shape
    print mean_change_GCM[m].shape
    map.pcolormesh(y,x,mean_change_GCM[m],cmap=cmap,\
            vmin=-20,vmax=20)
    #axes[j,0].set_xticks([])
    #axes[j,0].set_yticks([])
    #axes[j,0].xaxis.set_visible(False)
    plt.title(str(GCM_Model[m])+' '+GCM_resolution[m],fontsize= 6)
    #plt.annotate(str(GCM_Model[m])+' '+GCM_resolution[m], xy=(1, 0), fontsize=5,\
        #horizontalalignment='left', verticalalignment='bottom')

def Plot_RCM(m):
    map=Basemap(projection='cyl',llcrnrlat=lats.min(),urcrnrlat=lats.max(),\
        llcrnrlon=lons.min(),urcrnrlon=lons.max(),resolution='h')
    ctang.setMap(map)
    x,y=map(lats,lons)
    map.pcolormesh(y,x,mean_change[m],cmap=cmap,vmin=-20,vmax=20)
    #axes[m,1].set_xlabel('lon',fontsize=0.4)
    #axes[j,1].set_xticks([])
    #axes[j,1].set_yticks([])
    #axes[j,1].xaxis.set_visible(False)
    plt.title('Exp '+ str(m+1)+','+str(RCM_Model[m]),fontsize= 6)
    #plt.annotate('Exp '+ str(m+1)+','+str(RCM_Model[m]), xy=(1, 0), fontsize=5,\
        #horizontalalignment='left', verticalalignment='bottom')


def NotAvailable():
    t='NotAvailable'
    plt.axis([0, 10, 0, 10])
    plt.text(5, 10, t, fontsize=18, style='oblique', ha='center',\
                    va='center', wrap=True)
#=================================================== 

for m in range(11):
    for k in range(4):
        plt.sca(axes[m,k]) # active shis subplot for GCM
        if k == 0:
            if GCM_Model[m] == 'EC-EARTH':
                NotAvailable
            else:
                Plot_GCM(m)
        if k == 1:
            Plot_RCM(m)
        if k == 2:
            if m == 10:
                plt.colorbar(orientation='horizontal',shrink=1.2) # draw colorbar
            else:
                if GCM_Model[m+10] == 'EC-EARTH':
                    NotAvailable
                else:
                    Plot_GCM(m+10)
        if k == 3:
            if m == 10:
                print "done"
            else:
                Plot_RCM(m+10)





plt.suptitle(Title)

#plt.tight_layout()
#plt.savefig('map.gcm_vs_rcm.eps',format='eps')
#plt.savefig('map.gcm_vs_rcm.pdf')
#plt.savefig('map.gcm_vs_rcm.png')
plt.show()

quit()
