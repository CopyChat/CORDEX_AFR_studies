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
import Taylor
import ctang

#=================================================== Definitions
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/'
N_model = 21
VAR ='rsds' # ,'tas','sfcWind') #,'PVpot')
OBSvar = 'ssrd'
N_column = 4
N_row = 4
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
        'EAR40.ssrd.daymean.1970-1999.SA2.timmean.nc',\
        'EAR40.ssrd.daymean.1970-1999.SA2.detrended.maskannual.timstd.nc',\
        'EAR40.ssrd.daymean.1970-1999.SA2.monmean.detrended.maskannual.timstd.nc',\
        'EAR40.ssrd.daymean.1970-1999.SA2.yearmean.detrended.timstd.nc')
filefix=(\
        '.hist.day.1970-1999.SA.timmean.nc' ,\
        '.hist.day.1970-1999.SA.detrended.daymean.maskannual.timstd.nc',\
        '.hist.day.1970-1999.SA.monmean.detrended.maskannual.timstd.nc',\
        '.hist.day.1970-1999.SA.yearmean.detrended.masknoooon.timstd.nc')

filefix_Remap=(\
        '.hist.day.1970-1999.SA.timmean.remap.nc' ,\
        '.hist.day.1970-1999.SA.detrended.daymean.maskannual.timstd.remap.nc',\
        '.hist.day.1970-1999.SA.monmean.detrended.maskannual.timstd.remap.nc',\
        '.hist.day.1970-1999.SA.yearmean.detrended.masknoooon.timstd.remap.nc')
Remap_CORDEX=(\
        'rsds_AFR-44.hist.daymean.timmean.ensmean.remap.nc',\
        'rsds_AFR-44.hist.daymean.timstd.ensmean.remap.nc',\
        'rsds_AFR-44.hist.monmean.timstd.ensmean.remap.nc',\
        'rsds_AFR-44.hist.yearmean.timstd.ensmean.remap.nc')

# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf(\
        Data+VAR+'_AFR-44_'+RCM_Model[1]+'.hist.day.1970-1999.SA.timmean.nc')

# Read lon,lat for OBS
lonsOBS,latsOBS=ctang.read_lonlat_netcdf_1D(\
        Data+OBSfile[0])

print lonsOBS.shape
# Read Ensmean of timmean for CORDEX & OBS
timmean_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[0])\
        for i in range(N_model)])
timmean_CORDEX_remap=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix_Remap[0])\
        for i in range(N_model)])
timmean_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar,\
        Data+OBSfile[0]))

# Read Ensmean of daily std for CORDEX & OBS
dailystd_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[1]) for i in range(N_model)])
dailystd_CORDEX_remap=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix_Remap[1])\
        for i in range(N_model)])
dailystd_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar,\
        Data+OBSfile[1]))

# Read Ensmean of monthly std for CORDEX & OBS
monstd_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[2]) for i in range(N_model)])
monstd_CORDEX_remap=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix_Remap[2])\
        for i in range(N_model)])
monstd_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar,\
        Data+OBSfile[2]))

# Read Ensmean of annual std for CORDEX & OBS
annualstd_CORDEX=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix[3]) for i in range(N_model)])
annualstd_CORDEX_remap=np.array([ctang.read_lonlatmap_netcdf(VAR,\
        Data+VAR+'_AFR-44_'+RCM_Model[i]+filefix_Remap[3])\
        for i in range(N_model)])
annualstd_OBS=np.array(ctang.read_lonlatmap_netcdf(OBSvar,\
        Data+OBSfile[3]))
print annualstd_OBS.shape

#=================================================== cal
# remove the missing points
timmean_CORDEX_remap[timmean_CORDEX_remap > 100055] = np.nan
dailystd_CORDEX_remap[dailystd_CORDEX_remap > 100000] = np.nan
monstd_CORDEX_remap[monstd_CORDEX_remap > 100000] = np.nan
annualstd_CORDEX_remap[annualstd_CORDEX_remap > 100000] = np.nan

# Ensmean of CORDEX 
Ensmean_timmean_CORDEX=np.mean(timmean_CORDEX,axis=0)
Ensmean_dailystd_CORDEX=np.mean(dailystd_CORDEX,axis=0)*100/Ensmean_timmean_CORDEX
Ensmean_monstd_CORDEX=np.mean(monstd_CORDEX,axis=0)*100/Ensmean_timmean_CORDEX
Ensmean_annualstd_CORDEX=np.mean(annualstd_CORDEX,axis=0)*100/Ensmean_timmean_CORDEX

Ensmean_timmean_CORDEX_remap=np.mean(timmean_CORDEX_remap,axis=0)
Ensmean_dailystd_CORDEX_remap=np.mean( dailystd_CORDEX_remap,axis=0)
Ensmean_monstd_CORDEX_remap=np.mean( monstd_CORDEX_remap,axis=0)
Ensmean_annualstd_CORDEX_remap=np.mean( annualstd_CORDEX_remap,axis=0)

Ensmean_timmean_Bias=(Ensmean_timmean_CORDEX_remap-timmean_OBS)*100/Ensmean_timmean_CORDEX_remap
Ensmean_dailystd_Bias=(Ensmean_dailystd_CORDEX_remap-dailystd_OBS)*100/Ensmean_timmean_CORDEX_remap
Ensmean_monstd_Bias=(Ensmean_monstd_CORDEX_remap-monstd_OBS)*100/Ensmean_timmean_CORDEX_remap
Ensmean_annualstd_Bias=(Ensmean_annualstd_CORDEX_remap-annualstd_OBS)*100/Ensmean_timmean_CORDEX_remap

print Ensmean_timmean_Bias.shape

Climatology=np.array([ Ensmean_timmean_CORDEX, Ensmean_dailystd_CORDEX,\
                    Ensmean_monstd_CORDEX, Ensmean_annualstd_CORDEX])
OBSData=np.array([ timmean_OBS, dailystd_OBS, monstd_OBS, annualstd_OBS])
BiasData=np.array([ Ensmean_timmean_Bias, Ensmean_dailystd_Bias,\
                    Ensmean_monstd_Bias, Ensmean_annualstd_Bias])

print Climatology[2]
print Climatology[2].shape
print BiasData[2]
print BiasData[2].shape
#=================================================== input
## Reference dataset
x = np.linspace(0,4*np.pi,100)
data = np.sin(x)
refstd = data.std(ddof=1)           # Reference standard deviation
print refstd

RefStd = np.array([np.nanstd(var) for var in (\
        timmean_OBS.flatten(),\
        (dailystd_OBS/Ensmean_timmean_CORDEX_remap).flatten(),\
        (monstd_OBS/Ensmean_timmean_CORDEX_remap).flatten(),\
        (annualstd_OBS/Ensmean_timmean_CORDEX_remap).flatten())])
RefStd = np.array([ RefStd[0],RefStd[1]*100, RefStd[2]*100, RefStd[3]*100])

# to get std and corr: 1st, make nan = mask
Samples0 = np.array([[np.nanstd(m,ddof=1)/RefStd[0], \
        np.ma.corrcoef(timmean_OBS.flatten(),m,allow_masked='Ture')[0,1]]\
        for m in [np.ma.array(timmean_CORDEX_remap[i,:,:].flatten(), \
            mask=np.isnan(timmean_CORDEX_remap[i,:,:].flatten())) \
            for i in range(N_model)]])

# time variability are needed to be normalised by Ensmean_timmean_CORDEX_remap
Samples1 = np.array([[np.nanstd(m,ddof=1)*100/RefStd[1], \
        np.ma.corrcoef(dailystd_OBS.flatten(),m,allow_masked='Ture')[0,1]] for m in \
        [np.ma.array((dailystd_CORDEX_remap[i,:,:]/Ensmean_timmean_CORDEX_remap).flatten(), \
            mask=np.isnan((dailystd_CORDEX_remap[i,:,:]/Ensmean_timmean_CORDEX_remap).flatten())) \
            for i in range(N_model)]])

Samples2 = np.array([[np.nanstd(m,ddof=1)*100/RefStd[2], \
        np.ma.corrcoef(monstd_OBS.flatten(),m,allow_masked='Ture')[0,1]] for m in \
        [np.ma.array((monstd_CORDEX_remap[i,:,:]/Ensmean_timmean_CORDEX_remap).flatten(), \
            mask=np.isnan((monstd_CORDEX_remap[i,:,:]/Ensmean_timmean_CORDEX_remap).flatten())) \
            for i in range(N_model)]])

Samples3 = np.array([[np.nanstd(m,ddof=1)*100/RefStd[3], \
        np.ma.corrcoef(annualstd_OBS.flatten(),m,allow_masked='Ture')[0,1]] for m in \
        [np.ma.array((annualstd_CORDEX_remap[i,:,:]/Ensmean_timmean_CORDEX_remap).flatten(), \
            mask=np.isnan((annualstd_CORDEX_remap[i,:,:]/Ensmean_timmean_CORDEX_remap).flatten())) \
            for i in range(N_model)]])
SAMPLE=np.array([Samples0, Samples1, Samples2, Samples3])
print SAMPLE[0]
print RefStd[0]
#=================================================== end of cal
#=================================================== 
# Models
m1 = data + 0.2*np.random.randn(len(x))    # Model 1
m2 = 0.8*data + .1*np.random.randn(len(x)) # Model 2
m3 = np.sin(x-np.pi/10)                    # Model 3
m4 = np.sin(x-np.pi/12)                    # Model 3

# Compute stddev and correlation coefficient of models
#for i in range(N_model):
samples = np.array([ [m.std(ddof=1), np.corrcoef(data, m)[0,1]]
                    for m in (m1,m2,m3,m4)])
#=================================================== plot
Title='Evaluation of the simulated RSDS in the historical period'
TTT1=('Climatology','OBS','Bias','Taylor diagram')
TTT0=('Mean', 'Daily variability', 'Monthly variability', 'Annual variability')
#=================================================== 
Unit=(\
    ('(W/m2)','(%)','(%)','(%)'),\
    ('(%)','(%)','(%)','(%)'),\
    ('(%)','(%)','(%)','(%)'),\
    ('(%)','(%)','(%)','(%)'))
#=================================================== 
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
    #axx.axis('off')
    axx.xaxis.set_visible(False)
    axx.yaxis.set_visible(False)

    plt.title('RSDS '+TTT0[m]+' '+TTT1[k]+' '+Unit[m][k],fontsize= 8)
    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=1) 
    cb.ax.tick_params(labelsize=5) 
#=================================================== 
#=================================================== ploting
fig, axes = plt.subplots(nrows=N_row, ncols=N_column,\
        figsize=(30, 16),facecolor='w', edgecolor='k') # figsize=(w,h)
#fig.subplots_adjust(hspace=0.35,top=0.96,wspace=0.3)
#=================================================== 
LIMIT=np.array([\
        [[150,300],[150,300],[-20,20]],\
        [[10,60],[10,60],[0,30]],\
        [[0,30],[0,30],[-5,5]],\
        [[0,5],[0,50],[-20,-10]]])

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

#plt.tight_layout()
plt.savefig('evaluation.eps',format='eps')
plt.savefig('evaluation.pdf')
plt.savefig('evaluation.png')
plt.show()

quit()
