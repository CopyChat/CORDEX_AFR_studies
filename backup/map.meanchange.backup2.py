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
from scipy import stats
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.basemap import Basemap , addcyclic
import textwrap

# to load my functions
import sys 
sys.path.append('/Users/ctang/Code/My_Python_Code/')
import ctang

#=================================================== Definitions
N_region = 9
N_model = 21
VAR = ('rsds','tas','sfcWind') #,'PVpot')
N_var = 6 # var to be ploted
#---------------------------------------------------  data
# 21 * 4 table: 21 models vs 4 vars

    #               region_1 region_2 ... region_N_region
    #  model_1
    #  model_2
    #  ...
    #  model_N_model
    #  ens
#--------------------------------------------------- 
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/'
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
filefix=(\
        '.hist.day.1970-1999.SA.timmean.nc' ,\
        '.rcp85.day.2070-2099.SA.timmean.nc')
#each map is mean_ref[0,0,:,:].shape

# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf(\
        Data+VAR[0]+'_AFR-44_'+RCM_Model[1]+filefix[0])

# Read Ensmean of timmean for CORDEX 

mean_ref = np.zeros((N_var,N_model,lons.shape[0],lats.shape[1]))
mean_future = np.zeros((N_var,N_model,lons.shape[0],lats.shape[1]))

for v in range(len(VAR)):
    mean_ref[v*2]=np.array([ctang.read_lonlatmap_netcdf(VAR[v],\
        Data+VAR[v]+'_AFR-44_'+RCM_Model[i]+filefix[0])\
        for i in range(N_model)])
    mean_future[v*2]=np.array([ctang.read_lonlatmap_netcdf(VAR[v],\
            Data+VAR[v]+'_AFR-44_'+RCM_Model[i]+filefix[1]) \
            for i in range(N_model)])
#=================================================== end of read
#=================================================== Cal
mean_change=mean_future-mean_ref
# make vws in %
mean_change[4]=(mean_change[4]*100)/mean_ref[4]

# calculate PVpot changes and changes improsed by tas and sfcwind
# estimation follow Martin Wild 2015
Ref_PVpot=ctang.PVpot(mean_ref[0],mean_ref[1])
Future_PVpot=ctang.PVpot(mean_future[0],mean_future[1])

# PVpot changes in % by jerez 2015
#Ref_PVpot=ctang.PVpot2(mean_ref[0],mean_ref[1],mean_ref[2])
#Future_PVpot=ctang.PVpot2(mean_future[0],mean_future[1],mean_future[2])

mean_change[1]= (Future_PVpot-Ref_PVpot)*100/Ref_PVpot 

# TAS introduce PVpot changes in %
#Delta_PVpot introduced by TAS/VWS is a3/a4 * rsds * Delta_TAS/VWS
#mean_change[3]=-4.715e-6*mean_ref[0]*mean_change[2]*100/Ref_PVpot # Jerez 2015
mean_change[3]=-0.0045*mean_ref[0]*mean_change[2]*100/Ref_PVpot # Martin 2015
#mean_change[3]=-4.715e-6*mean_ref[0]*mean_change[2]*100/Ref_PVpot # Jerez 2015
mean_change[5]= 7.64e-6*mean_ref[0]*mean_ref[4]*100/Ref_PVpot # Jerez 2015
print mean_change[3]
Emean_change=np.mean(mean_change,axis=1)
print Emean_change.shape
#--------------------------------------------------- 
# for t test in each point
T=2.086 # dof=20
t_value=mean_change
#--------------------------------------------------- 
#t_value in 4D N_var,N_model,lat,lon
def significant_map(t_value):
    for v in range(N_var):
        for m in range(N_model):
            print 'calculating t test in dim: ',str(v),str(m)
            for lat in range(t_value.shape[2]):
                for lon in range(t_value.shape[3]):
                    if abs(stats.ttest_1samp(\
                            t_value[v,:,lat,lon],t_value[v,m,lat,lon])[0]) > T:
                        t_value[v,m,lat,lon]=t_value[v,m,lat,lon]
                    else:
                        t_value[v,m,lat,lon]=np.NaN
    return t_value
#--------------------------------------------------- 

def robust(t_value):
# The following criteria are used to
# consider them as either uncertain or negligible:
#   negligible: NaN >= N_model/2
#   uncertain:  at least 2 significant individual differ in sign
#   robust:     NaN < N_model/2, sum(abs(x_i)) = abs(sum(x_i))

    count=1
    for v in range(N_var):
        for lat in range(t_value.shape[2]):
            for lon in range(t_value.shape[3]):
                if np.isnan(t_value[v,:,lat,lon]).sum() >= N_model/2:     # negligible
                    t_value[v,:,lat,lon]=[ 99999 for t in range(N_model)] 
                else:
                    if (N_model-np.isnan(t_value[v,:,lat,lon]).sum()) > np.abs(np.nansum(np.sign(t_value[v,:,lat,lon]))):
                        t_value[v,:,lat,lon]=[-99999 for t in range(N_model)] # uncertain
                    else:
                        t_value[v,:,lat,lon]=t_value[v,:,lat,lon]
                        #print "get robust point",count
                        count += 1
    return t_value

#t_value=robust(significant_map(t_value))
#ctang.Save2mat('ttest',t_value) # save to txt file to save time

#t_value=ctang.Loadmat('ttest')

# if not significance
t_value=mean_change

print t_value

Ensmean_change_ttest=np.nanmean(t_value,axis=1)
print Ensmean_change_ttest.shape
#Ensmean_change_ttest=np.ma.masked_invalid(np.nanmean(t_value,axis=1))

#=================================================== plot
degree_sign= u'\N{DEGREE SIGN}'
Title='rsds change %'
Unit=(('(W/m2)','(%)'),\
    ('('+degree_sign+'C)','(%)'),\
    ('(%)','(%)'))

TTT=(('Mean RSDS changes RCP8.5 ', 'PVpot changes RCP8.5 '),\
        ('Mean TAS changes RCP8.5 ', 'TAS-induced PVpot changes RCP8.5 '),\
        ('Mean VWS changes RCP8.5 ', 'VWS-induced PVpot changes RCP8.5 '))



#--------------------------------------------------- 
def PlotMap(array2D,lons,lats,m,k,axx,vmin,vmax):
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(vmin,vmax,9)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    map=Basemap(projection='cyl',llcrnrlat=lats[:,0].min(),urcrnrlat=lats[:,0].max(),\
            llcrnrlon=lons[0,:].min(),urcrnrlon=lons[0,:].max(),resolution='l')
    ctang.setMap(map)

    x,y=map(lats,lons)

    map.pcolormesh(y,x,array2D,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
    #axx.axis('off')
    axx.xaxis.set_visible(False)
    axx.yaxis.set_visible(False)

    plt.title(TTT[m][k]+Unit[m][k],fontsize= 8)
    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=1) 
    cb.ax.tick_params(labelsize=5) 

Limit=np.array([\
        [[-10,10],[-5,5]],\
        [[2,6],[-2.5,2.5]],\
        [[-10,10],[-2,2]]])
#--------------------------------------------------- 
fig, axes = plt.subplots(nrows=3, ncols=2,\
        figsize=(6, 8),facecolor='w', edgecolor='k') # figsize=(w,h)
#fig.subplots_adjust(hspace=0.35,top=0.96,wspace=0.3)

print Ensmean_change_ttest[1].min()

for i in range(len(VAR)):
    for j in range(2):
        plt.sca(axes[i,j]) # active shis subplot for GCM
        ax=axes[i,j]
        m=i*2+j
        print m
        PlotMap(Ensmean_change_ttest[m],lons,lats,i,j,ax,Limit[i,j][0],Limit[i,j][1])


#plt.savefig('map.max.min.eps',format='eps')
plt.savefig('map.max.min.png')
plt.show()

quit()
