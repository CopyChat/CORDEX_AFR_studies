#!/usr/bin/env python
"""
========
Ctang, A map of RSDS over Southern Africa.
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
import ctang

#=================================================== Definitions
OBS_Dir='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/obs.clt/'
VAR ='rsds' # ,'tas','sfcWind') #,'PVpot')
OBS='ERA_Interim'
N_column = 2
N_row = 21
N_plot = N_column*N_row
degree_sign= u'\N{DEGREE SIGN}'
### if plot only the obs not bias:
BIAS=1
#=================================================== test
##
#=================================================== end of test
# use CanESM2 instead of all not avail model:

Season=('DJF','JJA')

OBS_name=(\
	'CM_SAF',\
        'ISCCP',\
	'ERA_Interim',\
        'NCEP-NCAR',\
        )
Resolution=(\
	str('0.25*0.25')+degree_sign,\
	str('2.5*2.5')+degree_sign,\
	str('0.75*0.75')+degree_sign,\
        str('1.875*1.875')+degree_sign,\
        )
VAR=(\
    'cfc',\
    'cltisccp',\
    'tcc',\
    'tcdc',\
    )

Unit=(\
    '(%)',\
    '(%)',\
    '(%)',\
    '(%)',\
    )


Temp_cover=(\
    '1983-2005',\
    '1984-2005',\
    '1979-2005',\
    '1979-2005',\
    )
N_obs=len(OBS_name)

#=================================================== function of plot

def PlotMap(array2D,lons,lats,m,k,axx,vmin,vmax,cmap):
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

    if k==0:
        plt.title(OBS_name[k]+'  '+Unit[k],fontsize=10)
    else:
        ### if plot the bias
        if BIAS==1:
            plt.title(OBS_name[k]+' - '+OBS_name[0]+'  '+Unit[k],fontsize= 10)
        else:
            plt.title(OBS_name[k]+'  '+Unit[k],fontsize=10)


    cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=0.8) 
    cb.ax.tick_params(['{:.0f}'.format(x) for x in bounds ],labelsize=6) 
    #cbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in np.arange(cbar_min, cbar_max+cbar_step, cbar_step)], fontsize=16, weight='bold')
    
    axx.text(0.9, 0.9,str(Season[m]), ha='center', va='center', transform=axx.transAxes)

#===================================================  to plot

fig, axes = plt.subplots(nrows=2, ncols=N_obs,\
        # sharex=True, sharey=True,\
        figsize=(20, 6),facecolor='w', edgecolor='k') # figsize=(w,h)
fig.subplots_adjust(hspace=0.2,top=0.90,wspace=0.2)
#=================================================== 
#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars
for S in range(2):

    Reference='CFCmm.CM_SAF.CLARA-A2.1983-2005.GL.SA.'+str(Season[S])+'.timmean.nc'

    OBSfile=(\
        'cltisccp.monmean.1984-2005.SA.'+str(Season[S])+'.timmean.nc',\
        'ERA_In.clt.mon.mean.1979-2005.SA.'+str(Season[S])+'.timmean.nc',\
        'tcdc.eatm.gauss.mon.mean.1979-2005.SA.'+str(Season[S])+'.timmean.nc',\
        )

    Ref_remap=(\
        'CFCmm.CM_SAF.CLARA-A2.1983-2005.GL.SA.'+str(Season[S])+'.timmean.remap.isccp.nc',\
        'CFCmm.CM_SAF.CLARA-A2.1983-2005.GL.SA.'+str(Season[S])+'.timmean.remap.era_in.nc',\
        'CFCmm.CM_SAF.CLARA-A2.1983-2005.GL.SA.'+str(Season[S])+'.timmean.remap.ncep.nc',\
        )

#=================================================== calculate bias

# Read Ref as CM_saf
    print OBS_Dir+Reference
    Ref_OBS=np.array(ctang.read_lonlatmap_netcdf(VAR[0],OBS_Dir+Reference))
    print Ref_OBS.shape

    lons,lats=ctang.read_lonlat_netcdf_1D(\
        OBS_Dir+Reference)

    k=0 # column number = N_obs
    plt.sca(axes[S,k]) # active shis subplot for RCM
    axx=axes[S,k]

    cmap = plt.cm.jet

    print S,k,"-----------"
    PlotMap(Ref_OBS,lons,lats,S,k,axx,0,100,cmap)

    ### if plot the bias
    if BIAS==1:
        cmap = plt.cm.seismic

    for i in range(N_obs-1):
        k=k+1
        Ref_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(VAR[0],OBS_Dir+Ref_remap[i]))


        print OBS_Dir+OBSfile[i]
        OBS=np.array(ctang.read_lonlatmap_netcdf(VAR[i+1],OBS_Dir+OBSfile[i]))
        print Ref_OBS_remap.shape,OBS.shape

        if OBS_name[k] == 'ERA_Interim':
            OBS=OBS*100

        # mask the missing value > 1000
        OBS[OBS>1000]=np.nan

        lonsOBS,latsOBS=ctang.read_lonlat_netcdf_1D(\
            OBS_Dir+OBSfile[i])

        # plot the bias
        print S,k,"-----------"
        plt.sca(axes[S,k]) # active shis subplot for RCM
        axx=axes[S,k]

        ### if plot the bias
        if BIAS==1:
            PlotMap(OBS-Ref_OBS_remap,lonsOBS,latsOBS,S,k,axx,-50,50,cmap)
        else:
            PlotMap(OBS,lonsOBS,latsOBS,S,k,axx,0,100,cmap)



# plt.suptitle(Title)

OutputImage='rsds.obs.seasonal'
#plt.savefig(OutputImage+'.eps',format='eps')
plt.savefig(OutputImage+'.png')
plt.show()
quit()


