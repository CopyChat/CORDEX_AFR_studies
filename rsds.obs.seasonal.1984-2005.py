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
OBS_Dir='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/rsds_validation/'
VAR ='rsds' # ,'tas','sfcWind') #,'PVpot')

degree_sign= u'\N{DEGREE SIGN}'
### if plot only the obs not bias:
BIAS=1
#=================================================== test
##
#=================================================== end of test
# use CanESM2 instead of all not avail model:

Season=('DJF','JJA')

OBS_name=(\
	'SARAH-2',\
        # 'CM_SAF_CDR',\
	'SRB',\
	'ERA_Interim',\
        # 'NCAR',\
        'CFSR',\
        )

### if want to plot less OBS, please comment some OBS out.

Resolution=(\
	str('0.5*0.5')+degree_sign,\
        # str('0.3*0.3')+degree_sign,\
	str('1.0*1.0')+degree_sign,\
	str('0.75*0.75')+degree_sign,\
	# str('1.875*1.875')+degree_sign,\
	str('0.3*0.3')+degree_sign,\
        )
VAR=(\
    'SIS',\
    # 'SIS',\
    'sw_sfc_dn',\
    'ssrd',\
    # 'dswrf',\
    'DSWRF_P8_L1_GGA0',\
    )

Unit=(\
    '(W/$m^{2}$)',\
    # '(W/$m^{2}$)',\
    '(W/$m^{2}$)',\
    '(W/$m^{2}$)',\
    # '(W/$m^{2}$)',\
    '(W/$m^{2}$)',\
    )


Temp_cover=(\
    '1984-2005',\
    # '1984-2005',\
    '1984-2005',\
    '1984-2005',\
    # '1984-2005',\
    '1984-2005',\
    )
N_obs=len(OBS_name)

#=================================================== function of plot

def PlotMap(array2D,lons,lats,m,k,axx,vmin,vmax,cmap):
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
        figsize=(6.5*len(OBS_name), 6),facecolor='w', edgecolor='k') # figsize=(w,h)
fig.subplots_adjust(hspace=0.2,top=0.90,wspace=0.2)
#=================================================== 
#=================================================== reading data
# 21 * 4 table: 21 models vs 4 vars
for S in range(2):

    Reference='SISmm.SARAH-E_v2.1984-2005.SA.'+str(Season[S])+'.timmean.nc'

    OBSfile=(\
        # 'SISmm.CDR.mon.mean.1984-2005.SA.'+str(Season[S])+'.timmean.nc',\
        'srb_rel3.0_shortwave_monthly_utc.1984-2005.SA.'+str(Season[S])+'.timmean.nc',\
        'ERA_In.ssrd.mon.mean.1984-2005.SA.'+str(Season[S])+'.timmean.nc',\
        # 'dswrf.sfc.gauss.monmean.1984-2005.SA.'+str(Season[S])+'.timmean.nc',\
        'cfsr.rsds.mon.1984.2005.SA.'+str(Season[S])+'.timmean.nc',\
        )

    Ref_remap=(\
        # 'SISmm.SARAH-E_v2.1984-2005.SA.'+str(Season[S])+'.timmean.remap.cdr.nc',\
        'SISmm.SARAH-E_v2.1984-2005.SA.'+str(Season[S])+'.timmean.remap.srb.nc',\
        'SISmm.SARAH-E_v2.1984-2005.SA.'+str(Season[S])+'.timmean.remap.era_in.nc',\
        # 'SISmm.SARAH-E_v2.1984-2005.SA.'+str(Season[S])+'.timmean.remap.ncep.nc',\
        'SISmm.SARAH-E_v2.1984-2005.SA.'+str(Season[S])+'.timmean.remap.cfsr.nc',\
        )

#=================================================== calculate bias

# Read Ref as CM_saf
    print OBS_Dir+Reference
    Ref_OBS=np.array(ctang.read_lonlatmap_netcdf(VAR[0],OBS_Dir+Reference))
    print Ref_OBS.shape

    Ref_OBS[Ref_OBS == -999] = np.nan

    lons,lats=ctang.read_lonlat_netcdf_1D(\
        OBS_Dir+Reference)


    k=0 # column number = N_obs
    plt.sca(axes[S,k]) # active shis subplot for RCM
    axx=axes[S,k]

    cmap = plt.cm.jet

    print S,k,"-----------"
    PlotMap(Ref_OBS,lons,lats,S,k,axx,60,360,cmap)

    ### if plot the bias
    if BIAS==1:
        cmap = plt.cm.seismic

    for i in range(N_obs-1):
        k=k+1
        Ref_OBS_remap=np.array(ctang.read_lonlatmap_netcdf(VAR[0],OBS_Dir+Ref_remap[i]))

        Ref_OBS_remap[Ref_OBS_remap == -999] = np.nan

        print i
        print OBS_Dir+OBSfile[i]

        OBS=np.array(ctang.read_lonlatmap_netcdf(VAR[i+1],OBS_Dir+OBSfile[i]))
        print Ref_OBS_remap.shape,OBS.shape

        lonsOBS,latsOBS=ctang.read_lonlat_netcdf_1D(\
            OBS_Dir+OBSfile[i])

        # plot the bias
        print S,k,"-----------"
        plt.sca(axes[S,k]) # active shis subplot for RCM
        axx=axes[S,k]

        ### if plot the bias
        if BIAS==1:
            PlotMap(OBS-Ref_OBS_remap,lonsOBS,latsOBS,S,k,axx,-100,100,cmap)
        else:
            PlotMap(OBS,lonsOBS,latsOBS,S,k,axx,60,360,cmap)



# plt.suptitle(Title)

OutputImage='rsds.obs.seasonal'
#plt.savefig(OutputImage+'.eps',format='eps')
plt.savefig(OutputImage+'.png')
plt.show()
quit()


