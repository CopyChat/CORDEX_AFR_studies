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
import numpy as np
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.basemap import Basemap 

# to load my functions
import sys 
sys.path.append('/Users/ctang/Code/My_Python_Code/')
import ctang

#=================================================== Definitions
OBS='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/'
#=================================================== ploting
Title='RSDS annual mean of 1983-2005 from CM_SAF SARAH-2'
#--------------------------------------------------- 
#fig,axes = plt.subplots(111,figsize=(12, 12),facecolor='w', edgecolor='k')
fig, ax = plt.subplots(figsize=(8,5))
#=================================================== read
ref=np.array(ctang.read_lonlatmap_netcdf('SIS',\
        OBS+'SISmm.SARAH-2.1983-2015.yearmean.SA.timmean.nc'))

# Read lon,lat for model
lons,lats=ctang.read_lonlat_netcdf_1D(\
        OBS+'SISmm.SARAH-2.1983-2015.yearmean.SA.timmean.nc')
print ref.shape
print lons.shape

#=================================================== 
vmin=180
vmax=280

cmap = plt.cm.OrRd
cmaplist = [cmap(i) for i in range(cmap.N)]
bounds = np.linspace(vmin,vmax,11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

map=Basemap(projection='cyl',llcrnrlat=lats[:,0].min(),urcrnrlat=lats[:,0].max(),\
    llcrnrlon=lons[0,:].min(),urcrnrlon=lons[0,:].max(),resolution='h')
ctang.setMap(map)
x,y=map(lats,lons)

map.pcolormesh(y,x,ref,cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
ax.grid(False)

patterns = ['-', '+', 'x', 'o', 'O', '.', '*']  # more patterns

for p in [\
        #                ( (x,y) of left bottom, width, height)
        #patches.Rectangle( (10, -10), 20, 10, hatch=patterns[0],fill=False,alpha=0.1,),\
        patches.Rectangle( (10, -10), 20, 10, fill=False,alpha=0.9,lw=2),\
        patches.Rectangle( (30, -10), 10, 10, fill=False,alpha=0.9,lw=2),\
        patches.Rectangle( (13, -15), 27, 5,  fill=False,alpha=0.9,lw=2),\
        patches.Rectangle( (13, -30), 7, 15,  fill=False,alpha=0.9,lw=2),\
        patches.Rectangle( (20, -35), 15, 20, fill=False,alpha=0.9,lw=2),\
        patches.Rectangle( (43, -25), 7, 10,  fill=False,alpha=0.9,lw=2),\
        patches.Rectangle( (55, -21.5), 3, 2, fill=False,alpha=0.9,lw=2), ]:
    ax.add_patch(p)

plt.text(20, -5, 'Reg. 1', fontsize=8, fontweight='bold', ha='center', va='center', wrap=True,color='white')
plt.text(35, -5, 'Reg. 2', fontsize=8, fontweight='bold', ha='center', va='center', wrap=True,color='white')
plt.text(25, -12.5,'Reg. 3',fontsize=8, fontweight='bold', ha='center', va='center', wrap=True,color='white')
plt.text(17,-22.5,'Reg. 4',fontsize=8, fontweight='bold', ha='center', va='center', wrap=True,color='white')
plt.text(27.5,-25,'Reg. 5',fontsize=8, fontweight='bold', ha='center', va='center', wrap=True,color='white')
plt.text(47,-20,'Reg. 6',fontsize=8, fontweight='bold', ha='center', va='center', wrap=True,color='white')
plt.text(56,-17,'Reg. 7',fontsize=8, fontweight='bold', ha='center', va='center', wrap=True,color='white')

cb=plt.colorbar(cmap=plt.cm.jet,orientation='horizontal',shrink=0.6) 
cb.ax.tick_params(labelsize=6) 





plt.suptitle(Title)

#plt.savefig('ref.rsds.map.eps',format='eps')
plt.savefig('ref.rsds.map.png')


plt.show()
quit()
