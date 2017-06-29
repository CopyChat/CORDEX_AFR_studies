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
#=================================================== ploting
Title='domain'
#--------------------------------------------------- 
#fig,axes = plt.subplots(111,figsize=(12, 12),facecolor='w', edgecolor='k')
fig, ax = plt.subplots()
#=================================================== 
map=Basemap(projection='cyl',llcrnrlat=-50,urcrnrlat=20,\
        llcrnrlon=0,urcrnrlon=60,resolution='h')
ctang.setMap(map)
plt.suptitle(Title)

plt.savefig('empty.map.eps',format='eps')
plt.savefig('empty.map.png')


plt.show()
quit()
