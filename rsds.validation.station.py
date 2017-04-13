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
import pandas as pd
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
Data='/Users/ctang/Code/CORDEX_AFR_studies/data/validation/'

OBS_Dir='/Users/ctang/Code/CORDEX_AFR_studies/data/OBS/'
obs='rsds.southern.africa.gt.5year.1970-1999.csv'

N_model = 21
VAR ='rsds' 

NYEAR=313

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
# reading GEBA

OBSfile=OBS_Dir+obs
GEBA = np.array(pd.read_csv(OBSfile,index_col=False))

GEBA_Data = GEBA[:,7:19]
STATION=GEBA[:,0]
Jan=GEBA[:,7]


# for cordex data:
CORDEX_file=[Data+VAR+'.AFR-44_'+j+'.csv' for j in RCM_Model]
CORDEX = np.zeros((N_model,NYEAR))
CORDEX=np.array([np.array(pd.read_csv(i,index_col=False)) for i in CORDEX_file])
print CORDEX.shape
print CORDEX[1,1,3:15]


STATION_name ='station_name_32'
station_name = np.array(pd.read_csv(OBS_Dir+STATION_name,header=None))[:,1]
station_ID = np.array(pd.read_csv(OBS_Dir+STATION_name,header=None))[:,0]


#--------------------------------------------------- 
# functino to plot CORDEX vs OBS
def VS(x,y,ax,i,title):
    vmin=np.min(GEBA_Data)
    vmax=np.max(GEBA_Data)
    
    ax.set_xlim(vmin,vmax)
    ax.set_ylim(vmin,vmax)

    ax.set_xticks(range(150,351,50))
    ax.set_yticks(range(150,351,50))

    ax.set_xlabel('GEBA',fontsize=7)
    ax.set_ylabel('RCM',fontsize=7)

    ax.tick_params(direction='in',length=2,width=1,labelsize=6)

    ax.set_title(str(i+1)+". "+title[i],fontsize=6)
    
    ax.set_axisbelow(True)

    ax.scatter(x,y,s=1,facecolors='blue',zorder=2)
    ax.set_aspect('equal')

    #print(type(x))
    #print(type(y))
    #meanbias=np.mean(y-x)
    #ax.text( 150,350,str(meanbias),ha='center', rotation=0)   # plot vs line

    # no. of records
    NO=len(list(x))
    ax.text( 380,165,'#:'+str(NO),ha='right', fontsize=8, rotation=0)   

    # ref line:
    k=np.linspace(100,400,301)
    ax.plot(k,k,'k-',zorder=5,color='black') # identity line

    # linear fitting
    m,b = np.polyfit(x, y, 1)
    yy=[t*m+b for t in range(400)]
    ax.plot(range(400),yy,'--',color='red',zorder=10,label='fitting')

    legend = ax.legend(loc='upper left', shadow=False,prop={'size':8})

    # add correof
    cof=np.ma.corrcoef(x,y)[0,1]
    ax.text( 380,135,'cof:'+str(format(cof,'.2f')),ha='right', fontsize=8, rotation=0)   # plot vs line
    return format(cof,'.2f')
#--------------------------------------------------- 


#=================================================== plot by 21 models

def plot_by_model(title):
    COF=np.zeros((N_model,len(station_ID)))

    for i in range(N_model):
    #for i in range(2):
        print("plotting in model",str(i+1))
        fig, axes = plt.subplots(nrows=5, ncols=7,\
            figsize=(14,10),facecolor='w', edgecolor='k') # (w,h)
        fig.subplots_adjust(left=0.03,bottom=0.03,right=0.99,top=0.94,wspace=0.3,hspace=0.34)
        axes = axes.flatten() # reshape plots to 1D if needed

        for j in range(len(station_ID)):
            sta=station_ID[j]

            # prepare modeldata
            Model_array=CORDEX[i]
            Model_sta=np.array(Model_array[np.where(Model_array[:,1]==sta)])
            Model_sta=Model_sta[:,3:15].flatten()
        
            # prepare obs
            GEBA_PlotData=np.array(GEBA[np.where(GEBA[:,0]==sta)])
            GEBA_PlotData=GEBA_PlotData[:,7:19].flatten()

            # check
            print("-------",sta,j,Model_sta.shape,GEBA_PlotData.shape)

            # to plot
            COF[i,j]=VS(np.array(np.float32(GEBA_PlotData)),np.array(np.float32(Model_sta)),\
                axes[j],j,title)

        plt.suptitle('CORDEX simulated monthly RSDS (W/m2) vs GEBA in 32 stations '+\
                '('+RCM_Model[i]+' NO. '+str(i+1)+")",fontsize=14)
        outfile='validation.rsds.model_NO.'+str(i+1)
        plt.savefig(outfile+'.png')
        plt.savefig(outfile+'.eps', format='eps')

#=================================================== save cof
    headers=['Sta_'+str(i+1) for i in range(len(station_ID))]
    with open('CORDEX.validation.GEBA.1970-1999.cof.csv', 'w') as fp:
        fp.write(','.join(headers) + '\n')
        np.savetxt(fp, COF, '%5.2f', ',')
#=================================================== end plot by model
#=================================================== plot by stations 32 plots

def plot_by_station(title):
    COF=np.zeros((len(station_ID),N_model))
    #for i in range(len(station_ID)):
    for i in range(2):
        print("plotting in station",str(i+1))
        fig, axes = plt.subplots(nrows=4, ncols=6,\
            figsize=(14,9),facecolor='w', edgecolor='k') # (w,h)
        fig.subplots_adjust(left=0.03,bottom=0.04,right=0.98,hspace=0.32,top=0.92,wspace=0.35)
        axes = axes.flatten() # reshape plots to 1D if needed

        sta=station_ID[i]
        for j in range(N_model):
            # prepare modeldata
            Model_array=CORDEX[j]
            Model_sta=np.array(Model_array[np.where(Model_array[:,1]==sta)])
            Model_sta=Model_sta[:,3:15].flatten()
        
            # prepare obs
            GEBA_PlotData=np.array(GEBA[np.where(GEBA[:,0]==sta)])
            GEBA_PlotData=GEBA_PlotData[:,7:19].flatten()

            # check
            print("-------",sta,j,Model_sta.shape,GEBA_PlotData.shape)

            # to plot
            COF[i,j]=VS(np.array(np.float32(GEBA_PlotData)),np.array(np.float32(Model_sta)),\
                axes[j],j,title)

        plt.suptitle('CORDEX simulated monthly RSDS (W/m2) vs GEBA in '+str(station_name[i])+' (NO. '+str(i+1)+')' ,fontsize=16)
        outfile='validation.rsds.station_NO.'+str(i+1)
        plt.savefig(outfile+'.png')
        plt.savefig(outfile+'.eps', format='eps')

#=================================================== save cof
    headers=['RCM_'+str(i+1) for i in range(N_model)]
    with open('CORDEX.validation.GEBA.1970-1999.cof.csv', 'w') as fp:
        fp.write(','.join(headers) + '\n')
        np.savetxt(fp, COF, '%5.2f', ',')
#=================================================== end plot by stations 32 plots

plot_by_station(RCM_Model)
#plot_by_model(station_name)

#=================================================== end
plt.show()
quit()

