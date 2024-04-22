#!/usr/bin/env python
# -*- coding: utf-8 -*-
## system libraries
from __future__ import print_function
import traceback
from os.path import exists
import glob
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import *
from scipy import interpolate
import pyproj
import pickle
import os,tempfile,sys,glob,subprocess,multiprocessing,time,random
from pkg_resources import parse_version
import  getopt 
from eccodes import *

## matplotlib related
import matplotlib
matplotlib.use('Agg')
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from matplotlib.colors import from_levels_and_colors, BoundaryNorm
from matplotlib.ticker import MaxNLocator
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings("ignore") #for pesky contour warning messages
import CallesEgna

class gribFile:
	def __init__(self, **kwds):
		self.__dict__.update(kwds)
class mapKey:
	def __init__(self, **kwds):
		self.__dict__.update(kwds)
class ReadKey:
	def __init__(self, **kwds):
		self.__dict__.update(kwds)
class keyList:
	def __init__(self, **kwds):
		self.__dict__.update(kwds)
class modelVariable:
	def __init__(self, **kwds):
		self.__dict__.update(kwds)
class attributes:
        def __init__(self, **kwds):
                self.__dict__.update(kwds)
class Location:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
class xVars:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
class Constants:
	def __init__(self, **kwds):
		self.__dict__.update(kwds)

def Combine(object1,object2,delay): #object 1 has to be from 1h newer cycle
    merged = object1

    merged.validity=merged.validity[:,:-1]
    merged.initial=merged.initial[:,:-1]

    merged.surfVars.pres = np.append(merged.surfVars.pres[:,:-1],object2.surfVars.pres[:,delay:],axis=0) 
    merged.surfVars.prtp = np.append(merged.surfVars.prtp[:,:-1],object2.surfVars.prtp[:,delay:],axis=0) 
    merged.surfVars.tcc = np.append(merged.surfVars.tcc[:,:-1],object2.surfVars.tcc[:,delay:],axis=0)
    merged.surfVars.vis = np.append(merged.surfVars.vis[:,:-1],object2.surfVars.vis[:,delay:],axis=0)
    
    merged.icingIndexVars.icei2 = np.append(merged.icingIndexVars.icei2[:,:-1],object2.icingIndexVars.icei2[:,delay:],axis=0)

    merged.integratedVars.lgt = np.append(merged.integratedVars.lgt[:,:-1],object2.integratedVars.lgt[:,delay:],axis=0)
    merged.integratedVars.cb = np.append(merged.integratedVars.cb[:,:-1],object2.integratedVars.cb[:,delay:],axis=0)
    merged.integratedVars.cat_max = np.append(merged.integratedVars.cat_max[:,:-1],object2.integratedVars.cat_max[:,delay:],axis=0)
    merged.integratedVars.cat_b = np.append(merged.integratedVars.cat_b[:,:-1],object2.integratedVars.cat_b[:,delay:],axis=0)
    merged.integratedVars.cat_maxlev = np.append(merged.integratedVars.cat_maxlev[:,:-1],object2.integratedVars.cat_maxlev[:,delay:],axis=0)
    merged.integratedVars.cat_t = np.append(merged.integratedVars.cat_t[:,:-1],object2.integratedVars.cat_t[:,delay:],axis=0)
    merged.integratedVars.lmxice = np.append(merged.integratedVars.lmxice[:,:-1],object2.integratedVars.lmxice[:,delay:],axis=0)
    
    merged.surfAccumVars.rain = np.append(merged.surfAccumVars.rain[:,:-1],object2.surfAccumVars.rain[:,delay:],axis=0)
    merged.surfAccumVars.snow = np.append(merged.surfAccumVars.snow[:,:-1],object2.surfAccumVars.snow[:,delay:],axis=0)
    merged.surfAccumVars.grpl = np.append(merged.surfAccumVars.grpl[:,:-1],object2.surfAccumVars.grpl[:,delay:],axis=0)

    merged.clsVars.t = np.append(merged.clsVars.t[:,:-1],object2.clsVars.t[:,delay:],axis=0)
    merged.clsVars.q = np.append(merged.clsVars.q[:,:-1],object2.clsVars.q[:,delay:],axis=0)
    merged.clsVars.u = np.append(merged.clsVars.u[:,:-1],object2.clsVars.u[:,delay:],axis=0)
    merged.clsVars.v = np.append(merged.clsVars.v[:,:-1],object2.clsVars.v[:,delay:],axis=0)
    merged.clsVars.ugst = np.append(merged.clsVars.ugst[:,:-1],object2.clsVars.ugst[:,delay:],axis=0)
    merged.clsVars.vgst = np.append(merged.clsVars.vgst[:,:-1],object2.clsVars.vgst[:,delay:],axis=0)

    return merged

def CountMatches(number,array):
    return (array==number).sum()

def CountClassValues(base,top,array):
    summa = ((array >= base) & (array < top)).sum()
    return ((array >= base) & (array < top)).sum()

def Facecolor(bp,color):
    for patch in bp['boxes']:
        patch.set_facecolor(color)
    return

def EPSGrams(location):

    subplts=5
    fig, axarr = plt.subplots(subplts,sharex=True, figsize=(25*0.3, 29.5*0.3))#,constrained_layout=True)
    time_steps = matplotlib.dates.date2num(location.validity[0])
    width = 0.03
    medianprops = dict(linewidth=0,linestyle=None)
    linewidth = 1.5

    ########## Total cloud cover #############
    tcc = location.surfVars.tcc
    bp_tcc=axarr[0].boxplot(tcc*100,positions=time_steps,widths=width,showfliers=False,patch_artist=True,medianprops=medianprops,zorder=2)
    Facecolor(bp_tcc,'orange')
    axarr[0].plot(time_steps,tcc[0]*100,color='blue',linewidth=linewidth,linestyle='--',zorder=5)
    axarr[0].set_ylim(0,100)
    axarr[0].invert_yaxis()
    axarr[0].grid(axis='y',linewidth=0.5,linestyle=':')
    axarr[0].set_ylabel('Total cloudiness (' + r'$\mathrm{\%}$' +')',size=9)
    axarr[0].tick_params(labelleft=True,labelright=True,left=True,right=True)

    ######### Precipitation ##############
    rain_increment = location.surfAccumVars.rain[:,1:]-location.surfAccumVars.rain[:,0:-1]
    snow_increment = location.surfAccumVars.snow[:,1:]-location.surfAccumVars.snow[:,0:-1]
    grpl_increment = location.surfAccumVars.grpl[:,1:]-location.surfAccumVars.grpl[:,0:-1]
    precipitation = rain_increment+snow_increment+grpl_increment
    bp_preci=axarr[1].boxplot(precipitation,positions=time_steps[:-1],widths=width,showfliers=False,patch_artist=True,medianprops=medianprops,zorder=2)
    Facecolor(bp_preci,'lawngreen')
    axarr[1].plot(time_steps[:-1],precipitation[0],color='blue',linewidth=linewidth,linestyle='--',zorder=5)
    axarr[1].grid(axis='y',linewidth=0.5,linestyle=':')
    axarr[1].set_ylabel('Precipitation rate (' + r'$\mathrm{mm/h}$' +')',size=9)
    axarr[1].tick_params(labelleft=True,labelright=True,left=True,right=True)
    axarr[1].set_ylim(0,None)
    axarr[1].invert_yaxis()
    
    ######## Precipitation type/probability ############
    ptype = location.surfVars.prtp
    lgt = location.integratedVars.lgt 

    lgt_list = np.array([])
    for i in range(len(ptype[0,:])):
        lgt20 = sum(x >= 20 for x in lgt[:,i]) #probability of thunder is defined when lgt>=20
        lgt_list = np.append(lgt_list,lgt20/len(lgt[:,i]))
    axarr[2].plot(time_steps,lgt_list*100,color='darkred',linewidth=linewidth+0.3,linestyle=':',zorder=5,label="thunder")

    ptype_codes = [('drizzle',0,'yellow'),('rain',1,'green'),('sleet',2,'orange'),('snow',3,'deepskyblue'),('frz drizzle',4,'turquoise'),('frz rain',5,'blue'),('graupel',6,'red'),('hail',7,'hotpink')]
    bottom_bars=np.zeros(len(ptype[0,:]))
    for code in ptype_codes:
        bar_list = np.array([])
        for i in range(len(ptype[0,:])):
            bar = CountMatches(code[1],ptype[:,i])
            bar_list=np.append(bar_list,bar)
        axarr[2].bar(time_steps,(bar_list/len(ptype)*100),width=width,color=code[2],bottom=(bottom_bars/len(ptype)*100),label=code[0],edgecolor='k',linewidth=0.5,zorder=2)
        bottom_bars = bottom_bars+bar_list

    axarr[2].legend(loc='center left', bbox_to_anchor=(1,0.5),prop={'size': 7}, labelspacing=0.3)
    
    axarr[2].set_ylim(0,100)
    axarr[2].invert_yaxis()
    axarr[2].grid(axis='y',linewidth=0.5,linestyle=':')
    axarr[2].set_axisbelow(True)
    axarr[2].set_ylabel('Probability (' + r'$\mathrm{\%}$' +')',size=9)


    #########10m Wind and 10m Gust #################
    WS10m=np.sqrt(location.clsVars.u**2+location.clsVars.v**2)
    GUST10m=np.sqrt(location.clsVars.ugst**2+location.clsVars.vgst**2)
    bp_wind=axarr[3].boxplot(WS10m,positions=time_steps,widths=width-0.015,showfliers=False,patch_artist=True,medianprops=medianprops,zorder=7)
    bp_gust=axarr[3].boxplot(GUST10m,positions=time_steps,widths=width,showfliers=False,patch_artist=True,medianprops=medianprops,zorder=4)
    Facecolor(bp_wind,'bisque')
    Facecolor(bp_gust,'lightgrey')
    axarr[3].plot(time_steps,WS10m[0],color='blue',linewidth=linewidth,linestyle='--',zorder=10)
    axarr[3].plot(time_steps,GUST10m[0],color='blue',linewidth=linewidth,linestyle='--',zorder=10)
    axarr[3].grid(axis='y',linewidth=0.5,linestyle=':')
    axarr[3].set_yticks(np.arange(0,np.amax(GUST10m)+1,2.0))
    axarr[3].set_ylabel('GUST10m, WS10m (' + r'$\mathrm{m/s}$' +')',size=9)
    axarr[3].tick_params(labelleft=True,labelright=True,left=True,right=True)
    axarr[3].invert_yaxis()


    ############# temperature and dew point ################
    q2e=location.surfVars.pres*location.clsVars.q/(0.622+location.clsVars.q*(1.-0.622))
    Td2m=243.5*np.log(q2e/(6.112*100.))/(17.67-np.log(q2e/(6.112*100.)))
    T2m = location.clsVars.t-constants.tZero
    bp_t2m=axarr[4].boxplot(T2m,positions=time_steps,widths=width-0.015,showfliers=False,patch_artist=True,medianprops=medianprops,zorder=7)
    bp_td2m=axarr[4].boxplot(Td2m,positions=time_steps,widths=width,showfliers=False,patch_artist=True,medianprops=medianprops,zorder=4)
    Facecolor(bp_t2m,'salmon')
    Facecolor(bp_td2m,'lightgreen')
    axarr[4].plot(time_steps,T2m[0],color='blue',linewidth=linewidth,linestyle='--',zorder=10)
    axarr[4].plot(time_steps,Td2m[0],color='blue',linewidth=linewidth,linestyle='--',zorder=10)
    axarr[4].grid(axis='y',linewidth=0.5,linestyle=':')
    #axarr[4].set_yticks(np.arange(int(np.amin(Td2m))-1,int(np.amax(T2m))+1,2.0))
    axarr[4].set_ylabel('T2m, Td2m (' + r'$\mathrm{\degree C}$' +')',size=9)
    axarr[4].tick_params(labelleft=True,labelright=True,left=True,right=True)
    axarr[4].invert_yaxis()
    axarr[4].axhline(0, color='black',linewidth=0.6)
    
    ##### END SUBPLOTS
    for subplt in range(subplts):
        axarr[subplt].xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
        axarr[subplt].xaxis.set_major_locator(mdates.DayLocator(interval=1))
        axarr[subplt].xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
        axarr[subplt].xaxis.set_minor_locator(mdates.HourLocator(interval=6))
        axarr[subplt].tick_params(axis='x', which='major',pad=12,labelsize=7)
        axarr[subplt].tick_params(axis='x', which='minor',labelsize=7)
        axarr[subplt].xaxis.set_ticks_position('both')
        axarr[subplt].set_xlim([location.initial[0][0], location.validity[0][-1]])
        axarr[subplt].grid(axis='x',which='minor',linewidth=0.4)
        axarr[subplt].grid(axis='x',which='major',linewidth=0.4,color='k')
        axarr[subplt].invert_yaxis()
    fig.tight_layout()
    plt.subplots_adjust(hspace=0.2,top=0.932)
    fig.suptitle(f'Location: {location.name} (Lon: {location.longitude:.2f}, Lat: {location.latitude:.2f}), Base time: {location.initial[0][0].strftime("%d/%m/%y %H UTC")}, Dashed line: 000Mbr',y=0.995,size=10)
    params = {'legend.fontsize': '7',
              'axes.labelsize': '7',
              'axes.titlesize':'8',
              'xtick.labelsize':'7',
              'ytick.labelsize':'7'}
    rcParams.update(params)
    axarr[0].tick_params(axis='x',which='both',labeltop=True)
    plt.savefig(savePath+'EPSgram_' + location.name)
    plt.close()

###############################################################################################
###############################################################################################
fcCycle=sys.argv[1]
try:
    readData=sys.argv[2]
    print("readData set to ",readData, " in call")
except:
    readData="readData"
    print("readData = ",readData, " by default")
fileKey='/metcoop/transfers/' + fcCycle + '/fc????????' + fcCycle + '+00?grib2*_*'
wrkdir='/data/hirlam2/Python/'
savePath=wrkdir+ fcCycle +'/EPSgrams/'

typesOfLevel=['heightAboveGround']
levels = ['any']
shortNames= ['lsm']
perturbationNumbers=['any']
stepTypes= ['instant']
dataDates='any'
dataTimes='any'
validityDates='any'
validityTimes='any'
mapResolution='i'
makeMap=True
cmap = 'viridis'       
gribFill = 9999
listAll =  False  

constants=Constants(gravity=9.81, tZero=273.15, TfC=-52, TcC=52, TfK=221.15, TcK=291.15, knotcoef=1.94384)


locations=pickle.load(open('/data/hirlam2/Python/'+ fcCycle +'/EPSgrams/EPSG_locations.pickle','rb'))

#getting data from previous cycles for 30 members ensemble.
fc_Hour = datetime.datetime.strptime(fcCycle,"%H")
previous1H_cycle = (fc_Hour-timedelta(hours=1)).strftime("%H")
previous2H_cycle = (fc_Hour-timedelta(hours=2)).strftime("%H")
previous3H_cycle = (fc_Hour-timedelta(hours=3)).strftime("%H")
previous4H_cycle = (fc_Hour-timedelta(hours=4)).strftime("%H")
previous5H_cycle = (fc_Hour-timedelta(hours=5)).strftime("%H")

previous1H = pickle.load(open(wrkdir+ previous1H_cycle +'/EPSgrams/EPSG_locations.pickle','rb'))
previous2H = pickle.load(open(wrkdir+ previous2H_cycle +'/EPSgrams/EPSG_locations.pickle','rb'))
previous3H = pickle.load(open(wrkdir+ previous3H_cycle +'/EPSgrams/EPSG_locations.pickle','rb'))
previous4H = pickle.load(open(wrkdir+ previous4H_cycle +'/EPSgrams/EPSG_locations.pickle','rb'))
previous5H = pickle.load(open(wrkdir+ previous5H_cycle +'/EPSgrams/EPSG_locations.pickle','rb'))

###### combinating cycles and plotting ########
for location,location1H,location2H,location3H,location4H,location5H in zip(locations,previous1H,previous2H,previous3H,previous4H,previous5H):
    combined = location
    combined = Combine(combined,location1H,1)
    combined = Combine(combined,location2H,2)
    combined = Combine(combined,location3H,3)
    combined = Combine(combined,location4H,4)
    combined = Combine(combined,location5H,5)
    try:
        print('Plotting:',combined.name,combined.longitude,combined.latitude,datetime.datetime.now())
        EPSGrams(combined)
    
    except:
        print('Unable to plot:',combined.name,combined.longitude,combined.latitude,', Skipping..',datetime.datetime.now())
