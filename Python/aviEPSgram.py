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

    merged.icingIndexVars.icei2 = np.append(merged.icingIndexVars.icei2[:,:,:-1],object2.icingIndexVars.icei2[:,:,delay:],axis=0)

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

def CountClassValuesCAT(base,top,array):
    summa = ((array > base) & (array <= top)).sum()
    return ((array > base) & (array <= top)).sum()

def Facecolor(bp,color):
    for patch in bp['boxes']:
        patch.set_facecolor(color)
    return

def AviEPSGrams(location):
    subplts=5
    fig, axarr = plt.subplots(subplts,sharex=True, figsize=(25*0.3, 29.5*0.3))#,constrained_layout=True)
    time_steps = matplotlib.dates.date2num(location.validity[0])
    width = 0.03
    medianprops = dict(linewidth=0,linestyle=None)
    linewidth = 1.5

    ########## Cloud base probability #############
    cb = location.integratedVars.cb
    cb_classes = [(100,"maroon"),(500,"tomato"),(1000,"orange"),(1500,"yellow"),(2000,"lightgreen"),(5000,"darkturquoise"),(100000,"blue")]
    bottom_bars=np.zeros(len(cb[0,:]))
    previous_class = 0
    for cb_class in cb_classes:
        bar_list = np.array([])
        for i in range(len(cb[0,:])):   
            bar = CountClassValues(previous_class,cb_class[0],cb[:,i]*(1/0.3048)) #m to ft
            bar_list=np.append(bar_list,bar)
        previous_class = cb_class[0]
        if cb_class[0] == 100000:
            axarr[0].bar(time_steps,(bar_list/len(cb)*100),width=width,color=cb_class[1],bottom=(bottom_bars/len(cb)*100),label=f'>{5000}ft',edgecolor='k',linewidth=0.5)
        else:
            axarr[0].bar(time_steps,(bar_list/len(cb)*100),width=width,color=cb_class[1],bottom=(bottom_bars/len(cb)*100),label=f"<{cb_class[0]}ft",edgecolor='k',linewidth=0.5)
        bottom_bars = bottom_bars+bar_list
        
    axarr[0].legend(loc='center left', bbox_to_anchor=(1,0.5),prop={'size': 7}, labelspacing=0.3)
    axarr[0].set_ylim(0,100)
    axarr[0].invert_yaxis()
    axarr[0].grid(axis='both',which='both')
    axarr[0].set_axisbelow(True)
    axarr[0].set_title('Cloud base probability (' + r'$\mathrm{\%}$' +')',size=9)

    ######### CAT and Mean CAT hight ##############
    
    cat_max = location.integratedVars.cat_max
    cat_classes = [('severe',12,'red'),('moderate',8,'yellow'),('light',4,'lightgreen')]
    bottom_bars=np.zeros(len(cat_max[0,:]))
    previous_class = 100
    for cat_class in cat_classes:
        bar_list = np.array([])
        for i in range(len(cat_max[0,:])):
            bar = CountClassValuesCAT(cat_class[1],previous_class,cat_max[:,i])
            bar_list=np.append(bar_list,bar)
        previous_class = cat_class[1]
        axarr[1].bar(time_steps,(bar_list/len(cat_max)*100),width=width,bottom=(bottom_bars/len(cat_max)*100),label=f"{cat_class[0]} (index > {cat_class[1]})",edgecolor='k',linewidth=0.5,color=cat_class[2])
        bottom_bars = bottom_bars+bar_list

    axarr[1].legend(loc='upper center',bbox_to_anchor=(0.5,-0.05),prop={'size': 7}, labelspacing=0.1,ncol=3)
    axarr[1].set_ylim(0,100)
    axarr[1].invert_yaxis()
    axarr[1].grid(axis='both',which='both')
    axarr[1].set_axisbelow(True)
    axarr[1].set_title('CAT probability (' + r'$\mathrm{\%}$' +')' + ' and height (' + r'$\mathrm{ft}$' +')',size=9)
    
    axarr1_2=axarr[1].twinx()
    cat_base = location.integratedVars.cat_b
    cat_base_mean = np.nanmean(cat_base,axis=0)
    cat_top = location.integratedVars.cat_t
    cat_top_mean = np.nanmean(cat_top,axis=0)
    cat_max = location.integratedVars.cat_maxlev
    cat_max_mean = np.nanmean(cat_max,axis=0)

    axarr1_2.plot(time_steps,cat_base_mean*(1/0.3048),color='k',linewidth=0.7,linestyle='--',label="Mean level of CAT base and top")
    axarr1_2.plot(time_steps,cat_top_mean*(1/0.3048),color='k',linewidth=0.7,linestyle='--')
    axarr1_2.plot(time_steps,cat_max_mean*(1/0.3048),color='k',linewidth=0.8,label="Mean level of max CAT")
    axarr1_2.set_ylim(0,max(cat_top_mean)*(1/0.3048)+2000)
    axarr1_2.legend(loc='upper right',prop={'size': 7}, labelspacing=0.2)
    axarr1_2.set_ylabel("ft")
    #axarr1_2.set_yticks(np.arange(0,max(cat_top_mean)*(1/0.3048)+2000,5000))
    
    ######## Fog probability ############
    vis = location.surfVars.vis
    vis_classes = [(50,"red"),(250,"orange"),(500,"yellow"),(750,"lightgreen"),(1000,"deepskyblue")]
    bottom_bars=np.zeros(len(cb[0,:]))
    previous_class = 0
    for vis_class in vis_classes:
        bar_list = np.array([])
        for i in range(len(cb[0,:])):   
            bar = CountClassValues(previous_class,vis_class[0],vis[:,i])
            bar_list=np.append(bar_list,bar)
        previous_class = vis_class[0]
        axarr[2].bar(time_steps,(bar_list/len(cb)*100),width=width,color=vis_class[1],bottom=(bottom_bars/len(cb)*100),label=f"visibility < {vis_class[0]}m",edgecolor='k',linewidth=0.5)
        bottom_bars = bottom_bars+bar_list
        
    axarr[2].legend(loc='upper center',bbox_to_anchor=(0.5,-0.05),prop={'size': 7}, labelspacing=0.1,ncol=5)
    axarr[2].set_ylim(0,100)
    axarr[2].invert_yaxis()
    axarr[2].grid(axis='both',which='both')
    axarr[2].set_axisbelow(True)
    axarr[2].set_title('Fog probability (' + r'$\mathrm{\%}$' +')',size=9)


    #########Icing probability #################
    icingLevels=[150,300,750,1500,2250,3000,3750,4500,5250,6000]
    icei2=location.icingIndexVars.icei2

    icing_max = location.integratedVars.lmxice
    icing_max_mean = (np.nanmean(icing_max,axis=0))*(1/0.3048)

    trace_list = np.array([])
    light_list = np.array([])
    moderate_list = np.array([])
    severe_list = np.array([])
    for timestep in range(len(icei2[0,0,:])):
        zero = 0
        trace = 0
        light = 0
        moderate = 0
        severe = 0
        for member in range(len(icei2[:,0,0])):
            if 4 in icei2[member,:,timestep]:
                severe = severe + 1
            elif 3 in icei2[member,:,timestep]:
                moderate = moderate + 1
            elif 2 in icei2[member,:,timestep]:
                light = light + 1
            elif 1 in icei2[member,:,timestep]:
                trace = trace + 1

        trace_list = np.append(trace_list,trace)
        light_list = np.append(light_list,light)
        moderate_list = np.append(moderate_list,moderate)
        severe_list = np.append(severe_list,severe)

    len_icing = len(icei2[:,0,0])

    axarr[3].bar(time_steps,(severe_list/len_icing*100),width=width,color="red",label="severe",edgecolor='k',linewidth=0.5)
    axarr[3].bar(time_steps,(moderate_list/len_icing*100),bottom=(severe_list)/len_icing*100,width=width,color="orange",label="moderate",edgecolor='k',linewidth=0.5)
    axarr[3].bar(time_steps,(light_list/len_icing*100),bottom=(severe_list+moderate_list)/len_icing*100,width=width,color="yellow",label="light",edgecolor='k',linewidth=0.5)
    axarr[3].bar(time_steps,(trace_list/len_icing*100),bottom=(severe_list+moderate_list+light_list)/len_icing*100,width=width,color="deepskyblue",label="trace",edgecolor='k',linewidth=0.5)

    axarr[3].legend(loc='upper center',bbox_to_anchor=(0.5,-0.05),prop={'size': 7}, labelspacing=0.1,ncol=5)
    axarr[3].set_ylim(0,100)
    axarr[3].invert_yaxis()
    axarr[3].grid(axis='both',which='both')
    axarr[3].set_axisbelow(True)
    axarr[3].set_title('Icing probability (' + r'$\mathrm{\%}$' +')',size=9)

    axarr3_2=axarr[3].twinx()
    axarr3_2.plot(time_steps,icing_max_mean,color='k',linewidth=1.0,label="Mean level of max icing")
    axarr3_2.set_yticks(np.arange(0,21000,3000))
    axarr3_2.legend(loc='upper right',prop={'size': 7}, labelspacing=0.3)
    axarr3_2.set_ylabel("ft")
    axarr3_2.set_axisbelow(True)


    ############# precipitation type and precipitation/thunder probability ################
    ptype = location.surfVars.prtp
    lgt = location.integratedVars.lgt 

    lgt_list = np.array([])
    for i in range(len(ptype[0,:])):
        lgt20 = sum(x >= 20 for x in lgt[:,i]) #probability of thunder is defined when lgt>=20
        lgt_list = np.append(lgt_list,lgt20/len(lgt[:,i]))
    axarr[4].plot(time_steps,lgt_list*100,color='darkred',linewidth=linewidth+0.3,linestyle=':',zorder=5,label="thunder")

    ptype_codes = [('drizzle',0,'yellow'),('rain',1,'green'),('sleet',2,'orange'),('snow',3,'deepskyblue'),('frz drizzle',4,'turquoise'),('frz rain',5,'blue'),('graupel',6,'red'),('hail',7,'hotpink')]
    bottom_bars=np.zeros(len(ptype[0,:]))
    for code in ptype_codes:
        bar_list = np.array([])
        for i in range(len(ptype[0,:])):
            bar = CountMatches(code[1],ptype[:,i])
            bar_list=np.append(bar_list,bar)
        axarr[4].bar(time_steps,(bar_list/len(ptype)*100),width=width,color=code[2],bottom=(bottom_bars/len(ptype)*100),label=code[0],edgecolor='k',linewidth=0.5,zorder=2)
        bottom_bars = bottom_bars+bar_list

    axarr[4].legend(loc='center left', bbox_to_anchor=(1,0.5),prop={'size': 7}, labelspacing=0.3)
    
    axarr[4].set_ylim(0,100)
    axarr[4].invert_yaxis()
    axarr[4].grid(axis='both',which='both')
    axarr[4].set_axisbelow(True)
    axarr[4].set_title('Precipitation type and probability (' + r'$\mathrm{\%}$' +')',size=9)
    
    ##### END SUBPLOTS
    for subplt in range(subplts):
        axarr[subplt].xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
        axarr[subplt].xaxis.set_major_locator(mdates.DayLocator(interval=1))
        axarr[subplt].xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
        axarr[subplt].xaxis.set_minor_locator(mdates.HourLocator(interval=6))
        axarr[subplt].tick_params(axis='x', which='major',pad=12)
        axarr[subplt].set_xlim([location.initial[0][0], location.validity[0][-1]])
        axarr[subplt].grid(axis='x')
        axarr[subplt].invert_yaxis()
    fig.tight_layout()
    fig.suptitle(f'Location: {location.name} (Lon: {location.longitude:.2f}, Lat: {location.latitude:.2f}), Base time: {location.initial[0][0].strftime("%d/%m/%y %H UTC")}',y=0.999,size=10)
    params = {'legend.fontsize': '7',
              'axes.labelsize': '7',
              'axes.titlesize':'8',
              'xtick.labelsize':'7',
              'ytick.labelsize':'7'}
    rcParams.update(params)
    plt.savefig(savePath+'aviEPSgram_' + location.name)
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
savePath=wrkdir+ fcCycle +'/aviEPSgrams/'

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

#getting data from previous cycles for the 30 members ensemble 
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
        AviEPSGrams(combined)
    
    except:
        print('Unable to plot:',combined.name,combined.longitude,combined.latitude,', Skipping..',datetime.datetime.now())
