#!/usr/bin/env python
# -*- coding: utf-8 -*-
## system libraries
from __future__ import print_function
import traceback
from os.path import exists
import glob
## numeric python
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import *
import pyproj
import pickle
#import epygram
#import unittest2,os,tempfile,sys,glob,subprocess,multiprocessing,time,random
import os,tempfile,sys,glob,subprocess,multiprocessing,time,random
from pkg_resources import parse_version
#from cdo import Cdo,CDOException,CdoTempfileStore
import  getopt 
from eccodes import *

## matplotlib related
import matplotlib as mpl
mpl.use('Agg')
#mpl.use('QT4Agg')
#from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
#import matplotlib
import matplotlib.pyplot as plt
from pylab import *
#import matplotlib.axes as ax
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.dates import DateFormatter
import matplotlib.colors as mcolors #new
from matplotlib.colors import from_levels_and_colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
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
class Constants:
	def __init__(self, **kwds):
		self.__dict__.update(kwds)

def mapCreator(template,mapPath='koe',makeMap=True, CustomMap=False,resolution='i'):
    import cartopy.crs as ccrs
    lons=codes_get_array(template,"longitudes")
    lons=np.where(lons>180., lons-360., lons).reshape(codes_get(template,"Nj"),codes_get(template,"Ni"))
    lats=codes_get_array(template,"latitudes").reshape(codes_get(template,"Nj"),codes_get(template,"Ni"))
    dlats= np.zeros(np.shape(lats))
    elats= np.zeros((np.shape(lats)[0]+1,np.shape(lats)[1]+1))
    dlats[:-1,:] =  lats[1:,:]-lats[0:-1,:]
    dlats[-1,:]=dlats[-2,:]
    elats[0:-1,:-1] = lats-dlats/2.
    elats[-1,:-1]=lats[-1,:]+dlats[-1,:]/2.
    elats[:,-1] = elats[:,-2]
    dlons= np.zeros(np.shape(lons))
    elons= np.zeros((np.shape(lons)[0]+1,np.shape(lons)[1]+1))
    dlons[:,:-1] =  lons[:,1:]-lons[:,0:-1]
    dlons[:,-1]=dlons[:,-2]
    elons[:-1,0:-1,] = lons-dlons/2.
    elons[:-1,-1]=lons[:,-1]+dlons[:,-1]/2.
    elons[-1,:] = elons[-2,:]
    if makeMap:
        if CustomMap:
                m3 = 0#Basemap(
                #projection='merc', lat_0=60, lon_0=25,
                #llcrnrlon= elons.min() ,urcrnrlon= elons.max(),
                #llcrnrlat= elats.min(),urcrnrlat= elats.max(),
        #       llcrnrlon=15,urcrnrlon=35,
        #       llcrnrlat=55,urcrnrlat=75,
        #        resolution=resolution)
        elif codes_get(template,"gridType")=='lambert': 
                m3 =  ccrs.LambertConformal(central_longitude=codes_get(template,"LoVInDegrees"), central_latitude=codes_get(template,"LaDInDegrees"), false_easting=0.0, false_northing=0.0, secant_latitudes=None, 
                        standard_parallels=[codes_get(template,"Latin1InDegrees"), codes_get(template,"Latin2InDegrees")], globe=None, cutoff=-30)
                
        elif codes_get(template,"gridType")=='rotated_ll':
                m3 = 0#Basemap(
                #projection='rotpole', lat_0=codes_get_array(template,"latitudes").mean(),
                #o_lon_p=-(180-codes_get(template,"polon")), o_lat_p=-codes_get(template,"polat"), #Basemap uses the north pole, of course
                #lon_0=codes_get_array(template,"polon"),
                #llcrnrlon=elons[0,0],urcrnrlon=elons[-1,-1],
                #llcrnrlat=elats[0,0],urcrnrlat=elats[-1,-1],
                #llcrnrlon=codes_get_array(template,longitudes[0,0],urcrnrlon=codes_get_array(template,longitudes[-1,-1],
                #llcrnrlat=codes_get_array(template,latitudes[0,0],urcrnrlat=codes_get_array(template,latitudes[-1,-1],
                #resolution=resolution)
        #pickle.dump(m3,open(mapPath,'wb'),-1)
    else:
        a=1
        # If a pickle exists, you can uncomment the next comand to load it.
        #m3=pickle.load(open(mapPath,'rb'))

    return m3,lons,lats, elons,elats, m3.transform_points(ccrs.PlateCarree(),  lons,lats)[:,:,0:2], m3.transform_points(ccrs.PlateCarree(),  elons,elats)[:,:,0:2] 

def mapGrib(gribFile, verbose=False):
     #maps the content of a grib file  later usei
     thisMap=[]
     #print(verbose)
     buf=(subprocess.check_output(['grib_ls','-p','count,offset:i,shortName,typeOfLevel,level,perturbationNumber,stepType,dataDate,dataTime,validityDate,validityTime,indicatorOfParameter,edition',gribFile]).splitlines())
     for line in buf:
         edition='nothing found'
         if verbose:
             print(line)
         sample=line.split()
         try:
             edition=int(sample[-1])
         except:
             print('Skipping: ',sample)
         if edition==1:
             try:
                 thisMap.append(attributes(
                 filepath=gribFile,
                 count=int(sample[0]),
                 offset=float(sample[1]),
                 shortName=bytes.decode(sample[2]),
                 typeOfLevel=bytes.decode(sample[3]),
                 level=float(bytes.decode(sample[4])),
                 perturbationNumber=bytes.decode(sample[5]),
                 stepType=bytes.decode(sample[6]),
                 dataDate=int(sample[7]), dataTime=int(sample[8]), validityDate=int(sample[9]), validityTime=int(sample[10]),
                 indicatorOfParameter=int(sample[11])))
             except:
                 print('Skipping: ',  sample)
         elif edition==2:
             try:
                 thisMap.append(attributes(
                 filepath=gribFile,
                 count=int(sample[0]),
                 offset=float(sample[1]),
                 shortName=bytes.decode(sample[2]),
                 typeOfLevel=bytes.decode(sample[3]),
                 level=float(bytes.decode(sample[4])),
                 perturbationNumber=bytes.decode(sample[5]),
                 stepType=bytes.decode(sample[6]),
                 dataDate=int(sample[7]), dataTime=int(sample[8]), validityDate=int(sample[9]), validityTime=int(sample[10]),
                 indicatorOfParameter=bytes.decode(sample[11])))
             except:
                 print('Skipping: ',  sample)
         #else:
         #    print ('Not a valid etition number: ',sample)
     return(thisMap)

def getPV(shortName,typeOfLevel,level,perturbationNumber,stepType):
    '''Manual selector for plotting'''
    #print("Getting field:",shortName,typeOfLevel,"level:",level,stepType,"perturbationNumber:",perturbationNumber,datetime.datetime.now())
    if type(level) == range or type(level) == list or type(level) == ndarray:  # checks if it's a list/range/array selection of levels
        field=grib.seek(int((CallesEgna.findMessage(thisFile.gribMap, shortNames=[shortName], typesOfLevel=[typeOfLevel], levels=level,  perturbationNumbers=[perturbationNumber], stepTypes=[stepType]))[0].offset))
    else:   # assumes it's single level
        field=grib.seek(int((CallesEgna.findMessage(thisFile.gribMap, shortNames=[shortName], typesOfLevel=[typeOfLevel], levels=[float(level)],  perturbationNumbers=[perturbationNumber], stepTypes=[stepType]))[0].offset))
    gid = codes_grib_new_from_file(grib)
    
    pv=np.ma.masked_equal(codes_get_array(gid,"values"),gribFill).reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni"))
    return pv,gid

def getPV2 (message,release=True):
    ''''''
    global gid,m3,lons,lats,elons,elats,xy,exy # will make variables accessible outside of this function
    print ("done globals for mapCreator", datetime.datetime.now())
    ''''''
    print('Doing ',message.shortName)
    grib.seek(int(message.offset))
    print ("done seek", datetime.datetime.now())
    gid = codes_grib_new_from_file(grib)
    print ("done gid", datetime.datetime.now())
    try:
        m3
        print("m3 found!")
    except:
        m3,lons,lats,elons,elats,xy,exy=mapCreator(gid,selected.mapPath,resolution=selected.resolution,makeMap=True)
        print ("done setup m3", datetime.datetime.now())
    pv=np.ma.masked_equal(codes_get_array(gid,"values"),gribFill).reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni"))
    print ("done defining pv", datetime.datetime.now())
    if release:
        codes_release(gid)
    return pv

def mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,vartype,areatype='full',colormap='viridis',SaveFig=True,picname=None):
    print('Plotting:',vartype)
    plt.figure(figsize=figureSize)
    #print ("Done plt.figure(figsize=figureSize)", datetime.datetime.now())
    ax = plt.axes(projection=m3)
    #print ("Done plt.axes",datetime.datetime.now())
    ax.set_extent([exy[0,0,0],exy[-1,-1,0],exy[0,0,1],exy[-1,-1,1]],crs=m3)
    #print ("Done set_extend",datetime.datetime.now())
    coastColor='grey'

    #========================== vartype decisions
    if 'temperature' in vartype:
        name=codes_get(gid,"name")
        units=u"\u2103"
        pvmin=round(pv.min()-constants.Tmelt,2)
        pvmax=round(pv.max()-constants.Tmelt,2)
        pv=pv-constants.Tmelt
        if vartype == 'temperature_zoom':
            boundst=np.linspace(-24,24,25)
            bounds=np.linspace(-24,24,49)
            pcontour=False
            pbarbs=False
            ptext=True
            clayers=None
            ticks=boundst.astype(int)
            cmap=plt.get_cmap(colormap)
            norm=mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            frac=0.15 #0.050
            pad=0.05 #0.001
            labelsize=8
            clim_min=bounds.min()
            clim_max=bounds.max()
            step=None
            barblength=None
            linewidth=None
            weight=None
            texts=[[0.5,-0.05,"min: {} \u2103 max: {} \u2103".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
            cbtitle=''
            cbcolor=None
            cbweight=None
            cbrotation=0
            cbva='bottom'
            cbha='center'
            cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
            layers=[
                    [pv,cmap,8,cbar1],
                    ];
        else:
            boundst=np.linspace(-52,52,27)
            bounds=np.linspace(-52,52,53)
            pcontour=False
            pbarbs=False
            ptext=True
            clayers=None
            ticks=boundst.astype(int)
            cmap=plt.get_cmap(colormap)
            norm=mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            frac=0.15 #0.050
            pad=0.05 #0.001
            labelsize=8
            clim_min=bounds.min()
            clim_max=bounds.max()
            step=None
            barblength=None
            linewidth=None
            weight=None
            texts=[[0.5,-0.05,"min: {} \u2103 max: {} \u2103".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
            cbtitle=''
            cbcolor=None
            cbweight=None
            cbrotation=0
            cbva='bottom'
            cbha='center'
            cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
            layers=[
                    [pv,cmap,8,cbar1],
                    ];
    elif 'wind' in vartype:
        #units=codes_get(gid,"units")
        units=r'${m\ s^{-1}}$'
        if vartype == 'max_winds':
            name=codes_get(gid,"name")
            uReg=pv[0]
            vReg=pv[1]
            wspeed=pv[2]
            bounds=np.array([5,10,20,30,50,70,90,110])*100
            colors=['grey', 'm', 'r', 'orange', 'yellow', 'palegreen', 'lightblue', 'skyblue','dodgerblue']
            boundst=bounds
            pcontour=False
            pbarbs=True
            ptext=True
            clayers=None
            ticks=boundst.astype(int)
            cmap, norm = from_levels_and_colors(bounds,colors,extend='both')
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            frac=0.15
            pad=0.05
            labelsize=8
            pvmin=round(wspeed.min(),1)
            pvmax=round(wspeed.max(),1)
            clim_min=bounds.min()
            clim_max=bounds.max()
            step=40
            barblength=5
            linewidth=0.4
            weight=None
            texts=[[0.5,-0.05,"min: {} m  max: {} m".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
            cbtitle=''
            cbcolor=None
            cbweight=None
            cbrotation=0
            cbva='bottom'
            cbha='center'
            cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
            layers=[
                    [wspeed,cmap,8,cbar1],
                    ];
        if vartype == 'winds_wmax':
            name='W Max'
            pv=np.ma.masked_array(pv,pv<0.25)
            bounds=[0, 0.25, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16]
            boundst=bounds
            pcontour=False
            pbarbs=False
            ptext=True
            clayers=None
            ticks=boundst
            norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            cmap = plt.get_cmap('cool')
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            frac=0.15
            pad=0.05
            labelsize=8
            pvmin=round(pv.min(),1)
            pvmax=round(pv.max(),1)
            clim_min=np.array(bounds).min()
            clim_max=np.array(bounds).max()
            step=None
            barblength=None
            linewidth=None
            weight=None
            texts=[[0.5,-0.05,"min: {} m/s  max: {} m/s".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
            cbtitle=''
            cbcolor=None
            cbweight=None
            cbrotation=0
            cbva='bottom'
            cbha='center'
            cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
            layers=[
                    [pv,cmap,8,cbar1],
                    ];
        elif vartype == 'windspeed':
            name='10m winds'
            uReg=pv[0]
            vReg=pv[1]
            wspeed=pv[2]
            bounds=np.linspace(1,37,37)
            boundst=np.array([1, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 37])
            pcontour=False
            pbarbs=True
            ptext=True
            clayers=None
            ticks=boundst.astype(int)
            norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            cmap = plt.get_cmap(colormap)
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            frac=0.15
            pad=0.05
            labelsize=8
            pvmin=round(wspeed.min(),1)
            pvmax=round(wspeed.max(),1)
            clim_min=bounds.min()
            clim_max=bounds.max()
            step=40
            barblength=5
            linewidth=0.4
            weight=None
            texts=[[0.5,-0.05,"max: {} ".format(pvmax)+units,'bottom','center','horizontal','anchor']]
            cbtitle=''
            cbcolor=None
            cbweight=None
            cbrotation=0
            cbva='bottom'
            cbha='center'
            cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
            layers=[
                    [wspeed,cmap,8,cbar1],
                    ];
    elif 'precipitation' in vartype:
        RAIN=pv[0]
        GRPL=pv[2]
        SNOW=pv[1]
        PMSL=pv[3]
        #units=codes_get(gid,"units")
        units=r'$kg\ m^{-2}$'
        if vartype == 'precipitation_inst':
            units=r'$kg\ m^{-2} h^{-1}$'
            SOLID=(SNOW + GRPL)*3600.
            TOTAL=RAIN*3600. + SOLID
            SOLID=np.ma.masked_array(SOLID,SOLID<0.1)
            TOTAL=np.ma.masked_array(TOTAL,TOTAL<0.1)
            name='Precipitation Intensity'
            bounds=[0.2, 0.5, 1., 2., 4., 8., 16., 32., 64.]
            boundst=bounds
            pcontour=True
            pbarbs=False
            ptext=True
            texts=[[0.5,-0.05, "Max RAIN: {} Max SNOW: {} Max GRAUPEL: {}".format(round(RAIN.max()*3600,2),round(SNOW.max()*3600,2),round(GRPL.max()*3600,2)),'bottom','center','horizontal','anchor']]
            clayers=[PMSL/100]
            colors='red'
            linewidths=0.25
            ticks=bounds
            norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            czorder=11
            crange=range(950,1080,5)
            alpha=1.0
            extend='both'
            labelsize=8
            #pvmin=pv.min()
            #pvmin=pv.max()
            clim_min=np.array(bounds).min()
            clim_max=np.array(bounds).max()
            frac=0.050
            pad=0.001
            cbar1=[extend,boundst,frac,pad,'RAIN','greenyellow','bold',0,'bottom','left',ticks,None]
            cbar2=[extend,boundst,frac,pad,'SOLID','cyan','bold',0,'top','center','',None]
            layers=[
                    [TOTAL,'summer_r',8,cbar1],
                    [SOLID,'cool',9,cbar2],
                    ];
        if vartype == 'precipitation_accum_total':
            SOLID=(SNOW + GRPL)
            TOTAL=RAIN + SOLID
            SOLID=np.ma.masked_array(SOLID,SOLID<0.1)
            TOTAL=np.ma.masked_array(TOTAL,TOTAL<0.1)
            name='Total acc. prec.'
            bounds=np.array([0.5, 1., 2., 4., 8., 16., 32., 64.,128])
            boundst=bounds
            pcontour=False
            pbarbs=False
            ptext=True
            texts=[[0.5,-0.05, "Max RAIN: {} Max SNOW: {} Max GRAUPEL: {}".format(round(RAIN.max(),2),round(SNOW.max(),2),round(GRPL.max(),2)),'bottom','center','horizontal','anchor']]
            clayers=[PMSL/100]
            colors='red'
            linewidths=0.25
            ticks=bounds #np.logspace(-14,0,11)
            norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            czorder=11
            crange=range(950,1080,5)
            alpha=1.0
            extend='both'
            labelsize=8
            #pvmin=pv.min()
            #pvmin=pv.max()
            clim_min=bounds.min()
            clim_max=bounds.max()
            frac=0.15
            pad=0.05
            cbar1=[extend,boundst,frac,pad,'TOTAL','greenyellow','bold',0,'bottom','left',ticks,None]
            layers=[
                    [TOTAL,'summer_r',8,cbar1],
                    ];
        if vartype == 'precipitation_accum_solid':
            SOLID=(SNOW + GRPL)
            TOTAL=RAIN + SOLID
            SOLID=np.ma.masked_array(SOLID,SOLID<0.1)
            TOTAL=np.ma.masked_array(TOTAL,TOTAL<0.1)
            name='Total acc. snow and graupel'
            bounds=np.array([0.5, 1., 2., 4., 8., 16., 32., 64.,128])
            boundst=bounds
            pcontour=False
            pbarbs=False
            ptext=True
            texts=[[0.5,-0.05, "Max RAIN: {} Max SNOW: {} Max GRAUPEL: {}".format(round(RAIN.max(),2),round(SNOW.max(),2),round(GRPL.max(),2)),'bottom','center','horizontal','anchor']]
            clayers=[PMSL/100]
            colors='red'
            linewidths=0.25
            ticks=bounds #np.logspace(-14,0,11)
            norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            czorder=11
            crange=range(950,1080,5)
            alpha=1.0
            extend='both'
            labelsize=8
            #pvmin=pv.min()
            #pvmin=pv.max()
            clim_min=bounds.min()
            clim_max=bounds.max()
            frac=0.15
            pad=0.05
            cbar1=[extend,boundst,frac,pad,'SOLID','cyan','bold',0,'bottom','left',ticks,None]
            layers=[
                    [SOLID,'cool',8,cbar1],
                    ];
    elif 'hail' in vartype:
        name=codes_get(gid,"name")
        units="Risk of Hail"
        pv=np.ma.masked_array(pv,pv<(pv.min()+17.0))
        bounds = np.array([17, 20, 70])
        boundst=bounds
        pcontour=False
        pbarbs=False
        ptext=False
        clayers=None
        #ticks=boundst.astype(int)
        ticks=[]
        norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256, clip=False)
        #cmap = plt.get_cmap(colormap)
        cmap = mcolors.ListedColormap(['orange','red'])
        zorder=8
        czorder=None
        alpha=1.0
        extend='neither'
        frac=0.15
        pad=0.05
        labelsize=8
        #pvmin=int(pv.min())
        #pvmax=int(pv.max())
        clim_min=bounds.min()
        clim_max=bounds.max()
        step=None
        barblength=None
        linewidth=None
        weight=None
        # texts=[[0.5,-0.05,"min: {} max: {}".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
        cbtitle='Orange=Hail, Red=Large Hail'
        cbcolor='k'
        cbweight=None
        cbrotation=0
        cbva='bottom'
        cbha='center'
        cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
        layers=[
               [pv,cmap,8,cbar1],
               ];
    elif 'MLH' == vartype:
        name=codes_get(gid,"name")
        units=codes_get(gid,"units")
        pv=pv
        bounds=np.linspace(0,5000,51)
        boundst=np.linspace(0,5000,11)
        pcontour=False
        pbarbs=False
        ptext=True
        clayers=None
        ticks=np.linspace(0,5000,11).astype(int)
        norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
        cmap = plt.get_cmap(colormap)
        zorder=8
        czorder=None
        alpha=1.0
        extend='max'
        frac=0.15
        pad=0.05
        labelsize=8
        pvmin=int(pv.min())
        pvmax=int(pv.max())
        clim_min=bounds.min()
        clim_max=bounds.max()
        step=None
        barblength=None
        linewidth=None
        weight=None
        texts=[[0.5,-0.05,"min: {} m  max: {} m".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
        cbtitle=''
        cbcolor=None
        cbweight=None
        cbrotation=0
        cbva='bottom'
        cbha='center'
        cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
        layers=[
               [pv,cmap,8,cbar1],
               ];
    elif vartype == 'lightning':
        name='Lightning Intensity' #codes_get(gid,"name")
        units=r'$\rm{strikes \backslash hour\ in\ 30 \times\ 30km^{2}}$'
        bounds = np.array([1, 10, 50, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750])
        boundst=bounds
        pcontour=False
        pbarbs=False
        ptext=True
        clayers=None
        ticks=boundst.astype(int)
        norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
        cmap = plt.get_cmap(colormap)
        zorder=8
        czorder=None
        alpha=1.0
        extend='both'
        frac=0.15
        pad=0.05
        labelsize=8
        pvmin=int(pv.min())
        pvmax=int(pv.max())
        pv=np.ma.masked_array(pv,pv<5)
        clim_min=bounds.min()
        clim_max=bounds.max()
        step=None
        barblength=None
        linewidth=None
        weight=None
        texts=[[0.5,-0.05,"min: {} max: {}".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
        cbtitle=''
        cbcolor=None
        cbweight=None
        cbrotation=0
        cbva='bottom'
        cbha='center'
        cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
        layers=[
               [pv,cmap,8,cbar1],
               ];
    elif vartype == 'HTI':
        HTI=pv[0]
        PMSL=pv[1]
        name=codes_get(gid,"name")
        units=codes_get(gid,"units")
        bounds = np.array([0.73, 0.815, 0.9, 0.9454, 0.999])
        boundst=bounds
        pcontour=True
        pbarbs=False
        ptext=True
        clayers=[PMSL/100]
        ticks=['0.73','low', '0.9', 'moderate','severe \n 0.99']
        tickLabelRotation=[0,-90,0,-90,0]
        colors=['white','burlywood','burlywood','darkorange','darkorange','black']
        cmap, norm = from_levels_and_colors(bounds,colors,extend='both')        
        zorder=8
        clayers=[PMSL/100]
        colors='blue'
        linewidths=0.25
        czorder=11
        crange=range(950,1080,5)
        alpha=1.0
        extend='both'
        frac=0.15
        pad=0.05
        labelsize=8
        pvmin=HTI.min()
        pvmax=HTI.max()
        clim_min=bounds.min()
        clim_max=bounds.max()
        step=None
        barblength=None
        linewidth=None
        weight=None
        texts=[[0.5,-0.05,"min: {} max: {}".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
        cbtitle=''
        cbcolor=None
        cbweight=None
        cbrotation=0
        cbva='bottom'
        cbha='center'
        cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
        layers=[
               [HTI,cmap,8,cbar1],
               ];
    elif 'CAPE' in vartype:
        if vartype == 'MU-CAPE':
            name='MU-CAPE' #codes_get(gid,"name")
            units='$J\ kg{-1}$'
            pv=np.ma.masked_array(pv,pv<50)
            bounds=np.array([2,3,4,5,6,8,10,12,14,16,18,20,22,24,28])*100.
            boundst=bounds
            pcontour=False
            pbarbs=False
            ptext=True
            clayers=None
            ticks=boundst.astype(int)
            norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            cmap = plt.get_cmap(colormap)
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            frac=0.15
            pad=0.05
            labelsize=8
            pvmin=round(pv.min(), 2) #int(pv.min())
            pvmax=round(pv.max(), 2) #int(pv.max())
            clim_min=bounds.min()
            clim_max=bounds.max()
            step=None
            barblength=None
            linewidth=None
            weight=None
            texts=[[0.5,-0.05,"max: {} ".format(pvmax)+units,'bottom','center','horizontal','anchor']]
            cbtitle=''
            cbcolor=None
            cbweight=None
            cbrotation=0
            cbva='bottom'
            cbha='center'
            cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
            layers=[
                    [pv,cmap,8,cbar1],
                    ];
    elif vartype == 'visibility':
        name=codes_get(gid,"name")
        units=codes_get(gid,"units")
        pv=np.ma.masked_array(pv,pv>5000)
        bounds=[50,100,150,350,600,800,1000,5000]
        boundst=bounds
        pcontour=False
        colors=['purple', 'red', 'orange', 'yellow', 'greenyellow', 'limegreen', 'cornflowerblue', 'lightblue', 'white']
        pbarbs=False
        ptext=True
        clayers=None
        ticks=boundst
        cmap,norm = from_levels_and_colors(bounds,colors,extend='both')
        zorder=8
        czorder=None
        alpha=1.0
        extend='both'
        frac=0.15
        pad=0.05
        labelsize=8
        pvmin=int(pv.min())
        pvmax=int(pv.max())
        clim_min=np.array(bounds).min()
        clim_max=np.array(bounds).max()
        step=None
        barblength=None
        linewidth=None
        weight=None
        texts=[[0.5,-0.05,"min: {} m".format(pvmin),'bottom','center','horizontal','anchor']]
        cbtitle=''
        cbcolor=None
        cbweight=None
        cbrotation=0
        cbva='bottom'
        cbha='center'
        cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
        layers=[
               [pv,cmap,8,cbar1],
               ];
    elif vartype == 'ptype':
        name=codes_get(gid,"name")
        units=codes_get(gid,"units")
        pv=pv
        bounds=range(8)
        boundst=bounds
        pcontour=False
        pbarbs=False
        ptext=False
        clayers=None
        ticks=['drizzle','rain','sleet','snow','fr.drizzle','fr.rain','graupel','hail']
        cmap = mcolors.ListedColormap(['greenyellow','green','orange','#069af3','hotpink','purple','red','darkred'])
        norm=None
        zorder=8
        czorder=None
        alpha=1.0
        extend='neither'
        frac=0.15
        pad=0.05
        labelsize=8
        pvmin=pv.min()
        pvmax=pv.max()
        clim_min=-0.5
        clim_max=8-0.5
        step=None
        barblength=None
        linewidth=None
        weight=None
        texts=None
        cbtitle=''
        cbcolor=None
        cbweight='bold'
        cbrotation=0
        cbva='bottom'
        cbha='center'
        cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
        layers=[
               [pv,cmap,8,cbar1],
               ];
    elif 'icing' in vartype:
        if vartype == 'icing_height':
            TRACE=pv[0]
            LIGHT=pv[1]
            MODERATE=pv[2]
            SEVERE=pv[3]
            icinglevs=pv[4]
            name='Height of maximum icing' #codes_get(gid,"name")
            units='m' #codes_get(gid,"units")
            bounds=icinglevs
            boundst=bounds
            pcontour=False
            pbarbs=False
            ptext=False
            clayers=None
            ticks=boundst.astype(int)
            norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            #--
            trace=plt.get_cmap('Blues_r')
            light=plt.get_cmap('Greens_r')
            moderate=plt.get_cmap('PuOr')
            severe=plt.get_cmap('Reds_r')
            colorT=trace(np.linspace(0.15,1,256))
            colorL=light(np.linspace(0.5,1,256))
            colorM=moderate(np.linspace(0.,0.43,256))
            colorS=severe(np.linspace(0.15,1,256))
            trace=mcolors.LinearSegmentedColormap.from_list('trace',colorT)
            light=mcolors.LinearSegmentedColormap.from_list('light',colorL)
            moderate=mcolors.LinearSegmentedColormap.from_list('moderate',colorM)
            severe=mcolors.LinearSegmentedColormap.from_list('severe',colorS)
            #--
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            labelsize=8
            #pvmin=pv.min()
            #pvmin=pv.max()
            clim_min=bounds.min()
            clim_max=bounds.max()
            frac=0.050
            pad=0.001
            cbar1=[extend,boundst,frac,pad,'trace','lightblue','bold',0,'bottom','left',ticks,None]
            cbar2=[extend,boundst,frac,pad,'light','#01a049','bold',0,'top','center','',None]
            cbar3=[extend,boundst,frac,pad,'moderate','#f0944d','bold',0,'bottom','center','',None]
            cbar4=[extend,boundst,frac,pad,'severe','#ec2d01','bold',0,'top','right','',None]
            layers=[
                    [TRACE,trace,8,cbar1],
                    [LIGHT,light,8,cbar2],
                    [MODERATE,moderate,9,cbar3],
                    [SEVERE,severe,11,cbar4]
                    ];
    elif 'cloud' in vartype:
        name=codes_get(gid,"name")
        if vartype == 'cloud_amount':
            hcc=8.*np.ma.masked_array(pv[0],pv[0]<0.1)
            mcc=8.*np.ma.masked_array(pv[1],pv[1]<0.1)
            lcc=8.*np.ma.masked_array(pv[2],pv[2]<0.1)
            PMSL=pv[3]
            name='Cloud Fraction'
            units=''
            bounds=np.array([0,1,2,3,4,5,6,7,8])
            boundst=bounds
            pcontour=True
            pbarbs=False
            ptext=False
            ticks=boundst.astype(int)
            norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            #--
            light=plt.get_cmap('Blues')
            colorL=light(np.linspace(0.2,0.6,256))
            light=mcolors.LinearSegmentedColormap.from_list('High',colorL)
            moderate=plt.get_cmap('Greens')
            colorM=moderate(np.linspace(0.2,0.6,256))
            moderate=mcolors.LinearSegmentedColormap.from_list('Mid',colorM)
            severe=plt.get_cmap('Greys')
            colorS=severe(np.linspace(0.2,0.6,256))
            severe=mcolors.LinearSegmentedColormap.from_list('Low',colorS)
            #--
            clayers=[PMSL/100]
            colors='red'
            linewidths=0.25
            czorder=11
            crange=range(950,1080,5)
            alpha=1.0
            extend='neither'
            labelsize=8
            #pvmin=pv.min()
            #pvmin=pv.max()
            clim_min=bounds.min()
            clim_max=bounds.max()
            frac=0.050
            pad=0.001
            cbar1=[extend,boundst,frac,pad,'High','lightblue','bold',0,'bottom','left',ticks,None]
            cbar2=[extend,boundst,frac,pad,'Mid','#01a049','bold',0,'top','center','',None]
            cbar3=[extend,boundst,frac,pad,'Low','#f0944d','bold',0,'bottom','center','',None]
            layers=[
                    [hcc,light,8,cbar1],
                    [mcc,moderate,9,cbar2],
                    [lcc,severe,10,cbar3],
                    ];
            coastColor='darkorange'
        if vartype == 'cloud_base':
            name=codes_get(gid,"name")
            units='ft / m'
            pv=pv
            bounds=np.array([100,200,500,1000,1500,5000,10000, 26000])
            boundst=bounds
            pcontour=False
            pbarbs=False
            ptext=True
            clayers=None
            ticks=['100 / 30', '200 / 60', '500 / 152', '1000 / 305', '1500 / 457', '5000 / 1524', '10000 / 3048', '26000 / 7924']
            colors=['dimgrey', 'purple', 'red', 'orange', 'yellow', 'greenyellow', 'limegreen', 'cornflowerblue', 'lightblue']
            cmap, norm = from_levels_and_colors(bounds,colors,extend='both')
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            frac=0.15
            pad=0.05
            labelsize=8
            pvmin=int(pv.min())
            pvmax=int(pv.max())
            clim_min=bounds.min()
            clim_max=bounds.max()
            step=None
            barblength=None
            linewidth=None
            weight=None
            texts=[[0.5,-0.05,"min: {} ft  max: {} ft".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
            cbtitle=''
            cbcolor=None
            cbweight=None
            cbrotation=0
            cbva='bottom'
            cbha='center'
            cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
            layers=[
                    [pv,cmap,8,cbar1],
                    ];
        if vartype == 'cloud_top':
            name=codes_get(gid,"name")
            units='ft / m'
            pv=pv
            bounds=np.array([100,500,1000,1500,5000,10000,20000,26000,30000,35000])
            boundst=bounds
            pcontour=False
            pbarbs=False
            ptext=True
            clayers=None
            ticks=['100 / 30', '500 / 152', '1000 / 305', '1500 / 457', '5000 / 1524', '10000 / 3048', '20000 / 6069', '26000 / 7924', '30000 / 9144', '35000 / 10668']
            colors=['dimgrey', 'red', 'orange', 'yellow', 'greenyellow', 'limegreen', 'royalblue', 'cornflowerblue', 'lightblue', 'plum', 'mediumorchid']
            cmap, norm = from_levels_and_colors(bounds,colors,extend='both')
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            frac=0.15
            pad=0.05
            labelsize=8
            pvmin=int(pv.min())
            pvmax=int(pv.max())
            clim_min=bounds.min()
            clim_max=bounds.max()
            step=None
            barblength=None
            linewidth=None
            weight=None
            texts=[[0.5,-0.05,"min: {} ft  max: {} ft".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
            cbtitle=''
            cbcolor=None
            cbweight=None
            cbrotation=0
            cbva='bottom'
            cbha='center'
            cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
            layers=[
                    [pv,cmap,8,cbar1],
                    ];
        if vartype == 'cloud_depth':
            name='Cloud Depth'
            units='m'
            pv=pv
            bounds=np.linspace(0,10000,51)
            boundst=np.linspace(0,10000,11)
            pcontour=False
            pbarbs=False
            ptext=False
            clayers=None
            ticks=boundst.astype(int)
            #cmap=plt.cm.gnuplot(np.linspace(0.,1,256))
            #cmap=plt.cm.get_cmap("gnuplot")
            #cmap.set_over('#c9ff27')
            norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            cmap = plt.get_cmap(colormap)
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            frac=0.15
            pad=0.05
            labelsize=8
            pvmin=int(pv.min())
            pvmax=int(pv.max())
            clim_min=bounds.min()
            clim_max=bounds.max()
            step=None
            barblength=None
            linewidth=None
            weight=None
            texts=[[0.5,-0.05,"min: {} max: {}".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
            cbtitle=''
            cbcolor=None
            cbweight=None
            cbrotation=0
            cbva='bottom'
            cbha='center'
            cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
            layers=[
                    [pv,cmap,8,cbar1],
                    ];
    elif 'CAT' in vartype:
        if vartype == 'maxCAT':
            maxcat=pv[0]
            maxcatlev=pv[1]
            lightcat=np.ma.masked_array(maxcatlev,maxcat<4)
            moderatecat=np.ma.masked_array(maxcatlev,maxcat<8)
            severecat=np.ma.masked_array(maxcatlev,maxcat<12)
            name=codes_get(gid,"name")
            units='m'
            bounds=np.array([1,3,5,6,7,8,9,10])*1000
            boundst=bounds
            pcontour=False
            pbarbs=False
            ptext=False
            clayers=None
            ticks=boundst.astype(int)
            norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            #--
            light=plt.get_cmap('Blues_r')
            colorL=light(np.linspace(0.5,1,256))
            light=mcolors.LinearSegmentedColormap.from_list('light',colorL)
            moderate=plt.get_cmap('Greens_r')
            colorM=moderate(np.linspace(0.,0.43,256))
            moderate=mcolors.LinearSegmentedColormap.from_list('moderate',colorM)
            severe=plt.get_cmap('Reds_r')
            colorS=severe(np.linspace(0.15,1,256))
            severe=mcolors.LinearSegmentedColormap.from_list('severe',colorS)
            #--
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            labelsize=8
            #pvmin=pv.min()
            #pvmin=pv.max()
            clim_min=bounds.min()
            clim_max=bounds.max()
            frac=0.050
            pad=0.001
            cbar1=[extend,boundst,frac,pad,'light','lightblue','bold',0,'bottom','left',ticks,None]
            cbar2=[extend,boundst,frac,pad,'moderate','#01a049','bold',0,'top','center','',None]
            cbar3=[extend,boundst,frac,pad,'severe','#f0944d','bold',0,'bottom','center','',None]
            layers=[
                    [lightcat,light,8,cbar1],
                    [moderatecat,moderate,9,cbar2],
                    [severecat,severe,10,cbar3]
                    ];
    elif 'Inversions' in vartype:
        if vartype == 'Inversions':
            LIGHT=np.ma.masked_array(pv[0], pv[0]<1)    #inHigh
            MODERATE=np.ma.masked_array(pv[1], pv[1]<1) #inMiddle
            SEVERE=np.ma.masked_array(pv[2], pv[2]<1)   #inLow
            name='Temp. inversion strength' #codes_get(gid,"name")
            units=u"\u2103" #codes_get(gid,"units")
            bounds=np.array([0,2,4,6,8,10])
            boundst=bounds
            pcontour=False
            pbarbs=False
            ptext=True
            clayers=None
            ticks=boundst.astype(int)
            norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
            #--
            light=plt.get_cmap('Blues')
            moderate=plt.get_cmap('Greens')
            severe=plt.get_cmap('Oranges')
            colorL=light(np.linspace(0.2,0.6,256))
            colorM=moderate(np.linspace(0.2,0.6,256))
            colorS=severe(np.linspace(0.2,0.6,256))
            light=mcolors.LinearSegmentedColormap.from_list('High',colorL)
            moderate=mcolors.LinearSegmentedColormap.from_list('Mid',colorM)
            severe=mcolors.LinearSegmentedColormap.from_list('Low',colorS)
            #--
            zorder=8
            czorder=None
            alpha=1.0
            extend='both'
            labelsize=7
            pvmin=round(np.min([LIGHT,MODERATE,SEVERE]),1)
            pvmax=round(np.max([LIGHT,MODERATE,SEVERE]),1)
            clim_min=bounds.min()
            clim_max=bounds.max()
            frac=0.070
            pad=0.01
            texts=[[0.5,-0.05,"min: {}  max: {}".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
            cbar1=[extend,boundst,frac,pad,'<300m','#f0944d','bold',0,'top','center',ticks,None]
            cbar2=[extend,boundst,frac,pad,'300-800m','#01a049','bold',0,'bottom','center','',None]
            cbar3=[extend,boundst,frac,pad,'>800m','lightblue','bold',0,'top','center','',None]
            layers=[
                    [SEVERE,severe,8,cbar1],
                    [MODERATE,moderate,9,cbar2],
                    [LIGHT,light,11,cbar3]
                    ];
    else:
        print('no settings found.. setting default plot parameters')
        name=codes_get(gid,"name")
        units=codes_get(gid,"units")
        pv=pv
        bounds=np.linspace(int(pv.min()),int(pv.max()),50)
        boundst=bounds
        pcontour=False
        pbarbs=False
        ptext=False
        clayers=None
        ticks=bounds.astype(int)
        cmap = cmap = plt.get_cmap('viridis')
        norm=None
        zorder=8
        czorder=None
        alpha=1.0
        extend='neither'
        frac=0.15
        pad=0.05
        labelsize=8
        pvmin=pv.min()
        pvmax=pv.max()
        clim_min=-0.5
        clim_max=8-0.5
        step=None
        barblength=None
        linewidth=None
        weight=None
        texts=None
        cbtitle=''
        cbcolor=None
        cbweight='bold'
        cbrotation=0
        cbva='bottom'
        cbha='center'
        cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
        layers=[
               [pv,cmap,zorder,cbar1],
               ];

    #=========================== plot area size, other plot areas can be added later
    if areatype == 'full':
        ax.coastlines(resolution="50m", color=coastColor, linewidth=0.65, alpha=0.65,zorder=1000)
        #print ("Done coastlines", datetime.datetime.now())
        ax.add_feature(cfeature.BORDERS.with_scale('50m'),facecolor='none',edgecolor=coastColor,linestyle='-.',linewidth=0.5,alpha=0.65,zorder=1000)
        #print ("Done borders", datetime.datetime.now())
        ax.add_feature(cfeature.LAKES.with_scale('10m'),facecolor='none',edgecolor=coastColor,linestyle='-',linewidth=0.5,alpha=0.75,zorder=1000)
        #print ("Done lakes", datetime.datetime.now())
        ax.gridlines(draw_labels=False, linewidth=0.5, color='black', alpha=0.75, linestyle=':',zorder=1000)
        #print ("Done gridlines", datetime.datetime.now())
    else:   #default behaviour
        print ("Unknown area size.. Doing full area")
        ax.coastlines(resolution="50m", color=coastColor, linewidth=0.65, alpha=0.65)
        print ("Done coastlines", datetime.datetime.now())
        ax.add_feature(cfeature.BORDERS.with_scale('50m'),facecolor='none',edgecolor=coastColor,linestyle='-.',linewidth=0.5,alpha=0.65)
        print ("Done borders", datetime.datetime.now())
        ax.add_feature(cfeature.LAKES.with_scale('10m'),facecolor='none',edgecolor=coastColor,linestyle='-',linewidth=0.5,alpha=0.75)
        print ("Done lakes", datetime.datetime.now())
        ax.gridlines(draw_labels=False, linewidth=0.5, color='black', alpha=0.75, linestyle=':')
        print ("Done gridlines", datetime.datetime.now())

    #=========================== begin plotting

    for layer in layers:
        im=pcolormesh(elons,elats,layer[0], transform=ccrs.PlateCarree(),cmap=layer[1],norm=norm,zorder=layer[2],alpha=alpha)
        #for cbar in cbars:
            #print (cbar)
        cbar=layer[-1]
        extend=cbar[0];bounds=cbar[1];frac=cbar[2];pad=cbar[3];cbtitle=cbar[4];cbcolor=cbar[5];cbweight=cbar[6];cbrotation=cbar[7];cbva=cbar[8];cbha=cbar[9];ticklabels=cbar[10];weight=cbar[11];
        cbar=plt.colorbar(im,extend=extend,ticks=bounds, fraction=frac, pad=pad)
        cbar.ax.set_title(cbtitle, color=cbcolor,weight=cbweight,rotation=cbrotation, va=cbva,ha=cbha)
        #print(ticklabels)
        cbar.set_ticklabels(ticklabels)
        #cbar.ax.set_yticklabels(ticklabels,weight=weight)
        plt.clim(clim_min,clim_max)
        #plt.setp(cbars[0].ax.get_yticklabels(), fontsize=labelsize)
    if pcontour:
        for cpv in clayers:
            imc=contour(lons,lats,cpv,crange,transform=ccrs.PlateCarree(),colors=colors,linewidths=linewidths,zorder=czorder)
    if pbarbs:
        imb=barbs(lons[::step,::step],lats[::step,::step],uReg[::step,::step]*constants.knotcoef,vReg[::step,::step]*constants.knotcoef,transform=ccrs.PlateCarree(),length=barblength,linewidth=linewidth,zorder=100)
    if ptext:
        for sets in texts:
            #print(sets)
            x=sets[0];y=sets[1];text=sets[2];va=sets[3];ha=sets[4];rotation=sets[5];rotation_mode=sets[6];
            ax.text(x,y,text,va=va, ha=ha, rotation=rotation, rotation_mode=rotation_mode, transform=ax.transAxes)

    plt.title(name+', '+units+', '+codes_get(gid,"typeOfLevel")+' '+ str(codes_get(gid,"level"))+'\n '+ CallesEgna.getDates(gid)[0].strftime("%Y-%m-%d %H:%M UTC")+', valid: '+CallesEgna.getDates(gid)[1].strftime("%Y-%m-%d %H:%M UTC")+' '+codes_get(gid,"stepType"),fontsize=8)
    #plt.show()
    if SaveFig:    #savePath+'MaxWind_'+Cycle+str(K).zfill(2)+'.png'
        seconds=(CallesEgna.getDates(gid)[1]-CallesEgna.getDates(gid)[0]).seconds
        leadTimeInMinutes=seconds // 60
        #hours=seconds // 3600
        #minutes = (seconds % 3600) // 60
        if picname != None:
            #savefig(savePath+picname+CallesEgna.getDates(gid)[0].strftime("%H")+CallesEgna.getDates(gid)[1].strftime("%H")+'_'+singleMember+'.png')
            #savefig(savePath+picname + '_'+ CallesEgna.getDates(gid)[0].strftime("%H").zfill(2)+str(codes_get(gid,"step")).zfill(2) + '.png')
            #savefig(savePath+picname + '_'+ CallesEgna.getDates(gid)[0].strftime("%H").zfill(2)+'+'+str(hours).zfill(2)+'h'+str(minutes).zfill(2)+'m' + '.png')
            savefig(savePath+picname + '_'+ CallesEgna.getDates(gid)[0].strftime("%H").zfill(2)+'+'+str(leadTimeInMinutes).zfill(4) + '.png')
        else:
            savefig(savePath+codes_get(gid,"shortName")+ CallesEgna.getDates(gid)[0].strftime("%H").zfill(2)+'+'+CallesEgna.getDates(gid)[1].strftime("+%Hh%Mm").zfill(2) +'_'+codes_get(gid,"stepType")+'_'+codes_get(gid,"typeOfLevel")+str(codes_get(gid,"level"))+'.png')
        plt.close()


###############################################################################################
###############################################################################################

#New Colormaps:

#creating a new colormap for MLH
Mcol1=plt.cm.Purples(np.linspace(0.2,1,26))
Mcol2=plt.cm.Blues(np.linspace(0.4,1,26))
Mcol3=plt.cm.Greens(np.linspace(0.4,1,26))
Mcol4=plt.cm.Reds(np.linspace(0.4,1,51))
#Mcol5=plt.cm.spring(np.linspace(0,1,76))
Mcol5=plt.cm.YlOrBr_r(linspace(0.5,1,127))
MLHcolors=np.vstack((Mcol1,Mcol2,Mcol3,Mcol4,Mcol5))
myMLHmap=mcolors.LinearSegmentedColormap.from_list('my_MLHmap',MLHcolors)
myMLHmap.set_over('#ffffff')
plt.register_cmap(cmap=myMLHmap)

#creating a new colormap for Wmax

windcolor1=plt.cm.rainbow(np.linspace(0.,1,255))
windcolor2=plt.cm.binary(0)
windcolor=np.vstack((windcolor2,windcolor1))
cmap_wind=mcolors.LinearSegmentedColormap.from_list('my_colormapwind',windcolor)
cmap_wind.set_over('#742802')
cmap_wind.set_under('#fe019a')
plt.register_cmap(cmap=cmap_wind)

#creating a new colormap for Temperature
color1=plt.cm.YlGn(np.linspace(0.05,1,40))
color2=plt.cm.Purples(np.linspace(0.2,1,44))
color3=plt.cm.Blues(np.linspace(0.2,1,44))
color4=plt.cm.Reds(np.linspace(0.,1,44))
color5=plt.cm.YlOrRd(np.linspace(0.2,0.65,44))
color6=plt.cm.copper_r(np.linspace(0.2,0.95,40))
colors=np.vstack((color1,color2,color3,color4,color5,color6))
mymap=mcolors.LinearSegmentedColormap.from_list('my_colormap',colors)
mymap.set_under('#ff08e8')
mymap.set_over('#000000')
plt.register_cmap(cmap=mymap)

#creating a new colormap for TemperatureZoom
color11=plt.cm.YlGn(np.linspace(0.05,1,42))
color22=plt.cm.Purples(np.linspace(0.2,1,43))
color33=plt.cm.Blues(np.linspace(0.2,1,43))
color44=plt.cm.Reds(np.linspace(0.,1,43))
color55=plt.cm.YlOrRd(np.linspace(0.2,0.65,43))
color66=plt.cm.copper_r(np.linspace(0.2,0.95,42))
colorsnew=np.vstack((color11,color22,color33,color44,color55,color66))
mymapnew=mcolors.LinearSegmentedColormap.from_list('my_colormapnew',colorsnew)
mymapnew.set_under('#ff08e8')
mymapnew.set_over('#03012d')
plt.register_cmap(cmap=mymapnew)

#creating a new colormap for Windspeed
Col1=plt.cm.RdBu(np.linspace(0.6,1,64))
Col2=plt.cm.YlOrBr(np.linspace(0.2,0.7,64))
Col3=plt.cm.Reds(np.linspace(0.3,1,64))
Col4=plt.cm.copper_r(np.linspace(0.35,0.8,64))
Cols=np.vstack((Col1,Col2,Col3,Col4))
cmapU=mcolors.LinearSegmentedColormap.from_list('my_colormapUwind',Cols)
cmapU.set_over('#000000')
cmapU.set_under('#61e160')
plt.register_cmap(cmap=cmapU)

#creating a new colormap for MU-CAPE
cmapCAPE=mcolors.LinearSegmentedColormap.from_list('my_capecolormap',plt.cm.cubehelix_r(np.linspace(0.1,1,15)))
plt.register_cmap(cmap=cmapCAPE)

#creating a new colormap for RAIN
rainColors=mcolors.LinearSegmentedColormap.from_list('my_rainColorMap',plt.cm.summer_r(np.linspace(0,1,18)))
snowColors=mcolors.LinearSegmentedColormap.from_list('my_snowColorMap',plt.cm.cool(np.linspace(0,1,18)))
plt.register_cmap(cmap=rainColors)
plt.register_cmap(cmap=snowColors)

#creating a new colormap for HAIL
hailColors=mcolors.LinearSegmentedColormap.from_list('my_hailcolormap',plt.cm.hot_r(np.linspace(0.1,1,15)))
plt.register_cmap(cmap=hailColors)

#creating a new colormap for CloudDepth
cdepthcolor1=plt.cm.Blues(np.linspace(0.,1,128))
cdepthcolor2=plt.cm.Purples(np.linspace(0.,1,128))
Cols=np.vstack((cdepthcolor1,cdepthcolor2))
cmap_cdepth=mcolors.LinearSegmentedColormap.from_list('my_colormapcdepth',Cols)
plt.register_cmap(cmap=cmap_cdepth)


#creating a new colormap for LIGHTNING
licolor1=plt.cm.jet(np.linspace(0.,1,236))
licolor2=plt.cm.rainbow(np.linspace(0.,0.1,20))
licolor=np.vstack((licolor2,licolor1))
cmap_light=mcolors.LinearSegmentedColormap.from_list('my_colormaplight',licolor)
cmap_light.set_over('k')
cmap_light.set_under('w')
plt.register_cmap(cmap=cmap_light)

#creating a new colormap for Helicopter triggered lightning index
hticolors=plt.cm.YlOrBr(np.linspace(0,2,255))
cmap_hti=mcolors.LinearSegmentedColormap.from_list('my_colormaphti',hticolors)
cmap_hti.set_under('w')
cmap_hti.set_over('k')
plt.register_cmap(cmap=cmap_hti)

###############################################################################################
###############################################################################################

#-CF-dev  fcCycle='10'
#-CF-dev fcStep='006h15m'
#fcCycle='12'
#fcStep='006h15m'
#singleMember='mbr000'
fcCycle=sys.argv[1]
fcStep=sys.argv[2]
#singleMember=sys.argv[3]

# fileKey='/metcoop/transfers/' + fcCycle + '/fc????????' + fcCycle + '+' + fcStep + 'grib2*_' + singleMember
#-FC-dev fileKey='/home/forteliu/metcoop/transfers_mnwc/production/' + fcCycle + '/fc????????' + fcCycle + '+' + fcStep + 'grib2*'
fileKey='/metcoop/transfers_mnwc/production/' + fcCycle + '/fc????????' + fcCycle + '+' + fcStep + 'grib2*'
print(fileKey)
#-FC-dev wrkdir='/home/forteliu/Documents/Python/' + fcCycle +'/'
wrkdir='/data/hirlam2/Python_MNWC/' + fcCycle +'/'
gribMapPath=wrkdir+'gribmaps/'
savePath=wrkdir+ 'weathermaps/'

#fileKey=sys.argv[1]     #'/local_disk/fc2021111203+003grib2*_mbr000'
typesOfLevel=['any']    #['hybrid','isobaricInhPa','heightAboveGround']
levels = ['any']        #[float(0),float(2)]
shortNames= ['any']     #['t','cape']
stepTypes= ['any']      #['instant']
perturbationNumbers= ['any']
dataDates=['any']
dataTimes=['any']
validityDates=['any']
validityTimes=['any']
selector= 'latlonplot'  #'gribMapping'
#saveGribMap=True
mapPath=wrkdir+'geography/METCOOP25D_i.map'

mapResolution='i'
makeMap=True
cmap = 'viridis'       
figureSize=(13.5,12)

gribFill = 9999
listAll =  False 
doAll = False    #Ahto, can choose whether to manually select fields for plotting or have all fields in previous selections done in a loop.

constants=Constants(gravity=9.81, Tmelt=273.15, TfC=-52, TcC=52, TfK=221.15, TcK=291.15, knotcoef=1.94384) #gravity constant, Celcius vs Kelvin-conversion, Min and Max values for Temperature colorbar in C and K

plt.ion()

selected=modelVariable(mapPath=mapPath,resolution=mapResolution, makeMap=makeMap,labellocs=[1,0,0,1],
            keys=ReadKey(typesOfLevel=typesOfLevel, shortNames=shortNames, levels=levels, stepTypes=stepTypes))

if selector == "gribMapping":
    gribFiles=[]
    for ifi,infile in enumerate(sort(glob.glob(fileKey))):
        print('mapping grib file ',infile)
        thisFile=gribFile(filePath=infile,gribMap=mapGrib(infile))
        gribFiles.append(thisFile)
        if saveGribMap:
            print ('dumping map as '+gribMapPath+infile.split('/')[-1]+'.map')
            pickle.dump(thisFile,open(gribMapPath+infile.split('/')[-1]+'.map','wb'),-1)

elif selector == "latlonplot":
    for ifi, infile in enumerate(sort(glob.glob(fileKey))[::-1]):
        #Need to go through files in refverse order ( [::-1]) to ensure that fp-file comes first, as needed to define Zhyb, see below
        print (datetime.datetime.now())
        print ("Mka-1")
        print (ifi)
        print (infile)
        print ("Mka-2")
        try:
            thisFile=pickle.load(open(gribMapPath+infile.split('/')[-1]+'.map','rb'))
            print('Loaded '+gribMapPath+infile.split('/')[-1]+'.map')
            print ("Done loading", datetime.datetime.now())
        except:
            print('Loading gribmap failed, mapping '+infile)
            thisFile=gribFile(filePath=infile,gribMap=mapGrib(infile))
            #pickle.dump(thisFile,open(gribMapPath+infile.split('/')[-1]+'.map','wb'),-1)
            print ("Done mapping", datetime.datetime.now())
        if infile != thisFile.filePath:
            print('Invalid file or directory!',infile)
            raise ValueError (infile, filePath)
        grib=open(thisFile.filePath,'rb')
        
        #====================================== 
        if doAll:
            matchingMessages=CallesEgna.findMessage(thisFile.gribMap, shortNames=shortNames,typesOfLevel=typesOfLevel, levels=levels,  perturbationNumbers=perturbationNumbers, stepTypes=stepTypes)
            print ('found ', len(matchingMessages),' matching messages')
            print ("Done matching", datetime.datetime.now())

            for message in matchingMessages:
                pv = getPV2(message,release=False)
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,vartype=(message.shortName+'_'+message.stepType+'_'+str(int(message.level))))
                codes_release(gid)
#======================================
        else:         #plotting of specific weather maps based on specific demands
            template=grib.seek(int((CallesEgna.findMessage(thisFile.gribMap, shortNames=['t'], typesOfLevel=['heightAboveGround'], levels=[float(2)],  perturbationNumbers=['any'], stepTypes=['instant']))[0].offset))
            gid = codes_grib_new_from_file(grib)
            try:
                m3
            except:
                m3,lons,lats,elons,elats,xy,exy=mapCreator(gid,selected.mapPath,resolution=selected.resolution,makeMap=True)
                print ("done setup m3", datetime.datetime.now())
 
#            if 'grib2_fp_mbr' in infile:
            if 'grib2_fp' in infile:
                print ("Reading file:",infile, datetime.datetime.now())
                ## CAPE
                try:
                   pv,gid = getPV('cape','heightAboveGround',0,'any','instant')
                   mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'MU-CAPE','full','my_capecolormap',picname='MU-CAPE')
                   codes_release(gid)
                except:
                   print("No CAPE found in this file")
                pv,gid = getPV('pres','heightAboveGround',0,'any','instant') #just to get a template as cape is missing in the grib
                Zhyb = np.zeros((1,pv.shape[0],pv.shape[1]))
                #Zhyb = np.zeros((1,pv.shape[0],pv.shape[1]))
                for lev in range(30,66,1):
                    pv,gid = getPV('z','hybrid',lev,'any','instant');codes_release(gid)
                    Zhyb = np.vstack((Zhyb,np.array([pv])))
                Zhyb=Zhyb[1:]

            #elif 'grib2_mbr' in infile:
            elif 'grib2' in infile:
                print ("Reading file:",infile, datetime.datetime.now())
                ## 2m Temperatures
                pv,gid = getPV('t','heightAboveGround',2,'any','instant')
            
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'temperature','full','my_colormap',picname='T_2mC')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'temperature_zoom','full','my_colormapnew',picname='T_2mCzoom')
                codes_release(gid)

                ## WINDS
                u10m,gid = getPV('u','heightAboveGround',10,'any','instant');codes_release(gid)
                v10m,gid = getPV('v','heightAboveGround',10,'any','instant')
                
                uvspeed=np.sqrt(u10m**2 + v10m**2)
                ureg,vreg=CallesEgna.RegularWind(u10m,v10m,lons,lats)
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,[ureg,vreg,uvspeed],'windspeed','full','my_colormapUwind',picname='WindSPEED')
                codes_release(gid)

                for hyblev in range(1,66,1):
                    if hyblev == 1:
                        wmax,gid = getPV('w','hybrid',hyblev,'any','instant') #vertical velocity
                    else:
                        pv,gid = getPV('w','hybrid',hyblev,'any','instant')
                        wmax=np.where(pv>wmax,pv,wmax)
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,wmax,'winds_wmax','full','BuPu',picname='Wmax')
                codes_release(gid)
                # maxWind
                umax,gid = getPV('maxucol','maxWind',0,'any','instant');codes_release(gid)
                vmax,gid = getPV('maxvcol','maxWind',0,'any','instant');codes_release(gid)
                #umax,gid = getPV('u','maxWind',0,'any','instant');codes_release(gid)
                #vmax,gid = getPV('v','maxWind',0,'any','instant');codes_release(gid)
                ureg,vreg=CallesEgna.RegularWind(umax,vmax,lons,lats)
                pv,gid = getPV('lmxws','maxWind',0,'any','instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,[umax,vmax,pv],'max_winds','full',picname='MaxWind')
                codes_release(gid)

                ## PRECIPITATIONS
                pmsl,gid = getPV('pres','heightAboveSea',0,'any','instant');codes_release(gid)
                rain_inst,gid = getPV('rain','heightAboveGround',0,'any','instant');codes_release(gid)
                snow_inst,gid = getPV('snow','heightAboveGround',0,'any','instant');codes_release(gid)
                grpl_inst,gid = getPV('grpl','heightAboveGround',0,'any','instant')
                
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,[rain_inst,snow_inst,grpl_inst,pmsl],'precipitation_inst','full',['my_rainColorMap','my_snowColorMap'],picname='PrecInst')
                codes_release(gid)

                snow_accum,gid = getPV('snow','heightAboveGround',0,'any','accum');codes_release(gid)
                grpl_accum,gid = getPV('grpl','heightAboveGround',0,'any','accum');codes_release(gid)
                rain_accum,gid = getPV('rain','heightAboveGround',0,'any','accum')
            
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,[rain_accum,snow_accum,grpl_accum,pmsl],'precipitation_accum_solid','full',['my_rainColorMap','my_snowColorMap'],picname='PrecSolid')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,[rain_accum,snow_accum,grpl_accum,pmsl],'precipitation_accum_total','full',['my_rainColorMap','my_snowColorMap'],picname='PrecAcc')
                codes_release(gid)
                
                pv,gid = getPV('xhail','heightAboveGround',0,'any','instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'hail','full','my_hailcolormap',picname='Hail')
                codes_release(gid)

                ## MIXED LAYER HEIGHT
                pv,gid = getPV('mld','heightAboveGround',0,'any','instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'MLH','full','my_MLHmap',picname='MLH')
                codes_release(gid)

                ## VISIBILITY
                pv,gid = getPV('vis','heightAboveGround',0,'any','instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'visibility','full',picname='Visibility')
                codes_release(gid)

                ## PTYPE
                pv,gid = getPV('prtp','heightAboveGround',0,'any','instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'ptype','full',picname='PrecType')
                codes_release(gid)

                ## LIGHTNING
                pv,gid = getPV('lgt','entireAtmosphere',0,'any','instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'lightning','full','my_colormaplight',picname='Lightning')
                codes_release(gid)

                ## Helicopter induced lightning
                pv,gid = getPV('hti','entireAtmosphere',0,'any','instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,[pv,pmsl],'HTI','full','my_colormaphti',picname='HTI')
                codes_release(gid)
                
                ## ICING INDEX
                iceheights=np.array([150, 300,750,1500,2250,3000, 3750,4500,5250,6000])
                
                for lev in iceheights:
                    if lev == 150:
                        icei2,gid = getPV('icei2','heightAboveGround',lev,'any','instant')
                        icei2all = np.array([icei2]) #save indexes from all heights
                        iceilevels = np.zeros((len(iceheights),icei2.shape[0],icei2.shape[1]))
                        for i, level in enumerate(iceheights):
                            iceilevels[i,:,:]=level

                        codes_release(gid)
                    else:
                        pv,gid = getPV('icei2','heightAboveGround',lev,'any','instant')
                        icei2=np.where(pv>icei2,pv,icei2)   #max icing
                        icei2all = np.vstack((icei2all,np.array([pv])))

                TRACE=np.ma.masked_array(iceilevels,icei2all!=1.0)     #Icing Index = 1
                LIGHT=np.ma.masked_array(iceilevels,icei2all!=2.0)     #Icing Index = 2
                MODERATE=np.ma.masked_array(iceilevels,icei2all!=3.0)  #Icing Index = 3
                SEVERE=np.ma.masked_array(iceilevels,icei2all!=4.0)    #Icing Index = 4
                TRACEmin=np.min(TRACE,axis=0)
                LIGHTmin=np.min(LIGHT,axis=0)
                MODERATEmin=np.min(MODERATE,axis=0)
                SEVEREmin=np.min(SEVERE,axis=0)

                mapDraw2(m3,lons,lats,elons,elats,xy,exy,[TRACEmin,LIGHTmin,MODERATEmin,SEVEREmin,iceheights],'icing_height','full',picname='IcingHeight')
                codes_release(gid)

                ## CLOUDS
                # CLOUD AMOUNT
                hcc,gid = getPV('hcc','heightAboveGround',0,'any','instant');codes_release(gid)
                mcc,gid = getPV('mcc','heightAboveGround',0,'any','instant');codes_release(gid)
                lcc,gid = getPV('lcc','heightAboveGround',0,'any','instant')
                
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,[hcc,mcc,lcc,pmsl],'cloud_amount','full',picname='CloudAmount')
                codes_release(gid)
                # CLOUD BASE
                cb,gid = getPV('cb','entireAtmosphere',0,'any','instant')
                cbf=cb/0.3048 # in feet [ft]
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,cbf,'cloud_base','full',picname='CloudBase')
                codes_release(gid)
                # CLOUD TOP
                ct,gid = getPV('ct','entireAtmosphere',0,'any','instant')
                ctf=ct/0.3048 # in feet [ft]
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,ctf,'cloud_top','full',picname='CloudTop')
                #codes_release(gid)
                # CLOUD DEPTH
                pv = ct - cb
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,cbf,'cloud_depth','full','my_colormapcdepth',picname='CloudDepth')
                codes_release(gid)
                
                ## CAT
                cat_max,gid = getPV('cat_max','entireAtmosphere',0,'any','instant');codes_release(gid)
                cat_maxlev,gid = getPV('cat_maxlev','entireAtmosphere',0,'any','instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,[cat_max,cat_maxlev],'maxCAT','full',picname='MaxCAT')
                codes_release(gid)

                ## INVERSIONS
                Thyb = np.zeros((1,pv.shape[0],pv.shape[1]))
                for lev in range(30,66,1):
                    pv,gid = getPV('t','hybrid',lev,'any','instant');codes_release(gid)
                    Thyb = np.vstack((Thyb,np.array([pv])))
                Zs,gid = getPV('z','heightAboveGround',0,'any','instant')
                Thyb=Thyb[1:]
                null = np.zeros((Zs.shape))
                inHigh = np.zeros((Zs.shape))
                inMiddle = np.zeros((Zs.shape))
                inLow = np.zeros((Zs.shape))
                highLimit=800*constants.gravity
                lowLimit=300*constants.gravity
                dt = Thyb[-2]-Thyb[-1]
                notInversionIndeces = dt <= 0; dt[notInversionIndeces]=0.
                for kz in range(Zhyb.shape[0]-1,-0,-1):
                        dt=Thyb[kz-1,:,:]-Thyb[kz,:,:]
                        notInversionIndeces = dt <= 0; dt[notInversionIndeces]=0.
                        altitude = Zhyb[kz,:,:]-Zs[:,:]
                        highSelector=null*0
                        highIndices = altitude >= highLimit
                        highSelector[highIndices]=1.
                        lowSelector=null*0
                        lowIndices = altitude <= lowLimit
                        lowSelector[lowIndices]=1.
                        middleSelector=null + 1.
                        middleSelector[highIndices]=0.
                        middleSelector[lowIndices]=0.

                        inHigh= inHigh + dt*highSelector
                        inLow= inLow + dt*lowSelector
                        inMiddle= inMiddle + dt*middleSelector
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,[inHigh,inMiddle,inLow],'Inversions','full',picname='Tinv')
                codes_release(gid)
                del Thyb,Zhyb
                print ("Plots done!", datetime.datetime.now())
            else:
                print('Skipping unknown file:',infile)
        grib.close()

