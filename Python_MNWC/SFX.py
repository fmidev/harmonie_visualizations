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

def cf_truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    #code for creating a new colormap from an existing colormap
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
    'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
    cmap(np.linspace(minval, maxval, n)))
    return new_cmap

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
                 perturbationNumber=int(sample[5]),
                 stepType=bytes.decode(sample[6]),
                 dataDate=int(sample[7]), dataTime=int(sample[8]), validityDate=int(sample[9]), validityTime=int(sample[10]),
                 indicatorOfParameter=bytes.decode(sample[11])))
             except:
                 print('Skipping: ',  sample)
         else:
             print ('Not a valid etition number: ',sample)
     return(thisMap)

def getPV(indicatorOfParameter,typeOfLevel,level,stepType):
    '''Manual selector for plotting'''
    print("Getting field:",indicatorOfParameter,typeOfLevel,"level:",level,stepType,datetime.datetime.now())
    if type(level) == range or type(level) == list or type(level) == ndarray:  # checks if it's a list/range/array selection of levels
        field=grib.seek(int((CallesEgna.findMessage(thisFile.gribMap, shortNames=['any'], indicatorsOfParameter=[indicatorOfParameter], typesOfLevel=[typeOfLevel], levels=level, stepTypes=[stepType],perturbationNumbers=['any']))[0].offset))
    else:   # assumes it's single level
        field=grib.seek(int((CallesEgna.findMessage(thisFile.gribMap, shortNames=['any'], indicatorsOfParameter=[indicatorOfParameter], typesOfLevel=[typeOfLevel], levels=[float(level)], stepTypes=[stepType],perturbationNumbers=['any']))[0].offset))
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
    #print('Plotting:',vartype)
    plt.figure(figsize=figureSize)
    #print ("Done plt.figure(figzize=figureSize)", datetime.datetime.now())
    ax = plt.axes(projection=m3)
    #print ("Done plt.axes",datetime.datetime.now())
    ax.set_extent([exy[0,0,0],exy[-1,-1,0],exy[0,0,1],exy[-1,-1,1]],crs=m3)
    #print ("Done set_extend",datetime.datetime.now())

    #========================== vartype decisions
    if vartype == 'T2ISBA' or 'TDeep' in vartype or 'TsPatch' in vartype:
        pv=pv-constants.Tmelt
        if 'TDeep' in vartype:
            bounds=np.linspace(-30,30,61)
            boundst=np.linspace(-30,30,31)
        else:
            bounds=np.linspace(-52,52,53)
            boundst=np.linspace(-52,52,27)
        name=picname #codes_get(gid,"name")
        units=u"\u2103"
        pvmin=round(pv.min(),2)
        pvmax=round(pv.max(),2)
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
    elif vartype == 'SST':
        name='Water surface temperature' #codes_get(gid,"name")
        units=u"\u2103"
        sst=pv[0]-constants.Tmelt
        tsw=pv[1]-constants.Tmelt
        bounds_old=np.linspace(0,30,11)
        bounds=np.insert(bounds_old,0,[-20,-10,-2])
        boundst=bounds
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
        texts=[[0.5,-0.1,"Sea min: {} \u2103 Sea max: {} \u2103 \nLake min: {} \u2103 Lake max: {} \u2103".format(round(sst.min(),2),round(sst.max(),2),round(tsw.min(),2),round(tsw.max(),2)),'bottom','center','horizontal','anchor']]
        cbtitle=''
        cbcolor=None
        cbweight=None
        cbrotation=0
        cbva='bottom'
        cbha='center'
        cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
        layers=[
               [sst,cmap,8,cbar1],
               [tsw,cmap,9,None]
               ];
    elif vartype == 'SIC':
        name=codes_get(gid,"name")
        units=codes_get(gid,"units") #u"\u2103"
        pvmin=round(pv.min(),2)
        pvmax=round(pv.max(),2)
        pv=pv
        bounds=np.arange(0.0,1.2,0.20)
        boundst=bounds
        pcontour=False
        pbarbs=False
        ptext=True
        clayers=None
        ticks=boundst
        cmap=plt.get_cmap(colormap)
        norm=None #mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
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
    elif 'Frac' in vartype:  #RootSoilIceFrac
        name=picname #codes_get(gid,"name")
        units='%' #codes_get(gid,"units") #u"\u2103"
        pv=pv*100
        pvmin=round(pv.min(),2)
        pvmax=round(pv.max(),2)
        bounds=np.array([5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95])
        boundst=bounds
        pcontour=False
        pbarbs=False
        ptext=True
        clayers=None
        ticks=boundst
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
    elif 'SnowDensity' in vartype:
        name=picname #codes_get(gid,"name")
        units=r' $\mathrm{kg\ m^{-3}}$' #codes_get(gid,"units") #u"\u2103"
        pv=pv
        pvmin=round(pv.min(),2)
        pvmax=round(pv.max(),2)
        bounds=np.linspace(100,275,8)
        boundst=bounds
        pcontour=False
        pbarbs=False
        ptext=True
        clayers=None
        ticks=boundst.astype(int)
        cmap=plt.get_cmap(colormap)
        new_cmap = cf_truncate_colormap(cmap, 0.2, 0.8)
        new_cmap.set_over('#d3494e')
        new_cmap.set_under('#5d06e9')
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
        texts=[[0.5,-0.05,"min: {} max: {}".format(pvmin, pvmax),'bottom','center','horizontal','anchor']]
        cbtitle=''
        cbcolor=None
        cbweight=None
        cbrotation=0
        cbva='bottom'
        cbha='center'
        cbar1=[extend,boundst,frac,pad,cbtitle,cbcolor,cbweight,cbrotation,cbva,cbha,ticks,weight]
        layers=[
               [pv,new_cmap,8,cbar1],
               ];
    elif 'SnowDepth' in vartype or 'Flake' in vartype:
        name=picname #codes_get(gid,"name")
        if 'Flake' in vartype:
            units=r' $\mathrm{m}$' #codes_get(gid,"units")
        else:
            units='m'
        pv=pv
        pvmin=round(pv.min(),2)
        pvmax=round(pv.max(),2)
        bounds=[     0.01,    0.03,    0.1,      0.3,        1,       3,      10]
        colors=[    'w', '#0485d1','#448ee4','#95d0fc','#d5ffff','#61e160','#c1fd95', '#ffffb6']
        boundst=bounds
        pcontour=False
        pbarbs=False
        ptext=True
        clayers=None
        ticks=boundst
        cmap, norm = from_levels_and_colors(bounds,colors,extend='both')
        zorder=8
        czorder=None
        alpha=1.0
        extend='both'
        frac=0.15 #0.050
        pad=0.05 #0.001
        labelsize=8
        clim_min=np.array(bounds).min()
        clim_max=np.array(bounds).max()
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
    elif 'SnowWater' in vartype:
        name=picname #codes_get(gid,"name")
        units=r' $\mathrm{kg\ m^{-2}}$' #codes_get(gid,"units")
        pv=pv
        pvmin=round(pv.min(),2)
        pvmax=round(pv.max(),2)
        bounds=[     0.01,      3,    10,      30,        100,       300,      1000]
        colors=[    'w', '#0485d1','#448ee4','#95d0fc','#d5ffff','#61e160','#c1fd95', '#ffffb6']
        boundst=bounds
        pcontour=False
        pbarbs=False
        ptext=True
        clayers=None
        ticks=boundst
        cmap, norm = from_levels_and_colors(bounds,colors,extend='both')
        zorder=8
        czorder=None
        alpha=1.0
        extend='both'
        frac=0.15 #0.050
        pad=0.05 #0.001
        labelsize=8
        clim_min=np.array(bounds).min()
        clim_max=np.array(bounds).max()
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
    elif 'SoilIceS' in vartype or 'SoilIceR' in vartype or 'SoilWater' in vartype:
        name=picname #codes_get(gid,"name")
        units=r' $\mathrm{m^{3}\ m^{-3}}$' #codes_get(gid,"units")
        pv=pv
        pvmin=round(pv.min(),2)
        pvmax=round(pv.max(),2)
        bounds=np.array([0.001, 0.01, 0.03, 0.05, 0.07, 0.09, 0.15, 0.2, 0.25, 0.3, 0.5, 0.75])
        boundst=bounds
        pcontour=False
        pbarbs=False
        ptext=True
        clayers=None
        ticks=boundst
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
        ax.coastlines(resolution="50m", color='grey', linewidth=0.65, alpha=0.65,zorder=1000)
        #print ("Done coastlines", datetime.datetime.now())
        ax.add_feature(cfeature.BORDERS.with_scale('50m'),facecolor='none',edgecolor='grey',linestyle='-.',linewidth=0.5,alpha=0.65,zorder=1000)
        #print ("Done borders", datetime.datetime.now())
        ax.add_feature(cfeature.LAKES.with_scale('10m'),facecolor='none',edgecolor='grey',linestyle='-',linewidth=0.5,alpha=0.75,zorder=1000)
        #print ("Done lakes", datetime.datetime.now())
        ax.gridlines(draw_labels=False, linewidth=0.5, color='black', alpha=0.75, linestyle=':',zorder=1000)
        #print ("Done gridlines", datetime.datetime.now())
    else:   #default behaviour
        print ("Unknown area size.. Doing full area")
        ax.coastlines(resolution="50m", color='grey', linewidth=0.65, alpha=0.65)
        ax.coastlines(resolution="50m", color='grey', linewidth=0.65, alpha=0.65)
        print ("Done coastlines", datetime.datetime.now())
        ax.add_feature(cfeature.BORDERS.with_scale('50m'),facecolor='none',edgecolor='grey',linestyle='-.',linewidth=0.5,alpha=0.65)
        print ("Done borders", datetime.datetime.now())
        ax.add_feature(cfeature.LAKES.with_scale('10m'),facecolor='none',edgecolor='grey',linestyle='-',linewidth=0.5,alpha=0.75)
        print ("Done lakes", datetime.datetime.now())
        ax.gridlines(draw_labels=False, linewidth=0.5, color='black', alpha=0.75, linestyle=':')
        print ("Done gridlines", datetime.datetime.now())

    #=========================== begin plotting

    for layer in layers:
        im=pcolormesh(elons,elats,layer[0], transform=ccrs.PlateCarree(),cmap=layer[1],norm=norm,zorder=layer[2],alpha=alpha)
        #for cbar in cbars:
            #print (cbar)
        cbar=layer[-1]
        if cbar != None:
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

    try:
        plt.title(name+', '+units+', '+codes_get(gid,"typeOfLevel")+' '+ str(codes_get(gid,"level"))+'\n '+ CallesEgna.getDates(gid)[0].strftime("%Y-%m-%d %H:%M UTC")+', valid: '+CallesEgna.getDates(gid)[1].strftime("%Y-%m-%d %H:%M UTC")+' '+codes_get(gid,"stepType")+', mbr'+str(codes_get(gid,"perturbationNumber")),fontsize=8)
    except:
        plt.title(name+', '+units+', '+codes_get(gid,"typeOfLevel")+' '+ str(codes_get(gid,"level"))+'\n '+ CallesEgna.getDates(gid)[0].strftime("%Y-%m-%d %H:%M UTC")+', valid: '+CallesEgna.getDates(gid)[1].strftime("%Y-%m-%d %H:%M UTC")+' '+codes_get(gid,"stepType"),fontsize=8)
    #plt.show()
    if SaveFig:    #savePath+'MaxWind_'+Cycle+str(K).zfill(2)+'.png'
        if picname != None:
            print('saving as: ' + savePath+picname+ '+' +str(codes_get(gid,"step")).zfill(3)+'_' + singleMember + '_' + CallesEgna.getDates(gid)[0].strftime("%H").zfill(2)+'.png')
            savefig(savePath+picname+ '+' +str(codes_get(gid,"step")).zfill(3)+'_' + singleMember + '_' + CallesEgna.getDates(gid)[0].strftime("%H").zfill(2)+'.png')
        else:
            savefig(savePath+codes_get(gid,"shortName")+str(codes_get(gid,"step"))+'_'+codes_get(gid,"stepType")+'_'+codes_get(gid,"typeOfLevel")+str(codes_get(gid,"level"))+'_mbr'+str(codes_get(gid,"perturbationNumber"))+'.png')
        plt.close()


###############################################################################################
###############################################################################################

#New Colormaps:

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

#creating a new colormap for Soil water content
Scolor1=plt.cm.hsv_r(np.linspace(0.25,0.65,140))
Scolor2=plt.cm.hsv_r(np.linspace(0.75,1,116))
Scolors=np.vstack((Scolor1,Scolor2))
Soil_cmap=mcolors.LinearSegmentedColormap.from_list('Soil_cmap',Scolors)
Soil_cmap.set_under('#bc13fe')
Soil_cmap.set_over('#fe019a')
plt.register_cmap(cmap=Soil_cmap)

SSTcol1=plt.cm.Blues_r(np.linspace(0.5,0.95,60))
SSTcol2=plt.cm.YlOrRd(np.linspace(0,1,155))
SSTcol3=plt.cm.RdPu(np.linspace(0.5,0.8,41))
SSTcolors=np.vstack((SSTcol1,SSTcol2,SSTcol3))
mySSTmap=mcolors.LinearSegmentedColormap.from_list('my_SSTmap',SSTcolors) #creating a new colormap
mySSTmap.set_over('#7e1e9c')
mySSTmap.set_under('#030aa7')
plt.register_cmap(cmap=mySSTmap)

###############################################################################################
###############################################################################################

typesOfLevel=['heightAboveGround','meanSea']    #['hybrid','isobaricInhPa','heightAboveGround']
levels = ['any']        #[float(0),float(2)]
shortNames= ['any']     #['t','cape']
stepTypes= ['any']      #['instant']
indicatorsOfParameter=['any']
perturbationNumbers= ['any']
dataDates=['any']
dataTimes=['any']
validityDates=['any']
validityTimes=['any']
selector= 'latlonplot'  #'gribMapping'
mapResolution='i'

fcCycle=sys.argv[1]
fcStep=sys.argv[2]
singleMember=sys.argv[3]

# fileKey='/metcoop/transfers/' + fcCycle + '/fc????????' + fcCycle + '+' + fcStep + 'grib_sfxs_' + singleMember
fileKey='/metcoop/transfers_mnwc/production/' + fcCycle + '/fc????????' + fcCycle + '+' + fcStep + 'h00m' + 'grib_sfxs'
wrkdir='/data/hirlam2/Python_MNWC/' + fcCycle +'/'
gribMapPath=wrkdir+'gribmaps/'
savePath=wrkdir+ 'sfxfields/'
mapPath=wrkdir+'geography/METCOOP25D_i.map'

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
    for ifi, infile in enumerate(sort(glob.glob(fileKey))):
        print (datetime.datetime.now())
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
            matchingMessages=CallesEgna.findMessage(thisFile.gribMap, indicatorsOfParameter=indicatorsOfParameter, shortNames=shortNames,typesOfLevel=typesOfLevel, levels=levels, stepTypes=stepTypes)
            print ('found ', len(matchingMessages),' matching messages')
            print ("Done matching", datetime.datetime.now())

            for message in matchingMessages:
                pv = getPV2(message,release=False)
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv)
                codes_release(gid)
#======================================
        else:         #plotting of specific weather maps based on specific demands
            template=grib.seek(int((CallesEgna.findMessage(thisFile.gribMap,indicatorsOfParameter=[11] , typesOfLevel=['heightAboveGround'], levels=[float(802)], stepTypes=['instant'],shortNames=['any'],perturbationNumbers=['any']))[0].offset))
            gid = codes_grib_new_from_file(grib)
            try:
                m3
            except:
                m3,lons,lats,elons,elats,xy,exy=mapCreator(gid,selected.mapPath,resolution=selected.resolution,makeMap=True)
                print ("done setup m3", datetime.datetime.now())
            
            if 'sfxs' in infile:
                print ("Reading file:",infile, datetime.datetime.now())
                ## T2ISBA
                pv,gid = getPV(11,'heightAboveGround',802,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'T2ISBA','full','my_colormap',picname='T2ISBA')
                codes_release(gid)

                ## SST
                tsw,gid = getPV(11,'heightAboveGround',770,'instant');codes_release(gid)
                sst,gid = getPV(11,'meanSea',0,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,[sst,tsw],'SST','full','my_SSTmap',picname='SST')
                codes_release(gid)

                ## SeaIceConcentration
                pv,gid = getPV(91,'meanSea',0,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'SIC','full','Blues_r',picname='SeaIceConcentration')
                codes_release(gid)

                ## RootSoilIceFracPatch1
                pv,gid = getPV(193,'heightAboveGround',831,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'RootSoilIceFrac','full','Soil_cmap',picname='RootSoilIceFracPatch1')
                codes_release(gid)

                ## RootSoilIceFracPatch2
                pv,gid = getPV(193,'heightAboveGround',841,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'RootSoilIceFrac','full','Soil_cmap',picname='RootSoilIceFracPatch2')
                codes_release(gid)

                ## SnowDensityPatch1
                sdp1,gid = getPV(191,'heightAboveGround',830,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,sdp1,'SnowDensity','full','jet',picname='SnowDensityPatch1')
                codes_release(gid)

                ## SnowDensityPatch2
                sdp2,gid = getPV(191,'heightAboveGround',840,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,sdp2,'SnowDensity','full','jet',picname='SnowDensityPatch2')
                codes_release(gid)

                ## SnowWaterEquivalentPatch1
                swep1,gid = getPV(13,'heightAboveGround',830,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,swep1,'SnowWater','full','none',picname='SnowWaterEquivalentPatch1')
                codes_release(gid)

                ## SnowWaterEquivalentPatch2
                swep2,gid = getPV(13,'heightAboveGround',840,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,swep2,'SnowWater','full','none',picname='SnowWaterEquivalentPatch2')
                #codes_release(gid)

                ## SnowDepthPatch1
                pv=swep1/sdp1
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'SnowDepth','full','none',picname='SnowDepthPatch1')
                #codes_release(gid)

                ## SnowDepthPatch2
                pv=swep2/sdp2
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'SnowDepth','full','none',picname='SnowDepthPatch2')
                codes_release(gid)

                ## SoilIceRootPatch1
                pv,gid = getPV(193,'heightAboveGround',831,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'SoilIceRoot','full','Soil_cmap',picname='SoilIceRootPatch1')
                codes_release(gid)

                ## SoilIceRootPatch2
                pv,gid = getPV(193,'heightAboveGround',841,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'SoilIceRoot','full','Soil_cmap',picname='SoilIceRootPatch2')
                codes_release(gid)

                ## SoilIceSurfPatch1
                sisp1,gid = getPV(193,'heightAboveGround',830,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,sisp1,'SoilIceSurf','full','Soil_cmap',picname='SoilIceSurfPatch1')
                codes_release(gid)

                ## SoilIceSurfPatch2
                sisp2,gid = getPV(193,'heightAboveGround',840,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,sisp2,'SoilIceSurf','full','Soil_cmap',picname='SoilIceSurfPatch2')
                codes_release(gid)

                ## SoilWaterRootPatch1
                pv,gid = getPV(86,'heightAboveGround',831,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'SoilWater','full','Soil_cmap',picname='SoilWaterRootPatch1')
                codes_release(gid)

                ## SoilWaterRootPatch2
                pv,gid = getPV(86,'heightAboveGround',841,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'SoilWater','full','Soil_cmap',picname='SoilWaterRootPatch2')
                codes_release(gid)

                ## SoilWaterSurfPatch1
                swsp1,gid = getPV(86,'heightAboveGround',830,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,swsp1,'SoilWater','full','Soil_cmap',picname='SoilWaterSurfPatch1')
                codes_release(gid)

                ## SoilWaterSurfPatch2
                swsp2,gid = getPV(86,'heightAboveGround',840,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,swsp1,'SoilWater','full','Soil_cmap',picname='SoilWaterSurfPatch2')
                #codes_release(gid)

                ## SurfaceSoilIceFracPatch1
                pv=sisp1/(swsp1+sisp1)
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'SurfaceSoilIceFrac','full','Soil_cmap',picname='SurfaceSoilIceFracPatch1')
                #codes_release(gid)

                ## SurfaceSoilIceFracPatch2
                pv=sisp2/(swsp2+sisp2)
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'SurfaceSoilIceFrac','full','Soil_cmap',picname='SurfaceSoilIceFracPatch2')
                codes_release(gid)

                ## TDeepPatch1
                pv,gid = getPV(11,'heightAboveGround',831,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'TDeep','full','my_colormap',picname='TDeepPatch1')
                codes_release(gid)

                ## TDeepPatch2
                pv,gid = getPV(11,'heightAboveGround',841,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'TDeep','full','my_colormap',picname='TDeepPatch2')
                codes_release(gid)

                ## TsPatch1
                pv,gid = getPV(11,'heightAboveGround',830,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'TsPatch','full','my_colormap',picname='TsPatch1')
                codes_release(gid)

                ## TsPatch2
                pv,gid = getPV(11,'heightAboveGround',840,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'TsPatch','full','my_colormap',picname='TsPatch2')
                codes_release(gid)

                ## IceThicknessFlake
                pv,gid = getPV(18,'heightAboveGround',780,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'IceThicknessFlake','full','none',picname='IceThicknessFlake')
                codes_release(gid)

                ## SnowDepthFlake
                pv,gid = getPV(17,'heightAboveGround',780,'instant')
                mapDraw2(m3,lons,lats,elons,elats,xy,exy,pv,'SnowDepthFlake','full','none',picname='SnowDepthFlake')
                codes_release(gid)

                print ("Plots done!", datetime.datetime.now())
            else:
                print('Skipping unknown file:',infile)
        grib.close()

