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
from scipy import interpolate
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
import matplotlib
matplotlib.use('Agg')
#from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
#import matplotlib
import matplotlib.pyplot as plt
from pylab import *
#import matplotlib.axes as ax
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

def RegMats(lon,lat):
# Calle Fortelius January 2017, July 2017
# Returns coefficients needed to compute the true north and eastward components of a vector given relative to some projected grid
# Method: The local angle of the grid north relative to true north is approximated by examining the direction
# between correspoinding grid points on adjacent latitude rows using pyproj-Geod. While approximate, the naethod is 
# independent of the projection, that, consequently, needs not to be known.
# lon and lat are arrays of geographical latitudes and longitudes of the grid.
        g = pyproj.Geod(ellps='WGS84')

        #;az21=np.zeros((template.longitudes.shape))
        #(az12[:-1,:], az21[:-1,:], dist)=g.inv(template.longitudes[0:-1,:], template.latitudes[0:-1,:], 
        #                       template.longitudes[1:,:], template.latitudes[1:,:])
        az12=np.zeros((lon.shape)) 
        az12[:-1,:]=g.inv(lon[0:-1,:], lat[0:-1,:], lon[1:,:], lat[1:,:])[0]
        az12[-1,:]=az12[-2,:]; az12=-az12
        sinalpha=np.sin(az12*np.pi/180)
        cosalpha=np.cos(az12*np.pi/180)
        return cosalpha,sinalpha

def RegWinds(uproj,vproj,cosalpha,sinalpha):
# Calle Fortelius January 2017, July 2017
# Returns the true north and eastward components of a vector given relative to some grid
# using coefficients computed by Regmats
        return uproj*cosalpha - vproj*sinalpha, vproj*cosalpha + uproj*sinalpha

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
            pass
            #print('Skipping: ',sample)
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
            #print ('Not a valid etition number: ',sample)
    return(thisMap)

def datestr_to_datetime(date,time):
    if time=='0':
        time='0000'
    elif len(time)==3:
        time='0'+time
    return datetime.datetime.strptime(date+time,'%Y%m%d%H%M')

def defineLocations():
    #Note: Here we aqssume all members to have the same initial time.
    locations=[
    #Western Finland
    Location(longitude= 24.6000, latitude= 60.2000, name='Espoo'),
    Location(longitude= 23.6100, latitude= 60.8200, name='Forssa'),
    Location(longitude= 24.8000, latitude= 61.8600, name='Halli'),
    Location(longitude= 23.0000, latitude= 59.5000, name='Hanko_sea'),
    Location(longitude= 25.0400, latitude= 60.0700, name='Helsinki_meri'),
    Location(longitude= 25.0800, latitude= 60.2200, name='Helsinki-Itakeskus'),
    Location(longitude= 24.9500, latitude= 60.1800, name='Helsinki-Kaisaniemi'),
    Location(longitude= 24.9600, latitude= 60.2000, name='Helsinki-Kumpula'),
    Location(longitude= 25.0440, latitude= 60.2540, name='Helsinki-Malmi-EFHF'),
    Location(longitude= 24.9000, latitude= 60.3200, name='Helsinki-Vantaa'),
    Location(longitude= 24.4000, latitude= 60.2000, name='Hirlamsmossen'),
    Location(longitude= 21.3720, latitude= 60.2230, name='Houtskari'),
    Location(longitude= 24.8500, latitude= 60.6300, name='Hyvinkaa'),
    Location(longitude= 24.4000, latitude= 61.0000, name='Hameenlinna'),
    Location(longitude= 25.1000, latitude= 60.5000, name='Jarvenpaa'),
    Location(longitude= 20.5000, latitude= 61.0000, name='Joppe'),
    Location(longitude= 24.3000, latitude= 59.8690, name='Kallbadan'),
    Location(longitude= 23.0300, latitude= 63.1000, name='Kauhava'),
    Location(longitude= 22.7280, latitude= 60.1640, name='Kemio'),
    Location(longitude= 19.1670, latitude= 61.4170, name='Kilo13'),
    Location(longitude= 23.1000, latitude= 63.8000, name='Kokkola'),
    Location(longitude= 24.6830, latitude= 64.5380, name='Kopsa'),
    Location(longitude= 20.7820, latitude= 60.2560, name='Kumlinge'),
    Location(longitude= 23.7500, latitude= 61.5700, name='Lake-Nasijarvi'),
    Location(longitude= 24.1000, latitude= 60.2000, name='Lohja'),
    Location(longitude= 19.8800, latitude= 60.1200, name='Maarianhamina'),
    Location(longitude= 19.1310, latitude= 60.3030, name='Market'),
    Location(longitude= 22.4700, latitude= 61.8400, name='Niinisalo'),
    Location(longitude= 21.8100, latitude= 61.4800, name='Pori'),
    Location(longitude= 20.6670, latitude= 61.4670, name='Pori_sea'),
    Location(longitude= 21.5200, latitude= 61.1300, name='Rauma'),
    Location(longitude= 24.8000, latitude= 60.8000, name='Riihimaki'),
    Location(longitude= 24.1100, latitude= 60.7400, name='Rayskala'),
    Location(longitude= 23.1000, latitude= 60.4000, name='Salo'),
    Location(longitude= 22.8500, latitude= 62.7900, name='Seinajoki'),
    Location(longitude= 20.8160, latitude= 59.4670, name='Suomen_leijona'),
    Location(longitude= 23.4000, latitude= 60.0000, name='Tammisaari'),
    Location(longitude= 23.7840, latitude= 61.4900, name='Tampere'),
    Location(longitude= 23.6100, latitude= 61.4170, name='Tampere-EFTP'),
    Location(longitude= 22.2720, latitude= 60.4550, name='Turku'),
    Location(longitude= 22.2620, latitude= 60.5150, name='Turku-EFTU'),
    Location(longitude= 21.6700, latitude= 63.1000, name='Vaasa'),
    #Eastern Finland
    Location(longitude= 27.2000, latitude= 60.6000, name='Hamina'),
    Location(longitude= 28.8000, latitude= 61.2000, name='Imatra'),
    Location(longitude= 29.7500, latitude= 62.6300, name='Joensuu'),
    Location(longitude= 25.7300, latitude= 62.2400, name='Jyvaskyla'),
    Location(longitude= 26.6090, latitude= 63.7030, name='Kiuruvesi'),
    Location(longitude= 29.8160, latitude= 62.7210, name='Kontiolahti'),
    Location(longitude= 26.9000, latitude= 60.5000, name='Kotka'),
    Location(longitude= 26.7000, latitude= 60.9000, name='Kouvola'),
    Location(longitude= 27.6580, latitude= 62.8940, name='Kuopio'),
    Location(longitude= 27.6509, latitude= 63.0353, name='Kuopio-EFKU'),
    Location(longitude= 25.7000, latitude= 61.0000, name='Lahti'),
    Location(longitude= 25.6340, latitude= 60.9830, name='Lahti-Ski'),
    Location(longitude= 25.4000, latitude= 61.5000, name='Lake-Paijanne'),
    Location(longitude= 28.3700, latitude= 61.2500, name='Lake-Saimaa'),
    Location(longitude= 28.1000, latitude= 61.0500, name='Lappeenranta'),
    Location(longitude= 29.6290, latitude= 63.5110, name='Lieksa-Nurmes'),
    Location(longitude= 26.5300, latitude= 61.2500, name='Mantyharju'),
    Location(longitude= 27.2000, latitude= 61.7000, name='Mikkeli'),
    Location(longitude= 25.7500, latitude= 62.2900, name='Palokka'),
    Location(longitude= 25.8300, latitude= 60.4500, name='Porvoo'),
    Location(longitude= 28.2000, latitude= 62.0000, name='Rantasalmi'),
    Location(longitude= 28.8780, latitude= 61.8680, name='Savonlinna'),
    Location(longitude= 28.9320, latitude= 61.9460, name='Savonlinna-EFSA'),
    Location(longitude= 25.6700, latitude= 62.4000, name='Tikkakoski-EFJY'),
    Location(longitude= 26.9300, latitude= 60.9000, name='Utti'),
    Location(longitude= 27.8480, latitude= 62.3150, name='Varkaus'),
    Location(longitude= 27.8600, latitude= 62.1740, name='Varkaus-EFVR'),
    #Northern Finland
    Location(longitude= 23.4270, latitude= 68.3620, name='Enontekio-Hetta'),
    Location(longitude= 21.2630, latitude= 69.3100, name='Halti'),
    Location(longitude= 27.4200, latitude= 68.6200, name='Ivalo'),
    Location(longitude= 27.6800, latitude= 64.2800, name='Kajaani'),
    Location(longitude= 20.8000, latitude= 69.0500, name='Kilpisjarvi'),
    Location(longitude= 24.5800, latitude= 65.7800, name='Kemi'),
    Location(longitude= 27.3900, latitude= 66.7200, name='Kemijarvi'),
    Location(longitude= 24.8700, latitude= 67.6800, name='Kittila'),
    Location(longitude= 29.4600, latitude= 64.1100, name='Kuhmo'),
    Location(longitude= 29.2300, latitude= 65.9800, name='Kuusamo'),
    Location(longitude= 23.6800, latitude= 67.9700, name='Muonio'),
    Location(longitude= 25.3700, latitude= 64.9300, name='Oulu'),
    Location(longitude= 26.9700, latitude= 65.4000, name='Pudasjarvi'),
    Location(longitude= 25.7150, latitude= 66.5000, name='Rovaniemi'),
    Location(longitude= 25.8300, latitude= 66.5700, name='Rovaniemi-EFRO'),
    Location(longitude= 29.1410, latitude= 66.1710, name='Ruka'),
    Location(longitude= 28.6800, latitude= 66.8300, name='Salla'),
    Location(longitude= 26.6500, latitude= 67.3700, name='Sodankyla'),
    Location(longitude= 27.0200, latitude= 69.7500, name='Utsjoki'),
    Location(longitude= 28.2850, latitude= 64.1390, name='Vuokatti'),
    #Other
    Location(longitude= 17.9000, latitude= 59.5700, name='Arlanda'),
    Location(longitude= 32.4000, latitude= 67.4000, name='Kuola_NPP'),
    Location(longitude= 22.1200, latitude= 65.5500, name='Lulea'),
    Location(longitude= 30.3000, latitude= 59.9000, name='St_Petersburg'),
    Location(longitude= 24.8000, latitude= 59.4000, name='Tallinn'),
    Location(longitude= 26.7200, latitude= 58.3800, name='Tartu'),
    Location(longitude= 26.4600, latitude= 58.2600, name='Toravere'),
    Location(longitude= 26.5030, latitude= 58.0530, name='Otepaa'),
    Location(longitude= 24.1100, latitude= 56.9470, name='Riga'),
    Location(longitude= 25.2800, latitude= 54.6850, name='Vilnius'),
    Location(longitude= 12.5410, latitude= 55.6850, name='Copenhagen'),
    Location(longitude= 10.7520, latitude= 59.9460, name='Oslo'),
    Location(longitude= 10.6600, latitude= 59.9800, name='Oslo_ski'),
    Location(longitude= 10.6700, latitude= 59.9650, name='Holmenkollen'),
    Location(longitude= 10.2040, latitude= 59.7440, name='Drammen'),
    Location(longitude= 10.3100, latitude= 63.3780, name='Trondheim'),
    Location(longitude= 10.5030, latitude= 61.1340, name='Lillehammer'),
    Location(longitude= 8.68000, latitude= 62.9650, name='Surnadal'),
    Location(longitude= 16.1800, latitude= 58.5850, name='Norrkoping'),
    Location(longitude= 18.0480, latitude= 59.3400, name='Stockholm'),
    Location(longitude= 15.6570, latitude= 60.6210, name='Falun'),
    Location(longitude= 14.6590, latitude= 63.1910, name='Ostersund')
    ]
    return(locations)

def getData(fileKey, locations):
    print ("Start", datetime.datetime.now())    
    hybList=['z','t','q','tke','cwat_cond','ciwc_cond','rain_cond','snow_cond','grpl_cond'] #list of shortnames for model level fields
    presList=['u','v'] #list of shortnames for isobaricInhPa instant variables
    surfList=['z','pres','t','mld','rain','snow','grpl'] #list of shortnames for heightAboGround instant variable at level 0
    cls2List=['t','q','r'] #list of shortnames for instant heightAboGround variable at level 2 
    cls10List=['u','v'] #list of shortnames for instant heightAboGround variable at level 10 
    clsmaxList=['ugst','vgst','wsmax'] #list of shortnames for max heightAboGround variable at level 10
    surfAccumList=['rain','snow','grpl','wevap','snsub'] #list of shortnames for heightAboGround accum variable at level 0

    fileList = sort(glob.glob(fileKey))

    #look for members:
    differentMembers=[]
    for f in filter(lambda x:('mbr' in x[-6:-2]),fileList):
        if f[-6:] not in differentMembers:
                differentMembers.append(f[-6:])
    differentMembers.sort()
    #-MKa
    differentMembers=['mbr000']   # no member info in MNWC filenames => loop above finds nothing
    #-MKa
    print("Found ",len(differentMembers), "members: ", differentMembers)

    #look for steps:
    differentSteps=[]
    for f in fileList:
        ts= f[f.index('+')-10:f.index('grib')]
        if ts not in differentSteps:
                differentSteps.append(ts)
    differentSteps.sort()
    print(differentSteps)
    print("Found ",len(differentSteps), "steps: ", differentSteps)

    nzMax=66  #nzMax is the maximum allowed number of vertical levels, should be <= NLEV+1 (CF)
    dimMultiLevel=(len(differentMembers), nzMax, len(differentSteps))
    dimSingleLevel=(len(differentMembers), len(differentSteps))
    steps=np.linspace(1,len(differentSteps),len(differentSteps))
    halfsteps=np.linspace(0.5,len(differentSteps)+0.5,len(differentSteps)+1)
    for location in locations:
        location.differentSteps=differentSteps
        location.differentMembers=differentMembers
        vl=[datetime.datetime(year=1900,month=1,day=1,hour=0,minute=0) for i in range(len(differentSteps)*len(differentMembers))]
        location.initial=reshape(vl,(len(differentMembers),len(differentSteps)))
        location.validity=reshape(vl,(len(differentMembers),len(differentSteps)))
        xv=xVars()
        for vname in hybList:
            setattr(xv,vname,np.zeros(dimMultiLevel))
     #       print(location.name, vname,vname,  getattr(location,vname)[0,0,0])
        location.hybVars=xv
        location.hybVars.zhalf=np.zeros(dimMultiLevel)
        location.hybVars.pfull=np.zeros(dimMultiLevel)
        location.hybVars.density=np.zeros(dimMultiLevel)
        xv=xVars()
        for vname in surfList:
            setattr(xv,vname,np.zeros(dimSingleLevel))
        location.surfVars=xv
        xv=xVars()
        for vname in presList:
            setattr(xv,vname,np.zeros(dimMultiLevel))
        location.presVars=xv
        xv=xVars()
        for vname in cls2List+cls10List+clsmaxList:
            setattr(xv,vname,np.zeros(dimSingleLevel))
        location.clsVars=xv
        xv=xVars()
        for vname in surfAccumList:
            setattr(xv,vname,np.zeros(dimSingleLevel))
        location.surfAccumVars=xv
    for ifi, infile in enumerate(fileList):
        try:
            thisFile=pickle.load(open(gribMapPath+infile.split('/')[-1]+'.map','rb'))
            print('Loaded '+gribMapPath+infile.split('/')[-1]+'.map')
        except:
            print('Loading gribmap failed, mapping '+infile)
            thisFile=gribFile(filePath=infile,gribMap=mapGrib(infile), verbose=False)
            #pickle.dump(thisFile,open(gribMapPath+infile.split('/')[-1]+'.map','wb'),-1)
        if infile != thisFile.filePath:
            raise ValueError (infile, filePath)
        grib=open(thisFile.filePath,'rb')
        if 'fp' in thisFile.filePath:
            matchingMessages=CallesEgna.findMessage(thisFile.gribMap, shortNames=[hybList[0]], typesOfLevel=['hybrid'], levels=['any'], 
                perturbationNumbers=perturbationNumbers, stepTypes=['instant'],verbose=False)
        else:
            matchingMessages=CallesEgna.findMessage(thisFile.gribMap, shortNames=hybList[1:], typesOfLevel=['hybrid'], levels=['any'], 
                perturbationNumbers=perturbationNumbers, stepTypes=['instant'],verbose=False)
            #matchingMessages=CallesEgna.findMessage(thisFile.gribMap, shortNames=hybList, typesOfLevel=['heightAboveGround'], levels=[float(2)], 
        print ('found ', len(matchingMessages),' matching messages')
        # iMember=differentMembers.index(infile[-6:]) 
        iMember=0 #'mbr000'  # CF: iMember is an index = 0 (no member info in MNWC filenames)
        iTime=differentSteps.index(infile[infile.index('+')-10:infile.index('grib')])
        for message in matchingMessages:
            #print('Doing ',message.shortName, int(message.level))
            grib.seek(int(message.offset))
            gid = codes_grib_new_from_file(grib)
            values=np.ma.masked_equal(codes_get_array(gid,"values"),gribFill).reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni"))
            if 'lons' not in locals():
                lons,lats=codes_get_array(gid,"longitudes").reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni")),codes_get_array(gid,"latitudes").reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni"))
                cosa,sina=RegMats(lons,lats)
            if 'tree' not in locals():
                xs, ys, zs = CallesEgna.lon_lat_to_cartesian( codes_get_array(gid,"longitudes"), codes_get_array(gid,"latitudes"),R=6371000)
                tree = cKDTree(list(zip(xs, ys, zs)))
            for location in locations:
                try:
                    location.jyp
                except:
                    print('finding indeces for ',location.name,' at ',location.longitude,' E ', location.latitude,' N')
                    xpt, ypt, zpt = CallesEgna.lon_lat_to_cartesian(location.longitude,location.latitude,R=6371000)#,R=RPlanet)
                    d, ind = tree.query((xpt, ypt, zpt), k = 1)
                    location.jyp,location.ixp=np.unravel_index(ind, (codes_get(gid,"Nj"), codes_get(gid,"Ni")))
                getattr(location.hybVars, (codes_get(gid,"shortName")))[iMember, int(message.level-1), iTime] =  values[location.jyp,location.ixp]
                location.cosa=cosa[location.jyp,location.ixp]
                location.sina=sina[location.jyp,location.ixp]
                try:
                    location.hybVars.ahalf
                except:
                    hybridCoefficients=codes_get_array(gid,"pv")
                    NLEV=int(len(hybridCoefficients)/2 -1) #Number of full levels
                    location.hybVars.ahalf=hybridCoefficients[0:NLEV]
                    location.hybVars.bhalf=hybridCoefficients[NLEV:]
                    location.hybVars.fullLevels=np.linspace(1,NLEV,NLEV)
                    location.hybVars.halfLevels=np.linspace(0.5,NLEV+0.5,NLEV+1)
            codes_release(gid)
#       Here do single level fields at the surface
        if 'fp' in thisFile.filePath:
            plevs=[50,100,150,200,250,300,400,500,700,800,850,925,1000]
            for message in CallesEgna.findMessage(thisFile.gribMap, shortNames=presList, typesOfLevel=['isobaricInhPa'], levels=[300,500,850,925],perturbationNumbers=perturbationNumbers, stepTypes=['instant'],verbose=False):
                grib.seek(int(message.offset))
                gid = codes_grib_new_from_file(grib)
                values=np.ma.masked_equal(codes_get_array(gid,"values"),gribFill).reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni"))
                for location in locations:
                    #print(codes_get(gid,"shortName"))
                    getattr(location.presVars, (codes_get(gid,"shortName")))[iMember, plevs.index(int(message.level)), iTime] =  values[location.jyp,location.ixp]
                    #location.validity[iMember,iTime] = datestr_to_datetime( str(codes_get(gid,"validityDate")),str(codes_get(gid,"validityTime")) )
                    #location.initial[iMember,iTime] = datestr_to_datetime( str(codes_get(gid,"dataDate")),str(codes_get(gid,"dataTime")) )
                codes_release(gid)
        else:
            for message in CallesEgna.findMessage(thisFile.gribMap, shortNames=surfList, typesOfLevel=['heightAboveGround'], levels=[0.], perturbationNumbers=perturbationNumbers, stepTypes=['instant'],verbose=False):
                grib.seek(int(message.offset))
                gid = codes_grib_new_from_file(grib)
                values=np.ma.masked_equal(codes_get_array(gid,"values"),gribFill).reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni"))
                for location in locations:
                    #print(codes_get(gid,"shortName"))
                    getattr(location.surfVars, (codes_get(gid,"shortName")))[iMember, iTime] =  values[location.jyp,location.ixp]
                    location.validity[iMember,iTime] = datestr_to_datetime( str(codes_get(gid,"validityDate")),str(codes_get(gid,"validityTime")) )
                    location.initial[iMember,iTime] = datestr_to_datetime( str(codes_get(gid,"dataDate")),str(codes_get(gid,"dataTime")) )
                codes_release(gid)
            for message in (CallesEgna.findMessage(thisFile.gribMap, shortNames=cls2List, typesOfLevel=['heightAboveGround'], levels=[2.], perturbationNumbers=perturbationNumbers, stepTypes=['instant'],verbose=False)+
                       CallesEgna.findMessage(thisFile.gribMap, shortNames=cls10List, typesOfLevel=['heightAboveGround'], levels=[10.], perturbationNumbers=perturbationNumbers, stepTypes=['instant'],verbose=False)+
                       CallesEgna.findMessage(thisFile.gribMap, shortNames=clsmaxList, typesOfLevel=['heightAboveGround'], levels=[10.], perturbationNumbers=perturbationNumbers, stepTypes=['max'],verbose=False)):
                grib.seek(int(message.offset))
                gid = codes_grib_new_from_file(grib)
                values=np.ma.masked_equal(codes_get_array(gid,"values"),gribFill).reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni"))
                for location in locations:
                    getattr(location.clsVars,(codes_get(gid,"shortName")))[iMember, iTime] =  values[location.jyp,location.ixp]
                codes_release(gid)
            for message in (CallesEgna.findMessage(thisFile.gribMap, shortNames=surfAccumList, typesOfLevel=['heightAboveGround'], levels=[0.], perturbationNumbers=perturbationNumbers, stepTypes=['accum'],verbose=False)):
                grib.seek(int(message.offset))
                gid = codes_grib_new_from_file(grib)
                values=np.ma.masked_equal(codes_get_array(gid,"values"),gribFill).reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni"))
                for location in locations:
                    #print(codes_get(gid,"shortName"))
                    getattr(location.surfAccumVars,(codes_get(gid,"shortName")))[iMember, iTime] =  values[location.jyp,location.ixp]
                codes_release(gid)
        grib.close()
#       Do the remaining hybrid level variables
        for location in locations:
            location.hybVars.zhalf[iMember, NLEV,iTime] =  location.surfVars.z[iMember,iTime] #values[location.jyp,location.ixp]
            for iz in range(NLEV-1,-1,-1):
                location.hybVars.zhalf[iMember, iz, iTime]=2*location.hybVars.z[iMember, iz, iTime] - location.hybVars.zhalf[iMember, iz, iTime]
            for iz,a in enumerate(hybridCoefficients[0:NLEV]):
                location.hybVars.pfull[iMember, iz, iTime]= 0.5*(hybridCoefficients[iz+1]+a + (hybridCoefficients[NLEV+1+iz+1]+hybridCoefficients[NLEV+1+iz])*location.surfVars.pres[iMember, iTime])
                location.hybVars.density[iMember, iz, iTime]=   (hybridCoefficients[iz+1]-a + (hybridCoefficients[NLEV+1+iz+1]-hybridCoefficients[NLEV+1+iz])*location.surfVars.pres[iMember, iTime]) /( (location.hybVars.zhalf[iMember,iz,iTime]-location.hybVars.zhalf[iMember,iz+1,iTime]) )
                #print(a, hybridCoefficients[NLEV+1+iz], location.surfVars.pres[iMember, iTime], location.hybVars.pfull[iMember, iz, iTime])
    return(locations)

def metGrams(location):
    #location=locations[0]
    halfTimes=[location.initial[0][0]+datetime.timedelta(hours=g-0.5) for g in range(len(location.validity[0])+1)]
    NLEV=len(location.hybVars.fullLevels)
    subplts=5
    ##### PLOT PREPS
    fig, axarr = plt.subplots(subplts,sharex=True, figsize=(25*0.3, 29.5*0.3))#,constrained_layout=True)
    plevs=np.array([50,100,150,200,250,300,400,500,700,800,850,925,1000]) #5,10,11
    ticklocs=[100, 300, 500, 1000, 2000, 4000, 7000, 10000] #, 20000, 30000]
    ticklocs_for_feet = [152,305,457,762,1524,3048,6096,9144] 
    zlabs_for_feet = ['500ft','1000ft','1500ft','2500ft','5000ft','10000ft','20000ft','30000ft']
    zlabs=['0.1','0.3','0.5','1','2','4','7','10'] #,'','30']
    top=11000*constants.gravity #approximate max height above ground in the plot
    topLev= len(location.hybVars.z[0,:,0][location.hybVars.z[0,:,0] >=top])
	
    lowTop=2000*constants.gravity #approximate max height above ground in the plot
    lowtopLev= len(location.hybVars.z[0,:,0][location.hybVars.z[0,:,0] >=lowTop])

    f=interpolate.interp1d(location.hybVars.z[0,topLev:NLEV,0]/constants.gravity,
                            location.hybVars.fullLevels[topLev:NLEV],fill_value="extrapolate")

    fp=interpolate.interp1d(location.hybVars.pfull[0,topLev:NLEV,0],
                            location.hybVars.fullLevels[topLev:NLEV],fill_value="extrapolate")

    barblength=5
    linewidth=0.4
    clabelsize=9

    fig.suptitle(location.name + f' (Lon: {location.longitude:.2f}, Lat: {location.latitude:.2f}), '+location.initial[0][0].strftime("%Y-%m-%d Cycle %HZ")+", MNWC prod",size=11,y=0.99)
    ##### PANEL 1 DEFINITIONS -- Clouds (g/m3) and wind in the troposphere
    feet_axis = axarr[0].twinx()
    clev=0.002 #[0.000000001*1000]
    bounds=np.array([0.002,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,1])
    norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)

    ciwc=location.hybVars.ciwc_cond[0,topLev:NLEV,:] * 1000.*location.hybVars.density[0,topLev:NLEV,:]
    cwat=location.hybVars.cwat_cond[0,topLev:NLEV,:] * 1000.*location.hybVars.density[0,topLev:NLEV,:]
    pv=np.ma.masked_less(ciwc+cwat, clev)
    panel1=axarr[0].pcolormesh(halfTimes,location.hybVars.halfLevels[topLev:NLEV+1],pv,cmap='Blues',norm=norm)
    position=fig.add_axes([0.87,0.78,0.015,0.14])
    cbar=fig.colorbar(panel1,ax=axarr[0],cax=position,orientation='vertical',fraction=0.03, pad=0.01,aspect=20,ticks=bounds)#,shrink=0.5)
    
    cbar.ax.set_yticklabels(['0.002','0.01','0.05','0.1','0.2','0.3','0.4','0.5','0.6','1'])
    vectorstep=6


    ureg,vreg=RegWinds(location.presVars.u[0],location.presVars.v[0],location.cosa,location.sina)
    axarr[0].barbs(location.validity[0][::vectorstep], fp(plevs*100.)[[5,7,10,11]] , ureg[[5,7,10,11],::vectorstep]*constants.knotcoef, vreg[[5,7,10,11],::vectorstep]*constants.knotcoef,color='#ed0dd9',length=barblength,linewidth=linewidth,zorder=100)
    t0=axarr[0].contour(location.validity[0],location.hybVars.fullLevels[topLev:NLEV],  location.hybVars.t[0,topLev:NLEV,:]-constants.tZero,levels=[0],colors='red',linewidths=0.7,linestyles='dashed')
    try:
        axarr[0].clabel(t0, [0], inline=True, fmt=r'%.0f$^\circ$', fontsize=clabelsize)  
    except:
        pass
    if cwat.max() > clev:
        axarr[0].contour(location.validity[0],location.hybVars.fullLevels[topLev:NLEV],  cwat, [clev],colors='green')
    if ciwc.max() > clev:
        axarr[0].contour(location.validity[0],location.hybVars.fullLevels[topLev:NLEV],  ciwc, [clev],colors='blue')
    cbar.ax.set_ylabel("Cloud condensates " + r'$\mathrm{(gm^{-3})}$',size=7)
    cbar.ax.set_yticklabels(bounds,fontsize=7)
    axarr[0].set_yticks(ticks=f(ticklocs))
    axarr[0].set_yticklabels(labels=zlabs,size=7)
    #axarr[0].legend(loc='center left',bbox_to_anchor=(1,0.5))
    axarr[0].set_ylabel("altitude above sea (km/ft)",size=7)
    #axarr[0].set_title("Cloud condensates in " + r'$\mathrm{gm^{-3} }$',size=9)
    axarr[0].grid(axis='y',linestyle=':')
    axarr[0].set_axisbelow(True)
    
    feet_axis.set_yticks(ticks=f(ticklocs_for_feet))
    feet_axis.set_yticklabels(labels=zlabs_for_feet,size=7)
    feet_axis.set_ylim(axarr[0].get_ylim()[0],axarr[0].get_ylim()[1])
    feet_axis.invert_yaxis()

    ##### PANEL 2 DEFINITIONS -- Tropospheric rain, snow, and graupel as g/m3
    feet_axis2 =axarr[1].twinx()
    clev=0.01 #0.01 #[pv.min()]    
    rain = location.hybVars.rain_cond[0,topLev:NLEV,:] * 1000.*location.hybVars.density[0,topLev:NLEV,:]
    grpl = location.hybVars.grpl_cond[0,topLev:NLEV,:] * 1000.*location.hybVars.density[0,topLev:NLEV,:]
    snow = location.hybVars.snow_cond[0,topLev:NLEV,:] * 1000.*location.hybVars.density[0,topLev:NLEV,:]
    pv = np.ma.masked_less(rain+grpl+snow, clev)
    
    bounds=np.array([0.01,0.05,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3])
    norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
    panel2=axarr[1].pcolormesh(halfTimes,location.hybVars.halfLevels[topLev:NLEV+1], pv, cmap='summer_r', norm=norm)
    position=fig.add_axes([0.87,0.595,0.015,0.15])
    cbar=fig.colorbar(panel2,ax=axarr[1],cax=position,orientation='vertical',fraction=0.03, pad=0.01,aspect=20,ticks=bounds)
    #cbar.ax.set_xticklabels(['0.01','0.05','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1','1.5','2','3'])

    axarr[1].contour(location.validity[0],location.hybVars.fullLevels[topLev:NLEV],  location.hybVars.t[0,topLev:NLEV,:],[constants.tZero],colors='red',linewidths=0.7,linestyles='dashed')
    if rain.max()>clev:
        axarr[1].contour(location.validity[0],location.hybVars.fullLevels[topLev:NLEV], rain, [clev],colors='green')
    if grpl.max()>clev:
        axarr[1].contour(location.validity[0],location.hybVars.fullLevels[topLev:NLEV], grpl, [clev],colors='purple')
    if snow.max()>clev:
        axarr[1].contour(location.validity[0],location.hybVars.fullLevels[topLev:NLEV], snow, [clev],colors='blue')
    cbar.ax.set_ylabel("Rain, snow, graupel flux in " + r'$\mathrm{g m^{-3}}$',size=7)
    cbar.ax.set_yticklabels(bounds,fontsize=7)
    axarr[1].set_yticks(ticks=f(ticklocs))
    axarr[1].set_yticklabels(labels=zlabs,size=7)
    axarr[1].set_ylabel("altitude above sea (km/ft)",size=7)
    axarr[1].grid(axis='y',linestyle=':')
    axarr[1].set_axisbelow(True)
    feet_axis2.set_yticks(ticks=f(ticklocs_for_feet))
    feet_axis2.set_yticklabels(labels=zlabs_for_feet,size=7)
    feet_axis2.set_ylim(axarr[1].get_ylim()[0],axarr[1].get_ylim()[1])
    feet_axis2.invert_yaxis()
    #axarr[1].set_title('Flux of Rain, snow, graupel in ' + r'$\mathrm{g m^{-3}}$',size=9)

    ##### PANEL 3 DEFINITIONS -- Wind, precipitation and evaporation
    ureg,vreg=RegWinds(location.clsVars.u,location.clsVars.v,location.cosa,location.sina)

    RAINincrement=location.surfAccumVars.rain[0,1:]-location.surfAccumVars.rain[0,0:-1]
    SNOWincrement=location.surfAccumVars.snow[0,1:]-location.surfAccumVars.snow[0,0:-1]
    GRAUPELincrement=location.surfAccumVars.grpl[0,1:]-location.surfAccumVars.grpl[0,0:-1]
    baralign='edge'
    barwidth=(location.validity[0,1]-location.validity[0,0]).seconds/86400 #0.04
    axarr[2].bar(location.validity[0,:-1],RAINincrement+SNOWincrement+GRAUPELincrement,barwidth,align=baralign,edgecolor='none',color='lightgreen', zorder=1,label="rain")
    axarr[2].bar(location.validity[0,:-1],SNOWincrement+GRAUPELincrement,barwidth,align=baralign,edgecolor='none', color='deepskyblue', zorder=2,label="snow")
    axarr[2].bar(location.validity[0,:-1],GRAUPELincrement,barwidth,align=baralign,edgecolor='none', color='magenta', zorder=3,label="graupel")
    axarr[2].set_ylim(0.0,(RAINincrement+SNOWincrement+GRAUPELincrement).max()+1.0)
    axarr[2].legend(loc='upper left', bbox_to_anchor=(1.082, 0.8),prop={'size': 6})
    par3 = axarr[2].twinx()
    par3.plot(location.validity[0],(location.surfAccumVars.rain+location.surfAccumVars.grpl+location.surfAccumVars.snow)[0], color='green',label="precipitation")
    par3.plot(location.validity[0],(-location.surfAccumVars.wevap-location.surfAccumVars.snsub)[0], color='gray',label="evaporation")

    par3.set_ylabel("accumulation (mm)",size=7)
    par3.tick_params(labelsize=7)
    axarr[2].grid(axis='y',linewidth=0.5,linestyle=':')
    axarr[2].invert_yaxis()
    axarr[2].tick_params(labelsize=7)
    axarr[2].set_ylabel("rate (mm/h)",size=7)
    par3.legend(loc='lower left', bbox_to_anchor=(1.082, 0.2),prop={'size': 6})
    #axarr[2].set_title('Wind and Precipitation and evaporation',size=9)

    ##### PANEL 4 DEFINITIONS -- Temperature, TKE and Mixed Layer Depth
    dummy_twinx = axarr[3].twinx()
    dummy_twinx.tick_params(axis='y',left=False, right=False,labelleft=False, labelright=False)
    f=interpolate.interp1d((location.hybVars.z[0,lowtopLev:NLEV,0]-location.surfVars.z[0,0])/constants.gravity,
                            location.hybVars.fullLevels[lowtopLev:NLEV],fill_value="extrapolate")

    panel4=axarr[3].pcolormesh(halfTimes,location.hybVars.halfLevels[lowtopLev:NLEV+1],location.hybVars.t[0,lowtopLev:NLEV,:]-constants.tZero,cmap='RdYlBu_r')
    
    tkeContours=np.linspace(0.5,10,20)
    pressureContours=[300, 500, 850, 925, 1000]
    CNE=axarr[3].contour(location.validity[0],location.hybVars.fullLevels[lowtopLev:NLEV],location.hybVars.tke[0,lowtopLev:NLEV,:],levels=tkeContours[1:],colors='black',linewidths=0.3,linestyles='solid')
    CNE=axarr[3].contour(location.validity[0],location.hybVars.fullLevels[lowtopLev:NLEV],location.hybVars.tke[0,lowtopLev:NLEV,:],levels=[tkeContours[0]],colors='black',linewidths=0.6,linestyles='solid')
    CNT=axarr[3].contour(location.validity[0],location.hybVars.fullLevels[lowtopLev:NLEV],location.hybVars.t[0,lowtopLev:NLEV,:]-constants.tZero,colors='cyan',linewidths=1,linestyles='dashed',levels=[0])
    CNP=axarr[3].contour(location.validity[0],location.hybVars.fullLevels[lowtopLev:NLEV],location.hybVars.pfull[0,lowtopLev:NLEV,:]*0.01,colors='grey',linewidths=0.8,linestyles='dotted',levels=pressureContours)
    try:
        axarr[3].clabel(CNE, [tkeContours[0]], inline=True, fmt=r'%.1f', fontsize=clabelsize)  
    except:
        pass
    try:
        axarr[3].clabel(CNT, [0], inline=True, fmt=r'%.0f$^\circ$', fontsize=clabelsize)  
    except:
        pass
    position=fig.add_axes([0.84,0.235,0.015,0.15])
    cbar=fig.colorbar(panel4,ax=axarr[3],cax=position,orientation='vertical',fraction=0.03, pad=0.01,aspect=20)#,shrink=0.5)
    cbar.ax.set_ylabel("Temperature " + r'$\mathrm{(\degree C)}$',size=7)
    cbar.ax.tick_params(labelsize=7)
    try:
        axarr[3].clabel(CNP, pressureContours, inline=True, fmt=r'%.0f', fontsize=clabelsize) 
    except:
        pass
    axarr[3].plot(location.validity[0],f(location.surfVars.mld[0,:]),color='magenta',linewidth=1)
    axarr[3].set_yticks(ticks=f(ticklocs))
    axarr[3].set_yticklabels(labels=zlabs)
    axarr[3].set_ylim(lowtopLev+0.5,NLEV+0.5)
    axarr[3].set_ylabel("MLH, TKE (black lines)",size=7)
    axarr[3].tick_params(labelsize=7)
    #axarr[3].set_title('Temperature, TKE and MLH',size=9)

    par5 = axarr[4].twinx()
    WS10m=np.sqrt(location.clsVars.u**2+location.clsVars.v**2)[0]
    GUST=np.sqrt(location.clsVars.ugst**2+location.clsVars.vgst**2)[0]
    q2e=location.surfVars.pres*location.clsVars.q/(0.622+location.clsVars.q*(1.-0.622))
    Td2m=273.16+243.5*np.log(q2e/(6.112*100.))/(17.67-np.log(q2e/(6.112*100.)))
    MaxTrigger=5
    barbheight=np.mean(Td2m[0,:]-constants.tZero)
    vectorstep=3    
    par5.fill_between(location.validity[0],WS10m, GUST, where=GUST-WS10m >0,facecolor='lightblue', edgecolor='lightblue', interpolate=False, zorder=1)
    par5.fill_between(location.validity[0],WS10m+MaxTrigger, GUST, where=GUST-WS10m >= MaxTrigger,facecolor='deepskyblue', edgecolor='deepskyblue',interpolate=False,zorder=1)
    axarr[4].plot(location.validity[0],location.clsVars.t[0,:]-constants.tZero, color='darkorange',zorder=5,label="T_2m")
    axarr[4].plot(location.validity[0],Td2m[0,:]-constants.tZero, color='darkcyan',zorder=5,label="Td_2m")
    axarr[4].plot(location.validity[0],location.surfVars.t[0,:]-constants.tZero, color='r',zorder=5,label="T_s",linewidth=0.5)
    axarr[4].plot(location.validity[0],location.hybVars.t[0,NLEV-1]-constants.tZero, color='grey',zorder=5,label="T_nlev",linewidth=0.5)
    axarr[4].barbs(location.validity[0][::vectorstep], barbheight, ureg[0,::vectorstep]*constants.knotcoef, vreg[0,::vectorstep]*constants.knotcoef,color='gray',length=barblength,linewidth=linewidth+0.5,zorder=100)
    #par5.set_ylabel("T2m(o),Td2m(b),Ts(r),Tnlev(g)")
    axarr[4].legend(loc='center left', bbox_to_anchor=(1.08, 0.5),prop={'size': 6})
    axarr[4].grid(axis='y',linewidth=0.5,linestyle=':')
    axarr[4].axhline(0, color='black',linewidth=0.6)
    axarr[4].invert_yaxis()
    par5.set_ylabel("WS10m, WSMAX10m (m/s)",size=7)
    axarr[4].set_zorder(10)
    axarr[4].patch.set_visible(False)
    axarr[4].set_ylabel("Temperature " + r'$\mathrm{(\degree C)}$',size=7)
    par5.tick_params(labelsize=7)
    axarr[4].tick_params(labelsize=7)
    #axarr[4].set_title('Near surface temperature and wind',size=9)
    ##### END SUBPLOTS

    for subplt in range(subplts):
        axarr[subplt].xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
        axarr[subplt].xaxis.set_major_locator(mdates.DayLocator(interval=1))
        axarr[subplt].xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
        axarr[subplt].xaxis.set_minor_locator(mdates.HourLocator(interval=2))
        axarr[subplt].tick_params(axis='x', which='major',pad=8,labelsize=7)
        axarr[subplt].tick_params(axis='x', which='minor',labelsize=8)
        axarr[subplt].xaxis.set_ticks_position('both')
        axarr[subplt].set_xlim([location.initial[0][0], location.validity[0][-1]])
        axarr[subplt].grid(axis='x',which='minor',linewidth=0.4)
        axarr[subplt].grid(axis='x',which='major',linewidth=0.4,color='k')
        axarr[subplt].invert_yaxis()
    fig.tight_layout()
    plt.subplots_adjust(hspace=0.07,top=0.932)
    params = {'legend.fontsize': '7',
              'axes.labelsize': '7',
              'axes.titlesize':'8',
              'xtick.labelsize':'7',
              'ytick.labelsize':'7'}
    rcParams.update(params)
    axarr[0].tick_params(axis='x',which='both',labeltop=True)
    # top_axis = axarr[0].twiny()
    # top_axis.set_xticklabels(axarr[0].xaxis.get_major_formatter())
    # top_axis.set_xlabel(axarr[0].get_xlabel())
    # top_axis.set_xlim(axarr[0].get_xlim())

    plt.savefig(savePath+'Metg_' + singleMember + location.name)
    plt.close()

###############################################################################################
###############################################################################################
fcCycle=sys.argv[1]
singleMember=sys.argv[2]
try:
    readData=sys.argv[3]
    print("readData set to ",readData, " in call")
except:
    readData="readData"
    print("readData = ",readData, " by default")
# fileKey='/metcoop/transfers/' + fcCycle + '/fc????????' + fcCycle + '+???grib2*_' + singleMember
fileKey='/metcoop/transfers_mnwc/production/' + fcCycle + '/fc????????' + fcCycle + '+???h??mgrib2*'
print (fileKey)
wrkdir='/data/hirlam2/Python_MNWC/' + fcCycle +'/'
gribMapPath=wrkdir+'gribmaps/'
savePath=wrkdir+ 'meteograms/'

typesOfLevel=['heightAboveGround']
levels = ['any']
shortNames= ['lsm']
perturbationNumbers=['any']
stepTypes= ['instant']
dataDates='any'
dataTimes='any'
validityDates='any'
validityTimes='any'
selector= 'workingarea' #'gribMapping'
mapResolution='i'
makeMap=True
cmap = 'viridis'       
gribFill = 9999
listAll =  False  

constants=Constants(gravity=9.81, tZero=273.15, TfC=-52, TcC=52, TfK=221.15, TcK=291.15, knotcoef=1.94384)

#plt.ion()

if selector == "gribMapping":
    gribFiles=[]
    for ifi,infile in enumerate(sort(glob.glob(fileKey))):
        print('mapping grib file ',infile)
        thisFile=gribFile(filePath=infile,gribMap=mapGrib(infile))
        gribFiles.append(thisFile)
        if True: #saveGribMap:
            print ('dumping map as '+gribMapPath+infile.split('/')[-1]+'.map')
            pickle.dump(thisFile,open(gribMapPath+infile.split('/')[-1]+'.map','wb'),-1)

if selector == "workingarea":
    print(readData)
    if readData=="readData":
        print(readData)
        locations=defineLocations()
        locations=getData(fileKey, locations)
        pickle.dump(locations,open(savePath + 'MTG_locations' + singleMember +'.pickle','wb'),-1)
    else:
        locations=pickle.load(open(savePath + 'MTG_locations' + singleMember +'.pickle','rb'))
    for location in locations:
        try:
            print('Plotting:',location.name,location.longitude,location.latitude,datetime.datetime.now())
            metGrams(location)
        except:
            print('Unable to plot:',location.name,location.longitude,location.latitude,', Skipping..',datetime.datetime.now())
