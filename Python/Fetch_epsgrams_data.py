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
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
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
    #buf=(subprocess.check_output(['grib_ls','-p','count,offset:i,shortName,typeOfLevel,level,perturbationNumber,stepType,dataDate,dataTime,validityDate,validityTime,indicatorOfParameter,edition',gribFile]).splitlines())
    
    buf=subprocess.run(['grib_ls','-p','count,offset:i,shortName,typeOfLevel,level,perturbationNumber,stepType,dataDate,dataTime,validityDate,validityTime,indicatorOfParameter,edition',gribFile],check=True,stdout=subprocess.PIPE).stdout.splitlines()
    
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
                perturbationNumber=int(sample[5]),
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
    #Northern Sweden
    Location(longitude= 18.8200, latitude= 68.3500, name='Abisko'),
    Location(longitude= 20.3300, latitude= 67.8200, name='Kiruna'),
    Location(longitude= 20.8200, latitude= 67.1300, name='Gallivare'),
    Location(longitude= 20.1400, latitude= 66.5000, name='Jokkmokk'),
    Location(longitude= 19.2200, latitude= 66.3200, name='Tjakaape'),
    Location(longitude= 22.3600, latitude= 66.5300, name='Overkalix'),
    Location(longitude= 20.1300, latitude= 65.8700, name='Vidsel'),
    Location(longitude= 19.2700, latitude= 65.5800, name='Arvidsjaur'),
    Location(longitude= 22.1200, latitude= 65.5400, name='Kallax'),
    Location(longitude= 17.6900, latitude= 64.9600, name='Gunnarn'),
    Location(longitude= 21.0700, latitude= 64.6200, name='Skelleftea'),
    Location(longitude= 18.7200, latitude= 64.5500, name='Lycksele'),
    Location(longitude= 20.2800, latitude= 63.7800, name='Umea'),
    Location(longitude= 18.9800, latitude= 63.4000, name='Ornskoldsvik'),
    Location(longitude= 14.6590, latitude= 63.1910, name='Ostersund'),
    Location(longitude= 14.5000, latitude= 63.1800, name='Are-Ostersunds'),
    Location(longitude= 17.7700, latitude= 63.0300, name='Kramfors-Solleftea'),
    Location(longitude= 17.5000, latitude= 62.5500, name='Sundsvall-Timra'),
    Location(longitude= 22.1200, latitude= 65.5500, name='Lulea'),
    Location(longitude= 14.4200, latitude= 62.0500, name='Sveg'),
    Location(longitude= 12.8300, latitude= 61.1600, name='Salen'),
    Location(longitude= 13.7300, latitude= 60.6600, name='Malung'),
    Location(longitude= 16.9500, latitude= 60.5900, name='Gavle-Sandviken'),
    Location(longitude= 15.5200, latitude= 60.4200, name='Borlange'),
    Location(longitude= 13.7200, latitude= 61.3800, name='Trangslet'),
    #Southern Sweden
    Location(longitude= 17.6000, latitude= 59.8800, name='Uppsala-Arna'),
    Location(longitude= 18.7000, latitude= 59.7300, name='Norrtalje'),
    Location(longitude= 17.9200, latitude= 59.6500, name='Arlanda'),
    Location(longitude= 16.6300, latitude= 59.5800, name='Vasteras'),
    Location(longitude= 13.3400, latitude= 59.4400, name='Karlstad'),
    Location(longitude= 17.9400, latitude= 59.3500, name='Bromma'),
    Location(longitude= 14.7800, latitude= 59.3500, name='Kilsbergen'),
    Location(longitude= 21.5000, latitude= 59.3300, name='Nordostra-Ostersjon'),
    Location(longitude= 15.0300, latitude= 59.2200, name='Orebro'),
    Location(longitude= 16.1800, latitude= 58.5850, name='Norrkoping'),
    Location(longitude= 18.0480, latitude= 59.3400, name='Stockholm'),
    Location(longitude= 15.6570, latitude= 60.6210, name='Falun'),
    Location(longitude= 13.5100, latitude= 58.9500, name='Flygsektor-S92'),
    Location(longitude= 15.7200, latitude= 58.8300, name='Flygsektor-M84'),
    Location(longitude= 16.9000, latitude= 58.7800, name='Norrkoping-Kungsangen'),
    Location(longitude= 16.2500, latitude= 58.5900, name='Skavsta'),
    Location(longitude= 14.5100, latitude= 58.5100, name='Karlsborgs'),
    Location(longitude= 18.3200, latitude= 58.5000, name='Flygsektor-M23'),
    Location(longitude= 13.9700, latitude= 58.4600, name='Skovde'),
    Location(longitude= 12.7100, latitude= 58.4300, name='Satenas'),
    Location(longitude= 15.6800, latitude= 58.4100, name='Linkoping-SAAB'),
    Location(longitude= 15.5300, latitude= 58.4000, name='Malmens'),
    Location(longitude= 12.3500, latitude= 58.3100, name='Trollhattan-Vanersborgs'),
    Location(longitude= 17.2200, latitude= 57.8000, name='Flygsektor-M21'),
    Location(longitude= 11.8700, latitude= 57.7700, name='Save'),
    Location(longitude= 14.0700, latitude= 57.7500, name='Jonkopings'),
    Location(longitude= 12.2800, latitude= 57.6600, name='Landvetter'),
    Location(longitude= 18.3500, latitude= 57.6600, name='Visby'),
    Location(longitude= 14.1300, latitude= 57.2900, name='Hagshults'),
    Location(longitude= 13.9200, latitude= 56.9500, name='Ljungby'),
    Location(longitude= 14.7200, latitude= 56.9200, name='Vaxjo'),
    Location(longitude= 15.4500, latitude= 56.8500, name='Kosta'),
    Location(longitude= 12.8200, latitude= 56.6800, name='Halmstads'),
    Location(longitude= 16.2800, latitude= 56.6800, name='Kalmar'),
    Location(longitude= 17.3300, latitude= 56.5000, name='Flygsektor-S15'),
    Location(longitude= 12.8500, latitude= 56.3000, name='Angelholm'),
    Location(longitude= 15.2700, latitude= 56.2700, name='Kallinge'),
    Location(longitude= 13.2300, latitude= 56.0800, name='Ljungbyheds'),
    Location(longitude= 14.0800, latitude= 55.9200, name='Kristianstads'),
    Location(longitude= 13.3800, latitude= 55.5400, name='Malmo-Sturup'),
    Location(longitude= 19.5000, latitude= 55.1700, name='Sydostra-Ostersjon'),
    Location(longitude= 20.0000, latitude= 58.9400, name='Flygsektor-M32'),
    Location(longitude= 16.5200, latitude= 57.7800, name='Vasterviks'),
    #Other
    Location(longitude= 32.4000, latitude= 67.4000, name='Kuola_NPP'),
    Location(longitude= 30.3000, latitude= 59.9000, name='St-Petersburg'),
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
    Location(longitude= 8.68000, latitude= 62.9650, name='Surnadal')
    ]
    return(locations)

def getData(fileKey, locations):
    print ("Start", datetime.datetime.now())    
    hybList=['z'] #list of shortnames for model level fields
    presList=[''] #list of shortnames for isobaricInhPa instant variables
    surfList=['z','pres','prtp','tcc','vis'] #list of shortnames for heightAboGround instant variable at level 0
    cls2List=['t','q'] #list of shortnames for instant heightAboGround variable at level 2 
    cls10List=['u','v'] #list of shortnames for instant heightAboGround variable at level 10 
    clsmaxList=['ugst','vgst'] #list of shortnames for max heightAboGround variable at level 10
    surfAccumList=['rain','snow','grpl'] #list of shortnames for heightAboGround accum variable at level 0
    integratedList=['lgt','cb','cat_max','cat_b','cat_t','cat_maxlev','lmxice']
    icingIndexList = ['icei2']
    fileList = sort(glob.glob(fileKey))

    #look for members:
    differentMembers=[]
    for f in filter(lambda x:('mbr' in x[-6:-2]),fileList):
        if f[-6:] not in differentMembers:
                differentMembers.append(f[-6:])
    differentMembers.sort()
    print("Found ",len(differentMembers), "members: ", differentMembers)

    #look for steps:
    differentSteps=[]
    for f in fileList:
        ts= f[f.index('+')-10:f.index('grib')]
        if ts not in differentSteps:
                differentSteps.append(ts)
    differentSteps.sort()
    print("Found ",len(differentSteps), "steps: ", differentSteps)

    nzMax=66
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
        xv=xVars()
        for vname in integratedList:
            setattr(xv,vname,np.zeros(dimSingleLevel))
        location.integratedVars=xv
        xv=xVars()
        for vname in icingIndexList:
            setattr(xv,vname,np.zeros((len(differentMembers),10,len(differentSteps))))
        location.icingIndexVars=xv       

    for ifi, infile in enumerate(fileList):
        print('mapping '+infile)
        thisFile=gribFile(filePath=infile,gribMap=mapGrib(infile), verbose=False)
        if infile != thisFile.filePath:
            raise ValueError (infile, filePath)
        grib=open(thisFile.filePath,'rb')
        if 'fp' in thisFile.filePath:
            matchingMessages=CallesEgna.findMessage(thisFile.gribMap, shortNames=[hybList[0]], typesOfLevel=['hybrid'], levels=['any'], 
                perturbationNumbers=perturbationNumbers, stepTypes=['instant'],verbose=False)
        else:
            matchingMessages=CallesEgna.findMessage(thisFile.gribMap, shortNames=hybList[1:], typesOfLevel=['hybrid'], levels=['any'], 
                perturbationNumbers=perturbationNumbers, stepTypes=['instant'],verbose=False)
        iMember=differentMembers.index(infile[-6:])
        iTime=differentSteps.index(infile[infile.index('+')-10:infile.index('grib')])
        for message in matchingMessages:
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
                    #print('finding indeces for ',location.name,' at ',location.longitude,' E ', location.latitude,' N')
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
            for location in locations:
                location.hybVars.zhalf[iMember, NLEV,iTime] =  location.surfVars.z[iMember,iTime] #values[location.jyp,location.ixp]
                for iz in range(NLEV-1,-1,-1):
                    location.hybVars.zhalf[iMember, iz, iTime]=2*location.hybVars.z[iMember, iz, iTime] - location.hybVars.zhalf[iMember, iz, iTime]
                for iz,a in enumerate(hybridCoefficients[0:NLEV]):
                    location.hybVars.pfull[iMember, iz, iTime]= 0.5*(hybridCoefficients[iz+1]+a + (hybridCoefficients[NLEV+1+iz+1]+hybridCoefficients[NLEV+1+iz])*location.surfVars.pres[iMember, iTime])
                    location.hybVars.density[iMember, iz, iTime]=   (hybridCoefficients[iz+1]-a + (hybridCoefficients[NLEV+1+iz+1]-hybridCoefficients[NLEV+1+iz])*location.surfVars.pres[iMember, iTime]) /( (location.hybVars.zhalf[iMember,iz,iTime]-location.hybVars.zhalf[iMember,iz+1,iTime]) )
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
            for message in (CallesEgna.findMessage(thisFile.gribMap, shortNames=integratedList, typesOfLevel=['entireAtmosphere'], levels=[0.], perturbationNumbers=perturbationNumbers, stepTypes=['instant'],verbose=False)):
                grib.seek(int(message.offset))
                gid = codes_grib_new_from_file(grib)
                values=np.ma.masked_equal(codes_get_array(gid,"values"),gribFill).reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni"))
                for location in locations:
                    #print(codes_get(gid,"shortName"))
                    getattr(location.integratedVars,(codes_get(gid,"shortName")))[iMember, iTime] =  values[location.jyp,location.ixp]
                codes_release(gid)
            icingLevels=[150,300,750,1500,2250,3000,3750,4500,5250,6000]
            for message in (CallesEgna.findMessage(thisFile.gribMap, shortNames=icingIndexList, typesOfLevel=['heightAboveGround'], levels=icingLevels, perturbationNumbers=perturbationNumbers, stepTypes=['instant'],verbose=False)):
                grib.seek(int(message.offset))
                gid = codes_grib_new_from_file(grib)
                values=np.ma.masked_equal(codes_get_array(gid,"values"),gribFill).reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni"))
                for location in locations:
                    getattr(location.icingIndexVars,(codes_get(gid,"shortName")))[iMember,icingLevels.index(int(message.level)),iTime] =  values[location.jyp,location.ixp]
                codes_release(gid)
        grib.close()
    return(locations)


###############################################################################################
###############################################################################################
fcCycle=sys.argv[1]
try:
    readData=sys.argv[2]
    print("readData set to ",readData, " in call")
except:
    readData="readData"
    print("readData = ",readData, " by default")
fileKey='/metcoop/transfers/' + fcCycle + '/fc????????' + fcCycle + '+???grib2*_*'
wrkdir='/data/hirlam2/Python/'
savePath=wrkdir + fcCycle +'/EPSgrams/'

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


locations=defineLocations()
locations=getData(fileKey, locations)
pickle.dump(locations,open(savePath + 'EPSG_locations.pickle','wb'),-1)
