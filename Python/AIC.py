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
##import pyproj
import pickle
#import epygram
#import unittest2,os,tempfile,sys,glob,subprocess,multiprocessing,time,random
import os,tempfile,sys,glob,subprocess,multiprocessing,time,random
from pkg_resources import parse_version
#from cdo import Cdo,CDOException,CdoTempfileStore
import  getopt 
from eccodes import *

## matplotlib related
## import matplotlib
## matplotlib.use('QT4Agg')
#from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
#import matplotlib
import matplotlib.pyplot as plt
from pylab import *
#import matplotlib.axes as ax
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.dates as mdates
from matplotlib.colors import from_levels_and_colors, BoundaryNorm
from matplotlib.ticker import MaxNLocator
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import CallesEgna
import matplotlib
matplotlib.use('Agg')   # Enable to prevent opening graphics in batch mode

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

class MyClass(object):
    ## No need for dot syntax
    class_var = 1

    def __init__(self, i_var):
        self.i_var = i_var

## Need dot syntax as we've left scope of class namespace
MyClass.class_var
## 1


class Location:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class xVars:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class Constants:
	def __init__(self, **kwds):
		self.__dict__.update(kwds)

# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()

# Recast levels to new class
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
    print(verbose)
    buf=(subprocess.check_output([
        #'grib_ls','-p','count:i,offset:i,shortName,typeOfLevel,level:s,perturbationNumber:i,stepType,dataDate:i,dataTime:i,validityDate:i,validityTime:i,indicatorOfParameter:i,edition:i', gribFile]).splitlines())
        'grib_ls','-p','count,offset:i,shortName,typeOfLevel,level,perturbationNumber,stepType,dataDate,dataTime,validityDate,validityTime,indicatorOfParameter,edition', gribFile]).splitlines())
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
                print('Skipping edition 1: ',  sample)
        elif edition==2:
            try:
                thisMap.append(attributes(
                filepath=gribFile,
                count=int(bytes.decode(sample[0])),
                offset=float(bytes.decode(sample[1])),
                shortName=bytes.decode(sample[2]),
                typeOfLevel=bytes.decode(sample[3]),
                level=float(bytes.decode(sample[4])),
                perturbationNumber=bytes.decode(sample[5]),
                stepType=bytes.decode(sample[6]),
                dataDate=int(bytes.decode(sample[7])), dataTime=int(bytes.decode(sample[8])), 
                validityDate=int(bytes.decode(sample[9])), validityTime=int(bytes.decode(sample[10])),
                indicatorOfParameter=bytes.decode(sample[11])))
            except:
                print('trying for multiple steps', bytes.decode(sample[6])+bytes.decode(sample[7]))
                try:
                    thisMap.append(attributes(
                    filepath=gribFile,
                    count=int(bytes.decode(sample[0])),
                    offset=float(bytes.decode(sample[1])),
                    shortName=bytes.decode(sample[2]),
                    typeOfLevel=bytes.decode(sample[3]),
                    level=float(bytes.decode(sample[4])),
                    perturbationNumber=bytes.decode(sample[5]),
                    stepType=bytes.decode(sample[6])+bytes.decode(sample[7]),
                    dataDate=int(bytes.decode(sample[8])), dataTime=int(bytes.decode(sample[9])), 
                    validityDate=int(bytes.decode(sample[10])), validityTime=int(bytes.decode(sample[11])),
                    indicatorOfParameter=bytes.decode(sample[12])))
                except:
                    print('trying for multiple steps failed')
                    print('Skipping edition 2: ',  sample)
        else:
            print ('Not a valid etition number: ',sample)
    return(thisMap)
def predefined_vmax(shortName, level):
    # Temperature 
    if shortName=='t':
        if int(level)==10:
            return float(2.5)
        else:
            return float(0.6)
    # u,v  on levels 10,20
    elif shortName in ('u','v'): 
        if int(level) in (10,20):
            return float(4)
        else:
            return float(2)
    elif shortName=='q':
        if int(level)==10:
            return 1e-6
        else:
            return 5e-4

def datestr_to_datetime(date,time):
    if time=='0':
        time='0000'
    elif len(time)==3:
        time='0'+time
    return datetime.datetime.strptime(date+time,'%Y%m%d%H%M')

def accessMessage( gribMap, shortName='pres', typeOfLevel='heightAboveGround', level='0.',perturbationNumber='any', stepType='any', verbose=False):
    import CallesEgna
    message=CallesEgna.findMessage(gribMap, shortNames=[shortName], typesOfLevel=[typeOfLevel], levels=[level], 
                perturbationNumbers=[perturbationNumber], stepTypes=[stepType],verbose=False)[0]
    grib.seek(int(message.offset))
    gid=codes_grib_new_from_file(grib)
    return(np.ma.masked_equal(codes_get_array(gid,"values"),gribFill).reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni")), gid)


fileKey=''    
typesOfLevel=''
levels = ''     
shortNames= ''  
stepTypes= ''
dataDates='any'
dataTimes='any'
validityDates='any'
validityTimes='any'
selector= '' 
mapPath= '/.m3'
mapResolution='i'
makeMap=True
cmap = 'viridis'       
figureSize=(9,8)
gribFill = 9999
listAll =  False  

#AI_ua_2022011806+000.grib
#AI_ua_2021121906-2022011806.grib
#fcCycle='06'
fcCycle=sys.argv[1]

#singleMember='mbr000'
fileKey='/metcoop/transfers/' + fcCycle + '/AI_*_*????????' + fcCycle + '*.grib*'
wrkdir='/data/hirlam2/Python/' + fcCycle +'/'
#fileKey=wrkdir + '/AI_*_*????????' + fcCycle + '*.grib' #just for testin
gribMapPath=wrkdir+'gribmaps/'
savePath=wrkdir+ 'aninc/'

#dateGroup='2021112200+000'
#dateGroup='2021102300-2021112200'

typesOfLevel=['any']
levels = ['any']          
indicatorsOfParameters= ['any']
shortNames= ['any']
perturbationNumbers=['any']
stepTypes= ['any']

#selector= 'vprofile' 
#selector='gribMapping'
#saveGribMap=True
selector='analysisIncrements'

constants=Constants(gravity=9.81, tZero=273.15)

## For figure plotting or saving figures
figuresize=(9,8) #chose you own figuresize o plots (if chosen to be shown)
SaveChart=True #choose whether the plots shall be saved (True) or not (False)

## Define the NAME of the saved weather charts(format VariableCycleTimeindex.png)
#Timeindex="%Y%m%d%H" #remove some timeparameters if you prefer shorter names for the .png file
# NOT USED HERE


CloseChart=True # Choose whether the plots shall be displayed (False) or not (True) after plotting and saving all figures
##  Set as False AND uncomment the [matplotlib.use('Agg')] at the beginning of the code to force all the figures NOT to be displayed to the user
##  NOTE: If you comment or uncomment matplotlib.use('Agg'), restart the python for defining the default settings!  

plt.ion()


if selector == "gribMapping":
    gribFiles=[]
    for ifi,infile in enumerate(sort(glob.glob(fileKey))):
        print('mapping grib file ',infile)
        thisFile=gribFile(filePath=infile,gribMap=mapGrib(infile))
        gribFiles.append(thisFile)
        if True: #saveGribMap:
            print ('dumping map as '+gribMapPath+infile.split('/')[-1]+'.map')
            pickle.dump(thisFile,open(gribMapPath+infile.split('/')[-1]+'.map','wb'),-1)

elif selector == "analysisIncrements":
    cmap='seismic'
    for ifi, infile in enumerate(sort(glob.glob(fileKey))):
        infile_split=infile.split('/')[-1]## in case there are more dash symbols in directory
        if '-' in infile_split: # average over period: 2021102300-2021112200
            dateGroup= infile_split[infile_split.index('-')-10:infile_split.index('.grib')]
        else: #instantaneous: 2021112200+000
            dateGroup= infile_split[infile_split.index('grib')-15:infile_split.index('.grib')]
        if '-' in dateGroup:
            ztype='average'
        else:
            ztype='instantaneous'

        try:
            thisFile=pickle.load(open(gribMapPath+infile.split('/')[-1]+'.map','rb'))
            print('Loaded '+gribMapPath+infile.split('/')[-1]+'.map')
        except:
            print('Loading gribmap failed, mapping '+infile)
            thisFile=gribFile(filePath=infile,gribMap=mapGrib(infile),verbose=True)
            #pickle.dump(thisFile,open(gribMapPath+infile.split('/')[-1]+'.map','wb'),-1)
        if infile != thisFile.filePath:
            raise ValueError (infile, filePath)
        grib=open(thisFile.filePath,'rb')
        if 'totua' in thisFile.filePath:
            task='TOTUA'
        elif 'ua' in thisFile.filePath:
            task='UA'
        elif 'lsmix' in thisFile.filePath:
            task='LSMIX'
        elif 'sfx' in thisFile.filePath:
            task='SFX'
        elif 'sfc' in thisFile.filePath:
            task='SFC'
        matchingMessages=CallesEgna.findMessage(thisFile.gribMap, shortNames=shortNames, 
               typesOfLevel=typesOfLevel, levels=levels,  perturbationNumbers=perturbationNumbers, stepTypes=stepTypes,verbose=False)
        print ('found ', len(matchingMessages),' matching messages')
        for message in matchingMessages:
            print('Doing ', message.shortName)
            patch=''
            layer=''
            plotThis=True
            if True: #try:
                grib.seek(int(message.offset))
                gid = codes_grib_new_from_file(grib)
                #special treatment of surfex variables:
                tplv= codes_get(gid,"typeOfLevel")
                lv= codes_get(gid,"level")
                sn= codes_get(gid,"shortName")
                if tplv=='heightAboveGround' and lv in [830,831]:
                    patch='patch 1,'
                elif tplv=='heightAboveGround' and lv in [840,841]:
                    patch='patch 2,'
                if tplv=='heightAboveGround' and lv in [830,840]:
                    layer='surf lr,'
                elif tplv=='heightAboveGround' and lv in [831,841]:
                    layer='ground,'
                if tplv=='heightAboveGround' and (lv in [760,770] or (codes_get(gid,"shortName")=='pt')): #Skip SST and LST and SWE
                    plotThis=False
                print(tplv, lv, plotThis)
                try:
                    m3
                except:
                    m3,lons,lats,elons,elats,xy,exy=mapCreator(gid,mapPath,resolution=mapResolution,makeMap=True)
                    print ("done setup m3", datetime.datetime.now())
                if plotThis:
                    pv=np.ma.masked_equal(codes_get_array(gid,"values"),gribFill).reshape(codes_get(gid,"Nj"),codes_get(gid,"Ni"))
                    #Hmax=codes_get(gid,"max")
                    Hmax=pv.max()
                    indices = np.where(pv == Hmax)
                    latmaxloc=lats[indices[0][0],indices[1][0]]
                    lonmaxloc=lons[indices[0][0],indices[1][0]]

                    #Hmin=codes_get(gid,"min")
                    Hmin=pv.min()
                    indices = np.where(pv == Hmin)
                    latminloc=lats[indices[0][0],indices[1][0]]
                    lonminloc=lons[indices[0][0],indices[1][0]]

                    #Hmean=codes_get(gid,"average")
                    #Hstd=codes_get(gid,"standardDeviation")
                    limit=np.max([abs(Hmin),abs(Hmax)])
                    if tplv=='hybrid':
                        limit=predefined_vmax(message.shortName, lv)

                    fig=plt.figure(figsize=figureSize)
                    ax = plt.axes(projection=m3)
                    ax.set_extent([exy[0,0,0],exy[-1,-1,0],exy[0,0,1],exy[-1,-1,1]],crs=m3)
                    im=pcolormesh(elons,elats, pv, transform=ccrs.PlateCarree(),vmin=-limit, vmax=limit, cmap=cmap)
                    ax.coastlines(resolution="50m", color='grey',alpha=0.6)
                    ax.add_feature(cfeature.LAKES.with_scale('10m'),facecolor='none', edgecolor='grey',alpha=0.6)
                    ax.gridlines(linestyle=':', linewidth=1, color='gray', alpha=0.6)#,draw_labels=True)
                    ax.text(0.5,-0.05, "min: {:.2E} at ({:.2F},{:.2F}) max: {:.2E} at ({:.2F},{:.2F})".format(Hmin,lonminloc,latmaxloc, Hmax,lonmaxloc,latmaxloc),
                            va='bottom', ha='center', rotation='horizontal', rotation_mode='anchor', 
                            transform=ax.transAxes, fontsize=8) 
                    plt.colorbar(im)
                    
                    try: #grib 1
                        plt.title(codes_get(gid,"shortName")+', '+str(codes_get(gid,"indicatorOfParameter"))+', '+codes_get(gid,"units")+', '+codes_get(gid,"typeOfLevel")+', '+str(codes_get(gid,"level"))+'\n '
                                +task+', '+patch+ ' '+layer+' '+ dateGroup, fontsize=8)
                    except: #grib 2
                        plt.title(codes_get(gid,"name")+', '+codes_get(gid,"units")+', '+codes_get(gid,"typeOfLevel")+' '+ str(codes_get(gid,"level")) +'\n '+
                                task+ ' ' + patch +' ' + layer + ' ' + dateGroup, fontsize=8)
                        #plt.title(codes_get(gid,"name")+', '+codes_get(gid,"units")+', '+codes_get(gid,"typeOfLevel")+' '+ str(codes_get(gid,"level"))+'\n '+
                        #   CallesEgna.getDates(gid)[0].strftime("%Y%m%d%H:%M")+
                        #    ', valid: '+CallesEgna.getDates(gid)[1].strftime("%Y-%m-%d %H:%M UTC")+' '+codes_get(gid,"stepType")+
                        #    ', mbr000',fontsize=8)
                    #plt.show()
                    codes_release(gid)
                    SaveFig=SaveChart
                    if SaveFig:
                        print('saving '+savePath+task+'_'+sn+'_'+str(lv)+'_'+ztype+'.png')
                        savefig(savePath+task+'_'+sn+'_'+str(lv)+'_'+ztype+'_'+fcCycle+'.png')
                    CloseFig=CloseChart
                    if CloseFig:
                        plt.close(fig)

            else:#except:
                print('Ecxepition handling message:', message.shortName, message.typeOfLevel, message.level, message.perturbationNumber, message.stepType)
        grib.close()

