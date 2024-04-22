#!/usr/bin/env python
from datetime import datetime, timedelta
import numpy as np
import math
import pyproj
import subprocess
from eccodes import *
import pickle

def lon_lat_to_cartesian(lon, lat, R = 1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z
    
def findMessage(gribmap, shortNames='', typesOfLevel='', levels='', perturbationNumbers='', indicatorsOfParameter=['any'], stepTypes='',verbose=False):
    matches=[] 
    for message in gribmap:
        if( ((message.shortName          in shortNames)          or shortNames==['any']) and    
            ((message.typeOfLevel        in typesOfLevel)        or typesOfLevel==['any']) and    
            ((message.level              in levels)              or levels==['any']) and
            ((message.perturbationNumber in perturbationNumbers) or perturbationNumbers==['any']) and   
            ((message.indicatorOfParameter              in indicatorsOfParameter)              or indicatorsOfParameter==['any']) and
            ((message.stepType           in stepTypes)           or stepTypes==['any'])):
            if verbose:
                print('MATCH       :', message.shortName, message.indicatorOfParameter, message.typeOfLevel, message.level, message.perturbationNumber, message.stepType)
            matches.append(message)
        #print(message.shortName, message.indicatorOfParameter, message.typeOfLevel, message.level, message.perturbationNumber, message.stepType)
    return(matches)

def getDates(gid):
    return (datetime(int(str(codes_get(gid,"dataDate"))[0:4]),
                             int(str(codes_get(gid,"dataDate"))[4:6]),
                             int(str(codes_get(gid,"dataDate"))[6:8]),
                             int(str(codes_get(gid,"dataTime")).zfill(4)[0:2]),
                             int(str(codes_get(gid,"dataTime")).zfill(4)[2:4])),
            datetime(int(str(codes_get(gid,"validityDate"))[0:4]),
                             int(str(codes_get(gid,"validityDate"))[4:6]),
                             int(str(codes_get(gid,"validityDate"))[6:8]),
                             int(str(codes_get(gid,"validityTime")).zfill(4)[0:2]),
                             int(str(codes_get(gid,"validityTime")).zfill(4)[2:4])))

def GreatCirclePath_fwd(centerlong,centerlat,azimuth,halflength,step):
# Returns longitude-latitude pairs of points along a great circle,
# with a spacing of approximately step m. The path extends equal distances
# in either direction from the given point along the given direction. 
	# calculate distance between points
	g = pyproj.Geod(ellps='WGS84')
	(endlong, endlat, az21) = g.fwd(centerlong, centerlat, azimuth, halflength)
	(startlong, startlat, az12) = g.fwd(centerlong, centerlat, azimuth+180, halflength)
	(az12, az21, dist) = g.inv(startlong, startlat, endlong, endlat)	
	# calculate line string along path with an odd number of segments ~ step
	steps=int(dist/step)
	if (steps % 2 == 0):
		steps=steps+1
	# The centered distance from the center (exact only on a spherical earth):
	distance=np.linspace( -dist/2, dist/2, steps, endpoint=True) 
	lonlats = g.npts(startlong, startlat, endlong, endlat,
	                 steps-2)
	# npts doesn't include start/end points, so prepend/append them
	lonlats.insert(0, (startlong, startlat))
	lonlats.append((endlong, endlat))
	return lonlats,distance

def GreatCirclePath_inv(startlong,startlat,endlong,endlat,step):
# Returns longitude-latitude pairs of points along a great circle,
# with a spacing of approximately step m.
	# calculate distance between points
	g = pyproj.Geod(ellps='WGS84')
	(az12, az21, dist) = g.inv(startlong, startlat, endlong, endlat)
	# calculate line string along path with segments <= 1 km
	lonlats = g.npts(startlong, startlat, endlong, endlat,
	                 1 + int(dist / step))
	# npts doesn't include start/end points, so prepend/append them
	lonlats.insert(0, (startlong, startlat))
	lonlats.append((endlong, endlat))
	return lonlats

def RegularWind(uproj,vproj,lon,lat):
# Returns true north and eastward components of a vector given relative to some projected grid
# Method: The local angle of the grid north relative to true north is approximated by examining the direction
# between correspoinding grid points on adjacent latitude rows using pyproj-Geod.
# uproj and vproj are of dimension [time, height, y, x] 
        g = pyproj.Geod(ellps='WGS84')

	#;az21=np.zeros((template.longitudes.shape))
        #(az12[:-1,:], az21[:-1,:], dist)=g.inv(template.longitudes[0:-1,:], template.latitudes[0:-1,:], 
        #                       template.longitudes[1:,:], template.latitudes[1:,:])
        az12=np.zeros((lon.shape)) 
        az12[:-1,:]=g.inv(lon[0:-1,:], lat[0:-1,:], lon[1:,:], lat[1:,:])[0]
        az12[-1,:]=az12[-2,:]; az12=-az12
        sinalpha=np.sin(az12*np.pi/180)
        cosalpha=np.cos(az12*np.pi/180)
        return uproj*cosalpha - vproj*sinalpha, vproj*cosalpha + uproj*sinalpha

def uv2ffdd(u,v):
# Works on flattened arrays!
	ff = np.sqrt(u**2 + v**2)
	dd = np.zeros(u.size)
	#for i in np.nditer(ff):
	for i in np.arange(ff.size):
		#ff[i] = np.sqrt(u[i]**2 + v[i]**2)
		dd[i]=math.atan2(u[i],v[i])*180./np.pi + 180.
		if ff[i] == 0.:
			dd[i] = 0.
	return ff,dd

def gribMap(gribs, IwantAdatelist=True):
# We assume that all grib files accord to the same stencil.
    DateList=[];Inventory=['padding']
    #grbs = pygrib.open(gribfile)
    for field in gribs[0]:
        if field.typeOfLevel == 'meanSea':
            ityp=102
        elif field.typeOfLevel == 'heightAboveSea':
            ityp=102
        elif field.typeOfLevel == 'heightAboveGround':
            ityp=105
        elif field.typeOfLevel == 'hybrid':
            ityp=109
        else:
            ityp=field.typeOfLevel
        Inventory.append( [field.indicatorOfParameter,ityp,field.level,field.timeRangeIndicator] )
        if IwantAdatelist:
            stepUnits=gribs[0][1].stepUnits
            for grib in gribs:
                analDate=grib[1].analDate
                P1=grib[1].P1
                print ('Examining ', analDate,'+',P1,' at', datetime.now().time())
                if stepUnits==1:
                    validDate=analDate + timedelta(hours=float(P1))
                elif stepUnits==0:
                    validDate=analDate + timedelta(minutes=float(P1))
                DateList.append([analDate, P1, validDate])
#                	keys=['padding']
#			for field in grib:
#                		if field.typeOfLevel == 'meanSea':
#                        		ityp=102
#                		elif field.typeOfLevel == 'heightAboveSea':
#                        		ityp=102
#                		elif field.typeOfLevel == 'heightAboveGround':
#                        		ityp=105
#                		elif field.typeOfLevel == 'hybrid':
#                        		ityp=109
#                		else:
#                        		ityp=field.typeOfLevel
#				if field.stepUnits==1:
#                			validDate=field.analDate + timedelta(hours=float(field.P1))
#				elif field.stepUnits==0:
#                			validDate=field.analDate + timedelta(minutes=float(field.P1))
#				keys.append( [field.indicatorOfParameter,ityp,field.level,field.timeRangeIndicator, field.analDate, field.P1, validDate] )
#        		Inventory.append(keys)
    return Inventory, DateList
def gribSelector(gribs, indicatorOfParameter, typeOfLevel, levelList=['*'], timeRangeIndicator=0):
        gribList=[]
        for g in gribs:
                #print datetime.now().time()
                llist=[]
                for l in g(indicatorOfParameter=indicatorOfParameter,typeOfLevel=typeOfLevel, timeRangeIndicator=timeRangeIndicator):
                        if l.level in levelList or levelList == ['*']:
                                llist.append(l)
                gribList.append(llist)
        return gribList

def MakeArray(gribList,indicatorOfParameter,typeOfLevel, levlist, timeRangeIndicator=0, jy1=0,jy2=-1,ix1=0,ix2=-1):
# parameters:
#       griblist: list of lists containing references to the desired variable
#       jy1,jy2:  start and end values of the first horizontal index (y),  (0...)
#       ix1,ix2:  start and end values of the second horizontal index (x), (0...)
#       Default values cause extraction of the entire array
    NTIME=len(gribList)
    NLEV=len(levlist)
    NLAT=gribList[0][0].Ny
    NLON=gribList[0][0].Nx
    if jy2==-1:
        jy2=NLAT-1
    if ix2==-1:
        ix2=NLON-1
    Result=np.zeros((NTIME, NLEV, jy2-jy1+1, ix2-ix1+1))
    iTime=0
    for grib in gribList:
        #print 'Time step ',iTime, datetime.now().time()
	#levlist=[]; [levlist.append(lev.level )for lev in grib]
	#levlist=sorted(levlist)
        #print 'done levlist', datetime.now().time()
        for ll in grib:
            if (ll.indicatorOfParameter==indicatorOfParameter and ll.typeOfLevel==typeOfLevel and ll.timeRangeIndicator==timeRangeIndicator and ll.level in levlist):
                Result[iTime, levlist.index(ll.level),:]=ll.values[jy1:jy2+1, ix1:ix2+1]
        iTime=iTime+1
    return Result
#def MakeArray(ipar,leveltype, gribs,levelValue='notgiven',kz1=1,kz2=-1,jy1=1,jy2=-1,ix1=1,ix2=-1,itr=0,DateList=[]):
## Extract one field from a set of grib files; return in 4-d (or 3_D) list of lists of arrays ( time, (level), y-dimension, x-dimension)
## Parameters:
##       ipar:  grib identifyer of requested variable
##       leveltype:  level type
##       gribs: list of handles on a set of grib files. It is assumed that each grib file contains data from
##               exactly one time level characterized by analysis time and forecast range. It is further 
##               assumed, that a given variable, when present,  is given on identical grids in each file.
##       levelvalue: mist be given for non-hybrid level data
##       kz1,kz2:  start and end values of the levels (1..NLEV) for hybrid levels
##       jy1,jy2:  start and end values of the first horizontal index (y),  (1..NLAT)
##       ix1,ix2:  start and end values of the second horizontal index (x), (1..NLON)
##               Note: default setting of kz2,jy2,ix2=-1 results in extracting [kz1:,jy1:,ix1:]
##       itr: grib time range indicator
##       Datelist: if present a list of dates to be treated 
#
#	ThisField=[]        
#	if DateList==[]:
#        	for grib in gribs:
#                	print grib, datetime.now().time()
#			if leveltype=='hybrid' and levelValue=='notgiven':	
#				these=grib(indicatorOfParameter=ipar,typeOfLevel=leveltype,timeRangeIndicator=itr)
#				NLEV=len(these)
#				NLAT=these[0].Ny
#				NLON=these[0].Nx
#				if kz2==-1:
#					kz2=NLEV
#				if jy2==-1:
#					jy2=NLAT
#				if ix2==-1:
#					ix2=NLON
#
#				array=np.zeros((NLEV,NLAT,NLON))
#                		for this in these:
#					array[this.level-1,:]=this.values
#                        	ThisField.append(array[kz1-1:kz2,jy1-1:jy2,ix1-1:ix2])
#			else:
#				this=grib(indicatorOfParameter=ipar,typeOfLevel=leveltype,level=levelValue,timeRangeIndicator=itr)
#				NLAT=this[0].Ny
#				NLON=this[0].Nx
#				if jy2==-1:
#					jy2=NLAT
#				if ix2==-1:
#					ix2=NLON
#				ThisField.append(this[0].values[jy1-1:jy2,ix1-1:ix2] )	
#
#	else:
#        	for date in DateList:
#                	print date, datetime.now().time()
#                	index=DateList.index(date)
#			if leveltype=='hybrid' and levelValue=='notgiven':	
#				these=gribs[index](indicatorOfParameter=ipar,typeOfLevel=leveltype,timeRangeIndicator=itr)
#				NLEV=len(these)
#				NLAT=these[0].Ny
#				NLON=these[0].Nx
#				if kz2==-1:
#					kz2=NLEV
#				if jy2==-1:
#					jy2=NLAT
#				if ix2==-1:
#					ix2=NLON
#
#				array=np.zeros((NLEV,NLAT,NLON))
#                		for this in these:
#					array[this.level-1,:]=this.values
#                        	ThisField.append(array[kz1-1:kz2,jy1-1:jy2,ix1-1:ix2])
#			else:
#				this=gribs[index](indicatorOfParameter=ipar,typeOfLevel=leveltype,level=levelValue,timeRangeIndicator=itr)
#				NLAT=this[0].Ny
#				NLON=this[0].Nx
#				if jy2==-1:
#					jy2=NLAT
#				if ix2==-1:
#					ix2=NLON
#				ThisField.append(this[0].values[jy1-1:jy2,ix1-1:ix2] )	
#        return ThisField
#def MakeArray(ipar,itype, gribs, gribMap, kz1,kz2,jy1=0,jy2=-1,ix1=0,ix2=-1,itr=0,DateList=[]):
## Extract one field from a set of grib files; return in 4-d list of lists of arrays (time, level, y-dimension, x-dimension)
## Parameters:
##       ipar:  grib identifyer of requested variable
##       itype: grib identifyer of level type
##       gribs: list of handles on a set of grib files. It is assumed that each grib file contains data from
##               exactly one time level characterized by analysis time and forecast range. It is further 
##               assumed, that a given variable, when present,  is given on identical grids in each file.
##       DateList: if present, list of date info for each file to provide input for the extraction (must be present in gribs) (analysis time, step, valid time)
##       gribMap: list of lists with contents for every file in in gribs, structre:[dates][keylist]
##                date info is given again for each field
##       kz1,kz2:  start and end values of the levels, python style [kz1,kz2[
##       jy1,jy2:  start and end values of the first horizontal index (y), python style [jy1,jy2[
##       ix1,ix2:  start and end values of the second horizontal index (x), python style [ix1,ix2[
##               Note: default setting of jy2,ix2=-1 results in extracting [jy1:,ix1:]
##       itr: grib time range indicator  
#
#	ThisField=[]
#	if DateList==[]:
#        	for grib in gribs:
#                	print grib, datetime.now().time()
#                	array=[]
#                	for level in np.arange(kz1,kz2,1):
#                        	if jy2!=-1 and ix2 !=-1:
#                                	array.append(grib[gribMap.index([ipar,itype,level,itr])].values[jy1:jy2,ix1:ix2])
#                        	elif jy2 == -1:
#                                	array.append(grib[gribMap.index([ipar,itype,level,itr])].values[jy1:,ix1:ix2])
#                        	elif ix2 == -1:
#                                	array.append(grib[gribMap.index([ipar,itype,level,itr])].values[jy1:jy2,ix1:])
#                        	else:
#                                	array.append(grib[gribMap.index([ipar,itype,level,itr])].values[jy1:,ix1:])
#                	ThisField.append(array)
#	else:
#        	for date in DateList:
#                	print date, datetime.now().time()
#                	index=DateList.index(date)
#                	array=[]
#                	for level in np.arange(kz1,kz2,1):
#                        	if jy2!=-1 and ix2 !=-1:
#                                	array.append(gribs[index][gribMap.index([ipar,itype,level,itr])].values[jy1:jy2,ix1:ix2])
#                        	elif jy2 == -1:
#                                	array.append(gribs[index][gribMap.index([ipar,itype,level,itr])].values[jy1:,ix1:ix2])
#                        	elif ix2 == -1:
#                                	array.append(gribs[index][gribMap.index([ipar,itype,level,itr])].values[jy1:jy2,ix1:])
#                        	else:
#                                	array.append(gribs[index][gribMap.index([ipar,itype,level,itr])].values[jy1:,ix1:])
##                	for level in np.arange(kz1,kz2,1):
##                        	if jy2!=-1 and ix2 !=-1:
##                                	array.append(gribs[index][gribMap[index].index([ipar,itype,level,itr,date[0],date[1],date[2]])].values[jy1:jy2,ix1:ix2])
##                        	elif jy2 == -1:
##                                	array.append(gribs[index][gribMap[index].index([ipar,itype,level,itr,date[0],date[1],date[2]])].values[jy1:,ix1:ix2])
##                        	elif ix2 == -1:
##                               	array.append(gribs[index][gribMap[index].index([ipar,itype,level,itr,date[0],date[1],date[2]])].values[jy1:jy2,ix1:])
##                        	else:
##                                	array.append(gribs[index][gribMap[index].index([ipar,itype,level,itr,date[0],date[1],date[2]])].values[jy1:,ix1:])
#
#                	ThisField.append(array)
#        return ThisField
def gribmap(grbs):
        keys=['padding']
        #grbs = pygrib.open(gribfile)
        for grb in grbs:
                if grb.typeOfLevel == 'meanSea':
                        ityp=102
                elif grb.typeOfLevel == 'heightAboveSea':
                        ityp=102
                elif grb.typeOfLevel == 'heightAboveGround':
                        ityp=105
                elif grb.typeOfLevel == 'hybrid':
                        ityp=109
                else:
                        ityp=grb.typeOfLevel
                keys.append( [grb.indicatorOfParameter,ityp,grb.level] )
        return keys


def destag(u,v,Arakawa='C'):
	
	# u,v at mass points:
	if Arakawa=='C':
		mu=u.copy(); mv=v.copy()
		mu[:,0:-1]=(u[:,0:-1]+u[:,1:])/2.
		mv[0:-1,:]=(v[0:-1,:]+v[1:,:])/2. 
	elif Arakawa=='A':
		mu=u; mv=v
	else:
		raise ValueError("Arakawa= "+Arakawa+" is not implemented, options are C and A")
	return mu, mv

def rotwind2regwind(urot,vrot,reglon,reglat,polon,polat,distinctLongitudes,distinctLatitudes,Arakawa='A'):
	# this is a script fortransformattion of staggered u and v 
	# relative to rotated latlons to regular u and v at mass points
        # Based on the hirlam routines regrot, 
	# For Python April 2015, Carl Fortelius
	
	#*&0 : matrice pa
	#*&1 : matrice pb
	#*&2 : matrice pc
	#*&3 : matrice pd
	#*&4:  regular longitude of rotated south pole
	#*&5:  regular latitude of rotated south pole
	#*&6:  regular longitudes (precomputed)
	#*&7:  regular latitudes (precomputed)


	
	rd=np.pi/180.
	zsyc=np.sin(rd*(polat+90.))
	zcyc=np.cos(rd*(polat+90.))
	
	zsxreg=np.sin(rd*reglon)
	zcxreg=np.cos(rd*reglon)
	zsyreg=np.sin(rd*reglat)
	zcyreg=np.cos(rd*reglat)
	
	zsxmxc=np.sin(rd*(reglon-polon))
	zcxmxc=np.cos(rd*(reglon-polon))
	
	zsxrot=np.sin(rd*distinctLongitudes).T * np.ones(reglon.shape) 
	zcxrot=np.cos(rd*distinctLongitudes).T * np.ones(reglon.shape)
	zsyrot=(np.sin(rd*distinctLatitudes).T * np.ones(reglat.shape).T).T 
	zcyrot=(np.cos(rd*distinctLatitudes).T * np.ones(reglat.shape).T).T 
	
	
	zpa = zcyc*zsxmxc*zsxrot + zcxmxc*zcxrot
	zpb = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot - zcxmxc*zsxrot*zsyrot
	zpc = -zsyc*zsxrot/zcyreg
	zpd = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg
	#* this is a script for rotating hirlam winds according
	#* to predefined transformatin matrices prepared by
	#* the scripr regmats or other.
	#
	#
	#*&0 : rotated u (input)
	#*&1 : rotated v (input)
	#*&2 : regular u (output)
	#*&3 : regular v (output)
	#*&4 : matrice pa
	#*&5 : matrice pb
	#*&6 : matrice pc
	#*&7 : matrice pd
	
	
	#regular u,v at mass points:
	if Arakawa=='C':
		murot=urot.copy(); mvrot=vrot.copy()
		murot[:,0:-1]=(urot[:,0:-1]+urot[:,1:])/2.
		mvrot[0:-1,:]=(vrot[0:-1,:]+vrot[1:,:])/2. 
	elif Arakawa=='A':
		murot=urot; mvrot=vrot
	else:
		raise ValueError("Arakawa= "+Arakawa+" is not implemented, options are C and A")
	return zpa*murot + zpb*mvrot, zpc*murot + zpd*mvrot



