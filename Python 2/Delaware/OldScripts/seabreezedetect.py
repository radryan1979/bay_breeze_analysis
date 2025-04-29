from siphon.radarserver import RadarServer
from datetime import datetime, timedelta
from siphon.cdmr import Dataset
from metpy.plots import ctables  # For NWS colortable
from scipy import spatial
from scipy.interpolate import griddata
from dateutil import tz
from time import mktime
import numpy as np
import warnings
import matplotlib as mpl
warnings.filterwarnings("ignore", category=mpl.cbook.MatplotlibDeprecationWarning)
import matplotlib.pyplot as plt
import pyart
import numpy.ma as ma
import time
import pyproj
import cartopy.feature as cfeature
import cartopy.crs as ccrs

###########################################################################
###########################################################################
#                           FUNCTIONS
###########################################################################
###########################################################################

#Function to trim radar data to desired bounding rectangle:
#Required are 1D arrays of latitudes, longitudes, and reflectivities
#and two tuples of latitudes and longitudes defining the boundary of 
#the desired boundary.
def trim_rad_data(lats, lons, ref, boundinglat=0, boundinglon=0):
    for i in range(lats.shape[0]):
        for j in range(lats.shape[1]):
            if (lats[i][j]>boundinglat[0] and lats[i][j]<boundinglat[1] and \
                lons[i][j]>boundinglon[0] and lons[i][j]<boundinglon[1]):
                pass
            else:
                # lats[i] = np.nan; lons[i] = np.nan
                ref[i][j] = np.nan
                
            if (ref[i][j]>loref and ref[i][j]<hiref):
                pass
            else:
                ref[i][j] = np.nan
    return ref
    
#Find index of element in the array nearest to input value
#This will be helpful for finding latitudes to tranect.
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
    
#Regrid data to cartesian grid with polar structure, in other words
#i and j of 2d array are scan angle and range, respectively.
def regrid_to_cartesian(radar):
    display = pyart.graph.RadarMapDisplay(radar)
    x,y = display._get_x_y(0,True,None)
    x = x*1000; y = y*1000
    lambert_aea = {'proj': 'laea',
              'lat_0':lat0, 
              'lon_0':lon0, 
              'x_0':0., 
              'y_0':0.,
              'ellps': 'WGS84',
              'datum': 'WGS84',
              'R':6378137.0}
    lons, lats = pyart.core.cartesian_to_geographic(x,y,
        projparams = lambert_aea)
        
    #Otherwise, lats and lons have 1 extra data point.
    lats = lats[0:720, 0:1832]
    lons = lons[0:720, 0:1832]
    return lons,lats

#Convert polar structure grid to rectilinear
def grid_to_rectilinear(my_ref,lons,lats,nlon,nlat):
    rav_lats = lats.ravel()
    rav_lons = lons.ravel()
    rav_ref = my_ref.ravel()
    
    #Grid Data using matplotlib
    grid_lons = np.linspace(boundinglon[0],boundinglon[1],nlon)
    grid_lats = np.linspace(boundinglat[0],boundinglat[1],nlat)
    glon,glat = np.meshgrid(grid_lons,grid_lats)
    
    # Interpolate data onto grid using linear interpolation
    gref = griddata((rav_lons,rav_lats),rav_ref,(glon,glat),method='linear')
    
    return gref, glon, glat

#Calculate sea breeze front coordinates, raise warning if 
#one is not found.
def sea_breeze_front_coords(gref, glon, glat):
    #For easier processing:
    onedimlat = glat[:,0]
    onedimlon = glon[0,:]
    
    #Initialize sea breeze indexer
    sbfcnt=0
    #Initialize coordinate arrays for sea breeze front
    sbflat=[]
    sbflon=[]
    
    #Loop through transect lines
    for i in range(np.size(cstlnlat)):
        tlat = cstlnlat[i]
        tlon0 = cstlnlon[i]+0.01 #Just a bit offshore in case
        
        #Find index of nearest latitude in grid: this is our i value
        idxlat = find_nearest(onedimlat,tlat)#represents row containing that latitude
        idxlon = find_nearest(onedimlon,tlon0)#represents starting column
        
        #Initialize maxsumref for keeping track of where clustering occurs
        maxsumref = 0.0
        
        #Loop through previous j values - essentially moving from west to east
        for j in range(idxlon,idxlon-300,-1):
            #Test reflectivities surrounding this point to eliminate data
            sumref = np.count_nonzero(~np.isnan(gref[idxlat-5:idxlat+5,j-5:j+5]))
            
            #If larger cluster found, replace previous.
            if sumref>=60:
                #Maybe do some more conditionals based on the distance from a previously
                #found point (see matlab version) ----for now this is okay.
                if sumref>maxsumref:
                    maxsumref = sumref
                    maxrefj = j
            
        #If clustering sufficient for SB was found, place SB point
        if maxsumref>0:
            #Saving coordinates
            sbflat.append(glat[idxlat,maxrefj])
            sbflon.append(glon[idxlat,maxrefj])
            sbfcnt+=1
        
        else:
            sbflat.append(glat[idxlat,0])
            sbflon.append(np.nan)
    
    #If fewer than 34% of the transects (10 out of 30), 
    #then we call this an incomprehensible SB or not found.
    if sbfcnt/len(cstlnlat)<0.34:
        sbfound = True
        #Create tuple of coordinates    
        sbf = list(zip(sbflon,sbflat))
    else:
        sbfound = False
        sbf = ()
        
    
    return sbfound, sbf
    
def sb_day_analysis(rs, date):
    
    #Unit conversion for velocity and penetration:
    #latitude = 38 used for this conversion.
    #Source: http://www.csgnetwork.com/degreelenllavcalc.html
    met_per_degree = 87832.43
    
    #For analyzing until 23:59UTC ~sundown.
    nextday = date + timedelta(days=1)
    
    #Usually will be at hours=13 to start around sunrise.
    dt = date + timedelta(hours=13)
    
    query=rs.query()
    query.stations('KDOX').time_range(dt,nextday)#dt+timedelta(minutes=15))
    cat = rs.get_catalog(query)
    ds = list(cat.datasets.values())
    
    #Number of sea breeze front successfully identified. Resets every
    #day/every call of this function.
    num_sb_found=0
    
    for numrad in range(len(ds)):
        
        #Create radar object
        radar = pyart.io.read_nexrad_cdm(ds[numrad].access_urls['OPENDAP'])
        
        #Calculate timestamp to decrease computing time.
        timestamp = radar.time['units'].split(' ')[-1].split('T')
        timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
        timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
        
        prevtime = datetime(2018,12,31)
        #If next iteration is within 15 mins of previous iteration, this scan is skipped.
        timediff = (timestamp-prevtime).seconds/60.0/24.0
        if numrad>0 and timediff<0.25:
            continue
        prevtime = timestamp
        
        #Establish location of radar
        loc = pyart.io.nexrad_common.get_nexrad_location('KDOX')
        lon0 = loc[1] ; lat0 = loc[0]
        
        
        
        
    ###########################################################################
    ###########################################################################
    #                               CREAT TIMESTAMP
    ###########################################################################
    ###########################################################################
        
        
        #Create timestamp, if needed.
        # timestamp = radar.time['units'].split(' ')[-1].split('T')
        # timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
        # timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
        # from_zone = tz.gettz('UTC')
        # to_zone = tz.gettz('America/New_York')
        # utc = timestamp.replace(tzinfo=from_zone)
        # local = utc.astimezone(to_zone)
        # lt = time.localtime()
        # dst = lt.tm_isdst
        # lt = time.localtime()
        # dst = lt.tm_isdst
        # if dst == 0:
        #     et = "EDT"
        # else:
        #     et = "EST"
            
            
        
    ###########################################################################
    ###########################################################################
    #                   RE-GRIDDING TO LAT/LON - POLAR GRID
    ###########################################################################
    ###########################################################################
            
        lons, lats = regrid_to_cartesian(radar)
        
    ###########################################################################
    ###########################################################################
    #                           HANDLING REF DATA
    ###########################################################################
    ###########################################################################
        
        
        #Organizing reflectivity data
        ref = radar.get_field(0, 'reflectivity')
        unmasked = ma.getdata(ref)  
        #Following should not be needed, used further down.
        # my_lats = lats.ravel()
        # my_lons = lons.ravel()
        # my_ref = my_ref.ravel()
        
        #Utilizing function at beginning to trim the ref data ----not sure if needed.
        my_ref = trim_rad_data(lats, lons, unmasked, boundinglat, boundinglon)
        
    ###########################################################################
    ###########################################################################
    #                       RE-GRIDDING DATA TO RECTILINEAR GRID
    ###########################################################################
    ###########################################################################
        
        nlon = 1000; nlat = 1000; #can be changed, but seems good.
        gref, glon, glat = grid_to_rectilinear(my_ref,lons,lats,nlon,nlat)
        
    ###########################################################################
    ###########################################################################
    #                       ANALYZE INDIVIDUAL SCANS
    ###########################################################################
    ###########################################################################
    
        sbfound, sbf = sea_breeze_front_coords(gref,glon,glat)
        
        if sbfound:
            
            #Creating list of lats and lons for each timestep
            #These will then be used to calculate inland pen, velocity, etc.
            #to throw back at outter loop or calling script.
            sbflon[numrad,:], sbflat[numrad,:] = zip(*sbf)
            
            sbflon_n[numrad] = np.nanmean(sbflon[numrad,0:6])
            sbflon_nc[numrad] = np.nanmean(sbflon[numrad,6:12])
            sbflon_c[numrad] = np.nanmean(sbflon[numrad,12:18])
            sbflon_sc[numrad] = np.nanmean(sbflon[numrad,18:24])
            sbflon_s[numrad] = np.nanmean(sbflon[numrad,24:30])
            
            if num_sb_found>=1:
                #Calculate penetration and velocity - two (more?) methods.
                #1) Track 'sections' (for now, 5) representing different
                #   branches of the sea breeze-north, north-central, central,
                #   south-central, and south.
                
                #Calculate velocity, index-1 because we are finite differencing.
                #in units [m/hr], negative is westward
                sb_vel_n[num_sb_found-1] = (sbflon_n[numrad]-sbflon_n[idx_prev_sb_found])\
                                            /timediff*met_per_degree
                sb_vel_nc[num_sb_found-1] = (sbflon_nc[numrad]-sbflon_nc[idx_prev_sb_found])\
                                            /timediff*met_per_degree
                sb_vel_c[num_sb_found-1] = (sbflon_c[numrad]-sbflon_c[idx_prev_sb_found])\
                                            /timediff*met_per_degree
                sb_vel_sc[num_sb_found-1] = (sbflon_sc[numrad]-sbflon_sc[idx_prev_sb_found])\
                                            /timediff*met_per_degree
                sb_vel_s[num_sb_found-1] = (sbflon_s[numrad]-sbflon_s[idx_prev_sb_found])\
                                            /timediff*met_per_degree
                
                #Calculate inland penetration in meters.
                sb_pen_n[num_sb_found-1] = (sbflon_n[numrad]-np.mean(cstlnlon[0:6]))*met_per_degree
                sb_pen_nc[num_sb_found-1] = (sbflon_nc[numrad]-np.mean(cstlnlon[0:6]))*met_per_degree
                sb_pen_c[num_sb_found-1] = (sbflon_c[numrad]-np.mean(cstlnlon[0:6]))*met_per_degree
                sb_pen_sc[num_sb_found-1] = (sbflon_sc[numrad]-np.mean(cstlnlon[0:6]))*met_per_degree
                sb_pen_s[num_sb_found-1] = (sbflon_s[numrad]-np.mean(cstlnlon[0:6]))*met_per_degree
            
            num_sb_found+=1
            idx_prev_sb_found = numrad
        else:
            sbflon[numrad,:] = np.nan
            sbflat[numrad,:] = np.nan
            
            #Portray that no sb was found on this timestep.
            sbflon_n[numrad] = np.nan
            sbflon_nc[numrad] = np.nan
            sbflon_c[numrad] = np.nan
            sbflon_sc[numrad] = np.nan
            sbflon_s[numrad] = np.nan

    ###########################################################################
    ###########################################################################
    #                           ANALYZE DAILY DATA
    ###########################################################################
    ###########################################################################
    
    #Eliminate days with insufficient data (total of ~45 scans per day
    #analyzed by this program).
    if numrad<20:
        day_sb_found = False
        error = "Insufficient data: Not enough scans to analyze this day."
        sb_vel_n[:] = np.nan
        sb_pen_n[:] = np.nan
        sb_vel_nc[:] = np.nan
        sb_pen_nc[:] = np.nan
        sb_vel_c[:] = np.nan
        sb_pen_c[:] = np.nan
        sb_vel_sc[:] = np.nan
        sb_pen_sc[:] = np.nan
        sb_vel_s[:] = np.nan
        sb_pen_s[:] = np.nan
        
        sb_vel = list(zip(sb_vel_n, sb_vel_nc, sb_vel_c, sb_vel_sc, sb_vel_s))
        sb_pen = list(zip(sb_pen_n, sb_pen_nc, sb_pen_c, sb_pen_sc, sb_pen_s))
        return day_sb_found, error, sb_vel, sb_pen
    
    elif num_sb_found<10:
        day_sb_found = False
        error = "Sufficient data, but sea breeze not found."
        sb_vel_n[:] = np.nan
        sb_pen_n[:] = np.nan
        sb_vel_nc[:] = np.nan
        sb_pen_nc[:] = np.nan
        sb_vel_c[:] = np.nan
        sb_pen_c[:] = np.nan
        sb_vel_sc[:] = np.nan
        sb_pen_sc[:] = np.nan
        sb_vel_s[:] = np.nan
        sb_pen_s[:] = np.nan
        
        sb_vel = list(zip(sb_vel_n, sb_vel_nc, sb_vel_c, sb_vel_sc, sb_vel_s))
        sb_pen = list(zip(sb_pen_n, sb_pen_nc, sb_pen_c, sb_pen_sc, sb_pen_s))
        return day_sb_found, error, sb_vel, sb_pen



    ###########################################################################
    ###########################################################################
    #                           Return Function
    ###########################################################################
    ###########################################################################
    
    sb_vel = list(zip(sb_vel_n, sb_vel_nc, sb_vel_c, sb_vel_sc, sb_vel_s))
    sb_pen = list(zip(sb_pen_n, sb_pen_nc, sb_pen_c, sb_pen_sc, sb_pen_s))
    error = "Good day, sir."
    return day_sb_found, error, sb_vel, sb_pen
