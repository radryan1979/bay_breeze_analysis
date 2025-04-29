"""
Coded by Dan Moore
Test day(s) for sea breeze. Returns velocity along 5 latitudinal sections
(problems to be diagnosed).
As of now, the program willl successfully iterate through radar data by day,
in 15 min increments.

"""
from siphon.radarserver import RadarServer
from datetime import datetime,timedelta
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
import pandas as pd


###########################################################################
###########################################################################
#                           FUNCTIONS
###########################################################################
###########################################################################

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
    query.stations('KDOX').time_range(dt,nextday)
    cat = rs.get_catalog(query)
    ds = list(cat.datasets.values())
    
    #Number of sea breeze front successfully identified. Resets every
    #day/every call of this function.
    num_sb_found=0
    
    #Maintain scan times that are actually read by the program.
    scan_time=[]
    
    #Initializing areal location arrays:
    # sbflon_n=[]
    # sbflon_nc=[]
    # sbflon_c=[]
    # sbflon_sc=[]
    # sbflon_s=[]
    sbf_lon_day = []
    
    #Velocity arrays
    # sb_vel_n=[]
    # sb_vel_nc=[]
    # sb_vel_c=[]
    # sb_vel_sc=[]
    # sb_vel_s=[]
    sb_vel = []

    #Penetration arrays
    # sb_pen_n=[]
    # sb_pen_nc=[]
    # sb_pen_c=[]
    # sb_pen_sc=[]
    # sb_pen_s=[]
    sb_pen = []
    
    #If no scans are found for the given date:
    if len(ds)==0:
        day_sb_found = False
        error = "No data found for selected date."
        
        scan_time = np.nan
        sb_vel = np.nan
        sb_pen = np.nan
        
        return day_sb_found, error, sbf_lon_day, scan_time
        
    for numrad in range(len(ds)):
        
        try: 
            #Create radar object
            radar = pyart.io.read_nexrad_cdm(ds[numrad].access_urls['OPENDAP'])
        except:
            continue
        
        #Calculate timestamp to decrease computing time.
        timestamp = radar.time['units'].split(' ')[-1].split('T')
        timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
        timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
        
        if numrad>0:
            timediff = (timestamp-prevtime).seconds/60.0/60.0
        
        #If next iteration is within 15 mins of previous iteration, this scan is skipped.
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
            
        lons, lats = regrid_to_cartesian(radar, lon0, lat0)
        
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
        
        nlon = 1000; nlat = 1000 #can be changed, but seems good.
        gref, glon, glat = grid_to_rectilinear(my_ref,lons,lats,nlon,nlat)
        
    ###########################################################################
    ###########################################################################
    #                       ANALYZE INDIVIDUAL SCANS
    ###########################################################################
    ###########################################################################
    
        sbfound, sbf = sea_breeze_front_coords(gref,glon,glat)
        
        #Maintain scan time to know 
        scan_time.append(timestamp)
        
        if sbfound:
            
            #Creating list of lats and lons for each timestep
            #These will then be used to calculate inland pen, velocity, etc.
            #to throw back at outter loop or calling script.
            sbflon, sbflat = zip(*sbf)

            # sbflon_n.append(np.nanmean(sbflon[0:6]))
            # sbflon_nc.append(np.nanmean(sbflon[6:12]))
            # sbflon_c.append(np.nanmean(sbflon[12:18]))
            # sbflon_sc.append(np.nanmean(sbflon[18:24]))
            # sbflon_s.append(np.nanmean(sbflon[24:30]))
            sbf_lon_day.append(list(sbflon))
            
            if num_sb_found>=1:
                pass
    ###########################################################################
    ###########################################################################
    #                       VELOCITY & PENETRATION INFO
    ###########################################################################
    ###########################################################################
                
                
                #Calculate penetration and velocity - two (more?) methods.
                #1) Track 'sections' (for now, 5) representing different
                #   branches of the sea breeze-north, north-central, central,
                #   south-central, and south.
                
                #Calculate velocity, index-1 because we are finite differencing.
                #in units [m/hr], negative is westward
                # sb_vel_n.append((sbflon_n[num_sb_found]-sbflon_n[num_sb_found-1])\
                #                             /timediff*met_per_degree)
                # sb_vel_nc.append((sbflon_nc[num_sb_found]-sbflon_nc[num_sb_found-1])\
                #                             /timediff*met_per_degree)
                # sb_vel_c.append((sbflon_c[num_sb_found]-sbflon_c[num_sb_found-1])\
                #                             /timediff*met_per_degree)
                # sb_vel_sc.append((sbflon_sc[num_sb_found]-sbflon_sc[num_sb_found-1])\
                #                             /timediff*met_per_degree)
                # sb_vel_s.append((sbflon_s[num_sb_found]-sbflon_s[num_sb_found-1])\
                #                             /timediff*met_per_degree)
                
                #We will do this outside of this function.
                # sb_vel.append((sbf_lon_day[:,num_sb_found] - sbf_lon_day[:,num_sb_found-1])\
                #                             /timediff * met_per_degree)
                
                #Calculate inland penetration in meters.
                # sb_pen_n.append((sbflon_n[num_sb_found]-np.mean(cstlnlon[0:6]))*met_per_degree)
                # sb_pen_nc.append((sbflon_nc[num_sb_found]-np.mean(cstlnlon[0:6]))*met_per_degree)
                # sb_pen_c.append((sbflon_c[num_sb_found]-np.mean(cstlnlon[0:6]))*met_per_degree)
                # sb_pen_sc.append((sbflon_sc[num_sb_found]-np.mean(cstlnlon[0:6]))*met_per_degree)
                # sb_pen_s.append((sbflon_s[num_sb_found]-np.mean(cstlnlon[0:6]))*met_per_degree)
                
                #Not sure if there is a simpler way to do this.
                # sb_pen_temp = []
                # for j in range(len(cstlnlon)):
                #     sb_pen_temp[i] = (sbf_lon_day[i,num_sb_found]-cstlnlon[i]) * met_per_degree
                # sb_pen.append(sb_pen_temp)
            
            num_sb_found+=1
        else:

            #Portray that no sb was found on this timestep.
            # sbflon_n.append(np.nan)
            # sbflon_nc.append(np.nan)
            # sbflon_c.append(np.nan)
            # sbflon_sc.append(np.nan)
            # sbflon_s.append(np.nan)
            
            #If none found, append an nan for each transect.
            sbf_lon_day.append([np.nan]*len(cstlnlon))

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
        
        # sb_vel = np.nan
        # sb_pen = np.nan
        return day_sb_found, error, sbf_lon_day, scan_time
    
    elif num_sb_found<10:
        day_sb_found = False
        error = "Sufficient data, but sea breeze not found."
        
        sb_vel = np.nan
        sb_pen = np.nan
        return day_sb_found, error, sbf_lon_day, scan_time

    ###########################################################################
    ###########################################################################
    #                           Return Function
    ###########################################################################
    ###########################################################################
    
    # sb_vel = list(zip(sb_vel_n, sb_vel_nc, sb_vel_c, sb_vel_sc, sb_vel_s))
    # sb_pen = list(zip(sb_pen_n, sb_pen_nc, sb_pen_c, sb_pen_sc, sb_pen_s))
    day_sb_found = True
    error = "Good day, sir."
    return day_sb_found, error, sbf_lon_day, scan_time


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
def regrid_to_cartesian(radar, lon0, lat0):
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
    if sbfcnt/len(cstlnlat)>0.34:
        sbfound = True
        #Create tuple of coordinates    
        sbf = list(zip(sbflon,sbflat))
    else:
        sbfound = False
        sbf = ()
        
    
    return sbfound, sbf
    






###########################################################################
###########################################################################
#                           INITIALIZATIONS
###########################################################################
###########################################################################


#Boundaries for desired values
boundinglat = (38.15, 39.65)
boundinglon = (-76.45, -74.6)

#Defining the DE Coastline
cstlnlat = np.array((39.175, 39.15, 39.125, 39.1, 39.075, 39.05, 39.025, 39.0, #Bay Coastline
                    38.975, 38.95, 38.925, 38.9, 38.875, 38.85, 38.825, 
                    38.8, 38.775,
                    38.75, 38.725, 38.7, 38.675, 38.65, 38.625, 38.6,#Ocean Coastline
                    38.575, 38.55, 38.525, 38.5, 38.475, 38.45)).T
                    
cstlnlon = np.array((-75.41, -75.41, -75.41, -75.40, -75.40, -75.39, -75.35, -75.33,#Bay Coastline
                    -75.32, -75.31, -75.31, -75.29, -75.26, -75.24, -75.21,
                    -75.18, -75.08,
                    -75.08, -75.08, -75.07, -75.07, -75.07, -75.06, -75.06,#Ocean Coastline
                    -75.06, -75.06, -75.05, -75.05, -75.05, -75.05)).T

#Reflectivity thresholds
loref = 8; hiref = 25 #dBz thresholds

###########################################################################
###########################################################################
#                               READ IN DATA
###########################################################################
###########################################################################

#Read in data from unidata server - option to loop through.
#Eventually write this is a function to be called every day of interest 
#to analyze for sea breeze and spit out information.

rs=RadarServer('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')

sb_found = []
error_log = []

#Enter conditional below for initializing pandas dataframe
first_found = False

#Write for loop to go through range of dates over time period:

for i in range(1,32):
    
    date = datetime(2015, 8, i)
    
    #Important variables returned. Latitudes of transects important for analysis
    sb_found_bool, error, sbf_lon, scan_time = sb_day_analysis(rs, date)
    
    sb_found.append(sb_found_bool)
    error_log.append(error)
    
    #Only enter this the first sb day found to initialize dataframe.
    if sb_found_bool and not first_found:
        #Prevents entering this again.
        first_found = True
        
        #Creating dataframe for overall analysis
        df = pd.DataFrame(sbf_lon)
        df = df.transpose()
        df.columns = scan_time
        df.index = cstlnlat
        
        print("Sea breeze detected on ", date)
        
        print(datetime.now())
        
        print(error)
        
        continue
    
    print(error)
    
    if sb_found_bool:
        
        #Can do analysis by latitudinal band.
        # for lat in transect_lats:
        
        #Or by time.
        # for i in range(len(scan_time)):
        
        #Maintain temp dataframe to add to accumulating data frame for time-range.
        df_temp = pd.DataFrame(sbf_lon)
        df_temp = df_temp.transpose()
        df_temp.columns = scan_time
        df_temp.index = cstlnlat
        
        df = pd.concat([df,df_temp], axis=1)
        
        print("Sea breeze detected on ", date)
        # print("Coordinates are ", sbf_lon)
        # print("On "+str(datelist[i]+"sb_found is "+sb_found[i]+"with sb_vel = "+sb_vel+\
        #     "and sb_pen = "+sb_pen))
        
        print(datetime.now())
    
    else:
        print("No sea breeze detected on ", date)
        
        print(datetime.now())

df.to_csv("/Volumes/LaCie/SeaBreeze/RadarDetectOutput/FirstTest.csv")

#Notification that the program has finished
print('Program Finished')

