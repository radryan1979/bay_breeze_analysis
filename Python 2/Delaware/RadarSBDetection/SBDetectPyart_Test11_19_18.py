"""
Coded by Dan Moore
Final program to ingest 5 years worth of radar data and output
coordinates of SBF for analysis.

Updated: 11_14_18
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
from scipy import ndimage
from skimage import measure
from skimage import filters
from skimage.measure import regionprops


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
    
    #Establish location of radar
    loc = pyart.io.nexrad_common.get_nexrad_location('KDOX')
    lon0 = loc[1] ; lat0 = loc[0]
    
    #Number of sea breeze front successfully identified. Resets every
    #day/every call of this function.
    num_sb_found=0

    #Number of scans actually analyzed.
    scan_num = 0
    
    #Maintain scan information for scans that are actually read by the program.
    scan_time=[]
    precip_found=[]
    sea_breeze_found=[]
    radar_mode=[]
    
    #Initializing areal location arrays:
    sbf_lon_day = []
    
    #If no scans are found for the given date:
    if len(ds)==0:
        day_sb_found = False
        error = "No data found for selected date."
        
        scans_df=pd.DataFrame()
        
        return day_sb_found, error, sbf_lon_day, scans_df
    #Eliminate days with insufficient data (total of ~45 scans per day
    #analyzed by this program).
    elif len(ds)<20:
        day_sb_found = False
        error = "Insufficient data: Not enough scans to analyze this day."
        
        scans_df=pd.DataFrame()

        return day_sb_found, error, sbf_lon_day, scans_df
        
    for numrad in range(len(ds)):
        
        try: 
            #Create radar object
            radar = pyart.io.read_nexrad_cdm(ds[numrad].access_urls['OPENDAP'])
        except:
            print("Line 84: read_nexrad_cdm error. Please check NEXRAD server if problem persists.")
            continue
        
        #Calculate timestamp to decrease computing time.
        timestamp = radar.time['units'].split(' ')[-1].split('T')
        timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
        timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
        
        #Determine mode of radar (nsweeps>10, precip mode)
        if radar.nsweeps > 10:
            rad_mode = 'precip'
        else:
            rad_mode = 'clear'
        
        if scan_num>0:
            timediff = (timestamp-prevtime).seconds/60.0/60.0
        
        #If next iteration is within 15 mins of previous iteration, this scan is skipped.
        if scan_num>0 and timediff<0.25:
            continue
            
            
        
    ###########################################################################
    ###########################################################################
    #                   ASSIGNING LAT/LON - POLAR GRID
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
        
        #Check if there is substantial precipitation. If amount of reflectivity
        #values above 25 dBZ is more than 15% of all nonzero values in the 
        #dataset, then skip the day.
        if np.count_nonzero(unmasked[unmasked>25.0])/\
        np.count_nonzero(unmasked) >= 0.15:
            sbf_lon_day.append([np.nan]*len(cstlnlon))
            prevtime = timestamp
            scan_num+=1
            continue
        
        if lats.shape != ref.shape:
            print("Error, data shape mismatch, scan skipped.")
            continue
        
        #If scan passes all the above if statements, we will evaluate it.
        prevtime = timestamp
        scan_num+=1
        
        #Utilizing function at beginning to trim the ref data.
        my_ref = trim_rad_data(lats, lons, unmasked, boundinglat, boundinglon)
        
    ###########################################################################
    ###########################################################################
    #                       RE-GRIDDING DATA TO RECTILINEAR GRID
    ###########################################################################
    ###########################################################################
        
        nlon = 1000; nlat = 1000 #can be changed, but seems good.
        gref, glon, glat = grid_to_rectilinear(unmasked,lons,lats,nlon,nlat)
        
        #Trim data further, and detect precip.
        gref, sb_precip = connected_areas(gref)
        
    ###########################################################################
    ###########################################################################
    #                       ANALYZE INDIVIDUAL SCANS
    ###########################################################################
    ###########################################################################
        if num_sb_found == 0:
            sbfound, sbf = first_pass_sea_breeze_front_coords(gref,glon,glat)
        else:
            #If SBF has been found previously on this day, use previous coordinates
            #to minimize error of erroneous SBF coordinates.
            sbfound, sbf = sea_breeze_front_coords(gref,glon,glat,sbflon)
        
        #Maintain scan data:
        #Maintain scan time to know 
        scan_time.append(timestamp)
        precip_found.append(sb_precip)
        sea_breeze_found.append(sbfound)
        radar_mode.append(rad_mode)
        
        if sbfound:
            
            #Creating list of lats and lons for each timestep
            #These will then be used to calculate inland pen, velocity, etc.
            #to throw back at outter loop or calling script.
            sbflon, sbflat = zip(*sbf)

            sbf_lon_day.append(list(sbflon))
            
            num_sb_found+=1
        else:
            
            #If none found, append an nan for each transect.
            sbf_lon_day.append([np.nan]*len(cstlnlon))

    ###########################################################################
    ###########################################################################
    #                           ANALYZE DAILY DATA
    ###########################################################################
    ###########################################################################
    
    if num_sb_found<8:
        day_sb_found = False
        error = "Sufficient data, but sea breeze not found."
        
        scans_df = pd.DataFrame({'Scan_Time':scan_time,
                                'Precipitation':precip_found,
                                'Sea_Breeze':sea_breeze_found,
                                'Radar_Mode':radar_mode
                                })
        scans_df.set_index('Scan_Time')
        
        return day_sb_found, error, sbf_lon_day, scans_df

    ###########################################################################
    ###########################################################################
    #                           Return Function
    ###########################################################################
    ###########################################################################

    day_sb_found = True
    error = "Good day, sir."
    
    scans_df = pd.DataFrame({'Scan_Time':scan_time,
                            'Precipitation':precip_found,
                            'Sea_Breeze':sea_breeze_found,
                            'Radar_Mode':radar_mode
                            })
    scans_df.set_index('Scan_Time')
        
    return day_sb_found, error, sbf_lon_day, scans_df


#Function to trim radar data to desired bounding rectangle:
#Required are 1D arrays of latitudes, longitudes, and reflectivities
#and two tuples of latitudes and longitudes defining the boundary of 
#the desired boundary.
def trim_rad_data(lats, lons, ref, boundinglat, boundinglon):
    for i in range(lats.shape[0]):
        for j in range(lats.shape[1]):
            if (lats[i][j]>boundinglat[0] and lats[i][j]<boundinglat[1] and \
                lons[i][j]>boundinglon[0] and lons[i][j]<boundinglon[1]):
                pass
            else:
                # lats[i] = np.nan; lons[i] = np.nan
                ref[i][j] = np.nan
                
            if ref[i][j] >= loref:# and ref[i][j]<=hiref):
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
    display = pyart.graph.RadarMapDisplayCartopy(radar)
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
def first_pass_sea_breeze_front_coords(gref, glon, glat):
    #For easier processing:
    onedimlat = glat[:,0]
    onedimlon = glon[0,:]
    
    #Initialize sea breeze indexer
    sbfcnt=0
    #Initialize coordinate arrays for sea breeze front
    sbflat=[]
    sbflon=[]
    
    #Index list for knowing where a SBF was detected
    sbf_idx_list=[]
    
    #Set end of detection
    idxlonend = find_nearest(onedimlon,-76.0)
    
    #First Pass Loop through transect lines
    for i in range(np.size(cstlnlat)):
        tlat = cstlnlat[i]
        tlon0 = cstlnlon[i]+0.01 #Just a bit offshore in case
        
        #Find index of nearest latitude in grid: this is our i value
        idxlat = find_nearest(onedimlat,tlat)#represents row containing that latitude
        idxlon = find_nearest(onedimlon,tlon0)#represents starting column
        
        #Initialize maxsumref for keeping track of where clustering occurs
        maxsumref = 0.0
        
        #Initialize maxsumref for keeping track of where clustering occurs
        maxsumnan = 0.0
        
        #Loop through previous j values - essentially moving from west to east
        for j in range(idxlon,idxlonend,-1):
            #Test reflectivities surrounding this point to eliminate data
            sumref = np.count_nonzero(~np.isnan(gref[idxlat-5:idxlat+5,j-5:j+5]))
            
            #If larger cluster found, replace previous.
            if sumref>=60:
                
                #Added conditional to be within 0.3degrees of coastline - empirical.
                if sumref>maxsumref:# and \
                    # (onedimlon[idxlon]-onedimlon[j])<0.3:
                        maxsumref = sumref
                        maxrefj = j
        
        #If clustering sufficient for SB was found, place SB point
        if maxsumref>0:
            #Saving coordinates
            sbflat.append(glat[idxlat,maxrefj])
            sbflon.append(glon[idxlat,maxrefj])
            sbf_idx_list.append(i)
            sbfcnt+=1
        
        else:
            sbflat.append(glat[idxlat,0])
            sbflon.append(np.nan)
                
        #Catch runs smaller than 4 SBF points (out of 8 lines) so that 
        #points are not mistaken as SBF points.

                
        
    
    #If fewer than 34% of the transects (10 out of 30), 
    #then we call this an incomprehensible SB or not found.
    if sbfcnt/len(cstlnlat)>0.25:
        
        #Place sbf coords on transects that did not find sbf.
        # for i in range(len(sbf_idx_list)-1):
        #     points = [(sbflon[sbf_idx_list[i]],sbflat[sbf_idx_list[i]]),
        #             (sbflon[sbf_idx_list[i+1]],sbflat[sbf_idx_list[i+1]])]
        #     x_coords, y_coords = zip(*points)
        #     A = np.vstack([x_coords,np.ones(len(x_coords))]).T
        #     m, b = np.linalg.lstsq(A, y_coords)[0]
            
            #Difference between indices of consecutively found SBF points.
            # idx_diff = sbf_idx_list[i+1] - sbf_idx_list[i]
            
            #If idx_diff==1, then they are consecutive transects, skip.
            #If idx_diff>3, then they are probably two separate disconnected branches.
            # if idx_diff == 1 or idx_diff>3:
            #     continue
            
            # for j in range(sbf_idx_list[i]+1,sbf_idx_list[i+1]):
            #     sbflon[j] = (sbflat[j]-b)/m
        
        sbfound = True
        #Create tuple of coordinates    
        sbf = list(zip(sbflon,sbflat))
    else:
        sbfound = False
        sbf = list(zip(sbflon,sbflat))
        # sbf = ()
    
    return sbfound, sbf
    
#Calculate sea breeze front coordinates, raise warning if 
#one is not found.
def sea_breeze_front_coords(gref, glon, glat, prev_lon):
    #For easier processing:
    onedimlat = glat[:,0]
    onedimlon = glon[0,:]
    
    #Initialize sea breeze indexer
    sbfcnt=0
    #Initialize coordinate arrays for sea breeze front
    sbflat=[]
    sbflon=[]
    
    #Index list for knowing where a SBF was detected
    sbf_idx_list=[]
    
    #Set end of detection
    idxlonend = find_nearest(onedimlon,-76.0)
    
    #First Pass Loop through transect lines
    for i in range(np.size(cstlnlat)):
        tlat = cstlnlat[i]
        tlon0 = cstlnlon[i]+0.01 #Just a bit offshore in case
        
        #Find index of nearest latitude in grid: this is our i value
        idxlat = find_nearest(onedimlat,tlat)#represents row containing that latitude
        idxlon = find_nearest(onedimlon,tlon0)#represents starting column
        
        #Initialize maxsumref for keeping track of where clustering occurs
        maxsumref = 0.0
        
        #Initialize maxsumref for keeping track of where clustering occurs
        maxsumnan = 0.0
        
        #Loop through previous j values - essentially moving from west to east
        for j in range(idxlon,idxlonend,-1):
            #Test reflectivities surrounding this point to eliminate data
            sumref = np.count_nonzero(~np.isnan(gref[idxlat-5:idxlat+5,j-5:j+5]))
            
            #If larger cluster found, replace previous.
            if sumref>=60:
                
                try:
                    #Added conditional to be within 0.3degrees of previous point - empirical.
                    if sumref>maxsumref and \
                        abs(prev_lon[j]-onedimlon[j])<0.25:
                            maxsumref = sumref
                            maxrefj = j
                            
                except:
                    #Added conditional to be within 0.3degrees of previous point - empirical.
                    if sumref>maxsumref: 
                        maxsumref = sumref
                        maxrefj = j                    
        
        #If clustering sufficient for SB was found, place SB point
        if maxsumref>0:
            #Saving coordinates
            sbflat.append(glat[idxlat,maxrefj])
            sbflon.append(glon[idxlat,maxrefj])
            sbf_idx_list.append(i)
            sbfcnt+=1
        
        else:
            sbflat.append(glat[idxlat,0])
            sbflon.append(np.nan)
                
        #Catch runs smaller than 4 SBF points (out of 8 lines) so that 
        #points are not mistaken as SBF points.

                
        
    
    #If fewer than 34% of the transects (10 out of 30), 
    #then we call this an incomprehensible SB or not found.
    if sbfcnt/len(cstlnlat)>0.25:
        
        #Place sbf coords on transects that did not find sbf.
        # for i in range(len(sbf_idx_list)-1):
        #     points = [(sbflon[sbf_idx_list[i]],sbflat[sbf_idx_list[i]]),
        #             (sbflon[sbf_idx_list[i+1]],sbflat[sbf_idx_list[i+1]])]
        #     x_coords, y_coords = zip(*points)
        #     A = np.vstack([x_coords,np.ones(len(x_coords))]).T
        #     m, b = np.linalg.lstsq(A, y_coords)[0]
            
            #Difference between indices of consecutively found SBF points.
            # idx_diff = sbf_idx_list[i+1] - sbf_idx_list[i]
            
            #If idx_diff==1, then they are consecutive transects, skip.
            #If idx_diff>3, then they are probably two separate disconnected branches.
            # if idx_diff == 1 or idx_diff>3:
            #     continue
            
            # for j in range(sbf_idx_list[i]+1,sbf_idx_list[i+1]):
            #     sbflon[j] = (sbflat[j]-b)/m
        
        sbfound = True
        #Create tuple of coordinates    
        sbf = list(zip(sbflon,sbflat))
    else:
        sbfound = False
        sbf = list(zip(sbflon,sbflat))
        # sbf = ()
    
    return sbfound, sbf


def connected_areas(gref2):
    
    gref2[np.isnan(gref2)]=0
    
    #Determine if there is precip in the scan.
    precip=False
    
    labeled, ncomponents = ndimage.measurements.label(gref2)

    means = ndimage.mean(gref2, labeled, range(ncomponents + 1))
    
    areas = [r.filled_area for r in regionprops(labeled)]
    perimeters = [r.perimeter for r in regionprops(labeled)]
    labels = [r.label for r in regionprops(labeled)]
    eccentricity = [r.eccentricity for r in regionprops(labeled)]
    orientation = [r.orientation for r in regionprops(labeled)]
    
    for i in range(len(areas)):
        
        p_a_ratio = float(perimeters[i])/float(areas[i])
        idx = labeled==labels[i]
        
        if areas[i]>100 and means[i]>20.0:
            precip=True
        else:
            pass
        
        if areas[i]>750 and p_a_ratio>0.20 and means[i]<20.0 and \
            (-np.pi/2.<orientation[i]<-np.pi/6. or np.pi/2.>orientation[i]>np.pi/(6.)):
            pass
        else:
            labeled[idx]=0
            gref2[idx]=np.nan
    
    gref2[gref2==0]=np.nan
    
    return gref2, precip



###########################################################################
###########################################################################
#                           INITIALIZATIONS
###########################################################################
###########################################################################


#Boundaries for desired values
boundinglat = (38.15, 39.65)
boundinglon = (-76.45, -74.6)


#Defining the DE Coastline
cstlnlat = np.array((39.25, 39.225, 39.2, 39.175, 39.15, 39.125, 39.1, 
                    39.075, 39.05, 39.025, 39.0, #Bay Coastline
                    38.975, 38.95, 38.925, 38.9, 38.875, 38.85, 38.825, 
                    38.8, 38.775,
                    38.75, 38.725, 38.7, 38.675, 38.65, 38.625, 38.6,#Ocean Coastline
                    38.575, 38.55, 38.525, 38.5, 38.475, 38.45, 38.425,
                    38.4, 38.375, 38.35, 38.325, 38.3, 38.275, 38.25)).T
                    
cstlnlon = np.array((-75.42, -75.42, -75.41, -75.41, -75.41, -75.41, -75.40, 
                    -75.40, -75.39, -75.35, -75.33,#Bay Coastline
                    -75.32, -75.31, -75.31, -75.29, -75.26, -75.24, -75.21,
                    -75.18, -75.08,
                    -75.08, -75.08, -75.07, -75.07, -75.07, -75.06, -75.06,#Ocean Coastline
                    -75.06, -75.06, -75.05, -75.05, -75.05, -75.05, -75.05,
                    -75.05, -75.06, -75.07, -75.08, -75.09, -75.1, -75.11)).T

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

w_dir = "/home/1898/RadSBDetect/Test_11_19_18/"

sb_found = []
error_log = []

cols = ['Scan_Time','Precipitation','Sea_Breeze','Radar_Mode']

scan_log = pd.DataFrame(index=None,columns=cols)

scan_log.set_index('Scan_Time')


#Write for loop to go through range of dates over time period:

for i in range(2015,2016):
    
    print("Processing Year: ", i)
    print("The Current time is: ", datetime.now())
    
    for j in range(8,9):
        
        #Initialize dataframe every new month.
        df = pd.DataFrame()
        
        #Enter conditional below for initializing pandas dataframe
        first_found = False
        
        for k in range(4,10):
            
            #Make sure it skips the 31st day if it is June or September.
            if (j==6 or j==9) and k==31:
                continue
    
            date = datetime(i, j, k)
            
            #Important variables returned. Latitudes of transects important for analysis
            sb_found_bool, error, sbf_lon, scans_df = sb_day_analysis(rs, date)
            
            sb_found.append(sb_found_bool)
            error_log.append(error)
            
            scan_log = pd.concat([scan_log,scans_df], axis=1)
            
            #Date string
            date_string = date.strftime("%Y_%m_%d")
            
            print(date_string,error)
            
            day_log = pd.DataFrame({'Date':[date_string],
                                    'Sea_Breeze':[sb_found_bool],
                                    'Error':[error]
                                    })
            
            day_log.set_index('Date')
            
            day_log.to_csv(w_dir+"{0}_{1}_{2}DayLog.csv".\
                format(i,j,k))
            
            #Only enter this the first sb day found to initialize dataframe.
            if sb_found_bool and not first_found:
                #Prevents entering this again.
                first_found = True
                
                #Creating dataframe for overall analysis
                df = pd.DataFrame(sbf_lon)
                df = df.transpose()
                df.columns = scan_time
                df.index = cstlnlat
                
                df.to_csv(w_dir+"{0}_{1}_{2}RadarDetect.csv".\
                    format(i,j,k))
                
                print("Sea breeze detected on ", date)
                
                print("The Current time is: ", datetime.now())
                
                continue
            
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
                
                df_temp.to_csv(w_dir+"{0}_{1}_{2}RadarDetect.csv".\
                    format(i,j,k))
                
                df = pd.concat([df,df_temp], axis=1)
                
                print("Sea breeze detected on ", date)
                # print("Coordinates are ", sbf_lon)
                # print("On "+str(datelist[i]+"sb_found is "+sb_found[i]+"with sb_vel = "+sb_vel+\
                #     "and sb_pen = "+sb_pen))
                
                print("The Current time is: ", datetime.now())
            
            else:
                print("No sea breeze detected on ", date)
                
                print("The Current time is: ", datetime.now())
                
        
        if i==2013 and j==5:
            #First month analyzed, create overall dataframe for all data.
            total_df = df
            
            #Print monthly dataframes to a csv file for analysis.
            df.to_csv(w_dir+"{0}_{1}_RadarDetect.csv".format(i,j))
            scan_log.to_csv(w_dir+"{0}_{1}_ScanLog.csv".format(i,j))
            
            continue
        
        #Print monthly dataframes to a csv file for analysis.
        df.to_csv(w_dir+"{0}_{1}_RadarDetect.csv".format(i,j))
        scan_log.to_csv(w_dir+"{0}_{1}_ScanLog.csv".format(i,j))
        
        #Concatenate to overall dataframe:
        total_df = pd.concat([total_df,df], axis=1)
        

total_df.to_csv(w_dir+"2013_2017RadarDetect.csv")

#Notification that the program has finished
print('Program Finished')

