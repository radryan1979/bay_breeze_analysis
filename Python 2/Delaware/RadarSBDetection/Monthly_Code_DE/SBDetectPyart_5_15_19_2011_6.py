"""
Coded by Dan Moore
Final program to ingest 10 years worth of radar data and output
coordinates of SBF for analysis.

Updated to eliminate false detections and increase accuracy.

This will run as a test every hour between 10AM and 6PM to 
create initial dataset to analyze for quality/accuracy.

Update 4_15_19: allow for different shapes of reflectivities.
This will allow for any year (potentially) to be analyzed.

Updated to analyze KMLB (Melbourne, FL) radar.

Update 4_18_19: correction to handle errors in lats/lons.

Update 4_30_19: no update to code but previous version caught
on 6_7_2018, need to re-analyze 2018.

Update 5_15_19: Updated to conservative re-mapping for 
sampling onto rectilinear geographic grid for detection.

Analyzes KDOX Radar in Delaware.

Updated: 5_15_19
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
import xarray as xr
import xesmf as xe


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
    
    #Initialize daily precip dataframe
    cols = ['Scan_Time','Location','Mean_Reflectivity','Pixel_Area','Max_Reflectivity']
    precip_log = pd.DataFrame(index=None,columns=cols)
    precip_log.set_index('Scan_Time')

    #Number of scans actually analyzed.
    scan_num = 0
    
    #Maintain scan information for scans that are actually read by the program.
    scan_time=[]
    precip_found=[]
    sea_breeze_found=[]
    radar_mode=[]
    
    #Initializing areal location arrays:
    sbf_lon_day = []
    areal_precip_error=False
    
    #If no scans are found for the given date:
    if len(ds)==0:
        day_sb_found = False
        error = "No data found for selected date."
        
        scans_df=pd.DataFrame()
        
        return day_sb_found, error, sbf_lon_day, scans_df, scan_time, areal_precip_error
    #Eliminate days with insufficient data (total of ~45 scans per day
    #analyzed by this program).
    elif len(ds)<20:
        day_sb_found = False
        error = "Insufficient data: Not enough scans to analyze this day."
        
        scans_df=pd.DataFrame()

        return day_sb_found, error, sbf_lon_day, scans_df, scan_time, areal_precip_error
    
    #For analyzing top or half of hour:
    top_of_hour=True
        
    for number in range(13,24):
        
        #To process every half hour (twice per hour)
        for half in range(2):
        
            try: 
                #Create radar object (one every half hour - first scan of hour and then 
                #first scan after half past that hour)
                for ijk in range(len(ds)):
                    if top_of_hour:
                        if str(ds[ijk])[13:15]==str(number):
                            ds_use = ds[ijk]
                            top_of_hour=False
                            break
                    else:
                        if str(ds[ijk])[13:15]==str(number) and \
                            str(ds[ijk])[15:16]=='3':
                                ds_use = ds[ijk]
                                top_of_hour=True
                                break
                    
                radar = pyart.io.read_nexrad_cdm(ds_use.access_urls['OPENDAP'])
            except:
                print("Line 84: read_nexrad_cdm error. Please check NEXRAD server if problem persists.")
                continue
            
            #Calculate timestamp to decrease computing time.
            timestamp = radar.time['units'].split(' ')[-1].split('T')
            timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
            timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
            
            print("Scanning: {0}".format(str(timestamp)))
            
            #Determine mode of radar (nsweeps>10, precip mode)
            if radar.nsweeps > 10:
                rad_mode = 'precip'
            else:
                rad_mode = 'clear'
                
        ###########################################################################
        ###########################################################################
        #                   CHECK FOR COMMON ERROR
        ###########################################################################
        ###########################################################################
            common_error = check_lat_lon_error(radar, lon0, lat0)
            
            if common_error:
                #Erroneous values in dataset, skip.
                print("Error, erroneous values in lats/lons. Scan Skipped.")
                continue
            
        ###########################################################################
        ###########################################################################
        #                   ASSIGNING LAT/LON - POLAR GRID
        ###########################################################################
        ###########################################################################
                
            gref, glon, glat = conserve_remap_to_cartesian(radar, lon0, lat0)
            
        ###########################################################################
        ###########################################################################
        #                           HANDLING REF DATA
        ###########################################################################
        ###########################################################################
    
            #Check if there is substantial precipitation. If amount of reflectivity
            #values above 25 dBZ is more than 33% of all nonzero values in the 
            #dataset, then skip the day.
            if np.count_nonzero(gref[gref>20.0])/\
            np.count_nonzero(gref) >= 0.33:
                areal_precip_error=True
                continue
            
            if glat.shape != gref.shape:
                print("Error, data shape mismatch, scan skipped.")
                continue
    
            
            #If scan passes all the above if statements, we will evaluate it.
            prevtime = timestamp
            scan_num+=1
            
            #Utilizing function at beginning to trim the ref data.
            my_ref = trim_rad_data(glat, glon, gref, boundinglat, boundinglon)
            
            
        ###########################################################################
        ###########################################################################
        #                       ANALYZE INDIVIDUAL SCANS
        ###########################################################################
        ###########################################################################
            #######################
            #Set background reflectivity for testing later:
            onedimlat = glat[:,0]
            onedimlon = glon[0,:]
            #Box for background reflectivities
            NLatbgbox = find_nearest(onedimlat, 39.25)
            SLatbgbox = find_nearest(onedimlat, 38.25)
            ELonbgbox = find_nearest(onedimlon, -75.1)
            WLonbgbox = find_nearest(onedimlon,-76.0)
            #Average of background reflectivities:
            ave_bg_ref = np.nanmean(gref[SLatbgbox:NLatbgbox,WLonbgbox:ELonbgbox])
            ave_bg_refmed = np.nanmedian(gref[SLatbgbox:NLatbgbox,WLonbgbox:ELonbgbox])
            #######################
            
            #Trim data further, and detect precip.
            gref, sb_precip, temp_precip_df = connected_areas(gref)
            
            #Convert row_col coords to lat_lon for processing purposes.
            #Add this df as output to keep track of precip events for each scan.
            if sb_precip:
                temp_precip_df['Location'] = \
                rowcol_to_latlon(temp_precip_df['Location'], glat, glon)
                
                temp_precip_df['Scan_Time'] = [timestamp]*len(temp_precip_df)
                temp_precip_df.set_index('Scan_Time')
                
                precip_log = pd.concat([precip_log,temp_precip_df], sort=True)
        
            if num_sb_found == 0:
                sbfound, sbf = first_pass_sea_breeze_front_coords(gref,glon,glat,ave_bg_ref)
            else:
                #If SBF has been found previously on this day, use previous coordinates
                #to minimize error of erroneous SBF coordinates.
                sbfound, sbf = sea_breeze_front_coords(gref,glon,glat,sbflon,ave_bg_ref)
            
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
    #                           PRINT DAILY PRECIP LOG
    ###########################################################################
    ###########################################################################
    
    day=timestamp.day
    mo=timestamp.month
    yr=timestamp.year
    
    precip_log.to_csv(w_dir+"{0}_{1}_{2}PrecipLog.csv".format(yr,mo,day))

    ###########################################################################
    ###########################################################################
    #                           ANALYZE DAILY DATA
    ###########################################################################
    ###########################################################################
    
    if scan_num<10:
        day_sb_found = False
        error = "Insufficient data."
        
        scans_df = pd.DataFrame({'Scan_Time':scan_time,
                                'Precipitation':precip_found,
                                'Sea_Breeze':sea_breeze_found,
                                'Radar_Mode':radar_mode
                                })
        scans_df.set_index('Scan_Time')
        
        return day_sb_found, error, sbf_lon_day, scans_df, scan_time, areal_precip_error
        
    
    if num_sb_found<4:
        day_sb_found = False
        error = "Sufficient data, but sea breeze not found."
        
        scans_df = pd.DataFrame({'Scan_Time':scan_time,
                                'Precipitation':precip_found,
                                'Sea_Breeze':sea_breeze_found,
                                'Radar_Mode':radar_mode
                                })
        scans_df.set_index('Scan_Time')
        
        return day_sb_found, error, sbf_lon_day, scans_df, scan_time, areal_precip_error

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
        
    return day_sb_found, error, sbf_lon_day, scans_df, scan_time, areal_precip_error


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

#Check to see if data is flawed: has lats/lons>100
def check_lat_lon_error(radar, lon0, lat0):
    display = pyart.graph.RadarMapDisplay(radar)
    
    x,y = display._get_x_y(0,False,None)
    x_b,y_b = display._get_x_y(0,True,None)
    x = x*1000; y = y*1000; x_b=x_b*1000; y_b=y_b*1000
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
    lons_b, lats_b = pyart.core.cartesian_to_geographic(x_b,y_b,
        projparams = lambert_aea)
    
    if lons.max()>100.0 or lats.max()>100.0 or\
        lons_b.max()>100.0 or lats_b.max()>100.0:
            error=True #Erroneous values in dataset, skip.
    else:
        error=False
    
    return error
    
#Regrid data to cartesian grid with polar structure, in other words
#i and j of 2d array are scan angle and range, respectively.
def conserve_remap_to_cartesian(radar, lon0, lat0):
    display = pyart.graph.RadarMapDisplay(radar)
    
    x,y = display._get_x_y(0,False,None)
    x_b,y_b = display._get_x_y(0,True,None)
    x = x*1000; y = y*1000; x_b=x_b*1000; y_b=y_b*1000
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
    lons_b, lats_b = pyart.core.cartesian_to_geographic(x_b,y_b,
        projparams = lambert_aea)
    
    ref = radar.get_field(0, 'reflectivity')
    
    ds=xr.Dataset({'reflectivity': (['i', 'j'],  ref)},
            coords={'lon': (['i', 'j'], lons),
            'lat': (['i', 'j'], lats)})
    
    grid_in = {'lon': lons, 'lat': lats,
        'lon_b': lons_b, 'lat_b': lats_b}
    
    regridder = xe.Regridder(grid_in, out_grid, 'conservative')
    data_out=regridder(ref)
    regridder.clean_weight_file()
    
    return data_out, out_grid.lon.values, out_grid.lat.values

#Convert polar structure grid to rectilinear - NO LONGER USED (5-15-19)
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
def first_pass_sea_breeze_front_coords(gref, glon, glat, ave_bg_ref):
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
    
    #Index for runs. will eliminate short runs (outliers)
    sbf_run_idx_list=[]
    skip=0
    run=0
    
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
        
        #Loop through previous j values - essentially moving from west to east
        for j in range(idxlon,idxlonend,-1):
            #Test reflectivities surrounding this point to eliminate data
            sumref = np.count_nonzero(~np.isnan(gref[idxlat-5:idxlat+6,j-5:j+6]))
            
            #If larger cluster found, replace previous.
            if sumref>=60:
                
                #Average nugget reflectivity:
                ave_nug_ref = np.nanmean(gref[idxlat-5:idxlat+6,j-5:j+6])
                
                ave_bg_refW = np.nanmean(gref[idxlat-5:idxlat+6,j-10:j-4])
                
                ave_bg_refE = np.nanmean(gref[idxlat-5:idxlat+6,j+5:j+11])
                
                #Added conditional to be within 0.3degrees of coastline - empirical.
                if sumref>maxsumref and\
                    ave_nug_ref>ave_bg_refW and ave_nug_ref>ave_bg_refE and\
                    ave_nug_ref>15.0:# and \
                    # (onedimlon[idxlon]-onedimlon[j])<0.3:
                        maxsumref = sumref
                        maxrefj = j
        
        #If clustering sufficient for SB was found, place SB point
        if maxsumref>0:
            #Saving coordinates
            sbflat.append(glat[idxlat,maxrefj])
            sbflon.append(glon[idxlat,maxrefj])
            sbf_idx_list.append(i)
            sbf_run_idx_list.append(i)
            run+=1
            skip=0
            sbfcnt+=1
        
        else:
            sbflat.append(glat[idxlat,0])
            sbflon.append(np.nan)
            skip+=1
            if skip>=3:
                if run>=4:
                    #keep points in this run because it is long enough
                    sbf_run_idx_list=[]
                    pass
                else:
                    #those are erroneous values, get rid of points
                    for erase in sbf_run_idx_list:
                        sbflon[erase] = np.nan
                    sbf_idx_list = list(set(sbf_idx_list)^set(sbf_run_idx_list))
                    sbf_run_idx_list=[] #reset this look for next run
                    pass
                run=0

                
        
    
    #If fewer than 34% of the transects (10 out of 30), 
    #then we call this an incomprehensible SB or not found.
    if sbfcnt/float(len(cstlnlat))>0.19:
        #Place sbf coords on transects that did not find sbf.
        for i in range(len(sbf_idx_list)-1):
            points = [(sbflon[sbf_idx_list[i]],sbflat[sbf_idx_list[i]]),
                    (sbflon[sbf_idx_list[i+1]],sbflat[sbf_idx_list[i+1]])]
            x_coords, y_coords = zip(*points)
            A = np.vstack([x_coords,np.ones(len(x_coords))]).T
            m, b = np.linalg.lstsq(A, y_coords)[0]
            
            #Difference between indices of consecutively found SBF points.
            idx_diff = sbf_idx_list[i+1] - sbf_idx_list[i]
            
            #If idx_diff==1, then they are consecutive transects, skip.
            #If idx_diff>3, then they are probably two separate disconnected branches.
            if idx_diff == 1 or idx_diff>3:
                continue
            
            for j in range(sbf_idx_list[i]+1,sbf_idx_list[i+1]):
                sbflon[j] = (sbflat[j]-b)/float(m)
        
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
def sea_breeze_front_coords(gref, glon, glat, prev_lon, ave_bg_ref):
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
    
    #Index for runs. will eliminate short runs (outliers)
    sbf_run_idx_list=[]
    skip=0
    run=0
    
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
        
        #Loop through previous j values - essentially moving from west to east
        for j in range(idxlon,idxlonend,-1):
            #Test reflectivities surrounding this point to eliminate data
            sumref = np.count_nonzero(~np.isnan(gref[idxlat-5:idxlat+6,j-5:j+6]))
            
            #If larger cluster found, replace previous.
            if sumref>=60:
                
                #Average nugget reflectivity:
                ave_nug_ref = np.nanmean(gref[idxlat-5:idxlat+6,j-5:j+6])
                
                ave_bg_refW = np.nanmean(gref[idxlat-5:idxlat+6,j-10:j-4])
                
                ave_bg_refE = np.nanmean(gref[idxlat-5:idxlat+6,j+5:j+11])
                
                #Added conditional to be within 0.3degrees of coastline - empirical.
                if sumref>maxsumref and\
                    ave_nug_ref>ave_bg_refW and ave_nug_ref>ave_bg_refE and\
                    ave_nug_ref>15.0:# and \
                    # (onedimlon[idxlon]-onedimlon[j])<0.3:

                        maxsumref = sumref
                        maxrefj = j
        
        #If clustering sufficient for SB was found, place SB point
        if maxsumref>0:
            #Saving coordinates
            sbflat.append(glat[idxlat,maxrefj])
            sbflon.append(glon[idxlat,maxrefj])
            sbf_idx_list.append(i)
            sbf_run_idx_list.append(i)
            run+=1
            skip=0
            sbfcnt+=1
        
        else:
            sbflat.append(glat[idxlat,0])
            sbflon.append(np.nan)
            skip+=1
            if skip>=3:
                if run>=4:
                    #keep points in this run because it is long enough
                    sbf_run_idx_list=[]
                    pass
                else:
                    #those are erroneous values, get rid of points
                    for erase in sbf_run_idx_list:
                        sbflon[erase] = np.nan
                    sbf_idx_list = list(set(sbf_idx_list)^set(sbf_run_idx_list))
                    sbf_run_idx_list=[] #reset this look for next run
                    pass
                run=0

                
        
    
    #If fewer than 34% of the transects (10 out of 30), 
    #then we call this an incomprehensible SB or not found.
    if sbfcnt/float(len(cstlnlat))>0.19:
        #Place sbf coords on transects that did not find sbf.
        for i in range(len(sbf_idx_list)-1):
            points = [(sbflon[sbf_idx_list[i]],sbflat[sbf_idx_list[i]]),
                    (sbflon[sbf_idx_list[i+1]],sbflat[sbf_idx_list[i+1]])]
            x_coords, y_coords = zip(*points)
            A = np.vstack([x_coords,np.ones(len(x_coords))]).T
            m, b = np.linalg.lstsq(A, y_coords)[0]
            
            #Difference between indices of consecutively found SBF points.
            idx_diff = sbf_idx_list[i+1] - sbf_idx_list[i]
            
            #If idx_diff==1, then they are consecutive transects, skip.
            #If idx_diff>3, then they are probably two separate disconnected branches.
            if idx_diff == 1 or idx_diff>3:
                continue
            
            for j in range(sbf_idx_list[i]+1,sbf_idx_list[i+1]):
                sbflon[j] = (sbflat[j]-b)/float(m)
        
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
    
    s = ndimage.generate_binary_structure(2,2)
    labeled, ncomponents = ndimage.measurements.label(gref2,structure=s)

    means = ndimage.mean(gref2, labeled, range(1,ncomponents + 1))
    maxs = ndimage.maximum(gref2, labeled, range(1,ncomponents + 1))
    
    areas = [r.filled_area for r in regionprops(labeled)]
    # perimeters = [r.perimeter for r in regionprops(labeled)]
    labels = [r.label for r in regionprops(labeled)]
    # eccentricity = [r.eccentricity for r in regionprops(labeled)]
    # orientation = [r.orientation for r in regionprops(labeled)]
    locs = [r.centroid for r in regionprops(labeled)]
    maj_ax_len = [r.major_axis_length for r in regionprops(labeled)]
    
    precip_area = []
    precip_ref = []
    precip_xyloc = []
    precip_max = []
    
    for i in range(len(areas)):
        
        # p_a_ratio = float(perimeters[i])/float(areas[i])
        idx = labeled==labels[i]
        
        if maj_ax_len[i]==0.0:
            labeled[idx]=0
            gref2[idx]=np.nan
            continue

        comp_ratio = 1.2732*areas[i]/maj_ax_len[i]**2
        
        if areas[i]>100 and means[i]>22.5:
            precip_area.append(areas[i])
            precip_ref.append(means[i])
            precip_xyloc.append(locs[i])
            precip_max.append(maxs[i])
            precip=True
            labeled[idx]=0
            gref2[idx]=np.nan
        elif means[i]>25.0:
            labeled[idx]=0
            gref2[idx]=np.nan
        elif areas[i]<10:
            labeled[idx]=0
            gref2[idx]=np.nan
        else:
            pass
        
        if areas[i]>200 and comp_ratio<0.30 and means[i]<22.5:
            #Potential Sea breeze front detected.
            pass
        else:
            labeled[idx]=0
            gref2[idx]=np.nan
    
    gref2[gref2==0]=np.nan
    

    # labels = np.unique(labeled)
    # labeled = np.searchsorted(labels, labeled)
    
    # plt.imsave('/Users/dpmoore2927/Desktop/test2.png',labeled)
    
    if precip:
        precip_df = pd.DataFrame({'Location':precip_xyloc,
                            'Mean_Reflectivity':precip_ref,
                            'Pixel_Area':precip_area,
                            'Max_Reflectivity':precip_max
                            })
    else:
        precip_df = pd.DataFrame()
    
    return gref2, precip, precip_df

def rowcol_to_latlon(rowcol,lats,lons):
    
    latlon = []
    for rc in rowcol:
        latlon.append((lats[int(rc[0]),0],lons[0,int(rc[1])]))
    
    return latlon



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

#Grid x-array to be used for convservative re-mapping:
out_grid = xe.util.grid_2d(-76.45, -74.6, 0.002, 38.15, 39.65, 0.002)

#Reflectivity thresholds
loref = 10.0; hiref = 25.0 #dBz thresholds

###########################################################################
###########################################################################
#                               READ IN DATA
###########################################################################
###########################################################################

#Read in data from unidata server - option to loop through.
#Eventually write this is a function to be called every day of interest 
#to analyze for sea breeze and spit out information.

#http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/
#http://149.165.168.53/thredds/radarServer/nexrad/level2/S3/
rs=RadarServer('http://149.165.168.53/thredds/radarServer/nexrad/level2/S3/')

w_dir = "/home/1898/RadSBDetect/SBDetect5_15_19/"

sb_found = []
error_log = []

cols = ['Scan_Time','Precipitation','Sea_Breeze','Radar_Mode']

scan_log = pd.DataFrame(index=None,columns=cols)

scan_log.set_index('Scan_Time')

first_month=True


#Write for loop to go through range of dates over time period:

for i in range(2011,2012):
    
    print("Processing Year: ", i)
    print("The Current time is: ", datetime.now())
    
    for j in range(6,7):
        
        #Initialize dataframe every new month.
        df = pd.DataFrame()
        
        #Enter conditional below for initializing pandas dataframe
        first_found = False
        
        for k in range(1,32):
            
            #Make sure it skips the 31st day if it is June or September.
            if (j==6 or j==9) and k==31:
                continue
    
            date = datetime(i, j, k)
            
            #Important variables returned. Latitudes of transects important for analysis
            sb_found_bool, error, sbf_lon, scans_df, scan_time,\
            areal_precip_error = sb_day_analysis(rs, date)
            
            sb_found.append(sb_found_bool)
            error_log.append(error)
            
            scan_log = pd.concat([scan_log,scans_df])
            
            #Date string
            date_string = date.strftime("%Y_%m_%d")
            
            print(date_string,error)
            
            day_log = pd.DataFrame({'Date':[date_string],
                                    'Sea_Breeze':[sb_found_bool],
                                    'Error':[error],
                                    'Areal_Precip_Error':[areal_precip_error]
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
                
        
        if first_month:
            first_month=False
            #First month analyzed, create overall dataframe for all data.
            total_df = df.copy()
            
            #Print monthly dataframes to a csv file for analysis.
            df.to_csv(w_dir+"{0}_{1}_RadarDetect.csv".format(i,j))
            scan_log.to_csv(w_dir+"{0}_{1}_ScanLog.csv".format(i,j))
            
            continue
        
        #Print monthly dataframes to a csv file for analysis.
        df.to_csv(w_dir+"{0}_{1}_RadarDetect.csv".format(i,j))
        scan_log.to_csv(w_dir+"{0}_{1}_ScanLog.csv".format(i,j))
        
        #Concatenate to overall dataframe:
        total_df = pd.concat([total_df,df], axis=1)
        
df.to_csv(w_dir+"{0}_{1}_RadarDetect.csv".format(i,j))
scan_log.to_csv(w_dir+"{0}_{1}_ScanLog.csv".format(i,j))

# total_df.to_csv(w_dir+"2009_2018FLRadarDetect.csv")

#Notification that the program has finished
print('Program Finished')

