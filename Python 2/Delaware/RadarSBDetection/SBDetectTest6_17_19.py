"""
Coded by Dan Moore
Program to show method used to eliminate presence of precip in radar
scan, in order to distinguish sea breeze front from precip. Both of 
these reflect areas of higher reflectivity values compared with the 
background reflectivity values.

The method is shown to work in both clear air mode and precip mode.
The file used is simply modified.

Update 4-15_19: updated to analyze KMLB, Florida

Update 5-14-19: Include conservative re-mapping for detection.

Update 6-17-19: Change gradient test to include eliminated reflectivities
Yet to be implemented into actual algorithm.

Also updated to higher reflectivities on the same transect take priority
as SBF points.

Updated: 6-17-19
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
import cartopy.crs as ccrs
import cartopy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
from scipy import spatial
from scipy.interpolate import griddata
import sys
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

def sb_day_analysis(rs, date, station):
    
    #If function is passed a blank 'rs', then it uses the local file here.
    #For the purpose of this assignment, only local files will be used.
    
    #Precip case: 2017_7_4_CaseStudy.nc
    #           2017_7_4_1820CaseStudy.nc
    #Clear case: 2016_8_19_CaseStudy.nc
    #           2016_8_19_2022CaseStudy.nc
    #False case: 2016_8_20_1523CaseStudy.nc
    #Bad Data: 2016_8_12CaseStudy.nc
    
    #FL CASES:
    #Clear Case: 2013_8_16FLCaseStudy.nc
    #Problem Case: 2018_6_7CaseStudy.nc 

    
    if rs == '':
        radar = pyart.io.read_cfradial(w_dir+'2016_8_19_CaseStudy.nc')
        
    else:
    
        query=rs.query()
        query.stations(station).time(date)
        catalog = rs.get_catalog(query)
        dataset = list(catalog.datasets.values())[0]
        
        # now that we have the data, let's go ahead and open it with pyart
        radar = pyart.io.read_nexrad_cdm(dataset.access_urls['OPENDAP'])
    
    #Establish location of radar
    loc = pyart.io.nexrad_common.get_nexrad_location(station)
    lon0 = loc[1] ; lat0 = loc[0]
    
    #Number of sea breeze front successfully identified. Resets every
    #day/every call of this function.
    num_sb_found=0

    #Number of sea breeze front successfully identified. Resets every
    #day/every call of this function.
    num_sb_found = 0
    
    #Maintain scan times that are actually read by the program.
    scan_time=[]
    
    #Initializing areal location arrays:
    sbf_lon_day = []
        
    #Determine mode of radar (nsweeps>10, precip mode)
    if radar.nsweeps > 10:
        rad_mode = 'precip'
    else:
        rad_mode = 'clear'

    #Calculate timestamp to decrease computing time.
    timestamp = radar.time['units'].split(' ')[-1].split('T')
    timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
    timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
    
    ###########################################################################
    ###########################################################################
    #                   CHECK FOR COMMON ERROR
    ###########################################################################
    ###########################################################################
    common_error = check_lat_lon_error(radar, lon0, lat0)
    
    if common_error:
        #Erroneous values in dataset, skip.
        print("Error, erroneous values in lats/lons. Scan Skipped.")
        sys.exit()
            
        
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
        
        
    #Organizing reflectivity data
    # ref = radar.get_field(0, 'reflectivity')
    # unmasked = ma.getdata(ref)  
    
    #Check if there is substantial precipitation. If amount of reflectivity
    #values above 25 dBZ is more than 10% of all nonzero values in the 
    #dataset, then skip the day.
    if np.count_nonzero(gref[gref>20.0])/\
    np.count_nonzero(gref) >= 0.10:
        sbf_lon_day.append([np.nan]*len(cstlnlon))
        print('too much precip, day not analyzed.')
    
    print(glat.shape, gref.shape)
    
    if glat.shape != gref.shape:
        print("Error, data shape mismatch, scan skipped.")
        return
    
    # og_unmasked = np.array(unmasked)
    plot_sbf(timestamp, gref, glon, glat,"plot1")
    
    bg_gref = gref.copy()
    
    #Utilizing function at beginning to trim the ref data.
    my_ref = trim_rad_data(glat, glon, gref, boundinglat, boundinglon)
        
    ###########################################################################
    ###########################################################################
    #                       RE-GRIDDING DATA TO RECTILINEAR GRID
    ###########################################################################
    ###########################################################################
        
    # og_gref, og_glon, og_glat = grid_to_rectilinear(og_unmasked,lons,lats,nlon,nlat)
    plot_sbf(timestamp, gref, glon, glat,"plot2")
    
    
    
    # gref, glon, glat = grid_to_rectilinear(unmasked,lons,lats,nlon,nlat)
    
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
    
    
        
    ###########################################################################
    ###########################################################################
    #                       ANALYZE INDIVIDUAL SCANS
    ###########################################################################
    ###########################################################################
    
    gref, sb_precip, temp_precip_df = connected_areas(gref)
    
    plot_sbf(timestamp, gref, glon, glat,"plot3")
    plot_sbf(timestamp, gref, glon, glat,"plot4") # with lines
    
    
    
    #Convert row_col coords to lat_lon for processing purposes.
    #Add this df as output to keep track of precip events for each scan.
    if sb_precip:
        temp_precip_df['Location'] = \
        rowcol_to_latlon(temp_precip_df['Location'], glat, glon)
        
        temp_precip_df['Scan_Time'] = [timestamp]*len(temp_precip_df)
        temp_precip_df.set_index('Scan_Time')

    ###########################################################################
    ###########################################################################
    #                           PRINT DAILY PRECIP LOG
    ###########################################################################
    ###########################################################################
    
    day=timestamp.day
    mo=timestamp.month
    yr=timestamp.year
    
    temp_precip_df.to_csv(w_dir+"{0}_{1}_{2}PrecipLog.csv".format(yr,mo,day))
    
    sbfound, sbf, sbfinterp = first_pass_sea_breeze_front_coords(gref,glon,glat,ave_bg_ref,
                                                        ave_bg_refmed,bg_gref)#,rad_mode)
                                                        
    plot_sbf_zoom(timestamp, gref, glon, glat,"plot5") # zoom to region surrounding sbf point
    
    print('sbfound =',sbfound)
    
    plot_sbf(timestamp, gref, glon, glat, "plot6", sbf, sbfinterp)
    
    return sbfound, sbf
    

        


#Function to trim radar data to desired bounding rectangle:
#Required are 1D arrays of latitudes, longitudes, and reflectivities
#and two tuples of latitudes and longitudes defining the boundary of 
#the desired boundary.
def trim_rad_data(lats, lons, ref, boundinglat, boundinglon):
    for i in range(lats.shape[0]):
        for j in range(lats.shape[1]):
            # if (lats[i][j]>boundinglat[0] and lats[i][j]<boundinglat[1] and \
            #     lons[i][j]>boundinglon[0] and lons[i][j]<boundinglon[1]):
            #     pass
            # else:
            #     # lats[i] = np.nan; lons[i] = np.nan
            #     ref[i][j] = np.nan
                
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

#Convert polar structure grid to rectilinear
# def grid_to_rectilinear(my_ref,lons,lats,nlon,nlat):
#     rav_lats = lats.ravel()
#     rav_lons = lons.ravel()
#     rav_ref = my_ref.ravel()
    
#     #Grid Data using matplotlib
#     grid_lons = np.linspace(boundinglon[0],boundinglon[1],nlon)
#     grid_lats = np.linspace(boundinglat[0],boundinglat[1],nlat)
#     glon,glat = np.meshgrid(grid_lons,grid_lats)
    
#     # Interpolate data onto grid using linear interpolation
#     gref = griddata((rav_lons,rav_lats),rav_ref,(glon,glat),method='linear')
    
#     return gref, glon, glat

#Calculate sea breeze front coordinates, raise warning if 
#one is not found.
def first_pass_sea_breeze_front_coords(gref, glon, glat,ave_bg_ref,ave_bg_refmed,bg_gref):
    #For easier processing:
    onedimlat = glat[:,0]
    onedimlon = glon[0,:]
    
    #TEST - make background reflectivities 0 if negative
    bg_gref[bg_gref<0.0] = 0.0
    
    #Initialize sea breeze indexer
    sbfcnt=0
    loop=0
    #Initialize coordinate arrays for sea breeze front
    sbflat=[]
    sbflon=[]
    
    #Index list for knowing where a SBF was detected
    sbf_idx_list=[]
    
    #Index for runs. will eliminate short runs (outliers)
    sbf_run_idx_list=[]######
    skip=0#####
    run=0######
    
    #Set end of detection
    idxlonend = find_nearest(onedimlon,-81.4)
    
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
        
        ave_nug_refFINAL=0.0
        
        #Loop through previous j values - essentially moving from west to east
        for j in range(idxlon,idxlonend,-1):
            #Test reflectivities surrounding this point to eliminate data
            sumref = np.count_nonzero(~np.isnan(gref[idxlat-5:idxlat+6,j-5:j+6]))

            
            #If larger cluster found, replace previous.
            if sumref>=60:
                
                #Average nugget reflectivity:
                ave_nug_ref = np.nanmean(gref[idxlat-5:idxlat+6,j-5:j+6])
                
                ave_bg_refW = np.nanmean(bg_gref[idxlat-5:idxlat+6,j-10:j-4])
                
                ave_bg_refE = np.nanmean(bg_gref[idxlat-5:idxlat+6,j+5:j+11])
                
                #Added conditional to be within 0.3degrees of coastline - empirical.
                if  sumref>maxsumref and\
                    ave_nug_ref>ave_bg_refW and ave_nug_ref>ave_bg_refE and\
                    ave_nug_ref>10.0 and ave_nug_ref>ave_nug_refFINAL:# and \
                    # (onedimlon[idxlon]-onedimlon[j])<0.3:
                        maxsumref = sumref
                        maxrefj = j
                        
                        ave_nug_refFINAL = ave_nug_ref
                        ave_bg_refWFINAL = ave_bg_refW
                        ave_bg_refEFINAL = ave_bg_refE
        
        #If clustering sufficient for SB was found, place SB point
        if maxsumref>0:
            print(loop)
            #Saving coordinates
            sbflat.append(glat[idxlat,maxrefj])
            sbflon.append(glon[idxlat,maxrefj])
            sbf_idx_list.append(i)
            sbf_run_idx_list.append(i)######
            run+=1#####
            skip=0#####
            sbfcnt+=1
            
            print("nug", ave_nug_refFINAL)
            print("West", ave_bg_refWFINAL)
            print("East", ave_bg_refEFINAL)
            print("Lat" , sbflat[loop])
            print("Lon" , sbflon[loop])
        
        else:
            print(loop)
            print("Not successful")
            sbflat.append(glat[idxlat,0])
            sbflon.append(np.nan)
            skip+=1#####
            if skip>=3:######
                if run>=4:#####
                    #keep points
                    sbf_run_idx_list=[]
                    pass####
                else:#####
                    #those are erroneous values, get rid of points
                    print("outlier",sbf_run_idx_list)
                    for erase in sbf_run_idx_list:
                        sbflon[erase] = np.nan
                    sbf_run_idx_list=[] #reset this
                    pass
                run=0#####
                
                
                
        loop+=1
        #Catch runs smaller than 4 SBF points (out of 8 lines) so that 
        #points are not mistaken as SBF points.

                
    print("ratio",sbfcnt/len(cstlnlat))
    print("ratio2", sbfcnt/float(len(cstlnlat)))
    
    sbfloninterp=[]
    sbflatinterp=[]
    
    #If fewer than 34% of the transects (10 out of 30), 
    #then we call this an incomprehensible SB or not found.
    if sbfcnt/len(cstlnlat)>0.19:
        
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
                sbflon[j] = (sbflat[j]-b)/m
                sbfloninterp.append(sbflon[j])
                sbflatinterp.append(sbflat[j])
        
        sbfound = True
        #Create tuple of coordinates    
        sbf = list(zip(sbflon,sbflat))
        sbfinterp = list(zip(sbfloninterp,sbflatinterp))
    else:
        sbfound = False
        sbf = list(zip(sbflon,sbflat))
        sbfinterp = list(zip(sbfloninterp,sbflatinterp))
        # sbf = ()
    
    return sbfound, sbf, sbfinterp
    
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
    
    #Box for background reflectivities
    NLatbgbox = find_nearest(onedimlat, 39.25)
    SLatbgbox = find_nearest(onedimlat, 38.25)
    ELonbgbox = find_nearest(onedimlon, -75.1)
    WLonbgbox = idxlonend
    
    #Average of background reflectivities:
    ave_bg_ref = np.nanmean(gref[SLatbgbox:NLatbgbox,WLonbgbox:ELonbgbox])
    
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
                
                #Average nugget reflectivity:
                ave_nug_ref = np.nanmean(gref[idxlat-5:idxlat+5,j-5:j+5])
                
                try:
                    #Added conditional to be within 0.3degrees of previous point - empirical.
                    if ave_nug_ref>ave_bg_ref and sumref>maxsumref and \
                        abs(prev_lon[j]-onedimlon[j])<0.25:
                            maxsumref = sumref
                            maxrefj = j
                            
                except:
                    if ave_nug_ref>ave_bg_ref and sumref>maxsumref: 
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
    if sbfcnt/len(cstlnlat)>0.20:
        
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
                sbflon[j] = (sbflat[j]-b)/m
        
        sbfound = True
        #Create tuple of coordinates    
        sbf = list(zip(sbflon,sbflat))
    else:
        sbfound = False
        sbf = list(zip(sbflon,sbflat))
        # sbf = ()
    
    return sbfound, sbf, sbf
    

def plot_sbf(date, gref, glon, glat, title, sbf=(), sbfinterp=()):
    
    if title == "plot1" or title == "plot2" or title == "plot3":
        lines=False
    else:
        lines=True
    
    #For easier processing:
    onedimlon = glon[0,:]
    onedimlat = glat[:,0]
    
    translonE = np.zeros(len(cstlnlon))
    
    #Set end of detection (-81.4 for FL, -76.0 for DE)
    idxlonend = find_nearest(onedimlon,-76.0)
    translonW = onedimlon[idxlonend]
    
    #Loop through transect lines
    for i in range(np.size(cstlnlat)):
        tlon0 = cstlnlon[i]+0.01 #Just a bit offshore in case
        
        #Find index of nearest longitude in grid: this is our i value
        idxlon = find_nearest(onedimlon,tlon0)#represents starting column
        
        translonE[i] = onedimlon[idxlon]
    
    
    fig = plt.figure(figsize=[8, 8], dpi=100)
    my_ax = plt.axes(projection = ccrs.PlateCarree())
    
    r = my_ax.pcolormesh(glon, glat, gref, cmap=pyart.graph.cm.NWSRef,
                        norm=mpl.colors.Normalize(vmin=0,vmax=60))
    cbar = fig.colorbar(r, orientation='vertical')
    cbar.set_label('Equivalent Reflectivity Factor (dBZ)', rotation=90)
    
    states = cartopy.feature.NaturalEarthFeature(category='cultural',
                                  name='admin_1_states_provinces_lines',
                                  scale='10m', facecolor='none')
    coast = cartopy.feature.NaturalEarthFeature(category='physical', scale='10m',
                                facecolor='none', name='coastline')
                                
    reader = shpreader.Reader('/Users/dpmoore2927/Python3/countyl010g_shp_nt00964/countyl010g.shp')
    counties = list(reader.geometries())
    COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())
    
    my_ax.add_feature(COUNTIES, facecolor='none', edgecolor='lightgray')                            
    my_ax.add_feature(states, linestyle='-', edgecolor='darkslategray',linewidth=1)
    my_ax.add_feature(coast, linestyle='-', edgecolor='darkslategray',linewidth=1)
    
    plt.title('KDOX NEXRAD {0}'.format(date))
    my_ax.set(aspect=1,
            xlim=(glon[0][0],glon[0][-1]),
            ylim=(glat[0][1],glat[-1][0]))
    
    if lines:
        for i in range(len(cstlnlat)):
            plt.plot([translonW, translonE[i]], [cstlnlat[i], cstlnlat[i]],
                 color='purple',
                 linewidth=0.6,
                 transform=ccrs.PlateCarree(),
                 )
    
    #Plot SBF
    for pt in range(len(sbf)):
        plt.plot(sbf[pt][0],sbf[pt][1],'rp',markersize=4)
        # for i in range(-10,11,1):
        #     lon1 = find_nearest(onedimlon,sbf[pt][0])
            
        #     plt.plot(onedimlon[lon1+i],sbf[pt][1], 'bp', markersize=1)
    
    #Plot interpolated SBF points
    # for pt in range(len(sbfinterp)):
    #     plt.plot(sbfinterp[pt][0],sbfinterp[pt][1],'bp',markersize=4)
    
    #FLORIDA
    # plt.xticks(np.arange(-79.8,-81.8,-0.25))
    # plt.yticks(np.arange(27.375,28.9,0.25))
    
    #DELAWARE
    plt.xticks(np.arange(-74.75,-76.50,-0.25))
    plt.yticks(np.arange(38.25,39.75,0.25))
    
    plt.savefig(w_dir+'{0}.png'.format(title),transparent=True,dpi=500)
    return

def plot_sbf_zoom(date, gref, glon, glat, title, sbf=(), sbfinterp=()):
    
    if title == "plot1" or title == "plot2" or title == "plot3":
        lines=False
    else:
        lines=True
    
    #For easier processing:
    onedimlon = glon[0,:]
    onedimlat = glat[:,0]
    
    translonE = np.zeros(len(cstlnlon))
    
    #Set end of detection (-81.4 for FL, -76.0 for DE)
    idxlonend = find_nearest(onedimlon,-76.0)
    translonW = onedimlon[idxlonend]
    
    #Loop through transect lines
    for i in range(np.size(cstlnlat)):
        tlon0 = cstlnlon[i]+0.01 #Just a bit offshore in case
        
        #Find index of nearest longitude in grid: this is our i value
        idxlon = find_nearest(onedimlon,tlon0)#represents starting column
        
        translonE[i] = onedimlon[idxlon]
    
    
    fig = plt.figure(figsize=[8, 8], dpi=100)
    my_ax = plt.axes(projection = ccrs.PlateCarree())
    
    r = my_ax.pcolormesh(glon, glat, gref, cmap=pyart.graph.cm.NWSRef,
                        norm=mpl.colors.Normalize(vmin=0,vmax=60))
    cbar = fig.colorbar(r, orientation='vertical')
    cbar.set_label('Equivalent Reflectivity Factor (dBZ)', rotation=90)
    
    states = cartopy.feature.NaturalEarthFeature(category='cultural',
                                  name='admin_1_states_provinces_lines',
                                  scale='10m', facecolor='none')
    coast = cartopy.feature.NaturalEarthFeature(category='physical', scale='10m',
                                facecolor='none', name='coastline')
                                
    reader = shpreader.Reader('/Users/dpmoore2927/Python3/countyl010g_shp_nt00964/countyl010g.shp')
    counties = list(reader.geometries())
    COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())
    
    my_ax.add_feature(COUNTIES, facecolor='none', edgecolor='lightgray')                            
    my_ax.add_feature(states, linestyle='-', edgecolor='darkslategray',linewidth=1)
    my_ax.add_feature(coast, linestyle='-', edgecolor='darkslategray',linewidth=1)
    
    #Plot rectangle around nug:
    lon1 = onedimlon[find_nearest(onedimlon,sbf[30][0])-5]
    lat1 = onedimlat[find_nearest(onedimlat,sbf[30][1])-5]
    nug = mpl.patches.Rectangle((lon1,lat1), 10*0.002, 10*0.002,ec='blue',fill=False)
    my_ax.add_patch(nug)
    
    #Plot rectangle of gradients:
    lon2 = onedimlon[find_nearest(onedimlon,sbf[30][0])-10]
    lat2 = onedimlat[find_nearest(onedimlat,sbf[30][1])-5]
    grad1 = mpl.patches.Rectangle((lon2,lat2), 5*0.002, 10*0.002,ec='red',fill=False)
    my_ax.add_patch(grad1)
    lon3 = onedimlon[find_nearest(onedimlon,sbf[30][0])+5]
    lat3 = onedimlat[find_nearest(onedimlat,sbf[30][1])-5]
    grad2 = mpl.patches.Rectangle((lon3,lat3), 5*0.002, 10*0.002,ec='red',fill=False)
    my_ax.add_patch(grad2)
    
    plt.title('KDOX NEXRAD {0}'.format(date))
    my_ax.set(aspect=1,
            xlim=(lon1-0.02,lon1+0.04),
            ylim=(lat1-0.02,lat1+0.04))
    
    if lines:
        for i in range(len(cstlnlat)):
            plt.plot([translonW, translonE[i]], [cstlnlat[i], cstlnlat[i]],
                 color='purple',
                 linewidth=0.6,
                 transform=ccrs.PlateCarree(),
                 )
    
    #Plot SBF
    # for pt in range(len(sbf)):
    #     # plt.plot(sbf[pt][0],sbf[pt][1],'rp',markersize=4)
    #     for i in range(-5,6,1):
    #         for j in range(-5,6,1):
    #             
    #             lon1 = find_nearest(onedimlon,sbf[pt][0])
                
    #             plt.plot(onedimlon[lon1+i],sbf[pt][1], 'bp', markersize=1)
    
    #Plot interpolated SBF points
    # for pt in range(len(sbfinterp)):
    #     plt.plot(sbfinterp[pt][0],sbfinterp[pt][1],'bp',markersize=4)
    
    #FLORIDA
    # plt.xticks(np.arange(-79.8,-81.8,-0.25))
    # plt.yticks(np.arange(27.375,28.9,0.25))
    
    #DELAWARE
    # plt.xticks(np.arange(-74.75,-76.50,-0.25))
    # plt.yticks(np.arange(38.25,39.75,0.25))
    
    plt.savefig(w_dir+'{0}.png'.format(title),transparent=True,dpi=500)
    return


def connected_areas(gref2):
    
    gref2[np.isnan(gref2)]=0
    
    #Determine if there is precip in the scan.
    precip=False
    
    s = ndimage.generate_binary_structure(2,2)
    labeled, ncomponents = ndimage.measurements.label(gref2,structure=s)
    

    means = ndimage.mean(gref2, labeled, range(1,ncomponents + 1))
    maxs = ndimage.maximum(gref2, labeled, range(1,ncomponents + 1))
    
    areas = [r.filled_area for r in regionprops(labeled)]
    perimeters = [r.perimeter for r in regionprops(labeled)]
    labels = [r.label for r in regionprops(labeled)]
    eccentricity = [r.eccentricity for r in regionprops(labeled)]
    orientation = [r.orientation for r in regionprops(labeled)]
    locs = [r.centroid for r in regionprops(labeled)]
    maj_ax_len = [r.major_axis_length for r in regionprops(labeled)]
    
    precip_area = []
    precip_ref = []
    precip_xyloc = []
    precip_max = []
    
    for i in range(len(areas)):
        
        p_a_ratio = float(perimeters[i])/float(areas[i])
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

station='KDOX'

if station=='KMLB':
    #Boundaries for desired values
    boundinglat = (27.375, 28.875)
    boundinglon = (-81.725, -79.875)
    
    #Defining the FL Coastline
    cstlnlat = np.arange(28.65,27.64,-0.025)
    
    cstlnlon = np.arange(-80.57,-80.14,0.0105)
    
    #Grid x-array to be used for convservative re-mapping:
    out_grid = xe.util.grid_2d(-81.725, -79.875, 0.002, 27.375, 28.875, 0.002)
    
elif station=='KDOX':
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

w_dir = "/Users/dpmoore2927/Desktop/TestNC/"

###########################################################################
###########################################################################
#                               READ IN DATA
###########################################################################
###########################################################################

#Read in data from unidata server - option to loop through.
#Eventually write this is a function to be called every day of interest 
#to analyze for sea breeze and spit out information.

# rs=RadarServer('http://149.165.168.53/thredds/radarServer/nexrad/level2/S3/')
rs=''

date = datetime(2018,7,2,22,00)

sb_found_bool, sbf = sb_day_analysis(rs, date, station)

