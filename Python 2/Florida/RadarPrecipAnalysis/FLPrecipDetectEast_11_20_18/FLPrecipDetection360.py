"""
Coded by: Dan Moore

This program will ingest a list of dates from the Hughes
detection algorithm in Florida, determine which radar to
use, then attempt to detect convective activity.

We will define convective activity as the presence of a
cluster of reflectivity values of mean>40dBZ remaining in
view of radar for at least 30 mins. Once these criteria are
met, the program will exit with 'success' for that date, 
or if no convective activity was found on the given date, 
a 'failure' will be returned.

The output will be the original list of dates where sea breeze
was detected, with a 'True' or 'False' associated with it.

Updated: 11-26-18
"""

from siphon.radarserver import RadarServer
from datetime import datetime,timedelta
from siphon.cdmr import Dataset
from scipy import spatial
from scipy.interpolate import griddata
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyart
import numpy.ma as ma
import pyproj
import cartopy.crs as ccrs
import pandas as pd
import sys
from time import sleep

def detect_convect(rad, date):
    
    #Number of regions within radar scope. Miami and tampa
    #have 3 regions, melbourne has 2.
    if rad == melbourne:
        num_reg = 2
    else:
        num_reg = 3
        
    #Keep track of number of scans actually evaluated.
    scan_num = 0
        
    #Initialize temp dataframe to house bins.
    #This will be sent back to calling program.
    df_t = pd.DataFrame(columns=['Date','Region','Time(UTC)','40_45dBZ','45_50dBZ',\
                                '50_55dBZ','55_60dBZ','60_65dBZ','65_70dBZ',\
                                '70_75dBZ','75_80dBZ','80+dBZ'])
    
    #Initialize convective boolean:
    convect_found = False

    #For analyzing until 23:59UTC ~sundown.
    nextday = date + timedelta(days=1)
    
    #Usually will be at hours=13 to start around sunrise.
    dt = date + timedelta(hours=13)
    
    #Define radar server location:
    rs=RadarServer('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')
    
    query=rs.query()
    query.stations(rad).time_range(dt,nextday)
    try:
        cat = rs.get_catalog(query)
    except ConnectionError:
        print("Connection Refused, trying again.")
        #Sleep for two minutes if this happens
        sleep(300)
        try:
            cat = rs.get_catalog(query)
        except:
            print("Day not cooperating. Day Skipped.")
            error = "Connection Error."
            return None, error, None
    ds = list(cat.datasets.values())
    
    del rs, query, cat

    #If no scans are found for the given date:
    if len(ds)==0:
        convect_found = False
        error = "No data found for selected date."
        time = np.nan
        
        return convect_found, error, df_t
    
    #Establish location of radar
    loc = pyart.io.nexrad_common.get_nexrad_location(rad)
    lon0 = loc[1] ; lat0 = loc[0]
    
    for numrad in range(0,len(ds),2):
        
        try: 
            #Create radar object
            radar = pyart.io.read_nexrad_cdm(ds[numrad].access_urls['OPENDAP'])
        except:
            print("Line 84: read_nexrad_cdm error. Please check NEXRAD server.")
            continue
        
        #Calculate timestamp to decrease computing time.
        time = radar.time['units'].split(' ')[-1].split('T')
        time = time[0] + ' ' + time[1][:-1]
        timestamp = datetime.strptime(time, '%Y-%m-%d %H:%M:%S')
        
        if scan_num>0:
            timediff = (timestamp-prevtime).seconds/60.0/60.0
        
        #If next iteration is within 15 mins of previous iteration, this scan is skipped.
        if scan_num>0 and timediff<0.25:
            continue
        
        prevtime = timestamp
        scan_num += 1
        
        #Regrid to cartesian coordinates:
        lons, lats = regrid_to_cartesian(radar, lon0, lat0)
        
        #Organizing reflectivity data
        ref = radar.get_field(0, 'reflectivity')
        unmasked = ma.getdata(ref)
        
        my_ref = trim_rad_data(lats, lons, unmasked, boundinglat, boundinglon)
        
        #If no convective pixels are found, skip scan.
        if np.count_nonzero(~np.isnan(my_ref)) == 0:
            #Fill array to show that scan was analyzed
            df_temp = rad_to_bins_fill_zero(rad,date,time)
            df_t = df_t.append(df_temp, sort=True)
            continue
        else:
            convect_found=True

        #Re-gridding to rectilinear grid
        nlon = 500; nlat = 500 #can be changed, but seems good.
        gref, glon, glat, gderror = grid_to_rectilinear(my_ref,lons,lats,nlon,nlat)
        
        if gderror:
            print("griddata error, scan skipped.")
            continue
        
        df_temp = rad_to_bins(gref,rad,date,time)
        
        df_t = df_t.append(df_temp, sort=True)


    #Eliminate days with insufficient data (total of ~45 scans per day
    #analyzed by this program).
    if numrad<20:
        convect_found = False
        error = "Insufficient data: Not enough scans to analyze this day."
        
        return convect_found, error, df_t

    elif not convect_found:
        error = "Sufficient data, but convective activity not detected."

        return convect_found, error, df_t
    
    else:
        error = "All good"
        
        return convect_found, error, df_t


#Function to trim radar data to desired bounding rectangle:
#Required are 1D arrays of latitudes, longitudes, and reflectivities
#and two tuples of latitudes and longitudes defining the boundary of 
#the desired boundary.
def trim_rad_data(lats, lons, ref, boundinglat=0, boundinglon=0):
    for i in range(ref.shape[0]-1):
        for j in range(ref.shape[1]-1):
            if (lats[i][j]>boundinglat[0] and lats[i][j]<boundinglat[1] and \
                lons[i][j]>boundinglon[0] and lons[i][j]<boundinglon[1]):
                pass
            else:
                # lats[i] = np.nan; lons[i] = np.nan
                ref[i][j] = np.nan
                
            if ref[i][j]>=conv_thresh:
                pass
            else:
                ref[i][j] = np.nan
    return ref

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
    
    try:
        # Interpolate data onto grid using linear interpolation
        gref = griddata((rav_lons,rav_lats),rav_ref,(glon,glat),method='linear')
        gderror=False
    except:
        gref = np.nan
        gderror=True
    
    return gref, glon, glat,gderror


def bin_refs(gref,mask):
    
    bin1=0; bin2=0; bin3=0; bin4=0
    bin5=0; bin6=0; bin7=0; bin8=0; bin9=0
    
    #Mask array by multiplying by zero/one out/in mask.
    temp_ref = gref * mask
    
    if np.count_nonzero(~np.isnan(temp_ref)) == 0:
        ref_bins = np.array((bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9))
        
        return ref_bins
    
    for i in range(np.shape(temp_ref)[0]):
        for j in range(np.shape(temp_ref)[1]):
            if temp_ref[i][j]==0:
                continue
            elif temp_ref[i][j]>=40 and temp_ref[i][j]<45:
                bin1+=1
            elif temp_ref[i][j]>=45 and temp_ref[i][j]<50:
                bin2+=1
            elif temp_ref[i][j]>=50 and temp_ref[i][j]<55:
                bin3+=1
            elif temp_ref[i][j]>=55 and temp_ref[i][j]<60:
                bin4+=1
            elif temp_ref[i][j]>=60 and temp_ref[i][j]<65:
                bin5+=1
            elif temp_ref[i][j]>=65 and temp_ref[i][j]<70:
                bin6+=1
            elif temp_ref[i][j]>=70 and temp_ref[i][j]<75:
                bin7+=1
            elif temp_ref[i][j]>=75 and temp_ref[i][j]<80:
                bin8+=1
            elif temp_ref[i][j]>=80:
                bin9+=1
    
    ref_bins = np.array((bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9))
    
    return ref_bins
    
def rad_to_bins(gref, rad, date, time):
    
    datestr = datetime.strftime(date, '%d-%b-%Y')

    if rad == miami:
        region_names = np.array(('miami_east_1','miami_east_2','miami_west_1'))
        data_array = np.empty((3,9))
        
        dates = []; times = []; regions= []
        for i in range(np.shape(miami_mask)[2]):
            temp_bins = bin_refs(gref, miami_mask[:,:,i])
            dates.append(datestr)
            times.append(time)
            regions.append(region_names[i])
            data_array[i,0] = temp_bins[0]
            data_array[i,1] = temp_bins[1]
            data_array[i,2] = temp_bins[2]
            data_array[i,3] = temp_bins[3]
            data_array[i,4] = temp_bins[4]
            data_array[i,5] = temp_bins[5]
            data_array[i,6] = temp_bins[6]
            data_array[i,7] = temp_bins[7]
            data_array[i,8] = temp_bins[8]
    
    elif rad == melbourne:
        region_names = np.array(('mlb_east_1','mlb_east_2'))
        data_array = np.empty((2,9))
        
        dates = []; times = []; regions= []
        for i in range(np.shape(mlb_mask)[2]):
            temp_bins = bin_refs(gref, mlb_mask[:,:,i])
            dates.append(datestr)
            times.append(time)
            regions.append(region_names[i])
            data_array[i,0] = temp_bins[0]
            data_array[i,1] = temp_bins[1]
            data_array[i,2] = temp_bins[2]
            data_array[i,3] = temp_bins[3]
            data_array[i,4] = temp_bins[4]
            data_array[i,5] = temp_bins[5]
            data_array[i,6] = temp_bins[6]
            data_array[i,7] = temp_bins[7]
            data_array[i,8] = temp_bins[8]
        
    elif rad == tampa:
        region_names = np.array(('tampa_west_1','tampa_west_2','tampa_west_3'))
        data_array = np.empty((3,9))
        
        dates = []; times = []; regions= []
        for i in range(np.shape(miami_mask)[2]):
            temp_bins = bin_refs(gref, tampa_mask[:,:,i])
            dates.append(datestr)
            times.append(time)
            regions.append(region_names[i])
            data_array[i,0] = temp_bins[0]
            data_array[i,1] = temp_bins[1]
            data_array[i,2] = temp_bins[2]
            data_array[i,3] = temp_bins[3]
            data_array[i,4] = temp_bins[4]
            data_array[i,5] = temp_bins[5]
            data_array[i,6] = temp_bins[6]
            data_array[i,7] = temp_bins[7]
            data_array[i,8] = temp_bins[8]
    
    df_t1 = pd.DataFrame({'Date': dates,
                       'Region': regions,
                       'Time(UTC)': times,
                       '40_45dBZ': data_array[:,0],
                       '45_50dBZ': data_array[:,1],
                       '50_55dBZ': data_array[:,2],
                       '55_60dBZ': data_array[:,3],
                       '60_65dBZ': data_array[:,4],
                       '65_70dBZ': data_array[:,5],
                       '70_75dBZ': data_array[:,6],
                       '75_80dBZ': data_array[:,7],
                       '80+dBZ':   data_array[:,8]
                       })
    
    return df_t1

#When scan has no convective activity, fill temp array
#with zeros so that we know the scan time was analyzed.
def rad_to_bins_fill_zero(rad, date, time):
    
    datestr = datetime.strftime(date, '%d-%b-%Y')

    if rad == miami:
        region_names = np.array(('miami_east_1','miami_east_2','miami_west_1'))
        dates = []; times = []; regions= []
        for i in range(np.shape(miami_mask)[2]):
            dates.append(datestr)
            times.append(time)
            regions.append(region_names[i])
    
    elif rad == melbourne:
        region_names = np.array(('mlb_east_1','mlb_east_2'))
        
        dates = []; times = []; regions= []
        for i in range(np.shape(mlb_mask)[2]):
            dates.append(datestr)
            times.append(time)
            regions.append(region_names[i])
        
    elif rad == tampa:
        region_names = np.array(('tampa_west_1','tampa_west_2','tampa_west_3'))
        
        dates = []; times = []; regions= []
        for i in range(np.shape(tampa_mask)[2]):
            dates.append(datestr)
            times.append(time)
            regions.append(region_names[i])
    
    df_t1 = pd.DataFrame({'Date': dates,
                       'Region': regions,
                       'Time(UTC)': times,
                       '40_45dBZ': 0,
                       '45_50dBZ': 0,
                       '50_55dBZ': 0,
                       '55_60dBZ': 0,
                       '60_65dBZ': 0,
                       '65_70dBZ': 0,
                       '70_75dBZ': 0,
                       '75_80dBZ': 0,
                       '80+dBZ':   0
                       })
    
    return df_t1


#Initializing constants:
conv_thresh = 40.0#dBZ
miami = 'KAMX' 
miami_st = np.array([410, 420, 425, 440, 450])
tampa = 'KTBW'
tampa_st = np.array([350, 360, 380, 480, 490])
melbourne = 'KMLB'
melbourne_st = np.array([340, 371, 435])


#Create mask arrays for each region
miami_mask = np.zeros((500,500,3))
#Set east_1 mask
for i in range(np.shape(miami_mask)[0]):
    for j in range(np.shape(miami_mask)[1]):
        if ((i>=205 and i<499) and (j>=382 and j<455)) or \
            ((i>=100 and i<205) and (j>=352 and j<425)) or \
            ((i>=50 and i<100) and (j>=295 and j<390)):
                miami_mask[i,j,0] = 1.0

#Set east_2 mask            
for i in range(np.shape(miami_mask)[0]):
    for j in range(np.shape(miami_mask)[1]):
        if ((i>=205 and i<499) and (j>=290 and j<382)) or \
            ((i>=100 and i<205) and (j>=290 and j<352)):
                miami_mask[i,j,1] = 1.0
#Set west_1 mask
    #bounding points (i,j):NW(460,10),SW(25,266),NE(460,105),SE(25,361)
#Slope of lines
m1 = -435.0/256.0
#"y-intercept" for each bounding line
b_west1 = 460.0 - m1 * 10.0
b_east1 = 460.0 - m1 * 105.0
for i in range(np.shape(miami_mask)[0]):
    for j in range(np.shape(miami_mask)[1]):
        if ((i>=25 and i<460) and (j>((i-b_west1)/m1) and \
            j<((i-b_east1)/m1))):
                miami_mask[i,j,2] = 1.0

#Set Tampa mask (West_1, West_2, West_3)
tampa_mask = np.zeros((500,500,3))
    #West_1
    #bounding points (i,j):NW(420,140),NE(420,210),SW(100,250),SE(100,320)
m7 = -320.0/100.0
b_west7 = 420.0 - m7 * 140.0
b_east7 = 420.0 - m7 * 210.0
for i in range(np.shape(tampa_mask)[0]):
    for j in range(np.shape(tampa_mask)[1]):
        if ((i>=100 and i<=420) and (j>((i-b_west7)/m7) and \
            j<((i-b_east7)/m7))):
                tampa_mask[i,j,0] = 1.0
        
    #West_2
    #bounding points (i,j):NW(420,210),NE(420,280),SW(100,320),SE(100,390)
b_west8 = 420.0 - m7 * 210.0
b_east8 = 420.0 - m7 * 280.0
for i in range(np.shape(tampa_mask)[0]):
    for j in range(np.shape(tampa_mask)[1]):
        if ((i>=100 and i<=420) and (j>((i-b_west8)/m7) and \
            j<((i-b_east8)/m7))):
                tampa_mask[i,j,1] = 1.0
        
    #West_3
    #bounding points (i,j):NW(420,280),NE(420,350),SW(100,390),SE(100,460)
b_west9 = 420.0 - m7 * 280.0
b_east9 = 420.0 - m7 * 350.0
for i in range(np.shape(tampa_mask)[0]):
    for j in range(np.shape(tampa_mask)[1]):
        if ((i>=100 and i<=420) and (j>((i-b_west9)/m7) and \
            j<((i-b_east9)/m7))):
                tampa_mask[i,j,2] = 1.0
        
#Set Melbourne mask (East_1, East_2)
#Create mask array
mlb_mask = np.zeros((500,500,2))
    #East_1
    #bounding points (i,j):NW(480,154),NE(480,221),SW(24,361),SE(24,428)
m3 = -456.0/207.0
b_west5 = 480.0 - m3 * 154.0
b_east5 = 480.0 - m3 * 221.0
for i in range(np.shape(mlb_mask)[0]):
    for j in range(np.shape(mlb_mask)[1]):
        if  ((i>=24 and i<=480) and (j>((i-b_west5)/m3) and \
            j<((i-b_east5)/m3))):
                mlb_mask[i,j,0] = 1.0

    #East_2
    #bounding points (i,j):NW(480,87),NE(480,154),SW(24,294),SE(24,361)     
b_east6 = 480.0 - m3 * 154.0
b_west6 = 480.0 - m3 * 87.0
for i in range(np.shape(mlb_mask)[0]):
    for j in range(np.shape(mlb_mask)[1]):
        if  ((i>=24 and i<=480) and (j>((i-b_west6)/m3) and \
            j<((i-b_east6)/m3))):
                mlb_mask[i,j,1] = 1.0        



#Initialize exmpty arrays:
convection = []
error_log = []
miami_dates = []
tampa_dates = []
melbourne_dates = []

#Initialize dataframe for reflectivity bins:
mo_df = pd.DataFrame(columns=['Date','Region','Time(UTC)','40_45dBZ','45_50dBZ',\
                                '50_55dBZ','55_60dBZ','60_65dBZ','65_70dBZ',\
                                '70_75dBZ','75_80dBZ','80+dBZ'])

#Initialize dataframe for reflectivity bins:
ref_df = pd.DataFrame(columns=['Date','Region','Time(UTC)','40_45dBZ','45_50dBZ',\
                                '50_55dBZ','55_60dBZ','60_65dBZ','65_70dBZ',\
                                '70_75dBZ','75_80dBZ','80+dBZ'])
                                
station = '360'
coast = 'West'

workingdir = "/home/1898/FLPrecip/2013_2017{0}/{1}/".format(coast,station)

#Ingest date file:
df = pd.read_csv(workingdir+'Hughes{0}Times{1}.csv'.format(coast,station))

date = datetime.strptime(df['Date'][0], '%d-%b-%Y')

mo = date.month
yr = date.year

for i in df.index:

    #Determine if date has already been processed by 
    #the given radar. If so, continue. If not, set radar
    #and boundaries for processing.
    if df['Station'][i] in miami_st:
        if df['Date'][i] in miami_dates:
            continue
        rad_st = miami
        miami_dates.append(df['Date'][i])
        boundinglon = (-81.9, -80.0)
        boundinglat = (25.1, 27.0)
        
    elif df['Station'][i] in tampa_st:
        if df['Date'][i] in tampa_dates:
            continue
        rad_st = tampa
        tampa_dates.append(df['Date'][i])
        boundinglon = (-83.2, -81.2)
        boundinglat = (24.5, 26.5)
    elif df['Station'][i] in melbourne_st:
        if df['Date'][i] in melbourne_dates:
            continue
        rad_st = melbourne
        melbourne_dates.append(df['Date'][i])
        boundinglon = (-82.0, -80.0)
        boundinglat = (27.0, 29.0)
        
    else:
        print("Station not in range. {0} not processed.".format(df['Date'][i]))
        continue

    date = datetime.strptime(df['Date'][i], '%d-%b-%Y')
    
    print("Processing {0}, station: {1}, at {2}.".format(\
    df['Station'][i],df['Date'][i],datetime.now()))
    
    if date.month == mo and date.year == yr:
        pass
    else:
        #If new month, write monthly dataframe
        mo_df.to_csv(workingdir+"{0}_{1}_RefBins.csv".format(mo,yr))
        
        mo = date.month
        
        yr = date.year
        
        ref_df = ref_df.append(mo_df,sort=True)
        
        #Initialize dataframe for reflectivity bins:
        mo_df = pd.DataFrame(columns=['Date','Region','Time(UTC)','40_45dBZ','45_50dBZ',\
                                        '50_55dBZ','55_60dBZ','60_65dBZ','65_70dBZ',\
                                        '70_75dBZ','75_80dBZ','80+dBZ'])
    
    convect_bool, error, df_temp = detect_convect(rad_st, date)
    
    convection.append(convect_bool); error_log.append(error)
        
    #Append temp dataframe to monthly data frame
    mo_df = mo_df.append(df_temp,sort=True)
    
    #Delete temperary dataframe
    del df_temp

mo_df.to_csv(workingdir+"{0}_{1}_RefBins.csv".format(mo,yr))

ref_df = ref_df.append(mo_df,sort=True)
                                
ref_df.to_csv(workingdir+"RefBins.csv")

