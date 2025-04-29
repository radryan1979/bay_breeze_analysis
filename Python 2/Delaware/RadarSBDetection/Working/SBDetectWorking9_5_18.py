"""
Coded by Dan Moore
9/5/18
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
    
###########################################################################
###########################################################################
#                           INITIALIZATIONS
###########################################################################
###########################################################################


#Boundaries for desired values
boundinglat = (38.15, 39.65)
boundinglon = (-76.45, -74.6)

#Defining the DE Coastline
cstlnlat = np.array((39.15, 39.125, 39.1, 39.075, 39.05, 39.025, 39.0, #Bay Coastline
                    38.975, 38.95, 38.925, 38.9, 38.875, 38.85, 38.825, 
                    38.8, 38.775,
                    38.75, 38.725, 38.7, 38.675, 38.65, 38.625, 38.6,#Ocean Coastline
                    38.575, 38.55, 38.525, 38.5, 38.475, 38.45)).T
                    
cstlnlon = np.array((-75.41, -75.41, -75.40, -75.40, -75.39, -75.35, -75.33,#Bay Coastline
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

#ALSO HAVE TO MAKE SURE ONLY EVERY 15-30 MINS IS USED TO PRESERVE COMPUTING.
rs=RadarServer('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')
query=rs.query()
dt=datetime(2018, 6, 1, 14)
query.stations('KDOX').time_range(dt,dt+timedelta(hours=1))
cat = rs.get_catalog(query)
ds = list(cat.datasets.values())[0]

# for numrad in range(len(ds)):
radar = pyart.io.read_nexrad_cdm(ds.access_urls['OPENDAP'])

#Calculate timestamp to decrease computing time.
# timestamp = radar.time['units'].split(' ')[-1].split('T')
# timestamp = timestamp[0] + ' ' + timestamp[1][:-1]
# timestamp = datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S')
# #If next iteration is within 15 mins of previous iteration, this scan is skipped.
# if numrad>0 and ((timestamp-prevtime).seconds/60.0)<15:
#     continue
# prevtime = timestamp


#Alternative if file is local
# radar = pyart.io.read_cfradial("/Users/dpmoore2927/Desktop/20180601cfrad")

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
    
# Grid data into geographic coordinates
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
    

###########################################################################
###########################################################################
#                           HANDLING REF DATA
###########################################################################
###########################################################################


#Organizing reflectivity data
ref = radar.get_field(0, 'reflectivity')
unmasked = ma.getdata(ref)  
#Otherwise, lats and lons have 1 extra data point.
lats = lats[0:720, 0:1832]
lons = lons[0:720, 0:1832]
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

rav_lats = lats.ravel()
rav_lons = lons.ravel()
rav_ref = my_ref.ravel()

#Grid Data using matplotlib
nlon = 1000; nlat = 1000; #can be changed, but seems good.
grid_lons = np.linspace(boundinglon[0],boundinglon[1],nlon)
grid_lats = np.linspace(boundinglat[0],boundinglat[1],nlat)
glon,glat = np.meshgrid(grid_lons,grid_lats)

# Interpolate data onto grid using linear interpolation
gref = griddata((rav_lons,rav_lats),rav_ref,(glon,glat),method='linear')


###########################################################################
###########################################################################
#                           ANALYZE DATA
###########################################################################
###########################################################################

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
    tlon0 = cstlnlon[i]+0.01
    
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
    
#Create tuple of coordinates    
sbf = list(zip(sbflon,sbflat))



###########################################################################
###########################################################################
#                           PLOTTING DATA
###########################################################################
###########################################################################


# #Parameters for simple plotting
# newproj = ccrs.LambertConformal(central_latitude=lat0,central_longitude=lon0)
# cmap = 'pyart_NWSRef'
# levs = np.linspace(0,80,41,endpoint=True)
# ticks = np.linspace(0,80,9,endpoint=True)
# label = 'Radar Reflectivity Factor ($\mathsf{dBZ}$)'
# #Normalize the colormap based on the levels provided above
# norm = mpl.colors.BoundaryNorm(levs,256)
# #Features on the map------not currently working
# # political_boundaries = cfeature.NaturalEarthFeature(category='cultural',
# #                               name='admin_0_boundary_lines_land',
# #                               scale='10m', facecolor='none')
# # states = cfeature.NaturalEarthFeature(category='cultural',
# #                               name='admin_1_states_provinces_lines',
# #                               scale='10m', facecolor='none')

# #Making a binary array for imaging, if necessary
# new_ref = ma.masked_where(np.isnan(my_ref),my_ref)
# bin_ref = (my_ref>0).astype(int)
# plt_ref = ma.masked_where(bin_ref==0,bin_ref)

#Plotting rectilinear grid
# figrect = plt.figure(figsize=[12,12], dpi=300)
# plt.pcolormesh(glon,glat,gref,norm=norm,cmap=cmap)
# plt.axis([boundinglon[0], boundinglon[1], boundinglat[0], boundinglat[1]])
# for pt in range(len(sbf)):
#     plt.plot(sbf[pt][0],sbf[pt][1],'rp',markersize=10)
# plt.show()


# #Plotting data in polar grid
# # figpol = plt.figure(figsize=[12, 12], dpi=300)
# # ax = fig.add_subplot(1,1,1, projection=newproj)
# # im = ax.pcolormesh(lons, lats, plt_ref)#, norm=norm, cmap=cmap)
# # ax.set_extent((boundinglon[0],boundinglon[1],boundinglat[0],boundinglat[1]),crs=newproj)
# #Adding boundaries, etc. ----Not currently working
# # ax.add_feature(states)
# # ax.add_feature(political_boundaries)


# #Plotting points to test locations, if needed.
# # plt.plot(lons[360][100],lats[360][100],'bp',markersize=10)
# # plt.plot(lons[360][150],lats[360][150],'bp',markersize=10)

# # plt.plot(lons[300][60],lats[300][60],'gp',markersize=10)
# # plt.plot(lons[300][92],lats[300][92],'gp',markersize=10)



# #Saving figures
# figrect.savefig("/Users/dpmoore2927/Desktop/gridfigure.png", dpi=300, bbox_inches='tight')
# # figpol.savefig("/Users/dpmoore2927/Desktop/ogfigure.png", dpi=300, bbox_inches='tight')




#Notification that the program has finished
print('done')