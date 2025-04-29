import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from siphon.catalog import TDSCatalog, get_latest_access_url
from netCDF4 import Dataset
import matplotlib as mpl
import time
import os.path
import numpy as np
from datetime import datetime, timedelta
import pyart
from siphon.radarserver import RadarServer, get_radarserver_datasets
import numpy.ma as ma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
from scipy import spatial
from scipy.interpolate import griddata
import sys


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
    
    #Set end of detection
    idxlonend = find_nearest(onedimlon,-75.85)
    
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
        for j in range(idxlon,idxlonend,-1):
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
cat = TDSCatalog('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/catalog.xml')
query=rs.query()
date=datetime(2018,7,2,19)
query.stations('KDOX').time(date)
cat = rs.get_catalog(query)
ds = list(cat.datasets.values())
rs.validate_query(query)
catalog = rs.get_catalog(query)
dataset = list(catalog.datasets.values())[0]

# now that we have the data, let's go ahead and open it with pyart
radar = pyart.io.read_nexrad_cdm(dataset.access_urls['OPENDAP'])
lats = radar.gate_latitude
lons = radar.gate_longitude
min_lon = lons['data'].min() + 3.75
min_lat = lats['data'].min() + 3.0
max_lat = lats['data'].max() - 3.0
max_lon = lons['data'].max() - 4.25
loc = pyart.io.nexrad_common.get_nexrad_location('KDOX')
lon0 = loc[1] ; lat0 = loc[0]


projection = ccrs.Mercator(
                central_longitude=lon0,
                min_latitude=min_lat, max_latitude=max_lat)

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

#Find SBF
sbfound, sbf = sea_breeze_front_coords(gref,glon,glat)

if sbfound:
    sbflon, sbflat = zip(*sbf)
else:
    print("SBF Not Found.")
    sys.exit()

###########################################################################
###########################################################################
#                       GET LINES FOR TRANSECTS
###########################################################################
###########################################################################

#For easier processing:
onedimlon = glon[0,:]

translonE = np.zeros(len(cstlnlon))

#Set end of detection
idxlonend = find_nearest(onedimlon,-75.85)
translonW = onedimlon[idxlonend]

#Loop through transect lines
for i in range(np.size(cstlnlat)):
    tlon0 = cstlnlon[i]+0.01 #Just a bit offshore in case
    
    #Find index of nearest longitude in grid: this is our i value
    idxlon = find_nearest(onedimlon,tlon0)#represents starting column
    
    translonE[i] = onedimlon[idxlon]

#######################################
#######################################
#PLOTTING
#######################################
#######################################
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

for i in range(len(cstlnlat)):
    plt.plot([translonW, translonE[i]], [cstlnlat[i], cstlnlat[i]],
         color='purple',
         linewidth=0.6,
         transform=ccrs.PlateCarree(),
         )

#Plot SBF
for pt in range(len(sbf)):
    plt.plot(sbf[pt][0],sbf[pt][1],'rp',markersize=4)

plt.xticks(np.arange(-74.75,-76.50,-0.25))
plt.yticks(np.arange(38.25,39.75,0.25))

plt.show()



