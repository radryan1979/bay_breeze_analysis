
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
import pandas as pd
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator


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
                
            # if ref[i][j]>conv_thresh:
            #     pass
            # else:
            #     ref[i][j] = np.nan
            
            if ref[i][j]>8.0:
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


###########################################################################
###########################################################################
#                               READ IN DATA
###########################################################################
###########################################################################

#Read in data from unidata server - option to loop through.
#Eventually write this is a function to be called every day of interest 
#to analyze for sea breeze and spit out information.

#Boundaries for desired values
boundinglat = (38.15, 39.65)
boundinglon = (-76.45, -74.6) 
stations=np.array(['DBRG'])

w_dir = "/Users/dpm/Desktop/"

rs=RadarServer('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')
cat = TDSCatalog('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/catalog.xml')
query=rs.query()
date1=datetime(2017,7,4,16,45)
date2=datetime(2017,7,4,19,0)
query.stations('KDOX').time(date1)
rs.validate_query(query)
catalog = rs.get_catalog(query)
dataset = list(catalog.datasets.values())[0]

# now that we have the data, let's go ahead and open it with pyart
radar = pyart.io.read_nexrad_cdm(dataset.access_urls['OPENDAP'])
lats = radar.gate_latitude
lons = radar.gate_longitude
min_lon = lons['data'].min() + 2.5
min_lat = lats['data'].min() + 2.0
max_lat = lats['data'].max() - 2.0
max_lon = lons['data'].max() - 2.5
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
#                       CONVECTIVE ACTIVITY
###########################################################################
###########################################################################


query.stations('KDOX').time(date2)
rs.validate_query(query)
catalog2 = rs.get_catalog(query)
dataset2 = list(catalog2.datasets.values())[0]

radar2 = pyart.io.read_nexrad_cdm(dataset2.access_urls['OPENDAP'])
lats = radar.gate_latitude
lons = radar.gate_longitude
min_lon = lons['data'].min() + 2.5
min_lat = lats['data'].min() + 2.0
max_lat = lats['data'].max() - 2.0
max_lon = lons['data'].max() - 2.5

###########################################################################
###########################################################################
#                   ASSIGNING LAT/LON - POLAR GRID
###########################################################################
###########################################################################
    
lons2, lats2 = regrid_to_cartesian(radar2, lon0, lat0)

###########################################################################
###########################################################################
#                           HANDLING REF DATA
###########################################################################
###########################################################################


#Organizing reflectivity data
ref2 = radar2.get_field(0, 'reflectivity')
unmasked2 = ma.getdata(ref2)  

#Utilizing function at beginning to trim the ref data ----not sure if needed.
my_ref2 = trim_rad_data(lats2, lons2, unmasked2, boundinglat, boundinglon)

###########################################################################
###########################################################################
#                       RE-GRIDDING DATA TO RECTILINEAR GRID
###########################################################################
###########################################################################

gref2, glon, glat = grid_to_rectilinear(my_ref2,lons2,lats2,nlon,nlat)


#######################################
#######################################
    #READING TIME-SERIES
#######################################
#######################################

df = pd.read_csv('/Users/dpm/Desktop/7_4_17CaseStudy.csv', 
    header=0,index_col=0,parse_dates=True)
df = df.rename(columns={"AirTemperature(degC)":"Temp (\N{DEGREE SIGN}C)",
                        "WindSpeed(ms-1)":"Wind Speed (m s-1)",
                        "WindDirection(deg)":"Wind Direction (\N{DEGREE SIGN})"
                        })



#######################################
#######################################
#PLOTTING
#######################################
#######################################
fig = plt.figure(figsize=[15, 10], dpi=100)
ax1 = plt.subplot(121,projection = ccrs.PlateCarree())

r = ax1.pcolormesh(glon, glat, gref, cmap=pyart.graph.cm.NWSRef,
                    norm=mpl.colors.Normalize(vmin=0,vmax=60))

states = cartopy.feature.NaturalEarthFeature(category='cultural',
                              name='admin_1_states_provinces_lines',
                              scale='10m', facecolor='none')
coast = cartopy.feature.NaturalEarthFeature(category='physical', scale='10m',
                            facecolor='none', name='coastline')
                            
reader = shpreader.Reader('/Users/dpm/Documents/Python/countyl010g_shp_nt00964/countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())

ax1.add_feature(COUNTIES, facecolor='none', edgecolor='lightgray')                            
ax1.add_feature(states, linestyle='-', edgecolor='darkslategray',linewidth=1)
ax1.add_feature(coast, linestyle='-', edgecolor='darkslategray',linewidth=1)

ax1.set_title('KDOX NEXRAD {0}'.format(date1))
ax1.set(aspect=1,
        xlim=(glon[0][0],glon[0][-1]),
        ylim=(glat[0][1],glat[-1][0]))

for st in stations:
    plt.plot(-75.59,38.72, '*',
    markersize=16, markerfacecolor=(0,1,1,1),
    markeredgewidth=1, markeredgecolor='k', label='DBRG')

plt.plot(-75.44,38.88, '*',\
markersize=16, markerfacecolor=(0,1,0,1),\
markeredgewidth=1, markeredgecolor='k', label='DMIL')

ax1.legend(loc='lower right')

plt.xticks(np.arange(-74.75,-76.50,-0.25))
plt.yticks(np.arange(38.25,39.75,0.25))

ax2 = plt.subplot(122,projection = ccrs.PlateCarree())

r2 = ax2.pcolormesh(glon, glat, gref2, cmap=pyart.graph.cm.NWSRef,
                    norm=mpl.colors.Normalize(vmin=0,vmax=60))

ax2.add_feature(COUNTIES, facecolor='none', edgecolor='lightgray')                            
ax2.add_feature(states, linestyle='-', edgecolor='darkslategray',linewidth=1)
ax2.add_feature(coast, linestyle='-', edgecolor='darkslategray',linewidth=1)

ax2.set_title('KDOX NEXRAD {0}'.format(date2))
ax2.set(aspect=1,
        xlim=(glon[0][0],glon[0][-1]),
        ylim=(glat[0][1],glat[-1][0]))

for st in stations:
    plt.plot(-75.59,38.72, '*',
    markersize=16, markerfacecolor=(0,1,1,1),
    markeredgewidth=1, markeredgecolor='k', label='DBRG')

plt.plot(-75.44,38.88, '*',\
markersize=16, markerfacecolor=(0,1,0,1),\
markeredgewidth=1, markeredgecolor='k', label='DMIL')
    
ax2.legend(loc='lower right')

plt.xticks(np.arange(-74.75,-76.50,-0.25))
plt.yticks(np.arange(38.25,39.75,0.25))

fig.subplots_adjust(right=0.8)

cbar_ax = fig.add_axes([0.85,0.2,0.02,0.6])
cbar = fig.colorbar(r2, orientation='vertical', cax=cbar_ax)
cbar.set_label('Equivalent Reflectivity Factor (dBZ)', rotation=90)

# ax_inset = fig.add_axes([0.115,0.26,0.14,0.14], projection=ccrs.PlateCarree())
# ax_inset.set_global()
# ax_inset.add_feature(COUNTIES, facecolor='none', edgecolor='lightgray')                            
# ax_inset.add_feature(states, linestyle='-', edgecolor='darkslategray',linewidth=1)
# ax_inset.add_feature(coast, linestyle='-', edgecolor='darkslategray',linewidth=1)
# ax_inset.add_patch(mpatches.Rectangle(xy=[-76.4, 38.2], width=1.8, height=1.45,
#                 facecolor='None', edgecolor='red'))
# ax_inset.set(aspect=1,
#         xlim = (-83.0, -68.0),
#         ylim = (32.0, 45.0))

fig.suptitle('July 4, 2017 Sea Breeze',fontsize=24,y=0.85)

plt.savefig('/Users/dpm/Desktop/7_4_17RadCS.png', transparent=True)


fig2 = plt.figure(figsize=[15, 9], dpi=100)
majloc = mpl.ticker.MultipleLocator(base=60.0)
minloc = mpl.ticker.MultipleLocator(base=30.0)
fmt = mpl.dates.DateFormatter('%h:%m')
my_ax1 = plt.subplot(311)
df['Temp (\N{DEGREE SIGN}C)'].plot(color='r')
# my_ax1.set_ylabel('Temperature (\N{DEGREE SIGN}C)',rotation='vertical')
my_ax1.set_xlabel('Time in UTC',fontsize=24)
my_ax1.xaxis.set_major_locator(majloc)
my_ax1.xaxis.set_minor_locator(minloc)
my_ax2 = plt.subplot(312)
df['Wind Speed (m s-1)'].plot(color='g',sharex=my_ax1,)
# my_ax2.set_ylabel('Wind Speed (m s-1)',rotation='vertical')
# my_ax2.set_xlabel('Time in UTC',fontsize=24)
my_ax2.xaxis.set_major_locator(majloc)
my_ax2.xaxis.set_minor_locator(minloc)
my_ax3 = plt.subplot(313)
df['Wind Direction (\N{DEGREE SIGN})'].plot(sharex=my_ax1)
my_ax3.set_yticks((0,90,180,270,360))
labels = [item.get_text() for item in my_ax3.get_yticklabels()]
labels[0] = 'N'
labels[1] = 'E'
labels[2] = 'S'
labels[3] = 'W'
labels[4] = 'N'
my_ax3.set_yticklabels(labels)
# my_ax3.set_ylabel('Wind Direction (\N{DEGREE SIGN})',rotation='vertical')
my_ax3.set_xlabel('Time in UTC',fontsize=16)
my_ax3.xaxis.set_major_locator(majloc)
my_ax3.xaxis.set_minor_locator(minloc)

fig2.patches.extend([plt.Rectangle((0.415,0.099),0.02,0.694,
                                  fill=True, color='r', alpha=0.4, zorder=1000,
                                  transform=fig.transFigure, figure=fig)])
legend_elements = [mpatches.Patch(facecolor='red', edgecolor='None',
                         alpha=0.4, label='Sea Breeze Passage')]
fig2.legend(handles=legend_elements,loc=(0.77,0.115))
                                  

plt.savefig('/Users/dpm/Desktop/7_4_17StationCS.png')



