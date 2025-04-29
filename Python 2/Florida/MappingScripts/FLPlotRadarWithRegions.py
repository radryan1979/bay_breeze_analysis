
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
                
            if ref[i][j]>conv_thresh:
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

# Turn on for Miami Radar:
boundinglon = (-81.9, -80.0)
boundinglat = (25.1, 27.0)  
stations=np.array([410,420,425,440,450])

# # Turn on for Melbourne Radar:
# boundinglon = (-82.0, -80.0)
# boundinglat = (27.0, 29.0) 
# stations=np.array([340,371,435])

# # Turn on for Tampa Radar:
# boundinglon = (-83.2, -81.2)
# boundinglat = (26.5, 28.5) 
# stations=np.array([350,360,380,480,490])

conv_thresh = 40.0#dBZ

# Read in station locations
locs=pd.read_csv("/Users/dpmoore2927/Desktop/UD/Hyperion/FloridaStationLocations.csv",
    header=None)
locs=locs.set_index(0)
                

rs=RadarServer('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')
cat = TDSCatalog('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/catalog.xml')
query=rs.query()
date=datetime(2018,7,2,20)
query.stations('KAMX').time(date)
cat = rs.get_catalog(query)
ds = list(cat.datasets.values())
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
loc = pyart.io.nexrad_common.get_nexrad_location('KAMX')
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

nlon = 500; nlat = 500 #can be changed, but seems good.
gref, glon, glat = grid_to_rectilinear(my_ref,lons,lats,nlon,nlat)

###########################################################################
###########################################################################
#                       MAKING POLYGONS FOR MASKS
###########################################################################
###########################################################################

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
        if ((i>=205 and i<499) and (j>=190 and j<382)) or \
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
        if ((i>=25 and i<460) and (j>((i-b_west1)/m1)) and \
            (j<(i-b_east1)/m1)):
                miami_mask[i,j,2] = 1.0
                
#Masking reflectivities
new_ref = gref*miami_mask[:,:,0]+gref*miami_mask[:,:,1]+gref*miami_mask[:,:,2]


#Create mask arrays for each region
#Set east_1 mask
miami_P1 = np.array([[glon[499][455]   ,    glat[499][455]] ,
                     [glon[499][382]  ,    glat[499][382]],
                     [glon[205][382]  ,    glat[205][382]],
                     [glon[205][352]  ,    glat[205][352]],
                     [glon[100][352]  ,    glat[100][352]],
                     [glon[100][295]  ,    glat[100][295]],
                     [glon[50][295]   ,    glat[50][295]] ,
                     [glon[50][390]   ,    glat[50][390]] ,
                     [glon[100][390]  ,    glat[100][390]],
                     [glon[100][425]   ,    glat[100][425]] ,
                     [glon[205][425]   ,    glat[205][425]] ,
                     [glon[205][455]   ,    glat[205][455]] 
                     ])

#Set east_2 mask            
miami_P2 = np.array([[glon[499][382]   ,    glat[499][382]] ,
                     [glon[499][290]   ,    glat[499][290]] ,
                     [glon[100][290]   ,    glat[100][290]] ,
                     [glon[100][352]   ,    glat[100][352]] ,
                     [glon[205][352]   ,    glat[205][352]] ,
                     [glon[205][382]   ,    glat[205][382]] 
                     ])

#bounding points (i,j):NW(460,10),SW(25,266),NE(460,101),SE(25,361)
#Set west_1 mask
miami_P3 = np.array([[glon[460][10]   ,    glat[460][10]] ,
                     [glon[460][105]   ,    glat[460][105]] ,
                     [glon[25][361]    ,    glat[25][361]]  ,
                     [glon[25][266]    ,    glat[25][266]] 
                     ])




# #Set Melbourne mask (East_1, East_2)
# #Create mask array
# mlb_mask = np.zeros((500,500,2))
#     #East_1
#     #bounding points (i,j):NW(480,154),NE(480,221),SW(24,361),SE(24,428)
# m3 = -456.0/207.0
# b_west5 = 480.0 - m3 * 154.0
# b_east5 = 480.0 - m3 * 221.0
# for i in range(np.shape(mlb_mask)[0]):
#     for j in range(np.shape(mlb_mask)[1]):
#         if  ((i>=24 and i<=480) and (j>((i-b_west5)/m3) and \
#             j<((i-b_east5)/m3))):
#                 mlb_mask[i,j,0] = 1.0

#     #East_2
#     #bounding points (i,j):NW(480,87),NE(480,154),SW(24,294),SE(24,361)     
# b_east6 = 480.0 - m3 * 154.0
# b_west6 = 480.0 - m3 * 87.0
# for i in range(np.shape(mlb_mask)[0]):
#     for j in range(np.shape(mlb_mask)[1]):
#         if  ((i>=24 and i<=480) and (j>((i-b_west6)/m3) and \
#             j<((i-b_east6)/m3))):
#                 mlb_mask[i,j,1] = 1.0   


# #Masking reflectivities
# new_ref = gref*mlb_mask[:,:,0]+gref*mlb_mask[:,:,1]

# #bounding points (i,j):NW(480,154),NE(480,221),SW(24,361),SE(24,428)
# mlb_P1 = np.array([[glon[480][221]    ,    glat[480][221]]  ,
#                    [glon[480][154]    ,    glat[480][154]]  ,
#                    [glon[24][361]     ,    glat[24][361]]   ,
#                    [glon[24][428]     ,    glat[24][428]] 
#                    ])

# #bounding points (i,j):NW(480,87),NE(480,154),SW(24,294),SE(24,361)  

# mlb_P2 = np.array([[glon[480][154]    ,    glat[480][154]] ,
#                    [glon[480][87]     ,    glat[480][87]]  ,
#                    [glon[24][294]     ,    glat[24][294]]  ,
#                    [glon[24][361]     ,    glat[24][261]] 
#                    ])



# #Set Tampa mask (West_1, West_2, West_3)
# #Create mask array
# tampa_mask = np.zeros((500,500,3))
#     #West_1
#     #bounding points (i,j):NW(420,140),NE(420,210),SW(100,250),SE(100,320)
# m7 = -320.0/100.0
# b_west7 = 420.0 - m7 * 140.0
# b_east7 = 420.0 - m7 * 210.0
# for i in range(np.shape(tampa_mask)[0]):
#     for j in range(np.shape(tampa_mask)[1]):
#         if  ((i>=100 and i<=420) and (j>((i-b_west7)/m7) and \
#             j<((i-b_east7)/m7))):
#                 tampa_mask[i,j,0] = 1.0

#     #West_2
#     #bounding points (i,j):NW(420,210),NE(420,280),SW(100,320),SE(100,390)     
# b_east8 = 420.0 - m7 * 280.0
# b_west8 = 420.0 - m7 * 210.0
# for i in range(np.shape(tampa_mask)[0]):
#     for j in range(np.shape(tampa_mask)[1]):
#         if  ((i>=100 and i<=420) and (j>((i-b_west8)/m7) and \
#             j<((i-b_east8)/m7))):
#                 tampa_mask[i,j,1] = 1.0   
                
#     #West_3
#     #bounding points (i,j):NW(420,280),NE(420,350),SW(100,390),SE(100,460)   
# b_east9 = 420.0 - m7 * 350.0
# b_west9 = 420.0 - m7 * 280.0
# for i in range(np.shape(tampa_mask)[0]):
#     for j in range(np.shape(tampa_mask)[1]):
#         if  ((i>=100 and i<=420) and (j>((i-b_west9)/m7) and \
#             j<((i-b_east9)/m7))):
#                 tampa_mask[i,j,2] = 1.0  


# #Masking reflectivities
# new_ref = gref*tampa_mask[:,:,0]+gref*tampa_mask[:,:,1]+gref*tampa_mask[:,:,2]

# #bounding points (i,j):NW(420,140),NE(420,210),SW(100,250),SE(100,320)
# tampa_P1 = np.array([[glon[420][210]   ,    glat[420][210]] ,
#                      [glon[420][140]   ,    glat[420][140]] ,
#                      [glon[100][250]   ,    glat[100][250]] ,
#                      [glon[100][320]   ,    glat[100][320]] 
#                      ])

# #bounding points (i,j):NW(420,210),NE(420,280),SW(100,320),SE(100,390)   
# tampa_P2 = np.array([[glon[420][210]   ,    glat[420][210]]  ,
#                      [glon[420][280]   ,    glat[420][280]]  ,
#                      [glon[100][390]   ,    glat[100][390]]  ,
#                      [glon[150][320]   ,    glat[100][320]] 
#                      ])

# #bounding points (i,j):NW(420,280),NE(420,350),SW(100,390),SE(100,460)  
# tampa_P3 = np.array([[glon[420][280]   ,    glat[420][280]]  ,
#                      [glon[420][350]   ,    glat[420][350]]  ,
#                      [glon[100][460]   ,    glat[100][460]]  ,
#                      [glon[100][390]   ,    glat[100][390]] 
#                      ])



#######################################
#######################################
#PLOTTING
#######################################
#######################################
fig = plt.figure(figsize=[8, 8], dpi=100)
my_ax = plt.axes(projection = ccrs.PlateCarree())

r = my_ax.pcolormesh(glon, glat, new_ref, cmap=pyart.graph.cm.NWSRef,
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

plt.title('KAMX NEXRAD {0}'.format(date))
my_ax.set(aspect=1,
        xlim=(glon[0][0],glon[0][-1]),
        ylim=(glat[0][1],glat[-1][0]))

P1 = matplotlib.patches.Polygon(miami_P1, \
    facecolor = 'none', edgecolor= 'green', closed=True)
P2 = matplotlib.patches.Polygon(miami_P2, \
    facecolor = 'none', edgecolor= 'blue', closed=True)
P3 = matplotlib.patches.Polygon(miami_P3, \
    facecolor = 'none', edgecolor= 'pink', closed=True)
# P4 = matplotlib.patches.Polygon(mlb_P1, \
#     facecolor = 'none', edgecolor= 'orange', closed=True)
# P5 = matplotlib.patches.Polygon(mlb_P2, \
#     facecolor = 'none', edgecolor= 'red', closed=True)
# P6 = matplotlib.patches.Polygon(tampa_P1, \
#     facecolor = 'none', edgecolor= 'purple', closed=True)
# P7 = matplotlib.patches.Polygon(tampa_P2, \
#     facecolor = 'none', edgecolor= 'yellow', closed=True)
# P8 = matplotlib.patches.Polygon(tampa_P3, \
#     facecolor = 'none', edgecolor= 'violet', closed=True)

my_ax.add_patch(P1)
my_ax.add_patch(P2)
my_ax.add_patch(P3)
# my_ax.add_patch(P4)
# my_ax.add_patch(P5)
# my_ax.add_patch(P6)
# my_ax.add_patch(P7)
# my_ax.add_patch(P8)

for st in stations:
    plt.plot(locs.loc[st][3],locs.loc[st][2],'gp',markersize=8)

# # For Tampa
# plt.xticks(np.arange(-81.20,-83.21,-0.25))
# plt.yticks(np.arange(26.5,28.6,0.25))

# For Miami
plt.xticks(np.arange(-80.0,-82.1,-0.25))
plt.yticks(np.arange(25.1,27.1,0.25))

# # Turn on for Melbourne Radar:
# plt.xticks(np.arange(-80.0,-82.1,-0.25))
# plt.yticks(np.arange(27.0,29.1,0.25))


plt.show()


