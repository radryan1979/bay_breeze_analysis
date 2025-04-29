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



###NEED TO INSTALL THESE ON FARBER
import xarray as xr
import xesmf as xe

station='KMLB'

#Boundaries for desired values
boundinglat = (27.375, 28.875)
boundinglon = (-81.725, -79.875)

boundinglon=(-105.19,76.159)
boundinglat = (-81.79,84.44)

#Defining the FL Coastline
cstlnlat = np.arange(28.65,27.64,-0.025)

cstlnlon = np.arange(-80.57,-80.14,0.0105)
    


    

w_dir = "/Users/dpmoore2927/Desktop/TestNC/"

radar = pyart.io.read_cfradial(w_dir+'2013_8_16FLCaseStudy.nc')

#Establish location of radar
loc = pyart.io.nexrad_common.get_nexrad_location(station)
lon0 = loc[1] ; lat0 = loc[0]



###DEPENDS ON WHICH RADAR - COULD TAKE PLACE OF BOUNDING LAT/LON???
ds_out = xe.util.grid_2d(-81.725, -79.875, 0.003, 27.375, 28.875, 0.003)
###I THINK 0.005 IS GOOD ENOUGH - changes number of 'pixels'
###to analyze
###BUT DON"T NEED TO TRIM DATA NOW, THIS DOES IT FOR US

###MUST BE INSIDE DETECT ALGORITHM
display = pyart.graph.RadarMapDisplayCartopy(radar)

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

gref=np.random.rand(np.shape(lons)[0],np.shape(lats)[1])
ref = radar.get_field(0, 'reflectivity')

ds=xr.Dataset({'reflectivity': (['i', 'j'],  ref)},
        coords={'lon': (['i', 'j'], lons),
        'lat': (['i', 'j'], lats)})

grid_in = {'lon': lons, 'lat': lats,
    'lon_b': lons_b, 'lat_b': lats_b}

regridder = xe.Regridder(grid_in, ds_out, 'conservative')
data_out=regridder(ref)
#return data_out, ds_out.lon, ds_out.lat
####################################



fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(111)
r = ax.pcolormesh(ds_out.lon_b, ds_out.lat_b, data_out, cmap=pyart.graph.cm.NWSRef,
                        norm=mpl.colors.Normalize(vmin=0,vmax=60))




fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(111)
t = ax.pcolormesh(lons, lats, ref)
display.set_limits(xlim=(-80.8, -80.75), ylim=(27.6, 27.65), ax=ax)
r = ax.pcolormesh(ds_out.lon_b, ds_out.lat_b, data_out,alpha=0.5,edgecolor='k')
plt.show()



