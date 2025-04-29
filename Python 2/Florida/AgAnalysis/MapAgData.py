"""
Coded by: Dan Moore

This program will map the Florida land Use data
to be used by ag analysis to determine amount 
of precipitation over agriculture areas in FL.

Updated: 3-12-19
"""

import pandas as pd
import numpy as np
from datetime import datetime,timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import os
from cartopy import config
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader
from shapely.ops import cascaded_union

fig = plt.figure(figsize=[15,15],dpi=100)
ax1 = plt.subplot(111,projection=ccrs.PlateCarree())

# states = cfeature.NaturalEarthFeature(category = 'cultural',
#                 name='admin_1_states_provinces_lines',
#                 scale='10m',facecolor='none')
# coast = cfeature.NaturalEarthFeature(category='physical',scale='10m',
#                 facecolor='none',name='coastline')

# fname='/Volumes/LaCie/GISData/FL_LandUse/FL_Agriculture_ProjectUTM17.shp'
geoms = Reader(fname).geometries()
ag_polygon = cascaded_union(list(geoms))
# ag_land=ShapelyFeature(geoms,ccrs.UTM(17), facecolor='none')

ax1.add_feature(ag_land, facecolor='blue', edgecolor='lightgray')
# ax1.add_feature(states, linestyle='-', edgecolor='darsklategray', linewidth=1)
# ax1.add_feature(coast, linestyle='-', edgecolor='darsklategray',linewidth=1)

ax1.set(aspect=1, xlim=(-83.5,-79.2),
                  ylim=(25.0,30.0))
                  
plt.show()
