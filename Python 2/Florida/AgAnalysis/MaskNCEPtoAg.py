"""
Coded by: Dan Moore

This program mask the NCEP arrays to only grids
inside agricultural land use polygons to determine amount 
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
from shapely.geometry import Point
                    

####Import NCEP grid
NCEPdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/"
NCEPstr = NCEPdirpath + "2017" + "/" + "2017" + "05" + \
                    "01" + "13" + ".nc"
dataset = netcdf_dataset(NCEPstr)
lats = dataset.variables['lats'][:,:]
lons = dataset.variables['lons'][:,:]



####Import binary ag NCEP grid
agBinarypath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/NCEP_FL_AgBinary.csv" #Created using ArcMap
NCEP_ag = pd.read_csv(agBinarypath)
NCEP_ag = NCEP_ag.drop(["FID","XCoord","YCoord"],axis=1)
ag_binary = np.array(NCEP_ag["AG_BINARY"]).reshape(lats.shape[0],lats.shape[1])


