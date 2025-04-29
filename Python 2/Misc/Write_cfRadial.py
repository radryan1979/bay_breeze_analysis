"""
Coded by Dan Moore
Final program to ingest 5 years worth of radar data and output
coordinates of SBF for analysis.

Updated: 11_14_18
"""
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


rs=RadarServer('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/nexrad/level2/S3/')
cat = TDSCatalog('http://thredds-aws.unidata.ucar.edu/thredds/radarServer/catalog.xml')
query=rs.query()
date=datetime(2017,7,4,19,54)
query.stations('KDOX').time(date)
cat = rs.get_catalog(query)
ds = list(cat.datasets.values())
rs.validate_query(query)
catalog = rs.get_catalog(query)
dataset = list(catalog.datasets.values())[0]

# now that we have the data, let's go ahead and open it with pyart
radar = pyart.io.read_nexrad_cdm(dataset.access_urls['OPENDAP'])

pyart.io.write_cfradial('/Users/dpm/Desktop/2017_7_4_CaseStudy.nc', radar)





