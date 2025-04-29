"""
Coded by Dan Moore
9/10/18

The purpose of this code is to analyze output from "SBDetectPyart.py",
which is in .csv format, containing coordinates of detected sea breeze
fronts. This program will take in that data, and plot it to assist in 
troubleshooting the detection algorithm.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyart
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import pandas as pd
from datetime import datetime
import cartopy.crs as ccrs


df = pd.read_csv("/Volumes/LaCie/SeaBreeze/RadarDetectOutput/FirstTest.csv", index_col=0)
df.columns = pd.to_datetime(df.columns)

#Parameters for simple plotting
# newproj = ccrs.LambertConformal(central_latitude=lat0,central_longitude=lon0)
# cmap = 'pyart_NWSRef'
# levs = np.linspace(0,80,41,endpoint=True)
# ticks = np.linspace(0,80,9,endpoint=True)
# label = 'Radar Reflectivity Factor ($\mathsf{dBZ}$)'
# #Normalize the colormap based on the levels provided above
# norm = mpl.colors.BoundaryNorm(levs,256)


#Boundaries for desired values
boundinglat = (38.15, 39.65)
boundinglon = (-76.45, -74.6)

unique_date = df.columns.map(pd.Timestamp.date).unique()

print(len(unique_date))

for date in unique_date:
    for column in df:
        lons = np.empty(len(df.index)); lats = np.empty(len(df.index))
        
        if column.date() == date:
            if np.nansum(df[column][:])>-1.0:
                continue
            for i in range(len(lons)):
                lons[i], lats[i] = df[column][df.index[i]], df.index[i]
            # Plotting rectilinear grid
            fig = plt.figure(figsize=[12,12], dpi=300)
            ax = fig.add_subplot(1,1,1, projection = ccrs.Mercator())
            ax.set_extent([boundinglon[0], boundinglon[1], boundinglat[0], boundinglat[1]])
            lon_g = list(lons)
            lat_g = list(lats)
            ax.scatter(lon_g,lat_g,transform=ccrs.Geodetic())
            ax.coastlines(resolution='50m')
            ax.gridlines()
            plt.title(str(column.strftime("%Y-%m-%d %H:%M:%S")))
            plt.show()
            
            input("Press Enter to Continue...")
            

        else:
            continue












