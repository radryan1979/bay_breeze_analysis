import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import numpy as np
import os
import math
import datetime as dt

from cartopy import config
import cartopy.crs as ccrs

in_path ='/Users/dpm/Desktop/VRES_Data/test/'

for filename in os.listdir(in_path):
    
    if filename.endswith(".nc"):
        # print (in_path+filename)

        fname = in_path +filename
        dataset = netcdf_dataset(fname)
        U_bot = dataset.variables['UBOT'][ 0, :, :]
        V_bot = dataset.variables['VBOT'][ 0, :, :]
        Temp = dataset.variables['TS'][0,:,:]
        phi_bot = (270-(np.arctan2(V_bot,U_bot)*(180/math.pi)))
        lats = dataset.variables['lat'][:]
        lons = dataset.variables['lon'][:]
        times = dataset.variables['time'][:]
        
        #Create julian date array
        start = dt.datetime(1984,1,1)      # This is the "days since" part
        julian_date = np.array([start + dt.timedelta(times[i]) for i in range(len(times))])
        
        
        
        
        
        #Plot attempt:
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        plt.plot(lons, lats, Temp, 60,
                     transform=ccrs.PlateCarree())
        
        plt.colorbar().set_label("degK")
        
        plt.plot(lons1[27,10], lats1[27,10],'ro', markersize=10)
        plt.plot(lons1[27,11], lats1[27,11],'ro', markersize=10)
        plt.plot(lons1[26,12], lats1[26,12],'ro', markersize=10)
        plt.plot(lons1[26,13], lats1[26,13],'ro', markersize=10)
        
        plt.legend()
        plt.title("NCEP Precipitation 1 May, 2018 18UTC")
        ax.coastlines()
        # plt.savefig("/Users/dpm/Desktop/Test1.png", dpi=500)
        plt.show()
        
    
