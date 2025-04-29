import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import numpy as np
import os

from cartopy import config
import cartopy.crs as ccrs

fname1='/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/2018/2018050120.nc'
dataset1 = netcdf_dataset(fname1)
precip1 = dataset1.variables['P'][ :, :]
lats1 = dataset1.variables['lats'][:,:]
lons1 = dataset1.variables['lons'][:,:]

fname2='/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/2018/2018050118.nc'
dataset2 = netcdf_dataset(fname2)
precip2 = dataset2.variables['P'][ :, :]
lats2 = dataset2.variables['lats'][:,:]
lons2 = dataset2.variables['lons'][:,:]

temp_precip1 = precip1>0.5
temp_precip1 = temp_precip1.astype(np.int)
temp_precip2 = precip2>0.5
temp_precip2 = temp_precip2.astype(np.int)

precip_cum = (temp_precip1+temp_precip2).astype(np.float16)
precip_cum = precip_cum/2.0*100.0

ax = plt.axes(projection=ccrs.PlateCarree())

plt.contourf(lons2, lats2, precip_cum, 60,
             transform=ccrs.PlateCarree())

plt.colorbar().set_label("mm/hr")

plt.plot(lons1[60,250], lats1[60,250],'ro', markersize=11, label="West Coast")
plt.plot(lons1[75,190], lats1[75,190],'bo', markersize=11, label="SE Coast")

plt.legend()

plt.title("NCEP Precipitation 1 May, 2018 18UTC")

ax.coastlines()

plt.savefig("/Users/dpm/Desktop/Test1.png", dpi=500)
plt.show()