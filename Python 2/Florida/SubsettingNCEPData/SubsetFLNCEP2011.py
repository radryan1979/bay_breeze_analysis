import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import numpy as np
import os

from cartopy import config
import cartopy.crs as ccrs

in_path ='/Volumes/LaCie/SeaBreeze/RawData/GCIP_EOP_GriddedPrecip/2011/'

for filename in os.listdir(in_path):
    
    if filename.endswith(".nc"):
        # print (in_path+filename)

        fname = in_path +filename
        dataset = netcdf_dataset(fname)
        precip = dataset.variables['A_PCP_GDS5_SFC_acc1h'][ :, :]
        lats = dataset.variables['g5_lat_0'][:,:]
        lons = dataset.variables['g5_lon_1'][:,:]
        
        latSubset = lats[120:330 , 800:1080] 
        lonSubset = lons[120:330 , 800:1080] 
        precipSubset = precip[120:330 , 800:1080] 
        
        out_path = '/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/2011/'
        file_name = filename[4:14]+".nc"
        file = out_path+file_name
        ncfile = netcdf_dataset(file,'w',format='NETCDF3_CLASSIC')
        
        #create dimensions
        ncfile.createDimension('lat',latSubset.shape[0])
        ncfile.createDimension('lon',lonSubset.shape[1])
        
        #define variables
        latitude = ncfile.createVariable('lats','d',('lat','lon'))
        longitude = ncfile.createVariable('lons','d',('lat','lon'))
        precipitation = ncfile.createVariable('P','d',('lat','lon'))
        longitude[:,:] = lonSubset
        latitude[:,:] = latSubset
        precipitation[:,:] = precipSubset
        
        #close ncfile
        ncfile.close()
