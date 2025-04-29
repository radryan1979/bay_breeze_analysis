import iris
import iris_grib
import xarray as xr
import matplotlib.pyplot as plt

 

cubes = iris.load_cube('/Volumes/LaCie/SeaBreeze/RawData/GCIP_EOP_GriddedPrecip/128.117.83.221/pub/download/data/dpmud121173/ST4.2009040100.06h')       # each variable in the netcdf file is a cube
iris.save(cubes[0],'/Users/dpm/Desktop/output1.nc')  # save a specific variable to grib

ds = xr.open_dataset('/Users/dpm/Desktop/output1.nc')

