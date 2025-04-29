import pyart
import numpy as np
import os
from pathlib import Path

for yr in range(2013,2017):

    pathlist = Path('/home/work/geog.rauscher/DATA/NEXRAD/%d/NetCDF/' %(yr)).glob('*.nc')
    outdirectory='/home/work/geog.rauscher/DATA/NEXRAD/%d/GriddedNetCDF/' %(yr)
    
    for path in pathlist:
        
        filename=str(path)
        
        outfile=filename[-26:]
    
        radardat=pyart.io.read(filename)
        
        z=1
        y=635
        x=800
        
        zlim=(1000.0,1000.0)
        ylim=(-250000.0,250000.0)
        xlim=(-250000.0,250000.0)
        
        griddedrad=pyart.map.grid_from_radars(radardat,grid_shape=(z,y,x),grid_limits=(zlim,ylim,xlim),write_point_lon_lat_alt=True)
        
        #lon,lat=griddedrad.get_point_longitude_latitude()
        
        lat=griddedrad.point_latitude
        lon=griddedrad.point_longitude
        
        griddedrad.add_field('lat',lat)
        griddedrad.add_field('lon',lon)
        
        pyart.io.write_grid(outdirectory+outfile,griddedrad,format='NETCDF4')
    
    
