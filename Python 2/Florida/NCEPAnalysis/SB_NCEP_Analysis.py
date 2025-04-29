"""
Coded by: Dan Moore

This program will ingest a list of dates from the Hughes
detection algorithm in Florida, and a list of Gridded Precip
Files from NCEP.

We will run analyses on the output files to statistically
analyze the timing of SBF passage, in conjunction with
convection intensity, coverage, and location.

The output will be several analyses of time-series post-
SBF passage, including rainfall coverage and intensity.

Updated: 2-13-19
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


times_dir = '/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results1_17_19/'
coasts = ['West']

NCEPdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/"


for coast in coasts:
    times_infile = times_dir+'Hughes{0}.csv'.format(coast)
    station_df = pd.read_csv(times_infile,header=0,index_col=0)
    stations = pd.unique(station_df.index)
    
    for st in stations:
        
        temp_df = station_df.loc[st]
        
        sbcount=0
            
        for _,row in temp_df.iterrows():

            if row["Date"][7:11]=='2008':
                continue
            
            SBF_hour = str(int(row['Hour']/1))
            SBF_min = str(int(row['Hour']%1*60))
            date_str = row['Date']+' '+SBF_hour+':'+SBF_min
            SBF_time = datetime.strptime(date_str,'%d-%b-%Y %H:%M')

            
            
            filenotfound=0
            
            accum_hours = 12 #X hour accumulation analysis. 
            event_threshold = 0.5 #(in mm) accumulation total for storms - probability of precipitation over this amount.
            
            for i in range(accum_hours):
                
                #Begin the next hour after SBF detection (this includes the SBF detection time) - past 1 hour.
                begtime = SBF_time+timedelta(hours=i)
                NCEPstr = NCEPdirpath + str(begtime.year) + "/" + str(begtime.year) + str(begtime.month).zfill(2) + \
                        str(begtime.day).zfill(2) + str(begtime.hour).zfill(2) + ".nc"
                
                if os.path.isfile(NCEPstr):
                
                    dataset = netcdf_dataset(NCEPstr)
                    precip_data = dataset.variables['P'][ :, :]
                    # iter_precip = precip_data>0.5 #0.5mm ~ 0.02in
                    # iter_precip = iter_precip.astype(np.int)
                    
                    if i==0:
                        event_precip = precip_data.copy()
                    else:
                        event_precip = event_precip+precip_data
                else:
                    print(NCEPstr, "File not found.")
                    filenotfound+=1
            
            if filenotfound==6:
                continue
            else:
            
                temp_POP = event_precip>event_threshold #Events more than 1 inch
                temp_POP = temp_POP.astype(np.int)
                
                if sbcount == 0:
                    POP = temp_POP.copy()
                else:
                    POP = POP+temp_POP

            sbcount+=1
        
        POP = POP/np.float(sbcount)*100.0
        lats = dataset.variables['lats'][:,:]
        lons = dataset.variables['lons'][:,:]
            
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        
        plt.contourf(lons, lats, POP, 60,
                     transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
                     
        plt.title("{0} 12Hr PoP Given SBC".format(st))
        
        plt.colorbar(ticks=[0,10,20,30,40,50,60,70,80,90,100],extend='both').set_label("PoP (%)")
        plt.clim(0,100)
        
        ax.coastlines()
        
        plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results1_17_19/SBC_PoP_{0}.png".format(st), dpi=500)
        plt.show()
        
    
    
    
    
