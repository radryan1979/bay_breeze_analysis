"""
Coded by: Dan Moore

This program will ingest a list of dates from the typing
program.

We will create maps of mean accumulation totals for each
regime over the period 12UTC-00UTC, in addition to 
standard deviation of accumulation totals.

Updated: 2-19-19
"""


import pandas as pd
import numpy as np
from datetime import datetime,timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import os
import sys

from cartopy import config
import cartopy.crs as ccrs


reg_infile = '/Volumes/LaCie/SeaBreeze/Florida/TampaRegimeTyping/TampaTypingUpdate.csv'

NCEPdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/"



df = pd.read_csv(reg_infile,header=0,index_col=0)
regimes = df[df.Regime.notna()]
regimes = regimes.Regime.unique()

for reg in range(1,10):
    
    temp_df = df.loc[df.Regime==reg]
    
    count=0
        
    for _,row in temp_df.iterrows():

        if row["Year"]==2008.0:
            continue
        if row["Month"]<4 or row["Month"]>10:
            continue
        
        # SBF_hour = str(int(row['Hour']/1))
        # SBF_min = str(int(row['Hour']%1*60))
        # date_str = row['Date']+' '+SBF_hour+':'+SBF_min
        # SBF_time = datetime.strptime(date_str,'%d-%b-%Y %H:%M')
        
        year=int(row["Year"])
        month=int(row["Month"])
        day=int(row["Day"])
        dtime=datetime(year,month,day,13)

        
        
        filenotfound=0
        
        accum_hours = 12 #X hour accumulation analysis. 
        event_threshold = 0.05 #(in mm) accumulation total for storms - probability of precipitation over this amount.
        
        for i in range(accum_hours):
            
            #Begin the next hour after SBF detection (this includes the SBF detection time) - past 1 hour.
            begtime = dtime+timedelta(hours=i)
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
                    event_precip = event_precip + precip_data
            else:
                print(NCEPstr, "File not found.")
                filenotfound+=1
        
        if filenotfound==accum_hours:
            continue

        if count == 0:
            accum_precip = event_precip.copy()
            sum_squares = np.square(event_precip)
        else:
            accum_precip = accum_precip+event_precip
            sum_squares = sum_squares + np.square(event_precip)

        count+=1
    
    mean_accum = np.divide(accum_precip,float(count))
    
    #Standard Deviation Calculation: https://www.strchr.com/standard_deviation_in_one_pass
    # sum_squares_2 = np.divide(sum_squares,float(count))
    # mean_squared = np.square(mean_accum)
    # stdev = np.sqrt(np.abs(sum_squares - mean_squared))
    # cor_var = np.divide(stdev,mean_accum)
    
    #Standard Deviation Calculation: http://mathcentral.uregina.ca/QQ/database/QQ.09.06/h/murtaza1.html
    sum_squared = np.square(accum_precip)
    sum_squares_2 = np.multiply(count,sum_squares)
    denom = float(count*(count-1))
    variance = np.divide((sum_squares_2 - sum_squared),denom)
    stdev = np.sqrt(variance)
    cor_var = np.divide(stdev,mean_accum)
    
    lats = dataset.variables['lats'][:,:]
    lons = dataset.variables['lons'][:,:]
        
        
        
    #PLOT Mean    
    ax = plt.axes(projection=ccrs.PlateCarree())
    # cmap = mpl.cm.gist_ncar
    # norm = mpl.colors.Normalize(vmin=0, vmax=15)
    plt.contourf(lons, lats, mean_accum,# 60,
                 transform=ccrs.PlateCarree())#, norm=norm, cmap=cmap)
    plt.title("Regime {0} Average 12Hr \nRainfall Accumulation".format(int(reg)))
    plt.colorbar()#ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]).set_label("Rainfall Accumulation (mm)")
    # plt.clim(0,100)
    ax.coastlines()
    
    plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results1_17_19/Figures/RegimeAccumPrecip/Regime{0}MeanAccumPrecip.png".format(int(reg)),
        dpi=500)
    plt.show()




    #PLOT Standard Deviation
    ax = plt.axes(projection=ccrs.PlateCarree())
    # cmap = mpl.cm.gist_ncar
    # norm = mpl.colors.Normalize(vmin=0, vmax=15)
    plt.contourf(lons, lats, stdev,# 60,
                 transform=ccrs.PlateCarree())#, norm=norm, cmap=cmap)
    plt.title("Regime {0} Standard Deviation 12Hr \nRainfall Accumulation".format(int(reg)))
    plt.colorbar()#ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]).set_label("Rainfall Accumulation (mm)")
    # plt.clim(0,100)
    ax.coastlines()
    
    plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results1_17_19/Figures/RegimeAccumPrecip/Regime{0}STDPrecip.png".format(int(reg)),
        dpi=500)
    plt.show()
    
    
    
    #PLOT Correlation of Variance
    ax = plt.axes(projection=ccrs.PlateCarree())
    # cmap = mpl.cm.gist_ncar
    # norm = mpl.colors.Normalize(vmin=0, vmax=15)
    plt.contourf(lons, lats, cor_var,# 60,
                 transform=ccrs.PlateCarree())#, norm=norm, cmap=cmap)
    plt.title("Regime {0} Correlation of Variance 12Hr \nRainfall Accumulation".format(int(reg)))
    plt.colorbar()#ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]).set_label("Rainfall Accumulation (mm)")
    # plt.clim(0,100)
    ax.coastlines()
    
    plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results1_17_19/Figures/RegimeAccumPrecip/Regime{0}CVPrecip.png".format(int(reg)),
        dpi=500)
    plt.show()

