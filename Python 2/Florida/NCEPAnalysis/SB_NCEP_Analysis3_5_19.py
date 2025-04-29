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

Updated: 3-27-19
"""


import pandas as pd
import numpy as np
from datetime import datetime,timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import os
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader


times_dir = '/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/'

NCEPdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/"


times_infile = times_dir+'Hughes.csv'
station_df = pd.read_csv(times_infile,header=0,index_col=0)
station_df["Date"] = pd.to_datetime(station_df["Date"])
station_df = station_df[(station_df["Date"].dt.month>=6) & (station_df["Date"].dt.month<=8)]

for cst in ["East","West"]:
    if cst=="East":
        st=[340,410,420,440]
    else:
        st=[360,350,380,450,480,490]
    
    temp_df = station_df.loc[st]
    temp_df=temp_df.drop_duplicates(subset="Date")
    temp_df=temp_df[temp_df["Type"]=="classic"]
    
    sbcount=0
        
    for _,row in temp_df.iterrows():
        
        SBF_hour = "13"#str(int(row['Hour']/1)) ###Just looking for SB days, not necessarily 12 hours after SB passage
        SBF_min = "00"#str(int(row['Hour']%1*60))
        date_str = str(row["Date"].day).zfill(2)+"-"+str(row["Date"].month).zfill(2)+"-"+str(row["Date"].year)+' '+SBF_hour+':'+SBF_min
        SBF_time = datetime.strptime(date_str,'%d-%m-%Y %H:%M')
        #SBF_time = SBF_time+timedelta(hours=4)#Local to UTC for NCEP Analysis

        
        
        filenotfound=0
        
        accum_hours = 15 #X hour accumulation analysis. 
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

                if np.any(precip_data>10000.0):
                    print(NCEPstr, "Error in data, skipping file.")
                    filenotfound+=1
                    continue
            
                if i==0:
                    event_precip = precip_data.copy()
                else:
                    event_precip = event_precip+precip_data
            else:
                print(NCEPstr, "File not found.")
                filenotfound+=1
        
        if filenotfound==accum_hours:
            continue
        else:
        
            temp_POP = event_precip>event_threshold #Events more than 1 inch
            temp_POP = temp_POP.astype(np.int)
            
            if sbcount == 0:
                POP = temp_POP.copy()
            else:
                POP = POP+temp_POP
        
        if sbcount == 0:
            accum_precip = event_precip.copy()
            # sum_squares = np.square(event_precip)
        else:
            accum_precip = accum_precip+event_precip

        sbcount+=1
    
    POP = POP/np.float(sbcount)*100.0
    mean_accum_cst = np.divide(accum_precip,float(11))
    lats = dataset.variables['lats'][:,:]
    lons = dataset.variables['lons'][:,:]
    
    
    #PLOT PoP
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.background_patch.set_facecolor('k')
    cmap = mpl.cm.gist_ncar
    norm = mpl.colors.Normalize(vmin=0, vmax=100)
    sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,100))
    sm._A = []
    plt.contourf(lons, lats, POP, 60,
                 transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
    states = cfeature.NaturalEarthFeature(category="cultural",
                        name="admin_1_states_provinces_lines",
                        scale="10m", facecolor="none")
    coast = cfeature.NaturalEarthFeature(category="physical",
                        name="coastline",
                        scale="10m", facecolor="none")
    ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
    ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
    plt.title("15Hr PoP Given \n {0} CoastClassic SB".format(cst))
    plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100],extend='both').set_label("PoP (%)")
    ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
    plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/15HrPrecip/SBC_PoP_{0}Coast.png".format(cst), dpi=500)
    plt.show()
    
    
    #PLOT mean
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.background_patch.set_facecolor('k')
    cmap = mpl.cm.gist_ncar
    norm = mpl.colors.Normalize(vmin=0, vmax=mean_accum.max())
    sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,mean_accum.max()))
    sm._A = []
    states = cfeature.NaturalEarthFeature(category="cultural",
                        name="admin_1_states_provinces_lines",
                        scale="10m", facecolor="none")
    coast = cfeature.NaturalEarthFeature(category="physical",
                        name="coastline",
                        scale="10m", facecolor="none")
    ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
    ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
    plt.contourf(lons, lats, mean_accum_cst, 60,
                 transform=ccrs.PlateCarree(), norm=norm, cmap=cmap)
    plt.title("Florida Average Summer Daytime \nRainfall Accumulation {0} Coast Classic SB".format(cst))
    plt.colorbar(sm).set_label("Rainfall Accumulation (mm)")#ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]).set_label("Rainfall Accumulation (mm)")
    ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
    plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/15HrPrecip/{0}Coast_ClassicSB_PrecipAccum.png".format(cst),
        dpi=500)
    plt.show()
    
    
    #Plot percentage of summertime rainfall:
    rain_perc = np.divide(mean_accum_cst,mean_accum)#mean_accum calculated in "NCEP_Rainfall_Amount.py"
    rain_perc = np.multiply(rain_perc,100.0)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.background_patch.set_facecolor('k')
    cmap = mpl.cm.gist_ncar
    norm = mpl.colors.Normalize(vmin=0, vmax=70)
    sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,70))
    sm._A = []
    states = cfeature.NaturalEarthFeature(category="cultural",
                        name="admin_1_states_provinces_lines",
                        scale="10m", facecolor="none")
    coast = cfeature.NaturalEarthFeature(category="physical",
                        name="coastline",
                        scale="10m", facecolor="none")
    ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
    ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
    plt.contourf(lons, lats, rain_perc, 70,
                 transform=ccrs.PlateCarree(), norm=norm, cmap=cmap)
    plt.title("Percentage of Summertime Rainfall Accumulation \n During {0} Coast Classic SB".format(cst))
    plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100]).set_label("Percent of Summertime Rainfall (%)")
    ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
    plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/15HrPrecip/{0}Coast_ClassicSB_PercentPrecip.png".format(cst),
        dpi=500)
    plt.show()




