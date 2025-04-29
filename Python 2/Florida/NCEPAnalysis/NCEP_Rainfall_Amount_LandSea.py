"""
Coded by: Dan Moore

This program will ingest a list of dates from the sea 
breeze detection and analyze sea-breeze vs. non-sea breeze
days and determine the distribution of daily precipitation
amounts in the ag-designated grid cells.

We will then determine a means for analyzing statistical 
meaning of the two sets of dates.

Update 5-28-19: Changed dates such that east coast and west
coast days do not contain dates in both coast sb dataframe.

Updated: 5-28-19
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
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

####FIRST: import NCEP and mask to ag data:
####Import NCEP grid
NCEPdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/"
NCEPstr = NCEPdirpath + "2017" + "/" + "2017" + "05" + \
                    "01" + "13" + ".nc"
dataset = netcdf_dataset(NCEPstr)
lats = dataset.variables['lats'][:,:]
lons = dataset.variables['lons'][:,:]

####Import binary ag NCEP grid
agBinarypath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/NCEP_LandSea_Binary.csv" #Created using ArcMap
NCEP_ag = pd.read_csv(agBinarypath)
NCEP_ag = NCEP_ag.drop(["XCoord","YCoord"],axis=1)
ag_binary = np.array(NCEP_ag["LANDBINARY"]).reshape(lats.shape[0],lats.shape[1])
for i in range(ag_binary.shape[0]):
    for j in range(ag_binary.shape[1]):
        if lats[i][j]>30.0:
            ag_binary[i][j]=0


mean_accum = np.loadtxt(NCEPdirpath+"mean_accum.csv",dtype=bytes,delimiter=',')
mean_accum = mean_accum.astype(float)
# print(type(mean_accum))
# print(mean_accum.shape)

ave_ag_precip = np.multiply(mean_accum,ag_binary)
summer_ag_precip_accum = np.sum(ave_ag_precip)

####Begin rest of program
times_dir = '/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/'
NCEPdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/"

times_infile = times_dir+'Hughes_SST_Regime.csv'
station_df = pd.read_csv(times_infile,header=0,index_col=0)
station_df["Date"] = pd.to_datetime(station_df["Date"])

#THIS SECTION CHANGES IF WE CHANGE DATE RANGE FROM SUMMER TO ALL:
full_dates = pd.date_range("01 April 2008","October 31 2018", freq="D")
dates_mask = (full_dates.month>=6) & (full_dates.month<=8)
full_dates=full_dates[dates_mask]

station_df = station_df[(station_df["Date"].dt.month>=6) & (station_df["Date"].dt.month<=8)]

#######################################################

ag_analysis=True

east_st=[340,410,420,440]
west_st=[360,350,380,450,480,490]
east_days = station_df.loc[east_st]
# east_days["Date"]=pd.to_datetime(east_days["Date"])
# east_days = east_days[(east_days["Type"]=="DPWS") | (east_days["Type"]=="classicWS")]
east_days = east_days[east_days["Type"]=="classic"]
east_days = east_days.drop_duplicates(subset="Date")

west_days = station_df.loc[west_st]
# west_days["Date"]=pd.to_datetime(west_days["Date"])
# west_days = west_days[(west_days["Type"]=="DPWS") | (west_days["Type"]=="classicWS")]
west_days = west_days[west_days["Type"]=="classic"]
west_days = west_days.drop_duplicates(subset="Date")
bothcst_days = station_df[(station_df.Date.isin(west_days["Date"])) & (station_df.Date.isin(east_days["Date"]))]
bothcst_days = bothcst_days.drop_duplicates(subset="Date")
east_days = east_days[~east_days.Date.isin(bothcst_days["Date"])]
west_days = west_days[~west_days.Date.isin(bothcst_days["Date"])]
total_sbdays = pd.concat([west_days["Date"],east_days["Date"],bothcst_days["Date"]])
total_sbdays = total_sbdays.drop_duplicates()
nonsb_days = full_dates[~full_dates.isin(total_sbdays)]
print(len(total_sbdays)/1012.0)

for cst in ["EC_SB","WC_SB","BC_SB","Non_SB"]:#

    sbcount=0
    
    cst_rf_array=np.array([])
    
    if cst == "EC_SB":
        dates = pd.DatetimeIndex(east_days["Date"].values)
    elif cst == "WC_SB":
        dates = pd.DatetimeIndex(west_days["Date"].values)
    elif cst == "BC_SB":
        dates = pd.DatetimeIndex(bothcst_days["Date"].values)
    else:
        dates = pd.DatetimeIndex(nonsb_days.values)
        
    for date in dates:
        
        SBF_hour = "13"#str(int(row['Hour']/1)) ###Just looking for SB days, not necessarily 12 hours after SB passage
        SBF_min = "00"#str(int(row['Hour']%1*60))
        date_str = str(date.day).zfill(2)+"-"+str(date.month).zfill(2)+"-"+str(date.year)+' '+SBF_hour+':'+SBF_min
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
                if ag_analysis:
                    precip_data = np.multiply(precip_data,ag_binary)

                if np.any(precip_data>10000.0) or np.any(precip_data.mask):
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
        
            temp_POP = event_precip>event_threshold
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

        cst_rf_array = np.append(cst_rf_array,event_precip[event_precip>0.5])

        sbcount+=1
    
    
    POP = POP/np.float(sbcount)*100.0
    POP[POP==0]=np.nan
    mean_accum_cst = np.divide(accum_precip,float(11))
    mean_accum_cst[mean_accum_cst==0] = np.nan
    lats = dataset.variables['lats'][:,:]
    lons = dataset.variables['lons'][:,:]
    
    #Percent of ag precip:
    if ag_analysis:
        cst_ag_precip_total = np.nansum(mean_accum_cst)
        perc_ag_precip = cst_ag_precip_total/summer_ag_precip_accum*100
        print("{0} days make up {1} of total days.\n The percent of daytime precipitation over land on {0} days is {2}.".format(cst,len(dates)/1012.0*100,perc_ag_precip))
    
    
    if cst == "EC_SB":
        title = "East Coast Classic SB"
    elif cst == "WC_SB":
        title = "West Coast Classic SB"
    elif cst == "BC_SB":
        title = "Classic SB Both Coasts"
    else:
        title = "No Classic SB"
    
    """
    if ag_analysis:
    
        #PLOT PoP
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor('k')
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,100))
        sm._A = []
        plt.contourf(lons, lats, POP, 60,
                     transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
                     
        plt.title("15Hr PoP Given \n {0}".format(title))
        plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100],extend='both').set_label("PoP (%)")
        states = cfeature.NaturalEarthFeature(category="cultural",
                            name="admin_1_states_provinces_lines",
                            scale="10m", facecolor="none")
        coast = cfeature.NaturalEarthFeature(category="physical",
                            name="coastline",
                            scale="10m", facecolor="none")
        ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
        # plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/LandSea/SBC_PoP_{0}_Classic.png".format(cst), 
        #     dpi=500)
        plt.show()
        
        
        #PLOT mean
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor('k')
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=mean_accum.max())
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,mean_accum.max()))
        sm._A = []
        plt.contourf(lons, lats, mean_accum_cst, 60,
                     transform=ccrs.PlateCarree(), norm=norm, cmap=cmap)
        plt.title("Average Summer Daytime Rainfall\n Accumulation {0}".format(title))
        plt.colorbar(sm).set_label("Rainfall Accumulation (mm)")#ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]).set_label("Rainfall Accumulation (mm)")
        states = cfeature.NaturalEarthFeature(category="cultural",
                            name="admin_1_states_provinces_lines",
                            scale="10m", facecolor="none")
        coast = cfeature.NaturalEarthFeature(category="physical",
                            name="coastline",
                            scale="10m", facecolor="none")
        ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
        # plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/LandSea/{0}_PrecipAccum.png".format(cst),
        #     dpi=500)
        plt.show()
        
        
        #Plot percentage of summertime rainfall:
        rain_perc = np.divide(mean_accum_cst,mean_accum)#mean_accum calculated in "NCEP_Rainfall_Amount.py"
        rain_perc = np.multiply(rain_perc,100.0)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor('k')
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=50)
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,50))
        sm._A = []
        ax.contourf(lons, lats, rain_perc, 50,
                     transform=ccrs.PlateCarree(), norm=norm, cmap=cmap)
        plt.title("Percentage of Summertime Rainfall Accumulation \n During Days With {0}".format(title))
        plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100]).set_label("Percent of Summertime Rainfall (%)")
        states = cfeature.NaturalEarthFeature(category="cultural",
                            name="admin_1_states_provinces_lines",
                            scale="10m", facecolor="none")
        coast = cfeature.NaturalEarthFeature(category="physical",
                            name="coastline",
                            scale="10m", facecolor="none")
        ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
        # plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/LandSea/{0}_Classic_PercentPrecip.png".format(cst),
        #     dpi=500)
        plt.show()
        
        
        #Distribution of daily rainfall totals
        bins = np.linspace(0, 100, 20)
        x=cst_rf_array
        #Calculating Interquartile Range as a measure of variance:
        q75, q25 = np.percentile(x, [75 ,25])
        iqr = q75 - q25
        #Plot distribution:
        n1,_,_=plt.hist(x,bins,alpha=1)
        plt.xlabel('Daytime Rainfall Accumulation [mm]')
        plt.ylabel('Frequency')
        plt.title('Daily Daytime Precipitation \n During Days With {0}'.format(title))
        plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/3.0, 'Med=%(Med)2.2f mm, \nIQR=%(IQR)2.2f mm' %\
                 {'Med':np.median(x),'IQR':iqr})
        plt.xlim(0)
        # plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/LandSea/DailyPrecip_{0}_Distribution.png".format(cst),
        #     dpi=500)
        plt.show()
        
    else:
        #PLOT PoP
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor('k')
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,100))
        sm._A = []
        plt.contourf(lons, lats, POP, 60,
                     transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
                     
        plt.title("15Hr PoP Given \n {0}".format(title))
        plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100],extend='both').set_label("PoP (%)")
        states = cfeature.NaturalEarthFeature(category="cultural",
                            name="admin_1_states_provinces_lines",
                            scale="10m", facecolor="none")
        coast = cfeature.NaturalEarthFeature(category="physical",
                            name="coastline",
                            scale="10m", facecolor="none")
        ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
        plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/15HrPrecip/SBC_PoP_{0}_Classic.png".format(cst), 
            dpi=500)
        plt.show()
        
        
        #PLOT mean
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor('k')
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=mean_accum.max())
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,mean_accum.max()))
        sm._A = []
        plt.contourf(lons, lats, mean_accum_cst, 60,
                     transform=ccrs.PlateCarree(), norm=norm, cmap=cmap)
        plt.title("Average Summer Daytime Rainfall\n Accumulation {0}".format(title))
        plt.colorbar(sm).set_label("Rainfall Accumulation (mm)")#ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]).set_label("Rainfall Accumulation (mm)")
        states = cfeature.NaturalEarthFeature(category="cultural",
                            name="admin_1_states_provinces_lines",
                            scale="10m", facecolor="none")
        coast = cfeature.NaturalEarthFeature(category="physical",
                            name="coastline",
                            scale="10m", facecolor="none")
        ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
        plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/15HrPrecip/{0}_PrecipAccum.png".format(cst),
            dpi=500)
        plt.show()
        
        
        #Plot percentage of summertime rainfall:
        rain_perc = np.divide(mean_accum_cst,mean_accum)#mean_accum calculated in "NCEP_Rainfall_Amount.py"
        rain_perc = np.multiply(rain_perc,100.0)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor('k')
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=50)
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,50))
        sm._A = []
        ax.contourf(lons, lats, rain_perc, 50,
                     transform=ccrs.PlateCarree(), norm=norm, cmap=cmap)
        plt.title("Percentage of Summertime Rainfall Accumulation \n During Days With {0}".format(title))
        plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100]).set_label("Percent of Summertime Rainfall (%)")
        states = cfeature.NaturalEarthFeature(category="cultural",
                            name="admin_1_states_provinces_lines",
                            scale="10m", facecolor="none")
        coast = cfeature.NaturalEarthFeature(category="physical",
                            name="coastline",
                            scale="10m", facecolor="none")
        ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
        plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/15HrPrecip/{0}_Classic_PercentPrecip.png".format(cst),
            dpi=500)
        plt.show()
        
    """
        
# print(chisquare([21.4393, 28.4774, 36.6627,14.7139], f_exp=[22.0356, 30.6324, 31.5217,15.8103]))



