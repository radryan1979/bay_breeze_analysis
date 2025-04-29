"""
Coded by: Dan Moore

This program will ingest a list of dates from the typing
program.

We will emulate the NWS Typing Probability of Precipitation
by Regime types.

The output will be several analyses of time-series post-
SBF passage, including rainfall coverage and intensity.

Updated: 4-3-19
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


reg_infile = '/Volumes/LaCie/SeaBreeze/Florida/TampaRegimeTyping/TampaTypingUpdate.csv'
times_dir = '/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/'
NCEPdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/"



times_infile = times_dir+'Hughes.csv'
station_df = pd.read_csv(times_infile,header=0,index_col=0)
station_df["Date"] = pd.to_datetime(station_df["Date"])



#THIS SECTION CHANGES IF WE CHANGE DATE RANGE FROM SUMMER TO ALL:
full_dates = pd.date_range("01 April 2008","October 31 2018", freq="D")
dates_mask = (full_dates.month>=6) & (full_dates.month<=8)
full_dates=full_dates[dates_mask]

station_df = station_df[(station_df["Date"].dt.month>=6) & (station_df["Date"].dt.month<=8)]


df = pd.read_csv(reg_infile,header=0,index_col=0)
df = df[(df["Month"]>=4) & (df["Month"]<=10)]
regimes = df[df.Regime.notna()]
regimes = regimes.Regime.unique()

#Separate regimes into date arrays for each:
reg1 = df.loc[df.Regime==1.0]
reg1_dates = np.array([])
for _,row in reg1.iterrows():
    year=int(row["Year"])
    month=int(row["Month"])
    day=int(row["Day"])
    dtime=datetime(year,month,day)
    reg1_dates = np.append(reg1_dates,dtime)
reg2 = df.loc[df.Regime==2.0]
reg2_dates = np.array([])
for _,row in reg2.iterrows():
    year=int(row["Year"])
    month=int(row["Month"])
    day=int(row["Day"])
    dtime=datetime(year,month,day)
    reg2_dates = np.append(reg2_dates,dtime)
reg3 = df.loc[df.Regime==3.0]
reg3_dates = np.array([])
for _,row in reg3.iterrows():
    year=int(row["Year"])
    month=int(row["Month"])
    day=int(row["Day"])
    dtime=datetime(year,month,day)
    reg3_dates = np.append(reg3_dates,dtime)
reg4 = df.loc[df.Regime==4.0]
reg4_dates = np.array([])
for _,row in reg4.iterrows():
    year=int(row["Year"])
    month=int(row["Month"])
    day=int(row["Day"])
    dtime=datetime(year,month,day)
    reg4_dates = np.append(reg4_dates,dtime)
reg5 = df.loc[df.Regime==5.0]
reg5_dates = np.array([])
for _,row in reg5.iterrows():
    year=int(row["Year"])
    month=int(row["Month"])
    day=int(row["Day"])
    dtime=datetime(year,month,day)
    reg5_dates = np.append(reg5_dates,dtime)
reg6 = df.loc[df.Regime==6.0]
reg6_dates = np.array([])
for _,row in reg6.iterrows():
    year=int(row["Year"])
    month=int(row["Month"])
    day=int(row["Day"])
    dtime=datetime(year,month,day)
    reg6_dates = np.append(reg6_dates,dtime)
reg7 = df.loc[df.Regime==7.0]
reg7_dates = np.array([])
for _,row in reg7.iterrows():
    year=int(row["Year"])
    month=int(row["Month"])
    day=int(row["Day"])
    dtime=datetime(year,month,day)
    reg7_dates = np.append(reg7_dates,dtime)
reg8 = df.loc[df.Regime==8.0]
reg8_dates = np.array([])
for _,row in reg8.iterrows():
    year=int(row["Year"])
    month=int(row["Month"])
    day=int(row["Day"])
    dtime=datetime(year,month,day)
    reg8_dates = np.append(reg8_dates,dtime)
reg9 = df.loc[df.Regime==9.0]
reg9_dates = np.array([])
for _,row in reg9.iterrows():
    year=int(row["Year"])
    month=int(row["Month"])
    day=int(row["Day"])
    dtime=datetime(year,month,day)
    reg9_dates = np.append(reg9_dates,dtime)


#Setting up sea breeze detection date arrays
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

for reg in regimes:
    
    if reg == 1.0:
        reg_dates = pd.DatetimeIndex(reg1_dates)
    elif reg == 2.0:
        reg_dates = pd.DatetimeIndex(reg2_dates)
    elif reg == 3.0:
        reg_dates = pd.DatetimeIndex(reg3_dates)
    elif reg == 4.0:
        reg_dates = pd.DatetimeIndex(reg4_dates)
    elif reg == 5.0:
        reg_dates = pd.DatetimeIndex(reg5_dates)
    elif reg == 6.0:
        reg_dates = pd.DatetimeIndex(reg6_dates)
    elif reg == 7.0:
        reg_dates = pd.DatetimeIndex(reg7_dates)
    elif reg == 8.0:
        reg_dates = pd.DatetimeIndex(reg8_dates)
    elif reg == 9.0:
        reg_dates = pd.DatetimeIndex(reg9_dates)
        

    
    for cst in ["All","EC_SB","WC_SB","BC_SB","Non_SB"]:

        sbcount=0
        
        if cst == "All":
            dates = pd.DatetimeIndex(full_dates)
        elif cst == "EC_SB":
            dates = pd.DatetimeIndex(east_days["Date"].values)
        elif cst == "WC_SB":
            dates = pd.DatetimeIndex(west_days["Date"].values)
        elif cst == "BC_SB":
            dates = pd.DatetimeIndex(bothcst_days["Date"].values)
        else:
            dates = pd.DatetimeIndex(nonsb_days.values)
        
        dates = dates[dates.isin(reg_dates)]
        print(reg,cst,len(dates))
        
            
        for date in dates:
            
            SBF_hour = "13"#str(int(row['Hour']/1)) ###Just looking for SB days, not necessarily 12 hours after SB passage
            SBF_min = "00"#str(int(row['Hour']%1*60))
            date_str = str(date.day).zfill(2)+"-"+str(date.month).zfill(2)+"-"+str(date.year)+' '+SBF_hour+':'+SBF_min
            SBF_time = datetime.strptime(date_str,'%d-%m-%Y %H:%M')
        
                
                
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
            
                temp_POP = event_precip>event_threshold 
                temp_POP = temp_POP.astype(np.int)
                
                if sbcount == 0:
                    POP = temp_POP.copy()
                else:
                    POP = POP+temp_POP
            
            
            if sbcount == 0:
                accum_precip = event_precip.copy()
            else:
                accum_precip = accum_precip+event_precip
            
            sbcount+=1
                
                
        POP = POP/np.float(sbcount)*100.0
        POP[POP==0]=np.nan
        mean_accum_cst = np.divide(accum_precip,float(11))
        mean_accum_cst[mean_accum_cst==0] = np.nan
        lats = dataset.variables['lats'][:,:]
        lons = dataset.variables['lons'][:,:]
        
        
        if cst == "All":
            title = "Regime {0} Overall 12-03UTC".format(int(reg))
        elif cst == "EC_SB":
            title = "Regime {0} East Coast \nClassic SB 12-03UTC".format(int(reg))
        elif cst == "WC_SB":
            title = "Regime {0} West Coast \nClassic SB 12-03UTC".format(int(reg))
        elif cst == "BC_SB":
            title = "Regime {0} Classic SB \nBoth Coasts 12-03UTC".format(int(reg))
        else:
            title = "Regime {0} No Classic SB \n12-03UTC".format(int(reg))
        
        
        #PLOT PoP
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor('k')
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,100))
        sm._A = []
        plt.contourf(lons, lats, POP, 100,
                     transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
        plt.title(title+" PoP")
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
        plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/PrecipByRegime15Hr/Reg_{0}_SBC_PoP_{1}.png".format(int(reg),cst), 
            dpi=500)
        plt.show()
        
        
        
        """
        #Plot percentage of summertime rainfall:
        rain_perc = np.divide(mean_accum_cst,mean_accum)#mean_accum calculated in "NCEP_Rainfall_Amount.py"
        rain_perc = np.multiply(rain_perc,100.0)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor('k')
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=30)
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,30))
        sm._A = []
        states = cfeature.NaturalEarthFeature(category="cultural",
                            name="admin_1_states_provinces_lines",
                            scale="10m", facecolor="none")
        coast = cfeature.NaturalEarthFeature(category="physical",
                            name="coastline",
                            scale="10m", facecolor="none")
        ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        plt.contourf(lons, lats, rain_perc, 30,
                     transform=ccrs.PlateCarree(), norm=norm, cmap=cmap)
        plt.title("Percentage of Summertime Rainfall \n During "+title)
        plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100]).set_label("Percent of Summertime Rainfall (%)")
        ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
        plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/PrecipByRegime15Hr/PercPrecip/{0}_Reg{1}_PercentPrecip.png".format(cst,int(reg)), 
            dpi=500)
        plt.show()
        """
        
        """
        #PLOT STUDY AREA
        area = np.empty_like(POP)
        area[:,:] = 20
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor('k')
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,100))
        sm._A = []
        plt.contourf(lons, lats, area, 100,
                     transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
        plt.title("Florida NCEP Subset Study Area")
        # plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100],extend='both').set_label("PoP (%)")
        states = cfeature.NaturalEarthFeature(category="cultural",
                            name="admin_1_states_provinces_lines",
                            scale="10m", facecolor="none")
        coast = cfeature.NaturalEarthFeature(category="physical",
                            name="coastline",
                            scale="10m", facecolor="none")
        ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.set(xlim=(-90,-76.5),ylim=(23.5,34.5))
        plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/NCEP_Study_Region.png".format(int(reg),cst), 
            dpi=500)
        plt.show()
        """
        """
        ####Import binary ag NCEP grid
        agBinarypath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/NCEP_FL_AgBinary.csv" #Created using ArcMap
        NCEP_ag = pd.read_csv(agBinarypath)
        NCEP_ag = NCEP_ag.drop(["FID","XCoord","YCoord"],axis=1)
        ag_binary = np.array(NCEP_ag["AG_BINARY"]).reshape(lats.shape[0],lats.shape[1])
        ag_mask = np.multiply(mean_accum,ag_binary)
        ag_mask[ag_mask==0]=np.nan
        ag_mask[ag_mask>1]=10#PICK COLOR
        
        landBinaryPath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/NCEP_LandSea_Binary.csv" #Created using ArcMap
        NCEP_land = pd.read_csv(landBinaryPath)
        NCEP_land = NCEP_land.drop(["XCoord","YCoord"],axis=1)
        land_binary = np.array(NCEP_land["LANDBINARY"]).reshape(lats.shape[0],lats.shape[1])
        land_mask = np.multiply(mean_accum,land_binary)
        land_mask[land_mask==0]=np.nan
        land_mask[land_mask>0]=80#PICK COLOR
        
        #PLOT STUDY AREA PLUS LAND SEA PLUS AG
        area = np.empty_like(POP)
        area[:,:] = 20
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor('k')
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,100))
        sm._A = []
        plt.contourf(lons, lats, area, 100,
                     transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
        plt.contourf(lons, lats, land_mask, 100,
                     transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
        plt.contourf(lons, lats, ag_mask, 100,
                     transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
        plt.title("Florida NCEP Subset Study Area")
        # plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100],extend='both').set_label("PoP (%)")
        states = cfeature.NaturalEarthFeature(category="cultural",
                            name="admin_1_states_provinces_lines",
                            scale="10m", facecolor="none")
        coast = cfeature.NaturalEarthFeature(category="physical",
                            name="coastline",
                            scale="10m", facecolor="none")
        ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.set(xlim=(-90,-76.5),ylim=(23.5,34.5))
        plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/NCEP_Study_Region_AgAndLand.png".format(int(reg),cst), 
            dpi=500)
        plt.show()
        
        
        #Plot Grid sizes of Ag
        area = np.empty_like(POP)
        area[:,:] = 20
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.background_patch.set_facecolor('k')
        cmap = mpl.cm.gist_ncar
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,100))
        sm._A = []
        ag_pcolor = np.ma.array(ag_mask,mask=np.isnan(ag_mask))
        land_pcolor = np.ma.array(land_mask,mask=np.isnan(land_mask))
        plt.pcolormesh(lons, lats, land_pcolor,
                     transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
        plt.pcolormesh(lons, lats, ag_pcolor,
                     transform=ccrs.PlateCarree(),cmap=cmap, norm=norm, edgecolor='k')
        # plt.title("Florida NCEP Subset Study Area")
        # plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100],extend='both').set_label("PoP (%)")
        states = cfeature.NaturalEarthFeature(category="cultural",
                            name="admin_1_states_provinces_lines",
                            scale="10m", facecolor="none")
        coast = cfeature.NaturalEarthFeature(category="physical",
                            name="coastline",
                            scale="10m", facecolor="none")
        ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        ax.set(xlim=(-81.5,-80.5),ylim=(26.65,27.4))
        plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/NCEP_Study_Region_AgInset.png".format(int(reg),cst), 
            dpi=500)
        plt.show()
        """
        
