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


times_dir = '/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/SST_Indices/'

NCEPdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/"

# overall_pop = POP.copy()

agBinarypath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/NCEP_LandSea_Binary.csv" #Created using ArcMap
NCEP_ag = pd.read_csv(agBinarypath)
NCEP_ag = NCEP_ag.drop(["XCoord","YCoord"],axis=1)
ag_binary = np.array(NCEP_ag["LANDBINARY"]).reshape(lats.shape[0],lats.shape[1])

mean_accum = np.loadtxt(NCEPdirpath+"mean_accum.csv",dtype=bytes,delimiter=',')
mean_accum = mean_accum.astype(float)

ave_ag_precip = np.multiply(mean_accum,ag_binary)
summer_ag_precip_accum = np.sum(ave_ag_precip)


NE_infile = times_dir+'41009_SST_Index.csv'
SE_infile = times_dir+'MLRF1_SST_Index.csv'
W_infile = times_dir+'42013_SST_Index.csv'
NE_df = pd.read_csv(NE_infile,header=0,index_col=0)
SE_df = pd.read_csv(SE_infile,header=0,index_col=0)
W_df = pd.read_csv(W_infile,header=0,index_col=0)
NE_df.index = pd.to_datetime(NE_df.index)
SE_df.index = pd.to_datetime(SE_df.index)
W_df.index = pd.to_datetime(W_df.index)
NE_df = NE_df[(NE_df.index.month>=6) & (NE_df.index.month<=8)]
SE_df = SE_df[(SE_df.index.month>=6) & (SE_df.index.month<=8)]
W_df = W_df[(W_df.index.month>=6) & (W_df.index.month<=8)]

for cst in ["NE","SE","W"]:
    if cst=="NE":
        temp_df = NE_df.copy()
    elif cst == "SE":
        temp_df = SE_df.copy()
    else:
        temp_df = W_df.copy()
    
    for anom in ["Positive", "Negative"]:
        if anom == "Positive":
            temp_df2 = temp_df[temp_df["Index"]=="Positive"]
        else:
            temp_df2 = temp_df[temp_df["Index"]=="Negative"]
    
        dates = temp_df2.index
        SEneg_dates = dates.copy()
        
        print(cst, anom, len(dates))
        
        sbcount=0
        
        cst_rf_array=np.array([])

            
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
                    # precip_data = np.multiply(precip_data,ag_binary)
    
                    if np.any(precip_data>10000.0) or np.any(precip_data.mask):
                        print(NCEPstr, "Error in data, skipping file.")
                        filenotfound+=1
                        continue
                    
                    precip_data = np.multiply(precip_data,ag_binary)
                
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
        
        
        if cst == "NE":
            title = "North East Coast {0} SST Anomaly".format(anom)
        elif cst == "SE":
            title = "South East Coast {0} SST Anomaly".format(anom)
        elif cst == "W":
            title = "West Coast {0} SST Anomaly".format(anom)

        cst_ag_precip_total = np.nansum(mean_accum_cst)
        perc_ag_precip = cst_ag_precip_total/summer_ag_precip_accum*100
        print("{0} {1} anomaly days make up {2} of total days.\n The percent of daytime precipitation over land on {0} {1} anomaly days is {3}.".format(cst,anom,len(dates)/1012.0*100,perc_ag_precip))
    
        
        #PLOT PoP
        # ax = plt.axes(projection=ccrs.PlateCarree())
        # ax.background_patch.set_facecolor('k')
        # cmap = mpl.cm.gist_ncar
        # norm = mpl.colors.Normalize(vmin=0, vmax=100)
        # sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,100))
        # sm._A = []
        # plt.contourf(lons, lats, POP, 60,
        #              transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
        # states = cfeature.NaturalEarthFeature(category="cultural",
        #                     name="admin_1_states_provinces_lines",
        #                     scale="10m", facecolor="none")
        # coast = cfeature.NaturalEarthFeature(category="physical",
        #                     name="coastline",
        #                     scale="10m", facecolor="none")
        # ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        # ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        # plt.title("15Hr PoP Given "+title)
        # plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100],extend='both').set_label("PoP (%)")
        # ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
        # plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/SSTPrecipSummer/{0}Coast_{1}SST_PoP.png".format(cst,anom), dpi=500)
        # plt.show()
        
        # #PLOT relative PoP
        # ax = plt.axes(projection=ccrs.PlateCarree())
        # ax.background_patch.set_facecolor('k')
        # cmap = mpl.cm.seismic
        # norm = mpl.colors.Normalize(vmin=-15, vmax=15)
        # sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(-15,15))
        # sm._A = []
        # plt.contourf(lons, lats, POP-overall_pop, 60,
        #              transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
        # states = cfeature.NaturalEarthFeature(category="cultural",
        #                     name="admin_1_states_provinces_lines",
        #                     scale="10m", facecolor="none")
        # coast = cfeature.NaturalEarthFeature(category="physical",
        #                     name="coastline",
        #                     scale="10m", facecolor="none")
        # ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        # ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        # plt.title("15Hr PoP Anomaly Given \n"+title)
        # plt.colorbar(sm,ticks=[-15,-10,-5,0,5,10,15],extend='both').set_label("PoP (%)")
        # ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
        # plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/SSTPrecipSummer/{0}Coast_{1}SST_RelPoP.png".format(cst,anom), dpi=500)
        # plt.show()
        
        
        # #PLOT mean
        # ax = plt.axes(projection=ccrs.PlateCarree())
        # ax.background_patch.set_facecolor('k')
        # cmap = mpl.cm.gist_ncar
        # norm = mpl.colors.Normalize(vmin=0, vmax=mean_accum.max())
        # sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,mean_accum.max()))
        # sm._A = []
        # states = cfeature.NaturalEarthFeature(category="cultural",
        #                     name="admin_1_states_provinces_lines",
        #                     scale="10m", facecolor="none")
        # coast = cfeature.NaturalEarthFeature(category="physical",
        #                     name="coastline",
        #                     scale="10m", facecolor="none")
        # ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        # ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        # plt.contourf(lons, lats, mean_accum_cst, 60,
        #              transform=ccrs.PlateCarree(), norm=norm, cmap=cmap)
        # plt.title("Florida Average Summer Daytime \nRainfall Accumulation\n"+title)
        # plt.colorbar(sm).set_label("Rainfall Accumulation (mm)")#ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]).set_label("Rainfall Accumulation (mm)")
        # ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
        # plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/SSTPrecipSummer/{0}Coast_{1}SST_PrecipAccum.png".format(cst,anom), 
        #     dpi=500)
        # plt.show()
        
        
        # #Plot percentage of summertime rainfall:
        # rain_perc = np.divide(mean_accum_cst,mean_accum)#mean_accum calculated in "NCEP_Rainfall_Amount.py"
        # rain_perc = np.multiply(rain_perc,100.0)
        # ax = plt.axes(projection=ccrs.PlateCarree())
        # ax.background_patch.set_facecolor('k')
        # cmap = mpl.cm.gist_ncar
        # norm = mpl.colors.Normalize(vmin=0, vmax=50)
        # sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,50))
        # sm._A = []
        # states = cfeature.NaturalEarthFeature(category="cultural",
        #                     name="admin_1_states_provinces_lines",
        #                     scale="10m", facecolor="none")
        # coast = cfeature.NaturalEarthFeature(category="physical",
        #                     name="coastline",
        #                     scale="10m", facecolor="none")
        # ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
        # ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
        # plt.contourf(lons, lats, rain_perc, 50,
        #              transform=ccrs.PlateCarree(), norm=norm, cmap=cmap)
        # plt.title("Percentage of Summertime Rainfall \n During "+title)
        # plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100]).set_label("Percent of Summertime Rainfall (%)")
        # ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
        # plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/SSTPrecipSummer/{0}Coast_{1}SST_PercentPrecip.png".format(cst,anom), 
        #     dpi=500)
        # plt.show()




