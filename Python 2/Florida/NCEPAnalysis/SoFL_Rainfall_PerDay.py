"""
Coded by: Dan Moore

This program will ingest a list of dates from the typing
program.

We will create maps of mean accumulation totals for the
entire summer over the period 12UTC-00UTC, in addition to 
standard deviation of accumulation totals.

Update 3/28/19: 15 hours instead of 12. goes until 3UTC

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
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from scipy.stats import wilcoxon,ttest_1samp,mannwhitneyu


reg_infile = '/Volumes/LaCie/SeaBreeze/Florida/TampaRegimeTyping/TampaTypingUpdate.csv'

NCEPdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/"

#Setting land-sea binary
NCEPstr = NCEPdirpath + "2017" + "/" + "2017" + "05" + \
                    "01" + "13" + ".nc"
dataset = netcdf_dataset(NCEPstr)
lats = dataset.variables['lats'][:,:]
lons = dataset.variables['lons'][:,:]

agBinarypath = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/NCEP_LandSea_Binary.csv" #Created using ArcMap
NCEP_ag = pd.read_csv(agBinarypath)
NCEP_ag = NCEP_ag.drop(["XCoord","YCoord"],axis=1)
ag_binary = np.array(NCEP_ag["LANDBINARY"]).reshape(lats.shape[0],lats.shape[1])
for i in range(ag_binary.shape[0]):
    for j in range(ag_binary.shape[1]):
        if lats[i][j]>30.0:
            ag_binary[i][j]=0


df = pd.read_csv(reg_infile,header=0,index_col=0)

    
temp_df = df.copy()

count=0
sbcount=0
mocnt=0
yrcnt=0

yr=2008
mo=6

dates = []
rf_arr = []

for _,row in temp_df.iterrows():
    
    

    if row["Month"]<6 or row["Month"]>8:
        continue
    
    # SBF_hour = str(int(row['Hour']/1))
    # SBF_min = str(int(row['Hour']%1*60))
    # date_str = row['Date']+' '+SBF_hour+':'+SBF_min
    # SBF_time = datetime.strptime(date_str,'%d-%b-%Y %H:%M')
    
    year=int(row["Year"])
    month=int(row["Month"])
    day=int(row["Day"])
    dtime=datetime(year,month,day,13)
    
    if month!=mo:
        # np.savetxt(NCEPdirpath+"total_yearly_accum_{0}_{1}.csv".format(mo,yr),monthly_precip,delimiter=',')
        mo=month
        del monthly_precip
        mocnt=0
        
        if (year != yr):
            
            if year==2009:
                print(yr)
                sum_squares = np.square(yearly_precip)
                yr = year
                print("Yearly",yearly_precip[73,220],yearly_precip[74,220])
                
                yrcnt=0
            else:
                print(yr)
                sum_squares = sum_squares + np.square(yearly_precip)
                yr = year
                print("Yearly",yearly_precip[73,220],yearly_precip[74,220])
                 
                yrcnt=0
    
    filenotfound=0
    
    accum_hours = 15 #X hour accumulation analysis. 
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
            
            if np.any(precip_data>10000.0) or np.any(precip_data.mask):
                print(NCEPstr, "Error in data, skipping file.")
                filenotfound+=1
                continue

            if i==0:
                event_precip = precip_data.copy()
            else:
                event_precip = event_precip + precip_data
        else:
            print(NCEPstr, "File not found.")
            filenotfound+=1
    
    if filenotfound==accum_hours:
        print("count # {0} no good.".format(count))
        count+=1
        continue

    temp_POP = event_precip>event_threshold
    temp_POP = temp_POP.astype(np.int)
    
    if sbcount == 0:
        POP = temp_POP.copy()
    else:
        POP = POP+temp_POP.copy()
    
    if count == 0:
        accum_precip = event_precip.copy()
        yearly_precip = event_precip.copy()
        # sum_squares = np.square(event_precip)
    else:
        accum_precip = accum_precip+event_precip
        yearly_precip = yearly_precip + event_precip
        # sum_squares = sum_squares + np.square(event_precip)
        
    #Troubleshooting:
    # print(event_precip[73,220],event_precip[74,220])
    
    if mocnt == 0:
        monthly_precip = event_precip.copy() 
    else:
        monthly_precip = monthly_precip + event_precip
    
    dates = np.append(dates,dtime.date())
    rf_arr = np.append(rf_arr,np.mean(np.multiply(event_precip.copy(),ag_binary)))
    
    mocnt+=1
    yrcnt+=1
    count+=1
    sbcount+=1

print("Done")
#2018
print(yr)
# np.savetxt(NCEPdirpath+"total_monthly_accum_{0}_{1}.csv".format(mo,yr),monthly_precip,delimiter=',')

daily_rf = pd.DataFrame({"Date": dates,
                        "Perc_Masdar_perday": rf_arr}#/np.nansum(np.multiply(accum_precip.copy(),ag_binary))},
                        )
daily_rf.Date=pd.to_datetime(daily_rf.Date)

total_sum=np.nanmean(np.multiply(accum_precip.copy(),ag_binary))

#Histogram:
bins = np.linspace(0, 6, 50)
x1=daily_rf[np.isin(daily_rf.Date,east_days.Date)]["Perc_Masdar_perday"]#["Perc_Masdar_perday"]#
x2=daily_rf[np.isin(daily_rf.Date,west_days.Date)]["Perc_Masdar_perday"]#["Perc_Masdar_perday"]#
x3=daily_rf[np.isin(daily_rf.Date,bothcst_days.Date)]["Perc_Masdar_perday"]#["Perc_Masdar_perday"]#
x4=daily_rf[np.isin(daily_rf.Date,nonsb_days)]["Perc_Masdar_perday"]#["Perc_Masdar_perday"]#

x_overall = daily_rf["Perc_Masdar_perday"]

print(mannwhitneyu(x2,x4,alternative='greater'))

data=[x1,x2,x3,x4]
fig4, ax4 = plt.subplots()
ax4.set_title('Precipitation Accumulation \non Sea Breeze Days',fontsize=18)
ax4.boxplot(data, showfliers=False,showmeans=True)
labels=["ECSB","WCSB","BCSB","NoSB"]
ax4.set_xticklabels(labels,fontsize=14)
ax4.set_ylabel('Areal Mean Precipitation [mm]',fontsize=12)
plt.tight_layout()
plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/15HrPrecip/BoxWhiskerSB.png",dpi=300)
plt.show()

x5=daily_rf[np.isin(daily_rf.Date,NEpos_dates)]["Perc_Masdar_perday"]#["Perc_Masdar_perday"]#
x5_rep=daily_rf[~np.isin(daily_rf.Date,NEpos_dates)]["Perc_Masdar_perday"]
x6=daily_rf[np.isin(daily_rf.Date,NEneg_dates)]["Perc_Masdar_perday"]#["Perc_Masdar_perday"]#
x6_rep=daily_rf[~np.isin(daily_rf.Date,NEneg_dates)]["Perc_Masdar_perday"]
x7=daily_rf[np.isin(daily_rf.Date,SEpos_dates)]["Perc_Masdar_perday"]#["Perc_Masdar_perday"]#
x7_rep=daily_rf[~np.isin(daily_rf.Date,SEpos_dates)]["Perc_Masdar_perday"]
x8=daily_rf[np.isin(daily_rf.Date,SEneg_dates)]["Perc_Masdar_perday"]#["Perc_Masdar_perday"]#
x8_rep=daily_rf[~np.isin(daily_rf.Date,SEneg_dates)]["Perc_Masdar_perday"]
x9=daily_rf[np.isin(daily_rf.Date,Wpos_dates)]["Perc_Masdar_perday"]#["Perc_Masdar_perday"]#
x9_rep=daily_rf[~np.isin(daily_rf.Date,Wpos_dates)]["Perc_Masdar_perday"]
x0=daily_rf[np.isin(daily_rf.Date,Wneg_dates)]["Perc_Masdar_perday"]#["Perc_Masdar_perday"]#
x0_rep=daily_rf[~np.isin(daily_rf.Date,Wneg_dates)]["Perc_Masdar_perday"]

print(mannwhitneyu(x0,x0_rep,alternative='greater'))
x8_rep.median()

x_use = x0
print(x_use.mean(),x_use.median(),x_use.max(),x_use.min())

data=[x5,x6,x7,x8,x9,x0]
fig4, ax4 = plt.subplots()
ax4.set_title('Precipitation Accumulation \nby SST Extremes',fontsize=18)
ax4.boxplot(data, showfliers=False,showmeans=True)
labels=["NE+","NE-","SE+","SE-","W+","W-"]
ax4.set_xticklabels(labels,fontsize=14)
ax4.set_ylabel('Areal Mean Precipitation [mm]',fontsize=12)
plt.tight_layout()
plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/15HrPrecip/BoxWhiskerSST.png",dpi=300)
plt.show()


#Calculating Interquartile Range as a measure of variance:
# q75, q25 = np.percentile(x, [75 ,25])
# iqr = q75 - q25
#Plot distribution:
n1,_,_=plt.hist(x,bins,alpha=1)
plt.xlabel('Daytime Rainfall Accumulation [mm]')
plt.ylabel('Frequency')
# plt.title('Daily Daytime Precipitation \n During Days With {0}'.format(title))
# plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/3.0, 'Med=%(Med)2.2f mm, \nIQR=%(IQR)2.2f mm' %\
#          {'Med':np.median(x),'IQR':iqr})
plt.xlim(0)
# plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/LandSea/DailyPrecip_{0}_Distribution.png".format(cst),
#     dpi=500)
plt.show()



z_statistic, p_value = wilcoxon((x - x.median()))
print(p_value)

######NEED TO GET VALUES FROM NCEP_Rainfall_Amount_LandSea.py FIRST
# y = daily_rf[np.isin(daily_rf["Date"],dates.values)]["Perc_Masdar_perday"]
# z_statistic, p_value = wilcoxon((y - x.median()))
# z_statistic, p_value = ttest_1samp(x,1/999)



# count=11
 

#Standard Deviation Calculation: https://www.strchr.com/standard_deviation_in_one_pass
# sum_squares_2 = np.divide(sum_squares,float(count))
# mean_squared = np.square(mean_accum)
# stdev = np.sqrt(np.abs(sum_squares - mean_squared))
# cor_var = np.divide(stdev,mean_accum)

#Standard Deviation Calculation: http://mathcentral.uregina.ca/QQ/database/QQ.09.06/h/murtaza1.html
# mean_accum = np.divide(accum_precip.data,float(11))
# np.savetxt(NCEPdirpath+"mean_accum.csv",mean_accum,delimiter=',')
# sum_squared = np.square(accum_precip)
# sum_squares_2 = np.multiply(count,sum_squares)
# denom = float(count*(count-1))
# variance = np.divide((sum_squares_2 - sum_squared),denom)
# stdev = np.sqrt(abs(variance))
# cor_var = np.divide(stdev,mean_accum)
# lats = dataset.variables['lats'][:,:]
# lons = dataset.variables['lons'][:,:]
# POP = POP/np.float(sbcount)*100.0
# POP[POP==0]=np.nan

print("Done")
    
    

# #PLOT Mean    
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.background_patch.set_facecolor('k')
# cmap = mpl.cm.gist_ncar
# norm = mpl.colors.Normalize(vmin=0, vmax=mean_accum.max())
# sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,mean_accum.max()))
# sm._A = []
# plt.contourf(lons, lats, mean_accum, 60,
#              transform=ccrs.PlateCarree(), norm=norm, cmap=cmap)
# plt.title("Florida Average Summer \n Daytime Rainfall Accumulation")
# plt.colorbar(sm).set_label("Rainfall Accumulation (mm)")
# states = cfeature.NaturalEarthFeature(category="cultural",
#                     name="admin_1_states_provinces_lines",
#                     scale="10m", facecolor="none")
# coast = cfeature.NaturalEarthFeature(category="physical",
#                     name="coastline",
#                     scale="10m", facecolor="none")
# ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
# ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
# ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
# # plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/15HrPrecip/MeanAccumPrecip.png",
# #     dpi=500)
# plt.show()




# #PLOT Standard Deviation
# ax = plt.axes(projection=ccrs.PlateCarree())
# # cmap = mpl.cm.gist_ncar
# # norm = mpl.colors.Normalize(vmin=0, vmax=15)
# plt.contourf(lons, lats, stdev,# 60,
#              transform=ccrs.PlateCarree())#, norm=norm, cmap=cmap)
# plt.title("Florida Summer 12Hr Rainfall\n  Accumulation Standard Deviation")
# plt.colorbar(sm).set_label("Rainfall Accumulation (mm)")#ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]).set_label("Rainfall Accumulation (mm)")
# # plt.clim(0,100)
# ax.coastlines()

# plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/STDPrecip.png",
#     dpi=500)
# plt.show()



# #PLOT Correlation of Variance
# ax = plt.axes(projection=ccrs.PlateCarree())
# cmap = mpl.cm.gist_ncar
# norm = mpl.colors.Normalize(vmin=0, vmax=15)
# plt.contourf(lons, lats, cor_var,# 60,
#              transform=ccrs.PlateCarree())#, norm=norm, cmap=cmap)
# plt.title("Florida 12Hr Rainfall Accumulation \nJJA Correlation of Variance ")
# plt.colorbar(sm).set_label("Rainfall Accumulation (mm)")#ticks=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]).set_label("Rainfall Accumulation (mm)")


# plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/CVPrecip.png",
#     dpi=500)
# plt.show()



# #PLOT POP    #PLOT PoP
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.background_patch.set_facecolor('k')
# cmap = mpl.cm.gist_ncar
# norm = mpl.colors.Normalize(vmin=0, vmax=100)
# sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,100))
# sm._A = []
# plt.contourf(lons, lats, POP, 100,
#              transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
             
# plt.title("JJA 15Hr Daytime\nRainfall Frequency")
# plt.colorbar(sm,ticks=[0,10,20,30,40,50,60,70,80,90,100]).set_label("PoP (%)")
# states = cfeature.NaturalEarthFeature(category="cultural",
#                     name="admin_1_states_provinces_lines",
#                     scale="10m", facecolor="none")
# coast = cfeature.NaturalEarthFeature(category="physical",
#                     name="coastline",
#                     scale="10m", facecolor="none")
# ax.add_feature(states, linestyle="-", edgecolor="darkslategray", linewidth=1)
# ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
# reader = shpreader.Reader("/Volumes/LaCie/DPM_Code/Python/countyl010g_shp_nt00964/countyl010g.shp")
# counties = list(reader.geometries())
# COUNTIES = cfeature.ShapelyFeature(counties,ccrs.PlateCarree())
# ax.add_feature(coast, linestyle="-", edgecolor="darkslategray", linewidth=1)
# ax.set(xlim=(-83.75,-79.75),ylim=(24.4,30.5))
# plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/Figures/JJA15HrPoP_New.png", 
#     dpi=500)
# plt.show()


