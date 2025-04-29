"""
Coded by: Dan Moore

Purpose: Create climatologies of variables by regime.
For instance, I will plot the average relative humidity
difference between a test and reference station
by regime to determine effect of sea breeze when wind
directions aren't sufficient to detect presence.

Updated: 2-18-19
"""

import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt
import sys

#Which stations do you want to look at?
test_st         =       340
ref_st          =       390

#Which Regime, if applicable, do you want to study?
reg             =       1

#Individual Month Analysis?
mo_analysis     =       True        #Set to true if you want to look at data by month
mo              =       4           #Set to month of desired analysis

#This dataset contains met data for entire study period
met_path_name       =       "/Volumes/LaCie/SeaBreeze/Florida/UFLFormattedData/"
type_path           =       "/Volumes/LaCie/SeaBreeze/Florida/TampaRegimeTyping/TampaTypingUpdate.csv"
out_path            =       "/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results1_17_19/Figures/"

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#Comparing two stations for the entire study period or by month, not just SBC days.

# test_df                         =       pd.read_csv(met_path_name+"{0}_formatted.csv".format(test_st),header=0,index_col=None)
# test_df["local_eastern_time"]   =       test_df["local_eastern_time"].astype("datetime64")
# ref_df                          =       pd.read_csv(met_path_name+"{0}_formatted.csv".format(ref_st),header=0,index_col=None)
# ref_df["local_eastern_time"]    =       ref_df["local_eastern_time"].astype("datetime64")


# if mo_analysis:
#     summer_test                 =       test_df[(test_df["local_eastern_time"].dt.month==mo)]
# else:
#     summer_test                 =       test_df[(test_df["local_eastern_time"].dt.month>3) & (test_df["local_eastern_time"].dt.month<11)]
# summer_hrly_wd_test             =       summer_test["wind_direction_10m_deg"].groupby(summer_test["local_eastern_time"].dt.hour).mean()
# summer_hrly_ws_test             =       summer_test["wind_speed_10m_mph"].groupby(summer_test["local_eastern_time"].dt.hour).mean()

# if mo_analysis:
#     summer_ref              =       ref_df[(ref_df["local_eastern_time"].dt.month==mo)]
# else:
#     summer_ref              =       ref_df[(ref_df["local_eastern_time"].dt.month>3) & (ref_df["local_eastern_time"].dt.month<11)]
# summer_hrly_wd_ref      =       summer_ref["wind_direction_10m_deg"].groupby(summer_ref["local_eastern_time"].dt.hour).mean()
# summer_hrly_ws_ref      =       summer_ref["wind_speed_10m_mph"].groupby(summer_ref["local_eastern_time"].dt.hour).mean()

# fig, ax = plt.subplots(1,1, figsize=(15,4))
# ax.quiver(summer_hrly_wd_ref.index, np.ones(len(summer_hrly_wd_ref)) * 1.5, 
#                                     summer_hrly_ws_ref.values*np.cos((90-summer_hrly_wd_ref.values)*np.pi/180.0)*-1.0, 
#                                     summer_hrly_ws_ref.values*np.sin((90-summer_hrly_wd_ref.values)*np.pi/180.0)*-1.0,
#                                     color='b', headwidth=4, headlength=4, width=0.0025, label="Reference")
# ax.quiver(summer_hrly_wd_test.index, np.ones(len(summer_hrly_wd_test)) * 1.5, 
#                                     summer_hrly_ws_test.values*np.cos((90-summer_hrly_wd_test.values)*np.pi/180.0)*-1.0, 
#                                     summer_hrly_ws_test.values*np.sin((90-summer_hrly_wd_test.values)*np.pi/180.0)*-1.0,
#                                     color='r', headwidth=4, headlength=4, width=0.0025, label="Test")
# plt.yticks([])
# plt.xticks(range(24))
# plt.xlabel("Eastern Time")
# plt.legend()
# if mo_analysis:
#     plt.title("Florida Wind Climatology for \nTest Station: {0} and Reference Station: {1} \nMonth: {2}".format(test_st, ref_st, mo))
#     # plt.savefig(out_path+"{0}vs{1}MonthlyWind_{2}.png".format(test_st,ref_st,mo),dpi=500)
# else:
#     plt.title("Florida Summer Wind Climatology for \nTest Station: {0} and Reference Station: {1}".format(test_st, ref_st))
#     # plt.savefig(out_path+"{0}vs{1}SummerWind.png".format(test_st,ref_st),dpi=500)
# plt.show()
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#Comparing days where SB occured to reference climatology

# #This dataset contains dates of SB detection
# sb_path_name       =       "/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results1_17_19/"

# if test_st in [350,360,380,450,480,490]:
#     coast="West"
# elif test_st in [340,410,420,440]:
#     coast="East"
# else:
#     print("Invalid Test Station.")
#     sys.exit()

# #Dataframe to contain dates of SBC detected at given test station
# sb_df                           =       pd.read_csv(sb_path_name+"Hughes{0}_SST_Regime.csv".format(coast),header=0,index_col=False)
# sb_df                           =       sb_df[sb_df["Station"]==test_st]
# sb_df["Date"]                   =       sb_df["Date"].astype("datetime64")
# sb_df                           =       sb_df["Date"].astype(str)


# #Find dates listed above in formatted dataset of meteorological data
# test_df                         =       pd.read_csv(met_path_name+"{0}_formatted.csv".format(test_st),header=0,index_col=None)
# test_df["local_eastern_time"]   =       test_df["local_eastern_time"].astype(str)
# ref_df                          =       pd.read_csv(met_path_name+"{0}_formatted.csv".format(ref_st),header=0,index_col=None)
# ref_df["local_eastern_time"]    =       ref_df["local_eastern_time"].astype("datetime64")

# #Test Station
# summer_test                         =       test_df[test_df["local_eastern_time"].astype(str).str[0:10].isin(sb_df)]
# summer_test["local_eastern_time"]   =       summer_test["local_eastern_time"].astype("datetime64")
# if mo_analysis:
#     summer_test                     =       summer_test[(summer_test["local_eastern_time"].dt.month==mo)]
# else:
#     pass

# #Reference Station
# summer_hrly_wd_test                 =       summer_test["wind_direction_10m_deg"].groupby(summer_test["local_eastern_time"].dt.hour).mean()
# summer_hrly_ws_test                 =       summer_test["wind_speed_10m_mph"].groupby(summer_test["local_eastern_time"].dt.hour).mean()
# if mo_analysis:
#     summer_ref              =       ref_df[(ref_df["local_eastern_time"].dt.month==mo)]
# else:
#     summer_ref              =       ref_df[(ref_df["local_eastern_time"].dt.month>3) & (ref_df["local_eastern_time"].dt.month<11)]
    
# summer_hrly_wd_ref      =       summer_ref["wind_direction_10m_deg"].groupby(summer_ref["local_eastern_time"].dt.hour).mean()
# summer_hrly_ws_ref      =       summer_ref["wind_speed_10m_mph"].groupby(summer_ref["local_eastern_time"].dt.hour).mean()


# #Plot
# fig, ax = plt.subplots(1,1, figsize=(15,4))
# ax.quiver(summer_hrly_wd_ref.index, np.ones(len(summer_hrly_wd_ref)) * 3.5, 
#                                     summer_hrly_ws_ref.values*np.cos((90-summer_hrly_wd_ref.values)*np.pi/180.0)*-1.0, 
#                                     summer_hrly_ws_ref.values*np.sin((90-summer_hrly_wd_ref.values)*np.pi/180.0)*-1.0,
#                                     color='b', headwidth=4, headlength=4, width=0.0025, label="Reference")
# ax.quiver(summer_hrly_wd_test.index, np.ones(len(summer_hrly_wd_test)) * 3.5, 
#                                     summer_hrly_ws_test.values*np.cos((90-summer_hrly_wd_test.values)*np.pi/180.0)*-1.0, 
#                                     summer_hrly_ws_test.values*np.sin((90-summer_hrly_wd_test.values)*np.pi/180.0)*-1.0,
#                                     color='r', headwidth=4, headlength=4, width=0.0025, label="Test")

# plt.yticks([])
# plt.xticks(range(24))
# plt.xlabel("Eastern Time")
# plt.legend()
# plt.title("Wind Climatology for Two Florida Stations: Test Station Limited to SBC Days")
# plt.show()

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#Comparing days wind climatologies for certain regimes
test_df                         =       pd.read_csv(met_path_name+"{0}_formatted.csv".format(test_st),header=0,index_col=None)
test_df["local_eastern_time"]   =       test_df["local_eastern_time"].astype("datetime64")
ref_df                          =       pd.read_csv(met_path_name+"{0}_formatted.csv".format(ref_st),header=0,index_col=None)
ref_df["local_eastern_time"]    =       ref_df["local_eastern_time"].astype("datetime64")
type_df                         =       pd.read_csv(type_path,header=0,index_col=0)

#Adding new column to type dataframe that is type datetime64
for index,row in type_df.iterrows():
    type_df.loc[index,"Datetime"]       =       str(row["Month"])+"-"+str(row["Day"])+"-"+str(int(row["Year"]))

type_df["Datetime"]             =       type_df["Datetime"].astype("datetime64")
type_df["Datetime"]             =       type_df["Datetime"].astype(str)

for reg in range(1,10):
    temp_type_df                    =       type_df[type_df["Regime"]==reg]
    
    
    regime_test                     =       test_df[test_df["local_eastern_time"].astype(str).str[0:10].isin(temp_type_df["Datetime"])]
    regime_ref                      =       ref_df[ref_df["local_eastern_time"].astype(str).str[0:10].isin(temp_type_df["Datetime"])]
    
    #Eliminate erroneous data points.
    regime_ref.loc[regime_ref["temp_dp_2m_C"] <-0.0,"temp_dp_2m_C"] = np.nan
    regime_test.loc[regime_test["temp_dp_2m_C"] <-0.0,"temp_dp_2m_C"] = np.nan
    
    #Create mask to eliminate instances of 'NaN'
    mask                            =       np.isfinite(regime_test["rh_2m_pct"]) & np.isfinite(regime_ref["rh_2m_pct"])
    regime_test["RH_Diff"]          =       regime_test["rh_2m_pct"].loc[mask]-regime_ref["rh_2m_pct"].loc[mask]
    regime_rh_diff                  =       regime_test["RH_Diff"].groupby(regime_test["local_eastern_time"].dt.hour).mean()
    
    
    
    
    
    
    
    regime_hrly_wd_test             =       regime_test["wind_direction_10m_deg"].groupby(regime_test["local_eastern_time"].dt.hour).mean()
    regime_hrly_ws_test             =       regime_test["wind_speed_10m_mph"].groupby(regime_test["local_eastern_time"].dt.hour).mean()
    
    regime_hrly_wd_ref      =       regime_ref["wind_direction_10m_deg"].groupby(regime_ref["local_eastern_time"].dt.hour).mean()
    regime_hrly_ws_ref      =       regime_ref["wind_speed_10m_mph"].groupby(regime_ref["local_eastern_time"].dt.hour).mean()
    
    fig, ax = plt.subplots(1,1, figsize=(15,4))
    
    plt.plot(regime_rh_diff, color='g', label="RH Difference")
    
    ax.quiver(regime_hrly_wd_ref.index, np.ones(len(regime_hrly_wd_ref)) * regime_rh_diff.mean(), 
                                        regime_hrly_ws_ref.values*np.cos((90-regime_hrly_wd_ref.values)*np.pi/180.0)*-1.0, 
                                        regime_hrly_ws_ref.values*np.sin((90-regime_hrly_wd_ref.values)*np.pi/180.0)*-1.0,
                                        color='b', headwidth=4, headlength=4, width=0.0025, label="Reference Wind")
    ax.quiver(regime_hrly_wd_test.index, np.ones(len(regime_hrly_wd_test)) * regime_rh_diff.mean(), 
                                        regime_hrly_ws_test.values*np.cos((90-regime_hrly_wd_test.values)*np.pi/180.0)*-1.0, 
                                        regime_hrly_ws_test.values*np.sin((90-regime_hrly_wd_test.values)*np.pi/180.0)*-1.0,
                                        color='r', headwidth=4, headlength=4, width=0.0025, label="Test Wind")
    
    plt.xticks(range(24))
    plt.xlabel("Eastern Time")
    plt.legend()
    plt.ylabel("Test-Ref RH (%)")
    plt.title("Florida Regime {0} Test-Reference Dew Point Temp Climatology \n Test Station: {1} and Reference Station: {2}".format(reg, test_st, ref_st))
    plt.savefig(out_path+"RH_Difference/{0}/{0}vs{1}Regime{2}RH.png".format(test_st,ref_st,reg),dpi=500)
    plt.show()
