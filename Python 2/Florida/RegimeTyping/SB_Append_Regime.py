"""
Coded by: Dan Moore

This program will ingest a list of dates from the Hughes
detection algorithm in Florida, and a list of Regime Data.

We are going to append the Synoptic Regime Type 
as new column to each detection. Then we can use UFLHistoEast.py
and UFLHistoWest.py to analyze the statistics of each by month/year.
We can also analyze 

Updated: 2-5-19
"""


import pandas as pd
import numpy as np
from datetime import datetime,timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import os


times_dir = '/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/'
# coasts = ['West']

REGdirpath = "/Volumes/LaCie/SeaBreeze/Florida/TampaRegimeTyping/"

REGpath = REGdirpath + 'TampaTypingUpdate.csv'

Reg_df = pd.read_csv(REGpath,header=0,index_col=0)


times_infile = times_dir+'Hughes_SST.csv'
station_df = pd.read_csv(times_infile,header=0,index_col=[0,1])

for index,row in station_df.iterrows():
    
    SBF_hour = str(int(row['Hour']/1))
    SBF_min = str(int(row['Hour']%1*60)).zfill(2)
    date_str = str(index[1])+' '+SBF_hour+':'+SBF_min
    SBF_time = datetime.strptime(date_str,'%Y-%m-%d %H:%M')
    
    if len(Reg_df.loc[(Reg_df["Year"]==SBF_time.year) & \
        (Reg_df["Month"]==SBF_time.month) & (Reg_df["Day"]==SBF_time.day)])==0:
            station_df.loc[index,'Regime'] = "NaN"
    else:
        station_df.loc[index,'Regime'] = Reg_df.loc[(Reg_df["Year"]==SBF_time.year) & \
            (Reg_df["Month"]==SBF_time.month) & (Reg_df["Day"]==SBF_time.day)]["Regime"].values[0]
        

station_df.to_csv(times_dir+"Hughes_SST_Regime.csv")
    
    
    

