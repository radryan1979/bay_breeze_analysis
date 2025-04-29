"""
Coded by: Dan Moore

This program will ingest a list of dates from the Hughes
detection algorithm in Florida, and a list of SST Indices.

We will append West, SE and NE coastline SST anomaly indices
for each SB case.

Update 5/24/19: Analyzing 2 and 3 standard deviations.

Updated: 5/24/19
"""


import pandas as pd
import numpy as np
from datetime import datetime,timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import os


times_dir = '/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/'

BUOYdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/SST_Indices/"

W_SSTpath = BUOYdirpath + '42013_SST_Index_2STD.csv'
NE_SSTpath = BUOYdirpath + '41113_SST_Index_2STD.csv'
SE_SSTpath = BUOYdirpath + 'MLRF1_SST_Index_2STD.csv'

W_SST_df = pd.read_csv(W_SSTpath,header=0,index_col=0, parse_dates=True)
NE_SST_df = pd.read_csv(NE_SSTpath,header=0,index_col=0, parse_dates=True)
SE_SST_df = pd.read_csv(SE_SSTpath,header=0,index_col=0, parse_dates=True)

times_infile = times_dir+'Hughes.csv'
station_df = pd.read_csv(times_infile,header=0,index_col=[0,1], parse_dates=True)

for index,row in station_df.iterrows():
    
    date_str = str(index[1])
    # SBF_time = datetime.strptime(date_str,'%d-%b-%Y')
    # date_str = str(SBF_time)

    if pd.isna(NE_SST_df.loc[date_str]["Index"]):
        station_df.loc[index,'SST_Anomaly_NE'] = "NaN"
    else:
        station_df.loc[index,'SST_Anomaly_NE'] = NE_SST_df.loc[date_str]["Index"]
    if pd.isna(SE_SST_df.loc[date_str]["Index"]):
        station_df.loc[index,'SST_Anomaly_SE'] = "NaN"
    else:
        station_df.loc[index,'SST_Anomaly_SE'] = SE_SST_df.loc[date_str]["Index"]
    if pd.isna(W_SST_df.loc[date_str]["Index"]):
        station_df.loc[index,'SST_Anomaly_W'] = "NaN"
    else:
        station_df.loc[index,'SST_Anomaly_W'] = W_SST_df.loc[date_str]["Index"]
        

station_df.to_csv(times_dir+"Hughes_SST_2STD.csv")
    
    
    

