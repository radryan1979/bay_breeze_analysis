"""
Coded by: Dan Moore

This program will ingest a list of dates from the Hughes
detection algorithm in Florida, and a list of Buoy Data.

We are going to append West, East and South Coast SST 
as new columns to each detection. Then we can use UFLHistoEast.py
and UFLHistoWest.pyto analyze the statistics of each by month/year. 
We can also use a new program, FL_SST_Histograms.py to characterize 
climatologies for each month by coast.

Updated: 1-28-19
"""


import pandas as pd
import numpy as np
from datetime import datetime,timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import os


times_dir = '/Volumes/LaCie/SeaBreeze/Florida/Results1_17_19/'
coasts = ['West']

BUOYdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/"

westSSTpath = BUOYdirpath + '42013_WestFL/42013COMBINED.csv'
eastSSTpath = BUOYdirpath + '41009_EastFL/41009COMBINED.csv'
southSST_insidepath = BUOYdirpath + 'PKYF1_SouthFL/PKYF1COMBINED.csv'
southSST_outsidepath = BUOYdirpath + 'MLRF1_SouthFL/MLRF1COMBINED.csv'

westSST_df = pd.read_csv(westSSTpath,header=0, skiprows=[1],index_col=0)
westSST_df = westSST_df.reset_index(drop=True)
westSST_df = westSST_df.replace(999,np.NaN)
eastSST_df = pd.read_csv(eastSSTpath,header=0,skiprows=[1],index_col=0)
eastSST_df = eastSST_df.reset_index(drop=True)
eastSST_df = eastSST_df.replace(999,np.NaN)
southSSTinside_df = pd.read_csv(southSST_insidepath,header=0,skiprows=[1],index_col=0)
southSSTinside_df = southSSTinside_df.reset_index(drop=True)
southSSTinside_df = southSSTinside_df.replace(999,np.NaN)
southSSToutside_df = pd.read_csv(southSST_outsidepath,header=0,skiprows=[1],index_col=0)
southSSToutside_df = southSSToutside_df.reset_index(drop=True)
southSSToutside_df = southSSToutside_df.replace(999,np.NaN)


for coast in coasts:
    times_infile = times_dir+'Hughes{0}.csv'.format(coast)
    station_df = pd.read_csv(times_infile,header=0,index_col=[0,1])
    
    for index,row in station_df.iterrows():
        
        #No Data for 2008 or 2018 in Buoy data
        if index[1][7:11]=='2008' or index[1][7:11]=='2018':
            continue
        
        SBF_hour = str(int(row['Hour']/1))
        SBF_min = str(int(row['Hour']%1*60))
        date_str = str(index[1])+' '+SBF_hour+':'+SBF_min
        SBF_time = datetime.strptime(date_str,'%d-%b-%Y %H:%M')

        station_df.loc[index,'SST_East'] = np.nanmean(eastSST_df.loc[(eastSST_df["#YY"]==SBF_time.year) & \
            (eastSST_df["MM"]==SBF_time.month) & (eastSST_df["hh"]==SBF_time.hour) & \
            (eastSST_df["DD"]==SBF_time.day)]["WTMP"])
        station_df.loc[index,'SST_West'] = np.nanmean(westSST_df.loc[(westSST_df["#YY"]==SBF_time.year) & \
            (westSST_df["MM"]==SBF_time.month) & (westSST_df["hh"]==SBF_time.hour) & \
            (westSST_df["DD"]==SBF_time.day)]["WTMP"])
        station_df.loc[index,'SST_South_InsideKeys'] = np.nanmean(southSSTinside_df.loc[(southSSTinside_df["#YY"]==SBF_time.year) & \
            (southSSTinside_df["MM"]==SBF_time.month) & (southSSTinside_df["hh"]==SBF_time.hour) & \
            (southSSTinside_df["DD"]==SBF_time.day)]["WTMP"])
        station_df.loc[index,'SST_South_OutsideKeys'] = np.nanmean(southSSToutside_df.loc[(southSSToutside_df["#YY"]==SBF_time.year) & \
            (southSSToutside_df["MM"]==SBF_time.month) & (southSSToutside_df["hh"]==SBF_time.hour) & \
            (southSSToutside_df["DD"]==SBF_time.day)]["WTMP"])
        

station_df.to_csv(times_dir+"Hughes{0}_SST.csv".format(coast))
    
    
    

