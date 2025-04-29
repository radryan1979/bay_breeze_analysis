"""
The purpose of this code is to calculate the probability
of sea breeze by regime and create a table of probability of
each type of sea breeze at each station by regime.

Coded by Dan Moore

Updated: 3/22/19
"""

import pandas as pd
import numpy as np
from datetime import datetime,timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import os

from cartopy import config
import cartopy.crs as ccrs


times_dir = '/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/'

regime_infile = "/Volumes/LaCie/SeaBreeze/Florida/TampaRegimeTyping/TampaTypingUpdate.csv"


times_infile = times_dir+'Hughes_SST_Regime.csv'
station_df = pd.read_csv(times_infile,header=0,index_col=0)
reg_df = pd.read_csv(regime_infile,header=0,index_col=0)

classic = station_df[station_df["Type"]=="classic"]
cdp = station_df[station_df["Type"]=="classicDP"]
wk = station_df[station_df["Type"]=="weak"]
dpws = station_df[station_df["Type"]=="DPWS"]
cws = station_df[station_df["Type"]=="classicWS"]

df = pd.DataFrame([])

cnt=0

for st in station_df.index.drop_duplicates():

    classic_count = classic.loc[st]["Regime"].groupby(classic.loc[st]["Regime"]).count()
    sb_count = station_df.loc[st]["Regime"].groupby(station_df.loc[st]["Regime"]).count()
    wk_count = wk.loc[st]["Regime"].groupby(wk.loc[st]["Regime"]).count()
    cdp_count = cdp.loc[st]["Regime"].groupby(cdp.loc[st]["Regime"]).count()
    dpws_count = dpws.loc[st]["Regime"].groupby(dpws.loc[st]["Regime"]).count()
    cws_count = cws.loc[st]["Regime"].groupby(cws.loc[st]["Regime"]).count()
    
    temp_df = pd.DataFrame({"{0} SB Count".format(st): sb_count,
                            "{0} Classic Count".format(st): classic_count,
                            "{0} Classic Count".format(st): wk_count,
                            "{0} Classic Count".format(st): cdp_count,
                            "{0} Classic Count".format(st): dpws_count,
                            "{0} Classic Count".format(st): cws_count
                            })
    if cnt==0:
        df=temp_df.copy()
    else:
        df = pd.concat([df,temp_df],axis=1)
    
    
    cnt+=1


regime_count = reg_df["Regime"].groupby(reg_df["Regime"]).count()

final_df = df.div(regime_count, axis='index')*100

final_df = final_df.T

final_df.to_csv("/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/SB_ProbabilityByRegime.csv")