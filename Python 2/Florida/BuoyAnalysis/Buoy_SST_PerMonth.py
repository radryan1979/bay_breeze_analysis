"""
Coded by: Dan Moore

This program will ingest Buoy Data.

We are going to average SSTs by month to compare with
Teleconnection indices and monthly precipitation or 
sea breeze patterns/frequencies.

Updated: 2-6-19
"""


import pandas as pd
import numpy as np
from datetime import datetime,timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import os


# BUOYdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/"
BUOYdirpath = "/Users/dpmoore2927/Desktop/"

westSSTpath = BUOYdirpath + '42013COMBINED.csv'
eastSSTpath = BUOYdirpath + '41009COMBINED.csv'
southSST_insidepath = BUOYdirpath + 'PKYF1COMBINED.csv'
southSST_outsidepath = BUOYdirpath + 'MLRF1COMBINED.csv'

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


westSST_monthly = westSST_df.groupby(["#YY","MM"])["WTMP"].mean()
westSST_monthly = westSST_monthly.reset_index()
westSST_monthly = westSST_monthly[(westSST_monthly["MM"]>3) & (westSST_monthly["MM"]<11)]
westSST_monthly = westSST_monthly.reset_index(drop=True)

eastSST_monthly = eastSST_df.groupby(["#YY","MM"])["WTMP"].mean()
eastSST_monthly = eastSST_monthly.reset_index()
eastSST_monthly = eastSST_monthly[(eastSST_monthly["MM"]>3) & (eastSST_monthly["MM"]<11)]
eastSST_monthly = eastSST_monthly.reset_index(drop=True)

southSSTinside_monthly = southSSTinside_df.groupby(["#YY","MM"])["WTMP"].mean()
southSSTinside_monthly = southSSTinside_monthly.reset_index()
southSSTinside_monthly = southSSTinside_monthly[(southSSTinside_monthly["MM"]>3) & (southSSTinside_monthly["MM"]<11)]
southSSTinside_monthly = southSSTinside_monthly.reset_index(drop=True)

southSSToutside_monthly = southSSToutside_df.groupby(["#YY","MM"])["WTMP"].mean()
southSSToutside_monthly = southSSToutside_monthly.reset_index()
southSSToutside_monthly = southSSToutside_monthly[(southSSToutside_monthly["MM"]>3) & (southSSToutside_monthly["MM"]<11)]
southSSToutside_monthly = southSSToutside_monthly.reset_index(drop=True)

# total_df = pd.DataFrame({   "WestSST":          westSST_monthly["WTMP"],
#                             "EastSST":          eastSST_monthly["WTMP"],
#                             "SouthSSTinside":   southSSTinside_monthly["WTMP"],
#                             "SouthSSToutside":  southSSToutside_monthly["WTMP"]
#                         })

eastSST_monthly.to_csv("/Users/dpmoore2927/Desktop/EastSST_byMonth.csv")
southSSTinside_monthly.to_csv("/Users/dpmoore2927/Desktop/SouthSSTInside_byMonth.csv")
southSSToutside_monthly.to_csv("/Users/dpmoore2927/Desktop/SouthSSTOutside_byMonth.csv")

