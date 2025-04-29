"""
Coded by: Daniel Moore

The purpose of this program is to ingest a file and spit out Sea Breeze instances detected using the DetAlg module with either the 15min data.

"""

import numpy as np
import math as m
import pandas as pd
import os.path
import datetime as d
import csv

Station_list=[110]
"""
,120,121,130,140,150,160,170,180,230,240,241,250,\
260,270,275,280,290,302,303,304,311,320,330,340,350,360,371,380,390,405,\
410,420,425,435,440,450,455,460,470,480,490]
"""
maxdata=(365*5+1)*24*4

numsecs=(365*5+1)*24*60*60
#number of seconds in the 5 year range

base = d.datetime(2013, 1, 1)#start date
Time = np.array([base + d.timedelta(seconds=i) \
                  for i in range(0,numsecs,900)])
#implicit loop to skip 900 seconds or 15 minutes


for station in Station_list: #Loop through stations
    path_name="/Volumes/LaCie/SeaBreeze/UFLFormattedData/"
    suffix=".csv"
    infile=open(os.path.join\
    (path_name,str(station)+suffix),"r")

    df=pd.read_csv(infile,sep=',',header=0,index_col=1)

    df.index=pd.to_datetime(df.index)

    df = df.loc[~df.index.duplicated(keep='first')]

    df=df.drop(["StationID"],axis=1)
        #erasing StationID column as it's in the filename

    #print(df.head(10))

    df=df.reindex(index=Time,fill_value='NaN')
        #Reindex data frame with all inclusive Time_stamp to
        #eliminate missing dates in data - fills blank spaces with
        #'NaN'

    df.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/UFLFormattedData/"+\
        str(station)+"_formatted.csv"))
        #Prints dataframe to csv.




