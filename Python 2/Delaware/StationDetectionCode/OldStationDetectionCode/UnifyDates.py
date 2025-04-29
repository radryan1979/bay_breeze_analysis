"""
Coded by: Daniel Moore

The purpose of this program is to ingest the final dataset files and run statistics on each event by station
"""

import numpy as np
import math as m
import pandas as pd
import os.path
import datetime as d
import csv
from copy import deepcopy

Station_list=["DRHB","DSMY","DADV","DBBB","DDFS",\
"DBNG","DGES","DGUM","DIRL","DJCR","DMIL","DSND",\
"DSEA","DLAU","DBRG"]

#Expected number of data points for all five years
maxdata=(365*5+1)*24*12


for station in Station_list: #Loop through years of interest
    path_name="/Volumes/LaCie/SeaBreeze/FinalDates/"
    suffix="Dates.csv"
    infile=open(os.path.join\
    (path_name,station+suffix),"r")

    dfdates=pd.read_csv(infile,header=0,index_col=0)

    dates=dfdates['0'].tolist()

    numSB=len(dates)

    for i in range(numSB):
        print(station,dates)
        """

        sbtime=d.datetime.strptime(dates[i], '%m/%d/%y %H:%M')
        sbtime=sbtime.strftime('%Y-%m-%d %H:%M:%S')
        dates[i]=sbtime
        """


    df=pd.DataFrame(dates)
    df.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/FinalDates/"\
                          ,station+"Dates.csv"))






