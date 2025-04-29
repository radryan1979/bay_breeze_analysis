"""
Coded by: Daniel Moore

The purpose of this program is to ingest the final dataset files and just spit out the dates of each SB event captured by each station and output a pandas dataframe file of just the dates.
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
    path_name="/Volumes/LaCie/SeaBreeze/FinalDataset/"
    suffix="final.csv"
    infile=open(os.path.join\
    (path_name,station+suffix),"r")

    #Read data
    obs_list=infile.readlines()
    n=len(obs_list)

    Time=np.empty(n,dtype=d.datetime)

    i=0;j=0

    for obs in obs_list[1:]:

        Time[i]=obs.split(",")[0]

        i+=1

    numdata=len(Time)

    dates=np.empty(n,dtype=d.datetime)

    i=36

    while i < numdata:

        dates[j]=Time[i]
        j+=1

        i+=72

    numSB=np.count_nonzero(dates)
        #Number of nonzero rows

    dates=np.resize(dates,numSB)
        #Trims trailing blank rows


    df=pd.DataFrame(dates)

    #print(df)


    df.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/FinalDates/"\
                              ,station+"Dates.csv"))




