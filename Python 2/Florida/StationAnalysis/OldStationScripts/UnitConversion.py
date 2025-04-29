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

for station in Station_list: #Loop through years of interest
    path_name="/Volumes/LaCie/SeaBreeze/UFLFormattedData/"
    suffix=".csv"
    infile=open(os.path.join\
    (path_name,str(station)+suffix),"r")

    df=pd.read_csv(infile,sep=',',header=0,index_col=0)

    df['rain_2m_inches']=df['rain_2m_inches']*25.4
    df.rename(columns={'rain_2m_inches':'rain_2m_mm'},inplace=True)

    df.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/UFLFormattedData/"+\
        str(station)+".csv"))
        #Prints dataframe to csv.




