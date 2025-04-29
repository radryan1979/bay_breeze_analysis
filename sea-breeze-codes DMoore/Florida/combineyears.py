"""
Coded by: Dan Moore

The purpose of this module is to take all separate files of each station per year, and combine into one csv per station

"""

import numpy as np
import math as m
import pandas as pd
import os.path
import datetime as d





Station_list=[110]
"""
,120,121,130,140,150,160,170,180,230,240,241,250,\
260,270,275,280,290,302,303,304,311,320,330,340,350,360,371,380,390,405,\
410,420,425,435,440,450,455,460,470,480,490]
"""
Year_list=[2013,2014,2015,2016,2017]




for station in Station_list:

    for year in Year_list:

        path_name="/Volumes/LaCie/SeaBreeze/UFLFormattedData/"
        suffix=".csv"
        infile=open(os.path.join\
        (path_name,str(station)+'_'+str(year)+suffix),"r")

        if year==2013:
            df=pd.read_csv(infile,sep=',',header=0,index_col=0)
        else:
            df_temp=pd.read_csv(infile,sep=',',header=0,index_col=0)
            df=df.append(df_temp)

    df.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/UFLFormattedData/"+\
        str(station)+".csv"))
        #Prints dataframe to csv.
