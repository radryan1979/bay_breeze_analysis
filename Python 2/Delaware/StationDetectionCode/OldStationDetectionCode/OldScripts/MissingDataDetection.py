"""
Coded by: Dan Moore

The purpose of this program is to read through the data and detect missing data and fill those gaps with 'nan'.


*********
Abandoned because completed all in FormatData.py.
*********
"""

import numpy as np
import math as m
import pandas as pd
import SBFilters
import os.path
import datetime as d

#Allows csv to save unlimited amount of rows
np.set_printoptions(threshold=np.inf)

numsecs=(365*5+1)*24*60*60
#number of seconds in the 5 year range

base = d.datetime(2013, 1, 1)#start date
Time = np.array([base + d.timedelta(seconds=i) \
                  for i in range(0,numsecs,300)])
#implicit loop to skip 300 seconds or 5 minutes

Station_list=["DRHB"]
"""
,"DBLK","DSMY","DADV","DBBB","DDFS","DELN",\
"DHAR","DBNG","DGES","DGUM","DIRL","DJCR","DMIL","DSJR","DSND",\
"DWAR","DSEA","DLAU","DBRG"]
"""

maxdata=len(Time)

for station in Station_list:

    path_name="/Volumes/LaCie/SeaBreeze/FormattedData/"
    suffix=".csv"
    infileDEOS=open(os.path.join\
    (path_name+station+suffix),"r")

    df=pd.read_csv(infileDEOS,sep=",",header=0,index_col=0)

    df=df.reindex(d.datetime(Time))

    #Print to csv
    df.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/CompleteData/"+\
            str(station)+".csv"),sep=",")

