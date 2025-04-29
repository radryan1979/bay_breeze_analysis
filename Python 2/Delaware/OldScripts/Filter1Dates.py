
import numpy as np
import math as m
import pandas as pd
import os.path
import datetime as d
import csv

count=0
myset=set()

Station_list=["DRHB","DSMY","DADV","DBBB","DDFS",\
"DBNG","DGES","DGUM","DIRL","DJCR","DMIL","DSND",\
"DSEA","DLAU","DBRG"]

for station in Station_list: #Loop through years of interest
    path_name="/Volumes/LaCie/SeaBreeze/Delaware/Filter1/"
    suffix=".csv"
    infile=path_name+station+suffix
    
    df=pd.read_csv(infile,header=None)
    
    for i in range(len(df.columns)):
        date=df[i].values[0][0:10]
        if int(date[0:4])>2017 or int(date[0:4])<2016:
            continue
        elif int(date[5:7])>9 or int(date[5:7])<5:
            continue
        elif date in myset:
            continue
        else:
            count+=1
            myset.add(date)
    
    print(station,count)
    
        