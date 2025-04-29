import numpy as np
import math as m
import pandas as pd
import SBFilters
import os.path
import datetime as d
station="DRHB"

year=2013

path_name="/Volumes/LaCie/SeaBreeze/OriginalData/StationData/"
suffix=".dat"
infileDEOS=open(os.path.join\
(path_name,str(year),station+str(year)+suffix),"r")

df=pd.read_csv(infileDEOS,sep=',',header=0,index_col=False)

df["Time_stamp"] = df["Date"].map(str) + " " + df[" Time"].map(str)

df["Time_stamp"]=df["Time_stamp"].str[:-10]

df=df.drop(columns=["Date"," Time"])

df=df.set_index("Time_stamp")

df=df.iloc[::-1]

df1=df1.append(df)

path_name="/Volumes/LaCie/SeaBreeze/OriginalData/StationData/"
suffix=".dat"
infileDEOS1=open(os.path.join\
(path_name,str(2014),station+str(2014)+suffix),"r")

df_temp=pd.read_csv(infileDEOS1,sep=',',header=0,index_col=False)

df_temp["Time_stamp"] = df_temp["Date"].map(str) + " " + df_temp[" Time"].map(str)

df_temp["Time_stamp"]=df_temp["Time_stamp"].str[:-10]

df_temp=df_temp.drop(columns=["Date"," Time"])

df_temp=df_temp.set_index("Time_stamp")

df_temp=df_temp.iloc[::-1]

df1=df1.append(df_temp)

df1.to_csv(os.path.join("/Users/dpmoore2927/Desktop/"+\
        str(station)+".csv"))
