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
times_infile = times_dir+'Hughes.csv'
station_df = pd.read_csv(times_infile,header=0,index_col=0)




temp_df = station_df.copy()
temp_df["Date"]=pd.to_datetime(temp_df["Date"])

temp_df=temp_df[temp_df["Type"]=="classic"]

temp_df=temp_df.drop_duplicates(subset="Date")

temp_df=temp_df[(temp_df["Date"].dt.month>3) & (temp_df["Date"].dt.month<11)]

temp_df=temp_df[(temp_df["Type"]=="classic") | (temp_df["Type"]=="classicDP")]


temp_df=temp_df.groupby("Date").size().reset_index(name='count')
temp_df.groupby("count").count()



#For summer sb percentage table:
temp_df = station_df.copy()
temp_df["Date"]=pd.to_datetime(temp_df["Date"])
temp_df=temp_df[(temp_df["Date"].dt.month>3) & (temp_df["Date"].dt.month<11)]
temp_df=temp_df.groupby(["Station","Type"]).size()#/2354.0*100
#Normalizing to how many dates are missing a substantial amount of data
temp_df[340] = temp_df[340]/(2354.0-10)*100
temp_df[350] = temp_df[350]/(2354.0-9)*100
temp_df[360] = temp_df[360]/(2354.0-46)*100
temp_df[380] = temp_df[380]/(2354.0-32)*100
temp_df[410] = temp_df[410]/(2354.0-10)*100
temp_df[420] = temp_df[420]/(2354.0-16)*100
temp_df[440] = temp_df[440]/(2354.0-10)*100
temp_df[450] = temp_df[450]/(2354.0-13)*100
temp_df[480] = temp_df[480]/(2354.0-13)*100
temp_df[490] = temp_df[490]/(2354.0-15)*100
print(temp_df)


#For analysis of multiple detections:
temp_df = station_df.copy()
temp_df["Date"]=pd.to_datetime(temp_df["Date"])
temp_df=temp_df[temp_df["Type"]=="classic"]
temp_df=temp_df[(temp_df["Date"].dt.month>5) & (temp_df["Date"].dt.month<9)]
temp_df=temp_df.groupby("Date").size().reset_index(name='count')
temp_df.groupby("count").count()

temp_df2 = station_df.loc[480]
temp_df2["Date"]=pd.to_datetime(temp_df2["Date"])
# temp_df2=temp_df2[temp_df2["Type"]=="classic"]
temp_df2=temp_df2[(temp_df2["Date"].dt.month>5) & (temp_df2["Date"].dt.month<9)]
st=[380,490]
temp_df = station_df.loc[st]
temp_df["Date"]=pd.to_datetime(temp_df["Date"])
# temp_df=temp_df[temp_df["Type"]=="classic"]
temp_df=temp_df[(temp_df["Date"].dt.month>5) & (temp_df["Date"].dt.month<9)]
temp_df[np.isin(temp_df["Date"],temp_df2["Date"])]

#Same as above but by coast:
cst="West"


if cst=="East":
    st=[340,410,420,440]
else:
    st=[360,350,380,450,480,490]

temp_df = station_df.loc[st]
temp_df["Date"]=pd.to_datetime(temp_df["Date"])
temp_df=temp_df[temp_df["Type"]=="classic"]
temp_df=temp_df[(temp_df["Date"].dt.month>5) & (temp_df["Date"].dt.month<9)]
temp_df=temp_df.groupby("Date").size().reset_index(name='count')
temp_df.groupby("count").count()/np.sum(temp_df.groupby("count").count())*100



#Determining how many data points are missing in temperature and wind direction
for st in [340,350,360,380,410,420,440,450,480,490]:
    file = "/Volumes/LaCie/SeaBreeze/Florida/UFLFormattedData/{0}_formatted.csv".format(st)
    data_df = pd.read_csv(file,header=0,index_col=None)
    data_df["local_eastern_time"] = pd.to_datetime(data_df["local_eastern_time"])
    temp_df = data_df[(data_df["local_eastern_time"].dt.month>3) & (data_df["local_eastern_time"].dt.month<11)]
    temp_df = temp_df.groupby(data_df["local_eastern_time"].dt.date).count()
    print(st)
    print(temp_df[(temp_df["wind_direction_10m_deg"]<72) | (temp_df["temp_air_2m_C"]<72)]["local_eastern_time"].count())
    temp_df = data_df[(data_df["local_eastern_time"].dt.month>5) & (data_df["local_eastern_time"].dt.month<9)]
    temp_df = temp_df.groupby(data_df["local_eastern_time"].dt.date).count()
    print(temp_df[(temp_df["wind_direction_10m_deg"]<72) | (temp_df["temp_air_2m_C"]<72)]["local_eastern_time"].count())


