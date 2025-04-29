"""
Coded by: Dan Moore

The purpose of this module is to convert original station data into a usable format for my purposes. This means establishing a timestamp column as the first for reference. Also fills in missing data with 'NaN'.

Input file columns and unit structure:
   Date      Time     airT wspd wdir pres rad   GST prec relH
yyyy-mm-dd  hr:mn:s   degC m/s  deg  mbar W.m-2 m/s mm    %
"""

import numpy as np
import math as m
import pandas as pd
import SBFilters
import os.path
import datetime as d

#Set years of interest
Year_list=[2013,2014,2015,2016,2017]

Station_list1=[110,120,121,130,140,150,160,170,180,230,240,241,250,\
260,270,275,280,290,302,303,304,311,320] #just '-1' csv's

Station_list2=[330,340,350,360,371,380,390,405,\
410,420,425,435,440,450,455,460,470,480,490] #just '-2' csv's

#Expected number of data points for all five years
maxdata=(365*5+1)*24*4

numsecs=(365*5+1)*24*60*60
#number of seconds in the 5 year range

base = d.datetime(2013, 1, 1)#start date
Time = np.array([base + d.timedelta(seconds=i) \
                  for i in range(0,numsecs,900)])
#implicit loop to skip 900 seconds or 15 minutes

for year in Year_list:

    for station in Station_list2: #Loop through years of interest
        path_name="/Volumes/LaCie/SeaBreeze/FloridaStationDataRaw/"
        suffix=".csv"
        infile=open(os.path.join\
        (path_name,str(year)+'-2'+suffix),"r")

        df=pd.read_csv(infile,sep=',',header=0,index_col=0)
                        #-creates pandas dataframe from infile
                        #-creates header from first row of data
                        #-I don't want an index column because I will
                        #make one

        df=df.loc[station]


        #df["Time_stamp"] = df["Date"].map(str) + " " + \
        #df[" Time"].map(str)    #-creates a Time_stamp column from
                                #Date and Time columns

        #df["Time_stamp"]=df["Time_stamp"].str[:-10]
                        #Cuts off "(UTC)" at end of string

        df["local_eastern_time"] = pd.to_datetime(\
                        df["local_eastern_time"])
                        #Converts "Time_stamp" to type datetime

        #df=df.drop(columns=["StationID"])
                        #erasing Date and Time columns - no longer
                        #needed because we have Time_stamp column

        #df=df.set_index("Time_stamp")
                        #Setting Time_stamp as index

        #df=df.append(df_temp)
                        #-Appending this temp dataframe to original
                        #dataframe for every year after 2013

    #df=df.iloc[::-1]
        #Flips dataframe upside down - originally in format:
        #12-12-2017 23:55, ... ,1-1-2013 00:00
        #Now in proper format: 1-1-2013 00:00, ... ,12-12-2017 23:55

        df=df.fillna('NaN')
            #Fills blanks with NaN

        #df=df.reindex(index=Time,fill_value='NaN')
        #Reindex data frame with all inclusive Time_stamp to
        #eliminate missing dates in data - fills blank spaces with
        #'NaN'


        df.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/UFLFormattedData/"+\
            str(station)+'_'+str(year)+".csv"))
            #Prints dataframe to csv.




