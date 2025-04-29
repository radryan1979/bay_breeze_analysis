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
Year_list=[2017,2016,2015,2014,2013]

Station_list=["DRHB","DBLK","DSMY","DADV","DBBB","DDFS","DELN",\
"DHAR","DBNG","DGES","DGUM","DIRL","DJCR","DMIL","DSJR","DSND",\
"DWAR","DSEA","DLAU","DBRG"]

#Expected number of data points for all five years
maxdata=(365*5+1)*24*12

numsecs=(365*5+1)*24*60*60
#number of seconds in the 5 year range

base = d.datetime(2013, 1, 1)#start date
Time = np.array([base + d.timedelta(seconds=i) \
                  for i in range(0,numsecs,300)])
#implicit loop to skip 300 seconds or 5 minutes

for station in Station_list:

    for year in Year_list: #Loop through years of interest
        path_name="/Volumes/LaCie/SeaBreeze/OriginalData/StationData/"
        suffix=".dat"
        infileDEOS=open(os.path.join\
        (path_name,str(year),station+str(year)+suffix),"r")

        if str(year)=='2017': # have to set up final data frame
                              # with correct columns

            df=pd.read_csv(infileDEOS,sep=',',header=0,index_col=False)
                            #-creates pandas dataframe from infile
                            #-creates header from first row of data
                            #-I don't want an index column because I will
                            #make one

            df["Time_stamp"] = df["Date"].map(str) + " " + \
            df[" Time"].map(str)    #-creates a Time_stamp column from
                                    #Date and Time columns

            df["Time_stamp"]=df["Time_stamp"].str[:-10]
                            #Cuts off "(UTC)" at end of string

            df["Time_stamp"] = pd.to_datetime(df["Time_stamp"])
                            #Converts "Time_stamp" to type datetime

            df=df.drop(columns=["Date"," Time"])
                            #erasing Date and Time columns - no longer
                            #needed because we have Time_stamp column

            df=df.set_index("Time_stamp")
                            #Setting Time_stamp as index

        else:

            df_temp=pd.read_csv(infileDEOS,sep=',',header=0,index_col=False)
                            #-creates temp pandas dataframe from infile
                            #-creates header from first row of data
                            #-I don't want an index column because I will
                            #make one
                            #-This dataframe will be appended to original
                            #dataframe for printing purposes

            df_temp["Time_stamp"] = df_temp["Date"].map(str) + " " + \
            df_temp[" Time"].map(str) #-creates a Time_stamp column from
                                      #Date and Time columns

            df_temp["Time_stamp"]=df_temp["Time_stamp"].str[:-10]
                            #See above

            df_temp["Time_stamp"] = pd.to_datetime(df_temp["Time_stamp"])
                            #See above

            df_temp=df_temp.drop(columns=["Date"," Time"])
                            #See above

            df_temp=df_temp.set_index("Time_stamp")
                            #See above

            df=df.append(df_temp)
                            #-Appending this temp dataframe to original
                            #dataframe for every year prior to 2017

    df=df.iloc[::-1]
        #Flips dataframe upside down - originally in format:
        #12-12-2017 23:55, ... ,1-1-2013 00:00
        #Now in proper format: 1-1-2013 00:00, ... ,12-12-2017 23:55

    print(station,df.isnull().sum())
    df=df.fillna('NaN')
        #Fills blanks with NaN

    df=df.reindex(index=Time,fill_value='NaN')
        #Reindex data frame with all inclusive Time_stamp to
        #eliminate missing dates in data - fills blank spaces with
        #'NaN'

    """
    df.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/FormattedData/"+\
        str(station)+".csv"))
        #Prints dataframe to csv.
    """



