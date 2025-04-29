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

#Expected number of data points for all five years
maxdata=(365*5+1)*24*4



for station in Station_list: #Loop through years of interest
    path_name="/Volumes/LaCie/SeaBreeze/UFLFormattedData/"
    suffix="_formatted.csv"
    infile=open(os.path.join\
    (path_name,str(station)+suffix),"r")

    #Read data
    df=pd.read_csv(infile,sep=',',header=0,index_col=False)
    df=df.set_index(pd.DatetimeIndex(df['local_eastern_time']))


    """
    ***********************************************
    BEGIN FILTER1
    ***********************************************
    """

    #Initialize values
    k=0;i=0;Winddir=np.zeros(maxdata)
    Temp=np.zeros(maxdata);Filt12=np.zeros(maxdata,dtype='datetime64[s]')
    Time=np.zeros(maxdata,dtype='datetime64[s]')


    Temp[:]=df['temp_air_2m_C']

    Winddir[:]=df['wind_direction_10m_deg']

    Time[:]=df.index

    for i in range(2,maxdata):

        nexthr=np.zeros(6)
        nexthr=Winddir[i:i+2]

        if Temp[i]!="NaN" and Winddir[i]!="NaN"\
        and Temp[i-2]!="NaN" and Winddir[i-2]!="NaN":

            #Filter 1: See module header for details
            if 180<=Winddir[i-2]<=360 and \
                all(nex<=180 for nex in nexthr) and \
                all(nex>=0 for nex in nexthr):

                #Filter 2: See module header for details
                if Temp[i]>5 and (Temp[i-2]-Temp[i])>=1:
                    Filt12[k]=Time[i]
                    k+=1

    #Return to top of file
    infile.seek(0)

    numSB=np.count_nonzero(Filt12)
        #Number of nonzero rows

    #diagnostic check
    #print(numSB)




    Filt12=np.resize(Filt12,numSB)
        #Trims trailing blank rows

    """
    #Sends resulting timestamps to new file for future analysis
    with open(os.path.join("/Volumes/LaCie/SeaBreeze/Filter1/",\
        station+".csv"),"w") as csvFile:
        Fileout = csv.writer(csvFile,delimiter=",")
        Fileout.writerow(Filt12)
    """

    """
    ***********************************************
    END FILTER12
    ***********************************************
    """

    """
    ***********************************************
    BEGIN FILTER3
    ***********************************************
    """


    print(df.loc[pd.date_range(start=pd.to_datetime(Filt12[0])-d.timedelta(hours=3),end=pd.to_datetime(Filt12[0]),freq='15min')]['rain_2m_mm'])
    """
    #Initialize values...again
    k=0;i=0;j=0;
    Filt3=np.empty_like(Filt12)

    for i in range(len(Filt12)):

        if sum(df['rain_2m_mm'][d.datetime(Filt12[i])-d.timedelta(seconds=180)\
                    :d.datetime(Filt12[i])])==0 \
            and sum(df['rain_2m_mm'][d.datetime(Filt12[i]):\
                    d.datetime(Filt12[i])+d.timedelta\
                    (seconds=60)])>1:
            Filt3[j]=Filt12[i]

            j+=1

    numSB=np.count_nonzero(Filt3)
        #Number of nonzero rows

    #diagnostic check
    print(numSB)

    Filt3=np.resize(Filt3,numSB)
        #Trims trailing blank rows
    """
