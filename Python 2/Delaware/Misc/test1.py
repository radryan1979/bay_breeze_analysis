"""
<<<<<<< HEAD
Coded by: Daniel Moore

The purpose of this program is to ingest a file and spit out Sea Breeze instances detected using the DetAlg module with either the 15min data.

"""

=======
Test document to test filters out one by one utilizing the same input files for the actual document.
"""


>>>>>>> 41cf07a708178b4f2a7e06bb78964b2a9df47970
import numpy as np
import math as m
import pandas as pd
import os.path
<<<<<<< HEAD
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
=======

import DetAlg


#changes when more/less data is introduced
Station_list=["DRHB","DBLK","DSMY"]

#set years of interest
Year_list=[2013,2014,2015,2016,2017]


yrcnt=0 #initialize year count
while(yrcnt<len(Year_list)): #loop through years of interest
    year=Year_list[yrcnt]
    if year==2013 or year==2014 or year==2015 or year==2017:
        numdays=365
    else:
        numdays=366 #2016 is a leap year

    #initialize
    #Sea_breeze=np.zeros((numdays,len(Station_list)))
    #Filter12=np.zeros((numdays,len(Station_list)))
    Filter3=np.zeros(numdays)
    #Date=np.zeros(numdays)

    stcnt=0 #initialize station count every time year changes
    while(stcnt<len(Station_list)): #loop thru stations
        station=Station_list[stcnt]

        path_name="/Volumes/LaCie/SeaBreeze/OriginalData/StationData/"
        suffix=".dat"
        infile=open(os.path.join\
        (path_name,str(year),station+str(year)+suffix),"r")

        #if stcnt==0: #first pass, fill Date array
            #Date=ReadDate.date(infile,numdays)

        #we send data to filt12 to test for sea breeze based on
        #winddir and temp drop over 30 min intervals - sea sub
        #Filter12[:,stcnt]=DetAlg.filt12(infile,station,numdays)

        #Tests for precipitation at any station for each day
        #replaces each 0 with 1 if there is a precip at each
        #individual station.
        Filter3=DetAlg.filt3(infile,Filter3,numdays)

        stcnt+=1
    yrcnt+=1

>>>>>>> 41cf07a708178b4f2a7e06bb78964b2a9df47970
