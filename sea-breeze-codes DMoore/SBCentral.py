"""
Coded by: Daniel Moore

The purpose of this program is to ingest a file and spit out Sea Breeze instances detected using the DetAlg module with either the 5min data.

This is an upgrade to the previous version in that it utilizes only a single pandas dataframe rather than large numpy array. It also creates a list of instancs of SB, not just if a day had one. It will be used to more efficiently evaluate the characteristics of what allows for sea breeze precipitation.
"""

import numpy as np
import math as m
import pandas as pd
import os.path
import datetime as d
import csv

Station_list=["DRHB","DBLK","DSMY","DADV","DBBB","DDFS","DELN",\
"DHAR","DBNG","DGES","DGUM","DIRL","DJCR","DMIL","DSJR","DSND",\
"DWAR","DSEA","DLAU","DBRG"]

#Expected number of data points for all five years
maxdata=(365*5+1)*24*12



for station in Station_list: #Loop through years of interest
    path_name="/Volumes/LaCie/SeaBreeze/FormattedData/"
    suffix=".csv"
    infileDEOS=open(os.path.join\
    (path_name,station+suffix),"r")

    #Read data
    obs_list=infileDEOS.readlines()
    n=len(obs_list)


    """
    ***********************************************
    BEGIN FILTER1
    ***********************************************
    """

    #Initialize values
    k=0;i=0;Winddir=np.zeros(n)
    Temp=np.zeros(n);Filt12=np.empty(n,dtype=d.datetime)
    Time=np.empty(n,dtype=d.datetime)

    #Loop reading the data
    for obs in obs_list[1:]:

        Temp[i]=float(obs.split(",")[1])

        Winddir[i]=float(obs.split(",")[3])

        Time[i]=obs.split(",")[0]

        i+=1

    i=0

    while i<len(Temp):
        if i>5: #Will only continue if we can refer back a half hour
                #in data

            nexthr=np.zeros(6)
            nexthr=Winddir[i:i+6]

            if Temp[i]!="NaN" and Winddir[i]!="NaN"\
            and Temp[i-6]!="NaN" and Winddir[i-6]!="NaN":

                #Filter 1: See module header for details
                if 180<=Winddir[i-6]<=360 and \
                    all(nex<=180 for nex in nexthr) and \
                    all(nex>=0 for nex in nexthr):

                    #Filter 2: See module header for details
                    if Temp[i]>5 and (Temp[i-6]-Temp[i])>=1:
                        Filt12[k]=Time[i]
                        k+=1

        i+=1

    #Return to top of file
    infileDEOS.seek(0)

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

    #Initialize values...again
    k=0;i=0;j=0;Precip=np.zeros(n)
    Filt3=np.empty_like(Filt12)
    Time=np.empty(n,dtype=d.datetime)

    for obs in obs_list[1:]:

        Time[i]=obs.split(",")[0]
        Precip[i]=obs.split(",")[7]
        i+=1

    i=0

    while k<len(Filt12):

        if Filt12[k]==Time[i]:

            #test for precip on either side

            if sum(Precip[i-36:i])=0 and sum(Precip[i:i+11])>1:
                Filt3[j]=Filt12[k]

                j+=1

            k+=1

        i+=1

    numSB=np.count_nonzero(Filt3)
        #Number of nonzero rows

    #diagnostic check
    print(numSB)

    Filt3=np.resize(Filt3,numSB)
        #Trims trailing blank rows
    """
    #Sends resulting timestamps to new file for future analysis
    with open(os.path.join("/Volumes/LaCie/SeaBreeze/Filter3/",\
        station+".csv"),"w") as csvFile:
        Fileout = csv.writer(csvFile,delimiter=",")
        Fileout.writerow(Filt3)
    """

    """
    ***********************************************
    END FILTER3
    ***********************************************
    """

    """
    ***********************************************
    Create spreadsheet of data on either side of
    time stamps found by the previous filters
    ***********************************************
    """

    #Make new array with new set of time stamps which are
    #the time stamps on either side of the ones found by Filt3
    #for analysis purposes

    TempTS=np.empty(n,dtype=d.datetime)
        #Temp array to house datetime objects before
        #converting them to strings for use against index
        #of dataframe.
    numdata=0 #Will re-evaluate this after every loop to determine
              #which elements of array to fill

    numsecs6hr=6*60*60 #Number of seconds in 6 hours which is our
                       #chosen time interval for data
                       #anslysis, this is used for loop purposes

    for time in Filt3:
        time=d.datetime.strptime(time, '%Y-%m-%d %H:%M:%S')
            #Convert string from Filt3 to object=datetime so that
            #manipulation via timedelta can be used
        if time not in TempTS:
            #In order to prevent duplicates in the final array
            start=time-d.timedelta(hours=3)
                #Begin list of time-stamps at three hours before
                #indicated time from Filt3
            TempTS[numdata:numdata+72]=[start+d.timedelta(seconds=i)
                for i in range(0,numsecs6hr,300)]
                    #This sets a range of times from the start time above
                    #to a time 6 hours later. This way we will have data
                    #three hours on either side of the time of sea breeze
            numdata=np.count_nonzero(TempTS)
                #Re-evaluation of how many items are in the array
                #this is vital to placing new items into the correct
                #place as seen in the array filling line.

    #print(station,numdata/72)

    TempTS=Filt3

    #TempTS=np.resize(TempTS,numdata)
        #Re-size the array to cut off empty elements

    df=pd.read_csv(infileDEOS,sep=',',header=0,index_col=0)
                    #-creates pandas dataframe from infile
                    #-creates header from first row of data

    FinTS=np.empty_like(TempTS)

    i=0
    while i<len(TempTS):
        FinTS[i]=TempTS[i].strftime('%Y-%m-%d %H:%M:%S')
            #Changing type datetime into type string for use
            #with the pandas dataframe - re-indexing.
        i+=1

    newdf=df[df.index.isin(FinTS)]

    print(newdf)
"""
    newdf.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/FinalDates/"\
                              ,station+"Dates.csv"))
"""




