"""
Coded by: Daniel Moore

This module houses all filters for sea breeze detection for
my thesis. The filters are as follows:
1 - Current wind direction has easterly component and wind
    direction 30 mins prior has westerly component
2 - Air temperature has dropped 1 degC over same 30 min interval
3 - Precipitation occurs at at least one station during the day
4 -
5 - No precipitation at the station in previous 2 hours
6 -

This is an update to previous versions in that it will be used hand in hand with "SBCentral.py" utilizing a centralized pandas dataframe and taking only instances of sea breeze. Initial sea breeze detected will be from a change in wind direction. Then we will eliminate dates systematically with subsequent filters as described above.

Input file columns and unit structure:

Time_stamp',' Air Temperature (deg. C)',' Wind Speed (m.s-1)',' Wind Direction (deg.)',' Barometric Pressure (mbar)',' Solar Radiation (W.m-2)',' Wind Gust Speed (5) (m.s-1)',' Gage Precipitation (5) (mm)',' Relative Humidity (%)''

"""
import numpy as np
import math as m
import pandas as pd
import datetime as d
import os.path
import csv


def filt1_2(infile,station):

    #Read data
    obs_list=infile.readlines()
    n=len(obs_list)

    #Initialize values
    k=0;i=0;Winddir=np.zeros(n)
    Temp=np.zeros(n);Filt12=np.empty(n,dtype=d.datetime)

    #Loop reading the data
    for obs in obs_list[1:]:

        Temp[i]=float(obs.split(",")[1])

        Winddir[i]=float(obs.split(",")[3])

        if i>5: #Will only continue if we can refer back a half hour
                #in data

            if Temp[i]!="NaN" and Winddir[i]!="NaN"\
            and Temp[i-6]!="NaN" and Winddir[i-6]!="NaN":

                #Filter 1: See module header for details
                if 180<=Winddir[i-6]<=360 and 0<=Winddir[i]<180:

                    #Filter 2: See module header for details
                    if (Temp[i]-Temp[i-6])<=-1:
                        Filt12[k]=obs.split(",")[0]
                        k+=1

        i+=1


#Print
    #Return to top of file
    infile.seek(0)

    numSB=np.count_nonzero(Filt12)
        #Number of nonzero rows

    Filt12=np.resize(Filt12,numSB)
        #Trims trailing blank rows

    #Sends resulting timestamps to new file for future analysis
    with open(os.path.join("/Volumes/LaCie/SeaBreeze/Filter1/",\
        station+".csv"),"w") as csvFile:
        Fileout = csv.writer(csvFile,delimiter=",")
        Fileout.writerow(Filt12)

    return



def filt3(infile,station):

    #Read station data
    obs_list1=infile.readlines()
    n=len(obs_list1)

    #Read list of time_stamps that passed Filter 1
    obs_list2=os.path.join("/Volumes/LaCie/SeaBreeze/Filter1/",\
        station+".csv").readlines()
    numtimes=len(obs_list2)

    #Initialize values
    k=0;i=0;Precip=np.zeros(n)
    Filt3=np.empty(numtimes,dtype=d.datetime)

    #Loop reading the data
    for obs in obs_list[1:]:

        if

        Temp[i]=float(obs.split(",")[1])

        Winddir[i]=float(obs.split(",")[3])

        if i>5: #Will only continue if we can refer back a half hour
                #in data

            if Temp[i]!="NaN" and Winddir[i]!="NaN"\
            and Temp[i-6]!="NaN" and Winddir[i-6]!="NaN":

                #Filter 1: See module header for details
                if 180<=Winddir[i-6]<=360 and 0<=Winddir[i]<180:

                    #Filter 2: See module header for details
                    if (Temp[i]-Temp[i-6])<=-1:
                        Filt12[k]=obs.split(",")[0]
                        k+=1

        i+=1


