"""
Daniel Moore

This module houses all filters for sea breeze detection for
my thesis. The filters are as follows:
1 - Current wind direction has easterly component and wind
    direction 30 mins prior has westerly component
2 - Air temperature has dropped 1 degC over same 30 min interval
3 - Precipitation occurs at at least one station during the day
4 -
5 -
6 -

"""

import numpy as np
import math as m
import csv


#Needs input file sent to module in format:
#infile = open ("/Users/dpmoore2927/Desktop/Test.csv", "r")
def filt12(infile,station,numdays):

    #Read data
    obs_list=infile.readlines()
    n=len(obs_list)

    #Initialize values
    k=0;i=0;rawdate="";date=0.0;
    Date=np.zeros(numdays);Winddir=np.zeros(n)
    Temp=np.zeros(n);Filt12=np.zeros(numdays)

    #Loop reading the data
    for obs in obs_list[1:]:
        rawdate=obs.split(",")[0]

        date=float(rawdate.replace("-","")) #Puts date into format: yyyymmdd

        #Check for/Set date
        if k==0:#If first pass, sets date
            Date[k]=date
            k+=1
        elif date!=Date[k-1]:#After first pass, checks to see if current date
            Date[k]=date     #is different from previous
            k+=1

        #Check for missing data. If not missing, fill arrays
        if obs.split(",")[2]=="NAN" or obs.split(",")[2]=="":
            Temp[i]=999
        else:
            Temp[i]=float(obs.split(",")[2])

        if obs.split(",")[4]=="" or obs.split(",")[4]=="NAN":
            Winddir[i]=999
        else:
            Winddir[i]=float(obs.split(",")[4])

        if i>5: #Will only continue if we can refer back a half hour
                #in data

            if Temp[i]!=999 and Winddir[i]!=999\
            and Temp[i-6]!=999 and Winddir[i-6]!=999:

                if Filt12[k-1]==1: #Skips calculations if date has already
                    i+=1       #been flagged as a SB day
                    continue

                #Filter 1: See module header for details
                if 180<=Winddir[i-6]<=360 and 0<=Winddir[i]<180:

                    #Filter 2: See module header for details
                    if (Temp[i-6]-Temp[i])<=-1:
                        Filt12[k-1]=1

        i+=1


#Print
    #Return to top of file
    infile.seek(0)

    #print(station,sum(Filt12)) #purely check to see if detecting any SBD
    return Filt12





def filt3(infile,Filt3,numdays):

    """
    takes in column of 365 zeros and makes it 1 if that day
    has precipitation, else if the date already has a 1, then
    the day is skipped - this way we don't overwrite a 1 with
    a 0.
    """

    #Read data
    obs_list=infile.readlines()

    #Initialize values
    """
    -k is counting variable to fill daily arrays
    -rawdate is used to take the date string and convert it to a readable date below.
    -minprecip is the lowest amount of precipitation we accept as a notification that there was precipitation on a given day (in mm)
    -date is readable version of rawdate
    -Date is array that houses the dates. we use this to compare to the date of the next row of data to see if we have changed days.
    """
    k=0;rawdate="";date=0.0;minprecip=0.0
    Date=np.zeros(numdays);

    #Loop reading the data
    for obs in obs_list[1:]:
        rawdate=obs.split(",")[0]
        date=float(rawdate.replace("-","")) #Puts date into format: yyyymmdd

        if k==0:#If first pass, sets date
            Date[k]=date
            k+=1
        elif date!=Date[k-1]:#After first pass, checks to see if current date
            Date[k]=date     #is different from previous
            k+=1

        #Check to see if date has already been flagged as SB.
        if Filt3[k-1]==1.0:
            continue

        #Filter 3: See module header for details.
        if obs.split(",") [8]=="": #Check for missing data
            Filt3[k-1]=0.0
        elif float(obs.split(",")[8])>minprecip:#See initialization for info.
            Filt3[k-1]=1.0
        else:
            Filt3[k-1]=0.0

    infile.seek(0)#Return to top of the file
    return Filt3




