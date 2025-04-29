"""
Coded by: Daniel Moore

The purpose of this program is to ingest the final dataset files and run statistics on each event by station
"""

import numpy as np
import math as m
import pandas as pd
import os.path
import datetime as d
import csv
from copy import deepcopy

Station_list=["DRHB","DSMY","DADV","DBBB","DDFS",\
"DBNG","DGES","DGUM","DIRL","DJCR","DMIL","DSND",\
"DSEA","DLAU","DBRG"]


totdays=79        #Number of total cases from all stations

index=0           #Use to keep track of where items go in stat arrays

#Stat Dataframe
AAGStats=pd.DataFrame({'Time_stamp':np.empty(totdays,dtype=str),
                      'Station':np.empty(totdays,dtype=str),
                      'WindDirBef':np.zeros((totdays)),
                      'WindDirAft':np.zeros((totdays)),
                      'WindChange':np.zeros((totdays)),
                      'WindSpeedBef':np.zeros((totdays)),
                      'WindSpeedAft':np.zeros((totdays)),
                      'Precip6hr':np.zeros((totdays)),
                      'Precip12hr':np.zeros((totdays)),
                      'Precip24hr':np.zeros((totdays)),
                      'TempBef':np.zeros((totdays)),
                      'TempAft':np.zeros((totdays)),
                      'TempChange':np.zeros((totdays)),
                      'RelHumBef':np.zeros((totdays)),
                      'RelHumAft':np.zeros((totdays)),
                      'PrecipTiming':np.zeros((totdays))})

"""
WindChange over half hour
WindSpeedBef averaged over hour prior SB passage
WindSpeedAft Averaged over hour after SB passage
PrecipXhr Summed for X hours after
TempChange Over half hour
"""


for station in Station_list: #Loop through years of interest
    path_name="/Volumes/LaCie/SeaBreeze/FinalDates/"
    suffix="Dates.csv"
    infile=open(os.path.join\
    (path_name,station+suffix),"r")

    path_name="/Volumes/LaCie/SeaBreeze/FormattedData/"
    suffix=".csv"
    infileDEOS=open(os.path.join\
    (path_name,station+suffix),"r")

    obs_list=infileDEOS.readlines()
    n=len(obs_list)

    time=np.empty(n,dtype=d.datetime);temp=np.zeros((n))
    precip=np.zeros((n));ws=np.zeros((n))
    wdir=np.zeros((n));srad=np.zeros((n))
    relhum=np.zeros((n))



    i=0

    for obs in obs_list[1:]:

        time[i]=obs.split(",")[0]
        precip[i]=float(obs.split(",")[7])
        temp[i]=float(obs.split(",")[1])
        ws[i]=float(obs.split(",")[2])
        wdir[i]=float(obs.split(",")[3])
        srad[i]=float(obs.split(",")[5])
        relhum[i]=float(obs.split(",")[8])

        i+=1



    dfdates=pd.read_csv(infile,header=0,index_col=0)

    dates=dfdates['0'].tolist()

    numSB=len(dates)

    AAGStats.loc[index:numSB+index-1,['Station']]=station


    j=0

    for i in range(n):
        if time[i]==dates[j]:
            AAGStats.loc[index,['Time_stamp']]=dates[j]
            for k in range(i,i+12):
                if precip[k]>0:
                    AAGStats.loc[index,['PrecipTiming']]=(k-i)*5.
                    break

            #Run stats
            AAGStats.loc[index,['WindChange']]=abs(wdir[i-6]-wdir[i])
            if abs(wdir[i-6]-wdir[i])>180:
                AAGStats.loc[index,['WindChange']]=abs(wdir[i-6]-wdir[i])-180

            AAGStats.loc[index,['WindDirBef']]=np.mean(wdir[i-24:i])
            AAGStats.loc[index,['WindDirAft']]=np.mean(wdir[i:i+12])

            AAGStats.loc[index,['WindSpeedBef']]=np.mean(ws[i-12:i])
            AAGStats.loc[index,['WindSpeedAft']]=np.mean(ws[i:i+12])
            AAGStats.loc[index,['TempBef']]=np.mean(temp[i-12:i])
            AAGStats.loc[index,['TempAft']]=np.mean(temp[i:i+12])
            AAGStats.loc[index,['TempChange']]=temp[i-6]-temp[i]
            AAGStats.loc[index,['Precip6hr']]=np.nansum(precip[i:i+72])
            AAGStats.loc[index,['Precip12hr']]=np.nansum(precip[i:i+144])
            AAGStats.loc[index,['Precip24hr']]=np.nansum(precip[i:i+288])
            AAGStats.loc[index,['RelHumBef']]=np.mean(relhum[i-12:i])
            AAGStats.loc[index,['RelHumAft']]=np.mean(relhum[i:i+12])





            index+=1        #Put stats in correct location
            j+=1
            if j==numSB:
                break


"""
AAGStats.set_index(['Station','Time_stamp'],inplace=True, drop=True)

AAGStats.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/Statistics/"\
                          ,"AAGStats.csv"))



statindex=['Mean','MAD']

dftotalstats=pd.DataFrame({'WindDirBef':np.zeros((2)),
                      'WindDirAft':np.zeros((2)),
                      #'WindChange':np.zeros((2)),
                      'WindSpeedBef':np.zeros((2)),
                      'WindSpeedAft':np.zeros((2)),
                      #'Precip6hr':np.zeros((2)),
                      'Precip12hr':np.zeros((2)),
                      #'Precip24hr':np.zeros((2)),
                      'TempBef':np.zeros((2)),
                      'TempAft':np.zeros((2)),
                      #'TempChange':np.zeros((2)),
                      'RelHumBef':np.zeros((2)),
                      'RelHumAft':np.zeros((2)),
                      'PrecipTiming':np.zeros((2))},
                      index=statindex)

dftotalstats.loc['Mean',['WindDirBef']]=AAGStats['WindDirBef'].mean()
dftotalstats.loc['MAD',['WindDirBef']]=AAGStats['WindDirBef'].mad()
dftotalstats.loc['Mean',['WindDirAft']]=AAGStats['WindDirAft'].mean()
dftotalstats.loc['MAD',['WindDirAft']]=AAGStats['WindDirAft'].mad()
dftotalstats.loc['Mean',['WindSpeedBef']]=AAGStats['WindSpeedBef'].mean()
dftotalstats.loc['MAD',['WindSpeedBef']]=AAGStats['WindSpeedBef'].mad()
dftotalstats.loc['Mean',['WindSpeedAft']]=AAGStats['WindSpeedAft'].mean()
dftotalstats.loc['MAD',['WindSpeedAft']]=AAGStats['WindSpeedAft'].mad()
dftotalstats.loc['Mean',['Precip12hr']]=AAGStats['Precip12hr'].mean()
dftotalstats.loc['MAD',['Precip12hr']]=AAGStats['Precip12hr'].mad()
dftotalstats.loc['Mean',['TempBef']]=AAGStats['TempBef'].mean()
dftotalstats.loc['MAD',['TempBef']]=AAGStats['TempBef'].mad()
dftotalstats.loc['Mean',['TempAft']]=AAGStats['TempAft'].mean()
dftotalstats.loc['MAD',['TempAft']]=AAGStats['TempAft'].mad()
dftotalstats.loc['Mean',['RelHumBef']]=AAGStats['RelHumBef'].mean()
dftotalstats.loc['MAD',['RelHumBef']]=AAGStats['RelHumBef'].mad()
dftotalstats.loc['Mean',['RelHumAft']]=AAGStats['RelHumAft'].mean()
dftotalstats.loc['MAD',['RelHumAft']]=AAGStats['RelHumAft'].mad()
dftotalstats.loc['Mean',['PrecipTiming']]=AAGStats['PrecipTiming'].mean()
dftotalstats.loc['MAD',['PrecipTiming']]=AAGStats['PrecipTiming'].mad()

dftotalstats.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/Statistics/"\
                          ,"OverallStats.csv"))




"""


