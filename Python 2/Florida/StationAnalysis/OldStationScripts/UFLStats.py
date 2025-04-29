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

Station_list=[270,290,302,304,320,340,371,455,435,410,425,420,\
              440] #These stations are on East Coast

"""
[110,120,121,130,140,150,160,170,180,230,240,241,250,\
260,270,275,280,290,302,303,304,311,320,330,340,350,360,371,380,390,405,\
410,420,425,435,440,450,455,460,470,480,490]
"""


totdays=1040  #19#number of cases for station 110 for testing
#1040        #Number of total cases from all stations

index=0           #Use to keep track of where items go in stat arrays

#Stat Dataframe
UFLStats=pd.DataFrame({'Time_stamp':np.empty(totdays,dtype=str),
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
    path_name="/Volumes/LaCie/SeaBreeze/Florida/FinalFlorida/EastCoast/"
    suffix="times.csv"
    infile=open(os.path.join\
    (path_name,str(station)+suffix),"r")

    path_name="/Volumes/LaCie/SeaBreeze/Florida/UFLFormattedData/"
    suffix="_formatted.csv"
    infileUFL=open(os.path.join\
    (path_name,str(station)+suffix),"r")

    #Read data
    df=pd.read_csv(infileUFL,sep=',',header=0,index_col=False)
    df=df.set_index(pd.DatetimeIndex(df['local_eastern_time']))

    dfdates=pd.read_csv(infile,header=0,index_col=0)

    dates=dfdates.index.tolist()

    numSB=len(dates)

    UFLStats.loc[index:numSB+index-1,['Station']]=str(station)

    totprecip=df['rain_2m_mm'].sum()



    j=0

    for i in range(numSB):

        UFLStats.loc[index,['Time_stamp']]=dates[i]

        #Run stats
        """
        UFLStats.loc[index,['WindChange']]=abs(wdir[i-6]-wdir[i])
        if abs(wdir[i-6]-wdir[i])>180:
            UFLStats.loc[index,['WindChange']]=abs(wdir[i-6]-wdir[i])-180
        """

        UFLStats.loc[index,['WindDirBef']]=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i])-d.timedelta(hours=2),end=\
            pd.to_datetime(dates[i]),freq='15min')]\
            ['wind_direction_10m_deg'].mean()

        UFLStats.loc[index,['WindDirAft']]=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i]),end=\
            pd.to_datetime(dates[i])+d.timedelta(hours=1),freq='15min')]\
            ['wind_direction_10m_deg'].mean()

        UFLStats.loc[index,['WindSpeedBef']]=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i])-d.timedelta(hours=1),end=\
            pd.to_datetime(dates[i]),freq='15min')]\
            ['wind_speed_10m_mph'].mean()

        UFLStats.loc[index,['WindSpeedAft']]=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i]),end=\
            pd.to_datetime(dates[i])+d.timedelta(hours=1),freq='15min')]\
            ['wind_speed_10m_mph'].mean()

        UFLStats.loc[index,['TempBef']]=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i])-d.timedelta(hours=1),end=\
            pd.to_datetime(dates[i]),freq='15min')]\
            ['temp_air_2m_C'].mean()

        UFLStats.loc[index,['TempAft']]=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i]),end=\
            pd.to_datetime(dates[i])+d.timedelta(hours=1),freq='15min')]\
            ['temp_air_2m_C'].mean()

        UFLStats.loc[index,['TempChange']]=df.loc[pd.to_datetime(\
            dates[i])-d.timedelta(minutes=30)]['temp_air_2m_C']\
            -df.loc[pd.to_datetime(\
            dates[i])]['temp_air_2m_C']

        UFLStats.loc[index,['Precip6hr']]=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i]),end=\
            pd.to_datetime(dates[i])+d.timedelta(hours=6),freq='15min')]\
            ['rain_2m_mm'].sum()

        UFLStats.loc[index,['Precip12hr']]=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i]),end=\
            pd.to_datetime(dates[i])+d.timedelta(hours=12),freq='15min')]\
            ['rain_2m_mm'].sum()

        stationprecip=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i]),end=\
            pd.to_datetime(dates[i])+d.timedelta(hours=12),freq='15min')]\
            ['rain_2m_mm'].sum()

        UFLStats.loc[index,['Precip24hr']]=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i]),end=\
            pd.to_datetime(dates[i])+d.timedelta(hours=24),freq='15min')]\
            ['rain_2m_mm'].sum()

        UFLStats.loc[index,['RelHumBef']]=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i])-d.timedelta(hours=1),end=\
            pd.to_datetime(dates[i]),freq='15min')]\
            ['rh_2m_pct'].mean()

        UFLStats.loc[index,['RelHumAft']]=df.loc[pd.date_range(start=\
            pd.to_datetime(dates[i]),end=\
            pd.to_datetime(dates[i])+d.timedelta(hours=1),freq='15min')]\
            ['rh_2m_pct'].mean()

        index+=1        #Put stats in correct location

        sbprecipperc=float(stationprecip/totprecip)*100
        print(sbprecipperc)

        j+=1
        if j==numSB:
            break
"""

UFLStats.set_index(['Station','Time_stamp'],inplace=True, drop=True)

UFLStats.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/Statistics/"\
                          ,"UFLStats.csv"))


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
                      'RelHumAft':np.zeros((2))},
                      #'PrecipTiming':np.zeros((2)),
                      index=statindex)

dftotalstats.loc['Mean',['WindDirBef']]=UFLStats['WindDirBef'].mean()
dftotalstats.loc['MAD',['WindDirBef']]=UFLStats['WindDirBef'].mad()
dftotalstats.loc['Mean',['WindDirAft']]=UFLStats['WindDirAft'].mean()
dftotalstats.loc['MAD',['WindDirAft']]=UFLStats['WindDirAft'].mad()
dftotalstats.loc['Mean',['WindSpeedBef']]=UFLStats['WindSpeedBef'].mean()
dftotalstats.loc['MAD',['WindSpeedBef']]=UFLStats['WindSpeedBef'].mad()
dftotalstats.loc['Mean',['WindSpeedAft']]=UFLStats['WindSpeedAft'].mean()
dftotalstats.loc['MAD',['WindSpeedAft']]=UFLStats['WindSpeedAft'].mad()
dftotalstats.loc['Mean',['Precip12hr']]=UFLStats['Precip12hr'].mean()
dftotalstats.loc['MAD',['Precip12hr']]=UFLStats['Precip12hr'].mad()
dftotalstats.loc['Mean',['TempBef']]=UFLStats['TempBef'].mean()
dftotalstats.loc['MAD',['TempBef']]=UFLStats['TempBef'].mad()
dftotalstats.loc['Mean',['TempAft']]=UFLStats['TempAft'].mean()
dftotalstats.loc['MAD',['TempAft']]=UFLStats['TempAft'].mad()
dftotalstats.loc['Mean',['RelHumBef']]=UFLStats['RelHumBef'].mean()
dftotalstats.loc['MAD',['RelHumBef']]=UFLStats['RelHumBef'].mad()
dftotalstats.loc['Mean',['RelHumAft']]=UFLStats['RelHumAft'].mean()
dftotalstats.loc['MAD',['RelHumAft']]=UFLStats['RelHumAft'].mad()
#dftotalstats.loc['Mean',['PrecipTiming']]=UFLStats['PrecipTiming'].mean()
#dftotalstats.loc['MAD',['PrecipTiming']]=UFLStats['PrecipTiming'].mad()

dftotalstats.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/Statistics/"\
                          ,"UFLOverallStats.csv"))


"""


