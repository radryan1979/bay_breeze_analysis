"""
Coded by: Daniel Moore

The purpose of this program is to ingest a user input date, and produce a graph following several variables (probably temperature, wind direction and precipitation before and after sea breeze passage)
"""

import numpy as np
import math as m
import pandas as pd
import os.path
import datetime as d
import csv
from matplotlib import pyplot as plt
from math import radians

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

Station_list=["DRHB","DSMY","DADV","DBBB","DDFS",\
"DBNG","DGES","DGUM","DIRL","DJCR","DMIL","DSND",\
"DSEA","DLAU","DBRG"]

numdata=36        #Number of total cases from all stations

csmonth='07'#input("Enter month in format MM: ")
csday='14'#input("Enter day in format DD: ")
csyear='2017'#input("Enter year in format YYYY: ")



cstime_stamp=csyear+'-'+csmonth+'-'+csday


#Graph Dataframe
FigureData=pd.DataFrame({'Time':np.zeros((numdata)),
                      'WindDir':np.zeros((numdata)),
                      'WindSpeed':np.zeros((numdata)),
                      'Precip':np.zeros((numdata)),
                      'Temp':np.zeros((numdata)),
                      'RelHum':np.zeros((numdata))})



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



    dfdates=pd.read_csv(infile,header=0,index_col=0)

    dates=dfdates['0'].tolist()

    numSB=len(dates)

    for k in range(numSB):
        if dates[k][:10]==cstime_stamp:

            tsofficial=dates[k]

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

            for i in range(n):
                if time[i]==tsofficial:
                    FigureData['WindSpeed']=ws[i-12:i+24]
                    FigureData['WindDir']=wdir[i-12:i+24]
                    FigureData['Temp']=temp[i-12:i+24]
                    FigureData['RelHum']=relhum[i-12:i+24]
                    FigureData['Time']=np.arange(-60,120,5)
                    for k in range(36):
                        FigureData.loc[k,['Precip']]=\
                            np.sum(precip[i-12:i-12+k])
                    break


            #Print figure

            fig = plt.figure()
            host = fig.add_subplot(111)

            par1 = host.twinx()

            host.set_xlim(-60, 120)
            host.set_ylim(0, 40)
            par1.set_ylim(0, 40)

            host.set_xlabel("Time Since Sea Breeze [min]")
            par1.set_ylabel("Cumulative Precipitation [mm]")
            host.set_ylabel("Temperature [\N{DEGREE SIGN}C]")

            color1 = plt.cm.viridis(0.5)
            color2 = plt.cm.viridis(0.8)

            p2, = par1.plot(FigureData['Time'],\
                            FigureData['Precip'], \
                            color=color1,label="Precipitation",\
                            linewidth=3)
            p1, = host.plot(FigureData['Time'],\
                            FigureData['Temp'], \
                            color=color2, label="Temperature",\
                            linewidth=3)

            lns = [p1, p2]
            host.legend(handles=lns, loc='best',fontsize=16)

            host.yaxis.label.set_color(p1.get_color())
            host.yaxis.label.set_size(16)
            par1.yaxis.label.set_color(p2.get_color())
            par1.yaxis.label.set_size(16)
            host.xaxis.label.set_size(16)

            Title=cstime_stamp+' Sea Breeze at '+station

            plt.title(Title,fontsize=20)

            plt.show()



            WDBef=np.zeros((12))
            WDAft=np.zeros((24))
            WSBef=np.zeros((12))
            WSAft=np.zeros((24))

            WDBef=FigureData['WindDir'].values[0:12]
            WSBef=FigureData['WindSpeed'].values[0:12]
            WDAft=FigureData['WindDir'].values[12:36]
            WSAft=FigureData['WindSpeed'].values[12:36]


            ax = plt.subplot(111, polar=True)
            ax.scatter(x=[radians(x) for x in WDBef], y=WSBef,\
                       label='Before Sea Breeze')
            ax.scatter(x=[radians(x) for x in WDAft], y=WSAft, \
                       label='After Sea Breeze')
            ax.set_theta_zero_location('N')
            ax.set_theta_direction(-1)
            plt.title(Title,fontsize=20)

            label_position=ax.get_rlabel_position()
            ax.text(np.radians(label_position-10),ax.get_rmax()/2.,\
                    'Wind Speed [m$s^-$$^1$]',rotation=69,\
                    ha='center',va='center',fontsize=12)

            ax.set_xlabel('Wind Direction',fontsize=14)



            plt.legend(loc='best',fontsize=12)

            plt.show()









