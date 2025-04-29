import numpy as np
import math as m
import pandas as pd
import os.path
import datetime as d
import csv

"""
def timerange(time):
    starttime=time-d.timedelta(seconds=3600)
    numsecs2hr=7200
    Rg=np.array()
    for i in range(starttime,numsecs2hr,300):



numsecs=(365*2)*24*60*60
#number of seconds in the 5 year range

numsecs2hr=7200

base = d.datetime(2013, 1, 1)#start date
TS = np.array([base + d.timedelta(seconds=i) \
                  for i in range(0,numsecs,14400)])

newTS=np.empty(365*2*24*60,dtype=d.datetime)
numdata=0


for time in TS:
    start=time - d.timedelta(hours=1)
    newTS[numdata:numdata+24]=[start+d.timedelta(seconds=i)
        for i in range(0,numsecs2hr,300)]
    numdata=np.count_nonzero(newTS)

newTS=np.resize(newTS,numdata)

print(newTS)

"""

"""
nexthr=np.zeros(12)
nexthr=[1,2,3,4,5,6,7,8,9,10,11,12]

print(nexthr)
print(all(nex>0 for nex in nexthr))

"""

Station_list=["DRHB","DBLK","DSMY","DADV","DBBB","DDFS","DELN",\
"DHAR","DBNG","DGES","DGUM","DIRL","DJCR","DMIL","DSJR","DSND",\
"DWAR","DSEA","DLAU","DBRG"]

for station in Station_list: #Loop through years of interest
    path_name="/Volumes/LaCie/SeaBreeze/FormattedData/"
    suffix=".csv"
    infileDEOS=open(os.path.join\
    (path_name,station+suffix),"r")

    #Read data
    obs_list=infileDEOS.readlines()
    n=len(obs_list)

    data=

    for obs in obs_list[1:]:



