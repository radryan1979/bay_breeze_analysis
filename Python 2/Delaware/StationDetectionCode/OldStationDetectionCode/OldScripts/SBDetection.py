"""
Daniel Moore

The purpose of this program is to ingest a file and spit out Sea Breeze
Days detected using the DetAlg module with either the 5min data
function or the 1hr data function

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

#Changes when more/less data is introduced
Station_list=["DADV","DBBB","DBLK","DBNG","DBRG","DDFS","DELN","DGES",\
"DGUM","DHAR","DIRL","DJCR","DLAU","DMIL","DRHB","DSEA","DSJR","DSMY",\
"DSND","DWAR"]

    #Initialize
"""
    -Sea_breeze array houses binary (yes or no/1 or 0) whether each day indicates a sea breeze or not, for each station.
    -Filter12 is also a binary array which indicates whether each station passes filters 1 and 2 of our detection algorithm each day (see DetAlg module)
    -Filter3 is binary array which indicates whether each day passes filter 3 of our detection algorithm (see DetAlg module). If there is precipitation at any station, the day is marked a success.
"""

numdata=len(Time)
Sea_breeze=np.zeros((numdata,len(Station_list)))
Filter12=np.zeros((numdays,len(Station_list)))
Filter3=np.zeros(numdays)
Date=np.zeros(numdays)

for station in Stateion_list: #Loop thru stations

    #creating path to individual DEOS statino data by year
    path_name="/Volumes/LaCie/SeaBreeze/FormattedData/"
    suffix=".csv"
    infileDEOS=open(os.path.join\
    (path_name,str(year),station+suffix),"r")

    #Send data to filt12 to test for sea breeze based on
    #winddir and temp drop over 30 min intervals - sea sub
    Filter12[:,stcnt]=SBFilters.filt12(infileDEOS,station,\
                    numdays,Time)

    #Tests for precipitation at any station for each day
    #replaces each 0 with 1 if there is a precip at each
    #individual station.
    Filter3=SBFilters.filt3(infileDEOS,Filter3,numdays)

#A loop to multiply all filters to create final SB determination
k=0#Counting variable to do-while loop
while(k<numdays):
    Sea_breeze[k,:]=Filter12[k,:]*Filter3[k]
    k+=1

#Create dataframe for printing purposes - columns with station as
#header and date as index(first column)
df=pd.DataFrame(Sea_breeze,index=Date,columns=Station_list)

#Create final column of summation of SBB totals from all stations
df['Sum']=df.iloc[:,0:len(Station_list)].sum(axis=1)

#Create csv for each year
df.to_csv(os.path.join("/Volumes/LaCie/SeaBreeze/Results/SBB"+\
    str(year)+".csv"),sep=",")



