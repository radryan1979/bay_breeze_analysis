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
import DetAlg
import ReadDate
import os.path

#Allows csv to save unlimited amount of rows
np.set_printoptions(threshold=np.inf)

#Changes when more/less data is introduced
Station_list=["DRHB","DBLK","DSMY"]

#Set years of interest
Year_list=[2013,2014,2015,2016,2017]


yrcnt=0 #Initialize year count
while(yrcnt<len(Year_list)): #Loop through years of interest
    year=Year_list[yrcnt]
    if year==2013 or year==2014 or year==2015 or year==2017:
        numdays=365
    else:
        numdays=366 #2016 is a leap year

    #Initialize
    """
    -Sea_breeze array houses binary (yes or no/1 or 0) whether each day indicates a sea breeze or not, for each station.
    -Filter12 is also a binary array which indicates whether each station passes filters 1 and 2 of our detection algorithm each day (see DetAlg module)
    -Filter3 is binary array which indicates whether each day passes filter 3 of our detection algorithm (see DetAlg module). If there is precipitation at any station, the day is marked a success.
    -Date is a 1D array that houses the dates for each year. This will act as the index for the final csv. See below.
    """
    Sea_breeze=np.zeros((numdays,len(Station_list)))
    Filter12=np.zeros((numdays,len(Station_list)))
    Filter3=np.zeros(numdays)
    Date=np.zeros(numdays)

    stcnt=0 #Initialize station count every time year changes
    while(stcnt<len(Station_list)): #Loop thru stations
        error1=0 #Notifies filters of missing dates
        station=Station_list[stcnt]

        path_name="/Volumes/LaCie/SeaBreeze/OriginalData/StationData/"
        suffix=".dat"
        infile=open(os.path.join\
        (path_name,str(year),station+str(year)+suffix),"r")

        if stcnt==0: #First pass, fill Date array
            Date=ReadDate.date(infile,numdays)
            if numdays!=np.count_nonzero(Date):
                print("Error: dates missing from",station,"in",year)
                error1=1

        #Send data to filt12 to test for sea breeze based on
        #winddir and temp drop over 30 min intervals - sea sub
        Filter12[:,stcnt]=DetAlg.filt12(infile,station,numdays)

        #Tests for precipitation at any station for each day
        #replaces each 0 with 1 if there is a precip at each
        #individual station.
        Filter3=DetAlg.filt3(infile,Filter3,numdays)

        stcnt+=1

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
    yrcnt+=1



