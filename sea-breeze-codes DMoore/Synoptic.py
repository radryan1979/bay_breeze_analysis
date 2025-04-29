"""
Daniel Moore

The purpose of this program will be to ingest radiosonde data and send it to
SBFilters - using filt4 - to determine synoptic conditions related to sea breeze
and spit out a .csv file that will have dates and a 1 or 0 based on detection of
SB or not. This can then be used in the SBDetection program to further break down
days that see SB.

Data Format:

      1          2          3          4          5          6           7
 LINTYP
                                header lines
    254        HOUR        DAY      MONTH       YEAR    (blank)     (blank)
      1       WBAN#       WMO#        LAT D      LON D     ELEV       RTIME
      2       HYDRO       MXWD      TROPL      LINES     TINDEX      SOURCE
      3     (blank)      STAID    (blank)    (blank)      SONDE     WSUNITS

                                data lines
      9    PRESSURE     HEIGHT       TEMP      DEWPT   WIND DIR    WIND SPD
      4
      5
      6
      7
      8
                                          LEGEND

LINTYP: type of identification line
        254 = indicates a new sounding in the output file
          1 = station identification line
          2 = sounding checks line
          3 = station identifier and other indicators line
          4 = mandatory level
          5 = significant level
          6 = wind level (PPBB) (GTS or merged data)
          7 = tropopause level (GTS or merged data)
          8 = maximum wind level (GTS or merged data)
          9 = surface level

HOUR:   time of report in UTC
LAT:    latitude in degrees and hundredths
LON:    longitude in degrees and hundredths

PRESSURE: in whole hectopascal/millibars (original format)
HEIGHT:   height in meters (m)
TEMP:     temperature in tenths of degrees Celsius
DEWPT:    dew point temperature in tenths of a degree Celsius
WIND DIR: wind direction in degrees
WIND SPD: tenths of a meter per second
"""


import numpy as np
import math as m
import pandas as pd
import SBFilters
import os.path

#Allows csv to save unlimited amount of rows
np.set_printoptions(threshold=np.inf)


#Set input file - radiosonde data from NOAA
infileRS=open("/Volumes/LaCie/SeaBreeze/Soundings/KWAL01012013_12312017.txt","r")

#Set years of interest
Year_list=[2013,2014,2015,2016,2017]

yrcnt=0 #Initialize year count
while(yrcnt<len(Year_list)): #Loop through years of interest
    year=Year_list[yrcnt]
    if year==2013 or year==2014 or year==2015 or year==2017:
        numdays=365
    else:
        numdays=366 #2016 is a leap year

    #Sends file to SBFilters to get a 1D array where the first row is 010120xx and the
    #365th row is 121220xx - this needs to be inversed for use in the detection algorithm
    Filter4=SBFilters.filt4(infileRS,year,numdays)

    #Create csv for each year
    np.savetxt(os.path.join("/Volumes/LaCie/SeaBreeze/Soundings/Results/Radiosonde"+\
        str(year)+".csv"),Filter4,delimiter=',')
    yrcnt+=1



























