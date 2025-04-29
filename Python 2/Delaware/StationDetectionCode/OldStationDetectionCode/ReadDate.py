"""
Daniel Moore

The purpose of this module will be to ingest daily data from DEOS
with 5 minute temporal scale and simply print out dates in a
decimal format.

"""

import numpy as np
import math as m

"""
Needs input file sent to module in format:
infile = open ("/Users/dpmoore2927/Desktop/Test.txt", "r")
OR .csv OR .dat
"""
def date(infile,numdays):

    #Read data
    obs_list=infile.readlines()

    #Initialize values
    """
    -k is counting variable to fill daily arrays
    -rawdate is used to take the date string and convert it to a readable date below.
    -date is readable version of rawdate
    -Date is array that houses the dates. we use this to compare to the date of the next row of data to see if we have changed days.
    """

    k=0;rawdate="";date=0.0
    Date=np.zeros(numdays)

    #Loop reading the data
    for obs in obs_list[1:]:
        rawdate=obs.split(",")[0]

        date=float(rawdate.replace("-","")) #Puts date into format: yyyymmdd

        #First pass, set date.
        if k==0:
            Date[k]=date
            k+=1
        #After first pass, check to see if date is same as previous
        #pass. If not, fill Date with new date.
        elif date!=Date[k-1]:
            Date[k]=date
            k+=1


    infile.seek(0)#Return to top of file
    return Date
