"""
Coded by: Dan Moore

The purpose of this program is to create a 5-minute time stamp for given years (for our purpose the years are 2013 through 2017), and outputs a .csv file with one column of those time stamps.
"""

import numpy as np
import datetime as d

#Allows csv to save unlimited amount of rows
np.set_printoptions(threshold=np.inf)

numsecs=(365*5+1)*24*60*60
#number of seconds in the 5 year range

base = d.datetime(2013, 1, 1)#start date
Time = np.array([base + d.timedelta(seconds=i) \
                  for i in xrange(0,numsecs,300)])
#implicit loop to skip 300 seconds or 5 minutes

np.savetxt("/Volumes/LaCie/SeaBreeze/TimeStamp.csv",\
            Time,fmt='%s',delimiter=",")



