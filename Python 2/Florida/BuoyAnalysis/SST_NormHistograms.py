"""
Coded by: Dan Moore

Purpose: Histograms of normalized SST to determine proper
anomoly values for each buoy

Updated: 3-7-19
"""

import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt

path_name="/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/SST_Indices/"
buoys=["MLRF1","41113","42013"]

for bu in buoys:
    df = pd.read_csv(path_name+'{0}_SST_Index.csv'.format(bu),header=0,index_col=0)
    x = df["Normalization"].dropna()
    bins = np.linspace(x.min(),x.max(),30)
    n1,_,_=plt.hist(x,bins,alpha=1)
    plt.xlabel('Normalized SST [\N{DEGREE SIGN}C]')
    plt.ylabel('Frequency')
    plt.title('Normalized SST \n Florida Buoy {0}'.format(bu))
    plt.text((bins.min())/4.0, (n1.max()-n1.min())/3.0, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nSTD=%(STD)2.2f\N{DEGREE SIGN}C' %\
             {'mu':x.mean(),'STD':x.std()})
    # plt.xlim(0)
    plt.savefig(path_name+"/Buoy{0}SST_Histogram.png".format(bu),dpi=500)
    plt.show()