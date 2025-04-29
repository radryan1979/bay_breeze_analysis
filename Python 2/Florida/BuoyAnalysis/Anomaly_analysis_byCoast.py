"""
Coded by: Dan Moore

This program will ingest a list of dates from the Hughes
detection algorithm in Florida, and a list of SST Indices.

We will analyze the statistics if there are more SB detected
during negative SST anomalies.

Update 6/4/19: Analyzing 2 and 3 standard deviations.

Updated: 6/4/19
"""

import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy import stats
from scipy.stats import chisquare

times_dir = '/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19/'

BUOYdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/SST_Indices/"

W_SSTpath = BUOYdirpath + '42013_SST_Index.csv'
NE_SSTpath = BUOYdirpath + '41009_SST_Index.csv'
SE_SSTpath = BUOYdirpath + 'MLRF1_SST_Index.csv'

W_SST_df = pd.read_csv(W_SSTpath,header=0,index_col=0, parse_dates=True)
NE_SST_df = pd.read_csv(NE_SSTpath,header=0,index_col=0, parse_dates=True)
SE_SST_df = pd.read_csv(SE_SSTpath,header=0,index_col=0, parse_dates=True)

times_infile = times_dir+'Hughes_SST.csv'
station_df = pd.read_csv(times_infile,header=0,index_col=[0], parse_dates=True)

E_df = station_df.loc[[340,410,420,440]]
E_df["Date"] = pd.to_datetime(E_df["Date"])
W_df = station_df.loc[[360,350,380,480,490,450]]
W_df["Date"] = pd.to_datetime(W_df["Date"])


print("Overall")



for i in range(3):
    if i==0:
        print("NE")
        temp_df = E_df.copy()
        temp_df = temp_df.drop_duplicates(subset="Date")
        temp_df["idx"]=temp_df["SST_Anomaly_NE"]
        SST_df = NE_SST_df.copy()
        # perc_neganom=len(temp_df[temp_df["SST_Anomaly_NE"]=="Negative"])/len(SST_df[SST_df["Index"]=="Negative"])*100
    elif i==1:
        print("SE")
        temp_df = E_df.copy()
        temp_df = temp_df.drop_duplicates(subset="Date")
        temp_df["idx"]=temp_df["SST_Anomaly_SE"]
        SST_df = SE_SST_df.copy()
        # perc_neganom=len(temp_df[temp_df["SST_Anomaly_SE"]=="Negative"])/len(SST_df[SST_df["Index"]=="Negative"])*100
    elif i==2:
        print("W")
        temp_df = W_df.copy()
        temp_df = temp_df.drop_duplicates(subset="Date")
        temp_df["idx"]=temp_df["SST_Anomaly_W"]
        SST_df = W_SST_df.copy()
        # perc_neganom=len(temp_df[temp_df["SST_Anomaly_W"]=="Negative"])/len(SST_df[SST_df["Index"]=="Negative"])*100
    
    
    # sb_detect_rate = len(temp_df)/len(SST_df)*100
    
    neg_exp = len(SST_df[SST_df["Index"]=="Negative"])/len(SST_df)*len(temp_df)
    pos_exp = len(SST_df[SST_df["Index"]=="Positive"])/len(SST_df)*len(temp_df)
    non_exp = len(SST_df[SST_df["Index"]=="None"])/len(SST_df)*len(temp_df)
    
    
    sb_neg_count = len(temp_df[temp_df["idx"]=="Negative"])
    sb_pos_count = len(temp_df[temp_df["idx"]=="Positive"])
    sb_non_count = len(temp_df[temp_df["idx"]=="None"])
    
    print(chisquare([sb_neg_count, sb_pos_count, sb_non_count], f_exp=[neg_exp, pos_exp, non_exp]))
    
    # print("Percent of sea breeze detection: {0}".format(sb_detect_rate))
    # print("Percent of sea breeze detection given neg anomaly: {0}".format(perc_neganom))
    

print("Summer")
for i in range(3):
    if i==0:
        print("NE")
        temp_df = E_df.copy()
        temp_df = temp_df.drop_duplicates(subset="Date")
        temp_df["idx"]=temp_df["SST_Anomaly_NE"]
        temp_df = temp_df[(temp_df["Date"].dt.month>=6) & (temp_df["Date"].dt.month<=8)]
        SST_df = NE_SST_df.copy()
        SST_df = SST_df[(SST_df.index.month>=6) & (SST_df.index.month<=8)]
        # try:
        #     perc_neganom=np.nan_to_num(len(temp_df[temp_df["SST_Anomaly_NE"]=="Negative"])/len(SST_df[SST_df["Index"]=="Negative"])*100)
        # except ZeroDivisionError:
        #     perc_neganom = 0
    elif i==1:
        print("SE")
        temp_df = E_df.copy()
        temp_df = temp_df.drop_duplicates(subset="Date")
        temp_df["idx"]=temp_df["SST_Anomaly_SE"]
        temp_df = temp_df[(temp_df["Date"].dt.month>=6) & (temp_df["Date"].dt.month<=8)]
        SST_df = SE_SST_df.copy()
        SST_df = SST_df[(SST_df.index.month>=6) & (SST_df.index.month<=8)]
        # try:
        #     perc_neganom=np.nan_to_num(len(temp_df[temp_df["SST_Anomaly_SE"]=="Negative"])/len(SST_df[SST_df["Index"]=="Negative"])*100)
        # except ZeroDivisionError:
        #     perc_neganom = 0
    elif i==2:
        print("W")
        temp_df = W_df.copy()
        temp_df = temp_df.drop_duplicates(subset="Date")
        temp_df["idx"]=temp_df["SST_Anomaly_W"]
        temp_df = temp_df[(temp_df["Date"].dt.month>=6) & (temp_df["Date"].dt.month<=8)]
        SST_df = W_SST_df.copy()
        SST_df = SST_df[(SST_df.index.month>=6) & (SST_df.index.month<=8)]
        # try:
        #     perc_neganom=np.nan_to_num(len(temp_df[temp_df["SST_Anomaly_W"]=="Negative"])/len(SST_df[SST_df["Index"]=="Negative"])*100)
        # except ZeroDivisionError:
        #     perc_neganom = 0
    
    
    # sb_detect_rate = len(temp_df)/len(SST_df)*100
    
    # print("Percent of sea breeze detection: {0}".format(sb_detect_rate))
    # print("Percent of sea breeze detection given neg anomaly: {0}".format(perc_neganom))
    
    neg_exp = len(SST_df[SST_df["Index"]=="Negative"])/len(SST_df)*len(temp_df)
    pos_exp = len(SST_df[SST_df["Index"]=="Positive"])/len(SST_df)*len(temp_df)
    non_exp = len(SST_df[SST_df["Index"]=="None"])/len(SST_df)*len(temp_df)
    
    
    sb_neg_count = len(temp_df[temp_df["idx"]=="Negative"])
    sb_pos_count = len(temp_df[temp_df["idx"]=="Positive"])
    sb_non_count = len(temp_df[temp_df["idx"]=="None"])

    
    print(chisquare([sb_neg_count, sb_pos_count, sb_non_count], f_exp=[neg_exp, pos_exp, non_exp]))
    
    
    
print("Summer Classic")
for i in range(3):
    if i==0:
        cst = "NE"
        print(cst)
        temp_df = E_df.copy()
        temp_df = temp_df[temp_df["Type"]=="classic"]
        temp_df = temp_df.drop_duplicates(subset="Date")
        temp_df["idx"]=temp_df["SST_Anomaly_NE"]
        temp_df = temp_df[(temp_df["Date"].dt.month>=6) & (temp_df["Date"].dt.month<=8)]
        SST_df = NE_SST_df.copy()
        SST_df = SST_df[(SST_df.index.month>=6) & (SST_df.index.month<=8)]
        # try:
        #     perc_neganom=np.nan_to_num(len(temp_df[temp_df["SST_Anomaly_NE"]=="Negative"])/len(SST_df[SST_df["Index"]=="Negative"])*100)
        # except ZeroDivisionError:
        #     perc_neganom = 0
    elif i==1:
        cst = "SE"
        print(cst)
        temp_df = E_df.copy()
        temp_df = temp_df[temp_df["Type"]=="classic"]
        temp_df = temp_df.drop_duplicates(subset="Date")
        temp_df["idx"]=temp_df["SST_Anomaly_SE"]
        temp_df = temp_df[(temp_df["Date"].dt.month>=6) & (temp_df["Date"].dt.month<=8)]
        SST_df = SE_SST_df.copy()
        SST_df = SST_df[(SST_df.index.month>=6) & (SST_df.index.month<=8)]
        # try:
        #     perc_neganom=np.nan_to_num(len(temp_df[temp_df["SST_Anomaly_SE"]=="Negative"])/len(SST_df[SST_df["Index"]=="Negative"])*100)
        # except ZeroDivisionError:
        #     perc_neganom = 0
    elif i==2:
        cst = "W"
        print(cst)
        temp_df = W_df.copy()
        temp_df = temp_df[temp_df["Type"]=="classic"]
        temp_df = temp_df.drop_duplicates(subset="Date")
        temp_df["idx"]=temp_df["SST_Anomaly_W"]
        temp_df = temp_df[(temp_df["Date"].dt.month>=6) & (temp_df["Date"].dt.month<=8)]
        SST_df = W_SST_df.copy()
        SST_df = SST_df[(SST_df.index.month>=6) & (SST_df.index.month<=8)]
        # try:
        #     perc_neganom=np.nan_to_num(len(temp_df[temp_df["SST_Anomaly_W"]=="Negative"])/len(SST_df[SST_df["Index"]=="Negative"])*100)
        # except ZeroDivisionError:
        #     perc_neganom = 0
    
    
    # sb_detect_rate = len(temp_df)/len(SST_df)*100
    
    # print("Percent of sea breeze detection: {0}".format(sb_detect_rate))
    # print("Percent of sea breeze detection given neg anomaly: {0}".format(perc_neganom))
    
    neg_exp = len(SST_df[SST_df["Index"]=="Negative"])/len(SST_df)*len(temp_df)
    pos_exp = len(SST_df[SST_df["Index"]=="Positive"])/len(SST_df)*len(temp_df)
    non_exp = len(SST_df[SST_df["Index"]=="None"])/len(SST_df)*len(temp_df)
    
    
    sb_neg_count = len(temp_df[temp_df["idx"]=="Negative"])
    sb_pos_count = len(temp_df[temp_df["idx"]=="Positive"])
    sb_non_count = len(temp_df[temp_df["idx"]=="None"])
    
    
    # print("Expected",neg_exp,pos_exp,non_exp)
    # print("Actual",sb_neg_count,sb_pos_count,sb_non_count)
    
    print(chisquare([sb_neg_count, sb_pos_count, sb_non_count], f_exp=[neg_exp, pos_exp, non_exp]))
    

    print("sb",SST_df.loc[temp_df["Date"].values]["Normalization"].mean(),
        "+/-",SST_df.loc[temp_df["Date"].values]["Normalization"].std())
    print("Ave",SST_df["Normalization"].mean(),"+/-",SST_df["Normalization"].std())
    
    x = SST_df.loc[temp_df["Date"].values]
    x = x[np.isfinite(x["Normalization"])]["Normalization"].values
    bins = np.linspace(np.floor(np.nanmin(x)),np.ceil(np.nanmax(x)), 30)
    plt.hist(x,bins,alpha=1)
    plt.xlabel('SST Deviation from Mean [\N{DEGREE SIGN}C]')
    plt.ylabel('Frequency')
    plt.title('Florida Buoy {0}'.format(cst))
    if i==0:
        plt.text(-5, 60, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nstd=%(std)2.2f\N{DEGREE SIGN}C' %\
                 {'mu':x.mean(),\
                 'std':x.std()})
    elif i==1:
        plt.text(1.5, 60, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nstd=%(std)2.2f\N{DEGREE SIGN}C' %\
                 {'mu':x.mean(),\
                 'std':x.std()})
    else: 
        plt.text(1.5, 35, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nstd=%(std)2.2f\N{DEGREE SIGN}C' %\
                 {'mu':x.mean(),\
                 'std':x.std()})
    plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/SST_Indices/SBDays_DetrendedHisto_Wstd_{0}.png".format(cst),dpi=500)
    plt.show()




