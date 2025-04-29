"""
Coded by: Dan Moore

This program will ingest a list of dates from the Hughes
detection algorithm in Florida, and a list of convection
detected by the 'FLPrecipDetect.py' algorithm.

We will run analyses on the output files to statistically
analyze the timing of SBF passage, in conjunction with
convection intensity, coverage, and location.

The output will be several analyses of time-series post-
SBF passage, including rainfall coverage and intensity.

Updated: 12-13-18
"""

# def count_pix(conv_count, temp_cum):
#     if type(conv_count) is list:
#         conv_count = temp_cum
#         return conv_count
#     else:
#         conv_count

import pandas as pd
import numpy as np
from datetime import datetime,timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
import os


times_dir = '/Volumes/LaCie/SeaBreeze/Florida/Results9_20_18/'
coasts = ['East','West']

conv_dir = ''


for coast in coasts:
    times_infile = times_dir+'Hughes{0}.csv'.format(coast)
    station_df = pd.read_csv(times_infile,header=0,index_col=0)
    stations = pd.unique(station_df.index)
    
    for st in stations:
        temp_df = station_df.loc[st]
        conv_path = '/Users/dpm/Desktop/Hyperion/Ref_Bins{0}/'.format(st)
        if os.path.isdir(conv_path):
            files = os.listdir(conv_path)
            conv_df = pd.DataFrame([])
            for file in files:
                if file=='.DS_Store':
                    continue
                conv_df = conv_df.append(pd.read_csv(conv_path+file,header=0,index_col=12,parse_dates=True))
            if st==350:
                conv_df = conv_df.drop(columns=['Unnamed: 0'])
            else:
                conv_df = conv_df.drop(columns=['Unnamed: 0','Time(UTC)'])
            
            precip_cum = np.zeros((len(temp_df),6))
            ref_cum = np.zeros((len(temp_df),6))
            j=0
            
            for _,row in temp_df.iterrows():
                SBF_hour = str(int(row['Hour']/1))
                SBF_min = str(int(row['Hour']%1*60))
                date_str = row['Date']+' '+SBF_hour+':'+SBF_min
                SBF_time = datetime.strptime(date_str,'%d-%b-%Y %H:%M')
                
                temp_cum = np.zeros((6,9))
                temp_int = np.zeros((6,9))
                
                for i in range(6):
                    begtime = SBF_time+timedelta(hours=i)
                    endtime = SBF_time+timedelta(hours=i+1)
                    temp_conv = conv_df.loc[(conv_df.index>begtime) & (conv_df.index<endtime)]
                    temp_conv = temp_conv.drop(columns=['Date','Region'])
                    conv_sum = temp_conv.sum(axis=0)
                    refs = np.arange(42.5,85.0,5.0)
                    int_ave = np.multiply(temp_conv.sum(axis=0),refs)
                    if temp_conv.empty:
                        continue
                    temp_cum[i][:] = conv_sum
                    temp_int[i][:] = int_ave
                    
                if temp_cum.sum() > 0:
                    precip_cum[j][:] = temp_cum.sum(axis=1)
                    ref_cum[j][:] = np.divide(temp_int.sum(axis=1),temp_cum.sum(axis=1))
                    j+=1
            precip = precip_cum[~np.all(precip_cum==0, axis=1)]
            precip_mask = precip!=0
            filtered_precip = [d[m] for d, m in zip(precip.T, precip_mask.T)]
            plt.xlabel('Hr Since SBF Passage')
            plt.ylabel('Precipitation Coverage [pixels]')
            plt.title('{0} Areal Coverage of Precipitation on SB Days'.format(st))
            plt.boxplot(filtered_precip)
            plt.show()
            # where_is_nan = np.isnan(ref_cum)
            # ref_cum[where_is_nan] = 0
            ref = ref_cum[~np.all(ref_cum==0, axis=1)]
            ref_mask = ~np.isnan(ref)
            filtered_ref = [d[m] for d, m in zip(ref.T, ref_mask.T)]
            plt.xlabel('Hr Since SBF Passage')
            plt.ylabel('Reflectivity [dBZ]')
            plt.title('{0} Precipitation Intensity on SB Days'.format(st))
            plt.boxplot(filtered_ref)
            plt.show()
    




