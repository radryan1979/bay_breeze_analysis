import numpy as np
import pandas as pd
import sys
import math

path_names =    ["/Volumes/LaCie/SeaBreeze/RawData/TampaSoundingData/TampaSounding08_14.csv",
                "/Volumes/LaCie/SeaBreeze/RawData/TampaSoundingData/TampaSounding15_18.csv"]

regime_df = pd.DataFrame()

i = 0

mean_dir = np.empty((4)); mean_spd = np.empty((4))
mean_dir[:] = np.nan; mean_spd[:] = np.nan

error = 0
skip = True

knots_per_ms = 1.94384

for path in path_names:
    snd_df = pd.read_csv(path,header=None,index_col=None)
    
    for _,row in snd_df.iterrows():
        if row[0][0]=="#":
            if row[4]=="12":
                skip=False
                year = row[1]
                month = row[2]
                day = row[3]
            else:
                skip=True
        if skip:
            continue
        elif row[2] == '100000':
            if len(row[4])>5:
                if row[4][4]=="-":
                    mean_dir[0] = float(row[6])
                    mean_spd[0] = float(row[7])*knots_per_ms/10.0
            else:
                mean_dir[0] = float(row[7])
                mean_spd[0] = float(row[8])*knots_per_ms/10.0
                
            if mean_spd[0]==-9999. or mean_spd[0]==-8888. or \
                mean_dir[0]==-9999. or mean_dir[0]==-8888.:
                    mean_spd[0] = np.nan; mean_dir[0] = np.nan

        elif row[2] == '92500':
            if len(row[4])>5:
                if row[4][4]=="-":
                    mean_dir[1] = float(row[6])
                    mean_spd[1] = float(row[7])*knots_per_ms/10.0
            else:
                mean_dir[1] = float(row[7])
                mean_spd[1] = float(row[8])*knots_per_ms/10.0
                
            if mean_spd[1]==-9999. or mean_spd[1]==-8888. or \
                mean_dir[1]==-9999. or mean_dir[1]==-8888.:
                    mean_spd[1] = np.nan; mean_dir[1] = np.nan
        elif row[2] == '85000':
            if len(row[4])>5:
                if row[4][4]=="-":
                    mean_dir[2] = float(row[6])
                    mean_spd[2] = float(row[7])*knots_per_ms/10.0
            else:
                mean_dir[2] = float(row[7])
                mean_spd[2] = float(row[8])*knots_per_ms/10.0
            
            if mean_spd[2]==-9999. or mean_spd[2]==-8888. or \
                mean_dir[2]==-9999. or mean_dir[2]==-8888.:
                    mean_spd[2] = np.nan; mean_dir[2] = np.nan
        elif row[2] == '70000':
            if len(row[4])>5:
                if row[4][4]=="-":
                    mean_dir[3] = float(row[6])
                    mean_spd[3] = float(row[7])*knots_per_ms/10.0
            else:
                mean_dir[3] = float(row[7])
                mean_spd[3] = float(row[8])*knots_per_ms/10.0
            
            if mean_spd[3]==-9999. or mean_spd[3]==-8888. or \
                mean_dir[3]==-9999. or mean_dir[3]==-8888.:
                    mean_spd[3] = np.nan; mean_dir[3] = np.nan
                
                
            mean_spd_total = np.nanmean(mean_spd)
            mean_dir_total = np.nanmean(mean_dir)
            
            if np.isnan(mean_spd_total) or np.isnan(mean_dir_total):
                error+=1
                regime="NaN"
                regime_df = regime_df.append({   "Year": year,
                                    "Month": month,
                                    "Day": day,
                                    "Regime": regime,
                                    "Mean_SPD": mean_spd_total,
                                    "Mean_DIR": mean_dir_total}, ignore_index=True)
                
                mean_dir = np.empty((4)); mean_spd = np.empty((4))
                mean_dir[:] = np.nan; mean_spd[:] = np.nan
                
                mean_spd_total = 0.0
                mean_dir_total = 0.0
                continue
            
            if mean_spd_total<=4.0:
                regime = "1"
            elif mean_dir_total>=0.0 and mean_dir_total<120.0:
                if mean_spd_total<=10.0:
                    regime = "2"
                elif mean_spd_total>10.0:
                    regime = "3"
            elif mean_dir_total>=120.0 and mean_dir_total<190.0:
                if mean_spd_total<=10.0:
                    regime = "4"
                elif mean_spd_total>10.0:
                    regime = "5"
            elif mean_dir_total>=190.0 and mean_dir_total<290.0:
                if mean_spd_total<=10.0:
                    regime = "6"
                elif mean_spd_total>10.0:
                    regime = "7"
            elif mean_dir_total>=290.0 and mean_dir_total<=360.0:
                if mean_spd_total<=10.0:
                    regime = "8"
                elif mean_spd_total>10.0:
                    regime = "9"
            
            
            if mean_spd_total<0 or mean_dir_total<0:
                regime="NaN"
            
            regime_df = regime_df.append({   "Year": year,
                                    "Month": month,
                                    "Day": day,
                                    "Regime": regime,
                                    "Mean_SPD": mean_spd,
                                    "Mean_DIR": mean_dir}, ignore_index=True)
            
            mean_dir = np.empty((4)); mean_spd = np.empty((4))
            mean_dir[:] = np.nan; mean_spd[:] = np.nan
            
            mean_spd_total = 0.0
            mean_dir_total = 0.0
            
            i+=1

print("error",error)
        
        
