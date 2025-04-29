import os
import pandas as pd
import numpy as np
import datetime as dt

# state='Delaware'
state='Florida'

directory = os.fsencode("/Volumes/LaCie/SeaBreeze/{0}/RadarDetection/Results_wo_conservremap/".format(state))

success=0
cnt=0
june_cnt=0
july_cnt=0
aug_cnt=0

first=True

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    
    if filename.endswith("ScanLog.csv"): 
        temp_df = pd.read_csv(os.path.join(directory,file).decode("utf-8"),index_col=0,header=0)
        if first:
            df=temp_df.copy()
            first=False
        else:
            df=pd.concat([df,temp_df],ignore_index=True)
    
    if filename.endswith("RadarDetect.csv"): 
        if filename[5] == "6":
            june_cnt+=1
        elif filename[5] == "7":
            july_cnt+=1
        elif filename[5] == "8":
            aug_cnt+=1
        success+=1
    elif filename.endswith("DayLog.csv"): 
        cnt+=1
    
df["Scan_Time"] = pd.to_datetime(df["Scan_Time"])
    
print("June: {0}".format(june_cnt))
print("July: {0}".format(july_cnt))
print("August: {0}".format(aug_cnt))