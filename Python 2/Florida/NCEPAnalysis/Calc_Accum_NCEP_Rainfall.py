"""
Programmed by: Dan Moore

The purpose of this program is to calculate accumulated
rainfall for one pixel in FL Study Area by month to correlate
with teleconections

Updated: 2-5-19
"""

import numpy as np
import pandas as pd
import os.path
from netCDF4 import Dataset as netcdf_dataset

dir_name = "/Volumes/LaCie/SeaBreeze/Florida/NCEP4kmGridSubset/"

cnt = 0

years=np.empty(120) # number of months
months=np.empty(120)
se_precip=np.empty(120); w_precip=np.empty(120)
se_precip[:]=np.nan; w_precip[:]=np.nan

days31 = [1,3,5,7,8,10,12]
days30 = [4,6,9,11]

for yr in range(2009,2019):
    for mo in range(1,13):
        print(yr,mo)
        if (yr==2012 or yr==2016) and mo==2:
            days=29
        elif mo == 2:
            days=28
        elif mo in days31:
            days=31
        elif mo in days30:
            days=30
        
        se_month_precip=0.0; w_month_precip=0.0
        
        for day in range(1,days+1):
            for hr in range(12,24):
                ncep_path = dir_name + "{0}/{0}{1}{2}{3}.nc".format(\
                yr,str(mo).zfill(2),str(day).zfill(2),str(hr).zfill(2))
                
                if os.path.isfile(ncep_path):
                    dataset=netcdf_dataset(ncep_path)
                    # print(dataset.variables["P"])
                    se_month_precip = se_month_precip+dataset.variables["P"][60,250]
                    w_month_precip = w_month_precip+dataset.variables["P"][75,190]
                    # print(ncep_path)
                else:
                    continue
                
        years[cnt] = yr
        months[cnt] = mo
        if mo<4 or mo>10:
            se_precip[cnt]=np.nan;              w_precip[cnt]=np.nan
        else:
            se_precip[cnt] = se_month_precip;   w_precip[cnt] = w_month_precip
        
        cnt+=1


accum_df = pd.DataFrame({   "Year":             years,
                            "Month":            months,
                            "AccumPrecip_SE":   se_precip,
                            "AccumPrecip_W":    w_precip
                        })

accum_df.to_csv(dir_name+"Test_Accum.csv")
                


