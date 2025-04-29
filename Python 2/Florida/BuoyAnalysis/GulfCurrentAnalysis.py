import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

BUOYdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/GulfCurrentPositionTest/"

yr=2017
buoy = "East"

EastSSTpath = BUOYdirpath + '180km_41010_{0}.csv'.format(yr)
MiddleSSTpath = BUOYdirpath + '30km_41009_{0}.csv'.format(yr)
WestSSTpath = BUOYdirpath + '1km_41113_{0}.csv'.format(yr)
# SouthSSTpath = BUOYdirpath + 'S_MLRF1_{0}.csv'.format(yr)

#For analyzing up and down coast:
# MiddleSSTpath = BUOYdirpath + 'Middle_LKWF1_{0}.csv'.format(yr)
# SouthSSTpath = BUOYdirpath + 'South_FWYF1_{0}.csv'.format(yr)

EastSST_df = pd.read_csv(EastSSTpath,header=0, skiprows=[1],index_col=None)
EastSST_df = EastSST_df.reset_index(drop=True)
EastSST_df = EastSST_df.replace(999,np.NaN)
EastSST_df['Datetime'] = pd.to_datetime(EastSST_df['#YY'].apply(str) + '-' + EastSST_df['MM'].apply(str) + '-' + EastSST_df["DD"].apply(str) + " " +\
                        EastSST_df['hh'].apply(str) + ':' + EastSST_df['mm'].apply(str))
EastSST_df = EastSST_df.set_index('Datetime')
# EastSST_df = EastSST_df["WTMP"]
EastSST_df = EastSST_df.rename(columns={"WTMP":"41010"})
EastSST_df = EastSST_df.resample('D').mean()

MiddleSST_df = pd.read_csv(MiddleSSTpath,header=0,skiprows=[1],index_col=None)
MiddleSST_df = MiddleSST_df.reset_index(drop=True)
MiddleSST_df = MiddleSST_df.replace(999,np.NaN)
MiddleSST_df['Datetime'] = pd.to_datetime(MiddleSST_df['#YY'].apply(str) + '-' + MiddleSST_df['MM'].apply(str) + '-' + MiddleSST_df["DD"].apply(str) + " " +\
                        MiddleSST_df['hh'].apply(str) + ':' + MiddleSST_df['mm'].apply(str))
MiddleSST_df = MiddleSST_df.set_index('Datetime')
# MiddleSST_df = MiddleSST_df["WTMP"]
MiddleSST_df = MiddleSST_df.rename(columns={"WTMP":"41009"})
MiddleSST_df = MiddleSST_df.resample('D').mean()

WestSST_df = pd.read_csv(WestSSTpath,header=0,skiprows=[1],index_col=None)
WestSST_df = WestSST_df.reset_index(drop=True)
WestSST_df = WestSST_df.replace(999,np.NaN)
WestSST_df['Datetime'] = pd.to_datetime(WestSST_df['#YY'].apply(str) + '-' + WestSST_df['MM'].apply(str) + '-' + WestSST_df["DD"].apply(str) + " " +\
                        WestSST_df['hh'].apply(str) + ':' + WestSST_df['mm'].apply(str))
WestSST_df = WestSST_df.set_index('Datetime')
# WestSST_df = WestSST_df["WTMP"]
WestSST_df = WestSST_df.rename(columns={"WTMP":"41113"})
WestSST_df = WestSST_df.resample('D').mean()

# SouthSST_df = pd.read_csv(SouthSSTpath,header=0, skiprows=[1],index_col=None)
# SouthSST_df = SouthSST_df.reset_index(drop=True)
# SouthSST_df = SouthSST_df.replace(999,np.NaN)
# SouthSST_df['Datetime'] = pd.to_datetime(SouthSST_df['#YY'].apply(str) + '-' + SouthSST_df['MM'].apply(str) + '-' + SouthSST_df["DD"].apply(str) + " " +\
#                         SouthSST_df['hh'].apply(str) + ':' + SouthSST_df['mm'].apply(str))
# SouthSST_df = SouthSST_df.set_index('Datetime')
# # SouthSST_df = SouthSST_df["WTMP"]
# SouthSST_df = SouthSST_df.rename(columns={"WTMP":"SouthSST"})
# SouthSST_df = SouthSST_df.resample('H').mean()

totalSST_df=pd.concat([EastSST_df["41010"],MiddleSST_df["41009"],WestSST_df["41113"]], axis=1)

#Create monthly mask for study period
mask = (totalSST_df.index.month >= 4) & (totalSST_df.index.month <= 10)


#Create mask for non-NAN values
# y_mask = np.isfinite(totalSST_df.loc[mask]["{0}SST".format(buoy)])
# y = totalSST_df["{0}SST".format(buoy)].loc[mask].loc[y_mask]



# #Create polyfit
# x = mdates.date2num(totalSST_df[mask][y_mask].index.tolist())

# z4 = np.polyfit(x, y, 3) #last number is degree of fit
# p4 = np.poly1d(z4) 

# xx = np.linspace(x.min(), x.max(), 100)
# dd = mdates.num2date(xx)




#PLOT
# plt.plot(totalSST_df["EastSST"].loc[mask].index,totalSST_df["EastSST"].loc[mask], linewidth=0.35, linestyle='-', color='b')
# plt.plot(totalSST_df["MiddleSST"].loc[mask].index,totalSST_df["MiddleSST"].loc[mask], linewidth=0.35, linestyle='-', color='k')
# plt.plot(totalSST_df["WestSST"].loc[mask].index,totalSST_df["WestSST"].loc[mask], linewidth=0.35, linestyle='-', color='g')
# plt.plot(totalSST_df["SouthSST"].loc[mask].index,totalSST_df["SouthSST"].loc[mask], linewidth=0.35, linestyle='-', color='m')
# plt.plot(dd,p4(xx),'-r',label="Trendline")
totalSST_df.loc[mask].plot(linewidth=0.35)
plt.legend()
plt.ylabel("SWT (\N{DEGREE SIGN}C)")
plt.title("2017 Florida East Coast SWT Analysis")
# plt.savefig("/Users/dpm/Desktop/Test69.png", dpi=500)
plt.tight_layout()
plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/GulfCurrentPositionTest/TestTrendline{0}.png".format(yr),dpi=500)
