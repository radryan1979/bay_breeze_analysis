import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

BUOYdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/LoopCurrentPositionTest/"

yr=2016
buoy = "East"

EastSSTpath = BUOYdirpath + '20km_42098_{0}.csv'.format(yr)
MiddleSSTpath = BUOYdirpath + '100km_42022_{0}.csv'.format(yr)
WestSSTpath = BUOYdirpath + '150km_42099_{0}.csv'.format(yr)

EastSST_df = pd.read_csv(EastSSTpath,header=0, skiprows=[1],index_col=None)
EastSST_df = EastSST_df.reset_index(drop=True)
EastSST_df = EastSST_df.replace(999,np.NaN)
EastSST_df['Datetime'] = pd.to_datetime(EastSST_df['#YY'].apply(str) + '-' + EastSST_df['MM'].apply(str) + '-' + EastSST_df["DD"].apply(str) + " " +\
                        EastSST_df['hh'].apply(str) + ':' + EastSST_df['mm'].apply(str))
EastSST_df = EastSST_df.set_index('Datetime')
# EastSST_df = EastSST_df["WTMP"]
EastSST_df = EastSST_df.rename(columns={"WTMP":"42098"})
EastSST_df = EastSST_df.resample('D').mean()

MiddleSST_df = pd.read_csv(MiddleSSTpath,header=0,skiprows=[1],index_col=None)
MiddleSST_df = MiddleSST_df.reset_index(drop=True)
MiddleSST_df = MiddleSST_df.replace(999,np.NaN)
MiddleSST_df['Datetime'] = pd.to_datetime(MiddleSST_df['#YY'].apply(str) + '-' + MiddleSST_df['MM'].apply(str) + '-' + MiddleSST_df["DD"].apply(str) + " " +\
                        MiddleSST_df['hh'].apply(str) + ':' + MiddleSST_df['mm'].apply(str))
MiddleSST_df = MiddleSST_df.set_index('Datetime')
# MiddleSST_df = MiddleSST_df["WTMP"]
MiddleSST_df = MiddleSST_df.rename(columns={"WTMP":"42022"})
MiddleSST_df = MiddleSST_df.resample('D').mean()

WestSST_df = pd.read_csv(WestSSTpath,header=0,skiprows=[1],index_col=None)
WestSST_df = WestSST_df.reset_index(drop=True)
WestSST_df = WestSST_df.replace(999,np.NaN)
WestSST_df['Datetime'] = pd.to_datetime(WestSST_df['#YY'].apply(str) + '-' + WestSST_df['MM'].apply(str) + '-' + WestSST_df["DD"].apply(str) + " " +\
                        WestSST_df['hh'].apply(str) + ':' + WestSST_df['mm'].apply(str))
WestSST_df = WestSST_df.set_index('Datetime')
# WestSST_df = WestSST_df["WTMP"]
WestSST_df = WestSST_df.rename(columns={"WTMP":"42099"})
WestSST_df = WestSST_df.resample('D').mean()

totalSST_df=pd.concat([EastSST_df["42098"],MiddleSST_df["42022"],WestSST_df["42099"]], axis=1)

#Create monthly mask for study period
mask = (totalSST_df.index.month >= 4) & (totalSST_df.index.month <= 10)


# #Create mask for non-NAN values
# y_mask = np.isfinite(totalSST_df.loc[mask]["{0}SST".format(buoy)])
# y = totalSST_df["{0}SST".format(buoy)].loc[mask].loc[y_mask]



# #Create polyfit
# x = mdates.date2num(totalSST_df[mask][y_mask].index.tolist())

# z4 = np.polyfit(x, y, 2) #last number is degree of fit
# p4 = np.poly1d(z4) 

# xx = np.linspace(x.min(), x.max(), 100)
# dd = mdates.num2date(xx)




#PLOT
totalSST_df.loc[mask].plot(linewidth=0.35)
# plt.plot(dd,p4(xx),'-r',label="Trendline")
plt.legend()
plt.ylabel("SWT (\N{DEGREE SIGN}C)")
plt.title("{0} Florida West Coast SWT Analysis".format(yr))
plt.tight_layout()
plt.savefig("/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/LoopCurrentPositionTest/TestTrendline{0}.png".format(yr),dpi=500)
