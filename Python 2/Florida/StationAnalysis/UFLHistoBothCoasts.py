"""
Coded by: Dan Moore

The purpose of this code is to manipulate statistics created by UFLStats.py and to create visuals such as histograms, etc.
"""

import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt

path_name="/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results1_17_19/"
suffix="HughesBothCoasts.csv"
infile=open(os.path.join\
(path_name,suffix),"r")

df=pd.read_csv(infile,sep=',',header=0,index_col=False)

infile_east = open(os.path.join(path_name,"HughesEast_SST.csv"),"r")
dfeast = pd.read_csv(infile_east,sep=',',header=0,index_col=False)

infile_west = open(os.path.join(path_name,"HughesWest_SST.csv"),"r")
dfwest = pd.read_csv(infile_west,sep=',',header=0,index_col=False)

"""
#Precipitation totaled over 6 hours after sea breeze passage
pd.DataFrame.hist(df,column='Precip6hr')
plt.xlabel('Precipitation [mm]')
plt.ylabel('Frequency')
plt.title('6 Hour Total Precipitation After Sea Breeze Passage')
plt.text(20, 17.5, u'\u03bc=%(mu)2.2f mm, MAD=%(MAD)2.2f mm' %\
         {'mu':df['Precip6hr'].mean(),'MAD':df['Precip6hr'].mad()})
plt.xlim(0)
plt.show()

#Precipitation totaled over 24 hours after sea breeze passage
pd.DataFrame.hist(df,column='Precip24hr')
plt.xlabel('Precipitation [mm]')
plt.ylabel('Frequency')
plt.title('24 Hour Total Precipitation After Sea Breeze Passage')
plt.text(40, 25, u'\u03bc=%(mu)2.2f mm, MAD=%(MAD)2.2f mm' %\
         {'mu':df['Precip24hr'].mean(),'MAD':df['Precip24hr'].mad()})
plt.xlim(0)
plt.show()
"""







####Use these two for standard deviation in leiu of MAD
#\u03c3=%(sd)2.2f mm,
#'sd':df['Precip6hr'].std(),\

# #Precipitation totaled 12 hours after sea breeze passage
# bins = np.linspace(0, 60, 20)
# x=df['TwelveHrPrecip'].dropna()

# plt.hist(x,bins,alpha=1)
# plt.xlabel('Precipitation [mm]')
# plt.ylabel('Frequency')
# plt.title('Hughes Algorithm Florida Both Coast')
# plt.text(30, 1250, u'\u03bc=%(mu)2.2f mm, \nMAD=%(MAD)2.2f mm' %\
#          {'mu':df['TwelveHrPrecip'].mean(),'MAD':df['TwelveHrPrecip'].mad()})
# plt.xlim(0)
# plt.savefig(path_name+"/Histograms/12HrPrecipBoth.png",dpi=1200)
# plt.show()


# #Wind directions before and after sea breeze passage
# # bins = np.linspace(0, 360, 45)


# # plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
# # plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
# # plt.legend(loc='upper right')
# # plt.xlabel('Wind Direction [\N{DEGREE SIGN}]')
# # plt.ylabel('Frequency')
# # plt.title('Hughes Algorithm Florida Both Coasts')
# # plt.text(10, 135, u'$\u03bc_b$=%(mu)2.2f\N{DEGREE SIGN}, \n$MAD_b$=%(MAD)2.2f\N{DEGREE SIGN}' %\
# #          {'mu':df['WDBefore'].mean(),'MAD':df['WDBefore'].mad()})
# # plt.text(10, 110, u'$\u03bc_a$=%(mu)2.2f\N{DEGREE SIGN}, \n$MAD_a$=%(MAD)2.2f\N{DEGREE SIGN}' %\
# #          {'mu':df['WDAfter'].mean(),'MAD':df['WDAfter'].mad()})
# # plt.xlim(0)
# # plt.show()

# #Wind Change before-after
# # bins = np.linspace(0, 180, 25)
# # y=df['WDBefore'].dropna()
# # x=df['WDAfter'].dropna()
# # z=abs(y-x)

# # plt.hist(x,bins,alpha=1)
# # plt.xlabel('Wind Direction Before - Wind Direction After [\N{DEGREE SIGN}]')
# # plt.ylabel('Frequency')
# # plt.title('Wind Direction Change \nDuring Sea Breeze Passage')
# # plt.text(100, 50, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}' %\
# #          {'mu':df['WindChange'].mean(),'MAD':df['WindChange'].mad()})
# # plt.xlim(0,180)
# # plt.savefig(path_name+"/Histograms/WSChangeBoth.png",dpi=1200)
# # plt.show()

# # #Temperature change before-after
# # bins = np.linspace(0, 8, 25)
# # x=df['TempChange']

# # plt.hist(x,bins,alpha=1)
# # plt.xlabel('Temperature Before - Temperature After [\N{DEGREE SIGN}C]')
# # plt.ylabel('Frequency')
# # plt.title('Temperature Change \nDuring Sea Breeze Passage')
# # plt.text(4, 110, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
# #          {'mu':df['TempChange'].mean(),'MAD':df['TempChange'].mad()})
# # plt.xlim(0)
# # plt.show()


# #Wind Speed before and after sea breeze passage
# bins = np.linspace(0, 16, 30)
# y=df['WSBefore'].dropna()
# x=df['WSAfter'].dropna()
# plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
# plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
# plt.legend(loc='upper left')
# plt.xlabel('Wind Speed [mi $hr^-$$^1$]')
# plt.ylabel('Frequency')
# plt.title('Hughes Algorithm Florida Both Coasts')
# plt.text(11, 500, u'$\u03bc_b$=%(mu)2.2f mi $hr^-$$^1$, \n$MAD_b$=%(MAD)2.2fmi $hr^-$$^1$' %\
#          {'mu':df['WSBefore'].mean(),'MAD':df['WSBefore'].mad()})
# plt.text(11, 350, u'$\u03bc_a$=%(mu)2.2f mi $hr^-$$^1$, \n$MAD_a$=%(MAD)2.2fmi $hr^-$$^1$' %\
#          {'mu':df['WSAfter'].mean(),'MAD':df['WSAfter'].mad()})
# plt.xlim(0)
# plt.savefig(path_name+"/Histograms/WSBoth.png",dpi=1200)
# plt.show()

# # #Timing of precipitation
# # bins = np.linspace(0, 50, 10)
# # x=df['PrecipTiming']
# # plt.hist(x,bins,alpha=1)
# # plt.xlabel('Precipitation Timing [min]')
# # plt.ylabel('Frequency')
# # plt.title('Time of Precipitation \nAfter Sea Breeze Passage')
# # plt.text(10, 40, u'\u03bc=%(mu)2.2f min, \nMAD=%(MAD)2.2f min' %\
# #          {'mu':df['PrecipTiming'].mean(),'MAD':df['PrecipTiming'].mad()})
# # plt.xlim(0)
# # plt.show()

# #Timing of SBF Passage
# bins = np.linspace(0, 24, 24)
# x=df['Hour']

# plt.hist(x,bins,alpha=1)
# plt.xlabel('Time of SBF Passage [LST]')
# plt.ylabel('Frequency')
# plt.title('Hughes Algorithm Florida Both Coasts')
# plt.text(1, 300, u'\u03bc=%(mu)2.2f LST, \nMAD=%(MAD)2.2f LST' %\
#          {'mu':df['Hour'].mean(),'MAD':df['Hour'].mad()})
# plt.xlim(0)
# plt.savefig(path_name+"/Histograms/SBFTimingBoth.png",dpi=1200)
# plt.show()


# #Relative Humidity before and after sea breeze passage
# bins = np.linspace(40, 100, 30)
# y=df['RHBefore'].dropna()
# x=df['RHAfter'].dropna()
# plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
# plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
# plt.legend(loc='upper left')
# plt.xlabel('Relative Humidity [%]')
# plt.ylabel('Frequency')
# plt.title('Hughes Algorithm Florida Both Coasts')
# plt.text(42, 150, u'$\u03bc_b$=%(mu)2.2f, \n$MAD_b$=%(MAD)2.2f' %\
#          {'mu':df['RHBefore'].mean(),'MAD':df['RHBefore'].mad()})
# plt.text(42, 100, u'$\u03bc_a$=%(mu)2.2f, \n$MAD_a$=%(MAD)2.2f' %\
#          {'mu':df['RHAfter'].mean(),'MAD':df['RHAfter'].mad()})
# plt.xlim(40)
# plt.savefig(path_name+"/Histograms/RHBoth.png",dpi=1200)
# plt.show()

# #Temperature before and after sea breeze passage
# bins = np.linspace(10, 35, 25)
# y=df['TBefore'].dropna()
# x=df['TAfter'].dropna()
# plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
# plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
# plt.legend(loc='upper left')
# plt.xlabel('Temperature [\N{DEGREE SIGN}C]')
# plt.ylabel('Frequency')
# plt.title('Hughes Algorithm Florida Both Coasts')
# plt.text(12, 200, u'$\u03bc_b$=%(mu)2.2f\N{DEGREE SIGN}C, \n$MAD_b$=%(MAD)2.2f\N{DEGREE SIGN}C' %\
#          {'mu':df['TBefore'].mean(),'MAD':df['TBefore'].mad()})
# plt.text(12, 100, u'$\u03bc_a$=%(mu)2.2f\N{DEGREE SIGN}C, \n$MAD_a$=%(MAD)2.2f\N{DEGREE SIGN}C' %\
#          {'mu':df['TAfter'].mean(),'MAD':df['TAfter'].mad()})
# plt.xlim(10)
# # plt.savefig(path_name+"/Histograms/TempBoth.png",dpi=1200)
# plt.show()









#Frequency by month
dfeast["Date"] = dfeast["Date"].astype("datetime64")
east_freq_month = dfeast["Date"].groupby(dfeast["Date"].dt.month).count()

dfwest["Date"] = dfwest["Date"].astype("datetime64")
west_freq_month = dfwest["Date"].groupby(dfwest["Date"].dt.month).count()

wid=0.25
plt.bar(west_freq_month.index-wid/2.0,west_freq_month, color = 'b', width = wid, label="West")
plt.bar(east_freq_month.index+wid/2.0,east_freq_month, color = 'r', width = wid, label="East")

plt.legend(loc='best')
plt.title("Sea Breeze Frequency by Month")
plt.xlabel("Month"); plt.ylabel("Frequency")

plt.savefig(path_name+"Figures/SBC_Frequency_Month.png", dpi=500)
plt.show()


#Frequency by year
dfeast["Date"] = dfeast["Date"].astype("datetime64")
east_freq_year = dfeast["Date"].groupby(dfeast["Date"].dt.year).count()

dfwest["Date"] = dfwest["Date"].astype("datetime64")
west_freq_year = dfwest["Date"].groupby(dfwest["Date"].dt.year).count()

wid=0.25
plt.bar(west_freq_year.index-wid/2.0,west_freq_year, color = 'b', width = wid, label="West")
plt.bar(east_freq_year.index+wid/2.0,east_freq_year, color = 'r', width = wid, label="East")

plt.legend(loc='best')
plt.title("Sea Breeze Frequency by Year")
plt.xlabel("Year"); plt.ylabel("Frequency")

plt.savefig(path_name+"Figures/SBC_Frequency_Year.png", dpi=500)
plt.show()


#Rainfall Frequency both coasts
df_east_intense = dfeast[dfeast["TwelveHrPrecip"]>3]
east_intense_precip = df_east_intense["Date"].groupby(df_east_intense["Date"].dt.month).count()

df_west_intense = dfwest[dfwest["TwelveHrPrecip"]>3]
west_intense_precip = df_west_intense["Date"].groupby(df_west_intense["Date"].dt.month).count()

df_east_less = dfeast[((dfeast["TwelveHrPrecip"]>0.05) & (dfeast["TwelveHrPrecip"]<3.0))]
east_less_precip = df_east_less["Date"].groupby(df_east_less["Date"].dt.month).count()

df_west_less = dfwest[((dfwest["TwelveHrPrecip"]>0.05) & (dfwest["TwelveHrPrecip"]<3.0))]
west_less_precip = df_west_less["Date"].groupby(df_west_less["Date"].dt.month).count()







fig, ax = plt.subplots()

bar_width = 0.25

ax.bar(west_intense_precip.index-bar_width/2.0,west_intense_precip,
                width = bar_width,
                color='b', label='West')

ax.bar(east_intense_precip.index+bar_width/2.0,east_intense_precip,
                width = bar_width,
                color='r', label='East')

ax.set_xlabel('Month')
ax.set_ylabel('Frequency')
ax.set_title('Intense (>3in) 12Hr Rainfall')
ax.legend()

fig.tight_layout()
plt.savefig(path_name+"Figures/StationIntensePrecip.png", dpi=500)
plt.show()



fig, ax = plt.subplots()

bar_width = 0.25

ax.bar(west_less_precip.index-bar_width/2.0,west_less_precip,
                width = bar_width,
                color='b', label='West')

ax.bar(east_less_precip.index+bar_width/2.0,east_less_precip,
                width = bar_width,
                color='r', label='East')

ax.set_xlabel('Month')
ax.set_ylabel('Frequency')
ax.set_title('Less (<3in) 12Hr Rainfall')
ax.legend()

fig.tight_layout()
plt.savefig(path_name+"Figures/StationLessPrecip.png", dpi=500)
plt.show()


