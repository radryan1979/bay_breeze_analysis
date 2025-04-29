"""
Coded by: Dan Moore

The purpose of this code is to manipulate statistics created by UFLStats.py and to create visuals such as histograms, etc.
"""

import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt

path_name="/Volumes/LaCie/SeaBreeze/Florida/Results1_17_19/"
suffix="HughesEast_SST.csv"
infile=open(os.path.join\
(path_name,suffix),"r")

df=pd.read_csv(infile,sep=',',header=0,index_col=False)

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

#Precipitation totaled 12 hours after sea breeze passage
bins = np.linspace(0, 60, 20)
x=df['TwelveHrPrecip'].dropna()

plt.hist(x,bins,alpha=1)
plt.xlabel('Precipitation [mm]')
plt.ylabel('Frequency')
plt.title('Hughes Algorithm Florida East Coast')
plt.text(30, 1250, u'\u03bc=%(mu)2.2f mm, \nMAD=%(MAD)2.2f mm' %\
         {'mu':df['TwelveHrPrecip'].mean(),'MAD':df['TwelveHrPrecip'].mad()})
plt.xlim(0)
plt.savefig(path_name+"/Histograms/12HrPrecipEast.png",dpi=1200)
plt.show()

#Land-Sea Temp gradient before SBC
bins = np.linspace(-10, 12, 25)
x=df['TBefore']-df["SST_East"]
x=x.dropna()
plt.hist(x,bins,alpha=1)
plt.xlabel('Land-Sea Temperature Difference [\N{DEGREE SIGN}C]')
plt.ylabel('Frequency')
plt.title('Hughes Algorithm Florida East Coast')
plt.text(-8, 400, u'\u03bc=%(mu)2.2f \N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f \N{DEGREE SIGN}C' %\
         {'mu':x.mean(),'MAD':x.mad()})
# plt.savefig(path_name+"/Histograms/AirSeaGradientEast.png",dpi=500)
plt.show()


#Wind directions before and after sea breeze passage
bins = np.linspace(0, 360, 45)
y=df['WDBefore'].dropna()
x=df['WDAfter'].dropna()

plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
plt.legend(loc='upper right')
plt.xlabel('Wind Direction [\N{DEGREE SIGN}]')
plt.ylabel('Frequency')
plt.title('Hughes Algorithm Florida East Coast')
plt.text(280, 135, u'$\u03bc_b$=%(mu)2.2f\N{DEGREE SIGN}, \n$MAD_b$=%(MAD)2.2f\N{DEGREE SIGN}' %\
         {'mu':df['WDBefore'].mean(),'MAD':df['WDBefore'].mad()})
plt.text(280, 90, u'$\u03bc_a$=%(mu)2.2f\N{DEGREE SIGN}, \n$MAD_a$=%(MAD)2.2f\N{DEGREE SIGN}' %\
         {'mu':df['WDAfter'].mean(),'MAD':df['WDAfter'].mad()})
plt.xlim(0)
plt.savefig(path_name+"/Histograms/WDEast.png",dpi=1200)
plt.show()

# bins = np.linspace(0, 180, 20)
# y=df['WDBefore'].dropna()
# x=df['WDAfter'].dropna()
# z=abs(x-y).dropna()

# plt.hist(z,bins,alpha=1)
# plt.xlabel('Wind Direction Before - Wind Direction After [\N{DEGREE SIGN}]')
# plt.ylabel('Frequency')
# plt.title('Wind Direction Change \nDuring Sea Breeze Passage')
# plt.text(100, 250, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}' %\
#          {'mu':z.mean(),'MAD':z.mad()})
# plt.xlim(0,180)
# plt.show()

# #Temperature change before-after
# bins = np.linspace(0, 8, 25)
# x=df['TempChange']

# plt.hist(x,bins,alpha=1)
# plt.xlabel('Temperature Before - Temperature After [\N{DEGREE SIGN}C]')
# plt.ylabel('Frequency')
# plt.title('Temperature Change \nDuring Sea Breeze Passage')
# plt.text(4, 110, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
#          {'mu':df['TempChange'].mean(),'MAD':df['TempChange'].mad()})
# plt.xlim(0)
# plt.show()


#Wind Speed before and after sea breeze passage
bins = np.linspace(0, 16, 35)
y=df['WSBefore'].dropna()
x=df['WSAfter'].dropna()
plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
plt.legend(loc='upper left')
plt.xlabel('Wind Speed [mi $hr^-$$^1$]')
plt.ylabel('Frequency')
plt.title('Hughes Algorithm Florida East Coast')
plt.text(0.25, 70, u'$\u03bc_b$=%(mu)2.2f mi $hr^-$$^1$, \n$MAD_b$=%(MAD)2.2fmi $hr^-$$^1$' %\
         {'mu':df['WSBefore'].mean(),'MAD':df['WSBefore'].mad()})
plt.text(0.25, 40, u'$\u03bc_a$=%(mu)2.2f mi $hr^-$$^1$, \n$MAD_a$=%(MAD)2.2fmi $hr^-$$^1$' %\
         {'mu':df['WSAfter'].mean(),'MAD':df['WSAfter'].mad()})
plt.xlim(0)
plt.savefig(path_name+"/Histograms/WSEast.png",dpi=1200)
plt.show()

# #Timing of precipitation
# bins = np.linspace(0, 50, 10)
# x=df['PrecipTiming']
# plt.hist(x,bins,alpha=1)
# plt.xlabel('Precipitation Timing [min]')
# plt.ylabel('Frequency')
# plt.title('Time of Precipitation \nAfter Sea Breeze Passage')
# plt.text(10, 40, u'\u03bc=%(mu)2.2f min, \nMAD=%(MAD)2.2f min' %\
#          {'mu':df['PrecipTiming'].mean(),'MAD':df['PrecipTiming'].mad()})
# plt.xlim(0)
# plt.show()

#Timing of SBF Passage
bins = np.linspace(0, 24, 24)
x=df['Hour']

plt.hist(x,bins,alpha=1)
plt.xlabel('Time of SBF Passage [LST]')
plt.ylabel('Frequency')
plt.title('Hughes Algorithm Florida East Coast')
plt.text(1, 300, u'\u03bc=%(mu)2.2f LST, \nMAD=%(MAD)2.2f LST' %\
         {'mu':df['Hour'].mean(),'MAD':df['Hour'].mad()})
plt.xlim(0)
plt.savefig(path_name+"/Histograms/SBFTimingEast.png",dpi=1200)
plt.show()

#Relative Humidity before and after sea breeze passage
bins = np.linspace(40, 100, 30)
y=df['RHBefore'].dropna()
x=df['RHAfter'].dropna()
plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
plt.legend(loc='upper left')
plt.xlabel('Relative Humidity [%]')
plt.ylabel('Frequency')
plt.title('Hughes Algorithm Florida East Coast')
plt.text(42, 160, u'$\u03bc_b$=%(mu)2.2f, \n$MAD_b$=%(MAD)2.2f' %\
         {'mu':df['RHBefore'].mean(),'MAD':df['RHBefore'].mad()})
plt.text(42, 100, u'$\u03bc_a$=%(mu)2.2f, \n$MAD_a$=%(MAD)2.2f' %\
         {'mu':df['RHAfter'].mean(),'MAD':df['RHAfter'].mad()})
plt.xlim(40)
plt.savefig(path_name+"/Histograms/RHEast.png",dpi=1200)
plt.show()

#Temperature before and after sea breeze passage
bins = np.linspace(10, 35, 25)
y=df['TBefore'].dropna()
x=df['TAfter'].dropna()
plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
plt.legend(loc='upper left')
plt.xlabel('Temperature [\N{DEGREE SIGN}C]')
plt.ylabel('Frequency')
plt.title('Hughes Algorithm Florida East Coast')
plt.text(12, 250, u'$\u03bc_b$=%(mu)2.2f\N{DEGREE SIGN}C, \n$MAD_b$=%(MAD)2.2f\N{DEGREE SIGN}C' %\
         {'mu':df['TBefore'].mean(),'MAD':df['TBefore'].mad()})
plt.text(12, 100, u'$\u03bc_a$=%(mu)2.2f\N{DEGREE SIGN}C, \n$MAD_a$=%(MAD)2.2f\N{DEGREE SIGN}C' %\
         {'mu':df['TAfter'].mean(),'MAD':df['TAfter'].mad()})
plt.xlim(10)
plt.savefig(path_name+"/Histograms/TempEast.png",dpi=1200)
plt.show()



df["Date"] = df["Date"].astype("datetime64")

#SBC SST by month
for i in range(4,11):
    bins = np.linspace(20, 32, 30)
    x=df.loc[df["Date"].dt.month==i]["SST_East"].dropna()
    plt.hist(x,bins,alpha=1)
    plt.xlabel('SST [\N{DEGREE SIGN}C]')
    plt.ylabel('Frequency')
    plt.title('East Hughes Algorithm Florida East Coast \nSST Month: {0}'.format(i))
    plt.text(20, 20, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
             {'mu':df.loc[df["Date"].dt.month==i]["SST_East"].mean(),\
             'MAD':df.loc[df["Date"].dt.month==i]["SST_East"].mad()})
    # plt.xlim(0)
    plt.savefig(path_name+"/Histograms/EastCoastSBC_SST/EastSBC_SSTEast_{0}.png".format(i),dpi=500)
    plt.show()


for i in range(4,11):
    bins = np.linspace(18, 33, 30)
    x=df.loc[df["Date"].dt.month==i]["SST_West"].dropna()
    plt.hist(x,bins,alpha=1)
    plt.xlabel('SST [\N{DEGREE SIGN}C]')
    plt.ylabel('Frequency')
    plt.title('East Hughes Algorithm Florida West Coast \nSST Month: {0}'.format(i))
    plt.text(20, 20, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
             {'mu':df.loc[df["Date"].dt.month==i]["SST_West"].mean(),\
             'MAD':df.loc[df["Date"].dt.month==i]["SST_West"].mad()})
    # plt.xlim(0)
    plt.savefig(path_name+"/Histograms/EastCoastSBC_SST/EastSBC_SSTWest_{0}.png".format(i),dpi=500)
    plt.show()

for i in range(4,11):
    bins = np.linspace(20, 36, 30)
    x=df.loc[df["Date"].dt.month==i]["SST_South_InsideKeys"].dropna()
    plt.hist(x,bins,alpha=1)
    plt.xlabel('SST [\N{DEGREE SIGN}C]')
    plt.ylabel('Frequency')
    plt.title('East Hughes Algorithm Florida South Coast \nInside Keys SST Month: {0}'.format(i))
    plt.text(20, 20, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
             {'mu':df.loc[df["Date"].dt.month==i]["SST_South_InsideKeys"].mean(),\
             'MAD':df.loc[df["Date"].dt.month==i]["SST_South_InsideKeys"].mad()})
    # plt.xlim(0)
    plt.savefig(path_name+"/Histograms/EastCoastSBC_SST/EastSBC_SSTSouthInsideKeys_{0}.png".format(i),dpi=500)
    plt.show()
    
for i in range(4,11):
    bins = np.linspace(20, 32, 30)
    x=df.loc[df["Date"].dt.month==i]["SST_South_OutsideKeys"].dropna()
    plt.hist(x,bins,alpha=1)
    plt.xlabel('SST [\N{DEGREE SIGN}C]')
    plt.ylabel('Frequency')
    plt.title('East Hughes Algorithm Florida South Coast \nOutside Keys SST Month: {0}'.format(i))
    plt.text(20, 20, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
             {'mu':df.loc[df["Date"].dt.month==i]["SST_South_OutsideKeys"].mean(),\
             'MAD':df.loc[df["Date"].dt.month==i]["SST_South_OutsideKeys"].mad()})
    # plt.xlim(0)
    plt.savefig(path_name+"/Histograms/EastCoastSBC_SST/EastSBC_SSTSouthOutsideKeys_{0}.png".format(i),dpi=500)
    plt.show()


