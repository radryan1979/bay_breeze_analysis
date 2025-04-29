"""
Coded by: Dan Moore

The purpose of this code is to manipulate statistics created by UFLStats.py and to create visuals such as histograms, etc.
"""

import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt

path_name="/Volumes/LaCie/SeaBreeze/Statistics/"
suffix="UFLStats.csv"
infile=open(os.path.join\
(path_name,suffix),"r")

df=pd.read_csv(infile,sep=',',header=0,index_col=[0,1])

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
bins = np.linspace(0, 60, 40)
x=df['Precip12hr']

plt.hist(x,bins,alpha=1)
plt.xlabel('Precipitation [mm]')
plt.ylabel('Frequency')
plt.title('12 Hour Total Precipitation \nAfter Sea Breeze Passage')
plt.text(30, 6, u'\u03bc=%(mu)2.2f mm, \nMAD=%(MAD)2.2f mm' %\
         {'mu':df['Precip12hr'].mean(),'MAD':df['Precip12hr'].mad()})
plt.xlim(0)
plt.show()


#Wind directions before and after sea breeze passage
bins = np.linspace(0, 360, 45)
x=df['WindDirBef']
y=df['WindDirAft']

plt.hist(x,bins,alpha=0.5,label='Before Sea \nBreeze')
plt.hist(y,bins,alpha=0.5,label='After Sea \nBreeze')
plt.legend(loc='upper right')
plt.xlabel('Wind Direction [\N{DEGREE SIGN}]')
plt.ylabel('Frequency')
plt.title('Wind Direction Before vs. After\nSea Breeze Passage')
plt.text(10, 8, u'$\u03bc_b$=%(mu)2.2f\N{DEGREE SIGN}, \n$MAD_b$=%(MAD)2.2f\N{DEGREE SIGN}' %\
         {'mu':df['WindDirBef'].mean(),'MAD':df['WindDirBef'].mad()})
plt.text(10, 6, u'$\u03bc_a$=%(mu)2.2f\N{DEGREE SIGN}, \n$MAD_a$=%(MAD)2.2f\N{DEGREE SIGN}' %\
         {'mu':df['WindDirAft'].mean(),'MAD':df['WindDirAft'].mad()})
plt.xlim(0)
plt.show()

#Wind Change before-after
bins = np.linspace(0, 180, 25)
x=df['WindChange']

plt.hist(x,bins,alpha=1)
plt.xlabel('Wind Direction Before - Wind Direction After [\N{DEGREE SIGN}]')
plt.ylabel('Frequency')
plt.title('Wind Direction Change \nDuring Sea Breeze Passage')
plt.text(100, 6, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}' %\
         {'mu':df['WindChange'].mean(),'MAD':df['WindChange'].mad()})
plt.xlim(0,180)
plt.show()

#Temperature change before-after
bins = np.linspace(0, 8, 25)
x=df['TempChange']

plt.hist(x,bins,alpha=1)
plt.xlabel('Temperature Before - Temperature After [\N{DEGREE SIGN}C]')
plt.ylabel('Frequency')
plt.title('Temperature Change \nDuring Sea Breeze Passage')
plt.text(4, 12, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
         {'mu':df['TempChange'].mean(),'MAD':df['TempChange'].mad()})
plt.xlim(0)
plt.show()


#Wind Speed before and after sea breeze passage
bins = np.linspace(0, 7, 35)
x=df['WindSpeedBef']
y=df['WindSpeedAft']
plt.hist(x,bins,alpha=0.5,label='Before Sea \nBreeze')
plt.hist(y,bins,alpha=0.5,label='After Sea \nBreeze')
plt.legend(loc='upper right')
plt.xlabel('Wind Speed [m$s^-$$^1$]')
plt.ylabel('Frequency')
plt.title('Wind Speed Before vs. After \nSea Breeze Passage')
plt.text(4, 7.5, u'$\u03bc_b$=%(mu)2.2f m$s^-$$^1$, \n$MAD_b$=%(MAD)2.2fm$s^-$$^1$' %\
         {'mu':df['WindSpeedBef'].mean(),'MAD':df['WindSpeedBef'].mad()})
plt.text(4, 5.5, u'$\u03bc_a$=%(mu)2.2f m$s^-$$^1$, \n$MAD_a$=%(MAD)2.2fm$s^-$$^1$' %\
         {'mu':df['WindSpeedAft'].mean(),'MAD':df['WindSpeedAft'].mad()})
plt.xlim(0)
plt.show()

#Timing of precipitation
bins = np.linspace(0, 50, 10)
x=df['PrecipTiming']
plt.hist(x,bins,alpha=1)
plt.xlabel('Precipitation Timing [min]')
plt.ylabel('Frequency')
plt.title('Time of Precipitation \nAfter Sea Breeze Passage')
plt.text(10, 12, u'\u03bc=%(mu)2.2f min, \nMAD=%(MAD)2.2f min' %\
         {'mu':df['PrecipTiming'].mean(),'MAD':df['PrecipTiming'].mad()})
plt.xlim(0)
plt.show()


#Relative Humidity before and after sea breeze passage
bins = np.linspace(40, 100, 30)
x=df['RelHumBef']
y=df['RelHumAft']
plt.hist(x,bins,alpha=0.5,label='Before Sea \nBreeze')
plt.hist(y,bins,alpha=0.5,label='After Sea \nBreeze')
plt.legend(loc='upper left')
plt.xlabel('Relative Humidity [%]')
plt.ylabel('Frequency')
plt.title('Relative Humidity Before vs. After \nSea Breeze Passage')
plt.text(70, 12, u'$\u03bc_b$=%(mu)2.2f, \n$MAD_b$=%(MAD)2.2f' %\
         {'mu':df['RelHumBef'].mean(),'MAD':df['RelHumBef'].mad()})
plt.text(70, 10, u'$\u03bc_a$=%(mu)2.2f, \n$MAD_a$=%(MAD)2.2f' %\
         {'mu':df['RelHumAft'].mean(),'MAD':df['RelHumAft'].mad()})
plt.xlim(40)
plt.show()

#Temperature before and after sea breeze passage
bins = np.linspace(10, 35, 25)
x=df['TempBef']
y=df['TempAft']
plt.hist(x,bins,alpha=0.5,label='Before Sea \nBreeze')
plt.hist(y,bins,alpha=0.5,label='After Sea \nBreeze')
plt.legend(loc='upper left')
plt.xlabel('Temperature [\N{DEGREE SIGN}C]')
plt.ylabel('Frequency')
plt.title('Temperature Before vs. After \nSea Breeze Passage')
plt.text(12, 8, u'$\u03bc_b$=%(mu)2.2f\N{DEGREE SIGN}C, \n$MAD_b$=%(MAD)2.2f\N{DEGREE SIGN}C' %\
         {'mu':df['TempBef'].mean(),'MAD':df['TempBef'].mad()})
plt.text(12, 6, u'$\u03bc_a$=%(mu)2.2f\N{DEGREE SIGN}C, \n$MAD_a$=%(MAD)2.2f\N{DEGREE SIGN}C' %\
         {'mu':df['TempAft'].mean(),'MAD':df['TempAft'].mad()})
plt.xlim(10)
plt.show()



