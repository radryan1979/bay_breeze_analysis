"""
Coded by: Dan Moore

The purpose of this code is to manipulate statistics for FL sea breeze detection.

Updated: 4-3-19
"""

import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt
import math
from scipy import stats

path_name="/Volumes/LaCie/SeaBreeze/Florida/StationDetection/Results3_6_19"
suffix="Hughes.csv"
infile=open(os.path.join\
(path_name,suffix),"r")

df=pd.read_csv(infile,sep=',',header=0,index_col=False)
df["Date"] = df["Date"].astype("datetime64")

classic = df[df["Type"]=="classic"]; weak = df[df["Type"]=="weak"]
classicDP = df[df["Type"]=="classicDP"]; DPWS = df[df["Type"]=="DPWS"]
classicWS = df[df["Type"]=="classicWS"];

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
types=["classic"]#,"classicDP","classicWS","weak","DPWS"]
east=[340,410,420,440]
west=[350,360,380,450,480,490]

for ty in types:
    for i in range(len(west)):
            
        temp_df = df[df["Type"]=="{0}".format(ty)]
        
        temp_df = temp_df[np.isin(temp_df["Station"],west[i])]
        
        # if i==0:
        #     coast="East"
        #     temp_df = temp_df[np.isin(temp_df["Station"],east)]
        # else:
        #     coast="West"
        #     temp_df = temp_df[np.isin(temp_df["Station"],west)]
        temp_df = temp_df.reset_index()
    
        # #Precipitation totaled 12 hours after sea breeze passage
        # bins = np.linspace(0, 60, 20)
        # x=temp_df['TwelveHrPrecip'].dropna()
        # n1,_,_=plt.hist(x,bins,alpha=1)
        # plt.xlabel('Precipitation [mm]')
        # plt.ylabel('Frequency')
        # plt.title('Precipitation SBC Type {0}\nFlorida {1} Coast'.format(ty,coast))
        # plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/3.0, u'\u03bc=%(mu)2.2f mm, \nMAD=%(MAD)2.2f mm' %\
        #          {'mu':temp_df['TwelveHrPrecip'].mean(),'MAD':temp_df['TwelveHrPrecip'].mad()})
        # plt.xlim(0)
        # plt.savefig(path_name+"/Figures/Histograms/{0}Coast/12HrPrecip_{1}.png".format(coast,ty),dpi=500)
        # plt.show()
        
        
        # Wind directions before and after sea breeze passage
        # numer_bef=0.0
        # denom_bef=0.0
        # numer_aft=0.0
        # denom_aft=0.0
        # for j in range(len(temp_df)):
            
        #     numer_bef += np.sin(temp_df["WDBefore"][j]*np.pi/180.0)
        #     denom_bef += np.cos(temp_df["WDBefore"][j]*np.pi/180.0)
        #     numer_aft += np.sin(temp_df["WDAfter"][j]*np.pi/180.0)
        #     denom_aft += np.cos(temp_df["WDAfter"][j]*np.pi/180.0)
            
        # numer_bef = numer_bef/float(len(temp_df))
        # denom_bef = denom_bef/float(len(temp_df))
        # numer_aft = numer_aft/float(len(temp_df))
        # denom_bef = denom_bef/float(len(temp_df))
        # r_bef = np.sqrt(np.square(numer_bef)+np.square(denom_bef))
        # r_aft = np.sqrt(np.square(numer_aft)+np.square(denom_aft))
        # numer_bef = numer_bef/r_bef
        # denom_bef = denom_bef/r_bef
        # numer_aft = numer_aft/r_aft
        # denom_aft = denom_aft/r_aft
        # mean_dir_bef = np.arctan(numer_bef/denom_bef)*180/np.pi
        # mean_dir_aft = np.arctan(numer_aft/denom_aft)*180/np.pi
        
        # if numer_bef>0 and denom_bef<0:
        #     mean_dir_bef = 180-mean_dir_bef
        # elif numer_bef<0 and denom_bef>0:
        #     mean_dir_bef = 360-mean_dir_bef
        # elif numer_bef<0 and denom_bef<0:
        #     mean_dir_bef = 180+mean_dir_bef
            
        # if numer_aft>0 and denom_aft<0:
        #     mean_dir_aft = 180-mean_dir_aft
        # elif numer_aft<0 and denom_aft>0:
        #     mean_dir_aft = 360-mean_dir_aft
        # elif numer_aft<0 and denom_aft<0:
        #     mean_dir_aft = 180+mean_dir_aft
            
        # cv_bef = 1.0 - (np.sqrt(np.square(numer_bef)+np.square(denom_bef))/float(len(temp_df)))
        # cv_aft = 1.0 - (np.sqrt(np.square(numer_aft)+np.square(denom_aft))/float(len(temp_df)))
        
        """
        cosSum_bef = 0.0  
        sinSum_bef = 0.0 
        cosSum_aft = 0.0  
        sinSum_aft = 0.0 
        angles_bef = np.array(temp_df["WDBefore"])
        angles_aft = np.array(temp_df["WDAfter"])
        for i in range(len(angles_aft)):  
            theCos_bef = math.cos(math.radians(float(angles_bef[i])))  
            theSin_bef = math.sin(math.radians(float(angles_bef[i])))  
            cosSum_bef += theCos_bef  
            sinSum_bef += theSin_bef  
            theCos_aft = math.cos(math.radians(float(angles_aft[i])))  
            theSin_aft = math.sin(math.radians(float(angles_aft[i])))  
            cosSum_aft += theCos_aft 
            sinSum_aft += theSin_aft  
        N = len(angles_aft)  
        C_bef = cosSum_bef/N  
        S_bef = sinSum_bef/N  
        C_aft = cosSum_aft/N  
        S_aft = sinSum_aft/N  
        theMean_bef = math.atan2(S_bef,C_bef)  
        theMean_aft = math.atan2(S_aft,C_aft)  
        if theMean_bef < 0.0:  
            theMean_bef += math.radians(360.0)  
        if theMean_aft < 0.0:  
            theMean_aft += math.radians(360.0)      
        theMean_bef = theMean_bef*180.0/np.pi
        theMean_aft = theMean_aft*180.0/np.pi
        
        cv_bef = 1.0 - (np.sqrt(np.square(sinSum_bef)+np.square(cosSum_bef))/N)
        cv_aft = 1.0 - (np.sqrt(np.square(sinSum_aft)+np.square(cosSum_aft))/N)
        
        bins = np.linspace(0, 360, 45)
        y=temp_df['WDBefore'].dropna()
        x=temp_df['WDAfter'].dropna()
        n1,_,_=plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
        n2,_,_=plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
        plt.legend(loc='best')
        plt.xlabel('Wind Direction [\N{DEGREE SIGN}]')
        plt.ylabel('Frequency')
        plt.title('Wind Direction SBC Type {0}\nFlorida {1} Coast'.format(ty,coast))
        plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/3.0, u'$\u03bc_b$=%(mu)2.2f\N{DEGREE SIGN}, \n$CV_b$=%(CV)2.2f' %\
                 {'mu':theMean_bef,'CV':cv_bef})
        plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/6.0, u'$\u03bc_a$=%(mu)2.2f\N{DEGREE SIGN}, \n$CV_a$=%(CV)2.2f' %\
                 {'mu':theMean_aft,'CV':cv_aft})
        plt.xlim(0)
        plt.savefig(path_name+"/Figures/Histograms/{0}Coast/WindDir_{1}.png".format(coast,ty),dpi=500)
        plt.show()
        
        
        # Wind Change before-after
        bins = np.linspace(0, 181, 25)
        y=df['WDBefore'].dropna()
        x=df['WDAfter'].dropna()
        z=abs(y-x)
        plt.hist(x,bins,alpha=1)
        plt.xlabel('Wind Direction Before - Wind Direction After [\N{DEGREE SIGN}]')
        plt.ylabel('Frequency')
        plt.title('Wind Direction Change \nDuring Sea Breeze Passage')
        plt.text(100, 50, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}' %\
                 {'mu':df['WindChange'].mean(),'MAD':df['WindChange'].mad()})
        plt.xlim(0,180)
        # plt.savefig(path_name+"/Figures/Histograms/{0}Coast/WSChange_{1}.png".format(coast,ty),dpi=500)
        plt.show()
        
        #Temperature change before-after
        bins = np.linspace(0, 8, 25)
        x=df['TempChange']
        plt.hist(x,bins,alpha=1)
        plt.xlabel('Temperature Before - Temperature After [\N{DEGREE SIGN}C]')
        plt.ylabel('Frequency')
        plt.title('Temperature Change \nDuring Sea Breeze Passage')
        plt.text(4, 110, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
                 {'mu':df['TempChange'].mean(),'MAD':df['TempChange'].mad()})
        plt.xlim(0)
        # plt.savefig(path_name+"/Figures/Histograms/{0}Coast/TempChange_{2}.png".format(coast,ty),dpi=500)
        plt.show()
        
        #Wind Speed before and after sea breeze passage
        bins = np.linspace(0, 16, 30)
        y=temp_df['WSBefore'].dropna()
        x=temp_df['WSAfter'].dropna()
        n1,_,_=plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
        n2,_,_=plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
        plt.legend(loc='best')
        plt.xlabel('Wind Speed [mi $hr^-$$^1$]')
        plt.ylabel('Frequency')
        plt.title('Wind Speed SBC Type {0}\nFlorida {1} Coast'.format(ty,coast))
        plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/3.0, u'$\u03bc_b$=%(mu)2.2f mi $hr^-$$^1$, \n$MAD_b$=%(MAD)2.2fmi $hr^-$$^1$' %\
                 {'mu':temp_df['WSBefore'].mean(),'MAD':temp_df['WSBefore'].mad()})
        plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/6.0, u'$\u03bc_a$=%(mu)2.2f mi $hr^-$$^1$, \n$MAD_a$=%(MAD)2.2fmi $hr^-$$^1$' %\
                 {'mu':temp_df['WSAfter'].mean(),'MAD':temp_df['WSAfter'].mad()})
        plt.xlim(0)
        plt.savefig(path_name+"/Figures/Histograms/{0}Coast/WS_{1}.png".format(coast,ty),dpi=500)
        plt.show()
        
        # #Timing of precipitation
        # bins = np.linspace(0, 50, 10)
        # x=temp_df['PrecipTiming']
        # plt.hist(x,bins,alpha=1)
        # plt.xlabel('Precipitation Timing [min]')
        # plt.ylabel('Frequency')
        # plt.title('Time of Precipitation \nAfter Sea Breeze Passage')
        # plt.text(10, 40, u'\u03bc=%(mu)2.2f min, \nMAD=%(MAD)2.2f min' %\
        #          {'mu':temp_df['PrecipTiming'].mean(),'MAD':temp_df['PrecipTiming'].mad()})
        # plt.xlim(0)
        # plt.show()
        """
        
        #Timing of SBF Passage
        bins = np.linspace(0, 24, 24)
        x=temp_df['Hour']
        # if i==0:
        #     rvs=x.copy()
        # else:
        #     cdf=x.copy()
        # n1,_,_=plt.hist(x,bins,alpha=1)
        # plt.xlabel('Time of SBF Passage [LST]')
        # plt.ylabel('Frequency')
        # plt.title('Timing SBC Type {0}\nFlorida {1} Coast'.format(ty,coast))
        # plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/3.0, u'\u03bc=%(mu)2.2f LST, \nMAD=%(MAD)2.2f LST' %\
        #          {'mu':temp_df['Hour'].mean(),'MAD':temp_df['Hour'].mad()})
        # plt.xlim(0)
        # # plt.savefig(path_name+"/Figures//Histograms/{0}Coast/SBF_Timing_{1}.png".format(coast,ty),dpi=500)
        # plt.show()
        
        print(west[i],x.mean())
        
        """
        #Relative Humidity before and after sea breeze passage
        bins = np.linspace(6, 31, 30)
        y=temp_df['DPTempBefore']
        y[y<-100]=np.nan
        y=y.dropna()
        x=temp_df['DPTempAfter'].dropna()
        x[x<-100]=np.nan
        x=x.dropna()
        n1,_,_=plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
        n2,_,_=plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
        plt.legend(loc='best')
        plt.xlabel('Dew Point Temperature [\N{DEGREE SIGN}C]')
        plt.ylabel('Frequency')
        plt.title('Dew Point Temperature SBC Type {0}\nFlorida {1} Coast'.format(ty,coast))
        plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/3.0, u'$\u03bc_b$=%(mu)2.2f, \n$MAD_b$=%(MAD)2.2f' %\
                 {'mu':y.mean(),'MAD':y.mad()})
        plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/6.0, u'$\u03bc_a$=%(mu)2.2f, \n$MAD_a$=%(MAD)2.2f' %\
                 {'mu':x.mean(),'MAD':x.mad()})
        plt.xlim(6)
        plt.savefig(path_name+"/Figures/Histograms/{0}Coast/DewPt_{1}.png".format(coast,ty),dpi=500)
        plt.show()
        
        #Temperature before and after sea breeze passage
        bins = np.linspace(0, 35, 35)
        y=temp_df['TBefore'].dropna()
        x=temp_df['TAfter'].dropna()
        n1,_,_=plt.hist(y,bins,alpha=0.5,label='Before Sea \nBreeze')
        n2,_,_=plt.hist(x,bins,alpha=0.5,label='After Sea \nBreeze')
        plt.legend(loc='best')
        plt.xlabel('Temperature [\N{DEGREE SIGN}C]')
        plt.ylabel('Frequency')
        plt.title('Temperature SBC Type {0}\nFlorida {1} Coast'.format(ty,coast))
        plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/3.0, u'$\u03bc_b$=%(mu)2.2f\N{DEGREE SIGN}C, \n$MAD_b$=%(MAD)2.2f\N{DEGREE SIGN}C' %\
                 {'mu':temp_df['TBefore'].mean(),'MAD':temp_df['TBefore'].mad()})
        plt.text((bins.max()-bins.min())/4.0, (n1.max()-n1.min())/6.0, u'$\u03bc_a$=%(mu)2.2f\N{DEGREE SIGN}C, \n$MAD_a$=%(MAD)2.2f\N{DEGREE SIGN}C' %\
                 {'mu':temp_df['TAfter'].mean(),'MAD':temp_df['TAfter'].mad()})
        plt.xlim(0)
        plt.savefig(path_name+"/Figures/Histograms/{0}Coast/Temp_{1}.png".format(coast,ty),dpi=500)
        plt.show()



#Frequency by month
freq_month = df["Date"].groupby(df["Date"].dt.month).count()

wid=0.25
plt.bar(freq_month.index,freq_month, color = 'r', width = wid)

plt.legend(loc='best')
plt.title("Sea Breeze Frequency by Month")
plt.xlabel("Month"); plt.ylabel("Frequency")

plt.savefig(path_name+"Figures/SBC_Frequency_Month.png", dpi=500)
plt.show()


#Frequency by year
df["Date"] = df["Date"].astype("datetime64")
freq_year = df["Date"].groupby(df["Date"].dt.year).count()

wid=0.25
plt.bar(freq_year.index,freq_year, color = 'r', width = wid)

plt.legend(loc='best')
plt.title("Sea Breeze Frequency by Year")
plt.xlabel("Year"); plt.ylabel("Frequency")

plt.savefig(path_name+"Figures/SBC_Frequency_Year.png", dpi=500)
plt.show()


#Rainfall Frequency both coasts
df_intense = df[df["TwelveHrPrecip"]>3]
intense_precip = df_intense["Date"].groupby(df_intense["Date"].dt.month).count()

df_less = df[((df["TwelveHrPrecip"]>0.05) & (df["TwelveHrPrecip"]<3.0))]
less_precip = df_less["Date"].groupby(df_less["Date"].dt.month).count()






###Plot more intense rainfall (>3in)
fig, ax = plt.subplots()

bar_width = 0.25

ax.bar(intense_precip.index,intense_precip,
                width = bar_width,
                color='r')

ax.set_xlabel('Month')
ax.set_ylabel('Frequency')
ax.set_title('Intense (>3in) 12Hr Rainfall')
ax.legend()

fig.tight_layout()
plt.savefig(path_name+"Figures/StationIntensePrecip.png", dpi=500)
plt.show()


###Plot less precipitation (less than 3 in)
fig, ax = plt.subplots()

bar_width = 0.25

ax.bar(less_precip.index,less_precip,
                width = bar_width,
                color='r')

ax.set_xlabel('Month')
ax.set_ylabel('Frequency')
ax.set_title('Less (<3in) 12Hr Rainfall')
ax.legend()

fig.tight_layout()
plt.savefig(path_name+"Figures/StationLessPrecip.png", dpi=500)
plt.show()






#Frequency by month by sbtype
classic_freq = classic["Date"].groupby(classic["Date"].dt.month).count()
classicDP_freq = classicDP["Date"].groupby(classicDP["Date"].dt.month).count()
classicWS_freq = classicWS["Date"].groupby(classicWS["Date"].dt.month).count()
weak_freq = weak["Date"].groupby(weak["Date"].dt.month).count()
DPWS_freq = DPWS["Date"].groupby(DPWS["Date"].dt.month).count()

fig, ax = plt.subplots()

wid=0.15
ax.bar(classic_freq.index-2.0*wid,classic_freq, color = 'b', width = wid, label="Classic")
ax.bar(classicDP_freq.index-1.0*wid,classicDP_freq, color = 'g', width = wid, label="Classic DP")
ax.bar(classicWS_freq.index+0.0*wid,classicWS_freq, color = 'r', width = wid, label="Classic WS")
ax.bar(weak_freq.index+1.0*wid,weak_freq, color = 'c', width = wid, label="Weak")
ax.bar(DPWS_freq.index+2.0*wid,DPWS_freq, color = 'm', width = wid, label="Dew Pt WS")

ax.legend(loc='best',bbox_to_anchor=(1, 0.8))
ax.set_title("Regime Frequency By Month By Type")
ax.set_xlabel("Month"); ax.set_ylabel("Frequency")

# fig.tight_layout()

plt.savefig(path_name+"Figures/SB_Type_Frequency_Month.png", dpi=500, bbox_inches='tight', pad_inches=0.5)
plt.show()

"""