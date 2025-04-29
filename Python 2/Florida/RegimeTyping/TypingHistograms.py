"""
Coded by: Dan Moore

The purpose of this code is to create a histogram of regimes.

Updated: 3-26-19
"""

import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt

path_name="/Volumes/LaCie/SeaBreeze/Florida/"
suffix="StationDetection/Results3_6_19/"
infile=open(os.path.join\
(path_name,suffix+"Hughes_SST_Regime.csv"),"r")

df=pd.read_csv(infile,sep=',',header=0,index_col=False)
# df=df[df["Station"]==station]


suffix2="TampaRegimeTyping/TampaTyping.csv"
infile=open(os.path.join\
(path_name,suffix2),"r")

type_df = pd.read_csv(infile,sep=',',header=0,index_col=0)
type_df = type_df[(type_df["Month"]>3) & (type_df["Month"]<11)]


#Histograms of Typing by coast
climo_df = type_df["Regime"].groupby(type_df["Regime"]).count()

for cst in ["East","West"]:
    if cst=="East":
        st=[340,410,420,440]
    else:
        st=[360,350,380,450,480,490]
    
    temp_df = df[np.isin(df["Station"],st)]

    detect_type_df = temp_df.drop_duplicates(subset="Date")
    detect_type_df = detect_type_df["Regime"].groupby(detect_type_df["Regime"]).count()
    
    
    plt.bar(left=detect_type_df.index,height=detect_type_df/climo_df*100.0)
    plt.xticks(range(1,10))
    plt.yticks(range(0,110,20))
    plt.title("Probability of SBC by Regime {0} Coast".format(cst))
    plt.ylabel("Probability of SBC")
    plt.xlabel("Regime")
    plt.savefig(path_name+suffix+"Figures/Prob_of_SBC_{0}.png".format(cst),dpi=500)
    plt.show()
    
    
    plt.bar(left=climo_df.index, height=climo_df.values)
    plt.xticks(range(1,10))
    # plt.yticks(range(0,110,20))
    plt.title("Frequency of Regime on {0} Coast".format(cst))
    plt.ylabel("Frequency")
    plt.xlabel("Regime")
    plt.savefig(path_name+suffix+"/Figures/Freq_of_Regimes{0}.png".format(cst),dpi=500)
    plt.show()


for index,row in type_df.iterrows():
    type_df.loc[index,"Datetime"]       =       str(int(row["Month"]))+"-"+str(int(row["Day"]))+"-"+str(int(row["Year"]))

type_df["Datetime"]             =       type_df["Datetime"].astype("datetime64")

regime1 = type_df[type_df["Regime"]==1.0]; regime2 = type_df[type_df["Regime"]==2.0]
regime3 = type_df[type_df["Regime"]==3.0]; regime4 = type_df[type_df["Regime"]==4.0]
regime5 = type_df[type_df["Regime"]==5.0]; regime6 = type_df[type_df["Regime"]==6.0]
regime7 = type_df[type_df["Regime"]==7.0]; regime8 = type_df[type_df["Regime"]==8.0]
regime9 = type_df[type_df["Regime"]==9.0]

reg_1_freq = regime1["Datetime"].groupby(regime1["Datetime"].dt.month).count()
reg_2_freq = regime2["Datetime"].groupby(regime2["Datetime"].dt.month).count()
reg_3_freq = regime3["Datetime"].groupby(regime3["Datetime"].dt.month).count()
reg_4_freq = regime4["Datetime"].groupby(regime4["Datetime"].dt.month).count()
reg_5_freq = regime5["Datetime"].groupby(regime5["Datetime"].dt.month).count()
reg_6_freq = regime6["Datetime"].groupby(regime6["Datetime"].dt.month).count()
reg_7_freq = regime7["Datetime"].groupby(regime7["Datetime"].dt.month).count()
reg_8_freq = regime8["Datetime"].groupby(regime8["Datetime"].dt.month).count()
reg_9_freq = regime9["Datetime"].groupby(regime9["Datetime"].dt.month).count()

fig, ax = plt.subplots()

wid=0.09
ax.bar(reg_1_freq.index-4.0*wid,reg_1_freq, color = 'b', width = wid, label="Reg 1")
ax.bar(reg_2_freq.index-3.0*wid,reg_2_freq, color = 'g', width = wid, label="Reg 2")
ax.bar(reg_3_freq.index-2.0*wid,reg_3_freq, color = 'r', width = wid, label="Reg 3")
ax.bar(reg_4_freq.index-1.0*wid,reg_4_freq, color = 'c', width = wid, label="Reg 4")
ax.bar(reg_5_freq.index+0.0*wid,reg_5_freq, color = 'm', width = wid, label="Reg 5")
ax.bar(reg_6_freq.index+1.0*wid,reg_6_freq, color = 'y', width = wid, label="Reg 6")
ax.bar(reg_7_freq.index+2.0*wid,reg_7_freq, color = 'k', width = wid, label="Reg 7")
ax.bar(reg_8_freq.index+3.0*wid,reg_8_freq, color = 'g', width = wid, label="Reg 8")
ax.bar(reg_9_freq.index+4.0*wid,reg_9_freq, color = 'r', width = wid, label="Reg 9")

ax.legend(loc='best',bbox_to_anchor=(1, 0.8))
ax.set_title("Regime Frequency By Month")
ax.set_xlabel("Month"); ax.set_ylabel("Frequency")

# fig.tight_layout()

plt.savefig(path_name+suffix+"Figures/Regime_Frequency_Month.png", dpi=500, bbox_inches='tight', pad_inches=0.5)
plt.show()


