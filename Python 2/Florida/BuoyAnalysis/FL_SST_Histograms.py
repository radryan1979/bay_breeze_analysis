import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt

BUOYdirpath = "/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/"

westSSTpath = BUOYdirpath + '42013_WestFL/42013COMBINED.csv'
eastSSTpath = BUOYdirpath + '41009_EastFL/41009COMBINED.csv'
southSST_insidepath = BUOYdirpath + 'PKYF1_SouthFL/PKYF1COMBINED.csv'
southSST_outsidepath = BUOYdirpath + 'MLRF1_SouthFL/MLRF1COMBINED.csv'

westSST_df = pd.read_csv(westSSTpath,header=0, skiprows=[1],index_col=0)
westSST_df = westSST_df.reset_index(drop=True)
westSST_df = westSST_df.replace(999,np.NaN)
eastSST_df = pd.read_csv(eastSSTpath,header=0,skiprows=[1],index_col=0)
eastSST_df = eastSST_df.reset_index(drop=True)
eastSST_df = eastSST_df.replace(999,np.NaN)
southSSTinside_df = pd.read_csv(southSST_insidepath,header=0,skiprows=[1],index_col=0)
southSSTinside_df = southSSTinside_df.reset_index(drop=True)
southSSTinside_df = southSSTinside_df.replace(999,np.NaN)
southSSToutside_df = pd.read_csv(southSST_outsidepath,header=0,skiprows=[1],index_col=0)
southSSToutside_df = southSSToutside_df.reset_index(drop=True)
southSSToutside_df = southSSToutside_df.replace(999,np.NaN)

#SBC SST by month
for i in range(4,11):
    bins = np.linspace(18, 33, 30)
    x=westSST_df.loc[westSST_df["MM"]==i]["WTMP"].dropna()
    plt.hist(x,bins,alpha=1)
    plt.xlabel('SST [\N{DEGREE SIGN}C]')
    plt.ylabel('Frequency')
    plt.title('Florida West Coast \nSST Month: {0}'.format(i))
    plt.text(20, 20, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
             {'mu':westSST_df.loc[westSST_df["MM"]==i]["WTMP"].dropna().mean(),\
             'MAD':westSST_df.loc[westSST_df["MM"]==i]["WTMP"].dropna().mad()})
    # plt.xlim(0)
    plt.savefig(path_name+"/Histograms/SST_Climatology/SSTWest_{0}.png".format(i),dpi=500)
    plt.show()


for i in range(4,11):
    bins = np.linspace(20, 32, 30)
    x=eastSST_df.loc[eastSST_df["MM"]==i]["WTMP"].dropna()
    plt.hist(x,bins,alpha=1)
    plt.xlabel('SST [\N{DEGREE SIGN}C]')
    plt.ylabel('Frequency')
    plt.title('Florida East Coast \nSST Month: {0}'.format(i))
    plt.text(20, 20, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
             {'mu':eastSST_df.loc[eastSST_df["MM"]==i]["WTMP"].dropna().mean(),\
             'MAD':eastSST_df.loc[eastSST_df["MM"]==i]["WTMP"].dropna().mad()})
    # plt.xlim(0)
    plt.savefig(path_name+"/Histograms/SST_Climatology/SSTEast_{0}.png".format(i),dpi=500)
    plt.show()

for i in range(4,11):
    bins = np.linspace(20, 36, 30)
    x=southSSTinside_df.loc[southSSTinside_df["MM"]==i]["WTMP"].dropna()
    plt.hist(x,bins,alpha=1)
    plt.xlabel('SST [\N{DEGREE SIGN}C]')
    plt.ylabel('Frequency')
    plt.title('Florida South Coast \nInside Keys SST Month: {0}'.format(i))
    plt.text(20, 20, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
             {'mu':southSSTinside_df.loc[southSSTinside_df["MM"]==i]["WTMP"].dropna().mean(),\
             'MAD':southSSTinside_df.loc[southSSTinside_df["MM"]==i]["WTMP"].dropna().mad()})
    # plt.xlim(0)
    plt.savefig(path_name+"/Histograms/SST_Climatology/SSTSouthInsideKeys_{0}.png".format(i),dpi=500)
    plt.show()
    
for i in range(4,11):
    bins = np.linspace(20, 32, 30)
    x=southSSToutside_df.loc[southSSToutside_df["MM"]==i]["WTMP"].dropna()
    plt.hist(x,bins,alpha=1)
    plt.xlabel('SST [\N{DEGREE SIGN}C]')
    plt.ylabel('Frequency')
    plt.title('Florida South Coast \nOutside Keys SST Month: {0}'.format(i))
    plt.text(20, 20, u'\u03bc=%(mu)2.2f\N{DEGREE SIGN}C, \nMAD=%(MAD)2.2f\N{DEGREE SIGN}C' %\
             {'mu':southSSToutside_df.loc[southSSToutside_df["MM"]==i]["WTMP"].dropna().mean(),\
             'MAD':southSSToutside_df.loc[southSSToutside_df["MM"]==i]["WTMP"].dropna().mad()})
    # plt.xlim(0)
    plt.savefig(path_name+"/Histograms/SST_Climatology/SSTSouthOutsideKeys_{0}.png".format(i),dpi=500)
    plt.show()
