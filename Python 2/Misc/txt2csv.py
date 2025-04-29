import pandas as pd
import glob, os

st=41113
 
os.chdir("/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/{0}_Data/".format(st))
results = pd.DataFrame([])
 
for counter, file in enumerate(glob.glob("{0}*.txt".format(st))):
    namedf = pd.read_table(file, delimiter="\s+", header=[0,1])#, skiprows=[1])
    results = results.append(namedf)
 
results.to_csv('/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/{0}_Data/{0}_08_18_COMBINED.csv'.format(st))