import csv


years = range(2008,2009)#2018)
buoy_dir = "/Volumes/LaCie/SeaBreeze/Florida/NDBC_Buoy/41010_GC_Control/"

for yr in years:
    txt_file = buoy_dir + "41010_{0}.txt".format(yr)
    csv_file = buoy_dir + "test.csv"
    
    
    in_txt = csv.reader(open(txt_file, "rb"), delimiter = '\t')
    out_csv = csv.writer(open(csv_file, 'wb'))
    
    out_csv.writerows(in_txt)