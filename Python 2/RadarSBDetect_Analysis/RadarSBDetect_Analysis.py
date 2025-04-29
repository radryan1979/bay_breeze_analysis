import os
import pandas as pd
import numpy as np
import datetime as dt
import pyproj
import cartopy.crs as ccrs
import cartopy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm


import numpy as np
import cartopy.crs as ccrs
import cartopy.geodesic as cgeo

def blank_axes(ax):
        """
        blank_axes:  blank the extraneous spines and tick marks for an axes

        Input:
        ax:  a matplotlib Axes object

        Output: None
        """


        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.set_ticks_position('none')
        ax.tick_params(labelbottom='off', labeltop='off', labelleft='off', labelright='off' ,\
                        bottom='off', top='off', left='off', right='off' )
    #end blank_axes


def _axes_to_lonlat(ax, coords):
    """(lon, lat) from axes coordinates."""
    display = ax.transAxes.transform(coords)
    data = ax.transData.inverted().transform(display)
    lonlat = ccrs.PlateCarree().transform_point(*data, ax.projection)

    return lonlat


def _upper_bound(start, direction, distance, dist_func):
    """A point farther than distance from start, in the given direction.

    It doesn't matter which coordinate system start is given in, as long
    as dist_func takes points in that coordinate system.

    Args:
        start:     Starting point for the line.
        direction  Nonzero (2, 1)-shaped array, a direction vector.
        distance:  Positive distance to go past.
        dist_func: A two-argument function which returns distance.

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    if distance <= 0:
        raise ValueError(f"Minimum distance is not positive: {distance}")

    if np.linalg.norm(direction) == 0:
        raise ValueError("Direction vector must not be zero.")

    # Exponential search until the distance between start and end is
    # greater than the given limit.
    length = 0.1
    end = start + length * direction

    while dist_func(start, end) < distance:
        length *= 2
        end = start + length * direction

    return end


def _distance_along_line(start, end, distance, dist_func, tol):
    """Point at a distance from start on the segment  from start to end.

    It doesn't matter which coordinate system start is given in, as long
    as dist_func takes points in that coordinate system.

    Args:
        start:     Starting point for the line.
        end:       Outer bound on point's location.
        distance:  Positive distance to travel.
        dist_func: Two-argument function which returns distance.
        tol:       Relative error in distance to allow.

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    initial_distance = dist_func(start, end)
    if initial_distance < distance:
        raise ValueError(f"End is closer to start ({initial_distance}) than "
                         f"given distance ({distance}).")

    if tol <= 0:
        raise ValueError(f"Tolerance is not positive: {tol}")

    # Binary search for a point at the given distance.
    left = start
    right = end

    while not np.isclose(dist_func(start, right), distance, rtol=tol):
        midpoint = (left + right) / 2

        # If midpoint is too close, search in second half.
        if dist_func(start, midpoint) < distance:
            left = midpoint
        # Otherwise the midpoint is too far, so search in first half.
        else:
            right = midpoint

    return right


def _point_along_line(ax, start, distance, angle=0, tol=0.01):
    """Point at a given distance from start at a given angle.

    Args:
        ax:       CartoPy axes.
        start:    Starting point for the line in axes coordinates.
        distance: Positive physical distance to travel.
        angle:    Anti-clockwise angle for the bar, in radians. Default: 0
        tol:      Relative error in distance to allow. Default: 0.01

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    # Direction vector of the line in axes coordinates.
    direction = np.array([np.cos(angle), np.sin(angle)])

    geodesic = cgeo.Geodesic()

    # Physical distance between points.
    def dist_func(a_axes, b_axes):
        a_phys = _axes_to_lonlat(ax, a_axes)
        b_phys = _axes_to_lonlat(ax, b_axes)

        # Geodesic().inverse returns a NumPy MemoryView like [[distance,
        # start azimuth, end azimuth]].
        return geodesic.inverse(a_phys, b_phys).base[0, 0]

    end = _upper_bound(start, direction, distance, dist_func)

    return _distance_along_line(start, end, distance, dist_func, tol)


def scale_bar(ax, location, length, metres_per_unit=1000, unit_name='km',
              tol=0.01, angle=0, color='black', linewidth=3, text_offset=0.005,
              ha='center', va='bottom', plot_kwargs=None, text_kwargs=None,
              **kwargs):
    """Add a scale bar to CartoPy axes.

    For angles between 0 and 90 the text and line may be plotted at
    slightly different angles for unknown reasons. To work around this,
    override the 'rotation' keyword argument with text_kwargs.

    Args:
        ax:              CartoPy axes.
        location:        Position of left-side of bar in axes coordinates.
        length:          Geodesic length of the scale bar.
        metres_per_unit: Number of metres in the given unit. Default: 1000
        unit_name:       Name of the given unit. Default: 'km'
        tol:             Allowed relative error in length of bar. Default: 0.01
        angle:           Anti-clockwise rotation of the bar.
        color:           Color of the bar and text. Default: 'black'
        linewidth:       Same argument as for plot.
        text_offset:     Perpendicular offset for text in axes coordinates.
                         Default: 0.005
        ha:              Horizontal alignment. Default: 'center'
        va:              Vertical alignment. Default: 'bottom'
        **plot_kwargs:   Keyword arguments for plot, overridden by **kwargs.
        **text_kwargs:   Keyword arguments for text, overridden by **kwargs.
        **kwargs:        Keyword arguments for both plot and text.
    """
    # Setup kwargs, update plot_kwargs and text_kwargs.
    if plot_kwargs is None:
        plot_kwargs = {}
    if text_kwargs is None:
        text_kwargs = {}

    plot_kwargs = {'linewidth': linewidth, 'color': color, **plot_kwargs,
                   **kwargs}
    text_kwargs = {'ha': ha, 'va': va, 'rotation': angle, 'color': color,
                   **text_kwargs, **kwargs}

    # Convert all units and types.
    location = np.asarray(location)  # For vector addition.
    length_metres = length * metres_per_unit
    angle_rad = angle * np.pi / 180

    # End-point of bar.
    end = _point_along_line(ax, location, length_metres, angle=angle_rad,
                            tol=tol)

    # Coordinates are currently in axes coordinates, so use transAxes to
    # put into data coordinates. *zip(a, b) produces a list of x-coords,
    # then a list of y-coords.
    ax.plot(*zip(location, end), transform=ax.transAxes, **plot_kwargs)

    # Push text away from bar in the perpendicular direction.
    midpoint = (location + end) / 2
    offset = text_offset * np.array([-np.sin(angle_rad), np.cos(angle_rad)])
    text_location = midpoint + offset

    # 'rotation' keyword argument is in text_kwargs.
    ax.text(*text_location, f"{length} {unit_name}", rotation_mode='anchor',
            transform=ax.transAxes, **text_kwargs)



state='Delaware' 
# state='Florida'

directory = os.fsencode("/Volumes/LaCie/SeaBreeze/RadSBDetect5_1_19/{0}/".format(state))

success=0
cnt=0
june_cnt=0
july_cnt=0
aug_cnt=0
june_suc_cnt=0
july_suc_cnt=0
aug_suc_cnt=0


conv_cnt=0

first=True
first2=True

conv_dates = []

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    
    if filename.endswith("ScanLog.csv"):
        temp_df = pd.read_csv(os.path.join(directory,file).decode("utf-8"),index_col=0,header=0)
        if first:
            df=temp_df.copy()
            first=False
        else:
            df=pd.concat([df,temp_df],ignore_index=True)
    
    if filename.endswith("PrecipLog.csv"):
        temp_precip_df = pd.read_csv(os.path.join(directory,file).decode("utf-8"),index_col=0,header=0)
        if np.any(temp_precip_df["Mean_Reflectivity"]>30.0) and np.any(temp_precip_df["Max_Reflectivity"]>50.0):
            conv_cnt+=1
            date123=filename.split("P")[0].split("_")
            
            conv_dates=np.append(conv_dates,dt.date(int(date123[0]),int(date123[1]),int(date123[2])))
            
    
    if filename.endswith("_RadarDetect.csv"):
        tempdetect_df = pd.read_csv(os.path.join(directory,file).decode("utf-8"),index_col=0,header=0)
        if first:
            detect_df=tempdetect_df.copy()
            first2=False
        else:
            detect_df=pd.concat([detect_df,tempdetect_df],axis=1)
    
    if filename.endswith("RadarDetect.csv"): 
        if filename[5] == "6":
            june_suc_cnt+=1
        elif filename[5] == "7":
            july_suc_cnt+=1
        elif filename[5] == "8":
            aug_suc_cnt+=1
        # success+=1
    elif filename.endswith("DayLog.csv"): 
        t_df = pd.read_csv(os.path.join(directory,file).decode("utf-8"),index_col=0,header=0)
        if t_df["Error"][0]=="Sufficient data, but sea breeze not found.":
            if filename[5] == "6":
                june_cnt+=1
            elif filename[5] == "7":
                july_cnt+=1
            elif filename[5] == "8":
                aug_cnt+=1
            cnt+=1
        if t_df["Error"][0]=="Good day, sir.":
            if filename[5] == "6":
                june_cnt+=1
            elif filename[5] == "7":
                july_cnt+=1
            elif filename[5] == "8":
                aug_cnt+=1
            cnt+=1
    
    if filename.endswith("DayLog.csv"):
        log_df = pd.read_csv(os.path.join(directory,file).decode("utf-8"),index_col=0,header=0)
        if log_df["Sea_Breeze"].any():
            success+=1
    

df["Scan_Time"] = pd.to_datetime(df["Scan_Time"])
df=df.drop_duplicates(["Scan_Time"])
df_daily=df.groupby(df["Scan_Time"].dt.date).sum()
print("Days with SB detection (at least one scan): {0}".format(df_daily[df_daily["Sea_Breeze"]>0.0]["Sea_Breeze"].count()))
print("Days with Precipitation (at least one scan): {0}".format(df_daily[df_daily["Precipitation"]>0.0]["Precipitation"].count()))
    
print("June: {0}".format(june_cnt))
print("July: {0}".format(july_cnt))
print("August: {0}".format(aug_cnt))





if state == "Delaware":
    cstlnlat = np.array((39.25, 39.225, 39.2, 39.175, 39.15, 39.125, 39.1, 
                        39.075, 39.05, 39.025, 39.0, #Bay Coastline
                        38.975, 38.95, 38.925, 38.9, 38.875, 38.85, 38.825, 
                        38.8, 38.775,
                        38.75, 38.725, 38.7, 38.675, 38.65, 38.625, 38.6,#Ocean Coastline
                        38.575, 38.55, 38.525, 38.5, 38.475, 38.45, 38.425,
                        38.4, 38.375, 38.35, 38.325, 38.3, 38.275, 38.25)).T
                        
    cstlnlon = np.array((-75.42, -75.42, -75.41, -75.41, -75.41, -75.41, -75.40, 
                        -75.40, -75.39, -75.35, -75.33,#Bay Coastline
                        -75.32, -75.31, -75.31, -75.29, -75.26, -75.24, -75.21,
                        -75.18, -75.08,
                        -75.08, -75.08, -75.07, -75.07, -75.07, -75.06, -75.06,#Ocean Coastline
                        -75.06, -75.06, -75.05, -75.05, -75.05, -75.05, -75.05,
                        -75.05, -75.06, -75.07, -75.08, -75.09, -75.1, -75.11)).T

                    
elif state=="Florida":
    cstlnlat = np.arange(28.65,27.64,-0.025)
                        
    cstlnlon = np.array((-80.63, -80.61, -80.59, -80.57, -80.56, -80.56, -80.55, 
                        -80.53, -80.53, -80.57, -80.59,#Bay Coastline
                        -80.60, -80.60, -80.61, -80.61, -80.60, -80.60, -80.60,
                        -80.60, -80.59,
                        -80.58, -80.57, -80.57, -80.56, -80.55, -80.53, -80.52,#Ocean Coastline
                        -80.51, -80.50, -80.48, -80.47, -80.46, -80.44, -80.43,
                        -80.42, -80.40, -80.49, -80.48, -80.47, -80.46, -80.45))







#SELECT IF WE ARE LOOKING AT INTENSE PRECIP
#TRUE = INTENSE PRECIP
#FALSE = ALL DAYS
Intense_Precip=False
if Intense_Precip:
    out_dir = "/Volumes/LaCie/SeaBreeze/RadSBDetect5_1_19/ConvArticleFigures/"
else:
    out_dir = "/Volumes/LaCie/SeaBreeze/RadSBDetect5_1_19/ArticleFigures/"





#Handle coordinate dataframe:
use_detect_df=detect_df.copy()
use_detect_df=use_detect_df.transpose()
use_detect_df.index = pd.to_datetime(use_detect_df.index)
success_dates = df_daily.index[df_daily["Sea_Breeze"]>1.0]
success_dates = np.array(success_dates)


if Intense_Precip:
    success_dates = success_dates[np.isin(success_dates,conv_dates)]
    
    
num_detected_scans = np.mean(df_daily["Sea_Breeze"][df_daily["Sea_Breeze"]>0.0])


j=0
for i in range(len(success_dates)):
    if success_dates[i].month>7:
        j+=1
onset_df = df.copy()
onset_df = onset_df[onset_df["Sea_Breeze"]]
onset_df = onset_df.reset_index(drop=True)
onset_df["Date"] = onset_df["Scan_Time"].dt.date

use_detect_df = use_detect_df[np.isin(np.array(use_detect_df.index.date),success_dates)]
cnt_by_lat = use_detect_df.groupby(use_detect_df.index.date).count()

first_coords = []
last_coords = []

for i in range(len(cstlnlat)):
    lat_first = []
    lat_last = []
    for j in range(len(cnt_by_lat)):
        if cnt_by_lat[round(cstlnlat[i],3)][j]>1:
            sub_use_detect_df = use_detect_df[np.isin(np.array(use_detect_df.index.date),cnt_by_lat.index[j])]
            first_detect_df = sub_use_detect_df[round(cstlnlat[i],3)].groupby(sub_use_detect_df.index.date).first()
            last_detect_df = sub_use_detect_df[round(cstlnlat[i],3)].groupby(sub_use_detect_df.index.date).last()
            lat_first = np.append(lat_first,first_detect_df[0])
            lat_last = np.append(lat_last,last_detect_df[0])
            
    first_coords = np.append(first_coords,np.median(lat_first))
    last_coords = np.append(last_coords,np.median(lat_last))
            









# onset_df = onset_df.groupby(onset_df["Scan_Time"].dt.date)
# onset_df = onset_df[np.isin(onset_df.index,success_dates)]
# onset_times = onset_df["Scan_Time"]
# onset_times = pd.DataFrame({"Date": onset_times.index,
#                             "Scan_Time": onset_times.values
#                             })
# onset_times = onset_times.set_index("Scan_Time")
# use_detect_df.index = pd.to_datetime(use_detect_df.index)
# use_detect_df = use_detect_df[np.isin(use_detect_df.index,onset_times.index)]
# new_df = pd.merge(left=onset_times,right=use_detect_df, left_index=True, left_on=None, right_index=True)
# first_coords = new_df.mean(axis=0, skipna=True)#####coords of average onset


# offset_df = df.copy()
# offset_df = offset_df[offset_df["Sea_Breeze"]]
# offset_df = offset_df.groupby(offset_df["Scan_Time"].dt.date).last()
# offset_df = offset_df[np.isin(offset_df.index,success_dates)]
# offset_times = offset_df["Scan_Time"]
# offset_times = pd.DataFrame({"Date": offset_times.index,
#                             "Scan_Time": offset_times.values
#                             })
# offset_times = offset_times.set_index("Scan_Time")

# use_detect_df=detect_df.copy()
# use_detect_df=use_detect_df.transpose()
# use_detect_df.index = pd.to_datetime(use_detect_df.index)
# use_detect_df = use_detect_df[np.isin(use_detect_df.index,offset_times.index)]

# new2_df = pd.merge(left=offset_times,right=use_detect_df, left_index=True, left_on=None, right_index=True)
# last_coords = new2_df.mean(axis=0, skipna=True)


#Analyzing distance from coast:
#Delaware:
if state == "Delaware":
    #Calculate in km (different distance between lon lines)
    dist_to_cst_first = (cstlnlon - first_coords)*87.0
    dist_to_cst_last = (cstlnlon - last_coords)*87.0
    spd_by_trans = (first_coords-last_coords)*87/3.6*1000/3600.0
    
    bay_dist_first = dist_to_cst_first[:19].mean()
    ocn_dist_first = dist_to_cst_first[19:].mean()
    bay_dist_last = dist_to_cst_last[:19].mean()
    ocn_dist_last = dist_to_cst_last[19:].mean()
    
                    
elif state=="Florida":
    #Calculate in km (different distance between lon lines)
    dist_to_cst_first = (cstlnlon - first_coords)*98.0
    dist_to_cst_last = (cstlnlon - last_coords)*98.0
    spd_by_trans = (first_coords-last_coords)*98.0/3.6*1000/3600.0
    ave_dist_first = dist_to_cst_first.mean()
    ave_dist_last = dist_to_cst_last.mean()


if Intense_Precip:
    spd_by_trans_int = spd_by_trans
else:
    spd_by_trans_non = spd_by_trans


#PLOT speed by transect
fig, ax = plt.subplots()
plt.plot(spd_by_trans_non,cstlnlat,  'rp', markersize=4, label='All Days')
plt.plot(spd_by_trans_int,cstlnlat,  'bp', markersize=4, label='Intense Precip')
plt.title("{0} Sea Breeze Velocity".format(state))
plt.xlabel("Velocity (m s$^-$$^1$)"); plt.ylabel("Latitude")
ax.set(aspect=1,
        xlim=(spd_by_trans_int.min()-0.1,spd_by_trans_int.max()+0.1),
        ylim=(cstlnlat[-1]-0.1,cstlnlat[0]+0.1))
plt.legend(loc='best')
plt.savefig(out_dir+'VelocityByTransect{0}.png'.format(state), dpi=500)
plt.show()









#Average times:
offset_times["Hour"]=offset_times.index.hour + offset_times.index.minute/60.0
onset_times["Hour"]=onset_times.index.hour + onset_times.index.minute/60.0

#Duration by month:
mo_duration=0;j=0
for i in range(len(onset_times)):
    if onset_times.index[i].month==6:
        mo_duration += offset_times.iloc[i]["Hour"]-onset_times.iloc[i]["Hour"]
        j+=1
mo_duration = mo_duration/j

duration = offset_times["Hour"]-onset_times["Hour"]
ave_dur = np.mean(duration)
std_dur = np.std(duration)




#Analyzing precip per event:
#SHOULD LOOK AT IF THERE IS PRECIP, WHERE IS IT? 
#AVE REF, MAX REF, WHERE IS IT IN RELATION TO SBF
tdf = df.copy()
tdf = tdf[np.isin(tdf["Scan_Time"].dt.date,success_dates)]
tdf = tdf.reset_index()
event = tdf["Scan_Time"][0].day
look4precip=False
skip2nxtday=False
event_precip=0
event_conv_precip=0
cell_cnt=0

size=[]
num_conv_cell=[]

event_mean_ref = []
overall_mean_ref = []
hr_aft_sb = []
hrly_size=[]
event_size=0


j=0
for i in range(1,len(tdf)):
    if tdf["Scan_Time"][i].day == event:
        if skip2nxtday:
            j+=1
            continue
            
        
    else:    
        skip2nxtday=False
        look4precip=False
                
        event = tdf["Scan_Time"][i].day
        
    if tdf["Sea_Breeze"][i]:
        look4precip=True
        time_of_p = tdf["Scan_Time"][i]
        y=time_of_p.year;m=time_of_p.month;d=time_of_p.day;h=time_of_p.hour
        file_str = "{0}_{1}_{2}PrecipLog.csv".format(y,m,d)

        prec_df = pd.read_csv(directory.decode("utf-8")+file_str,index_col=0,header=0)
        prec_df = prec_df.reset_index()
        prec_df["Scan_Time"]=pd.to_datetime(prec_df["Scan_Time"])
        
        first_precip=True
        
        for j in range(len(prec_df)):
            
            if prec_df["Scan_Time"][j].hour < h:
                continue
            else:
                if prec_df["Max_Reflectivity"][j]>50.0 and prec_df["Mean_Reflectivity"][j]>30.0:
                    
                    
                    
                    if first_precip:
                        precip_hr = prec_df["Scan_Time"][j].hour - h
                        first_precip=False
                        
                    
                    if precip_hr!=prec_df["Scan_Time"][j].hour - h:
                        print(precip_hr,event_size)
                        hr_aft_sb = np.append(hr_aft_sb,precip_hr)
                        overall_mean_ref = np.append(overall_mean_ref, sum(event_mean_ref)/float(event_size))
                        hrly_size = np.append(hrly_size,event_size)
                        event_size=0
                        event_mean_ref=[]
                        
                        precip_hr = prec_df["Scan_Time"][j].hour - h
                        
                    
                    event_size+=prec_df["Pixel_Area"][j]
                    event_mean_ref = np.append(event_mean_ref,prec_df["Mean_Reflectivity"][j]*prec_df["Pixel_Area"][j])
                    

                    
                    size=np.append(size,prec_df["Pixel_Area"][j])
                    
                    
                    cell_cnt+=1
                #     if skip2nxtday:
                #         continue
                #     else:
                #         event_conv_precip+=1
                #         # skip2nxtday=True
                #     # break
                # elif prec_df["Max_Reflectivity"][j]>0.0:
                #     if skip2nxtday:
                #         continue
                #     else:
                #         event_precip+=1
                    # skip2nxtday=True
                    # break
        
        if not first_precip:    
            hr_aft_sb = np.append(hr_aft_sb,precip_hr)
            overall_mean_ref = np.append(overall_mean_ref, sum(event_mean_ref)/float(event_size))
            hrly_size = np.append(hrly_size,event_size)
            event_size=0
            event_mean_ref=[]
                    
        skip2nxtday=True
        if cell_cnt>0:
            num_conv_cell=np.append(num_conv_cell,cell_cnt)
            cell_cnt=0

size_pixel_FL = 0.044 #km^2
size_pixel_DE = 0.0374 #km^2

a=300.0
b=1.4

if state=="Delaware":
    mean_ref_DE = pd.DataFrame({"Hour_aft_sb": hr_aft_sb,
                                "Norm_Mean_Ref": overall_mean_ref,
                                "Pixel_Size": hrly_size,
                                "Area": hrly_size*size_pixel_DE
                                })
    norm_ref_DE = mean_ref_DE["Norm_Mean_Ref"].groupby(mean_ref_DE["Hour_aft_sb"]).mean()
    hrly_size_DE = mean_ref_DE["Area"].groupby(mean_ref_DE["Hour_aft_sb"]).mean()
elif state=="Florida":
    mean_ref_FL = pd.DataFrame({"Hour_aft_sb": hr_aft_sb,
                                "Norm_Mean_Ref": overall_mean_ref,
                                "Pixel_Size": hrly_size,
                                "Area": hrly_size*size_pixel_FL
                                })
    norm_ref_FL = mean_ref_FL["Norm_Mean_Ref"].groupby(mean_ref_FL["Hour_aft_sb"]).mean()
    hrly_size_FL = mean_ref_FL["Area"].groupby(mean_ref_FL["Hour_aft_sb"]).mean()

norm_rainrate_DE = np.empty_like(norm_ref_DE)
for i in range(len(norm_ref_DE)):
    norm_rainrate_DE[i]=(10**(norm_ref_DE[i]/10.0)/200)**(5.0/8.0)
    
norm_rainrate_FL = np.empty_like(norm_ref_FL)
for i in range(len(norm_ref_FL)):
    norm_rainrate_FL[i]=(10**(norm_ref_FL[i]/10.0)/200)**(5.0/8.0)


fig, ax = plt.subplots()
wid=0.25
plt.bar(norm_ref_DE.index-0.5*wid,norm_rainrate_DE, color = 'b', width = wid,label="Delaware")
plt.bar(norm_ref_FL.index+0.5*wid,norm_rainrate_FL, color = 'g', width = wid,label="Florida")
plt.legend(loc='best')
# plt.yticks((3,3.5,4,4.5,5))
plt.xticks((0,1,2,3,4,5,6,7,8,9,10))
ax.set_ylim(2,5)
plt.title("Rainfall Intensity After SB Detection")
plt.xlabel("Hour After Sea Breeze Detection"); plt.ylabel("Normalized Rainfall Rate (mm $hr^-$$^1$)")
plt.savefig(out_dir+'RainfallIntensity.png', dpi=500)
plt.show()

fig, ax = plt.subplots()
wid=0.25
plt.bar(hrly_size_DE.index-0.5*wid,hrly_size_DE.values, color = 'b', width = wid,label="Delaware")
plt.bar(hrly_size_FL.index+0.5*wid,hrly_size_FL.values, color = 'g', width = wid,label="Florida")
plt.legend(loc='best')
# plt.yticks((30,31,32,33,34,35))
plt.xticks((0,1,2,3,4,5,6,7,8,9,10))
ax.set_ylim(0,2000)
plt.title("Intense Precipitation Coverage")
plt.xlabel("Hour After Sea Breeze Detection"); plt.ylabel("Area ($km^2$)")
plt.savefig(out_dir+'RainfallCoverage.png', dpi=500)
plt.show()


"""
Plots
"""
#Frequency by month by month/year
df_temp=df_daily[df_daily["Sea_Breeze"]>1.0]
df_temp.index = pd.to_datetime(df_temp.index)
mo_freq = df_temp["Sea_Breeze"].groupby(df_temp.index.month).count()
yr_freq = df_temp["Sea_Breeze"].groupby(df_temp.index.year).count()


fig, ax = plt.subplots()
wid=0.25
plt.bar(mo_freq.index,mo_freq, color = 'r', width = wid)
plt.legend(loc='best')
plt.title("{0} Sea Breeze Frequency by Month".format(state))
plt.xlabel("Month"); plt.ylabel("Frequency")
# plt.savefig(out_dir+'DetectionByMonth{0}.png'.format(state), dpi=500)
plt.show()

fig, ax = plt.subplots()
wid=0.5
plt.bar(yr_freq.index,yr_freq, color = 'r', width = wid)
plt.legend(loc='best')
plt.title("{0} Sea Breeze Frequency by Year".format(state))
plt.xlabel("Year"); plt.ylabel("Frequency")
# plt.savefig(out_dir+'DetectionByYear{0}.png'.format(state), dpi=500)
plt.show()



##SB detections by # of successful scans required for SB
df_one=df_daily["Sea_Breeze"][df_daily["Sea_Breeze"]>0.0].count()
df_two=df_daily["Sea_Breeze"][df_daily["Sea_Breeze"]>1.0].count()
df_three=df_daily["Sea_Breeze"][df_daily["Sea_Breeze"]>2.0].count()
df_four=df_daily["Sea_Breeze"][df_daily["Sea_Breeze"]>3.0].count()
df_five=df_daily["Sea_Breeze"][df_daily["Sea_Breeze"]>4.0].count()
df_six=df_daily["Sea_Breeze"][df_daily["Sea_Breeze"]>5.0].count()
df_seven=df_daily["Sea_Breeze"][df_daily["Sea_Breeze"]>6.0].count()

values = [df_one,df_two,df_three,df_four,df_five,df_six,df_seven]
nums = [1,2,3,4,5,6,7]

fig, ax = plt.subplots()
plt.plot(nums, values, '-o')
plt.title("{0} Sea Breeze Detection".format(state))
plt.xlabel("Required Number of Scans"); plt.ylabel("Successful Detections")
# plt.savefig(out_dir+'NumberOfScans{0}.png'.format(state), dpi=500)
plt.show()



##Time of first onset
onset_df = df.copy()
onset_df = onset_df[onset_df["Sea_Breeze"]]
onset_df = onset_df.groupby(onset_df["Scan_Time"].dt.date).first()
onset_df["Scan_Time"] = pd.to_datetime(onset_df["Scan_Time"])
onset_df["onset_time"] = onset_df['Scan_Time'].dt.hour + onset_df['Scan_Time'].dt.minute/60.0
onset_df["onset_time"] = onset_df["onset_time"].round(0)
onset_freq = onset_df["Sea_Breeze"].groupby([onset_df["Scan_Time"].dt.month,onset_df["onset_time"]]).count()

fig, ax = plt.subplots()
plt.plot(onset_freq[6].index, onset_freq[6], '-o',label="June")
plt.plot(onset_freq[7].index, onset_freq[7], '-*',label="July")
plt.plot(onset_freq[8].index, onset_freq[8], '-^',label="August")
plt.legend(loc='best')
plt.title("{0} Sea Breeze Onset".format(state))
plt.xlabel("Time of First Detection (UTC)"); plt.ylabel("Frequency")
# plt.savefig(out_dir+'OnsetTimeByMonth{0}.png'.format(state), dpi=500)
plt.show()


##Time of last detection
offset_df = df.copy()
offset_df = offset_df[offset_df["Sea_Breeze"]]
offset_df = offset_df.groupby(offset_df["Scan_Time"].dt.date).last()
offset_df["Scan_Time"] = pd.to_datetime(offset_df["Scan_Time"])
offset_df["onset_time"] = offset_df['Scan_Time'].dt.hour + offset_df['Scan_Time'].dt.minute/60.0
offset_df["onset_time"] = offset_df["onset_time"].round(0)
onset_freq = offset_df["Sea_Breeze"].groupby([offset_df["Scan_Time"].dt.month,offset_df["onset_time"]]).count()

fig, ax = plt.subplots()
plt.plot(onset_freq[6].index, onset_freq[6], '-o',label="June")
plt.plot(onset_freq[7].index, onset_freq[7], '-*',label="July")
plt.plot(onset_freq[8].index, onset_freq[8], '-^',label="August")
plt.legend(loc='best')
plt.title("{0} Sea Breeze Final Detection".format(state))
plt.xlabel("Time of Final Detection (UTC)"); plt.ylabel("Frequency")
# plt.savefig(out_dir+'OffsetTimeByMonth{0}.png'.format(state), dpi=500)
plt.show()



##PLOT First and last coordinates:
if state=="Delaware":
    xlimW=-76.45 ;xlimE=-74.6
    ylimS=38.15 ;ylimN=39.65
elif state=="Florida":
    xlimW=-81.725;xlimE=-79.875
    ylimS=27.375 ;ylimN=28.875


fig = plt.figure(figsize=[8, 8], dpi=500)
my_ax = plt.axes(projection = ccrs.PlateCarree())

states = cartopy.feature.NaturalEarthFeature(category='cultural',
                              name='admin_1_states_provinces_lines',
                              scale='10m', facecolor='none')
coast = cartopy.feature.NaturalEarthFeature(category='physical', scale='10m',
                            facecolor='none', name='coastline')
                            
reader = shpreader.Reader('/Volumes/LaCie/DPM_Code/Python/countyl010g_shp_nt00964/countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())

my_ax.add_feature(COUNTIES, facecolor='none', edgecolor='lightgray')                            
my_ax.add_feature(states, linestyle='-', edgecolor='darkslategray',linewidth=1)
my_ax.add_feature(coast, linestyle='-', edgecolor='darkslategray',linewidth=1)

plt.title('{0} NEXRAD'.format(state))
my_ax.set(aspect=1,
        xlim=(xlimW,xlimE),
        ylim=(ylimS,ylimN))
        
plt.plot(first_coords,cstlnlat,'rp',markersize=4, label="Average First Detection")
plt.plot(last_coords,cstlnlat,'bp',markersize=4, label="Average Final Detection")
if state=="Florida":
    plt.plot(-80.6356,28.0997,'gh',markersize=8,label="KMLB NEXRAD")
elif state=="Delaware":
    plt.plot(-75.435,38.821,'gh',markersize=8,label="KDOX NEXRAD")
plt.legend(loc='best')
scale_bar(my_ax, (0.1, 0.05), 50)

left = 0.72
bottom = 0.25
width = 0.16
height = 0.2
rect = [left,bottom,width,height]
rect = [left,bottom,width,height]
ax4 = plt.axes(rect)

# need a font that support enough Unicode to draw up arrow. need space after Unicode to allow wide char to be drawm?
ax4.text(0.5, 0.0,u'\u25B2 \nN ', ha='center', fontsize=30, family='Arial', rotation = 0)
blank_axes(ax4)

plt.savefig(out_dir+"{0}AveCoordinateMap.png".format(state), dpi=500)
plt.show()
