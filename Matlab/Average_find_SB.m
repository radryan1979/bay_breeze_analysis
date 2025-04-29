% Average_find_SB.m

% modified from original by D. Veron - 25 April 2013

% This code takes the data, averages the data to hourly values, stores the
% data in an hourly array and then apply seabreeze detection criteria.

CM2012= importdata('E:\LodiseSeabreeze_Data\capemay4h2012edit.txt'); %Location

% read in data from file and assign to various arrays

startdata = 1;
enddata = size(CM2012.data,1);

cmyear = CM2012.data(:,1);
cmmonth = CM2012.data(:,2);
cmday = CM2012.data(:,3);
cmhour = CM2012.data(:,4);
cmminute = CM2012.data(:,5);
cmwdir = CM2012.data(:,6);

cmwspd = CM2012.data(:,7);
cmwgst = CM2012.data(:,8);
    
cmpres = CM2012.data(:,9);
cmatmp = CM2012.data(:,10);
cmwtmp = CM2012.data(:,11);

aa = find(cmday==day & cmmonth ==month);
cmtime = cmday(aa)+((cmhour(aa)-1)/24) + (cmminute(aa)/1440)
cmJwdir = cmwdir(aa);
NN = find(cmJwdir == 999.0);
cmJwdir(NN) = NaN;
cmJ_wspd = cmwspd(aa);
cmJ_atmp = cmatmp(aa);
cmJ_pres = cmpres(aa);

% apply breeze criteria before averaging

zz=find(cmyear==2012);
for ll=1:length(zz)
   yy=find(cmwdir < 180 & cmwdir > 90);  % criteria for bay breeze
    cmbreeze=CM2012.data(yy,:);  %pulls out all data for these times
end

days=cmbreeze(:,3) %prints days in 2012 when this criteria was passed

% example routine to average data to a certain time range and then provide
% a criteria.

r=1;
s=1;
t=1;

for dd=1:12 %month loop
   month_index=find(CM2012.data(:,2)==dd) % finds data by month
    mean_mon_spd=mean(CM2012.data(month_index,7)); %creates monthly array
    mean_mon_dir=mean(CM2012.data(month_index,6));
    mean_year=mean(CM2012.data(month_index,1));
    mon_mean_array(r)=[mean_year dd mean_mon_spd mean_mon_dir];
    %if this doesn't work you need to add back in (r) indices for
    %everything
    day_count=[31 28 31 30 31 30 31 31 30 31 30 31]; %number of days per month
    mm=day_count(dd);
    r=r+1;
    for bb=1:mm
        day_index=find(CM2012.data(:,3)==bb & CM2012.data(:,2)==dd)
        daily_wind_spd=mean(CM2012.data(day_index,7));
        daily_wind_dir=mean(CM2012.data(day_index,6));
        day_mean_array(s)=[mean_year dd bb daily_wind_spd daily_wind_dir];
        s=s+1;
        for cc=1:24
            hour_index=find(CM2012.data(:,4)==cc & CM2012.data(:,3)==bb & CM2012.data(:,2)==dd);
            hourly_wind_spd=mean(CM2012.data(day_index,7));
            hourly_wind_dir=mean(CM2012.data(day_index,6));
            hourly_temp=mean(CM2012.data(day_index,10));
            hourly_mean_array(t)=[mean_year dd bb cc hourly_wind_spd hourly_wind_dir hourly_temp];
            t=t+1;
        end
    end
end

% save as .mat file
save mon_mean_array day_mean_array hourly_mean_array means.mat

% save as text files
save cm2012_m.txt mon_mean_array -ascii
save cm2012_d.txt day_mean_array -ascii
save cm2012_h.txt hourly_mean_array -ascii

% Now can apply criteria to averaged data.  This might be best for a change
% of wind direction.

%hourly_wind_dir=hourly_wind_dir+360;  %add 360 offset to avoid negative numbers
%change_wind=diff(hourly_wind_dir);  %takes the difference of wind direction with hour prior
%frontal_passage=find(change_wind>45); % looks for shifts greater than 45 deg per hour
%change_temp = diff(hourly_temp);
%temp_drop=find(change_temp>2); 
%SB = intersect(temp_drop,frontal_passage); 



% SB_count(1:12,1:31,1:24) = 0

% you could do the same thing for a temperature variable and combine the
% criteria. or create a wind direction criteria and a wind direction change
% criteria - you decide and we'll talk.
% 
% note that to make this really useful, you should save in the monthly, daily,
% and hourly files a time.  I suggest the best way is to  either use date2num or calculate a decimail julian day
%  something like jday=(sum of daycount from 1:aa + bb/mm + cc/24);

% clean up memory

clear
pack

