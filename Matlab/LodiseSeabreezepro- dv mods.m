%John Lodise
%Seabreeze Research

clc
clear

CM2012= importdata('E:\LodiseSeabreeze_Data\capemay4h2012edit.txt') %Location


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

aa = find(cmday==3 & cmmonth ==7);
cmtime = cmhour(aa) + (cmminute(aa)/60);
cmJwdir = cmwdir(aa);
cmJ_wspd = cmwspd(aa);
cmJ_atmp = cmatmp(aa);
cmJ_pres = cmpres(aa);

L2012 = importdata('E:\LodiseSeabreeze_Data\lewesd1h2012edit.txt');

lyear = L2012.data(:,1);
lmonth = L2012.data(:,2);
lday = L2012.data(:,3);
lhour = L2012.data(:,4);
lminute = L2012.data(:,5);
lwdir = L2012.data(:,6);

lwspd = L2012.data(:,7);
lwgst = L2012.data(:,8);
    
lpres = L2012.data(:,9);
latmp = L2012.data(:,10);
lwtmp = L2012.data(:,11);

aa = find(lday==3 & lmonth ==7);
ltime = lhour(aa) + (lminute(aa)/60);
lJwdir = lwdir(aa);
lJ_wspd = lwspd(aa);
lJ_atmp = latmp(aa);
lJ_pres = lpres(aa);

SJ2012 = importdata('E:\LodiseSeabreeze_Data\sjsn4h2012edit.txt');

sjyear = SJ2012.data(:,1);
sjmonth = SJ2012.data(:,2);
sjday = SJ2012.data(:,3);
sjhour = SJ2012.data(:,4);
sjminute = SJ2012.data(:,5);
sjwdir = SJ2012.data(:,6);

sjwspd = SJ2012.data(:,7);
sjwgst = SJ2012.data(:,8);
    
sjpres = SJ2012.data(:,9);
sjatmp = SJ2012.data(:,10);
sjwtmp = SJ2012.data(:,11);

aa = find(sjday==3 & sjmonth ==7);
sjtime = sjhour(aa) + (sjminute(aa)/60);
sjJwdir = sjwdir(aa);
sjJ_wspd = sjwspd(aa);
sjJ_atmp = sjatmp(aa);
sjJ_pres = sjpres(aa);

figure(1)
plot(sjtime,sjJwdir,'b',ltime,lJwdir,'g',cmtime,cmJwdir,'r')
hold on
title('Wind Direction on July 3rd')
xlabel('Time')
ylabel('Wind Direction ( 360 Degrees)')
Legend('Sjsn','Lewes','Cape may')

figure(2)
plot(sjtime,sjJ_atmp,'b',ltime,lJ_atmp,'g',cmtime,cmJ_atmp,'r')
hold on
title('Air Temperatures on July 3rd')
xlabel('Time')
ylabel('Degrees Celcius')
Legend('Sjsn','Lewes','Cape may')

%plot(sjtime,sjJ_pres,'b',ltime,lJ_pres,'g',cmtime,cmJ_pres,'r')

yy = find(cmyear==2012);

for k = 1:enddata
    ww = find(cmwdir == 90)
    end

% Good start - I've put another option below

zz=find(cmyear==2012);
for ll=1:length(zz)
    yy=find(cmwdir >= 180 & cmwdir < 360);
    cmbreeze=CM2012.data(yy,:);  %pulls out all data for these times
end

days=cmbreeze(:,3); %prints days in 2012 when this criteria was passed

% example routine to average data to a certain time range and then provide
% a criteria.

r=1;
s=1;
t=1;

for aa=1:12 %month loop
    month_index=find(cmbreeze(:,2)==aa) % finds data by month
    mean_mon_wind(r)=mean(cmbreeze(month_index,7); %creates monthly array
    mean_mon_dir(r)=mean(cmbreeze(month_index,8);
    day_count=[31 28 31 30 31 30 31 31 30 31 30 31]; %number of days per month
    mm=daycount(aa);
    r=r+1;
    for bb=1:mm
        day_index=find(cmbreeze(:,3)==bb & cmbreeze(:,2)==aa)
        daily_wind_spd(s)=mean(cmbreeze(day_index,7));
        daily_wind_dir(s)=mean(cmbreeze(day_index,8);
        s=s+1;
        for cc=1:24
            hour_index(t)=find(cmbreeze(:,4)==cc & cmbreeze(:,3)==bb & cmbreeze(:,2)==aa);
            hourly_wind_spd(t)=mean(cmbreeze(day_index,7));
            hourly_wind_dir(t)=mean(cmbreeze(day_index,8);
            t=t+1;
        end
    end
end
hourly_wind_dir=hourly_wind_dir+360;  %add 360 offset to avoid negative numbers
change_wind=diff(hourly_wind_dir);  %takes the difference of wind direction with hour prior
frontal_passge=find(change_wind>45); % looks for shifts greater than 45 deg per hour

% you could do the same thing for a temperature variable and combine the
% criteria. or create a wind direction criteria and a wind direction change
% criteria - you decide and we'll talk.
% 
% note that to make this really useful, you should save in the monthly, daily,
% and hourly files a time.  I suggest the best way is to  either use date2num or calculate a decimail julian day
%  something like jday=(sum of daycount from 1:aa + bb/mm + cc/24);
