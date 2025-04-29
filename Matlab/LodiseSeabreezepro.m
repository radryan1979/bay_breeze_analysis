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

%plot(sjtime,sjJwdir,'b',ltime,lJwdir,'g',cmtime,cmJwdir,'r')
title('Wind Direction on July 3rd')
xlabel('Time')
ylabel('Wind Direction ( 360 Degrees)')
Legend('Sjsn','Lewes','Cape may')

%plot(sjtime,sjJ_atmp,'b',ltime,lJ_atmp,'g',cmtime,cmJ_atmp,'r')
title('Air Temperatures on July 3rd')
xlabel('Time')
ylabel('Degrees Celcius')
Legend('Sjsn','Lewes','Cape may')

%plot(sjtime,sjJ_pres,'b',ltime,lJ_pres,'g',cmtime,cmJ_pres,'r')



yy = find(cmyear==2012);

for k = 1:enddata
    ww = find(cmwdir == 90)
end



