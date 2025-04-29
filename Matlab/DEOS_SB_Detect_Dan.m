%Coded by Dan Moore
%Edited 7/1/19
%The purpose of this code is to take in Station data and detect Sea Breeze
%presence at each  'test' stations, per Gilchrist, 2013. The following 
%stations will be used as test stations:
%DBBB - Bethany Beach, DE - Boardwalk
%DBNG - Bethany Beach, DE
%DSJR - Kitts Hummock, DE
%DWAR - Harbeson, DE



clc
clear
close all


%Read in folder and included files.
cd C:\Users\dpm\Documents\MATLAB\RadSBDEOSData\
files=dir('*.csv');

teststationlist = ['DBBB'; 'DBNG'; 'DSJR'; 'DWAR'];

%Define onshore direction for each station; they correlate with each test 
%station in list above.
onshore = [89.15; 87.22; 84.41; 82.06];

j=1;
seabreezetime(1,length(teststationlist))=datetime;
seabreezestn=char([]);
seabreezedate(1,length(teststationlist))=datetime;

for file = files'
    f=file.name;
    stationnum=f(1:4);
    
    M                    =      readtable(f);
    station(j).timeEST   =      M(:,1);
    station(j).temp       =     M(:,2);
    station(j).num       =      stationnum;
%     station(j).temp      =      M(:,4);
%     station(j).relhum    =      M(:,6);
    station(j).dewtemp   =      M(:,3);
    station(j).windspd   =      M(:,4);
    station(j).winddir   =      M(:,6);
    station(j).precip    =      M(:,8);
%     station(j).solrad    =      M(:,11);
    
    j=j+1;
end

stnums={station.num};

sbevent=1; %Number of SB across all stations
for i=1:length(teststationlist)
    sbcnt=1; %Number of SB each station begins at 1,
%     if isempty(find(strcmp(stnums,teststationlist(i,:))))
%         txt=sprintf('Test Station %s not found', teststationlist(i,:));
%         disp(txt)
%         continue
%     end
    onshoreWD = onshore(i);
    testwinddir = ...
        table2array(station(strcmp(stnums,teststationlist(i,:))).winddir);
    testwindspd = ...
        table2array(station(strcmp(stnums,teststationlist(i,:))).windspd);
    testtemp    = ...
        table2array(station(strcmp(stnums,teststationlist(i,:))).temp);
    testdp  = ...
        table2array(station(strcmp(stnums,teststationlist(i,:))).dewtemp);
    time        = ...
        table2array(station(strcmp(stnums,teststationlist(i,:))).timeEST);
    testprecip  = ...
        table2array(station(strcmp(stnums,teststationlist(i,:))).precip);
    
    
    %TEST FOR JUST CLASSIC SB%%%%%
    for k=9:size(testwinddir,1)
        if time(k).Month<6 || time(k).Month>8 || time(k).Hour<8 || ...
                time(k).Hour>19
            continue
        end
        
        testtempchange=testtemp(k-6)-testtemp(k);

        testWDchange=180.0-abs(abs(testwinddir(k)-...
            testwinddir(k-6))-180.0);

        
        %Wind Orientation.Is wind onshore or offshore. If this value is 
        %less than 90 degrees, it is onshore. If greater than 90, it is 
        %offshore:
        TestWO=180.0-abs(abs(testwinddir(k)-onshoreWD)-180.0);
        TestWObefore=180.0-abs(abs(testwinddir(k-6)-onshoreWD)-180.0);
        WDDiff=180.0-abs(abs(testwinddir(k-6)-testwinddir(k))-180.0);

        
        
        
        
        
        sbbool = false; %reset everytime
        
        %This 'if' section is from the FL_SB_Detection_FlowChart.docx
        %VERY INEFFICIENT
        if TestWO<90.0 & testwindspd(k)>2.25%Test Wind Direction Onshore
            if TestWObefore>90.0 || (TestWObefore<90.0 & ...
                    testwindspd(k-6)<2.25)
                if testtemp(k-6)-testtemp(k)>1.0
                    if WDDiff>45.0
                        %CLASSIC SEA BREEZE
                        sbbool=true;
                    else
                        sbbool=false;
                    end
                    

                else
                    sbbool=false;

                end
            end
            
        else
            sbbool=false;
            
        end

        
        if sbbool
            date=datetime(time(k).Year,time(k).Month,time(k).Day);
            if  ~ismember(date,seabreezedate(:,i))
                seabreezedate(sbcnt,i)=date;
                seabreezetime(sbcnt,i)=time(k);
                sbcnt=sbcnt+1;

                %Build Stats Here:
                tempbefore(sbevent)=nanmean(testtemp(k-12:k));
                tempafter(sbevent)=nanmean(testtemp(k:k+12));
                winddirbefore(sbevent)=nanmean(testwinddir(k-12:k));
                winddirafter(sbevent)=nanmean(testwinddir(k:k+12));
                windspdbefore(sbevent)=nanmean(testwindspd(k-12:k));
                windspdafter(sbevent)=nanmean(testwindspd(k:k+12));
                precip12hr(sbevent)=nansum(testprecip(k:k+144));
                dptempbefore(sbevent)=nanmean(testdp(k-12:k));
                dptempafter(sbevent)=nanmean(testdp(k:k+12));
                timing(sbevent)=time(k).Hour + time(k).Minute/60.0;

                %For reference:
                datestr(sbevent)=string(date);
                timestr(sbevent)=string(time(k).Hour + time(k).Minute/60.0);
                sbstation(sbevent)=string(teststationlist(i,:));


                sbevent=sbevent+1;
            end
        end
    end
end
    
%WRITE DATA
dataframe = table(sbstation(:),datestr(:),tempbefore(:),tempafter(:),...
    winddirbefore(:),winddirafter(:),windspdbefore(:),windspdafter(:),...
    precip12hr(:),dptempbefore(:), ...
    dptempafter(:),timing(:),...
    'VariableNames',{'Station','Date','TBefore','TAfter','WDBefore',...
    'WDAfter','WSBefore','WSAfter','TwelveHrPrecip',...
    'DPTempBefore','DPTempAfter','Hour'});

writetable(dataframe,...
    "C:/Users/dpm/Documents/MATLAB/RadSBDEOSData/Results7_1_19/SBResults.csv");

dataframetimes = table(sbstation(:),datestr(:),timestr(:),...
    'VariableNames',{'Station','Date','Time'});

writetable(dataframetimes,...
    "C:/Users/dpm/Documents/MATLAB/RadSBDEOSData/Results7_1_19/SBDates.csv");
        