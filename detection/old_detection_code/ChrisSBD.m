
> % Programmed by: Chris Hughes
>
> % Updated: Feb. 2013
>
> % seabreeze_detection.m
>
> clc
> clear
>
>
> %% Loading the databases needed
> load DBBB_ALL.mat;
> load DBNG_ALL.mat;
> % load DWAR.mat;
> % load DHAR.mat;
> % load DBRG.mat;
> load DLAU_ALL.mat;
> % load DBNG.mat;
> % load DELN.mat;
> % load DJCR.mat;
> % load DRHB.mat;
> % load DSBY.mat;
> % load DSTK.mat;
> % load DSJR.mat;
> % load DIRL.mat;
> % load LWSD1.mat;
> % load UNIVDE.mat;
> % load OCIM.mat;
> %% Initilizing Variables
> CAT1(1:31,1:12,1:7,1:13) = 0; %Classic SB condition: day,month,year,station
> CAT1_count(1:12,1:7,1:13) = 0;
>
> CAT2(1:31,1:12,1:7,1:13) = 0; %West Dominant Winds
> CAT2_count(1:12,1:7,1:13) = 0;
>
> CAT3(1:31,1:12,1:7,1:13) = 0; %East Dominant Winds
> CAT3_count(1:12,1:7,1:13) = 0;
>
> CAT4(1:31,1:12,1:7,1:13) = 0; %Variable Winds
> CAT4_count(1:12,1:7,1:13) = 0;
>
> CAT5(1:31,1:12,1:7,1:13) = 0; %Dew SB condition
> CAT5_count(1:12,1:7,1:13) = 0;
>
> CAT6(1:31,1:12,1:7,1:13) = 0; %Both Classic and Dew condition
> CAT6_count(1:12,1:7,1:13) = 0;
>
> CAT7(1:31,1:12,1:7,1:13) = 0; %Missing DATA
> CAT7_count(1:12,1:7,1:13) = 0;
>
> CAT8(1:31,1:12,1:7,1:13) = 0;  %East/West Split
> CAT8_count(1:12,1:7,1:13) = 0;
>
> CAT9(1:31,1:12,1:7,1:13) = 0; %South Feature
> CAT9_count(1:12,1:7,1:13) = 0;
>
> CAT10(1:31,1:12,1:7,1:13) = 0; %SB Ending
> CAT10_count(1:12,1:7,1:13) = 0;
>
> CAT11(1:31,1:12,1:7,1:13) = 0; %2nd SB
> CAT11_count(1:12,1:7,1:13) = 0;
>
> CAT12(1:31,1:12,1:7,1:13) = 0; %2nd SB reversal
> CAT12_count(1:12,1:7,1:13) = 0;
>
> CAT13(1:31,1:12,1:7,1:13) = 0; %Const East Temp Gradiant 1
> CAT13_count(1:12,1:7,1:13) = 0;
>
> CAT14(1:31,1:12,1:7,1:13) = 0; %Const East Temp Gradiant 2
> CAT14_count(1:12,1:7,1:13) = 0;
>
> for station = 4    %DBBB
> for temp_year = 2005:2011
>     the_year = temp_year - 2004;
>     for the_month = 1:12
>
>     if the_month == 1
>         days_in_month = 31;
>     elseif the_month == 2
>         if the_year == (4 || 8)    %leap year
>             days_in_month = 29;
>         else
>             days_in_month = 28;
>         end
>     elseif the_month == 3
>         days_in_month = 31;
>     elseif the_month == 4
>         days_in_month = 30;
>     elseif the_month == 5
>         days_in_month = 31;
>     elseif the_month == 6
>         days_in_month = 30;
>     elseif the_month == 7
>         days_in_month = 31;
>     elseif the_month == 8
>         days_in_month = 31;
>     elseif the_month == 9
>         days_in_month = 30;
>     elseif the_month == 10
>         days_in_month = 31;
>     elseif the_month == 11
>         days_in_month = 30;
>     elseif the_month == 12
>         if the_year == 7
>             days_in_month = 30;
>         else
>         days_in_month = 31;
>         end
>     end
>
>     %% Loop Through the Days of the Month
>
>
>     for j = 1:days_in_month %1:days_in_month % Days of the Month
>
>         [the_month j the_year+2004]
>
>         %Finding the starting time
>         year_finder = find(DBBB_ALL.year==temp_year);
>         month_finder = find(DBBB_ALL.month==the_month);
>         day_finder = find(DBBB_ALL.day==j); %Day
>         hour_finder = find(DBBB_ALL.hour==12); %Hour
>         minute_finder = find(DBBB_ALL.minute==0); %Minute
>
>
>         aa1 = intersect(day_finder,minute_finder); %For starting time of day
>         aa2 = intersect(month_finder,year_finder);
>         aa3 = intersect(hour_finder,aa1);
>
>         critical_time(j) = intersect(aa2,aa3); %Starting time
>
>         %Loops through 12 hours (12 readings per hour)
>         number_of_readings = 144; %8AM-8PM
>
>         for i = 1:number_of_readings
>
>                 %Opening Reference Files (DLAU)
>                 reference_temp(i) = DLAU_ALL.temp(critical_time(j)+i-1);
>                 reference_east_wdir(i)= DLAU_ALL.east_wdir(critical_time(j)+i-1);
>                 reference_wind_speed(i) = DLAU_ALL.wspeed(critical_time(j)+i-1);
>
>
>
>            if station == 1
>                 test_temp(i) = LWSD1.avg_temp(critical_time(j)+i-1); %Temp.
>                 test_dew(i) = LWSD1.avg_dew(critical_time(j)+i-1);    %Dew
>                 test_east_wdir(i) = LWSD1.avg_east_wdir(critical_time(j)+i-1); %East Wind Dir.
>                 test_wind_speed(i) = LWSD1.avg_wspeed(critical_time(j)+i-1); %Wind Speed
>                 test_wdir(i) = LWSD1.wdir(critical_time(j)+i-1); %Wind Dir.
>                 test_rh(i) = 50; %filler
>            elseif station == 2
>                 test_temp(i) = DRHB.avg_temp(critical_time(j)+i-1); %Temp.
>                 test_dew(i) = DRHB.avg_dew(critical_time(j)+i-1);    %Dew
>                 test_east_wdir(i) = DRHB.avg_east_wdir(critical_time(j)+i-1); %East Wind Dir.
>                 test_wind_speed(i) = DRHB.avg_wspeed(critical_time(j)+i-1); %Wind Speed
>                 test_wdir(i) = DRHB.wdir(critical_time(j)+i-1); %Wind Dir.
>                 test_rh(i) = DRHB.avg_rh(critical_time(j)+i-1);
>             elseif station == 3
>                 test_temp(i) = DIRL.avg_temp(critical_time(j)+i-1); %Temp.
>                 test_dew(i) = DIRL.avg_dew(critical_time(j)+i-1);    %Dew
>                 test_east_wdir(i) = DIRL.avg_east_wdir(critical_time(j)+i-1); %East Wind Dir.
>                 test_wind_speed(i) = DIRL.avg_wspeed(critical_time(j)+i-1); %Wind Speed
>                 test_wdir(i) = DIRL.wdir(critical_time(j)+i-1); %Wind Dir.
>                 test_rh(i) = DIRL.avg_rh(critical_time(j)+i-1);
>             elseif station == 4
>                 test_temp(i) = DBBB_ALL.temp(critical_time(j)+i-1); %Temp.
>                 test_east_wdir(i) = DBBB_ALL.east_wdir(critical_time(j)+i-1); %East Wind Dir.
>                 test_wind_speed(i) = DBBB_ALL.wspeed(critical_time(j)+i-1); %Wind Speed
>                 test_wdir(i) = DBBB_ALL.wdir(critical_time(j)+i-1); %Wind Dir.
>             elseif station == 5
>                 test_temp(i) = DBNG.avg_temp(critical_time(j)+i-1); %Temp.
>                 test_dew(i) = DBNG.avg_dew(critical_time(j)+i-1);    %Dew
>                 test_east_wdir(i) = DBNG.avg_east_wdir(critical_time(j)+i-1); %East Wind Dir.
>                 test_wind_speed(i) = DBNG.avg_wspeed(critical_time(j)+i-1); %Wind Speed
>                 test_wdir(i) = DBNG.wdir(critical_time(j)+i-1); %Wind Dir.
>                 test_rh(i) = DBNG.avg_rh(critical_time(j)+i-1);
>               elseif station == 6
>                 test_temp(i) = DSJR.avg_temp(critical_time(j)+i-1); %Temp.
>                 test_dew(i) = DSJR.avg_dew(critical_time(j)+i-1);    %Dew
>                 test_east_wdir(i) = DSJR.avg_east_wdir(critical_time(j)+i-1); %East Wind Dir.
>                 test_wind_speed(i) = DSJR.avg_wspeed(critical_time(j)+i-1); %Wind Speed
>                 test_wdir(i) = DSJR.wdir(critical_time(j)+i-1); %Wind Dir.
>                 test_rh(i) = DSJR.avg_rh(critical_time(j)+i-1);
>               elseif station == 7
>                 test_temp(i) = DWAR.avg_temp(critical_time(j)+i-1); %Temp.
>                 test_dew(i) = DWAR.avg_dew(critical_time(j)+i-1);    %Dew
>                 test_east_wdir(i) = DWAR.avg_east_wdir(critical_time(j)+i-1); %East Wind Dir.
>                 test_wind_speed(i) = DWAR.avg_wspeed(critical_time(j)+i-1); %Wind Speed
>                 test_wdir(i) = DWAR.wdir(critical_time(j)+i-1); %Wind Dir.
>                 test_rh(i) = DWAR.avg_rh(critical_time(j)+i-1);
>               elseif station == 8
>                 test_temp(i) = DELN.avg_temp(critical_time(j)+i-1); %Temp.
>                 test_dew(i) = DELN.avg_dew(critical_time(j)+i-1);    %Dew
>                 test_east_wdir(i) = DELN.avg_east_wdir(critical_time(j)+i-1); %East Wind Dir.
>                 test_wind_speed(i) = DELN.avg_wspeed(critical_time(j)+i-1); %Wind Speed
>                 test_wdir(i) = DELN.wdir(critical_time(j)+i-1); %Wind Dir.
>                 test_rh(i) = DELN.avg_rh(critical_time(j)+i-1);
>               elseif station == 9
>                 test_temp(i) = DJCR.avg_temp(critical_time(j)+i-1); %Temp.
>                 test_dew(i) = DJCR.avg_dew(critical_time(j)+i-1);    %Dew
>                 test_east_wdir(i) = DJCR.avg_east_wdir(critical_time(j)+i-1); %East Wind Dir.
>                 test_wind_speed(i) = DJCR.avg_wspeed(critical_time(j)+i-1); %Wind Speed
>                 test_wdir(i) = DJCR.wdir(critical_time(j)+i-1); %Wind Dir.
>                 test_rh(i) = DJCR.avg_rh(critical_time(j)+i-1);
>
>                  elseif station == 10
>                 test_temp(i) = DSTK.avg_temp(critical_time(j)+i-1); %Temp.
>                 test_dew(i) = DSTK.avg_dew(critical_time(j)+i-1);    %Dew
>                 test_east_wdir(i) = DSTK.avg_east_wdir(critical_time(j)+i-1); %East Wind Dir.
>                 test_wind_speed(i) = DSTK.avg_wspeed(critical_time(j)+i-1); %Wind Speed
>                 test_wdir(i) = DSTK.wdir(critical_time(j)+i-1); %Wind Dir.
>                 test_rh(i) = DSTK.avg_rh(critical_time(j)+i-1);
>
>                  elseif station == 11
>                 test_temp(i) = DSBY.avg_temp(critical_time(j)+i-1); %Temp.
>                 test_dew(i) = DSBY.avg_dew(critical_time(j)+i-1);    %Dew
>                 test_east_wdir(i) = DSBY.avg_east_wdir(critical_time(j)+i-1); %East Wind Dir.
>                 test_wind_speed(i) = DSBY.avg_wspeed(critical_time(j)+i-1); %Wind Speed
>                 test_wdir(i) = DSBY.wdir(critical_time(j)+i-1); %Wind Dir.
>                 test_rh(i) = DSBY.avg_rh(critical_time(j)+i-1);
>
>
>             end
>
>             %% Test for Classic Sea Breeze
>
>
>             if i>12
>                 SB(i) = 0;
>                 test_temp(i)
>                 if (test_temp(i-12) > test_temp(i)+2)
>                     if (test_east_wdir(i)>0)
>                         if(test_wind_speed(i)>1.0)
>                             if(test_east_wdir(i-12)<-15 || (test_east_wdir(i-12)>-90 && test_wind_speed(i-12)<1))
>                                 if (reference_east_wdir(i)<0 || (reference_east_wdir(i)>0 && reference_wind_speed(i)<0.5))
>                                     SB(i) = 1; %Collects all of the times that this condition passes
>
>                                 end
>                             end
>                         end
>                     end
>
>
>
>                 end
>             end
>
>             %% Test for Dew Sea Breeze
>            SB2(i) = 0;
>
>             %% East / West Split
>             if i>12
>                 if (test_east_wdir(i)>0 && test_wind_speed(i)>1.0 && (reference_east_wdir(i)<0 || (reference_east_wdir(i)>0 && reference_wind_speed(i)<0.5)) )%&& (reference_east_wdir(i-12)<0 || (reference_east_wdir(i-12)>0 && reference_wind_speed(i-12)<1.0))); %&& (test_east_wdir(i-12)<0 || (test_east_wdir(i-12)>0 && test_wind_speed(i-12)<2))
>
>
>
>                     SB4(i) = 1; %Collects all of the times that this condition passes
>
>
>                 else
>                     SB4(i) = 0;
>                 end
>             end
>
>              %% East / West Split 2
>             if i>12
>                 if ((reference_temp(i) > test_temp(i)+2) && test_east_wdir(i)>0 && test_wind_speed(i)>1.0 && (reference_east_wdir(i)<0 || (reference_east_wdir(i)>0 && reference_wind_speed(i)<0.5)) )%&& (reference_east_wdir(i-12)<0 || (reference_east_wdir(i-12)>0 && reference_wind_speed(i-12)<1.0))); %&& (test_east_wdir(i-12)<0 || (test_east_wdir(i-12)>0 && test_wind_speed(i-12)<2))
>
>
>                     SB5(i) = 1; %Collects all of the times that this condition passes
>
>                 else
>                     SB5(i) = 0;
>                 end
>             end
>
>
>             %% SB Reversal
>             if i>12
>                 if (test_east_wdir(i)<0 || test_temp(i) > reference_temp(i));
>                     SB6(i) = 1;
>                 else
>                     SB6(i) = 0;
>                 end
>             end
>
>
>              %% Temp Gradiant
>             if i>12
>                 if (test_temp(i) + 3 < reference_temp(i));
>                     SB7(i) = 1;
>                 else
>                     SB7(i) = 0;
>                 end
>             end
>
>              %% East/West Wind
>             if i>12
>                 if (test_east_wdir(i)>0 && reference_east_wdir(i)<0);
>                     SB8(i) = 1;
>                 else
>                     SB8(i) = 0;
>                 end
>             end
>
>             %%
>
>             if i>12 %After atleast 1 Hour has passed (7AM+)
>                 if (SB(i) == 1 || SB2(i) == 1 ||SB4(i) == 1)
>                     SB3(i) = 1;  %Shows when either SB condition has passed
>                 else
>                     SB3(i) = 0;
>                 end
>             end
>
>
>             %WRF SB
> %             if i>12 %After atleast 1 Hour has passed (7AM+)
> %                 if (test_temp(i)+1 < reference_temp(i)&& test_east_wdir(i)>0 && reference_east_wdir(i)<0 && test_wind_speed(i).*sin(test_east_wdir(i)*pi/180)-reference_wind_speed(i).*sin(reference_east_wdir(i).*pi/180));
> %                     SB9(i) = 1;  %Shows when either SB condition has passed
> %                 else
> %                     SB9(i) = 0;
> %                 end
> %             end
>
>
>
>             %% Looks for missing data from the test station
>             if test_temp(i) > -50
>                 temp_present(i) = 1;
>             else
>                 temp_present(i) = 0;
>             end
>
>             if reference_temp(i) > -50   %TEMP NOT RH --OK!!!
>                 ref_temp_present(i) = 1;
>             else
>                 ref_temp_present(i) = 0;
>             end
>
>             if test_wdir(i) == 0   %TEMP NOT RH --OK!!!
>                 wdir_present(i) = 0;
>             else
>                 wdir_present(i) = 1;
>             end
>
>             if reference_east_wdir(i) == 0   %TEMP NOT RH --OK!!!
>                 ref_wdir_present(i) = 0;
>             else
>                 ref_wdir_present(i) = 1;
>             end
>
>
>             %% Looks for West Winds from the test station
>             if (test_east_wdir(i)<0)
>                 NO_SB(i) = 1;
>             else
>                 NO_SB(i) = 0;
>             end
>
>             %% Looks for East Winds from the reference station
>             if (reference_east_wdir(i)>0)
>                 REF_EAST(i) = 1;
>             else
>                 REF_EAST(i) = 0;
>             end
>
>             %% Looks for East Winds from the test station
>             if (test_east_wdir(i)>0)
>                 TEST_EAST(i) = 1;
>             else
>                 TEST_EAST(i) = 0;
>             end
>
>             %% Looks for West Winds from the test station
>             if (test_east_wdir(i)<0)
>                 TEST_WEST(i) = 1;
>             else
>                 TEST_WEST(i) = 0;
>             end
>
>             %% Gets timing of first instance of classic SB (UTC)
>             if i>12
>                 if (sum(SB(i)) == 1 && CAT1(j,the_month,the_year,station) == 0)
>                     CAT1(j,the_month,the_year,station) = 8 + ((i-1)/12)  ;
>                     CAT1_count(the_month,the_year,station) = CAT1_count(the_month,the_year,station) + 1;
>                 end
>             end
>
>             %% Gets timing of first instance of Dew SB (UTC)
>             if i>12
>                 if sum(SB2(i)) == 1 && CAT5(j,the_month,the_year,station) == 0
>                     CAT5(j,the_month,the_year,station) = 0;
>                     CAT5_count(the_month,the_year,station) = 0;
>                 end
>             end
>
>
>             %% Gets timing of first instance of either SB condition (UTC)
>             if i>12
>                 if sum(SB3(i)) == 1 && CAT6(j,the_month,the_year,station) == 0
>                   CAT6(j,the_month,the_year,station) =  8 + ((i-1)/12)  ;
>                     CAT6_count(the_month,the_year,station) = CAT6_count(the_month,the_year,station) + 1;
>                 end
>             end
>
>
>             %% EAST/WEST Split
>             if i>12
>                 if sum(SB4(i)) == 1 && CAT8(j,the_month,the_year,station) == 0
>                    CAT8(j,the_month,the_year,station) =  8 + ((i-1)/12)  ;
>                     CAT8_count(the_month,the_year,station) = CAT8_count(the_month,the_year,station) + 1;
>                 end
>             end
>
>             %% EAST/WEST Part 2
>             if i>12
>                 if sum(SB5(i)) == 1 && CAT9(j,the_month,the_year,station) == 0
>                    CAT9(j,the_month,the_year,station) =  8 + ((i-1)/12)  ;
>                     CAT9_count(the_month,the_year,station) = CAT9_count(the_month,the_year,station) + 1;
>                 end
>             end
>
>             %% East Wind Temp Gradiant 1
>             if i>12
>                 if sum(SB7(i)) == 1 && CAT13(j,the_month,the_year,station) == 0
>                  CAT13(j,the_month,the_year,station) = 8 + (i/12);
>                     CAT13_count(the_month,the_year,station) = CAT13_count(the_month,the_year,station) + 1;
>                 end
>             end
>
>             %% East Wind Temp Gradiant 2
>             if i>12
>                 if sum(SB8(i)) == 1 && CAT14(j,the_month,the_year,station) == 0
>                   CAT14(j,the_month,the_year,station) = 8 + (i/12);
>                     CAT14_count(the_month,the_year,station) = CAT14_count(the_month,the_year,station) + 1;
>                 end
>             end
>             %% SB reversal
>             if i>12
>                 if sum(SB6(i)) == 1 && max(SB(1:i)) == 1 && CAT10(j,the_month,the_year,station) == 0 %COULD USE SB3
>                  CAT10(j,the_month,the_year,station) = 8 + (i/12);
>                     CAT10_count(the_month,the_year,station) = CAT10_count(the_month,the_year,station) + 1;
>                 end
>             end
>
>
>
>             %% Gets timing of second instance of either SB condition (UTC)
>             if i>12
>                 if sum(SB(i)) == 1 && CAT10(j,the_month,the_year,station) > 1 && CAT11(j,the_month,the_year,station) == 0
>                     CAT11(j,the_month,the_year,station) = 8 + (i/12);
>                     CAT11_count(the_month,the_year,station) = CAT11_count(the_month,the_year,station) + 1;
>                 end
>             end
>
>             %% Gets timing of reversal of 2nd SB (UTC) (TEMP)
>             if i>12
>                 if sum(SB6(i)) == 1 && CAT12(j,the_month,the_year,station) == 0 && CAT11(j,the_month,the_year,station) > 0
>                  CAT12(j,the_month,the_year,station) = 8 + (i/12);
>                     CAT12_count(the_month,the_year,station) = CAT12_count(the_month,the_year,station) + 1;
>                 end
>             end
>
>
>
>
>
>             %% Searching for prevailing conditions
>             if i == number_of_readings
>                 if sum(TEST_EAST(1:144))>134 && CAT6(j,the_month,the_year,station) == 0
>                     % East Wind
>                     CAT3(j,the_month,the_year,station) = 8 + (i/12);
>                     CAT3_count(the_month,the_year,station) = CAT3_count(the_month,the_year,station) + 1;
>                 elseif sum(TEST_WEST(1:144))>134 && CAT6(j,the_month,the_year,station) == 0
>                      % West Wind
>                     CAT2(j,the_month,the_year,station) = 8 + (i/12);
>                     CAT2_count(the_month,the_year,station) = CAT2_count(the_month,the_year,station) + 1;
>                 elseif (CAT6(j,the_month,the_year,station) == 0)
>                     % Variable Wind
>                    CAT4(j,the_month,the_year,station) = 8 + (i/12);
>                     CAT4_count(the_month,the_year,station) = CAT4_count(the_month,the_year,station) + 1;
>                 end
>
>                 % Changing condition to 'Missing' if more than
>                 % 10% of the data is missing for the day
>                 if (sum(temp_present(1:144))<116 || sum(ref_temp_present(1:144))<116 || sum(wdir_present(1:144))<116 || sum(ref_wdir_present(1:144))<116)  %%%% ||sum(ref_rh_present(1:144))<130 ||sum(test_rh_present(1:144))<130
>                     CAT1(j,the_month,the_year,station) = 0;
>                     CAT2(j,the_month,the_year,station) = 0;
>                     CAT3(j,the_month,the_year,station) = 0;
>                     CAT4(j,the_month,the_year,station) = 0;
>                     CAT5(j,the_month,the_year,station) = 0;
>                     CAT6(j,the_month,the_year,station) = 0;
>                     CAT7(j,the_month,the_year,station) = 1;
>                     CAT8(j,the_month,the_year,station) = 0;
>                     CAT9(j,the_month,the_year,station) = 0;
>                     CAT10(j,the_month,the_year,station) = 0;
>                     CAT11(j,the_month,the_year,station) = 0;
>                     CAT12(j,the_month,the_year,station) = 0;
>                     CAT13(j,the_month,the_year,station) = 0;
>                     CAT14(j,the_month,the_year,station) = 0;
>                 end
>             end
>         end %each reading
>     end %day
> end %station
> end
>
> end
>
> %% Summary
> % Renaming Categories using an a(i,j) matrix
>     for i = 1:days_in_month
>        for i2 = 1:12;
>            for i3 = 1:7
>                for i4 = 1:11
>
>
>         if CAT1(i,i2,i3,i4) > 0
>             a(i,i2,i3,i4) = 1; %Classic SB
>             if CAT5(i,i2,i3,i4) > 1
>                 a(i,i2,i3,i4) = 1; %Classic + Dew SB
>             end
>         elseif CAT5(i,i2,i3,i4) > 0
>             a(i,i2,i3,i4) = 2; %Dew SB
>         elseif CAT2(i,i2,i3,i4) > 0
>             a(i,i2,i3,i4) = 4; %West Wind
>         elseif CAT3(i,i2,i3,i4) > 0
>             a(i,i2,i3,i4) = 5; %East Wind
>         elseif CAT4(i,i2,i3,i4) > 0
>             a(i,i2,i3,i4) = 6; %Variable Wind
>         elseif CAT8(i,i2,i3,i4) > 0
>             a(i,i2,i3,i4) = 7; %East/West Split
>         elseif CAT9(i,i2,i3,i4) > 0
>             a(i,i2,i3,i4) = 8; %South Feature
>         elseif CAT10(i,i2,i3,i4) > 0
>             a(i,i2,i3,i4) = 9; %SB Reversal
>         elseif CAT7(i,i2,i3,i4) > 0
>             a(i,i2,i3,i4) = -99; %Missing DATA
>         end
>                end
>            end
>        end
>
>     end
>
>
> % save a
> % save CAT1
> % save CAT2
> % save CAT3
> % save CAT4
> % save CAT5
> % save CAT6
> % save CAT7
> % save CAT8
> % save CAT9
> % save CAT10
>
>