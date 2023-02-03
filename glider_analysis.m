%% Main script for analyzing OOI Irminger glider data
%Mirrors analysis in aircal_all, with additions for further processing

%% Add paths to data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\METBKA') %SUMO met buoy data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\From_Roo') %Processed glider data

%% Identify METBKA files to use for air calibration
% See metbka_assess for analysis of METBKA data to determine coverage of recovered vs telemetered files

Yr5.met = 'deployment0005_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20180608T172109.234000-20190809T080351.721000.nc';
Yr6.met = 'deployment0006_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20190805T152914.351000-20200826T104006.497000.nc';
Yr7.met = 'deployment0007_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20200817T173428.501000-20210819T064814.736000.nc';
Yr8.met = 'deployment0008_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20210812T170400.699000-20211103T161130.307000.nc';

%% Year 5 (2018 deployment) - initial air calibration & corrections for S & pressure
load G453.mat
load G363.mat
G363 = G363R;

%Air calibration
rhcorr = 1; mindateplot = datenum(2018,6,1);
[Yr5.T_363, Yr5.med_gain_363] = aircalfun(G363, 'Glider 363', -1, Yr5.met, mindateplot, rhcorr, mindateplot, datenum(2019,6,1));
[Yr5.T_453, Yr5.med_gain_453] = aircalfun(G453, 'Glider 453', 1, Yr5.met, mindateplot, rhcorr, mindateplot, datenum(2019,2,1));

%Interpolation and S & P corrections - Note that internal salinity setting is 35 for 363 and 0 for 453
Yr5.G363 = glider_interpCorrFun(G363, 35);
Yr5.G453 = glider_interpCorrFun(G453, 0);

clear G363 G363R G453

%% Year 6 (2019 deployment) - initial air calibration & corrections for S & pressure
load G525.mat %data until 3 April 2020
load G560.mat %data until Oct. 2019 (O2 sensor failed - other data ok through)

rhcorr = 1; mindateplot = datenum(2019,8,1);

[Yr6.T_560, Yr6.med_gain_560] = aircalfun(G560, 'Glider 560', 1, Yr6.met, mindateplot, rhcorr, mindateplot, datenum(2019,10,15));
[Yr6.T_525, Yr6.med_gain_525] = aircalfun(G525, 'Glider 525', -1, Yr6.met, mindateplot, rhcorr, mindateplot, datenum(2020,4,15));

%Interpolation and S & P corrections - Internal salinity setting is 0
Yr6.G525 = glider_interpCorrFun(G525, 0);
Yr6.G560 = glider_interpCorrFun(G560, 0);

clear G525 G560

%% Year 7 (2020 deployment) - initial air calibration & corrections for S & pressure
load G515.mat %data through 15 June 2021
load G365.mat %data through 21 Nov 2020

rhcorr = 1; mindateplot = datenum(2020,8,1);

[Yr7.T_515, Yr7.med_gain_515] = aircalfun(G515, 'Glider 515', 1, Yr7.met, mindateplot, rhcorr, mindateplot, datenum(2021,6,30));
[Yr7.T_365, Yr7.med_gain_365] = aircalfun(G365, 'Glider 365', 1, Yr7.met, mindateplot, rhcorr, mindateplot, datenum(2020,11,30));

%Interpolation and S & P corrections - Internal salinity setting should be 0
Yr7.G515 = glider_interpCorrFun(G515, 0);
Yr7.G365 = glider_interpCorrFun(G365, 0);

clear G515 G365

%% Year 8 (2021 deployment) - initial air calibration & corrections for S & pressure
load GL469.mat; G469 = T; %11 Aug to 18 Oct 2021
load GL537.mat; G537 = T; %11 Aug to 5 Sept 2021
load PG565.mat; G565 = T; %30 July 2021 to 10 Jan 2022 - no air cal

rhcorr = 1; mindateplot = datenum(2021,7,30);
[Yr8.T_469, Yr7.med_gain_515] = aircalfun(G469, 'Glider 469', 1, Yr8.met, mindateplot, rhcorr, mindateplot, datenum(2021,10,20));
[Yr8.T_537, Yr7.med_gain_537] = aircalfun(G537, 'Glider 537', 1, Yr8.met, mindateplot, rhcorr, datenum(2021,8,11), datenum(2021,8,24));
%[Yr8.T_565, Yr7.med_gain_565] = aircalfun(G565, 'Glider 565', 1, Yr8.met, mindateplot, rhcorr, mindateplot, datenum(2022,1,12)); %doesn't actuallyhave air cal, so this just shows that air & surf water match

%Interpolation and S & P corrections - Internal salinity setting should be 0
Yr8.G469 = glider_interpCorrFun(G469, 0);
Yr8.G537 = glider_interpCorrFun(G537, 0);
Yr8.G565 = glider_interpCorrFun(G565, 0);

clear G469 G537 G565 T

%% Tried to save output of above b/c takes 10 min to run, but had a glitch

%% Lag correction
% Correct for historesis effects between up and down profiles using Gordon et al. 2020 approach, as in wfp_lag
addpath(genpath('C:\Users\Palevsky\Documents\GitHub\optode-response-time'))

[Yr5.Lag453] = glider_lagAssessFun(Yr5.G453, 26000, 39400, 1010, 'Glider 453, Year 5', [43:51]);

[Yr5.Lag363] = glider_lagAssessFun(Yr5.G363, 10000, 120000, 1010, 'Glider 363, Year 5', [13:20]); %just first set

[Yr6.Lag525] = glider_lagAssessFun(Yr6.G525, 800, 68200, 1010, 'Glider 525, Year 6', [7:17]);

%Gliders 469, 537, and 565 Year 8 - no initial lag assessment

%% Plotting used in function and in manually finding ranges for each glider

GLin = Yr7.G365;
beg = 1000;
stop = 90000;
jump = 2;
jump2 = 2;
name = 'Glider 365, Year 7';

figure(100); clf
plot(GLin.daten(beg:jump2:stop), GLin.depth_interp(beg:jump2:stop), 'k.'); hold on;
scatter(GLin.daten(beg:jump:stop), GLin.depth_interp(beg:jump:stop), [], GLin.oxygen_saturation(beg:jump:stop),'filled'); colorbar; caxis([85 100])
set(gca,'YDir','reverse'); 
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
ylim([-10 1010])
title({['Glider ' name ', Initial paired up and down profiles (oxygen % saturation)']})



%% Next steps

% 2) Identify cross-calibration opportunities
    % a) WFP aligned profiles
    % b) Turn-around cruise casts (Winkler + SBE43 on CTD)
    % c) Intercalibration with FLMA and FLMB

% %% Visualize data diagnostics
% glider_aircal %air calibration
% glider_updown_compare %historesis effects between up and down profiles
% glider_winklercal %compare with Winkler calibration cast values
% 
% %% Compare with wire-following profiler
% wfp_analysis %Load wfp data and do initial processings
%     dist_compare = 4; %Select radius around mooring to use for comparison (kilometers)
%     time_compare = 1; %Select time difference between mooring and glider profiles to keep for comparison (days)
%     % Need to select profile direction (-1 == up 1 == down); -1 for 363 and 1 for 453
% [G363_profile_summary, wfp_profile_summary_363] = glider_wfp_compare(G363, Yr5_wfp, dist_compare, time_compare, -1);
% [G453_profile_summary, wfp_profile_summary_453] = glider_wfp_compare(G453, Yr5_wfp, dist_compare, time_compare, 1);

