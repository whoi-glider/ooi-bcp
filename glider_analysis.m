%% Main script for analyzing OOI Irminger glider data

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

%% Close plots from air cal function
close all

%% Calculate lag correction
%glider_lagAssess
%Can avoid running full lag assessment to save time - output boundary layer thickness below
tau_in = 66.8878;

%% Apply lag correction to glider data
%For each glider dataset, reshape data into format for Gordon functions, then apply glider_lagCorrectFun
addpath(genpath('C:\Users\palevsky\Documents\GitHub\optode-response-time'))

[Yr5.glg363] = glider_reshape(Yr5.G363);
[Yr5.glg363.doxy_lagcorr] = glider_lagCorrectFun(Yr5.glg363, tau_in);

[Yr5.glg453] = glider_reshape(Yr5.G453);
[Yr5.glg453.doxy_lagcorr] = glider_lagCorrectFun(Yr5.glg453, tau_in);

[Yr6.glg525] = glider_reshape(Yr6.G525);
[Yr6.glg525.doxy_lagcorr] = glider_lagCorrectFun(Yr6.glg525, tau_in);

[Yr6.glg560] = glider_reshape(Yr6.G560);
[Yr6.glg560.doxy_lagcorr] = glider_lagCorrectFun(Yr6.glg560, tau_in);

[Yr7.glg515] = glider_reshape(Yr7.G515);
[Yr7.glg515.doxy_lagcorr] = glider_lagCorrectFun(Yr7.glg515, tau_in);

[Yr7.glg365] = glider_reshape(Yr7.G365);
[Yr7.glg365.doxy_lagcorr] = glider_lagCorrectFun(Yr7.glg365, tau_in);

[Yr8.glg469] = glider_reshape(Yr8.G469);
[Yr8.glg469.doxy_lagcorr] = glider_lagCorrectFun(Yr8.glg469, tau_in);

[Yr8.glg537] = glider_reshape(Yr8.G537);
[Yr8.glg537.doxy_lagcorr] = glider_lagCorrectFun(Yr8.glg537, tau_in);

[Yr8.glg565] = glider_reshape(Yr8.G565);
[Yr8.glg565.doxy_lagcorr] = glider_lagCorrectFun(Yr8.glg565, tau_in);


%% Identify and flag outliers and spikes in data, and regrid on pressure surfaces (following same approach as wfp_analysis)

%tolerances for spikes (2 x those chosen for WFP for oxy and 20x for
%temp and sal, given higher variability near surface)
oxyspike = 6;
tempspike = 1;
salspike = 0.1;
oxymin = 240;
oxymax = 360;
tempmin = 3;
tempmax = 12;
salmin = 32;
salmax = 35.1;


%Select depth resolution and smoothing - current setting is 1 m resolution
%w/ 5-m smoothing
pres_grid = [1:1:1000];
smval = 5; %points to smooth over

[Yr5.glg363,Yr5.glg363.fract_flag] = glider_flag_regrid(Yr5.glg363, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr5.glg453,Yr5.glg453.fract_flag] = glider_flag_regrid(Yr5.glg453, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr6.glg525,Yr6.glg525.fract_flag] = glider_flag_regrid(Yr6.glg525, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr6.glg560,Yr6.glg560.fract_flag] = glider_flag_regrid(Yr6.glg560, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr7.glg515,Yr7.glg515.fract_flag] = glider_flag_regrid(Yr7.glg515, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr7.glg365,Yr7.glg365.fract_flag] = glider_flag_regrid(Yr7.glg365, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr8.glg469,Yr8.glg469.fract_flag] = glider_flag_regrid(Yr8.glg469, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr8.glg537,Yr8.glg537.fract_flag] = glider_flag_regrid(Yr8.glg537, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr8.glg565,Yr8.glg565.fract_flag] = glider_flag_regrid(Yr8.glg565, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);

%Consider whether to also regrid on isotherms?

%% Combine glider datasets
glgmerge{1} = Yr5.glg363;
glgmerge{2} = Yr5.glg453;
glgmerge{3} = Yr6.glg525;
glgmerge{4} = Yr6.glg560;
glgmerge{5} = Yr7.glg515;
glgmerge{6} = Yr7.glg365;
glgmerge{7} = Yr8.glg469;
%glgmerge{8} = Yr8.glg537; %Removing from further analysis because all
%salinity data flagged as outlier, and only lasted a few weeks anyway
glgmerge{8} = Yr8.glg565;

glidertitles = [{'Glider 363, Year 5','Glider 453, Year 5','Glider 560, Year 6','Glider 515, Year 6',...
    'Glider 515, Year 7','Glider 365, Year 7','Glider 469, Year 8','Glider 565, Year 8'}];

%% Add air cal output to combined glider datasets
glgmerge{1}.Taircal = Yr5.T_363;
glgmerge{2}.Taircal = Yr5.T_453;
glgmerge{3}.Taircal = Yr6.T_525;
glgmerge{4}.Taircal = Yr6.T_560;
glgmerge{5}.Taircal = Yr7.T_515;
glgmerge{6}.Taircal = Yr7.T_365;
glgmerge{7}.Taircal = Yr8.T_469;

%% Assessment histograms
figure(1); clf
subplot(311)
for i = 1:length(glgmerge)
    histogram(glgmerge{i}.doxy_lagcorr_grid(:)); hold on;
end
legend(glidertitles,'location','NW')
xlabel('Dissolved oxygen, uncalibrated (\mumol/kg)');
xlim([oxymin oxymax])

subplot(312)
for i = 1:length(glgmerge)
    histogram(glgmerge{i}.temp_grid(:)); hold on;
end
legend(glidertitles,'location','NW')
xlabel('Temperature (^oC)');
xlim([tempmin tempmax])

subplot(313)
for i = 1:length(glgmerge)
    histogram(glgmerge{i}.sal_grid(:)); hold on;
end
legend(glidertitles,'location','NW')
xlabel('Salinity (PSU)');
xlim([salmin + 2.5 salmax])

