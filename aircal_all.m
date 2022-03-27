%% Apply air calibration method over entire deployment for all gliders

%% Add paths to data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\METBKA') %SUMO met buoy data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\From_Roo') %Processed glider data

%% Analysis of METBKA data to determine coverage of recovered vs telemetered files:
%metbka_assess

%Using all recovered files:
Yr5.met = 'deployment0005_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20180608T172109.234000-20190809T080351.721000.nc';
Yr6.met = 'deployment0006_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20190805T152914.351000-20200826T104006.497000.nc';
Yr7.met = 'deployment0007_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20200817T173428.501000-20210819T064814.736000.nc';

%% Calculate air calibration - Year 5 (2018 deployment)
load G453.mat
load G363.mat
G363 = G363R;

rhcorr = 1; mindateplot = datenum(2018,6,1);

[Yr5.T_363, Yr5.med_gain_363] = aircalfun(G363, 'Glider 363', -1, Yr5.met, mindateplot, rhcorr, mindateplot, datenum(2019,6,1));
[Yr5.T_453, Yr5.med_gain_453] = aircalfun(G453, 'Glider 453', 1, Yr5.met, mindateplot, rhcorr, mindateplot, datenum(2019,2,1));

clear G363 G363R G453
%% Calculate air calibration - Year 6 (2019 deployment)
load G525.mat %data until 3 April 2020
load G560.mat %data until Oct. 2019 (O2 sensor failed - other data ok through)

rhcorr = 1; mindateplot = datenum(2019,8,1);

[Yr6.T_560, Yr6.med_gain_560] = aircalfun(G560, 'Glider 560', 1, Yr6.met, mindateplot, rhcorr, mindateplot, datenum(2019,10,15));
[Yr6.T_525, Yr6.med_gain_525] = aircalfun(G525, 'Glider 525', -1, Yr6.met, mindateplot, rhcorr, mindateplot, datenum(2020,4,15));

clear G525 G560
%% Calculate air calibration - Year 7 (2020 deployment)
load G515.mat %data through 15 June 2021
load G365.mat %data through 21 Nov 2020

rhcorr = 1; mindateplot = datenum(2020,8,1);

[Yr7.T_515, Yr7.med_gain_515] = aircalfun(G515, 'Glider 515', 1, Yr7.met, mindateplot, rhcorr, mindateplot, datenum(2021,6,30));
[Yr7.T_365, Yr7.med_gain_365] = aircalfun(G365, 'Glider 365', 1, Yr7.met, mindateplot, rhcorr, mindateplot, datenum(2020,11,30));

clear G515 G560

%% Scatter plots (to inspect before clearing)
% 
% beg = 10000;
% jump = 30;
% jump2 = 2;
% 
% figure(100); clf
% plot(G525.daten(beg:jump2:end/2), G525.depth_interp(beg:jump2:end/2), 'k.'); hold on;
% scatter(G525.daten(beg:jump:end/2), G525.depth_interp(beg:jump:end/2), [], G525.oxygen_saturation(beg:jump:end/2),'filled'); colorbar; caxis([85 100])
% set(gca,'YDir','reverse'); 
% datetick('x',2,'keeplimits')
% ylabel('Depth (m)')
% title('Glider 525, oxygen saturation')
% 
% figure(103); clf
% plot(G365.daten, G365.depth_interp, 'k.'); hold on;
% scatter(G365.daten, G365.depth_interp, [], G365.oxygen_saturation,'filled'); colorbar; caxis([85 100])
% set(gca,'YDir','reverse'); 
% datetick('x',2,'keeplimits')
% ylabel('Depth (m)')
% title('Glider 365, oxygen saturation')
