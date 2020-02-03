%% Post-recovery analysis of glider data from Irminger5 (June 2018 - July 2019)

%% Load data
addpath('C:/Users/palevsky/Dropbox/Irminger5/glider')
load G363.mat %this is the recovered data file (just for 363, only glider recovered)
load I5_final_telemetered.mat %final telemetered data for 363 and 453

%% Add path for functions
addpath('C:/Users/palevsky/Dropbox/MATLAB/OOI data processing/OOI_Irminger')

%% Interpolate GPS and CTD data onto all time points

G363R.lon_interp = naninterp1(G363R.time, G363R.longitude, G363R.time);
G363R.lat_interp = naninterp1(G363R.time, G363R.latitude, G363R.time);
G363R.salinity_interp = naninterp1(G363R.time, G363R.salinity, G363R.time);
G363R.pressure_interp = naninterp1(G363R.time, G363R.pressure, G363R.time);
G363R.temperature_interp = naninterp1(G363R.time, G363R.temperature, G363R.time);

G453.lon_interp = naninterp1(G453.time, G453.longitude, G453.time);
G453.lat_interp = naninterp1(G453.time, G453.latitude, G453.time);
G453.salinity_interp = naninterp1(G453.time, G453.salinity, G453.time);
G453.pressure_interp = naninterp1(G453.time, G453.pressure, G453.time);
G453.temperature_interp = naninterp1(G453.time, G453.temperature, G453.time);

%% Plot full depth data over entire deployments
figure(1); clf
    subplot(221)
plot(G363R.daten, G363R.depth_interp, 'k.'); hold on;
scatter(G363R.daten, G363R.depth_interp, [], G363R.oxygen_saturation,'filled'); colorbar; caxis([85 100])
set(gca,'YDir','reverse'); 
xlim([datenum(2018,6,8,12,0,0) datenum(2019,5,14,0,0,0)])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 363, oxygen saturation')
    subplot(222)
plot(G363R.daten, G363R.depth_interp, 'k.'); hold on;
scatter(G363R.daten, G363R.depth_interp, [], G363R.temperature,'filled'); colorbar; caxis([2.5 6.5])
set(gca,'YDir','reverse'); 
xlim([datenum(2018,6,8,12,0,0) datenum(2019,5,14,0,0,0)])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 363, temperature')
    subplot(223)
plot(G453.daten, G453.depth_interp, 'k.'); hold on;
scatter(G453.daten, G453.depth_interp, [], G453.oxygen_saturation,'filled'); colorbar; caxis([85 100])
set(gca,'YDir','reverse'); 
xlim([datenum(2018,6,8,12,0,0) datenum(2019,5,14,0,0,0)])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 453, oxygen saturation')
    subplot(224)
plot(G453.daten, G453.depth_interp, 'k.'); hold on;
scatter(G453.daten, G453.depth_interp, [], G453.temperature,'filled'); colorbar; caxis([2.5 6.5])
set(gca,'YDir','reverse'); 
xlim([datenum(2018,6,8,12,0,0) datenum(2019,5,14,0,0,0)])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 453, temperature')

%% Overwrite telemetered 363 with recovered data
G363 = G363R;

%% Determine historesis effect between climb and dive profiles
% This is used as basis for determining the lag correction
% Note that this is calculated before any P or S corrections are applied
glider_updown_compare

%% Apply the lag correction throughout the full glider deployment
    secinday = 60*60*24;
    tau = 50/secinday;
    timetol = 600*10;

    indnonan = find(~isnan(G363.oxygen_saturation));
G363.oxygen_saturation_lagcorr = NaN(height(G363),1);
G363.oxygen_concentration_lagcorr = NaN(height(G363),1);
G363.oxygen_saturation_lagcorr(indnonan) = lagCorr(G363.oxygen_saturation(indnonan), G363.time(indnonan), tau, timetol);
G363.oxygen_concentration_lagcorr(indnonan) = lagCorr(G363.oxygen_concentration(indnonan), G363.time(indnonan), tau, timetol);

    indnonan = find(~isnan(G453.oxygen_saturation));
G453.oxygen_saturation_lagcorr = NaN(height(G453),1);
G453.oxygen_concentration_lagcorr = NaN(height(G453),1);
G453.oxygen_saturation_lagcorr(indnonan) = lagCorr(G453.oxygen_saturation(indnonan), G453.time(indnonan), tau, timetol);
G453.oxygen_concentration_lagcorr(indnonan) = lagCorr(G453.oxygen_concentration(indnonan), G453.time(indnonan), tau, timetol);

%% Apply pressure and salinity corrections to lag-corrected oxygen data
%Note that internal salinity setting is 35 for 363 and 0 for 453

G363.O2_corr = aaoptode_salpresscorr(G363.oxygen_concentration_lagcorr, G363.temperature_interp, G363.salinity_interp, G363.pressure_interp, 35);
G363.O2sat_corr = G363.oxygen_saturation_lagcorr.*(1+G363.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

G453.O2_corr = aaoptode_salpresscorr(G453.oxygen_concentration_lagcorr, G453.temperature_interp, G453.salinity_interp, G453.pressure_interp, 0);
G453.O2sat_corr = G453.oxygen_saturation_lagcorr.*(1+G453.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

%% Use Winkler values to calculate initial gain corrections at the time of deployment
glider_winklercal %compare with Winkler calibration cast values

%% Apply air calibration method over entire deployment for both gliders

Yr5_met = 'deployment0005_GI01SUMO-SBD12-06-METBKA000-telemetered-metbk_a_dcl_instrument_20180608T172154.969000-20190809T080328.237000.nc';
rhcorr = 1; mindateplot = datenum(2018,6,1); maxdateplot = datenum(2019,6,1);

[T_363, med_gain_363] = aircalfun(G363, 'Glider 363', -1, Yr5_met, mindateplot, rhcorr, mindateplot, maxdateplot);
[T_453, med_gain_453] = aircalfun(G453, 'Glider 453', 1, Yr5_met, mindateplot, rhcorr, mindateplot, datenum(2019,2,1));