%% Post-recovery analysis of glider data from Irminger5 (June 2018 - July 2019)

%% Load data
addpath('C:/Users/palevsky/Dropbox/Irminger5/glider')
load G363.mat %this is the recovered data file (just for 363, only glider recovered)
load I5_final_telemetered.mat %final telemetered data for 363 and 453

%% Add path for functions
addpath('C:/Users/palevsky/Dropbox/MATLAB/OOI data processing/OOI_Irminger')

%% Correct glider data for S and pressure
%Note that internal salinity setting is 35 for 363 and 0 for 453

G363R.lon_interp = naninterp1(G363R.time, G363R.longitude, G363R.time);
G363R.lat_interp = naninterp1(G363R.time, G363R.latitude, G363R.time);
G363R.salinity_interp = naninterp1(G363R.time, G363R.salinity, G363R.time);
G363R.pressure_interp = naninterp1(G363R.time, G363R.pressure, G363R.time);
G363R.temperature_interp = naninterp1(G363R.time, G363R.temperature, G363R.time);
G363R.O2_corr = aaoptode_salpresscorr(G363R.oxygen_concentration, G363R.temperature_interp, G363R.salinity_interp, G363R.pressure_interp, 35);
G363R.O2sat_corr = G363R.oxygen_saturation.*(1+G363R.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

G453.lon_interp = naninterp1(G453.time, G453.longitude, G453.time);
G453.lat_interp = naninterp1(G453.time, G453.latitude, G453.time);
G453.salinity_interp = naninterp1(G453.time, G453.salinity, G453.time);
G453.pressure_interp = naninterp1(G453.time, G453.pressure, G453.time);
G453.temperature_interp = naninterp1(G453.time, G453.temperature, G453.time);
G453.O2_corr = aaoptode_salpresscorr(G453.oxygen_concentration, G453.temperature_interp, G453.salinity_interp, G453.pressure_interp, 0);
G453.O2sat_corr = G453.oxygen_saturation.*(1+G453.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

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


%% Overwrite telemetered 363 with recovered to run original scripts
G363 = G363R;
glider_updown_compare %historesis effects between up and down profiles
glider_winklercal %compare with Winkler calibration cast values