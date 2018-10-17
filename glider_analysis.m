%Main script for analyzing OOI Irminger-5 glider data
load ('latest');

%% Correct glider data for S and pressure
%Note that internal salinity setting is 35 for 363 and 0 for 453

G363.lon_interp = naninterp1(G363.time, G363.longitude, G363.time);
G363.lat_interp = naninterp1(G363.time, G363.latitude, G363.time);
G363.salinity_interp = naninterp1(G363.time, G363.salinity, G363.time);
G363.pressure_interp = naninterp1(G363.time, G363.pressure, G363.time);
G363.temperature_interp = naninterp1(G363.time, G363.temperature, G363.time);
G363.O2_corr = aaoptode_salpresscorr(G363.oxygen_concentration, G363.temperature_interp, G363.salinity_interp, G363.pressure_interp, 35);
G363.O2sat_corr = G363.oxygen_saturation.*(1+G363.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

G453.lon_interp = naninterp1(G453.time, G453.longitude, G453.time);
G453.lat_interp = naninterp1(G453.time, G453.latitude, G453.time);
G453.salinity_interp = naninterp1(G453.time, G453.salinity, G453.time);
G453.pressure_interp = naninterp1(G453.time, G453.pressure, G453.time);
G453.temperature_interp = naninterp1(G453.time, G453.temperature, G453.time);
G453.O2_corr = aaoptode_salpresscorr(G453.oxygen_concentration, G453.temperature_interp, G453.salinity_interp, G453.pressure_interp, 0);
G453.O2sat_corr = G453.oxygen_saturation.*(1+G453.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

%% Visualize data diagnostics
glider_aircal %air calibration
glider_updown_compare %historesis effects between up and down profiles
glider_winklercal %compare with Winkler calibration cast values

%% Compare with wire-following profiler
wfp_analysis %Load wfp data and do initial processings
    dist_compare = 4; %Select radius around mooring to use for comparison (kilometers)
    time_compare = 1; %Select time difference between mooring and glider profiles to keep for comparison (days)
    % Need to select profile direction (-1 == up 1 == down); -1 for 363 and 1 for 453
[G363_profile_summary, wfp_profile_summary_363] = glider_wfp_compare(G363, Yr5_wfp, dist_compare, time_compare, -1);
[G453_profile_summary, wfp_profile_summary_453] = glider_wfp_compare(G453, Yr5_wfp, dist_compare, time_compare, 1);