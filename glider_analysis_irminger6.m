%Main script for analyzing OOI Irminger-6 glider data
load ('latest');

%% Correct glider data for S and pressure
%Note that internal salinity setting is 35 for 525 and 0 for 560

G525.lon_interp = naninterp1(G525.time, G525.longitude, G525.time);
G525.lat_interp = naninterp1(G525.time, G525.latitude, G525.time);
G525.salinity_interp = naninterp1(G525.time, G525.salinity, G525.time);
G525.pressure_interp = naninterp1(G525.time, G525.pressure, G525.time);
G525.temperature_interp = naninterp1(G525.time, G525.temperature, G525.time);
G525.O2_corr = aaoptode_salpresscorr(G525.oxygen_concentration, G525.temperature_interp, G525.salinity_interp, G525.pressure_interp, 0);
G525.O2sat_corr = G525.oxygen_saturation.*(1+G525.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

G560.lon_interp = naninterp1(G560.time, G560.longitude, G560.time);
G560.lat_interp = naninterp1(G560.time, G560.latitude, G560.time);
G560.salinity_interp = naninterp1(G560.time, G560.salinity, G560.time);
G560.pressure_interp = naninterp1(G560.time, G560.pressure, G560.time);
G560.temperature_interp = naninterp1(G560.time, G560.temperature, G560.time);
G560.O2_corr = aaoptode_salpresscorr(G560.oxygen_concentration, G560.temperature_interp, G560.salinity_interp, G560.pressure_interp, 0);
G560.O2sat_corr = G560.oxygen_saturation.*(1+G560.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

%% Visualize data diagnostics
ctd_Irminger6 %process all CTD cast data, and then close plot windows to move on to glider analysis
    close all
glider_aircal_irminger6 %air calibration
glider_updown_compare_irminger6 %historesis effects between up and down profiles
glider_plotmap_irminger6 %plot glider locations during cruise
%glider_winklercal %compare with Winkler calibration cast values

% %% Compare with wire-following profiler - cut this b/c didn't get high
% enough res data in Yr 5, unlikely to do along way in Yr 6
% wfp_analysis %Load wfp data and do initial processings
%     dist_compare = 4; %Select radius around mooring to use for comparison (kilometers)
%     time_compare = 1; %Select time difference between mooring and glider profiles to keep for comparison (days)
%     % Need to select profile direction (-1 == up 1 == down); -1 for 525 and 1 for 560
% [G525_profile_summary, wfp_profile_summary_525] = glider_wfp_compare(G525, Yr5_wfp, dist_compare, time_compare, -1);
% [G560_profile_summary, wfp_profile_summary_560] = glider_wfp_compare(G560, Yr5_wfp, dist_compare, time_compare, 1);