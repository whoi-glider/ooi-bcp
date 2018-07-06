%Main script for analyzing OOI Irminger-5 glider data
load 'C:/Users/Hilary/Google Drive/2018 Irminger/Data_Gliders/latest.mat'

%% Correct glider data for S and pressure
%Note that internal salinity setting is 35 for 363 and 0 for 453

G363.salinity_interp = naninterp1(G363.time, G363.salinity, G363.time);
G363.pressure_interp = naninterp1(G363.time, G363.pressure, G363.time);
G363.temperature_interp = naninterp1(G363.time, G363.temperature, G363.time);
G363.O2_corr = aaoptode_salpresscorr(G363.oxygen_concentration, G363.temperature_interp, G363.salinity_interp, G363.pressure_interp, 35);

G453.salinity_interp = naninterp1(G453.time, G453.salinity, G453.time);
G453.pressure_interp = naninterp1(G453.time, G453.pressure, G453.time);
G453.temperature_interp = naninterp1(G453.time, G453.temperature, G453.time);
G453.O2_corr = aaoptode_salpresscorr(G453.oxygen_concentration, G453.temperature_interp, G453.salinity_interp, G453.pressure_interp, 0);

%% Visualize data diagnostics
glider_aircal %air calibration
glider_updown_compare %historesis effects between up and down profiles