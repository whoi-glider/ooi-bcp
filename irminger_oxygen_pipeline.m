%% Overview script for full oxygen calibration pipeline
addpath('C:/Users/palevsky/Dropbox/MATLAB/OOI data processing/OOI_Irminger')

%% Read in and process turn-around cruise data for deep isotherm characterization
ISO = 3.1; %Select deep isotherm to use for salinity and oxygen correction
cruise_casts

%% Read in and process wire-following profiler data
%Applies lag correction for oxygen, backscatter spike analysis, and outlier filtering
wfp_analysis

%% Apply salinity correction on selected deep isotherm
wfp_deepisotherm_salcorr

%% Recalculate oxygen using corrected salinity, and then apply deep isotherm correction
wfp_deepisotherm_O2corr

%% Merge corrected output into a single structure for each of optode and fluorometer data
wfp_mergeyears

%Output for further analysis is in wggmerge and wggmerge_fl - different
%number of profiles removed based on anomalous data, so different sizes

%% Save wfp merged output
save wfpmerge_output.mat wggmerge wggmerge_fl HYPMlat HYPMlon

%% Read in and process glider data (all deployments, years 2-8)
%Includes air calibration calculations and assessment when gliders were configured for air oxygen measurements (year 5 on)
glider_analysis

%%
glider_wfp_comparison




%% Gain corrections with cruise data processed by Kristen
% cruise_oxygen
% close all;
% 
% %% Analysis on deep isotherms
% wfp_deepisotherms


%save cruise_oxygen_output.mat castsum btlsum