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
%save('wfpmerge_output.mat', 'wggmerge', 'wggmerge_fl', 'HYPMlat', 'HYPMlon', 'pres_grid_hypm', '-v7.3');

%% Read in and process glider data (all deployments, years 2-8)
%Includes air calibration calculations and assessment when gliders were configured for air oxygen measurements (year 5 on)
glider_analysis

%% Calculate glider gain from intercomparison with deep isotherm corrected WFP
%Plots gain corrections from WFP swim-by points as compared with air-cal, where available
glider_wfp_comparison

%% Create merged glider products
%Synthesizes gridded, deep isotherm oxygen-calibrated data, picking best glider from each year
glider_mergeyears

%Create gridded, deep isotherm oxygen-calibrated dataset for all gliders
glider_gridall

%% Save glider merged output
%save('glidermerge_output.mat', 'glidermerge', 'pres_grid_glider', '-v7.3');
%save('glider_all_output.mat', 'glgmerge', '-v7.3');
%save('glider_griddall.mat', 'glidergrid', 'pres_grid_glider', '-v7.3');

%% Plot output merging glider, wfp, and fixed depth data
irminger_presentation_plotting

%% Compare all measurements in the mixed layer
ml_oxygen_compare

%% Earlier scripts no longer used here but with pieces to use in finalizing pipeline
% cruise_oxygen
% wfp_deepisotherms
% WFP_Pc_evalplots - checks pressure correction by looking at O2 uniformity over deeply convecting layers, as in Wolf et al. 2018
