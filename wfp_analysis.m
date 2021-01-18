%% Analyzes wire-following profiler (HYPM) data
%Note: completely modified from original version written in Sept. 2018 to
%now have separate scripts for each step modeled on Lucy's thesis pipeline

%% Load data
%load_HYPM_Yr5

depth_grid = [150:5:2600];
therm_grid = [1.1:0.05:5.6];

addpath('C:/Users/palevsky/Dropbox/Irminger5/mooring_data')
filename_Yr5 = ['deployment0005_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20180610T000207-20190630T071452.nc'];

[Yr5_wfp, Yr5_wfpgrid, Yr5_wfpgrid_therm] = load_HYPM_DOSTA_fun(filename_Yr5, depth_grid, therm_grid);

%% Calculate gain correction based on Winkler data
wfp_Irminger_winklercalibration_Yr5
    
%% Perform deep isotherm drift correction
wfp_deepIsotherm_driftCorrection_Yr5 %note - still needs to be updated to make better correction for drift at depth