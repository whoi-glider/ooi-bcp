%% Analyzes wire-following profiler (HYPM) data
%Note: completely modified from original version written in Sept. 2018 to
%now have separate scripts for each step modeled on Lucy's thesis pipeline

%% Load data
load_HYPM_Yr5

%% Calculate gain correction based on Winkler data
wfp_Irminger_winklercalibration_Yr5
    
%% Perform deep isotherm drift correction
wfp_deepIsotherm_driftCorrection_Yr5 %note - still needs to be updated to make better correction for drift at depth