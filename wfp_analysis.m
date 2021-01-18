%% Analyzes wire-following profiler (HYPM) data
%Note: completely modified from original version written in Sept. 2018 to
%now have separate scripts for each step modeled on Lucy's thesis pipeline

%% Load data from all deployments and interpolate onto even grid
    depth_grid = [150:5:2600];
    therm_grid = [1.1:0.05:5.6];

addpath('C:/Users/palevsky/Dropbox/OOI Irminger Sea/OOI_downloads/THREDDS_updated/HYPM')
filenames = {'deployment0001_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20140912T000208-20150812T103930.nc',...
    'deployment0002_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20150817T030206-20160628T060527.nc',...
    'deployment0003_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20160712T000207-20170712T072809.nc',...
    'deployment0004_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20170807T000204-20180615T185737.nc',...
    'deployment0005_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20180610T000207-20190630T071452.nc',...
    'deployment0006_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20190807T000205-20191231T225556.nc'};

for i = 1:6
    [wfp{i}, wfpgrid{i}, wfpgrid_therm{i}] = load_HYPM_DOSTA_fun(filenames{i}, depth_grid, therm_grid);
end

%% Load Winkler data for calibrations
addpath('C:/Users/palevsky/Dropbox/Wellesley/OOI_Irminger_students/CruiseData_Yrs1to4')
loadWinklerIrmingerYrs1to5

%% Calculate Year 1-5 gain corrections based on Winkler data
wfp_Irminger_winklercalibration_Yrs1to5
    
%% Perform deep isotherm drift correction
wfp_deepIsotherm_driftCorrection_Yr5 %note - still needs to be updated to make better correction for drift at depth