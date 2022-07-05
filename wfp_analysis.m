%% Analyzes wire-following profiler (HYPM) data
%Note: completely modified from original version written in Sept. 2018 to
%now have separate scripts for each step modeled on Lucy's thesis pipeline

%% Load oxygen data from all deployments and interpolate onto even grid
    depth_grid = [150:5:2600];
    therm_grid = [1.1:0.05:5.6];

%https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered/catalog.html
addpath('C:/Users/palevsky/Dropbox/OOI Irminger Sea/OOI_downloads/THREDDS_updated/HYPM')
filenames = {'deployment0001_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20140912T000208-20150812T103930.nc',...
    'deployment0002_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20150817T030206-20160628T060527.nc',...
    'deployment0003_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20160712T000207-20170712T072809.nc',...
    'deployment0004_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20170807T000204-20180615T185737.nc',...
    'deployment0005_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20180610T000207-20190630T071452.nc',...
    'deployment0006_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20190807T000205-20200509T172910.nc'...
    'deployment0007_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20200825T200205-20210819T060646.nc'};

addpath('C:/Users/palevsky/Dropbox/MATLAB/OOI data processing/OOI_Irminger_students/common')
for i = 1:7
    [wfp{i}, wfpgrid{i}, wfpgrid_therm{i}] = load_HYPM_DOSTA_fun(filenames{i}, depth_grid, therm_grid);
end

%% Read in fluorometer data
filenames_flord = {'deployment0001_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20140912T000208-20150812T103930.nc',...
    'deployment0002_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20150817T030206-20160628T060527.nc',...
    'deployment0003_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20160712T000207-20170712T072809.nc',...
    'deployment0004_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20170807T000204-20180615T185737.nc',...
    'deployment0005_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20180610T000207-20190630T071452.nc',...
    'deployment0006_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20190807T000205-20200509T172910.nc'};

filterres = 5; %choose number of points to use in running min/max filter for backscatter spike analysis
for i = 1:6
    [wfp_flord{i}, wfpgrid_flord{i}] = load_HYPM_FLORD_fun(filenames_flord{i}, depth_grid, filterres);
end

%% Load Winkler data for calibrations
addpath('C:/Users/palevsky/Dropbox/Wellesley/OOI_Irminger_students/CruiseData_Yrs1to4')
loadWinklerIrmingerYrs1to5

%% Calculate Year 1-5 gain corrections based on Winkler data
wfp_Irminger_winklercalibration_Yrs1to5

%% Apply initial gain corrections to Year 1-5 data

for i = 1:5
    wfp{i}.oxygen_gaincorr = wfp{i}.oxygen * gain_hypm(i);
    wfpgrid{i}.oxygen_gaincorr = wfpgrid{i}.O2conc * gain_hypm(i);
    wfpgrid_therm{i}.oxygen_gaincorr = wfpgrid_therm{i}.O2conc * gain_hypm(i);
end

%% Test plots to ensure that all is calculated appropriately
wfp_plotting

%% Perform deep isotherm drift correction
%wfp_deepIsotherm_driftCorrection_Yr5 %note - still needs to be updated to make better correction for drift at depth