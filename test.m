addpath('C:/Users/palevsky/Dropbox/Irminger5/mooring_data')

filename_DOSTA = ['deployment0005_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20180610T000207-20190630T071452.nc'];
depth_grid = [150:5:2600];
therm_grid = [1.1:0.05:5.6];

tic
[wfp_t, wfpgrid_t, wfptherm_t] = load_HYPM_DOSTA_fun(filename_DOSTA, depth_grid, therm_grid);
toc