function [wfp, wfpgrid] = load_HYPM_FLORD_fun(filename_FLORD, depth_grid, filterres)
%Function to load OOI profiler data (based on load_HYPM_Yr5)
% INPUTS:
%   filename_FLORD: name of the netcdf file downloaded from OOI Data Portal
%   (make sure this file is in the path prior to calling function)
%   depth_grid: depths (m) on which to grid output
        % Example: depth_grid = [150:5:2600];
%   filterres: resolution (number of points) used for spike filtering
% OUTPUTS:
%   (all outputs are structures with multiple variables)
%   wfp - extracted data from filename_FLORD
%   wfpgrid - gridded on depth_grid

%Load FLORD data
   wfp.time_flord = ncread(filename_FLORD,'time');
   wfp.lon_flord = ncread(filename_FLORD,'lon');
   wfp.lat_flord = ncread(filename_FLORD,'lat');
   wfp.temperature_flord = ncread(filename_FLORD,'ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
   wfp.pracsal_flord = ncread(filename_FLORD,'practical_salinity'); %standard_name = 'sea_water_practical_salinity'
   wfp.pressure_flord = ncread(filename_FLORD,'int_ctd_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
     %Fluorometer data
   wfp.backscatter = ncread(filename_FLORD,'optical_backscatter'); %long_name = 'Optical Backscatter' units = 'm-1'
   wfp.scat_total = ncread(filename_FLORD,'seawater_scattering_coefficient'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
   wfp.chla = ncread(filename_FLORD,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'
    %Convert to matlab time
   wfp.time_flord_mat = convertTime(wfp.time_flord);
    
%% Assign profile indices prior to gridding
wfp.depth_flord = -gsw_z_from_p(wfp.pressure_flord,wfp.lat_flord);
[wfp.profile_index,wfp.updown_index] = profileIndex(wfp.depth_flord);

%% Calculate backscatter spikes - after Briggs et al. 2011, from Jose

for u = 1:wfp.profile_index(end)
   [profileindexindex] = find (wfp.profile_index == u);
   minfilter = movmin (wfp.backscatter(profileindexindex), filterres);
   minmaxfilter = movmax (minfilter, filterres);
   wfp.filteredspikes (profileindexindex) = wfp.backscatter(profileindexindex) - minmaxfilter;
   minfilter = movmin (wfp.chla(profileindexindex), filterres);
   minmaxfilter = movmax (minfilter, filterres);
   wfp.filteredchlspikes (profileindexindex) = wfp.chla(profileindexindex) - minmaxfilter;
end

%% Grid data to consistent depth intervals for each profile
secinday = 60*60*24;

%Grid on depth intervals
scivars = [wfp.temperature_flord, wfp.pracsal_flord, wfp.backscatter, wfp.scat_total, wfp.chla, wfp.filteredspikes', wfp.filteredchlspikes'];
[wfpgrid] = glider_grid(wfp.time_flord, wfp.lat_flord, wfp.lon_flord, wfp.depth_flord, wfp.profile_index, wfp.updown_index', scivars, depth_grid);
    wfpgrid.depth_grid = depth_grid;
wfpgrid.time_start = convertTime(wfpgrid.time_start);
wfpgrid.duration = wfpgrid.duration/secinday;
wfpgrid.updown = wfpgrid.profile_direction;

%% Unpack scivars in gridded form

wfpgrid.T = squeeze(wfpgrid.scivars(:,1,:));
wfpgrid.S = squeeze(wfpgrid.scivars(:,2,:));
wfpgrid.backscatter = squeeze(wfpgrid.scivars(:,3,:));
wfpgrid.scat_total = squeeze(wfpgrid.scivars(:,4,:));
wfpgrid.chla = squeeze(wfpgrid.scivars(:,5,:));
wfpgrid.backscatter_spikes = squeeze(wfpgrid.scivars(:,6,:));
wfpgrid.chla_spikes = squeeze(wfpgrid.scivars(:,7,:));
wfpgrid.pdens = gsw_sigma0(wfpgrid.S,wfpgrid.T)+1000; 
wfpgrid.press = gsw_p_from_z(repmat(-wfpgrid.depth_grid,length(wfpgrid.profile_ind),1)',...
        repmat(wfpgrid.lat,1,length(wfpgrid.depth_grid))');

end

