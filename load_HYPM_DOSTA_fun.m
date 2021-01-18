function [wfp, wfpgrid, wfpgrid_therm] = load_HYPM_DOSTA_fun(filename_DOSTA, depth_grid, therm_grid)
%Function to load OOI profiler data (based on load_HYPM_Yr5)
% INPUTS:
%   filename_DOSTA: name of the netcdf file downloaded from OOI Data Portal
%   (make sure this file is in the path prior to calling function)
%   depth_grid: depths (m) on which to grid output
        % Example: depth_grid = [150:5:2600];
%   therm_grid: isotherms (C) on which to grid output
        % Example: therm_grid = [1.1:0.05:5.6];
% OUTPUTS:
%   (all outputs are structures with multiple variables)
%   wfp - extracted data from filename_DOSTA
%   wfpgrid - gridded on depth_grid
%   wfpgrid_therm - gridded on therm_grid

%Load DOSTA data
   wfp.time_dosta = ncread(filename_DOSTA,'time');
   wfp.lon_dosta = ncread(filename_DOSTA,'lon');
   wfp.lat_dosta = ncread(filename_DOSTA,'lat');
   wfp.temperature_dosta = ncread(filename_DOSTA,'ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
   wfp.pracsal_dosta = ncread(filename_DOSTA,'practical_salinity'); %standard_name = 'sea_water_practical_salinity'
   wfp.pressure_dosta = ncread(filename_DOSTA,'int_ctd_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
       %Optode data
   wfp.oxygen = ncread(filename_DOSTA,'dissolved_oxygen');%standard_name = 'moles_of_oxygen_per_unit_mass_in_sea_water' units = 'umol kg-1'
   wfp.optode_temperature = ncread(filename_DOSTA,'optode_temperature'); %long_name = 'Optode Temperature' units = 'deg_C'
    %Convert to matlab time
   wfp.time_dosta_mat = convertTime(wfp.time_dosta);
    
%% Assign profile indices prior to gridding
wfp.depth_dosta = -gsw_z_from_p(wfp.pressure_dosta,wfp.lat_dosta);
[wfp.profile_index,wfp.updown_index] = profileIndex(wfp.depth_dosta);

%% Calculate density in raw profiles to enable gridding on density surfaces
[wfp.SA_dosta, in_ocean] = gsw_SA_from_SP(wfp.pracsal_dosta, wfp.pressure_dosta, wfp.lon_dosta, wfp.lat_dosta); %absolute salinity from practical salinity - [SA, ~] = gsw_SA_from_SP(SP,p,long,lat)
wfp.CT_dosta = gsw_CT_from_t(wfp.SA_dosta, wfp.temperature_dosta, wfp.pressure_dosta); %Conservative Temperature from in-situ temperature - CT = gsw_CT_from_t(SA,t,p)
wfp.pdens = gsw_rho(wfp.SA_dosta, wfp.CT_dosta, 0); %calculate potential density at reference pressure of 0 (surface)

%% Grid data to consistent depth intervals/isotherms for each profile
secinday = 60*60*24;

%Grid on depth intervals
scivars = [wfp.temperature_dosta, wfp.pracsal_dosta, wfp.oxygen, wfp.optode_temperature];
[wfpgrid] = glider_grid(wfp.time_dosta, wfp.lat_dosta, wfp.lon_dosta, wfp.depth_dosta, wfp.profile_index, wfp.updown_index', scivars, depth_grid);
    wfpgrid.depth_grid = depth_grid;
wfpgrid.time_start = convertTime(wfpgrid.time_start);
wfpgrid.duration = wfpgrid.duration/secinday;
wfpgrid.updown = wfpgrid.profile_direction;
% Grid on isotherms
    [~,deepind] = unique(wfp.temperature_dosta); %remove duplicates so can use interp1 function
[wfpgrid_therm] = glider_grid_dens(wfp.time_dosta(deepind), wfp.lat_dosta(deepind), wfp.lon_dosta(deepind),...
    wfp.temperature_dosta(deepind), wfp.profile_index(deepind), wfp.updown_index(deepind)', [scivars(deepind,:), wfp.depth_dosta(deepind)], therm_grid);
    wfpgrid_therm.therm_grid = therm_grid;
wfpgrid_therm.time_start = convertTime(wfpgrid_therm.time_start);
wfpgrid_therm.duration = wfpgrid_therm.duration/secinday;
wfpgrid_therm.updown = wfpgrid_therm.profile_direction;    

%% Take mean of paired up and down profiles
    tol = 1; %only combine profiles where time_start is < 1 day apart
[wfpgrid.scivars_pair, wfpgrid.ind_pair] = profilePairMean(wfpgrid,tol);
[wfpgrid_therm.scivars_pair, wfpgrid_therm.ind_pair] = profilePairMean(wfpgrid_therm,tol);

%% Unpack scivars in gridded form

wfpgrid.T = squeeze(wfpgrid.scivars_pair(:,1,:));
wfpgrid.S = squeeze(wfpgrid.scivars_pair(:,2,:));
wfpgrid.O2conc = squeeze(wfpgrid.scivars_pair(:,3,:));
wfpgrid.optode_temperature = squeeze(wfpgrid.scivars_pair(:,4,:));
wfpgrid.pdens = gsw_sigma0(wfpgrid.S,wfpgrid.T)+1000; 
wfpgrid.press = gsw_p_from_z(repmat(-wfpgrid.depth_grid,length(wfpgrid.profile_ind),1)',...
        repmat(wfpgrid.lat,1,length(wfpgrid.depth_grid))');
    
wfpgrid_therm.T = squeeze(wfpgrid_therm.scivars_pair(:,1,:));
wfpgrid_therm.S = squeeze(wfpgrid_therm.scivars_pair(:,2,:));
wfpgrid_therm.O2conc = squeeze(wfpgrid_therm.scivars_pair(:,3,:));
wfpgrid_therm.optode_temperature = squeeze(wfpgrid_therm.scivars_pair(:,4,:));
wfpgrid_therm.depth = squeeze(wfpgrid_therm.scivars_pair(:,5,:));

%Calculate O2 saturation    
    wfpgrid.O2equil = gsw_O2sol_SP_pt(wfpgrid.S, wfpgrid.T);
    wfpgrid.O2sat = (wfpgrid.O2conc./wfpgrid.O2equil - 1)*100;   

%Calculate O2 saturation for thermgrid 
    wfpgrid_therm.O2equil = gsw_O2sol_SP_pt(wfpgrid_therm.S, wfpgrid_therm.T);
    wfpgrid_therm.O2sat = (wfpgrid_therm.O2conc./wfpgrid_therm.O2equil - 1)*100;

end

