%% Year 5 profiler data
%Note that this is a continuation of the Yrs 1-4 data processed by Lucy
%Wanzer in load_HYPM

%Wire-following profiler, Year 5, DOSTA
addpath('C:/Users/palevsky/Dropbox/Irminger5/mooring_data')
filename = ['deployment0005_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20180610T000207-20190630T071452.nc']; ncdisp(filename)
    Yr5_wfp.time_dosta = ncread(filename,'time');
    Yr5_wfp.lon_dosta = ncread(filename,'lon');
    Yr5_wfp.lat_dosta = ncread(filename,'lat');
    Yr5_wfp.temperature_dosta = ncread(filename,'ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr5_wfp.pracsal_dosta = ncread(filename,'practical_salinity'); %standard_name = 'sea_water_practical_salinity'
    Yr5_wfp.pressure_dosta = ncread(filename,'int_ctd_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
       %Optode data
    Yr5_wfp.oxygen = ncread(filename,'dissolved_oxygen');%standard_name = 'moles_of_oxygen_per_unit_mass_in_sea_water' units = 'umol kg-1'
    Yr5_wfp.optode_temperature = ncread(filename,'optode_temperature'); %long_name = 'Optode Temperature' units = 'deg_C'
    %Convert to matlab time
    Yr5_wfp.time_dosta_mat = convertTime(Yr5_wfp.time_dosta);
    
%Wire-following profiler, Year 5, Fluorometer    
filename = ['deployment0005_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20180610T000207-20190630T071452.nc']; ncdisp(filename)
    Yr5_wfp.time_flord = ncread(filename,'time');
    Yr5_wfp.lon_flord = ncread(filename,'lon');
    Yr5_wfp.lat_flord = ncread(filename,'lat');
    Yr5_wfp.temperature_flord = ncread(filename,'ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr5_wfp.pracsal_flord = ncread(filename,'practical_salinity'); %standard_name = 'sea_water_practical_salinity'
    Yr5_wfp.pressure_flord = ncread(filename,'int_ctd_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Fluorometer data
    Yr5_wfp.backscatter = ncread(filename,'optical_backscatter'); %long_name = 'Optical Backscatter' units = 'm-1'
    Yr5_wfp.scat_total = ncread(filename,'seawater_scattering_coefficient'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
    Yr5_wfp.chla = ncread(filename,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'
    %Convert to matlab time
    Yr5_wfp.time_flord_mat = convertTime(Yr5_wfp.time_flord);
    
%% Assign profile indices prior to gridding
Yr5_wfp.depth_dosta = -gsw_z_from_p(Yr5_wfp.pressure_dosta,Yr5_wfp.lat_dosta);
    [Yr5_wfp.profile_index,Yr5_wfp.updown_index] = profileIndex(Yr5_wfp.depth_dosta);

%% Calculate density in raw profiles to enable gridding on density surfaces
[Yr5_wfp.SA_dosta, in_ocean] = gsw_SA_from_SP(Yr5_wfp.pracsal_dosta, Yr5_wfp.pressure_dosta, Yr5_wfp.lon_dosta, Yr5_wfp.lat_dosta); %absolute salinity from practical salinity - [SA, ~] = gsw_SA_from_SP(SP,p,long,lat)
Yr5_wfp.CT_dosta = gsw_CT_from_t(Yr5_wfp.SA_dosta, Yr5_wfp.temperature_dosta, Yr5_wfp.pressure_dosta); %Conservative Temperature from in-situ temperature - CT = gsw_CT_from_t(SA,t,p)
Yr5_wfp.pdens = gsw_rho(Yr5_wfp.SA_dosta, Yr5_wfp.CT_dosta, 0); %calculate potential density at reference pressure of 0 (surface)

%% Grid data to consistent depth intervals for each profile
depth_grid = [150:5:2600];
therm_grid = [1.1:0.05:5.6];
secinday = 60*60*24;

%All profiles for year 5
scivars = [Yr5_wfp.temperature_dosta, Yr5_wfp.pracsal_dosta, Yr5_wfp.oxygen, Yr5_wfp.optode_temperature...
        Yr5_wfp.backscatter, Yr5_wfp.scat_total, Yr5_wfp.chla];
[Yr5_wfpgrid] = glider_grid(Yr5_wfp.time_dosta,Yr5_wfp.lat_dosta,Yr5_wfp.lon_dosta,Yr5_wfp.depth_dosta,Yr5_wfp.profile_index,Yr5_wfp.updown_index',scivars,depth_grid);
    Yr5_wfpgrid.depth_grid = depth_grid;
Yr5_wfpgrid.time_start = convertTime(Yr5_wfpgrid.time_start);
Yr5_wfpgrid.duration = Yr5_wfpgrid.duration/secinday;
Yr5_wfpgrid.updown = Yr5_wfpgrid.profile_direction;
% Grid on isotherms
    [~,deepind] = unique(Yr5_wfp.temperature_dosta); %remove duplicates so can use interp1 function
[Yr5_wfpgrid_therm] = glider_grid_dens(Yr5_wfp.time_dosta(deepind),Yr5_wfp.lat_dosta(deepind),Yr5_wfp.lon_dosta(deepind),...
    Yr5_wfp.temperature_dosta(deepind),Yr5_wfp.profile_index(deepind),Yr5_wfp.updown_index(deepind)',[scivars(deepind,:), Yr5_wfp.depth_dosta(deepind)],therm_grid);
    Yr5_wfpgrid_therm.therm_grid = therm_grid;
Yr5_wfpgrid_therm.time_start = convertTime(Yr5_wfpgrid_therm.time_start);
Yr5_wfpgrid_therm.duration = Yr5_wfpgrid_therm.duration/secinday;
Yr5_wfpgrid_therm.updown = Yr5_wfpgrid_therm.profile_direction;    

%% Take mean of paired up and down profiles
    tol = 1; %only combine profiles where time_start is < 1 day apart
[Yr5_wfpgrid.scivars_pair,Yr5_wfpgrid.ind_pair] = profilePairMean(Yr5_wfpgrid,tol);
[Yr5_wfpgrid_therm.scivars_pair,Yr5_wfpgrid_therm.ind_pair] = profilePairMean(Yr5_wfpgrid_therm,tol);

%% Unpack scivars in gridded form
%When using scivars, gets all profiles (both up and down)
%When using scivars_pair, takes mean of paired up and down profiles

% Year 5
Yr5_wfpgrid.T = squeeze(Yr5_wfpgrid.scivars_pair(:,1,:));
Yr5_wfpgrid.S = squeeze(Yr5_wfpgrid.scivars_pair(:,2,:));
Yr5_wfpgrid.O2conc = squeeze(Yr5_wfpgrid.scivars_pair(:,3,:));
Yr5_wfpgrid.optode_temperature = squeeze(Yr5_wfpgrid.scivars_pair(:,4,:));
Yr5_wfpgrid.backscatter = squeeze(Yr5_wfpgrid.scivars_pair(:,5,:));
Yr5_wfpgrid.scat_total = squeeze(Yr5_wfpgrid.scivars_pair(:,6,:));
Yr5_wfpgrid.chla = squeeze(Yr5_wfpgrid.scivars_pair(:,7,:));
Yr5_wfpgrid.pdens = gsw_sigma0(Yr5_wfpgrid.S,Yr5_wfpgrid.T)+1000; 
Yr5_wfpgrid.press = gsw_p_from_z(repmat(-Yr5_wfpgrid.depth_grid,length(Yr5_wfpgrid.profile_ind),1)',...
        repmat(Yr5_wfpgrid.lat,1,length(Yr5_wfpgrid.depth_grid))');
    
Yr5_wfpgrid_therm.T = squeeze(Yr5_wfpgrid_therm.scivars_pair(:,1,:));
Yr5_wfpgrid_therm.S = squeeze(Yr5_wfpgrid_therm.scivars_pair(:,2,:));
Yr5_wfpgrid_therm.O2conc = squeeze(Yr5_wfpgrid_therm.scivars_pair(:,3,:));
Yr5_wfpgrid_therm.optode_temperature = squeeze(Yr5_wfpgrid_therm.scivars_pair(:,4,:));
Yr5_wfpgrid_therm.backscatter = squeeze(Yr5_wfpgrid_therm.scivars_pair(:,5,:));
Yr5_wfpgrid_therm.scat_total = squeeze(Yr5_wfpgrid_therm.scivars_pair(:,6,:));
Yr5_wfpgrid_therm.chla = squeeze(Yr5_wfpgrid_therm.scivars_pair(:,7,:));  
Yr5_wfpgrid_therm.depth = squeeze(Yr5_wfpgrid_therm.scivars_pair(:,8,:));

%Calculate O2 saturation    
    Yr5_wfpgrid.O2equil = gsw_O2sol_SP_pt(Yr5_wfpgrid.S,Yr5_wfpgrid.T);
Yr5_wfpgrid.O2sat = (Yr5_wfpgrid.O2conc./Yr5_wfpgrid.O2equil - 1)*100;   

%Calculate O2 saturation for thermgrid 
 Yr5_wfpgrid_therm.O2equil = gsw_O2sol_SP_pt(Yr5_wfpgrid_therm.S,Yr5_wfpgrid_therm.T);
Yr5_wfpgrid_therm.O2sat = (Yr5_wfpgrid_therm.O2conc./Yr5_wfpgrid_therm.O2equil - 1)*100;

    