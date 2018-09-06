%% Analyze profiler data telemetered from current deployment

%Wire-following profiler, Year 5, DOSTA - telemetered data    
filename = ['deployment0005_GI02HYPM-WFP02-03-DOSTAL000-telemetered-dosta_ln_wfp_instrument_20180610T000207-20180706T023027.nc']; ncdisp(filename)
    Yr5_wfp.time = ncread(filename,'time');
    Yr5_wfp.lon = ncread(filename,'lon');
    Yr5_wfp.lat = ncread(filename,'lat');
    %CTD data
    Yr5_wfp.temperature = ncread(filename,'ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr5_wfp.pracsal = ncread(filename,'practical_salinity'); %standard_name = 'sea_water_practical_salinity'
    Yr5_wfp.pressure = ncread(filename,'int_ctd_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Optode data
    Yr5_wfp.oxygen = ncread(filename,'dissolved_oxygen'); %standard_name = 'moles_of_oxygen_per_unit_mass_in_sea_water' units = 'umol kg-1'
    Yr5_wfp.optode_temperature = ncread(filename,'optode_temperature'); %long_name = 'Optode Temperature' units = 'deg_C'
    %Convert to matlab time
    Yr5_wfp.time_mat = convertTime(Yr5_wfp.time);
    
%% Assign profile indices
Yr5_wfp.depth = -gsw_z_from_p(Yr5_wfp.pressure, Yr5_wfp.lat);
    [Yr5_wfp.profile_index, Yr5_wfp.updown_index] = profileIndex(Yr5_wfp.depth);

%% Calculate potential density
Yr5_wfp.pdens = gsw_p_from_z(-Yr5_wfp.depth, Yr5_wfp.lat);

%% Visualize all wfp data without any gridding

figure; clf;
    cmin = 250; cmax = 300; %manually set min and max
    mindepth = 150; maxdepth = 2700;
scatter(Yr5_wfp.time_mat, Yr5_wfp.depth,[], Yr5_wfp.oxygen, 'filled');
    axis([min(Yr5_wfp.time_mat) max(Yr5_wfp.time_mat) mindepth maxdepth]); caxis([cmin cmax]);
    datetick('x'); set(gca,'YDir','reverse'); ylabel('Depth (m)');
    hcb = colorbar; set(hcb,'location','eastoutside')
    datetick('x',2,'keeplimits');
    title('WFP oxygen concentration (raw)')