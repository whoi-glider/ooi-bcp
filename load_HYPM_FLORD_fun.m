function [wfp] = load_HYPM_FLORD_fun(filename_FLORD, filterres, max_tol)
%Function to load OOI profiler data (based on load_HYPM_Yr5)
% INPUTS:
%   filename_FLORD: name of the netcdf file downloaded from OOI Data Portal
%   (make sure this file is in the path prior to calling function)
%   filterres: resolution (number of points) used for spike filtering
%   max_tol: tolerance for backscatter outliers (remove all values above this)
% OUTPUTS:
%   (all outputs are structures with multiple variables)
%   wfp - extracted data from filename_FLORD, with backscatter outliers
%   filtered and Briggs 2011 minmax spike filtering applied at filterres

%Load FLORD data
   wfp.time_flord = ncread(filename_FLORD,'time');
   wfp.lon_flord = ncread(filename_FLORD,'lon');
   wfp.lat_flord = ncread(filename_FLORD,'lat');
   try
       wfp.temperature_flord = ncread(filename_FLORD,'sea_water_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
   catch
       wfp.temperature_flord = ncread(filename_FLORD,'ctdpf_ckl_seawater_temperature'); 
   end
   try
       wfp.pracsal_flord = ncread(filename_FLORD,'sea_water_practical_salinity'); %standard_name = 'sea_water_practical_salinity'
   catch
       wfp.pracsal_flord = ncread(filename_FLORD,'practical_salinity'); %standard_name = 'sea_water_practical_salinity'
   end
   wfp.pressure_flord = ncread(filename_FLORD,'int_ctd_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
     %Fluorometer data
   wfp.backscatter = ncread(filename_FLORD,'optical_backscatter'); %long_name = 'Optical Backscatter' units = 'm-1'
   wfp.scat_total = ncread(filename_FLORD,'seawater_scattering_coefficient'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
   wfp.chla = ncread(filename_FLORD,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'
    %Convert to matlab time
   wfp.time_flord_mat = convertTime(wfp.time_flord);
    
%% Assign profile indices
wfp.depth_flord = -gsw_z_from_p(wfp.pressure_flord,wfp.lat_flord);
[wfp.profile_index,wfp.updown_index] = profileIndex(wfp.depth_flord);

%% Find and index outlier values for backscatter
wfp.backscatter_w_outliers = wfp.backscatter; %create variable to preserve outliers
outindex = find (wfp.backscatter >= max_tol); %find all points greater that outlier tolerance
outprofileindex = unique (wfp.profile_index(outindex)); %identify all profiles that include any outliers
for i = 1:length(outprofileindex)
    ind = find (wfp.profile_index == outprofileindex (i));
    wfp.backscatter(ind) = NaN; %set all values to NaN in profiles containing outliers
end

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

end

