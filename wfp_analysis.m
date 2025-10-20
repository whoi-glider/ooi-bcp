%% Loads, does lag correction, and removes outliers from wire-following profiler (HYPM) data

%% Load oxygen data from all deployments and interpolate onto even grid
%This step reads in the THREDDS Gold Copy data from all cruises, unpacks
%variables, identifies profile numbers, and then grids paired up/down profiles.

%https://thredds.dataexplorer.oceanobservatories.org/thredds/catalog/ooigoldcopy/public/GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered/catalog.html
addpath('C:/Users/palevsky/Dropbox/OOI Irminger Sea/OOI_downloads/THREDDS_updated/HYPM')
filenames = {'deployment0001_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20140912T000208-20150812T103930.nc',...
    'deployment0002_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20150817T030206-20160628T060527.nc',...
    'deployment0003_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20160712T000207-20170712T072809.nc',...
    'deployment0004_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20170807T000204-20180615T185737.nc',...
    'deployment0005_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20180610T000207-20190630T071452.nc',...
    'deployment0006_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20190807T000205-20200509T172910.nc'...
    'deployment0007_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20200825T200205-20210819T060646.nc'...
    'deployment0008_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20210820T000204-20220710T045503.nc'...
    'deployment0009_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20220702T000208-20230611T104311.nc'...
    'deployment0010_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20230905T000209-20240625T090554.nc'};

for i = 1:10
    [wfp{i}] = load_HYPM_DOSTA_fun(filenames{i});
end

%% Read in fluorometer data and complete initial processing
filenames_flord = {'deployment0001_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20140912T000208-20150812T103930.nc',...
    'deployment0002_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20150817T030206-20160628T060527.nc',...
    'deployment0003_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20160712T000207-20170712T072809.nc',...
    'deployment0004_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20170807T000204-20180615T185737.nc',...
    'deployment0005_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20180610T000207-20190630T071452.nc',...
    'deployment0006_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20190807T000205-20200509T172910.nc'...
    'deployment0007_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20200825T200205-20210819T060646.nc'...
    'deployment0008_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20210814T000210-20220710T045503.nc'...
    'deployment0009_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20220702T000208-20230611T104311.nc' ...
    'deployment0010_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20230905T000209-20240625T090554.nc'};

filterres = 11; %choose number of points to use in running min/max filter for backscatter spike analysis
max_tol = 0.005;
for i = 1:10
    [wfp_flord{i}] = load_HYPM_FLORD_fun(filenames_flord{i}, filterres, max_tol);
end

%% Lag correction
%wfp_lag.m %this takes a long time to run, so output is saved
%load wfp_lag_output_2Sept2023.mat %"wgg" output from wfp_lag.m
load wfp_lag_output_20Oct2025.mat %updated version with years 9 and 10

%% Calculate depth intervals of data

figure; clf
for yr = [10,9,8,7,6,5,4,3,2,1]
    A = abs(diff(wgg{yr}.pres'));
    histogram(A(:)); hold on;
    a(yr) = nanmean(A(:));
    b(yr) = nanstd(A(:));
end
xlim([0 4.5])
legend('Year 10','Year 9','Year 8','Year 7','Year 6', 'Year 5', 'Year 4', 'Year 3', 'Year 2','Year 1')
title('Histogram all depth intervals in raw WFP measurements')
xlabel('Presure (dbar)')

%% Check lag-corr output for outliers

%Check for outliers/range test
figure; clf
for yr = 1:10
    A = wgg{yr}.doxy_lagcorr(:);
    histogram(A); hold on;
    r(yr,1) = nanmean(A) - 4*nanstd(A);
    r(yr,2) = nanmean(A) + 4*nanstd(A);
    r(yr,3) = length(find(A < r(yr,1)));
    r(yr,4) = length(find(A > r(yr,2)));
end
xlim([280 390])
legend('Year 1','Year 2','Year 3', 'Year 4', 'Year 5', 'Year 6', 'Year 7','Year 8','Year 9','Year 10')
title('Histogram all OOI Irminger WFP L1-oxygen, lag corrected only')


% Check for spikess
figure; clf
subplot(311)
for yr = 1:10
    A = diff(wgg{yr}.doxy_lagcorr');
    histogram((A(:))); hold on;
end
xlim([-4 4])
legend('Year 1','Year 2','Year 3', 'Year 4', 'Year 5', 'Year 6', 'Year 7','Year 8', 'Year 9','Year 10')
title('Histogram of paired sample difference, all OOI Irminger WFP L1-oxygen, lag corrected only')

subplot(312)
for yr = 1:10
    A = diff(wgg{yr}.temp');
    histogram((A(:))); hold on;
end
xlim([-0.05 0.05])
legend('Year 1','Year 2','Year 3', 'Year 4', 'Year 5', 'Year 6', 'Year 7','Year 8', 'Year 9','Year 10')
title('Histogram of paired sample difference, all OOI Irminger temperature')

subplot(313)
for yr = 1:10
    A = diff(wgg{yr}.pracsal');
    histogram((A(:))); hold on;
end
xlim([-0.01 0.01])
legend('Year 1','Year 2','Year 3', 'Year 4', 'Year 5', 'Year 6', 'Year 7','Year 8', 'Year 9','Year 10')
title('Histogram of paired sample difference, all OOI Irminger salinity')

%% Flag outlier/spike samples

%tolerances for spikes chosen based on histograms
oxyspike = 3;
tempspike = 0.05;
salspike = 0.005;
c = 0; %counter to keep track of # datapoints flagged

for yr = 1:10

    %Number of profiles
    num_profiles = length(wgg{yr}.updown);

    %Create array to hold flag
    wgg{yr}.flag = 0*ones(size(wgg{yr}.mtime));

    for i = 1:num_profiles
        O = diff(wgg{yr}.doxy_lagcorr(i,:));
        Oout = find(abs(O) > oxyspike);
        T = diff(wgg{yr}.temp(i,:));
        Tout = find(abs(T) > tempspike);
        S = diff(wgg{yr}.pracsal(i,:));
        Sout = find(abs(S) > salspike);
        if yr == 5 & i == 141
            wgg{yr}.flag(i,100:end) = 1; %flag entire section of profile for range issue
        else
            wgg{yr}.flag(i,Oout) = wgg{yr}.flag(i,Oout) + 1;
            wgg{yr}.flag(i,Tout) = wgg{yr}.flag(i,Tout) + 10;
            wgg{yr}.flag(i,Sout) = wgg{yr}.flag(i,Sout) + 100;
        end
        c = c + length(find(wgg{yr}.flag(i,:) > 0));
    end
end
 
%% Reformat fluorometer data following same structure as oxygen data
tol = 50; %only use profiles with at least 50 points

for yr = 1:10
numprofiles = max(wfp_flord{yr}.profile_index);
    wgg_flord{yr}.mtime = NaN(numprofiles,2000);
    wgg_flord{yr}.lon = NaN(numprofiles,2000);
    wgg_flord{yr}.lat = NaN(numprofiles,2000);
    wgg_flord{yr}.pracsal = NaN(numprofiles,2000);
    wgg_flord{yr}.pres = NaN(numprofiles,2000);
    wgg_flord{yr}.temp = NaN(numprofiles,2000);
    wgg_flord{yr}.chla = NaN(numprofiles,2000);
    wgg_flord{yr}.backscatter = NaN(numprofiles,2000);
    wgg_flord{yr}.backscatter_w_outliers = NaN(numprofiles,2000);
    wgg_flord{yr}.filteredspikes = NaN(numprofiles,2000);
    wgg_flord{yr}.filteredchlspikes = NaN(numprofiles,2000);
    wgg_flord{yr}.updown = NaN(numprofiles,1);
for i = 1:numprofiles
    indp = find(wfp{yr}.profile_index == i);
    if length(indp) > tol
        wgg_flord{yr}.mtime(i,1:length(indp)) = wfp_flord{yr}.time_flord_mat(indp);
        wgg_flord{yr}.lon(i,1:length(indp)) = wfp_flord{yr}.lon_flord(indp);
        wgg_flord{yr}.lat(i,1:length(indp)) = wfp_flord{yr}.lat_flord(indp);
        wgg_flord{yr}.pracsal(i,1:length(indp)) = wfp_flord{yr}.pracsal_flord(indp);
        wgg_flord{yr}.pres(i,1:length(indp)) = wfp_flord{yr}.pressure_flord(indp);
        wgg_flord{yr}.temp(i,1:length(indp)) = wfp_flord{yr}.temperature_flord(indp);
        wgg_flord{yr}.chla(i,1:length(indp)) = wfp_flord{yr}.chla(indp);
        wgg_flord{yr}.backscatter(i,1:length(indp)) = wfp_flord{yr}.backscatter(indp);
        wgg_flord{yr}.backscatter_w_outliers(i,1:length(indp)) = wfp_flord{yr}.backscatter_w_outliers(indp);
        wgg_flord{yr}.filteredspikes(i,1:length(indp)) = wfp_flord{yr}.filteredspikes(indp);
        wgg_flord{yr}.filteredchlspikes(i,1:length(indp)) = wfp_flord{yr}.filteredchlspikes(indp);
        wgg_flord{yr}.updown(i) = wfp_flord{yr}.updown_index(indp(1+tol));
    end
    indzero = find(wgg_flord{yr}.backscatter(i,:) == 0);
    wgg_flord{yr}.backscatter(i,indzero) = NaN;
end

% Cut out rows and columns of all NaNs
ind_good = find(~all(isnan(wgg_flord{yr}.backscatter')) & ~isnan(wgg_flord{yr}.updown)' == 1); %includes > tol pts and has up/down index
    wgg_flord{yr}.mtime = wgg_flord{yr}.mtime(ind_good,~all(isnan(wgg_flord{yr}.backscatter)));
    wgg_flord{yr}.lon = wgg_flord{yr}.lon(ind_good,~all(isnan(wgg_flord{yr}.backscatter)));
    wgg_flord{yr}.lat = wgg_flord{yr}.lat(ind_good,~all(isnan(wgg_flord{yr}.backscatter)));
    wgg_flord{yr}.pracsal = wgg_flord{yr}.pracsal(ind_good,~all(isnan(wgg_flord{yr}.backscatter)));
    wgg_flord{yr}.pres = wgg_flord{yr}.pres(ind_good,~all(isnan(wgg_flord{yr}.backscatter)));
    wgg_flord{yr}.temp = wgg_flord{yr}.temp(ind_good,~all(isnan(wgg_flord{yr}.backscatter)));
    wgg_flord{yr}.chla = wgg_flord{yr}.chla(ind_good,~all(isnan(wgg_flord{yr}.backscatter)));
    wgg_flord{yr}.backscatter = wgg_flord{yr}.backscatter(ind_good,~all(isnan(wgg_flord{yr}.backscatter)));
    wgg_flord{yr}.filteredspikes = wgg_flord{yr}.filteredspikes(ind_good,~all(isnan(wgg_flord{yr}.backscatter)));
    wgg_flord{yr}.filteredchlspikes = wgg_flord{yr}.filteredchlspikes(ind_good,~all(isnan(wgg_flord{yr}.backscatter)));
    wgg_flord{yr}.updown = wgg_flord{yr}.updown(ind_good);
    
end

%% Regrid reformatted flord data

%Select depth resolution and smoothing - current setting is 1 m resolution
%w/ 5-m smoothing
pres_grid = [150:1:2600];
S = 5; %points to smooth over

for yr = 1:10

    %Number of profile indices
    num_profiles = length(wgg_flord{yr}.updown);

    wgg_flord{yr}.time_start = NaN*ones(num_profiles,1);
    wgg_flord{yr}.duration = NaN*ones(num_profiles,1);
    wgg_flord{yr}.lat_profile = NaN*ones(num_profiles,1);
    wgg_flord{yr}.lon_profile = NaN*ones(num_profiles,1);
    wgg_flord{yr}.chla_grid = NaN*ones(length(pres_grid),num_profiles);
    wgg_flord{yr}.backscatter_grid = NaN*ones(length(pres_grid),num_profiles);
    wgg_flord{yr}.spikes_grid = NaN*ones(length(pres_grid),num_profiles);
    wgg_flord{yr}.chlspikes_grid = NaN*ones(length(pres_grid),num_profiles);
    wgg_flord{yr}.pracsal_grid = NaN*ones(length(pres_grid),num_profiles);
    wgg_flord{yr}.temp_grid = NaN*ones(length(pres_grid),num_profiles);

    for i = 1:num_profiles
        try
            ind = find(~isnan(wgg_flord{yr}.pres(i,:)) & ~isnan(wgg_flord{yr}.backscatter(i,:))); %no nan values for depth or backscatter
            wgg_flord{yr}.chla_grid(:,i) = movmean(interp1(wgg_flord{yr}.pres(i,ind), wgg_flord{yr}.chla(i,ind), pres_grid),S);
            wgg_flord{yr}.backscatter_grid(:,i) = movmean(interp1(wgg_flord{yr}.pres(i,ind), wgg_flord{yr}.backscatter(i,ind), pres_grid),S);
            wgg_flord{yr}.spikes_grid(:,i) = movmean(interp1(wgg_flord{yr}.pres(i,ind), wgg_flord{yr}.filteredspikes(i,ind), pres_grid),S);
            wgg_flord{yr}.chlspikes_grid(:,i) = movmean(interp1(wgg_flord{yr}.pres(i,ind), wgg_flord{yr}.filteredchlspikes(i,ind), pres_grid),S);
            wgg_flord{yr}.pracsal_grid(:,i) = movmean(interp1(wgg_flord{yr}.pres(i,ind), wgg_flord{yr}.pracsal(i,ind), pres_grid),S);
            wgg_flord{yr}.temp_grid(:,i) = movmean(interp1(wgg_flord{yr}.pres(i,ind), wgg_flord{yr}.temp(i,ind), pres_grid),S);
            wgg_flord{yr}.time_start(i) = nanmin(wgg_flord{yr}.mtime(i,:));
            wgg_flord{yr}.duration(i) = nanmax(wgg_flord{yr}.mtime(i,:)) - nanmin(wgg_flord{yr}.mtime(i,:));
            wgg_flord{yr}.lat_profile(i) = nanmean(wgg_flord{yr}.lat(i,:));
            wgg_flord{yr}.lon_profile(i) = nanmean(wgg_flord{yr}.lon(i,:)); 
        end
    end
    
end

%% Regrid oxygen data reformatted for lag correction

%Select depth resolution and smoothing - current setting is 1 m resolution
%w/ 5-m smoothing
pres_grid = [150:1:2600];
pres_grid_hypm = pres_grid;
S = 5; %points to smooth over

for yr = 1:10

    %Number of profile indices
    num_profiles = length(wgg{yr}.updown);

    wgg{yr}.time_start = NaN*ones(num_profiles,1);
    wgg{yr}.duration = NaN*ones(num_profiles,1);
    wgg{yr}.lat_profile = NaN*ones(num_profiles,1);
    wgg{yr}.lon_profile = NaN*ones(num_profiles,1);
    wgg{yr}.doxy_lagcorr_grid = NaN*ones(length(pres_grid),num_profiles);
    wgg{yr}.SA_grid = NaN*ones(length(pres_grid),num_profiles);
    wgg{yr}.CT_grid = NaN*ones(length(pres_grid),num_profiles);
    wgg{yr}.pracsal_grid = NaN*ones(length(pres_grid),num_profiles);
    wgg{yr}.temp_grid = NaN*ones(length(pres_grid),num_profiles);
    wgg{yr}.pdens_grid = NaN*ones(length(pres_grid),num_profiles);

    for i = 1:num_profiles
        ind = find(~isnan(wgg{yr}.pres(i,:)) & ~isnan(wgg{yr}.doxy_lagcorr(i,:)) & wgg{yr}.flag(i,:) == 0); %no nan values for depth or oxygen and no range or spike flags
        wgg{yr}.doxy_lagcorr_grid(:,i) = movmean(interp1(wgg{yr}.pres(i,ind), wgg{yr}.doxy_lagcorr(i,ind), pres_grid),S);
        wgg{yr}.SA_grid(:,i) = movmean(interp1(wgg{yr}.pres(i,ind), wgg{yr}.SA(i,ind), pres_grid),S);
        wgg{yr}.CT_grid(:,i) = movmean(interp1(wgg{yr}.pres(i,ind), wgg{yr}.CT(i,ind), pres_grid),S);
        wgg{yr}.pracsal_grid(:,i) = movmean(interp1(wgg{yr}.pres(i,ind), wgg{yr}.pracsal(i,ind), pres_grid),S);
        wgg{yr}.temp_grid(:,i) = movmean(interp1(wgg{yr}.pres(i,ind), wgg{yr}.temp(i,ind), pres_grid),S);
        wgg{yr}.pdens_grid(:,i) = movmean(interp1(wgg{yr}.pres(i,ind), wgg{yr}.pdens(i,ind), pres_grid),S);
            ind_pos = find(wgg{yr}.mtime(i,:) > 0);
        wgg{yr}.time_start(i) = nanmin(wgg{yr}.mtime(i,ind_pos));
        wgg{yr}.duration(i) = nanmax(wgg{yr}.mtime(i,:)) - nanmin(wgg{yr}.mtime(i,:));
        wgg{yr}.lat_profile(i) = nanmean(wgg{yr}.lat(i,:));
        wgg{yr}.lon_profile(i) = nanmean(wgg{yr}.lon(i,:)); 
    end
    
end

%% Grid on isotherms

pt_grid = [1.5:0.02:5]; %Izi did 1.5-2.9 at intervals of 0.1, as potential temperature
S = 5; %points to smooth over

for yr = 1:10
    
    %Calculate potential temperature
    wgg{yr}.ptemp = gsw_pt0_from_t(wgg{yr}.SA,wgg{yr}.temp,wgg{yr}.pres);

    %Initialize arrays to hold regridded data
    num_profiles = length(wgg{yr}.updown);
    wgg{yr}.doxy_lagcorr_ptgrid = NaN*ones(length(pt_grid),num_profiles);
    wgg{yr}.SA_ptgrid = NaN*ones(length(pt_grid),num_profiles);
    wgg{yr}.pracsal_ptgrid = NaN*ones(length(pt_grid),num_profiles);
    wgg{yr}.pres_ptgrid = NaN*ones(length(pt_grid),num_profiles);
    wgg{yr}.temp_ptgrid = NaN*ones(length(pt_grid),num_profiles);
    wgg{yr}.pdens_ptgrid = NaN*ones(length(pt_grid),num_profiles);

    for i = 1:num_profiles
        ind = find(~isnan(wgg{yr}.pres(i,:)) & ~isnan(wgg{yr}.doxy_lagcorr(i,:)) & wgg{yr}.flag(i,:) == 0); %no nan values for depth or oxygen and no range or spike flags
        wgg{yr}.doxy_lagcorr_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.doxy_lagcorr(i,ind), pt_grid),S);
        wgg{yr}.SA_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.SA(i,ind), pt_grid),S);
        wgg{yr}.pracsal_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.pracsal(i,ind), pt_grid),S);
        wgg{yr}.pres_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.pres(i,ind), pt_grid),S);
        wgg{yr}.temp_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.temp(i,ind), pt_grid),S);
        wgg{yr}.pdens_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.pdens(i,ind), pt_grid),S);
    end
    
end



