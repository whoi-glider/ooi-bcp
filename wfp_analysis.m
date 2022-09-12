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
    'deployment0007_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20200825T200205-20210819T060646.nc'...
    'deployment0008_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20210820T000204-20220710T045503.nc'};

addpath('C:/Users/palevsky/Dropbox/MATLAB/OOI data processing/OOI_Irminger_students/common')
for i = 1:8
    [wfp{i}, wfpgrid{i}, wfpgrid_therm{i}] = load_HYPM_DOSTA_fun(filenames{i}, depth_grid, therm_grid);
end

%% Lag correction
load lagyr1to8.mat %output from wfp_lag.mat

%% Initial look at lag-corrected data

profilerng = [1:327];
yr = 6;

xscat = wgg{yr}.mtime(profilerng,:);
yscat = wgg{yr}.pres(profilerng,:);
doxy_scat = wgg{yr}.doxy_lagcorr(profilerng,:);

figure(200); clf
scatter(xscat(:),yscat(:),[],doxy_scat(:),'filled')
datetick('x')
axis ij
colorbar

%Note that there is indeed good data to be salvaged that doesn't show up in
%paired profiles!

%% Calculate depth intervals of data

figure(201); clf
for yr = 1:8
    A = abs(diff(wgg{yr}.pres'));
    histogram(A(:)); hold on;
    a(yr) = nanmean(A(:));
    b(yr) = nanstd(A(:));
end

%% Check lag-corr output for outliers

%Check for outliers/range test
figure(101); clf
for yr = 1:8
    A = wgg{yr}.doxy_lagcorr(:);
    histogram(A); hold on;
    r(yr,1) = nanmean(A) - 4*nanstd(A);
    r(yr,2) = nanmean(A) + 4*nanstd(A);
    r(yr,3) = length(find(A < r(yr,1)));
    r(yr,4) = length(find(A > r(yr,2)));
end
xlim([220 320])
legend('Year 1','Year 2','Year 3', 'Year 4', 'Year 5', 'Year 6', 'Year 7','Year 8')
title('Histogram all OOI Irminger WFP L2-oxygen, lag corrected only')


% Check for spikess
figure(100); clf
subplot(311)
for yr = 1:8
    A = diff(wgg{yr}.doxy_lagcorr');
    histogram((A(:))); hold on;
end
xlim([-4 4])
legend('Year 1','Year 2','Year 3', 'Year 4', 'Year 5', 'Year 6', 'Year 7','Year 8')
title('Histogram of paired sample difference, all OOI Irminger WFP L2-oxygen, lag corrected only')

subplot(312)
for yr = 1:8
    A = diff(wgg{yr}.temp');
    histogram((A(:))); hold on;
end
xlim([-0.05 0.05])
legend('Year 1','Year 2','Year 3', 'Year 4', 'Year 5', 'Year 6', 'Year 7','Year 8')
title('Histogram of paired sample difference, all OOI Irminger temperature')

subplot(313)
for yr = 1:8
    A = diff(wgg{yr}.pracsal');
    histogram((A(:))); hold on;
end
xlim([-0.01 0.01])
legend('Year 1','Year 2','Year 3', 'Year 4', 'Year 5', 'Year 6', 'Year 7')
title('Histogram of paired sample difference, all OOI Irminger salinity')

%% Flag outlier/spike samples

%tolerances for spikes chosen based on histograms
oxyspike = 3;
tempspike = 0.05;
salspike = 0.005;
c = 0; %counter to keep track of # datapoints flagged

for yr = 1:7

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
    %c
end
 

%% Regrid data reformatted for lag correction

%Select depth resolution and smoothing - current setting is 1 m resolution
%w/ 5-m smoothing
depth_grid = [150:1:2600];
S = 5; %points to smooth over

for yr = 1:7

    %Number of profile indices
    num_profiles = length(wgg{yr}.updown);

    wgg{yr}.time_start = NaN*ones(num_profiles,1);
    wgg{yr}.duration = NaN*ones(num_profiles,1);
    wgg{yr}.lat_profile = NaN*ones(num_profiles,1);
    wgg{yr}.lon_profile = NaN*ones(num_profiles,1);
    wgg{yr}.doxy_lagcorr_grid = NaN*ones(length(depth_grid),num_profiles);
    wgg{yr}.SA_grid = NaN*ones(length(depth_grid),num_profiles);
    wgg{yr}.CT_grid = NaN*ones(length(depth_grid),num_profiles);
    wgg{yr}.pracsal_grid = NaN*ones(length(depth_grid),num_profiles);
    wgg{yr}.temp_grid = NaN*ones(length(depth_grid),num_profiles);
    wgg{yr}.pdens_grid = NaN*ones(length(depth_grid),num_profiles);

    for i = 1:num_profiles
        ind = find(~isnan(wgg{yr}.depth(i,:)) & ~isnan(wgg{yr}.doxy_lagcorr(i,:)) & wgg{yr}.flag(i,:) == 0); %no nan values for depth or oxygen and no range or spike flags
        wgg{yr}.doxy_lagcorr_grid(:,i) = movmean(interp1(wgg{yr}.depth(i,ind), wgg{yr}.doxy_lagcorr(i,ind), depth_grid),S);
        wgg{yr}.SA_grid(:,i) = movmean(interp1(wgg{yr}.depth(i,ind), wgg{yr}.SA(i,ind), depth_grid),S);
        wgg{yr}.CT_grid(:,i) = movmean(interp1(wgg{yr}.depth(i,ind), wgg{yr}.CT(i,ind), depth_grid),S);
        wgg{yr}.pracsal_grid(:,i) = movmean(interp1(wgg{yr}.depth(i,ind), wgg{yr}.pracsal(i,ind), depth_grid),S);
        wgg{yr}.temp_grid(:,i) = movmean(interp1(wgg{yr}.depth(i,ind), wgg{yr}.temp(i,ind), depth_grid),S);
        wgg{yr}.pdens_grid(:,i) = movmean(interp1(wgg{yr}.depth(i,ind), wgg{yr}.pdens(i,ind), depth_grid),S);
        wgg{yr}.time_start(i) = nanmin(wgg{yr}.mtime(i,:));
        wgg{yr}.duration(i) = nanmax(wgg{yr}.mtime(i,:)) - nanmin(wgg{yr}.mtime(i,:));
        wgg{yr}.lat_profile(i) = nanmean(wgg{yr}.lat(i,:));
        wgg{yr}.lon_profile(i) = nanmean(wgg{yr}.lon(i,:)); 
    end
    
end

%%
mergetoplot = [wgg{1}.doxy_lagcorr_grid wgg{2}.doxy_lagcorr_grid wgg{3}.doxy_lagcorr_grid wgg{4}.doxy_lagcorr_grid...
    wgg{5}.doxy_lagcorr_grid wgg{6}.doxy_lagcorr_grid wgg{7}.doxy_lagcorr_grid];

figure(3); clf
imagesc(mergetoplot); caxis([240 300]); colorbar


%% Read in fluorometer data
% filenames_flord = {'deployment0001_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20140912T000208-20150812T103930.nc',...
%     'deployment0002_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20150817T030206-20160628T060527.nc',...
%     'deployment0003_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20160712T000207-20170712T072809.nc',...
%     'deployment0004_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20170807T000204-20180615T185737.nc',...
%     'deployment0005_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20180610T000207-20190630T071452.nc',...
%     'deployment0006_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20190807T000205-20200509T172910.nc'};
% 
% filterres = 5; %choose number of points to use in running min/max filter for backscatter spike analysis
% for i = 1:6
%     [wfp_flord{i}, wfpgrid_flord{i}] = load_HYPM_FLORD_fun(filenames_flord{i}, depth_grid, filterres);
% end

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