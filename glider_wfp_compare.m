function [glider_profile_summary, wfp_profile_summary] = glider_wfp_compare(glider, mooring, dist_compare, time_compare, updownindex)

%%% Inputs
% glider = table of glider data (i.e. G453 or G363)
% mooring = structure of data (i.e. Yr5_wfp)
% dist_compare = tolerance for radius around mooring from which to take glider profiles for comparison - in km
% time_compare = tolerance for time difference between location-aligned glider and mooring profiles to compare - in days
% updownindex = determine whether to use glider profiles with up or down data

%%% Outputs
% glider_profile_summary = array with all profile numbers, and associated lat, lon, time, and distance from the mooring
% mooring_profile_summary = array with all profile numbers, and associated time, minimum time gap from a location-aligned glider profile, and min and max indices of associated glider profiles

% **Also plots a map of glider profiles over time, along with mooring location and highlights gliders within dist_compare of mooring

%% Find all unique integer values in glider and wfp profile_index
glider_profile_list = unique(glider.profile_index(mod(glider.profile_index,1)==0)); 
wfp_profile_list = unique(mooring.profile_index(mod(mooring.profile_index,1)==0));

%% Make a summary list of key data for each profile
glider_profile_summary = NaN*ones(length(glider_profile_list), 6);
glider_profile_summary(:,1) = glider_profile_list;

%Loop over all glider profiles
for i = 1:length(glider_profile_list)
    ind_prof = find(glider.profile_index == glider_profile_list(i));
    glider_profile_summary(i,2) = nanmean(glider.daten(ind_prof)); %date
    glider_profile_summary(i,3) = nanmean(glider.lat_interp(ind_prof)); %lat
    glider_profile_summary(i,4) = nanmean(glider.lon_interp(ind_prof)); %lon
    glider_profile_summary(i,5) = distlatlon(glider_profile_summary(i,3), nanmean(mooring.lat), glider_profile_summary(i,4), nanmean(mooring.lon)); %distance from mooring in km
    glider_profile_summary(i,6) = nanmean(glider.profile_direction(ind_prof)); %updown index
end

%% Determine when glider is close enough to mooring for comparison and is correct up vs down to have data
ind_compare = find(glider_profile_summary(:,5) < dist_compare & glider_profile_summary(:,6) == updownindex);

%% Find wfp profiles that align in time with glider swim-by times
wfp_profile_summary = NaN*ones(length(wfp_profile_list), 5);
wfp_profile_summary(:,1) = wfp_profile_list;

%Loop over all glider profiles
for i = 1:length(wfp_profile_list)
    ind_prof = find(mooring.profile_index == wfp_profile_list(i));
    wfp_profile_summary(i,2) = nanmean(mooring.time_mat(ind_prof)); %date
    %%% Minimum time gap between this profile and any glider profile in range, along with the index of the closest time-aligned glider profile
    wfp_profile_summary(i,3) = min(abs(wfp_profile_summary(i,2) - glider_profile_summary(ind_compare,2))); %minimum time gap between this profile and any glider profile in range
    %%% Need to instead find all within range rather than just minimum because some profiles don't have any data in them
    time_matches = find(abs(wfp_profile_summary(i,2) - glider_profile_summary(ind_compare,2)) < time_compare); %find profiles with time alignment within the range time_compare
    if length(time_matches) > 0
        wfp_profile_summary(i,4) = min(time_matches);
        wfp_profile_summary(i,5) = max(time_matches);
    end
end

%% Determine when glider and mooring profiles are close enough in time as well as space
ind_aligned = find(isnan(wfp_profile_summary(:,4)) == 0);

%% Plot aligned profiles from glider and mooring
% 1 subplot for each in ind_aligned (might need to be changed once have longer wfp record)
% Start with just oxygen for now, but should be easy to add later
num_aligned = length(ind_aligned);
figure; clf
for i = 1:num_aligned
    subplot(1, num_aligned, i)
        ind_wfp = find(mooring.profile_index == wfp_profile_summary(ind_aligned(i),1));
    plot(mooring.oxygen(ind_wfp), mooring.depth(ind_wfp),'k.'); hold on;
         ind_glider = find(glider.profile_index >= glider_profile_summary(ind_compare(wfp_profile_summary(ind_aligned(i),4)),1) & ...
             glider.profile_index <= glider_profile_summary(ind_compare(wfp_profile_summary(ind_aligned(i),5)),1));
    plot(glider.O2_corr(ind_glider), glider.depth_interp(ind_glider),'m.'); hold on;
    set(gca,'YDir','reverse'); ylabel('Depth (m)');
    axis([270 320 150 1000])
end

%% Map locations of each glider profile over time since deployment, compared to wfp location
figure; clf
scatter(glider_profile_summary(:,4), glider_profile_summary(:,3), [], glider_profile_summary(:,2) - min(glider_profile_summary(:,2)), 'filled'); hold on;
plot(glider_profile_summary(ind_compare,4), glider_profile_summary(ind_compare,3), 'k.'); hold on;
plot(nanmean(mooring.lon), nanmean(mooring.lat), 'ko', 'markersize', 15); hold on;
colorbar;
axis([-40.1 -38.9 59.6 60.05])

end