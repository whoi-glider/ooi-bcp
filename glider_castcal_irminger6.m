%Plot CTD cast data along with profile data from gliders collected within
%at set time and distance tolerance

%% Find all unique integer values in glider profile_index
G560_profile_list = unique(G560.profile_index(mod(G560.profile_index,1)==0)); 
G525_profile_list = unique(G525.profile_index(mod(G525.profile_index,1)==0)); 

%% Make a summary list of key data for each profile
G560_profile_summary = NaN*ones(length(G560_profile_list), 6);
G560_profile_summary(:,1) = G560_profile_list;

%Loop over all glider profiles
for i = 1:length(G560_profile_list)
    ind_prof = find(G560.profile_index == G560_profile_list(i));
    G560_profile_summary(i,2) = nanmean(G560.daten(ind_prof)); %date
    G560_profile_summary(i,3) = nanmean(G560.lat_interp(ind_prof)); %lat
    G560_profile_summary(i,4) = nanmean(G560.lon_interp(ind_prof)); %lon
    G560_profile_summary(i,5) = isnan(nanmean(G560.oxygen_saturation(ind_prof))); %1 == no O2 data, 0 = has data
    G560_profile_summary(i,6) = max(G560.depth(ind_prof)); %max depth
end

G525_profile_summary = NaN*ones(length(G525_profile_list), 6);
G525_profile_summary(:,1) = G525_profile_list;

%Loop over all glider profiles
for i = 1:length(G525_profile_list)
    ind_prof = find(G525.profile_index == G525_profile_list(i));
    G525_profile_summary(i,2) = nanmean(G525.daten(ind_prof)); %date
    G525_profile_summary(i,3) = nanmean(G525.lat_interp(ind_prof)); %lat
    G525_profile_summary(i,4) = nanmean(G525.lon_interp(ind_prof)); %lon
    G525_profile_summary(i,5) = isnan(nanmean(G525.oxygen_saturation(ind_prof))); % + nanmean(G525.conductivity(ind_prof))); %1 == no O2 data, 0 = has data
    G525_profile_summary(i,6) = max(G525.depth(ind_prof)); %max depth
end

%% Determine when glider has oxygen data and is producing profiles > chosen tolerance
    tol_profiledepth = 600; %only use glider profiles that go > 400 m
ind_oxy_G525 = find(G525_profile_summary(:,5) == 0 & G525_profile_summary(:,6) > tol_profiledepth);
ind_oxy_G560 = find(G560_profile_summary(:,5) == 0 & G560_profile_summary(:,6) > tol_profiledepth);

%% Find aligned glider profiles for each CTD cast near gliders
timetol = 0.5;
disttol = 4;
for i = 1:length(alignedcasts)
    ctd_time = timealigned(i);
    ctd_lat = castmeta_irminger6.lat(alignedind(i));
    ctd_lon = -castmeta_irminger6.lon(alignedind(i));
    %Find glider profiles within +/- timetol of the CTD cast
    G525_timealign = find(abs(G525_profile_summary(ind_oxy_G525,2) - ctd_time) < timetol);
    G560_timealign = find(abs(G560_profile_summary(ind_oxy_G560,2) - ctd_time) < timetol);
    %%%%%%% Find time-aligned glider profiles within disttol of the CTD cast for G525
    clear distance_calc distalign
    distance_calc = NaN; distalign = NaN;
    for j = 1:length(G525_timealign)
        distance_calc(j) = distlatlon(G525_profile_summary(ind_oxy_G525(G525_timealign(j)),3), ctd_lat,...
            G525_profile_summary(ind_oxy_G525(G525_timealign(j)),4), ctd_lon);
    end
    distalign = find(distance_calc < disttol);
    cast_glider(alignedcasts(i)).G525_profile_summary = NaN*ones(length(distalign), 7);
    %Create profile summary for aligned profiles
    cast_glider(alignedcasts(i)).G525_profile_summary(:,1:5) = G525_profile_summary((ind_oxy_G525(G525_timealign(distalign))),[1:4,6]);
    cast_glider(alignedcasts(i)).G525_profile_summary(:,6) = cast_glider(alignedcasts(i)).G525_profile_summary(:,2) - ctd_time;
    cast_glider(alignedcasts(i)).G525_profile_summary(:,7) = distance_calc(distalign);
    %%%%%%% Find time-aligned glider profiles within disttol of the CTD cast for G560
    clear distance_calc distalign
    distance_calc = NaN; distalign = NaN;
    for j = 1:length(G560_timealign)
        distance_calc(j) = distlatlon(G560_profile_summary(ind_oxy_G560(G560_timealign(j)),3), ctd_lat,...
            G560_profile_summary(ind_oxy_G560(G560_timealign(j)),4), ctd_lon);
    end
    distalign = find(distance_calc < disttol);
    cast_glider(alignedcasts(i)).G560_profile_summary = NaN*ones(length(distalign), 7);
    %Create profile summary for aligned profiles
    cast_glider(alignedcasts(i)).G560_profile_summary(:,1:5) = G560_profile_summary((ind_oxy_G560(G560_timealign(distalign))),[1:4,6]);
    cast_glider(alignedcasts(i)).G560_profile_summary(:,6) = cast_glider(alignedcasts(i)).G560_profile_summary(:,2) - ctd_time;
    cast_glider(alignedcasts(i)).G560_profile_summary(:,7) = distance_calc(distalign);
end

%% Plot aligned glider and CTD cast data with associated Winklers
% Note that correction based on saturation works better than concentration
% because of missing CTD data on GL525 at time of glider reballasting
% calibration cast (C008)

%Calculate O2 saturation for Winklers
Winkler6.O2sol = gsw_O2sol_SP_pt(Winkler6.S, Winkler6.T);

for i = 1:length(alignedcasts)
    figure; clf;
hd = plot((cast{alignedcasts(i)}.O2corr./cast{alignedcasts(i)}.O2sol)*100, cast{alignedcasts(i)}.D, 'k.'); hold on;
hu = plot(cast{alignedcasts(i)}.O2corr(cast{alignedcasts(i)}.maxindex:end)./cast{alignedcasts(i)}.O2sol(cast{alignedcasts(i)}.maxindex:end)*100,...
    cast{alignedcasts(i)}.D(cast{alignedcasts(i)}.maxindex:end), 'r.'); hold on;
    indWink = find(Winkler6.cast == alignedcasts(i));
hw = plot(Winkler6.O2_bcp(intersect(indWink, ind_bcp))./Winkler6.O2sol(intersect(indWink, ind_bcp))*100, Winkler6.depth(intersect(indWink, ind_bcp)), 'm.','markersize',15); hold on;
for j = 1:length(cast_glider(alignedcasts(i)).G525_profile_summary(:,1))
    ind_prof = find(G525.profile_index == cast_glider(alignedcasts(i)).G525_profile_summary(j,1));
    h525 = plot(G525.O2sat_corr(ind_prof), G525.depth_interp(ind_prof), 'b.'); hold on;
end
for j = 1:length(cast_glider(alignedcasts(i)).G560_profile_summary(:,1))
    ind_prof = find(G560.profile_index == cast_glider(alignedcasts(i)).G560_profile_summary(j,1));
    h560 = plot(G560.O2sat_corr(ind_prof), G560.depth_interp(ind_prof), 'c.'); hold on;
end
axis ij
ylim([0 1000])
xlabel('Oxygen % saturation')
ylabel('Depth (m)')
legend([hd, hu, hw, h525, h560],'CTD Downcast','CTD Upcast','Winkler O_2','GL525','GL560','location','southeast')
title(['Irminger6 glider calibration, Cast ' num2str(alignedcasts(i))])
end


