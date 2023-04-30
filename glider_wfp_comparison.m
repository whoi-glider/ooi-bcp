%% Run after wfp_analysis and glider_analysis

%% Plot glider locations
M1 = 8; M2 = 5;

figure(2); clf
for i = 1:length(glgmerge)
    indnonan = find((isnan(glgmerge{i}.lon_profile) + isnan(glgmerge{i}.lat_profile) + isnan(glgmerge{i}.time_start)) == 0);
    plot(glgmerge{i}.lon_profile(indnonan), glgmerge{i}.lat_profile(indnonan), '-'); hold on;
end
plot(HYPMlon, HYPMlat, 'ko','markerfacecolor','k','markersize',M1);
for yr = [5,9]
    castsumyr = castsum{yr};
    for i = 1:length(castsumyr)
        if height(castsumyr{i}) > 0
            plot(castsumyr{i}.lon, castsumyr{i}.lat, 'ko','markerfacecolor',nicecolor('wkk'),'markersize',M2);
        end
    end
end
legend([glidertitles,'WFP mooring sites','Turn-around cruise casts'],'location','SW')
axis([-40.1 -38.9 59.5 60.1])

%% Find WFP matchups with glider profiles

dist_tol = 4;
time_tol = 1;

%For each glider profile within dist_tol of HYPM, find nearest time-aligned
%HYPM profile (but only keep if within time_tol)
for i = 1:length(glgmerge)
    glgmerge{i}.HYPMdist_align = NaN(length(glgmerge{i}.profilelist),1);
    glgmerge{i}.HYPMdist_align_ind = NaN(length(glgmerge{i}.profilelist),1);
    glgmerge{i}.HYPMdist_align_time = NaN(length(glgmerge{i}.profilelist),1);
    for j = 1:length(glgmerge{i}.profilelist)
        ind_talign = find(abs(glgmerge{i}.time_start(j) - wggmerge.time) < time_tol);
        if length(ind_talign) > 0
            glgmerge{i}.HYPMdist_align(j) = distlatlon(glgmerge{i}.lat_profile(j), nanmean(wggmerge.lat(ind_talign)), glgmerge{i}.lon_profile(j), nanmean(wggmerge.lon(ind_talign)));
            if glgmerge{i}.HYPMdist_align(j) < dist_tol
                [glgmerge{i}.HYPMdist_align_time(j), indt] = min(abs(glgmerge{i}.time_start(j) - wggmerge.time(ind_talign)));
                glgmerge{i}.HYPMdist_align_ind(j) = ind_talign(indt);
            end
        end
    end
end

%% For glider profiles with WFP matchups, regrid on isotherms
%Then find matching depths and isotherms, and match up
pt_grid_glider = [1.5:0.02:11];
[pt_grid_overlap, ind_glider_ptgrid, ind_hypm_ptgrid] = intersect(pt_grid_glider, pt_grid, 'stable');
    pres_grid_hypm = [150:1:2600];
    pres_grid_glider = [1:1:1000];
[pres_grid_overlap, ind_glider_presgrid, ind_hypm_presgrid] = intersect(pres_grid_glider, pres_grid_hypm, 'stable');

for i = 1:length(glgmerge)
    ind = find(isnan(glgmerge{i}.HYPMdist_align_ind) == 0); %indices for glider profiles that have a HYPM matchup
    %Create arrays to hold oxygen and salinity regridded on isotherms
    glgmerge{i}.HYPMalign_doxy_lagcorr_pt = NaN(length(pt_grid_glider),length(ind));
    glgmerge{i}.HYPMalign_sal_pt = NaN(length(pt_grid_glider),length(ind));
    %Create table to hold stats of alignment comparisons
    vars = {'O2_presA_mean','O2_presA_median','O2_presA_std','T_presA_mean','T_presA_median','T_presA_std','S_presA_mean','S_presA_median','S_presA_std',...
        'O2_thermA_mean','O2_thermA_median','O2_thermA_std','S_thermA_mean','S_thermA_median','S_thermA_std','O2_presA_deepcor_mean','O2_presA_deepcor_median',...
        'O2_presA_deepcor_std','O2_thermA_deepcor_mean','O2_thermA_deepcor_median','O2_thermA_deepcor_std'};
    glgmerge{i}.HYPMalign_stats = array2table(nan(length(ind),length(vars)));
    glgmerge{i}.HYPMalign_stats.Properties.VariableNames = vars;
  
    for j = 1:length(ind)
        ind_nonan = find(isnan(glgmerge{i}.doxy_lagcorr_grid(:,ind(j))) + isnan(glgmerge{i}.sal_grid(:,ind(j))) == 0);
        if length(ind_nonan) > 10
            %Regrid each profile on isotherms
            try
                glgmerge{i}.HYPMalign_doxy_lagcorr_pt(:,j) = interp1(glgmerge{i}.temp_grid(ind_nonan,ind(j)), glgmerge{i}.doxy_lagcorr_grid(ind_nonan,ind(j)), pt_grid_glider);
                glgmerge{i}.HYPMalign_sal_pt(:,j) = interp1(glgmerge{i}.temp_grid(ind_nonan,ind(j)), glgmerge{i}.sal_grid(ind_nonan,ind(j)), pt_grid_glider);
            end
            %Calculate HYPM-glider align factors for:
            % Lag-corrected O2, pracsal, and potential T for each depth
            O2depthComp = glgmerge{i}.doxy_lagcorr_grid(ind_glider_presgrid,ind(j))./wggmerge.doxy_lagcorr(ind_hypm_presgrid,glgmerge{i}.HYPMdist_align_ind(ind(j)));
            O2depthComp_wfpdeepcorr = glgmerge{i}.doxy_lagcorr_grid(ind_glider_presgrid,ind(j))./...
                (wggmerge.doxy_lagcorr(ind_hypm_presgrid,glgmerge{i}.HYPMdist_align_ind(ind(j))).*wfp_profilegain_interp(glgmerge{i}.HYPMdist_align_ind(ind(j))));
            TdepthComp = glgmerge{i}.temp_grid(ind_glider_presgrid,ind(j))./wggmerge.temp(ind_hypm_presgrid,glgmerge{i}.HYPMdist_align_ind(ind(j)));
            SdepthComp = glgmerge{i}.sal_grid(ind_glider_presgrid,ind(j))./wggmerge.pracsal(ind_hypm_presgrid,glgmerge{i}.HYPMdist_align_ind(ind(j)));
            % Lag-corrected O2 and pracsal for each isotherm
            O2thermComp = glgmerge{i}.HYPMalign_doxy_lagcorr_pt(ind_glider_ptgrid,j)./wggmerge.doxy_lagcorr_pt(ind_hypm_ptgrid,glgmerge{i}.HYPMdist_align_ind(ind(j)));
            O2thermComp_wfpdeepcorr = glgmerge{i}.HYPMalign_doxy_lagcorr_pt(ind_glider_ptgrid,j)./...
                (wggmerge.doxy_lagcorr_pt(ind_hypm_ptgrid,glgmerge{i}.HYPMdist_align_ind(ind(j))).*wfp_profilegain_interp(glgmerge{i}.HYPMdist_align_ind(ind(j))));
            SthermComp = glgmerge{i}.HYPMalign_sal_pt(ind_glider_ptgrid,j)./wggmerge.pracsal_pt(ind_hypm_ptgrid,glgmerge{i}.HYPMdist_align_ind(ind(j)));
            % Mean, median, and standard deviation for each of those - put in table
            glgmerge{i}.HYPMalign_stats.O2_presA_mean(j) = nanmean(O2depthComp);
            glgmerge{i}.HYPMalign_stats.O2_presA_median(j) = nanmedian(O2depthComp);
            glgmerge{i}.HYPMalign_stats.O2_presA_std(j) = nanstd(O2depthComp);
            glgmerge{i}.HYPMalign_stats.T_presA_mean(j) = nanmean(TdepthComp);
            glgmerge{i}.HYPMalign_stats.T_presA_median(j) = nanmedian(TdepthComp);
            glgmerge{i}.HYPMalign_stats.T_presA_std(j) = nanstd(TdepthComp);
            glgmerge{i}.HYPMalign_stats.S_presA_mean(j) = nanmean(TdepthComp);
            glgmerge{i}.HYPMalign_stats.S_presA_median(j) = nanmedian(TdepthComp);
            glgmerge{i}.HYPMalign_stats.S_presA_std(j) = nanstd(TdepthComp);
            glgmerge{i}.HYPMalign_stats.O2_thermA_mean(j) = nanmean(O2thermComp);
            glgmerge{i}.HYPMalign_stats.O2_thermA_median(j) = nanmedian(O2thermComp);
            glgmerge{i}.HYPMalign_stats.O2_thermA_std(j) = nanstd(O2thermComp);
            glgmerge{i}.HYPMalign_stats.S_thermA_mean(j) = nanmean(SthermComp);
            glgmerge{i}.HYPMalign_stats.S_thermA_median(j) = nanmedian(SthermComp);
            glgmerge{i}.HYPMalign_stats.S_thermA_std(j) = nanstd(SthermComp);
            glgmerge{i}.HYPMalign_stats.O2_presA_deepcor_mean(j) = nanmean(O2depthComp_wfpdeepcorr);
            glgmerge{i}.HYPMalign_stats.O2_presA_deepcor_median(j) = nanmedian(O2depthComp_wfpdeepcorr);
            glgmerge{i}.HYPMalign_stats.O2_presA_deepcor_std(j) = nanstd(O2depthComp_wfpdeepcorr);
            glgmerge{i}.HYPMalign_stats.O2_thermA_deepcor_mean(j) = nanmean(O2thermComp_wfpdeepcorr);
            glgmerge{i}.HYPMalign_stats.O2_thermA_deepcor_median(j) = nanmedian(O2thermComp_wfpdeepcorr);
            glgmerge{i}.HYPMalign_stats.O2_thermA_deepcor_std(j) = nanstd(O2thermComp_wfpdeepcorr);
        end
    end
end

%% Plot histograms of alignments

figure(3); clf
subplot(231)
for i = 1:length(glgmerge)
    histogram(glgmerge{i}.HYPMalign_stats.O2_thermA_mean,[0.9:0.01:1.2]); hold on;
end
%legend(glidertitles,'location','NE')
xlabel('Glider: HYPM oxygen ratio in matchup profiles, aligned on isotherms');

subplot(232)
for i = 1:length(glgmerge)
    histogram(glgmerge{i}.HYPMalign_stats.O2_presA_mean,[0.5:0.01:1.]); hold on;
end
%legend(glidertitles,'location','NE')
xlabel('Glider: HYPM oxygen ratio in matchup profiles, aligned on pressure surfaces');

subplot(234)
for i = 1:length(glgmerge)
    ind = find(glgmerge{i}.HYPMalign_stats.O2_thermA_mean > 0.9 & glgmerge{i}.HYPMalign_stats.O2_thermA_mean < 1.2);
    histogram(glgmerge{i}.HYPMalign_stats.O2_thermA_std(ind),[0:0.002:0.05]); hold on;
end
%legend(glidertitles,'location','NE')
xlabel('Stdev of Glider: HYPM matchup profiles, aligned on isotherms');

subplot(235)
for i = 1:length(glgmerge)
    ind = find(glgmerge{i}.HYPMalign_stats.O2_presA_mean > 0.9 & glgmerge{i}.HYPMalign_stats.O2_presA_mean < 1.2);
    histogram(glgmerge{i}.HYPMalign_stats.O2_presA_std(ind),[0:0.002:0.05]); hold on;
end
%legend(glidertitles,'location','NE')
xlabel('Stdev of Glider: HYPM matchup profiles, aligned on pressure surfaces');

subplot(233)
for i = 1:length(glgmerge)
    histogram(glgmerge{i}.HYPMalign_stats.T_presA_mean,[0.9:0.01:1.2]); hold on;
end
%legend(glidertitles,'location','NE')
xlabel('Glider: HYPM temperature ratio in matchup profiles, aligned on pressure surfaces');

subplot(236)
for i = 1:length(glgmerge)
    ind = find(glgmerge{i}.HYPMalign_stats.T_presA_mean > 0.9 & glgmerge{i}.HYPMalign_stats.T_presA_mean < 1.2);
    histogram(glgmerge{i}.HYPMalign_stats.T_presA_std(ind),[0:0.002:0.05]); hold on;
end
legend([{glidertitles{5} glidertitles{8} glidertitles{1} glidertitles{2} glidertitles{3} ...
    glidertitles{4} glidertitles{6} glidertitles{7}}],'location','NE')
xlabel('Stdev of Glider: HYPM temp matchup profiles, aligned on pressure surfaces');

%% Assess "good-quality" matchups based on
% 1) Stdev of oxygen matchup ratio for profile
O2std_cutoff = 0.01;
% 2) Stdev of temperature matchup ratio for profile
Tstd_cutoff = 0.015;
% 3) Absolute value of temperature matchup ratio for profile
Tabs_cutoff = 0.05;

for i = 1:length(glgmerge)
    glgmerge{i}.HYPMalign_stats.flag = zeros(height(glgmerge{i}.HYPMalign_stats),1);
    indO = find(glgmerge{i}.HYPMalign_stats.O2_presA_std > O2std_cutoff);
    glgmerge{i}.HYPMalign_stats.flag(indO) = glgmerge{i}.HYPMalign_stats.flag(indO) + 1;
    indT = find(glgmerge{i}.HYPMalign_stats.T_presA_std > Tstd_cutoff);
    glgmerge{i}.HYPMalign_stats.flag(indT) = glgmerge{i}.HYPMalign_stats.flag(indT) + 10;
    indTabs = find(abs(glgmerge{i}.HYPMalign_stats.T_presA_mean - 1) > Tabs_cutoff);
    glgmerge{i}.HYPMalign_stats.flag(indTabs) = glgmerge{i}.HYPMalign_stats.flag(indTabs) + 100;
end

%% Plot time series of alignments
C_gl = cmocean('phase',18);
slope_pick = 0.25;
glg_reorder = [1:8];
C_gl = C_gl(glg_reorder,:);

figure(4); clf
for i = glg_reorder
    indnoflag = find(glgmerge{i}.HYPMalign_stats.flag == 0 & isnan(glgmerge{i}.HYPMalign_stats.O2_presA_deepcor_mean) == 0);
    if length(indnoflag) > 0
        indlist = glgmerge{i}.HYPMdist_align_ind(~isnan(glgmerge{i}.HYPMdist_align_ind));
        plot(wggmerge.time(indlist), 1./glgmerge{i}.HYPMalign_stats.O2_presA_deepcor_mean,'.k','markersize',2); hold on;
        plot(wggmerge.time(indlist(indnoflag)), 1./glgmerge{i}.HYPMalign_stats.O2_presA_deepcor_mean(indnoflag),'ko','markerfacecolor',C_gl(i,:),'markersize',3); hold on;
        t_start = wggmerge.time(indlist(indnoflag(1)));
        [P(i,:),Sfit{i}] = polyfit(wggmerge.time(indlist(indnoflag)) - t_start, 1./glgmerge{i}.HYPMalign_stats.O2_presA_deepcor_mean(indnoflag),1);
        [y_fit,delta] = polyval(P(i,:),wggmerge.time(indlist(indnoflag)) - t_start, Sfit{i});
        try
            air_corr_slopeset = (glgmerge{i}.Taircal.air_meas_dist(:,10)-slope_pick.*glgmerge{i}.Taircal.ml_o2sat)./(1-slope_pick);
%             plot(glgmerge{i}.Taircal.ml_daten, movmean(glgmerge{i}.Taircal.met_o2sat./glgmerge{i}.Taircal.air_corr, 60),...
%                 '-','linewidth',3,'color',C_gl(i,:).*[0.9 0.9 0.9]); hold on; %empirical slope
            plot(glgmerge{i}.Taircal.ml_daten, movmean(glgmerge{i}.Taircal.met_o2sat./air_corr_slopeset, 60),...
                '-','linewidth',5,'color',C_gl(i,:).*[0.7 0.7 0.7]); hold on; %slope set by slope pick
        end
        h(i) = plot(wggmerge.time(indlist(indnoflag)), y_fit, '-','linewidth',2.5,'color',C_gl(i,:)); hold on;
        plot(wggmerge.time(indlist(indnoflag)), y_fit - delta, '--','linewidth',1,'color',C_gl(i,:)); hold on;
        plot(wggmerge.time(indlist(indnoflag)), y_fit + delta, '--','linewidth',1,'color',C_gl(i,:)); hold on;
    else
        h(i) = h(i-1);
    end
end
%ylim([0.92 1.15])
xlim([min(wgg{1}.time_start)-20 max(wgg{8}.time_start - 150)])
datetick('x',2,'keeplimits')
legend(h(glg_reorder), glidertitles(glg_reorder),'location','SW')
title('Glider gain corrections: Comparison of deep isotherm (points & linear fits) and air calibration (thick line) approaches');


