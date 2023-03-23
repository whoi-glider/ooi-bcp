%% Create a merged dataset from the gridded lag-corrected O2 data

%Overarching data
wggmerge.time = wgg{1}.time_start;
wggmerge.duration = wgg{1}.duration;
wggmerge.lat = wgg{1}.lat_profile;
wggmerge.lon = wgg{1}.lon_profile;
wggmerge.deploy_yr = ones(length(wgg{1}.time_start),1);

%Data gridded on pressure surfaces
wggmerge.temp = wgg{1}.temp_grid;
wggmerge.pracsal = wgg{1}.pracsal_grid;
wggmerge.SA = wgg{1}.SA_grid;
wggmerge.CT = wgg{1}.CT_grid;
wggmerge.pdens = wgg{1}.pdens_grid;
wggmerge.doxy_lagcorr = wgg{1}.doxy_lagcorr_grid;

%Data gridded on isotherms (potential temperature)
wggmerge.doxy_lagcorr_pt = wgg{1}.doxy_lagcorr_ptgrid;
wggmerge.SA_pt = wgg{1}.SA_ptgrid;
wggmerge.pracsal_pt = wgg{1}.pracsal_ptgrid;
wggmerge.pres_pt = wgg{1}.pres_ptgrid;
wggmerge.temp_pt = wgg{1}.temp_ptgrid;
wggmerge.pdens_pt = wgg{1}.pdens_ptgrid;

for yr = 2:8
    wggmerge.time = [wggmerge.time; wgg{yr}.time_start];
    wggmerge.duration = [wggmerge.duration; wgg{yr}.duration];
    wggmerge.lat = [wggmerge.lat; wgg{yr}.lat_profile];
    wggmerge.lon = [wggmerge.lon; wgg{yr}.lon_profile];
    wggmerge.temp = [wggmerge.temp wgg{yr}.temp_grid];
    wggmerge.pracsal = [wggmerge.pracsal wgg{yr}.pracsal_grid];
    wggmerge.SA = [wggmerge.SA wgg{yr}.SA_grid];
    wggmerge.CT = [wggmerge.CT wgg{yr}.CT_grid];
    wggmerge.pdens = [wggmerge.pdens wgg{yr}.pdens_grid];
    wggmerge.doxy_lagcorr = [wggmerge.doxy_lagcorr wgg{yr}.doxy_lagcorr_grid];
    wggmerge.deploy_yr = [wggmerge.deploy_yr; yr*ones(length(wgg{yr}.time_start),1)];
    wggmerge.doxy_lagcorr_pt = [wggmerge.doxy_lagcorr_pt wgg{yr}.doxy_lagcorr_ptgrid];
    wggmerge.SA_pt = [wggmerge.SA_pt wgg{yr}.SA_ptgrid];
    wggmerge.pracsal_pt = [wggmerge.pracsal_pt wgg{yr}.pracsal_ptgrid];
    wggmerge.pres_pt = [wggmerge.pres_pt wgg{yr}.pres_ptgrid];
    wggmerge.temp_pt = [wggmerge.temp_pt wgg{yr}.temp_ptgrid];
    wggmerge.pdens_pt = [wggmerge.pdens_pt wgg{yr}.pdens_ptgrid];
end

%% Create a merged dataset from the gridded fluorometer data

%Overarching data
wggmerge_fl.time = wgg_flord{1}.time_start;
wggmerge_fl.duration = wgg_flord{1}.duration;
wggmerge_fl.lat = wgg_flord{1}.lat_profile;
wggmerge_fl.lon = wgg_flord{1}.lon_profile;
wggmerge_fl.deploy_yr = ones(length(wgg_flord{1}.time_start),1);

%Data gridded on pressure surfaces
wggmerge_fl.temp = wgg_flord{1}.temp_grid;
wggmerge_fl.pracsal = wgg_flord{1}.pracsal_grid;
wggmerge_fl.chla = wgg_flord{1}.chla_grid;
wggmerge_fl.backscatter = wgg_flord{1}.backscatter_grid;
wggmerge_fl.spikes = wgg_flord{1}.spikes_grid;
wggmerge_fl.chlspikes = wgg_flord{1}.chlspikes_grid;

for yr = 2:7
    wggmerge_fl.time = [wggmerge_fl.time; wgg_flord{yr}.time_start];
    wggmerge_fl.duration = [wggmerge_fl.duration; wgg_flord{yr}.duration];
    wggmerge_fl.lat = [wggmerge_fl.lat; wgg_flord{yr}.lat_profile];
    wggmerge_fl.lon = [wggmerge_fl.lon; wgg_flord{yr}.lon_profile];
    wggmerge_fl.temp = [wggmerge_fl.temp wgg_flord{yr}.temp_grid];
    wggmerge_fl.pracsal = [wggmerge_fl.pracsal wgg_flord{yr}.pracsal_grid];
    wggmerge_fl.deploy_yr = [wggmerge_fl.deploy_yr; yr*ones(length(wgg_flord{yr}.time_start),1)];
    wggmerge_fl.chla = [wggmerge_fl.chla wgg_flord{yr}.chla_grid];
    wggmerge_fl.backscatter = [wggmerge_fl.backscatter wgg_flord{yr}.backscatter_grid];
    wggmerge_fl.spikes = [wggmerge_fl.spikes wgg_flord{yr}.spikes_grid];
    wggmerge_fl.chlspikes = [wggmerge_fl.chlspikes wgg_flord{yr}.chlspikes_grid];
end

%% Initial look at lag-corrected, gridded data

figure(4); clf
imagesc(wggmerge.doxy_lagcorr); caxis([240 300]); colorbar
%imagesc(wggmerge.doxy_lagcorr_pt); caxis([240 300]); colorbar

%% Load MLD data from Izi
%Note that these are outdated - need to switch to newer MLD data
addpath('C:/Users/palevsky/Dropbox/OOI Irminger Sea/Files from Izi')
load OOI_HYPM_ML_Feb2021.mat

%% Scatter plot with subset of data
profilerng = [1:5:3371];
sz = 1;
C = cmocean('Dense'); %set colormap

doxy_scat = wggmerge.doxy_lagcorr(:,profilerng);
[X,Y] = meshgrid(wggmerge.time(profilerng), pres_grid_hypm);

figure(5); clf
scatter(X(:),Y(:),5,doxy_scat(:),'filled'); hold on;
plot(dt,MLD,'k.','markersize',5); hold on;
axis ij; axis tight
colormap(C); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Oxygen concentration, lag corrected (\mumol/kg)', 'Fontsize', 12)

%% Scatter plot on isotherms

[X,Y] = meshgrid(wggmerge.time, pt_grid);

figure(6); clf
    subplot(211)
scatter(X(:),Y(:),5, wggmerge.doxy_lagcorr_pt(:),'filled'); hold on;
axis tight
colormap(C); ylabel('Potential temperature', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Oxygen concentration, lag corrected (\mumol/kg)', 'Fontsize', 12)

	subplot(212)
scatter(X(:),Y(:),5, wggmerge.temp_pt(:),'filled'); hold on;
axis tight
colormap(C); ylabel('Potential temperature', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('In situ temperature', 'Fontsize', 12)


