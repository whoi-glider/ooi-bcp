%% Make test plots of wfp oxygen gain-corrected data

% %Create a merged data set for oxygen
% wfpmerge.time = wfpgrid{1}.time_start(wfpgrid{1}.ind_pair);
% wfpmerge.pdens = wfpgrid{1}.pdens;
% wfpmerge.T = wfpgrid{1}.T;
% wfpmerge.S = wfpgrid{1}.S;
% wfpmerge.O2conc = wfpgrid{1}.O2conc;
% wfpmerge.oxygen_gaincorr = wfpgrid{1}.oxygen_gaincorr;
% for i = 2:7
%     wfpmerge.time = [wfpmerge.time; wfpgrid{i}.time_start(wfpgrid{i}.ind_pair)];
%     wfpmerge.pdens = [wfpmerge.pdens wfpgrid{i}.pdens];
%     wfpmerge.T = [wfpmerge.T wfpgrid{i}.T];
%     wfpmerge.S = [wfpmerge.S wfpgrid{i}.S];
%     wfpmerge.O2conc = [wfpmerge.O2conc wfpgrid{i}.O2conc];
%     if i < 6
%         wfpmerge.oxygen_gaincorr = [wfpmerge.oxygen_gaincorr wfpgrid{i}.oxygen_gaincorr];
%     end
%     if i >= 6
%         wfpmerge.oxygen_gaincorr = [wfpmerge.oxygen_gaincorr wfpgrid{i}.O2conc];
%     end
% end
% 
% %Create a merged data set for the fluorometer
% wfpmerge_flord.time = wfpgrid_flord{1}.time_start;
% wfpmerge_flord.pdens = wfpgrid_flord{1}.pdens;
% wfpmerge_flord.T = wfpgrid_flord{1}.T;
% wfpmerge_flord.S = wfpgrid_flord{1}.S;
% wfpmerge_flord.backscatter = wfpgrid_flord{1}.backscatter;
% wfpmerge_flord.chla = wfpgrid_flord{1}.chla;
% wfpmerge_flord.backscatter_spikes = wfpgrid_flord{1}.backscatter_spikes;
% wfpmerge_flord.chla_spikes = wfpgrid_flord{1}.chla_spikes;
% for i = 2:6
%     wfpmerge_flord.time = [wfpmerge_flord.time; wfpgrid_flord{i}.time_start];
%     wfpmerge_flord.pdens = [wfpmerge_flord.pdens wfpgrid_flord{i}.pdens];
%     wfpmerge_flord.T = [wfpmerge_flord.T wfpgrid_flord{i}.T];
%     wfpmerge_flord.S = [wfpmerge_flord.S wfpgrid_flord{i}.S];
%     wfpmerge_flord.backscatter = [wfpmerge_flord.backscatter wfpgrid_flord{i}.backscatter];
%     wfpmerge_flord.chla = [wfpmerge_flord.chla wfpgrid_flord{i}.chla];
%     wfpmerge_flord.backscatter_spikes = [wfpmerge_flord.backscatter_spikes wfpgrid_flord{i}.backscatter_spikes];
%     wfpmerge_flord.chla_spikes = [wfpmerge_flord.chla_spikes wfpgrid_flord{i}.chla_spikes];
% end

%% Create a merged dataset from the gridded lag-corrected O2 data

wggmerge.time = wgg{1}.time_start;
wggmerge.duration = wgg{1}.duration;
wggmerge.lat = wgg{1}.lat_profile;
wggmerge.lon = wgg{1}.lon_profile;
wggmerge.temp = wgg{1}.temp_grid;
wggmerge.pracsal = wgg{1}.pracsal_grid;
wggmerge.SA = wgg{1}.SA_grid;
wggmerge.CT = wgg{1}.CT_grid;
wggmerge.pdens = wgg{1}.pdens_grid;
wggmerge.doxy_lagcorr = wgg{1}.doxy_lagcorr_grid;

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
end

%% Initial look at lag-corrected, gridded data

figure(4); clf
imagesc(wggmerge.doxy_lagcorr); caxis([240 300]); colorbar

%% Load MLD data from Izi
%Note that these are outdated - need to switch to newer MLD data
addpath('C:/Users/palevsky/Dropbox/OOI Irminger Sea/Files from Izi')
load OOI_HYPM_ML_Feb2021.mat

%% Scatter plot with subset of data
profilerng = [1:5:3371];
sz = 1;
C = cmocean('Dense'); %set colormap

doxy_scat = wggmerge.doxy_lagcorr(:,profilerng);
[X,Y] = meshgrid(wggmerge.time(profilerng), depth_grid);

figure(5); clf
scatter(X(:),Y(:),5,doxy_scat(:),'filled'); hold on;
plot(dt,MLD,'k.','markersize',5); hold on;
axis ij; axis tight
colormap(C); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Oxygen concentration, lag corrected (\mumol/kg)', 'Fontsize', 12)


%% Make contour plots
%Note that contouring doesn't make optimal use of all data points

% %Adjustable parameters for plotting
%     mindepth = 150; maxdepth = 2600;
%     cints = 60; %number of contour intervals
%     C = cmocean('Dense'); %set colormap
%     C2 = cmocean('Algae'); 
% 
% %Make plotting grid
% [X,Y] = meshgrid(wggmerge.time, depth_grid);
% %[X2,Y2] = meshgrid(wfpmerge_flord.time, depth_grid);

% figure(1); clf;
%     subplot(311) %Density
% cmin = 27.5; cmax = 27.8; %manually set min and max
%     cvec = [cmin:(cmax-cmin)/cints:cmax];
% contourf(X,Y,wggmerge.pdens - 1000,cvec,'linecolor','none'); hold on;
% axis([min(wggmerge.time) max(wggmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
% colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
% datetick('x',2,'keeplimits');
% title('\sigma_\theta', 'Fontsize', 12)
% 
%     subplot(312) %Temperature
% cmin = 2; cmax = 6; %manually set min and max
%     cvec = [cmin:(cmax-cmin)/cints:cmax];
% contourf(X,Y,wggmerge.temp,cvec,'linecolor','none'); hold on;
% axis([min(wggmerge.time) max(wggmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
% colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
% datetick('x',2,'keeplimits');
% title('Temperature (deg C)', 'Fontsize', 12)

% subplot(313) %Oxygen_corr concentration
% cmin = 220; cmax = 320; %manually set min and max
%     cvec = [cmin:(cmax-cmin)/cints:cmax];
% contourf(X,Y,wggmerge.O2conc,cvec,'linecolor','none'); hold on;
% %plot(dt,MLD,'k.','markersize',5); hold on;
% axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth - 200]); caxis([cmin cmax]);
% colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
% datetick('x',2,'keeplimits');
% title('Oxygen concentration, uncorrected (\mumol/kg)', 'Fontsize', 12)

%subplot(313) %Oxygen_corr concentration
% cmin = 220; cmax = 320; %manually set min and max
%     cvec = [cmin:(cmax-cmin)/cints:cmax];
% contourf(X,Y,wggmerge.doxy_lagcorr,cvec,'linecolor','none'); hold on;
% %plot(dt,MLD,'k.','markersize',5); hold on;
% axis([min(wggmerge.time) max(wggmerge.time) mindepth maxdepth - 200]); caxis([cmin cmax]);
% colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
% datetick('x',2,'keeplimits');
% title('Oxygen concentration, lag corrected (\mumol/kg)', 'Fontsize', 12)


% figure(2); clf
%     %subplot(212) %Chla
% cmin = 0; cmax = 0.2; %manually set min and max
%     cvec = [cmin:(cmax-cmin)/cints:cmax];
% contourf(X2,Y2,wfpmerge_flord.chla,cvec,'linecolor','none'); hold on;
% plot(dt,MLD,'k.','markersize',5); hold on;
% axis([min(wfpmerge_flord.time) max(wfpmerge_flord.time) mindepth maxdepth - 800]); caxis([cmin cmax]);
% colormap(C2); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
% datetick('x',2,'keeplimits');
% title('Chlorophyll concentration (\mug/L)', 'Fontsize', 12)
% %%
% figure(3); clf
%     C3 = [gray(256/4); cmocean('Matter')];
% imagesc(wfpmerge_flord.backscatter_spikes)
% colormap(C3); colorbar
% caxis([0 2E-3])