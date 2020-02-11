%% Post-recovery analysis of glider data from Irminger5 (June 2018 - July 2019)

%% Load data
addpath('C:/Users/palevsky/Dropbox/Irminger5/glider')
load G363.mat %this is the recovered data file (just for 363, only glider recovered)
load I5_final_telemetered.mat %final telemetered data for 363 and 453

%% Add path for functions
addpath('C:/Users/palevsky/Dropbox/MATLAB/OOI data processing/OOI_Irminger')

%% Interpolate GPS and CTD data onto all time points

G363R.lon_interp = naninterp1(G363R.time, G363R.longitude, G363R.time);
G363R.lat_interp = naninterp1(G363R.time, G363R.latitude, G363R.time);
G363R.salinity_interp = naninterp1(G363R.time, G363R.salinity, G363R.time);
G363R.pressure_interp = naninterp1(G363R.time, G363R.pressure, G363R.time);
G363R.temperature_interp = naninterp1(G363R.time, G363R.temperature, G363R.time);

G453.lon_interp = naninterp1(G453.time, G453.longitude, G453.time);
G453.lat_interp = naninterp1(G453.time, G453.latitude, G453.time);
G453.salinity_interp = naninterp1(G453.time, G453.salinity, G453.time);
G453.pressure_interp = naninterp1(G453.time, G453.pressure, G453.time);
G453.temperature_interp = naninterp1(G453.time, G453.temperature, G453.time);

%% Plot full depth data over entire deployments
figure(1); clf
    subplot(221)
plot(G363R.daten, G363R.depth_interp, 'k.'); hold on;
scatter(G363R.daten, G363R.depth_interp, [], G363R.oxygen_saturation,'filled'); colorbar; caxis([85 100])
set(gca,'YDir','reverse'); 
xlim([datenum(2018,6,8,12,0,0) datenum(2019,5,14,0,0,0)])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 363, oxygen saturation')
    subplot(222)
plot(G363R.daten, G363R.depth_interp, 'k.'); hold on;
scatter(G363R.daten, G363R.depth_interp, [], G363R.temperature,'filled'); colorbar; caxis([2.5 6.5])
set(gca,'YDir','reverse'); 
xlim([datenum(2018,6,8,12,0,0) datenum(2019,5,14,0,0,0)])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 363, temperature')
    subplot(223)
plot(G453.daten, G453.depth_interp, 'k.'); hold on;
scatter(G453.daten, G453.depth_interp, [], G453.oxygen_saturation,'filled'); colorbar; caxis([85 100])
set(gca,'YDir','reverse'); 
xlim([datenum(2018,6,8,12,0,0) datenum(2019,5,14,0,0,0)])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 453, oxygen saturation')
    subplot(224)
plot(G453.daten, G453.depth_interp, 'k.'); hold on;
scatter(G453.daten, G453.depth_interp, [], G453.temperature,'filled'); colorbar; caxis([2.5 6.5])
set(gca,'YDir','reverse'); 
xlim([datenum(2018,6,8,12,0,0) datenum(2019,5,14,0,0,0)])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 453, temperature')

%% Overwrite telemetered 363 with recovered data
G363 = G363R;

%% Determine historesis effect between climb and dive profiles
% This is used as basis for determining the lag correction
% Note that this is calculated before any P or S corrections are applied
glider_updown_compare

%% Apply the lag correction throughout the full glider deployment
    secinday = 60*60*24;
    tau = 50/secinday;
    timetol = 600*10;

    indnonan = find(~isnan(G363.oxygen_saturation));
G363.oxygen_saturation_lagcorr = NaN(height(G363),1);
G363.oxygen_concentration_lagcorr = NaN(height(G363),1);
G363.oxygen_saturation_lagcorr(indnonan) = lagCorr(G363.oxygen_saturation(indnonan), G363.time(indnonan), tau, timetol);
G363.oxygen_concentration_lagcorr(indnonan) = lagCorr(G363.oxygen_concentration(indnonan), G363.time(indnonan), tau, timetol);

    indnonan = find(~isnan(G453.oxygen_saturation));
G453.oxygen_saturation_lagcorr = NaN(height(G453),1);
G453.oxygen_concentration_lagcorr = NaN(height(G453),1);
G453.oxygen_saturation_lagcorr(indnonan) = lagCorr(G453.oxygen_saturation(indnonan), G453.time(indnonan), tau, timetol);
G453.oxygen_concentration_lagcorr(indnonan) = lagCorr(G453.oxygen_concentration(indnonan), G453.time(indnonan), tau, timetol);

%% Apply pressure and salinity corrections to lag-corrected oxygen data
%Note that internal salinity setting is 35 for 363 and 0 for 453

G363.O2_corr = aaoptode_salpresscorr(G363.oxygen_concentration_lagcorr, G363.temperature_interp, G363.salinity_interp, G363.pressure_interp, 35);
G363.O2sat_corr = G363.oxygen_saturation_lagcorr.*(1+G363.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

G453.O2_corr = aaoptode_salpresscorr(G453.oxygen_concentration_lagcorr, G453.temperature_interp, G453.salinity_interp, G453.pressure_interp, 0);
G453.O2sat_corr = G453.oxygen_saturation_lagcorr.*(1+G453.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

%% Use Winkler values to calculate initial gain corrections at the time of deployment
glider_winklercal %compare with Winkler calibration cast values

%% Apply air calibration method over entire deployment for both gliders

Yr5_met = 'deployment0005_GI01SUMO-SBD12-06-METBKA000-telemetered-metbk_a_dcl_instrument_20180608T172154.969000-20190809T080328.237000.nc';
rhcorr = 1; mindateplot = datenum(2018,6,1); maxdateplot = datenum(2019,6,1);

[T_363, med_gain_363] = aircalfun(G363, 'Glider 363', -1, Yr5_met, mindateplot, rhcorr, mindateplot, maxdateplot);
[T_453, med_gain_453] = aircalfun(G453, 'Glider 453', 1, Yr5_met, mindateplot, rhcorr, mindateplot, datenum(2019,2,1));

%% Grid data to consistent depth intervals for each profile
depth_grid = [0:5:1000];
secinday = 60*60*24;

%All profiles for year 5
    G363.temperature_interp(G363.temperature_interp == 0) = NaN; %remove zero values
scivars = [G363.temperature_interp, G363.salinity_interp, G363.O2_corr, G363.O2sat_corr...
        G363.backscatter, G363.chlorophyll];
[G363grid] = glider_grid(G363.daten, G363.lat_interp, G363.lon_interp, ...
    G363.depth_interp, G363.profile_index, G363.profile_direction, ...
    scivars,depth_grid);
G363grid.depth_grid = depth_grid;

    G453.temperature_interp(G453.temperature_interp == 0) = NaN; %remove zero values
scivars = [G453.temperature_interp, G453.salinity_interp, G453.O2_corr, G453.O2sat_corr...
        G453.backscatter, G453.chlorophyll];
[G453grid] = glider_grid(G453.daten, G453.lat_interp, G453.lon_interp, ...
    G453.depth_interp, G453.profile_index, G453.profile_direction, ...
    scivars,depth_grid);
G453grid.depth_grid = depth_grid;

%% Plot glider data over full deployment
    tol = 0.25; %only use profiles with at least 50% good data
ind_data_453 = find(sum(~isnan(squeeze(G453grid.scivars(:,4,:)))) > tol*length(depth_grid));
ind_data_363 = find(sum(~isnan(squeeze(G363grid.scivars(:,4,:)))) > tol*length(depth_grid));

%% Plot map of both gliders

%ind = find(Winkler_casts(:,1) >=5 & Winkler_casts(:,1) <=6);
M = 10;
    figure(10); clf
plot(-G453grid.lon(ind_data_453), G453grid.lat(ind_data_453),'k.','markersize',M/2); hold on;
plot(-G363grid.lon(ind_data_363), G363grid.lat(ind_data_363),'w.','markersize',M/2); hold on;
scatter(-G453grid.lon(ind_data_453), G453grid.lat(ind_data_453), M, G453grid.time_start(ind_data_453)); hold on;
scatter(-G363grid.lon(ind_data_363), G363grid.lat(ind_data_363), M, G363grid.time_start(ind_data_363)); hold on;
%plot(Winkler_casts(ind,4), Winkler_casts(ind,3),'r.','markersize',M*2); hold on;

plot(-OOImoorings.SUMO5(2), OOImoorings.SUMO5(1),'^k','markersize',M*2,'markerfacecolor','k'); hold on;
plot(-OOImoorings.HYPM5(2), OOImoorings.HYPM5(1),'^k','markersize',M*2,'markerfacecolor','k'); hold on;
plot(-OOImoorings.FLMA5(2), OOImoorings.FLMA5(1),'^k','markersize',M*2,'markerfacecolor','k'); hold on;
plot(-OOImoorings.FLMB5(2), OOImoorings.FLMB5(1),'^k','markersize',M*2,'markerfacecolor','k'); hold on;

%% Calculate times that each glider is in vicinity of HYPM
tol = 5; %number of kilometers from HYPM

for i = 1:length(ind_data_363)
    dist_HYPM_363(i) = distlatlon(G363grid.lat(ind_data_363(i)), OOImoorings.HYPM5(1), ...
        -G363grid.lon(ind_data_363(i)), -OOImoorings.HYPM5(2));
end
ind_HYPM_363 = ind_data_363(find(dist_HYPM_363 < tol));


for i = 1:length(ind_data_453)
    dist_HYPM_453(i) = distlatlon(G453grid.lat(ind_data_453(i)), OOImoorings.HYPM5(1), ...
        -G453grid.lon(ind_data_453(i)), -OOImoorings.HYPM5(2));
end
ind_HYPM_453 = ind_data_453(find(dist_HYPM_453 < tol));

%% Depth profiles over time
%Adjustable parameters for plotting
    mindepth = 0; maxdepth = 1000;
    cints = 60; %number of contour intervals
    C = cmocean('Dense'); %set colormap
    C2 = cmocean('Algae'); 
    
%Plot for 453
figure(1); clf
[X,Y] = meshgrid(G453grid.time_start(ind_data_453), G453grid.depth_grid);
cmin = 85; cmax = 110; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,squeeze(G453grid.scivars(:,4,ind_data_453)*med_gain_453),cvec,'linecolor','none'); hold on;
axis([min(G453grid.time_start(ind_data_453)) max(G453grid.time_start(ind_data_453)) mindepth maxdepth]); caxis([cmin cmax]);
colormap(cmocean('Balance','pivot',100)); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
plot(G453grid.time_start(ind_HYPM_453), 5*ones(size(ind_HYPM_453)), '.','markersize',15,'color',nicecolor('k'));
datetick('x',2,'keeplimits');
title('Glider 453 oxygen saturation', 'Fontsize', 12)

%Plot for 363
figure(2); clf
[X,Y] = meshgrid(G363grid.time_start(ind_data_363), G363grid.depth_grid);
cmin = 85; cmax = 110; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,squeeze(G363grid.scivars(:,4,ind_data_363)*med_gain_363),cvec,'linecolor','none'); hold on;
axis([min(G363grid.time_start(ind_data_363)) max(G363grid.time_start(ind_data_363)) mindepth maxdepth]); caxis([cmin cmax]);
colormap(cmocean('Balance','pivot',100)); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
plot(G363grid.time_start(ind_HYPM_363), 5*ones(size(ind_HYPM_363)), '.','markersize',15,'color',nicecolor('k'));
datetick('x',2,'keeplimits');
title('Glider 363 oxygen saturation', 'Fontsize', 12)

%% Extract HYPM profiles aligned with each glider profile

figure(10); clf
for i = 1:length(ind_HYPM_453)
    [minval(i), minind(i)] = min(abs(G453grid.time_start(ind_HYPM_453(i)) - Yr5_wfpgrid.time_start(Yr5_wfpgrid.ind_pair)));
    if round((i)/10) == ((i)/10)
        subplot(2,4,((i)/10))
        plot(Yr5_wfpgrid.oxygen_corr(:,minind(i))*gain_hypm, Yr5_wfpgrid.depth_grid, 'k.'); hold on;
        plot(squeeze(G453grid.scivars(:,3,ind_HYPM_453(i)))*med_gain_453, G453grid.depth_grid, 'r.'); hold on;
        axis ij
        xlim([260 310])
        ylim([50 1400])
        title(datestr(G453grid.time_start(ind_HYPM_453(i)),1))
    end
end
%%
figure(11); clf
for i = 1:length(ind_HYPM_363)
    [minval(i), minind(i)] = min(abs(G363grid.time_start(ind_HYPM_363(i)) - Yr5_wfpgrid.time_start(Yr5_wfpgrid.ind_pair)));
    if round((i)/10) == ((i)/10)
        subplot(3,4,((i)/10))
        plot(Yr5_wfpgrid.oxygen_corr(:,minind(i))*gain_hypm, Yr5_wfpgrid.depth_grid, 'k.'); hold on;
        plot(squeeze(G363grid.scivars(:,3,ind_HYPM_363(i)))*med_gain_363, G363grid.depth_grid, 'r.'); hold on;
        axis ij
        xlim([260 310])
        ylim([50 1400])
        title(datestr(G363grid.time_start(ind_HYPM_363(i)),1))
    end
end




