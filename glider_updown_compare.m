% Look at paired up and down profiles during deployment period to check for lag and resulting bias

%% Plot to identify when gliders were measuring both up and down in sequence
% figure(3); clf
%     subplot(211)
% plot(G363.daten, G363.depth_interp, 'k.'); hold on;
% scatter(G363.daten, G363.depth_interp, [], G363.oxygen_saturation,'filled'); colorbar; caxis([85 100])
% set(gca,'YDir','reverse'); 
% xlim([datenum(2018,6,9,12,0,0) datenum(2019,3,14,0,0,0)])
% ylim([0 400])
% datetick('x',2,'keeplimits')
% ylabel('Depth (m)')
% title('Glider 363, oxygen saturation')
%     subplot(212)
% plot(G363.daten, G363.depth_interp, 'k.'); hold on;
% scatter(G363.daten, G363.depth_interp, [], G363.temperature,'filled'); colorbar; caxis([2.5 6.5])
% set(gca,'YDir','reverse'); 
% xlim([datenum(2018,6,9,12,0,0) datenum(2019,3,14,0,0,0)])
% ylim([0 400])
% datetick('x',2,'keeplimits')
% ylabel('Depth (m)')
% title('Glider 363, temperature')
%     subplot(212)
% plot(G453.daten, G453.depth_interp, 'k.'); hold on;
% scatter(G453.daten, G453.depth_interp, [], G453.oxygen_saturation,'filled'); colorbar; caxis([80 110])
% set(gca,'YDir','reverse'); 
% xlim([datenum(2018,6,9,12,0,0) datenum(2019,1,30,0,0,0)])
% ylim([0 400])
% datetick('x',2,'keeplimits')
% ylabel('Depth (m)')
% title('Glider 453, oxygen saturation')

%Time to check historesis effect on 453 = [datenum(2018,6,12,18,0,0) datenum(2018,6,13,16,0,0)]
%Time to check historesis effect on 363 = [datenum(2018,6,9,21,30,0) datenum(2018,6,10,18,30,0)] and 
% [datenum(2018,6,12,17,0,0) datenum(2018,6,13,6,0,0)]

%Note that there is an additional time interval June 19-20 to pull out as
%well

%% Pull out paired time intervals

ind453paired = find(G453.daten > datenum(2018,6,12,18,0,0) & G453.daten < datenum(2018,6,13,16,0,0));
ind453up = find(G453.profile_direction == -1);
ind453down = find(G453.profile_direction == 1);
ind453pairedup = intersect(ind453paired, ind453up);
ind453paireddown = intersect(ind453paired, ind453down);

ind363paired_a = find(G363.daten > datenum(2018,6,9,21,30,0) & G363.daten < datenum(2018,6,10,18,30,0));
ind363paired_b = find(G363.daten > datenum(2018,6,12,17,0,0) & G363.daten < datenum(2018,6,13,6,0,0));
%ind363paired = union(ind363paired_a, ind363paired_b);
ind363paired = ind363paired_b;
ind363up = find(G363.profile_direction == -1);
ind363down = find(G363.profile_direction == 1);
ind363pairedup = intersect(ind363paired, ind363up);
ind363paireddown = intersect(ind363paired, ind363down);

%% Apply lag correction, testing range of choices for tau
    secinday = 60*60*24;
    tau = [15:5:90]/secinday;
    timetol = 600*10;

    inddata = find(isnan(G363.oxygen_saturation(ind363pairedup)) + isnan(G363.depth_interp(ind363pairedup)) + isnan(G363.daten(ind363pairedup)) == 0);
G363_up_lagcorr = NaN*zeros(length(inddata),3);
G363_up_lagcorr(:,1) = G363.daten(ind363pairedup(inddata));
G363_up_lagcorr(:,2) = G363.depth_interp(ind363pairedup(inddata));
for i = 1:length(tau)
    G363_up_lagcorr(:,2+i) = lagCorr(G363.oxygen_saturation(ind363pairedup(inddata)), G363_up_lagcorr(:,1), tau(i), timetol);
end
    inddata = find(isnan(G363.oxygen_saturation(ind363paireddown)) + isnan(G363.depth_interp(ind363paireddown)) + isnan(G363.daten(ind363paireddown)) == 0);
G363_down_lagcorr = NaN*zeros(length(inddata),3);
G363_down_lagcorr(:,1) = G363.daten(ind363paireddown(inddata));
G363_down_lagcorr(:,2) = G363.depth_interp(ind363paireddown(inddata));
for i = 1:length(tau)
    G363_down_lagcorr(:,2+i) = lagCorr(G363.oxygen_saturation(ind363paireddown(inddata)), G363_down_lagcorr(:,1), tau(i), timetol);
end

%Repeat for 453
    inddata = find(isnan(G453.oxygen_saturation(ind453pairedup)) + isnan(G453.depth_interp(ind453pairedup)) + isnan(G453.daten(ind453pairedup)) == 0);
G453_up_lagcorr = NaN*zeros(length(inddata),3);
G453_up_lagcorr(:,1) = G453.daten(ind453pairedup(inddata));
G453_up_lagcorr(:,2) = G453.depth_interp(ind453pairedup(inddata));
for i = 1:length(tau)
    G453_up_lagcorr(:,2+i) = lagCorr(G453.oxygen_saturation(ind453pairedup(inddata)), G453_up_lagcorr(:,1), tau(i), timetol);
end

    inddata = find(isnan(G453.oxygen_saturation(ind453paireddown)) + isnan(G453.depth_interp(ind453paireddown)) + isnan(G453.daten(ind453paireddown)) == 0);
G453_down_lagcorr = NaN*zeros(length(inddata),3);
G453_down_lagcorr(:,1) = G453.daten(ind453paireddown(inddata));
G453_down_lagcorr(:,2) = G453.depth_interp(ind453paireddown(inddata));
for i = 1:length(tau)
    G453_down_lagcorr(:,2+i) = lagCorr(G453.oxygen_saturation(ind453paireddown(inddata)), G453_down_lagcorr(:,1), tau(i), timetol);
end

%% Grid on even depth intervals 
depthgrid = [5:5:995];
G453_upgrid = NaN*ones(length(depthgrid),2);
G453_downgrid = NaN*ones(length(depthgrid),2);
G453_upgrid_lagcorr = NaN*ones(length(depthgrid),length(tau));
G453_downgrid_lagcorr = NaN*ones(length(depthgrid),length(tau));
G363_upgrid = NaN*ones(length(depthgrid),2);
G363_downgrid = NaN*ones(length(depthgrid),2);
G363_upgrid_lagcorr = NaN*ones(length(depthgrid),length(tau));
G363_downgrid_lagcorr = NaN*ones(length(depthgrid),length(tau));
for i = 1:length(depthgrid)
    ind_453 = find(G453.depth_interp >= depthgrid(i) - 2.5 & G453.depth_interp < depthgrid(i) + 2.5);
    ind_363 = find(G363.depth_interp >= depthgrid(i) - 2.5 & G363.depth_interp < depthgrid(i) + 2.5);
    G453_upgrid(i,1) = nanmean(G453.oxygen_saturation(intersect(ind_453,ind453pairedup)));
    G453_upgrid(i,2) = nanstd(G453.oxygen_saturation(intersect(ind_453,ind453pairedup)));
    G453_downgrid(i,1) = nanmean(G453.oxygen_saturation(intersect(ind_453,ind453paireddown)));
    G453_downgrid(i,2) = nanstd(G453.oxygen_saturation(intersect(ind_453,ind453paireddown)));
    G363_upgrid(i,1) = nanmean(G363.oxygen_saturation(intersect(ind_363,ind363pairedup)));
    G363_upgrid(i,2) = nanstd(G363.oxygen_saturation(intersect(ind_363,ind363pairedup)));
    G363_downgrid(i,1) = nanmean(G363.oxygen_saturation(intersect(ind_363,ind363paireddown)));
    G363_downgrid(i,2) = nanstd(G363.oxygen_saturation(intersect(ind_363,ind363paireddown)));
    
    ind_363_downgrid_lagcorr = find(G363_down_lagcorr(:,2) >= depthgrid(i) - 2.5 & G363_down_lagcorr(:,2) < depthgrid(i) + 2.5);
    ind_363_upgrid_lagcorr = find(G363_up_lagcorr(:,2) >= depthgrid(i) - 2.5 & G363_up_lagcorr(:,2) < depthgrid(i) + 2.5);
    for j = 1:length(tau)
        G363_downgrid_lagcorr(i,j) = nanmean(G363_down_lagcorr(ind_363_downgrid_lagcorr,j+2));
        G363_upgrid_lagcorr(i,j) = nanmean(G363_up_lagcorr(ind_363_upgrid_lagcorr,j+2));
    end
    
    ind_453_downgrid_lagcorr = find(G453_down_lagcorr(:,2) >= depthgrid(i) - 2.5 & G453_down_lagcorr(:,2) < depthgrid(i) + 2.5);
    ind_453_upgrid_lagcorr = find(G453_up_lagcorr(:,2) >= depthgrid(i) - 2.5 & G453_up_lagcorr(:,2) < depthgrid(i) + 2.5);
    for j = 1:length(tau)
        G453_downgrid_lagcorr(i,j) = nanmean(G453_down_lagcorr(ind_453_downgrid_lagcorr,j+2));
        G453_upgrid_lagcorr(i,j) = nanmean(G453_up_lagcorr(ind_453_upgrid_lagcorr,j+2));
    end
end

%% Plot paired time intervals together
figure(4); clf
L = 0.75;
M = 3;
    subplot(121)
h6 = plot(G453.oxygen_saturation(ind453pairedup), G453.depth_interp(ind453pairedup), '.', 'color', nicecolor('ryw'),'markersize',M); hold on;
h5 = plot(G453.oxygen_saturation(ind453paireddown), G453.depth_interp(ind453paireddown), '.', 'color', nicecolor('kwwkryww'),'markersize',M); hold on;
h1 = plot(G453_downgrid(:,1), depthgrid, '-', 'color', nicecolor('kwwkry'),'linewidth', L); hold on;
h3 = plot(G453_downgrid_lagcorr(:,8), depthgrid, '-', 'color', nicecolor('kwwkrykk'),'linewidth', L+1); hold on;
    %plot(G453_downgrid(:,1) + G453_downgrid(:,2), depthgrid, '--', 'color', nicecolor('kwwkry'),'linewidth', L/2); hold on;
    %plot(G453_downgrid(:,1) - G453_downgrid(:,2), depthgrid, '--', 'color', nicecolor('kwwkry'),'linewidth', L/2); hold on;
h2 = plot(G453_upgrid(:,1), depthgrid, '-', 'color', nicecolor('ry'),'linewidth', L); hold on;
h4 = plot(G453_upgrid_lagcorr(:,8), depthgrid, '-', 'color', nicecolor('ryk'),'linewidth', L+1); hold on;
    %plot(G453_upgrid(:,1) + G453_upgrid(:,2), depthgrid, '--', 'color', nicecolor('ry'),'linewidth', L/2); hold on;
    %plot(G453_upgrid(:,1) - G453_upgrid(:,2), depthgrid, '--', 'color', nicecolor('ry'),'linewidth', L/2); hold on;
set(gca,'YDir','reverse'); 
axis([88 113 5 250])
ylabel('Depth (m)'); xlabel('O_2 saturation')
legend([h5 h1 h3 h6 h2 h4],'All dive measurements','Mean dive profile','Lag-corrected mean dive profile',...
    'All climb measurements','Mean climb profile','Lag-corrected mean climb profile',...
    'location','southeast')
title('Glider 453, 12-13 June 2018')

    subplot(122)
h6 = plot(G363.oxygen_saturation(ind363pairedup), G363.depth_interp(ind363pairedup), '.', 'color', nicecolor('bw'),'markersize',M); hold on;
h5 = plot(G363.oxygen_saturation(ind363paireddown), G363.depth_interp(ind363paireddown), '.', 'color', nicecolor('bcwwkkww'),'markersize',M); hold on;
h1 = plot(G363_downgrid(:,1), depthgrid, '-', 'color', nicecolor('bcwwkk'),'linewidth', L); hold on;
h3 = plot(G363_downgrid_lagcorr(:,8), depthgrid, '-', 'color', nicecolor('bbwwkkkk'),'linewidth', L+1); hold on;
h2 = plot(G363_upgrid(:,1), depthgrid, '-', 'color', nicecolor('b'),'linewidth', L); hold on;
h4 = plot(G363_upgrid_lagcorr(:,8), depthgrid, '-', 'color', nicecolor('bbk'),'linewidth', L+1); hold on;
set(gca,'YDir','reverse'); 
axis([90 116 5 250])
ylabel('Depth (m)'); xlabel('O_2 saturation')
legend([h5 h1 h3 h6 h2 h4],'All dive measurements','Mean dive profile','Lag-corrected mean dive profile',...
    'All climb measurements','Mean climb profile','Lag-corrected mean climb profile',...
    'location','southeast')
title('Glider 363, 12-13 June 2018')

%% Quantitative estimate of difference between up- and down- grid

    indmaxsurf = 60;
[v,i] = min(abs(nanmean(G363_upgrid_lagcorr - G363_downgrid_lagcorr)))
[v,i] = min(abs(nanmean(G363_upgrid_lagcorr(1:indmaxsurf,:) - G363_downgrid_lagcorr(1:indmaxsurf,:))))
[v,i] = min(nanmean(abs(G363_upgrid_lagcorr - G363_downgrid_lagcorr)))
[v,i] = min(nanmean(abs(G363_upgrid_lagcorr(1:indmaxsurf,:) - G363_downgrid_lagcorr(1:indmaxsurf,:))))

[v,i] = min(abs(nanmean(G453_upgrid_lagcorr - G453_downgrid_lagcorr)))
[v,i] = min(abs(nanmean(G453_upgrid_lagcorr(1:indmaxsurf,:) - G453_downgrid_lagcorr(1:indmaxsurf,:))))
[v,i] = min(nanmean(abs(G453_upgrid_lagcorr - G453_downgrid_lagcorr)))
[v,i] = min(nanmean(abs(G453_upgrid_lagcorr(1:indmaxsurf,:) - G453_downgrid_lagcorr(1:indmaxsurf,:))))