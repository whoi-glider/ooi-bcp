% Look at paired up and down profiles during deployment period to check for lag and resulting bias

%% Plot to identify when gliders were measuring both up and down in sequence

%Pull out times of aligned casts
alignedcasts = [4,8,11];
[~, alignedind, ~] = intersect(castmeta_irminger6.castnum, alignedcasts);
timealigned = castmeta_irminger6.daytime(alignedind);

%Set beginning and end times
begtime = datenum(2019,8,5,0,0,0);
endtime = datenum(now); %datenum(2019,8,19);

figure(3); clf
    subplot(411)
plot(G525.daten, G525.depth_interp, 'k.'); hold on;
scatter(G525.daten, G525.depth_interp, [], G525.oxygen_saturation,'filled'); colorbar; caxis([85 100])
plot(timealigned, zeros(size(timealigned)), 'r.','markersize',20); hold on;
plot(timealigned-1, zeros(size(timealigned)), 'r.','markersize',12); hold on;
plot(timealigned+1, zeros(size(timealigned)), 'r.','markersize',12); hold on;
set(gca,'YDir','reverse'); 
xlim([begtime endtime])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 525, oxygen saturation')
    subplot(412)
plot(G525.daten, G525.depth_interp, 'k.'); hold on;
scatter(G525.daten, G525.depth_interp, [], G525.temperature,'filled'); colorbar; caxis([2.5 6.5])
plot(timealigned, zeros(size(timealigned)), 'r.','markersize',20); hold on;
plot(timealigned-1, zeros(size(timealigned)), 'r.','markersize',12); hold on;
plot(timealigned+1, zeros(size(timealigned)), 'r.','markersize',12); hold on;
set(gca,'YDir','reverse'); 
xlim([begtime endtime])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 525, temperature')
%     subplot(613)
% plot(G525.daten, G525.depth_interp, 'k.'); hold on;
% scatter(G525.daten, G525.depth_interp, [], G525.chlorophyll,'filled'); colorbar; caxis([0 0.2])
% set(gca,'YDir','reverse'); 
% xlim([datenum(2019,8,5,12,0,0) datenum(now)])
% ylim([0 1000])
% datetick('x',2,'keeplimits')
% ylabel('Depth (m)')
% title('Glider 525, chlorophyll')
    subplot(413)
plot(G560.daten, G560.depth_interp, 'k.'); hold on;
scatter(G560.daten, G560.depth_interp, [], G560.oxygen_saturation,'filled'); colorbar; caxis([85 110])
plot(timealigned, zeros(size(timealigned)), 'r.','markersize',20); hold on;
plot(timealigned-1, zeros(size(timealigned)), 'r.','markersize',12); hold on;
plot(timealigned+1, zeros(size(timealigned)), 'r.','markersize',12); hold on;
set(gca,'YDir','reverse'); 
xlim([begtime endtime])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 560, oxygen saturation')
    subplot(414)
plot(G560.daten, G560.depth_interp, 'k.'); hold on;
scatter(G560.daten, G560.depth_interp, [], G560.temperature,'filled'); colorbar; caxis([2.5 6.5])
plot(timealigned, zeros(size(timealigned)), 'r.','markersize',20); hold on;
plot(timealigned-1, zeros(size(timealigned)), 'r.','markersize',12); hold on;
plot(timealigned+1, zeros(size(timealigned)), 'r.','markersize',12); hold on;
set(gca,'YDir','reverse'); 
xlim([begtime endtime])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 560, temperature')
%     subplot(616)
% plot(G560.daten, G560.depth_interp, 'k.'); hold on;
% scatter(G560.daten, G560.depth_interp, [], G560.chlorophyll,'filled'); colorbar; caxis([0 0.2])
% set(gca,'YDir','reverse'); 
% xlim([datenum(2019,8,5,12,0,0) datenum(now)])
% ylim([0 1000])
% datetick('x',2,'keeplimits')
% ylabel('Depth (m)')
% title('Glider 560, chlorophyll')

%Time to check historesis effect on 560 = [datenum(2018,6,12,18,0,0) datenum(2018,6,13,16,0,0)]
%Time to check historesis effect on 525 = [datenum(2018,6,9,21,30,0) datenum(2018,6,10,18,30,0)] and 
% [datenum(2018,6,12,17,0,0) datenum(2018,6,13,6,0,0)]

%Note that there is an additional time interval June 19-20 to pull out as
%well

%% Pull out paired time intervals

ind560paired = find(G560.daten > datenum(2019,8,11,16,0,0) & G560.daten < datenum(2019,8,12,0,0,0));
ind560up = find(G560.profile_direction == -1);
ind560down = find(G560.profile_direction == 1);
ind560pairedup = intersect(ind560paired, ind560up);
ind560paireddown = intersect(ind560paired, ind560down);

%ind525paired_a = find(G525.daten > datenum(2018,6,9,21,30,0) & G525.daten < datenum(2018,6,10,18,30,0));
ind525paired_b = find(G525.daten > datenum(2019,8,11,16,0,0) & G525.daten < datenum(2019,8,12,0,0,0));
%ind525paired = union(ind525paired_a, ind525paired_b);
ind525paired = ind525paired_b;
ind525up = find(G525.profile_direction == -1);
ind525down = find(G525.profile_direction == 1);
ind525pairedup = intersect(ind525paired, ind525up);
ind525paireddown = intersect(ind525paired, ind525down);

depthgrid = [5:5:995];
G560_upgrid = NaN*ones(length(depthgrid),2);
G560_downgrid = NaN*ones(length(depthgrid),2);
for i = 1:length(depthgrid)
    ind_560 = find(G560.depth_interp >= depthgrid(i) - 2.5 & G560.depth_interp < depthgrid(i) + 2.5);
    ind_525 = find(G525.depth_interp >= depthgrid(i) - 2.5 & G525.depth_interp < depthgrid(i) + 2.5);
    G560_upgrid(i,1) = nanmean(G560.oxygen_saturation(intersect(ind_560,ind560pairedup)));
    G560_upgrid(i,2) = nanstd(G560.oxygen_saturation(intersect(ind_560,ind560pairedup)));
    G560_downgrid(i,1) = nanmean(G560.oxygen_saturation(intersect(ind_560,ind560paireddown)));
    G560_downgrid(i,2) = nanstd(G560.oxygen_saturation(intersect(ind_560,ind560paireddown)));
    G525_upgrid(i,1) = nanmean(G525.oxygen_saturation(intersect(ind_525,ind525pairedup)));
    G525_upgrid(i,2) = nanstd(G525.oxygen_saturation(intersect(ind_525,ind525pairedup)));
    G525_downgrid(i,1) = nanmean(G525.oxygen_saturation(intersect(ind_525,ind525paireddown)));
    G525_downgrid(i,2) = nanstd(G525.oxygen_saturation(intersect(ind_525,ind525paireddown)));
end

%% Plot paired time intervals together
figure(4); clf
L = 2.5;
M = 10;
    subplot(121)
plot(G560.oxygen_saturation(ind560pairedup), G560.depth_interp(ind560pairedup), '.', 'color', nicecolor('ry'),'markersize',M); hold on;
plot(G560.oxygen_saturation(ind560paireddown), G560.depth_interp(ind560paireddown), '.', 'color', nicecolor('kwwkry'),'markersize',M); hold on;
h1 = plot(G560_downgrid(:,1), depthgrid, '-', 'color', nicecolor('kwwkry'),'linewidth', L); hold on;
    plot(G560_downgrid(:,1) + G560_downgrid(:,2), depthgrid, '--', 'color', nicecolor('kwwkry'),'linewidth', L/2); hold on;
    plot(G560_downgrid(:,1) - G560_downgrid(:,2), depthgrid, '--', 'color', nicecolor('kwwkry'),'linewidth', L/2); hold on;
h2 = plot(G560_upgrid(:,1), depthgrid, '-', 'color', nicecolor('ry'),'linewidth', L); hold on;
    plot(G560_upgrid(:,1) + G560_upgrid(:,2), depthgrid, '--', 'color', nicecolor('ry'),'linewidth', L/2); hold on;
    plot(G560_upgrid(:,1) - G560_upgrid(:,2), depthgrid, '--', 'color', nicecolor('ry'),'linewidth', L/2); hold on;
set(gca,'YDir','reverse'); 
axis([85 110 5 1000])
ylabel('Depth (m)'); xlabel('O_2 saturation')
legend([h1 h2],'Glider 560, down','Glider 560, up','location','southeast')
title('Paired up and down profiles, 11 August 2019')
    subplot(122)
plot(G525.oxygen_saturation(ind525pairedup), G525.depth_interp(ind525pairedup), '.', 'color', nicecolor('b'),'markersize',M); hold on;
plot(G525.oxygen_saturation(ind525paireddown), G525.depth_interp(ind525paireddown), '.', 'color', nicecolor('bcwwkk'),'markersize',M); hold on;
h1 = plot(G525_downgrid(:,1), depthgrid, '-', 'color', nicecolor('bcwwkk'),'linewidth', L); hold on;
    plot(G525_downgrid(:,1) + G525_downgrid(:,2), depthgrid, '--', 'color', nicecolor('bcwwkk'),'linewidth', L/2); hold on;
    plot(G525_downgrid(:,1) - G525_downgrid(:,2), depthgrid, '--', 'color', nicecolor('bcwwkk'),'linewidth', L/2); hold on;
h2 = plot(G525_upgrid(:,1), depthgrid, '-', 'color', nicecolor('b'),'linewidth', L); hold on;
    plot(G525_upgrid(:,1) + G525_upgrid(:,2), depthgrid, '--', 'color', nicecolor('b'),'linewidth', L/2); hold on;
    plot(G525_upgrid(:,1) - G525_upgrid(:,2), depthgrid, '--', 'color', nicecolor('b'),'linewidth', L/2); hold on;
set(gca,'YDir','reverse'); 
axis([85 110 5 1000])
ylabel('Depth (m)'); xlabel('O_2 saturation')
legend([h1 h2],'Glider 525, down','Glider 525, up','location','southeast')
title('Paired up and down profiles, 11 August 2019')

%%
figure(5); clf
%Other variables from GL560 on August 11
L = 2.5;
M = 10;
    subplot(141)
plot(G560.temperature(ind560pairedup), G560.depth_interp(ind560pairedup), '.', 'color', nicecolor('rrk'),'markersize',M); hold on;
set(gca,'YDir','reverse'); 
ylabel('Depth (m)'); xlabel('Temperature')
title('Temperature, 8/11/19, GL560')
    subplot(142)
plot(G560.salinity(ind560pairedup), G560.depth_interp(ind560pairedup), '.', 'color', nicecolor('ccbk'),'markersize',M); hold on;
set(gca,'YDir','reverse'); 
ylabel('Depth (m)'); xlabel('Salinity')
title('Salinity, 8/11/19, GL560')
    subplot(143)
plot(G560.chlorophyll(ind560pairedup), G560.depth_interp(ind560pairedup), '.', 'color', nicecolor('ggy'),'markersize',M); hold on;
set(gca,'YDir','reverse'); 
ylabel('Depth (m)'); xlabel('Chlorophyll')
title('Chlorophyll, 8/11/19, GL560')
    subplot(144)
plot(G560.backscatter(ind560pairedup), G560.depth_interp(ind560pairedup), '.', 'color', nicecolor('ry'),'markersize',M); hold on;
set(gca,'YDir','reverse'); 
ylabel('Depth (m)'); xlabel('Backscatter')
title('Backscatter, 8/11/19, GL560')