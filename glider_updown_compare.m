% Look at paired up and down profiles during deployment period to check for lag and resulting bias

%% Plot to identify when gliders were measuring both up and down in sequence
figure(3); clf
    subplot(211)
plot(G363.daten, G363.depth_interp, 'k.'); hold on;
scatter(G363.daten, G363.depth_interp, [], G363.oxygen_saturation,'filled'); colorbar; caxis([80 110])
set(gca,'YDir','reverse'); 
xlim([datenum(2018,6,9,12,0,0) datenum(2018,6,22,0,0,0)])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 363, oxygen saturation')
    subplot(212)
plot(G453.daten, G453.depth_interp, 'k.'); hold on;
scatter(G453.daten, G453.depth_interp, [], G453.oxygen_saturation,'filled'); colorbar; caxis([80 110])
set(gca,'YDir','reverse'); 
xlim([datenum(2018,6,9,12,0,0) datenum(2018,6,22,0,0,0)])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 453, oxygen saturation')

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

depthgrid = [5:5:995];
G453_upgrid = NaN*ones(length(depthgrid),2);
G453_downgrid = NaN*ones(length(depthgrid),2);
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
end

%% Plot paired time intervals together
figure(4); clf
L = 2.5;
M = 10;
    subplot(121)
plot(G453.oxygen_saturation(ind453pairedup), G453.depth_interp(ind453pairedup), '.', 'color', nicecolor('ry'),'markersize',M); hold on;
plot(G453.oxygen_saturation(ind453paireddown), G453.depth_interp(ind453paireddown), '.', 'color', nicecolor('kwwkry'),'markersize',M); hold on;
h1 = plot(G453_downgrid(:,1), depthgrid, '-', 'color', nicecolor('kwwkry'),'linewidth', L); hold on;
    plot(G453_downgrid(:,1) + G453_downgrid(:,2), depthgrid, '--', 'color', nicecolor('kwwkry'),'linewidth', L/2); hold on;
    plot(G453_downgrid(:,1) - G453_downgrid(:,2), depthgrid, '--', 'color', nicecolor('kwwkry'),'linewidth', L/2); hold on;
h2 = plot(G453_upgrid(:,1), depthgrid, '-', 'color', nicecolor('ry'),'linewidth', L); hold on;
    plot(G453_upgrid(:,1) + G453_upgrid(:,2), depthgrid, '--', 'color', nicecolor('ry'),'linewidth', L/2); hold on;
    plot(G453_upgrid(:,1) - G453_upgrid(:,2), depthgrid, '--', 'color', nicecolor('ry'),'linewidth', L/2); hold on;
set(gca,'YDir','reverse'); 
axis([85 110 5 1000])
ylabel('Depth (m)'); xlabel('O_2 saturation')
legend([h1 h2],'Glider 453, down','Glider 453, up','location','southeast')
title('Paired up and down profiles, 12-13 June 2018')
    subplot(122)
plot(G363.oxygen_saturation(ind363pairedup), G363.depth_interp(ind363pairedup), '.', 'color', nicecolor('b'),'markersize',M); hold on;
plot(G363.oxygen_saturation(ind363paireddown), G363.depth_interp(ind363paireddown), '.', 'color', nicecolor('bcwwkk'),'markersize',M); hold on;
h1 = plot(G363_downgrid(:,1), depthgrid, '-', 'color', nicecolor('bcwwkk'),'linewidth', L); hold on;
    plot(G363_downgrid(:,1) + G363_downgrid(:,2), depthgrid, '--', 'color', nicecolor('bcwwkk'),'linewidth', L/2); hold on;
    plot(G363_downgrid(:,1) - G363_downgrid(:,2), depthgrid, '--', 'color', nicecolor('bcwwkk'),'linewidth', L/2); hold on;
h2 = plot(G363_upgrid(:,1), depthgrid, '-', 'color', nicecolor('b'),'linewidth', L); hold on;
    plot(G363_upgrid(:,1) + G363_upgrid(:,2), depthgrid, '--', 'color', nicecolor('b'),'linewidth', L/2); hold on;
    plot(G363_upgrid(:,1) - G363_upgrid(:,2), depthgrid, '--', 'color', nicecolor('b'),'linewidth', L/2); hold on;
set(gca,'YDir','reverse'); 
axis([85 110 5 1000])
ylabel('Depth (m)'); xlabel('O_2 saturation')
legend([h1 h2],'Glider 363, down','Glider 363, up','location','southeast')
title('Paired up and down profiles, 12-13 June 2018')