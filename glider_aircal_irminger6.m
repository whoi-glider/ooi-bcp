%%% Calculations and visualization of glider oxygen measurements during air
%%% calibration intervals between profiles

%% Calculate statistics for each air interval and keep only those with low stdev
    spacer = 10; %number of air measurements at beginning of each surfacing to cut
    tol = 1; %only keep air data if standard deviation of oxygen saturation measurements after cutting spacer is less than this value
[stats_560] = glideraircal_stats(G560, spacer, tol);
[stats_525] = glideraircal_stats(G525, spacer, tol);

%% Plot histogram of all air calibration data for full deployment
figure(1); clf
    hold on;
subplot(411)
    d = rem(G560.profile_index,1) == 0.5;
histogram(G560.oxygen_saturation(d),[99:0.2:110]); hold on;
    d = rem(G525.profile_index,1) == 0.5;
histogram(G525.oxygen_saturation(d),[99:0.2:110]); hold on;
title('Histogram of all raw oxygen saturation data during air intervals')
legend('Glider 560','Glider 525')

subplot(412)
    ind = find(stats_560(:,7) == 1);
histogram(stats_560(ind,3),[99:0.2:110]); hold on;
    ind = find(stats_525(:,7) == 1);
histogram(stats_525(ind,3),[99:0.2:110]); hold on;
title('Histogram of mean oxygen saturation data during "good" air intervals')
legend('Glider 560','Glider 525')

subplot(413)
    d = rem(G560.profile_index,1) == 0.5;
histogram(real(G560.O2_corr(d)),[300:2:440]); hold on;
    d = rem(G525.profile_index,1) == 0.5;
histogram(real(G525.O2_corr(d)),[300:2:440]); hold on;
title('Histogram of all raw, S-corrected oxygen concentration data during air intervals')
legend('Glider 560','Glider 525')

subplot(414)
    ind = find(stats_560(:,7) == 1);
histogram(stats_560(ind,1),[300:2:440]); hold on;
    ind = find(stats_525(:,7) == 1);
histogram(stats_525(ind,1),[300:2:440]); hold on;
title('Histogram of mean oxygen concentration data during "good" air intervals')
legend('Glider 560','Glider 525')


%% Look over time at air cal
figure(2); clf
    d = rem(G560.profile_index,1) == 0.5;
plot(G560.daten(d), G560.oxygen_saturation(d),'b.'); hold on;
    ind = find(stats_560(:,7) == 1);
plot(stats_560(ind,6), stats_560(ind,3), 'co','markerfacecolor','c'); hold on;
    d = rem(G525.profile_index,1) == 0.5;
plot(G525.daten(d), G525.oxygen_saturation(d),'.', 'color', 'r');
    ind = find(stats_525(:,7) == 1);
plot(stats_525(ind,6), stats_525(ind,3), 'ro','markerfacecolor','r'); hold on;
ylim([90 105])
xlim([datenum(2019,8,6,12,0,0) datenum(now)])
datetick('x',2,'keeplimits')
legend('Glider 560, all','Glider 560, "good" interval mean','Glider 560','Glider 525, "good" interval mean','location','northeast')
ylabel('Raw oxygen saturation (%)')
title('Oxygen measurements during glider surface intervals measuring air')
