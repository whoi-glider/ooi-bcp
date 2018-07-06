%%% Calculations and visualization of glider oxygen measurements during air
%%% calibration intervals between profiles

%% Calculate statistics for each air interval and keep only those with low stdev
    spacer = 10; %number of air measurements at beginning of each surfacing to cut
    tol = 1; %only keep air data if standard deviation of oxygen saturation measurements after cutting spacer is less than this value
[stats_453] = glideraircal_stats(G453, spacer, tol);
[stats_363] = glideraircal_stats(G363, spacer, tol);

%% Plot histogram of all air calibration data for full deployment
figure(1); clf
    hold on;
subplot(411)
    d = rem(G453.profile_index,1) == 0.5;
histogram(G453.oxygen_saturation(d),[99:0.2:110]); hold on;
    d = rem(G363.profile_index,1) == 0.5;
histogram(G363.oxygen_saturation(d),[99:0.2:110]); hold on;
title('Histogram of all raw oxygen saturation data during air intervals')
legend('Glider 453','Glider 363')

subplot(412)
    ind = find(stats_453(:,7) == 1);
histogram(stats_453(ind,3),[99:0.2:110]); hold on;
    ind = find(stats_363(:,7) == 1);
histogram(stats_363(ind,3),[99:0.2:110]); hold on;
title('Histogram of mean oxygen saturation data during "good" air intervals')
legend('Glider 453','Glider 363')

subplot(413)
    d = rem(G453.profile_index,1) == 0.5;
histogram(real(G453.O2_corr(d)),[300:2:440]); hold on;
    d = rem(G363.profile_index,1) == 0.5;
histogram(real(G363.O2_corr(d)),[300:2:440]); hold on;
title('Histogram of all raw, S-corrected oxygen concentration data during air intervals')
legend('Glider 453','Glider 363')

subplot(414)
    ind = find(stats_453(:,7) == 1);
histogram(stats_453(ind,1),[300:2:440]); hold on;
    ind = find(stats_363(:,7) == 1);
histogram(stats_363(ind,1),[300:2:440]); hold on;
title('Histogram of mean oxygen concentration data during "good" air intervals')
legend('Glider 453','Glider 363')


%% Look over time at air cal
figure(2); clf
    d = rem(G453.profile_index,1) == 0.5;
plot(G453.daten(d), G453.oxygen_saturation(d),'b.'); hold on;
    ind = find(stats_453(:,7) == 1);
plot(stats_453(ind,6), stats_453(ind,3), 'co','markerfacecolor','c'); hold on;
    d = rem(G363.profile_index,1) == 0.5;
plot(G363.daten(d), G363.oxygen_saturation(d),'.', 'color', nicecolor('ry'));
    ind = find(stats_363(:,7) == 1);
plot(stats_363(ind,6), stats_363(ind,3), 'ro','markerfacecolor','r'); hold on;
ylim([97 115])
xlim([datenum(2018,6,9,12,0,0) datenum(now)])
datetick('x',2,'keeplimits')
legend('Glider 453, all','Glider 453, "good" interval mean','Glider 453','Glider 363, "good" interval mean','location','northeast')
ylabel('Raw oxygen saturation (%)')
title('Oxygen measurements during glider surface intervals measuring air')
