%%%%%% This is still just starter code to figure out lag correction, but
%%%%%% seems to work on test data and example up-down profiles from
%%%%%% Irminger-6. The trick here is that the lag time scale (tau) is ~30s
%%%%%% and in telemtered glider data, that is the resolution of
%%%%%% measurements, limiting the ability to apply this method (i.e. the
%%%%%% lag calculation is only incorporating 2 measurements). This should
%%%%%% be checked using recovered data from 363 from Irminger-5 before
%%%%%% applying broadly, but seems likely to be something we will want to
%%%%%% use.

test.oxygensat = G560.oxygen_saturation(ind560pairedup);
test.depth_interp = G560.depth_interp(ind560pairedup);
test.daten = G560.daten(ind560pairedup);

inddata = find(isnan(test.oxygensat) + isnan(test.depth_interp) + isnan(test.daten) == 0);
    secinday = 60*60*24;
    tau = 40/secinday;
    timetol = 600*10;
[O2_lagcorr] = lagCorr(test.oxygensat(inddata), test.daten(inddata), tau, timetol);

figure(10); clf
plot(test.daten, test.oxygensat, 'k.'); hold on;
plot(test.daten(inddata), O2_lagcorr, 'r.'); hold on;
datetick('x')

figure(11); clf
plot(test.oxygensat, test.depth_interp, 'k.'); hold on;
plot(O2_lagcorr, test.depth_interp(inddata), 'r.'); hold on;
axis ij

%% Plot paired time intervals together

    tau = 40/secinday;
    timetol = 600*10;

    inddata = find(isnan(G560.oxygen_saturation(ind560pairedup)) + isnan(G560.depth_interp(ind560pairedup)) + isnan(G560.daten(ind560pairedup)) == 0);
G560_upgrid_lagcorr = NaN*zeros(length(inddata),3);
G560_upgrid_lagcorr(:,1) = G560.daten(ind560pairedup(inddata));
G560_upgrid_lagcorr(:,2) = G560.depth_interp(ind560pairedup(inddata));
G560_upgrid_lagcorr(:,3) = lagCorr(G560.oxygen_saturation(ind560pairedup(inddata)), G560_upgrid_lagcorr(:,1), tau, timetol);
    inddata = find(isnan(G560.oxygen_saturation(ind560paireddown)) + isnan(G560.depth_interp(ind560paireddown)) + isnan(G560.daten(ind560paireddown)) == 0);
G560_downgrid_lagcorr = NaN*zeros(length(inddata),3);
G560_downgrid_lagcorr(:,1) = G560.daten(ind560paireddown(inddata));
G560_downgrid_lagcorr(:,2) = G560.depth_interp(ind560paireddown(inddata));
G560_downgrid_lagcorr(:,3) = lagCorr(G560.oxygen_saturation(ind560paireddown(inddata)), G560_downgrid_lagcorr(:,1), tau, timetol);



figure(12); clf
L = 1.5;
M = 10;
    subplot(121)
plot(G560.oxygen_saturation(ind560pairedup), G560.depth_interp(ind560pairedup), '.', 'color', nicecolor('ry'),'markersize',M); hold on;
plot(G560.oxygen_saturation(ind560paireddown), G560.depth_interp(ind560paireddown), '.', 'color', nicecolor('kwwkry'),'markersize',M); hold on;
h1 = plot(G560_downgrid(:,1), depthgrid, '-', 'color', nicecolor('kwwkry'),'linewidth', L); hold on;
h2 = plot(G560_upgrid(:,1), depthgrid, '-', 'color', nicecolor('ry'),'linewidth', L); hold on;
plot(G560_upgrid_lagcorr(:,3), G560_upgrid_lagcorr(:,2), 'r-','markersize',M);
plot(G560_downgrid_lagcorr(:,3), G560_downgrid_lagcorr(:,2), 'm-','markersize',M);
set(gca,'YDir','reverse'); 
axis([85 110 5 1000])
ylabel('Depth (m)'); xlabel('O_2 saturation')
legend([h1 h2],'Glider 560, down','Glider 560, up','location','southeast')
title('Paired up and down profiles, 11 August 2019')
    subplot(122)
plot(G525.oxygen_saturation(ind525pairedup), G525.depth_interp(ind525pairedup), '.', 'color', nicecolor('b'),'markersize',M); hold on;
plot(G525.oxygen_saturation(ind525paireddown), G525.depth_interp(ind525paireddown), '.', 'color', nicecolor('bcwwkk'),'markersize',M); hold on;
h1 = plot(G525_downgrid(:,1), depthgrid, '-', 'color', nicecolor('bcwwkk'),'linewidth', L); hold on;
h2 = plot(G525_upgrid(:,1), depthgrid, '-', 'color', nicecolor('b'),'linewidth', L); hold on;

set(gca,'YDir','reverse'); 
axis([85 110 5 1000])
ylabel('Depth (m)'); xlabel('O_2 saturation')
legend([h1 h2],'Glider 525, down','Glider 525, up','location','southeast')
title('Paired up and down profiles, 11 August 2019')

%% Test of lag correction in CTD data
%%%% Note that since lag scale incorporates many more points, artefacts are
%%%% added to the data by doing this, likely leading to better results just
%%%% by averaging the up & down casts. This will work for casts without
%%%% bottle samples (no stops on way up), but will need a different
%%%% approach on casts where the upcast included long enough soaks to
%%%% eliminate lag issue.

tau = 30; %note that function is written for time, but this is in depth space --> know that winch went ~60 m/min
timetol = 100;

A = lagCorr(cast{castnums(i)}.O2corr(1:cast{castnums(i)}.maxindex), cast{castnums(i)}.D(1:cast{castnums(i)}.maxindex), tau, timetol);

figure(30); clf
    i = 3;
plot(cast{castnums(i)}.O2corr(1:cast{castnums(i)}.maxindex), cast{castnums(i)}.D(1:cast{castnums(i)}.maxindex), 'k.'); hold on;
plot(lagCorr(cast{castnums(i)}.O2corr(1:cast{castnums(i)}.maxindex), cast{castnums(i)}.D(1:cast{castnums(i)}.maxindex), tau, timetol),...
    cast{castnums(i)}.D(1:cast{castnums(i)}.maxindex), 'm.'); hold on;
plot(cast{castnums(i)}.O2corr(cast{castnums(i)}.maxindex:end), cast{castnums(i)}.D(cast{castnums(i)}.maxindex:end), 'r.'); hold on;
%Reversed sign on depth to account for the fact that this is the up rather
%than down cast
plot(lagCorr(cast{castnums(i)}.O2corr(cast{castnums(i)}.maxindex:end), -cast{castnums(i)}.D(cast{castnums(i)}.maxindex:end), tau, timetol),...
    cast{castnums(i)}.D(cast{castnums(i)}.maxindex:end), 'g.'); hold on;
axis ij
xlabel('Oxygen (\muM) from Aanderaa optode')
ylabel('Depth (m)')
title(castnames(i))