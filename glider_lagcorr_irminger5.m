%%%%%% This is still just starter code to figure out lag correction, but
%%%%%% seems to work on test data and example up-down profiles from
%%%%%% Irminger-6. The trick here is that the lag time scale (tau) is ~30s
%%%%%% and in telemtered glider data, that is the resolution of
%%%%%% measurements, limiting the ability to apply this method (i.e. the
%%%%%% lag calculation is only incorporating 2 measurements). This should
%%%%%% be checked using recovered data from 363 from Irminger-5 before
%%%%%% applying broadly, but seems likely to be something we will want to
%%%%%% use.

% testdown.oxygensat = G363.oxygen_saturation(ind363paireddown);
% testdown.depth_interp = G363.depth_interp(ind363paireddown);
% testdown.daten = G363.daten(ind363paireddown);
% 
% testup.oxygensat = G363.oxygen_saturation(ind363pairedup);
% testup.depth_interp = G363.depth_interp(ind363pairedup);
% testup.daten = G363.daten(ind363pairedup);
% 
% inddata_down = find(isnan(testdown.oxygensat) + isnan(testdown.depth_interp) + isnan(testdown.daten) == 0);
% inddata_up = find(isnan(testup.oxygensat) + isnan(testup.depth_interp) + isnan(testup.daten) == 0);
    secinday = 60*60*24;
    tau = 30/secinday;
    timetol = 600*10;
% [O2_lagcorr_down] = lagCorr(testdown.oxygensat(inddata_down), testdown.daten(inddata_down), tau, timetol);
% [O2_lagcorr_up] = lagCorr(testup.oxygensat(inddata_up), testup.daten(inddata_up), tau, timetol);
% 
% figure(100); clf
% plot(testup.daten, testup.oxygensat, 'k.'); hold on;
% plot(testdown.daten, testdown.oxygensat, '.','color',nicecolor('kkw')); hold on;
% plot(testup.daten(inddata_up), O2_lagcorr_up, 'r.'); hold on;
% plot(testdown.daten(inddata_down), O2_lagcorr_down, 'b.'); hold on;
% datetick('x')
% 
% figure(110); clf
% plot(testup.oxygensat, testup.depth_interp, 'k.'); hold on;
% plot(testdown.oxygensat, testdown.depth_interp, '.','color',nicecolor('kkw')); hold on;
% plot(O2_lagcorr_up, testup.depth_interp(inddata_up), 'r.'); hold on;
% plot(O2_lagcorr_down, testdown.depth_interp(inddata_down), 'b.'); hold on;
% axis ij

%% Plot paired time intervals together



figure(12); clf
L = 1.5;
M = 10;
M1 = 4;
M2 = 4;
    %subplot(121)
plot(G363.oxygen_saturation(ind363pairedup), G363.depth_interp(ind363pairedup), '.', 'color', nicecolor('ry'),'markersize',M1); hold on;
plot(G363.oxygen_saturation(ind363paireddown), G363.depth_interp(ind363paireddown), '.', 'color', nicecolor('kwwkry'),'markersize',M1); hold on;
plot(G363_upgrid_lagcorr(:,3), G363_upgrid_lagcorr(:,2), 'r.','markersize',M2);
plot(G363_downgrid_lagcorr(:,3), G363_downgrid_lagcorr(:,2), 'm.','markersize',M2);
set(gca,'YDir','reverse'); 
axis([85 110 5 1000])
ylabel('Depth (m)'); xlabel('O_2 saturation')
%legend([h1 h2],'Glider 363, down','Glider 363, up','location','southeast')
title('Paired up and down profiles, 12-13 June 2018')
%     subplot(122)
% plot(G363.oxygen_saturation(ind363pairedup), G363.depth_interp(ind363pairedup), '.', 'color', nicecolor('b'),'markersize',M); hold on;
% plot(G363.oxygen_saturation(ind363paireddown), G363.depth_interp(ind363paireddown), '.', 'color', nicecolor('bcwwkk'),'markersize',M); hold on;
% h1 = plot(G363_downgrid(:,1), depthgrid, '-', 'color', nicecolor('bcwwkk'),'linewidth', L); hold on;
% h2 = plot(G363_upgrid(:,1), depthgrid, '-', 'color', nicecolor('b'),'linewidth', L); hold on;
% 
% set(gca,'YDir','reverse'); 
% axis([85 110 5 1000])
% ylabel('Depth (m)'); xlabel('O_2 saturation')
% legend([h1 h2],'Glider 525, down','Glider 525, up','location','southeast')
% title('Paired up and down profiles, 11 August 2019')

