
figure; clf
    M1 = 10; M2 = 14;
latminplot = 59.6; latmaxplot = 60.1; lonminplot = -40.1; lonmaxplot = -38.9;
plot(OOImoorings.SUMO4(2), OOImoorings.SUMO4(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
plot(OOImoorings.HYPM4(2), OOImoorings.HYPM4(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
plot(OOImoorings.FLMA4(2), OOImoorings.FLMA4(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
plot(OOImoorings.FLMB4(2), OOImoorings.FLMB4(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
plot(OOImoorings.SUMO5(2), OOImoorings.SUMO5(1),'^','markerfacecolor','k','markeredgecolor','k','markersize',M2); hold on;
plot(OOImoorings.HYPM5(2), OOImoorings.HYPM5(1),'^','markerfacecolor','k','markeredgecolor','k','markersize',M2); hold on;
plot(OOImoorings.FLMA5(2), OOImoorings.FLMA5(1),'^','markerfacecolor','k','markeredgecolor','k','markersize',M2); hold on;
plot(OOImoorings.FLMB5(2), OOImoorings.FLMB5(1),'^','markerfacecolor','k','markeredgecolor','k','markersize',M2); hold on;
scatter(G363.longitude, G363.latitude, [], G363.daten,'filled','markeredgecolor','k'); hold on;
scatter(G453.longitude, G453.latitude, [], G453.daten,'filled','markeredgecolor','r'); hold on;
h = colorbar; caxis([datenum(2018,6,9) datenum(now)]);
datetick(h,'y',6);
axis([lonminplot lonmaxplot latminplot latmaxplot])
title('Irminger-5 Glider surfacings and mooring locations'); xlabel('Longitude (deg W)'); ylabel('Latitude (deg N)')