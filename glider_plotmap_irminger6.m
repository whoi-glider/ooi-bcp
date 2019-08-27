load OOImooringLocations.mat
load castmeta_irminger6.mat

%find indices within given date range
    datemin = datenum(2019,8,4,12,0,0); datemax = datenum(now);
ind_525 = find(G525.daten > datemin & G525.daten < datemax);
ind_560 = find(G560.daten > datemin & G560.daten < datemax);

%plot maps showing locations of both gliders during given date range
figure(6); clf
    M1 = 14; M2 = 10;
latminplot = 59.6; latmaxplot = 60.1; lonminplot = -40.1; lonmaxplot = -38.9;
for i = 1:2
    subplot(1,2,i)
plot(OOImoorings.SUMO6(2), OOImoorings.SUMO6(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
plot(OOImoorings.HYPM6(2), OOImoorings.HYPM6(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
plot(OOImoorings.FLMA6(2), OOImoorings.FLMA6(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
plot(OOImoorings.FLMB6(2), OOImoorings.FLMB6(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
plot(OOImoorings.SUMO5(2), OOImoorings.SUMO5(1),'^','markerfacecolor','k','markeredgecolor','k','markersize',M2); hold on;
plot(OOImoorings.HYPM5(2), OOImoorings.HYPM5(1),'^','markerfacecolor','k','markeredgecolor','k','markersize',M2); hold on;
plot(OOImoorings.FLMA5(2), OOImoorings.FLMA5(1),'^','markerfacecolor','k','markeredgecolor','k','markersize',M2); hold on;
plot(OOImoorings.FLMB5(2), OOImoorings.FLMB5(1),'^','markerfacecolor','k','markeredgecolor','k','markersize',M2); hold on;
if i == 1
    scatter(G525.longitude(ind_525), G525.latitude(ind_525), [], G525.daten(ind_525),'filled','markeredgecolor','k'); hold on;
elseif i == 2
    scatter(G560.longitude(ind_560), G560.latitude(ind_560), [], G560.daten(ind_560),'filled','markeredgecolor','r'); hold on;
end
plot(-castmeta_irminger6.lon, castmeta_irminger6.lat,'o','markerfacecolor','m','markeredgecolor','k','markersize',M1/2); hold on;
h = colorbar; caxis([datenum(2019,8,4) datenum(now)]);
datetick(h,'y',6);
axis([lonminplot lonmaxplot latminplot latmaxplot])
title('Irminger-6 Glider surfacings and mooring locations'); xlabel('Longitude (deg W)'); ylabel('Latitude (deg N)')
end