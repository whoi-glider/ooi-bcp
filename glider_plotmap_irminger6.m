%% Make maps showing glider tracks, mooring locations, and CTD coast locations
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
    scatter(G560.longitude(ind_560), G560.latitude(ind_560), [], G560.daten(ind_560),'filled','markeredgecolor','k'); hold on;
end
plot(-castmeta_irminger6.lon, castmeta_irminger6.lat,'o','markerfacecolor','r','markeredgecolor','k','markersize',M1/2); hold on;
h = colorbar; caxis([datenum(2019,8,4) datenum(now)]);
datetick(h,'y',6);
axis([lonminplot lonmaxplot latminplot latmaxplot])
if i == 1
    title('GL525 surfacings with CTD cast and mooring locations');
elseif i == 2
    title('GL560 surfacings with CTD cast and mooring locations');
end
xlabel('Longitude (deg W)'); ylabel('Latitude (deg N)')
end
%% Plot focusing on glider locations corresponding in space/time to CTD casts
timetol = 1; %one day before or after cast

    figure(7); clf
for i = 1:length(alignedcasts)
    %Determine time range for each cast
ind_525 = find(G525.daten > timealigned(i) - timetol & G525.daten < timealigned(i) + timetol);
ind_560 = find(G560.daten > timealigned(i) - timetol & G560.daten < timealigned(i) + timetol);
    %Plot glider surfacings and cast location    
    subplot(1,3,i)
plot(OOImoorings.SUMO6(2), OOImoorings.SUMO6(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
plot(OOImoorings.HYPM6(2), OOImoorings.HYPM6(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
plot(OOImoorings.FLMA6(2), OOImoorings.FLMA6(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
plot(OOImoorings.FLMB6(2), OOImoorings.FLMB6(1),'v','markerfacecolor','k','markeredgecolor','k','markersize',M1); hold on;
    scatter(G525.longitude(ind_525), G525.latitude(ind_525), [], G525.daten(ind_525) - timealigned(i),'filled','markeredgecolor','k'); hold on;
    scatter(G560.longitude(ind_560), G560.latitude(ind_560), [], G560.daten(ind_560) - timealigned(i),'filled','markeredgecolor','w'); hold on;
h = colorbar; colormap(cmocean('balance')); caxis([-1 1]);
plot(-castmeta_irminger6.lon(alignedind(i)), castmeta_irminger6.lat(alignedind(i)),'o','markerfacecolor','y','markeredgecolor','k','markersize',M1/2); hold on;
title(['Glider surfacings within 24 hours of Cast ' num2str(castmeta_irminger6.castnum(alignedind(i)))])
end
