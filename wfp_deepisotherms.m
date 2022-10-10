%% Reproduce Izi's plot of # measurements at each isotherm

for yr = 1:8
    yrind = find(wggmerge.deploy_yr == yr);
    for j = 1:length(pt_grid)
        iso_sum(yr,j) = sum(~isnan(wggmerge.doxy_lagcorr_pt(j,yrind)));
    end
end

figure(7); clf
M = 10;
for yr = 1:8
    plot(pt_grid,iso_sum(yr,:),'.','markersize',M); hold on;
end
xlabel('Potential temperature')
ylabel('Number of measurements')
legend('Year 1: 2014-2015','Year 2: 2015-2016','Year 3: 2016-2017','Year 4: 2017-2018',...
    'Year 5: 2018-2019','Year 6: 2019-2020','Year 7: 2020-2021','Year 8: 2021-2022')

%% Plot example along a deep isotherm
%Use merged record gridded on isotherms

ISO = 2.5; %isotherm to plot
iso_ind = find(pt_grid == ISO);

figure(8); clf
    subplot(311)
for yr = 1:8
    yrind = find(wggmerge.deploy_yr == yr)
    plot(wggmerge.time(yrind), wggmerge.doxy_lagcorr_pt(iso_ind,yrind),'.'); hold on;
end
datetick('x'); ylabel('L2 oxygen, \mumol/kg'); title(['OOI Irminger WFP data on the ' num2str(ISO) ' \theta isotherm'])
    subplot(312)
for yr = 1:8
    yrind = find(wggmerge.deploy_yr == yr)
    plot(wggmerge.time(yrind), wggmerge.pracsal_pt(iso_ind,yrind),'.'); hold on;
end
datetick('x'); ylabel('Prac salinity, uncorr')
    subplot(313)
for yr = 1:8
    yrind = find(wggmerge.deploy_yr == yr)
    plot(wggmerge.time(yrind), wggmerge.pres_pt(iso_ind,yrind),'.'); hold on;
end
datetick('x'); ylabel('Pressure, db')

