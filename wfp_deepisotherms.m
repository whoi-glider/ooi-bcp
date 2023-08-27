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

ISO = 3.1; %isotherm to plot
iso_ind = find(pt_grid == ISO);

figure(8); clf
    subplot(311)
for yr = 1:8
    yrind = find(wggmerge.deploy_yr == yr);
    plot(wggmerge.time(yrind), wggmerge.doxy_lagcorr_pt(iso_ind,yrind),'.'); hold on;
    try
        plot(castalign{yr}.time, castalign{yr}.SBE_interp(iso_ind)./castalign{yr}.prho_interp(iso_ind)*1000,'ko','markerfacecolor','k'); hold on;
    end
end
datetick('x'); ylabel('L2 oxygen, \mumol/kg'); title(['OOI Irminger WFP data on the ' num2str(ISO) ' \theta isotherm'])
ylim([230 290])
    subplot(312)
for yr = 1:8
    yrind = find(wggmerge.deploy_yr == yr);
    plot(wggmerge.time(yrind), wggmerge.pracsal_pt(iso_ind,yrind),'.'); hold on;
    try
        plot(castalign{yr}.time, castalign{yr}.SP_interp(iso_ind),'ko','markerfacecolor','k'); hold on;
    end
end
datetick('x'); ylabel('Prac salinity, uncorr')
    subplot(313)
for yr = 1:8
    yrind = find(wggmerge.deploy_yr == yr);
    plot(wggmerge.time(yrind), wggmerge.pres_pt(iso_ind,yrind),'.'); hold on;
end
datetick('x'); ylabel('Pressure, db')

%%
figure(9); clf
for yr = [1:9]
    yrind = find(wggmerge.deploy_yr == yr);
    if yr == 4 | yr > 7
        plot(wggmerge.time(yrind), wggmerge.doxy_lagcorr_pt(iso_ind,yrind).*gain_yr_pt(yr,1),'.'); hold on;
    else
        plot(wggmerge.time(yrind), wggmerge.doxy_lagcorr_pt(iso_ind,yrind).*gain_yr_pt(yr,3),'.'); hold on;
    end
    try
        plot(castalign{yr}.time, castalign{yr}.SBE_interp(iso_ind)./castalign{yr}.prho_interp(iso_ind)*1000,'ko','markerfacecolor','k'); hold on;
    end
end
ylim([240 300])
datetick('x',2); ylabel('Gain-corrected L2 oxygen, \mumol/kg'); title(['OOI Irminger WFP data on the ' num2str(ISO) ' \theta isotherm'])

%% Apply a rough deep isotherm correction

iso_test = [2.4, 2.6, 2.8, 3.0, 3.1];
C_gl = cmocean('phase',length(iso_test)+2);
clear leg h
mksz = 8;
ftsz = 12;

figure(10); clf
for i = 1:length(iso_test)
    %Determine value of deep isotherm for pinning from cruise casts
    ISO = iso_test(i); %isotherm to plot
    iso_ind = find(pt_grid == ISO);
    deep_iso_pin = NaN(9,1);
    leg{i} = num2str(iso_test(i));
    for yr = 1:9
        try
            deep_iso_pin(yr) = castalign{yr}.SBE_interp(iso_ind)./castalign{yr}.prho_interp(iso_ind)*1000;
        end
    end
    
    %Calculate profile-specific gain correction to pin to deep isotherm
    wfp_profilegain_nosmooth_plot = nanmedian(deep_iso_pin)./wggmerge.doxy_lagcorr_pt(iso_ind,:);
    wfp_profilegain_plot = movmean(nanmedian(deep_iso_pin)./wggmerge.doxy_lagcorr_pt(iso_ind,:),30,'omitnan'); %30 profile smoothing
    %try
        %plot(wggmerge.time, wfp_profilegain_nosmooth, '.','color',C_gl(i,:)); hold on;
        h(i) = plot(wggmerge.time, wfp_profilegain_plot, '.', 'color',C_gl(i,:),'markersize',mksz); hold on;
    %end
    if ISO == 3.1
        %Interpolate onto times with missing values - use 'previous'
        ind = find(isnan(wfp_profilegain_plot) == 0);
        wfp_profilegain_interp = interp1(wggmerge.time(ind), wfp_profilegain_plot(ind), wggmerge.time, 'previous');
        wfp_profilegain = wfp_profilegain_plot;
    end
    clear wfp_profilegain_nosmooth_plot wfp_profilegain_plot
end

    xlim([min(wggmerge.time)-30 max(wggmerge.time)-100])
    datetick('x',2,'keeplimits')
    title('OOI Irminger WFP oxygen, deep isotherm gain correction (30-profile smoothing across range of isotherms)','fontsize',ftsz)
    %legend('Individual profile values','30-profile moving mean')
    legend(h,leg,'fontsize',ftsz)
    ylabel('Gain')

%% Pull all data on selected deep isotherm across all casts, not just HYPM-aligned
ISO = 3.1;
%Loop over all casts in castsum to find values interpolated at pt == ISO
%Put output values in table "casts"

for i = 1:height(casts)
    casttbl = castsum{casts.year(i)}(casts.castnum(i));
    casts.ISOdepth(i) = interp1(casttbl{1}.pt, casttbl{1}.depth, ISO);
    casts.ISO_SP(i) = interp1(casttbl{1}.pt, casttbl{1}.SP, ISO);
    casts.ISO_prho(i) = interp1(casttbl{1}.pt, casttbl{1}.prho, ISO);
    casts.ISO_DOcorr_umolkg(i) = interp1(casttbl{1}.pt, casttbl{1}.DOcorr_umolkg, ISO);
end
%%
%Identify all values within +/- 2 stdev of mean for depth, salinity, and density
ind_ISO_depth = find(casts.ISOdepth < nanmean(casts.ISOdepth) + 2*nanstd(casts.ISOdepth) & casts.ISOdepth > nanmean(casts.ISOdepth) - 2*nanstd(casts.ISOdepth));
ind_ISO_SP = find(casts.ISO_SP < nanmean(casts.ISO_SP) + 2*nanstd(casts.ISO_SP) & casts.ISO_SP > nanmean(casts.ISO_SP) - 2*nanstd(casts.ISO_SP));
ind_ISO_prho = find(casts.ISO_prho < nanmean(casts.ISO_prho) + 2*nanstd(casts.ISO_prho) & casts.ISO_prho > nanmean(casts.ISO_prho) - 2*nanstd(casts.ISO_prho));
ind_phys = intersect(intersect(ind_ISO_depth, ind_ISO_SP),ind_ISO_prho);

%Within the non-outlier values for physical properties, identify DO values within +/- 2 stdev of mean
ind_ISO_DO = find(casts.ISO_DOcorr_umolkg(ind_phys) < nanmean(casts.ISO_DOcorr_umolkg(ind_phys)) + 2*nanstd(casts.ISO_DOcorr_umolkg(ind_phys)) &...
    casts.ISO_DOcorr_umolkg(ind_phys) > nanmean(casts.ISO_DOcorr_umolkg(ind_phys)) - 2*nanstd(casts.ISO_DOcorr_umolkg(ind_phys)));

%Final index with outliers removed
ind_ISO = ind_phys(ind_ISO_DO);

%% Plot variability of data on deep isotherm
ind_nearHYPM = find(casts.HYPMdist < 50);
ind_ISOnearHYPM = intersect(ind_ISO, ind_nearHYPM); %doesn't seem to make a difference, tried values from 50-200 for distance
ind_no2020 = find(casts.time < datenum(2020,1,1) | casts.time > datenum(2021,1,1));
ind_ISOno2020 = intersect(ind_ISO, ind_no2020);
ind_ISOno2020nearHYPM = intersect(ind_ISOno2020, ind_ISOnearHYPM); %still doesn't seem to make a difference
ind_noDepthOutliers = find(casts.ISOdepth > 1800 & casts.ISOdepth < 2050);
ind_ISOno2020noDepthOutliers = intersect(ind_noDepthOutliers, ind_ISOno2020); %doesn't seem to make a difference
ind_noSalOutliers = find(casts.ISO_SP < (mean(casts.ISO_SP(ind_ISO))+ std(casts.ISO_SP(ind_ISO))) &...
    casts.ISO_SP > (mean(casts.ISO_SP(ind_ISO))- std(casts.ISO_SP(ind_ISO))));
ind_ISOno2020noSalOutliers = intersect(ind_noSalOutliers, ind_ISOno2020);

figure(11); clf
    subplot(4,1,1)
histogram(casts.ISOdepth(ind_ISO))
    subplot(4,1,2)
histogram(casts.ISO_SP(ind_ISO))
    subplot(4,1,3)
histogram(casts.ISO_prho(ind_ISO))
    subplot(4,1,4)
histogram(casts.ISO_DOcorr_umolkg(ind_ISO))

    M = 16;
figure(12); clf
%scatter(casts.time(ind_ISOno2020), casts.ISO_DOcorr_umolkg(ind_ISOno2020),[],casts.ISO_SP(ind_ISOno2020),'filled');
plot(casts.time(ind_ISOno2020noSalOutliers), casts.ISO_DOcorr_umolkg(ind_ISOno2020noSalOutliers),'k.','markersize',M);
text(datenum(2014,3,1), 278.5, ['mean = ' num2str(mean(casts.ISO_DOcorr_umolkg(ind_ISOno2020noSalOutliers)),4) char(177) num2str(std(casts.ISO_DOcorr_umolkg(ind_ISOno2020noSalOutliers)),2) ' \mumol/kg'])
datetick('x','keeplimits')
ylabel('\mumol/kg')
ylim([268 281])
title('Dissolved oxygen on 3.1 deg isotherm from calibrated oxygen profiles on cruise casts')
%c = colorbar;
%ylabel(c, 'Salinity')

%% Uncertainty quantification: how well do we know the deep isotherm?
iso_percenterr = std(casts.ISO_DOcorr_umolkg(ind_ISOno2020noSalOutliers))./mean(casts.ISO_DOcorr_umolkg(ind_ISOno2020noSalOutliers))*100;
iso_rmse = sqrt(mean((casts.ISO_DOcorr_umolkg(ind_ISOno2020noSalOutliers) - mean(casts.ISO_DOcorr_umolkg(ind_ISOno2020noSalOutliers))).^2));
