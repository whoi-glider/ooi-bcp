%% Calculate deep isotherm salinity gain correction
sal_iso = mean(casts.ISO_SP(ind_ISOno2020)); %Salinity to pin to, from cruise casts interpolated on deep isotherm
iso_ind = find(pt_grid == ISO); %Index of deep isotherm in gridded wfp data

%Apply correction to data prior to creating merged data product to avoid seams across years
for yr = 1:10
    wgg{yr}.pracsal_gain = sal_iso./wgg{yr}.pracsal_ptgrid(iso_ind,:);
end

%% Plot to check correction
figure; clf
    subplot(411)
for yr = 1:10
    plot(wgg{yr}.time_start, wgg{yr}.doxy_lagcorr_ptgrid(iso_ind,:),'.'); hold on;
end
datetick('x'); ylabel('L1 oxygen, \mumol/kg'); title(['OOI Irminger WFP data on the ' num2str(ISO) ' \theta isotherm'])
ylim([280 350])
    subplot(412)
for yr = 1:10
    plot(wgg{yr}.time_start, wgg{yr}.pracsal_ptgrid(iso_ind,:),'.'); hold on;
    plot(wgg{yr}.time_start, wgg{yr}.pracsal_ptgrid(iso_ind,:).*wgg{yr}.pracsal_gain,'k.'); hold on;
end
datetick('x'); ylabel('Prac. salinity')
    subplot(414)
for yr = 1:10
    plot(wgg{yr}.time_start, wgg{yr}.pres_ptgrid(iso_ind,:),'.'); hold on;
end
datetick('x'); ylabel('Pressure, db')
    subplot(413)
for yr = 1:10
    plot(wgg{yr}.time_start, wgg{yr}.pracsal_gain,'k.'); hold on;
end
datetick('x'); 
ylabel('Salinity gain')

%% Apply salinity correction to all wfp data
for yr = 1:10
    [~,len] = size(wgg{yr}.mtime);
    % Calculate corrected practical salinity based on gain
    wgg{yr}.pracsal_corr = wgg{yr}.pracsal.*repmat(wgg{yr}.pracsal_gain, len, 1)';
    wgg{yr}.pracsal_grid_corr = wgg{yr}.pracsal_grid.*repmat(wgg{yr}.pracsal_gain, length(pres_grid_hypm), 1);
    wgg{yr}.pracsal_ptgrid_corr = wgg{yr}.pracsal_ptgrid.*repmat(wgg{yr}.pracsal_gain, length(pt_grid), 1);
    % Re-calculate derived properties using corrected salinity  
    [wgg{yr}.SA, ~] = gsw_SA_from_SP(wgg{yr}.pracsal_corr, wgg{yr}.pres, wgg{yr}.lon, wgg{yr}.lat); %absolute salinity from practical salinity - [SA, ~] = gsw_SA_from_SP(SP,p,long,lat)
    wgg{yr}.CT = gsw_CT_from_t(wgg{yr}.pracsal_corr, wgg{yr}.temp, wgg{yr}.pres); %Conservative Temperature from in-situ temperature - CT = gsw_CT_from_t(SA,t,p)
    wgg{yr}.pdens = gsw_rho(wgg{yr}.SA, wgg{yr}.CT, 0); %calculate potential density at reference pressure of 0 (surface)
end

