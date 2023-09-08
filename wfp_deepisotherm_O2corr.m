%% Recalculate oxygen using corrected salinity

for yr = 1:8
    wgg{yr}.doxy_lagcorr_salcorr_uM = aaoptode_salpresscorr(wgg{yr}.doxy_lagcorr, wgg{yr}.temp, wgg{yr}.pracsal_corr, wgg{yr}.pres, 0); %Oxygen concentration in uM
    wgg{yr}.doxy_lagcorr_salcorr_umolkg = wgg{yr}.doxy_lagcorr_salcorr_uM./(wgg{yr}.pdens/1000); %Divide by potential density to get oxygen in umol/kg
end

%% Redo oxygen and physical property gridding using corrected values

pres_grid = pres_grid_hypm;
S = 5; %points to smooth over

for yr = 1:8
    %Number of profile indices
    num_profiles = length(wgg{yr}.updown);
    %Calculate potential temperature
    wgg{yr}.ptemp = gsw_pt0_from_t(wgg{yr}.SA,wgg{yr}.temp,wgg{yr}.pres);
    for i = 1:num_profiles
        ind = find(~isnan(wgg{yr}.pres(i,:)) & ~isnan(wgg{yr}.doxy_lagcorr(i,:)) & wgg{yr}.flag(i,:) == 0); %no nan values for depth or oxygen and no range or spike flags
        wgg{yr}.doxy_lagcorr_grid(:,i) = movmean(interp1(wgg{yr}.pres(i,ind), wgg{yr}.doxy_lagcorr_salcorr_umolkg(i,ind), pres_grid),S);
        wgg{yr}.SA_grid(:,i) = movmean(interp1(wgg{yr}.pres(i,ind), wgg{yr}.SA(i,ind), pres_grid),S);
        wgg{yr}.CT_grid(:,i) = movmean(interp1(wgg{yr}.pres(i,ind), wgg{yr}.CT(i,ind), pres_grid),S);
        wgg{yr}.pracsal_grid(:,i) = movmean(interp1(wgg{yr}.pres(i,ind), wgg{yr}.pracsal_corr(i,ind), pres_grid),S);
        wgg{yr}.pdens_grid(:,i) = movmean(interp1(wgg{yr}.pres(i,ind), wgg{yr}.pdens(i,ind), pres_grid),S);
        ind = find(~isnan(wgg{yr}.ptemp(i,:)) & ~isnan(wgg{yr}.doxy_lagcorr(i,:)) & wgg{yr}.flag(i,:) == 0); %no nan values for depth or oxygen and no range or spike flags
        if length(ind) > 0
            wgg{yr}.doxy_lagcorr_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.doxy_lagcorr_salcorr_umolkg(i,ind), pt_grid),S);
            wgg{yr}.SA_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.SA(i,ind), pt_grid),S);
            wgg{yr}.pracsal_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.pracsal_corr(i,ind), pt_grid),S);
            wgg{yr}.pres_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.pres(i,ind), pt_grid),S);
            wgg{yr}.temp_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.temp(i,ind), pt_grid),S);
            wgg{yr}.pdens_ptgrid(:,i) = movmean(interp1(wgg{yr}.ptemp(i,ind), wgg{yr}.pdens(i,ind), pt_grid),S);
        else
            wgg{yr}.doxy_lagcorr_ptgrid(:,i) = NaN;
        end
    end
end

%% Calculate deep isotherm oxygen gain correction
oxy_iso = mean(casts.ISO_DOcorr_umolkg(ind_ISOno2020)); %Oxygen concentration to pin to, from cruise casts interpolated on deep isotherm
iso_ind = find(pt_grid == ISO); %Index of deep isotherm in gridded wfp data
smoothval = 5; %moving median over 5 profiles

%Apply correction to data prior to creating merged data product to avoid seams across years
for yr = 1:8
    wgg{yr}.oxy_gain = movmedian(oxy_iso./wgg{yr}.doxy_lagcorr_ptgrid(iso_ind,:),smoothval,'omitnan');
end

%% Plot to check correction
figure; clf
    subplot(411)
for yr = 1:8
    plot(wgg{yr}.time_start, wgg{yr}.doxy_lagcorr_ptgrid(iso_ind,:),'.'); hold on;
    plot(wgg{yr}.time_start, wgg{yr}.doxy_lagcorr_ptgrid(iso_ind,:).*wgg{yr}.oxy_gain,'k.'); hold on;
end
datetick('x'); ylabel('Oxygen, \mumol/kg'); title(['OOI Irminger WFP data on the ' num2str(ISO) ' \theta isotherm'])
ylim([230 290])
    subplot(413)
for yr = 1:8
    plot(wgg{yr}.time_start, wgg{yr}.pracsal_ptgrid(iso_ind,:),'.'); hold on;
end
datetick('x'); ylabel('Prac. salinity')
    subplot(414)
for yr = 1:8
    plot(wgg{yr}.time_start, wgg{yr}.pres_ptgrid(iso_ind,:),'.'); hold on;
end
datetick('x'); ylabel('Pressure, db')
    subplot(412)
for yr = 1:8
    plot(wgg{yr}.time_start, wgg{yr}.oxy_gain,'k.'); hold on;
end
datetick('x'); 
ylabel('Oxygen gain')
ylim([0.98 1.2])

%% Apply oxygen correction to all wfp data
for yr = 1:8
    [~,len] = size(wgg{yr}.mtime);
    % Calculate corrected oxygen concentration based on gain
    wgg{yr}.doxy_gaincorr = wgg{yr}.doxy_lagcorr_salcorr_umolkg.*repmat(wgg{yr}.oxy_gain, len, 1)';
    wgg{yr}.doxy_grid_gaincorr = wgg{yr}.doxy_lagcorr_grid.*repmat(wgg{yr}.oxy_gain, length(pres_grid_hypm), 1);
    wgg{yr}.doxy_ptgrid_gaincorr = wgg{yr}.doxy_lagcorr_ptgrid.*repmat(wgg{yr}.oxy_gain, length(pt_grid), 1);
end

