%% Recalculate oxygen using corrected salinity
%Also calculate optode-specific pressure compensation coefficient

figure(100); clf
C_yrs = cmocean('phase',9);
D_emp = NaN*ones(8,1);
num2align = [2 2 2 2 3 2 2 8]; %select earliest full depth WFP upcast in dataset to align with calcast (down for #7 to get deep enough)
calcasts = [5 13 28 41 49 75 96 116]; %hand selection casts for comparison
for yr = 1:8
    %Pull out cal cast closest to HYPM at deployment
    %yrind = find(casts.year == yr);
    %[~,hypmind] = min(casts.HYPMdist(yrind));
    %casttbl = castsum{casts.year(yrind(hypmind))}(casts.castnum(yrind(hypmind)));
    casttbl = castsum{casts.year(calcasts(yr))}(casts.castnum(calcasts(yr)));
    [wgg{yr}.doxy_lagcorr_salcorr_uM, D_emp(yr), optodecalcast{yr}.O2salcorr, minerrval(yr)] = aaoptode_salpresscorr_empiricalD(wgg{yr}.doxy_lagcorr, wgg{yr}.temp, wgg{yr}.pracsal_corr, wgg{yr}.pres, 0, ...
        casttbl{1}.prs, casttbl{1}.DOcorr_umolkg./(casttbl{1}.prho/1000), num2align(yr)); %Oxygen concentration in uM

%Plot calcast along with 1st cast in dataset
plot(casttbl{1}.DOcorr_umolkg./(casttbl{1}.prho/1000), casttbl{1}.prs, '-','linewidth',2,'color',C_yrs(yr,:)); hold on;

end
%close all

axis ij
ylabel('Pressure (db)')
xlabel('Oxygen (\muM)')
legend('2014','2015','2016','2017','2018','2019','2020','2021')
xlim([260 310])
title('Calibrated SBE43 cal-cast profiles for P compensation')

figure(101); clf
for yr = 1:8
    casttbl = castsum{casts.year(calcasts(yr))}(casts.castnum(calcasts(yr)));
    plot(optodecalcast{yr}.O2salcorr,casttbl{1}.prs,'-','linewidth',2,'color',C_yrs(yr,:)); hold on;
end
axis ij
ylabel('Pressure (db)')
xlabel('Oxygen (\muM)')
legend('2014','2015','2016','2017','2018','2019','2020','2021','location','SE')
%xlim([260 310])
title('S-corrected WFP optode cal-cast profiles prior to P compensation')

%Calculate overall mean pressure coefficient from years with good data
indgoodD = find(D_emp > 0.01);
D_overall = mean(D_emp(indgoodD));

for yr = 1:8
    if yr == 1 | yr == 4
        wgg{yr}.doxy_lagcorr_salcorr_uM = aaoptode_salpresscorr_fixedD(wgg{yr}.doxy_lagcorr, wgg{yr}.temp, wgg{yr}.pracsal_corr, wgg{yr}.pres, 0,D_overall);
        yr
    end
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
    %subplot(411)
for yr = 1:8
    plot(wgg{yr}.time_start, wgg{yr}.doxy_lagcorr_ptgrid(iso_ind,:),'.'); hold on;
    %plot(wgg{yr}.time_start, wgg{yr}.doxy_lagcorr_ptgrid(iso_ind,:).*wgg{yr}.oxy_gain,'k.'); hold on;
end
plot(casts.time(ind_ISOno2020), casts.ISO_DOcorr_umolkg(ind_ISOno2020),'k.','markersize',20); hold on;
datetick('x'); ylabel('Oxygen, \mumol/kg'); title(['Wire-following profiler L2-equivalent data on the ' num2str(ISO) ' \theta isotherm'])
xlim([datenum(2014,5,1) datenum(2022,10,1)])
ylim([230 290])
%     subplot(413)
% for yr = 1:8
%     plot(wgg{yr}.time_start, wgg{yr}.pracsal_ptgrid(iso_ind,:),'.'); hold on;
% end
% datetick('x'); ylabel('Prac. salinity')
%    subplot(414)
for yr = 1:8
    %plot(wgg{yr}.time_start, wgg{yr}.pres_ptgrid(iso_ind,:),'.'); hold on;
    if yr == 1
        A = wgg{yr}.pres_ptgrid(iso_ind,:);
    else
        A = [A, wgg{yr}.pres_ptgrid(iso_ind,:)];
    end
end
%datetick('x'); ylabel('Pressure, db')
%    subplot(412)
% for yr = 1:8
%     plot(wgg{yr}.time_start, wgg{yr}.oxy_gain,'k.'); hold on;
% end
% datetick('x'); 
% ylabel('Oxygen gain')
% ylim([0.98 1.2])

%% Apply oxygen correction to all wfp data
for yr = 1:8
    [~,len] = size(wgg{yr}.mtime);
    % Calculate corrected oxygen concentration based on gain
    wgg{yr}.doxy_gaincorr = wgg{yr}.doxy_lagcorr_salcorr_umolkg.*repmat(wgg{yr}.oxy_gain, len, 1)';
    wgg{yr}.doxy_grid_gaincorr = wgg{yr}.doxy_lagcorr_grid.*repmat(wgg{yr}.oxy_gain, length(pres_grid_hypm), 1);
    wgg{yr}.doxy_ptgrid_gaincorr = wgg{yr}.doxy_lagcorr_ptgrid.*repmat(wgg{yr}.oxy_gain, length(pt_grid), 1);
end

