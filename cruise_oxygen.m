%% Read in corrected cruise oxygen data
% CTD SBE43 data are corrected by Kristen using Winklers

mltoumol = 44.661;

%% Read in processed data and reformat
addpath(genpath('G:\Shared drives\NSF_Irminger\OOI Cruises CTD Casts\CTD_Data\Processed'))
load('AllYears_Processed_KF.mat')
clear castlist

load Year1_Processed_KF.mat
Winkler_umolkg = [btlsum_tbl_yr1.Winkler1_umolkg; btlsum_tbl_yr1.Winkler2_umolkg];
btlsum_tbl_yr1 = [btlsum_tbl_yr1; btlsum_tbl_yr1];
btlsum_tbl_yr1.Winkler_umolkg = Winkler_umolkg;
castlist{1} = [2:9];

load Year2_Processed_KF.mat
castlist{2} = [1:13];

load Year3_Processed_KF.mat
castlist{3} = [1:10];

load Year4_Processed_KF.mat
castlist{4} = [9:12];

load Year5_Processed_KF.mat
Winkler_umolkg = [btlsum_tbl_yr5.Winkler_OOI_umolkg; btlsum_tbl_yr5.Winkler1_HIP_umolkg; btlsum_tbl_yr5.Winkler2_HIP_umolkg];
btlsum_tbl_yr5 = [btlsum_tbl_yr5; btlsum_tbl_yr5; btlsum_tbl_yr5];
btlsum_tbl_yr5.Winkler_umolkg = Winkler_umolkg;
castlist{5} = [1:23];

castlist{6} = cast_num{6};
castlist{7} = cast_num{7}; 

load Year8_Processed_KF.mat
castlist{8} = [1:12];

load Year9_Processed_KF.mat
Winkler_umolkg = [btlsum_tbl_yr9.Winkler1_umolkg; btlsum_tbl_yr9.Winkler2_umolkg; btlsum_tbl_yr9.Winkler3_umolkg];
btlsum_tbl_yr9 = [btlsum_tbl_yr9; btlsum_tbl_yr9; btlsum_tbl_yr9];
btlsum_tbl_yr9.Winkler_umolkg = Winkler_umolkg;
castlist{9} = [6:23];

%Merge bottle summary tables
btlsum{1} = btlsum_tbl_yr1;
btlsum{2} = btlsum_tbl_yr2;
btlsum{3} = btlsum_tbl_yr3;
btlsum{4} = btlsum_tbl_yr4;
btlsum{5} = btlsum_tbl_yr5;
btlsum{8} = btlsum_tbl_yr8;
btlsum{9} = btlsum_tbl_yr9;

castsum = downcasts;  

%Over-write with old file output - eventually update to use finalized data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\From_Hilary')
load cruise_oxygen_output.mat

%% Loop over data from all years
for yr = [1:5,8,9]
    castsumyr = castsum{yr};
    castlistyr = castlist{yr};
    
    %% Plot all casts and interpolate at Winkler depth
    figure(100*yr); clf
    for i = 1:length(castlistyr)
        if castlistyr(i) < 17 %set for subplot max
        subplot(4,4,i)
            btlid = find(btlsum{yr}.Cast == i);
        plot(castsumyr{castlistyr(i)}.DOcorr_umolkg, castsumyr{castlistyr(i)}.prs, 'k-'); hold on;
        if length(btlid) > 0
            plot(btlsum{yr}.Winkler_umolkg(btlid), btlsum{yr}.prs(btlid), 'ro'); hold on;
        end
        axis ij
        xlabel('SBE43 calibrated O_2, \mumol/kg')
        ylabel('Pressure (db)')
        end
    end

    %% Check goodness of calibration

    %Regression plot
    figure(100*yr + 1); clf
    scatter(btlsum{yr}.Winkler_umolkg, btlsum{yr}.DOcorr_umolkg, [], btlsum{yr}.prs, 'filled'); hold on;
    plot([140:350],[140:350],'k--')
    c = colorbar; ylabel(c, 'Depth (m)')
    axis([140 345 140 345])
    xlabel('Winkler oxygen, \mumol/L'); ylabel('Calibrated SBE43 oxygen, \mumol/kg')

    %Calculate uncertainty from residuals
    U = nanmean(abs(btlsum{yr}.Winkler_umolkg - btlsum{yr}.DOcorr_umolkg));
        deep_cutoff = 200; %subset deep values only
        deep = find(btlsum{yr}.prs > deep_cutoff);
    U_deep = nanmean(abs(btlsum{yr}.Winkler_umolkg(deep) - btlsum{yr}.DOcorr_umolkg(deep)));  

    %Plot residuals
    figure(100*yr + 2); clf
    plot(btlsum{yr}.prs, (btlsum{yr}.Winkler_umolkg - btlsum{yr}.DOcorr_umolkg), 'k.'); hold on;
    plot([0 3000], 0*[0 3000], 'k--')
    ylabel('Residual, Winkler - SBE O_2, \mumol/kg')
    xlabel('Pressure (db)')
    ylim([-10 10])
    text(500, 4, ['Mean residual = ' num2str(U) ' \mumol/kg'])
    text(500, 3.4, ['Mean residual > ' num2str(deep_cutoff) 'm = ' num2str(U_deep) ' \mumol/kg'])

end

%% Create table to hold information about casts and plot locations

load OOImooringLocations.mat %lat-lon coordinates of each year's moorings
clear casts
casts = table;
casts.year = [1*ones(length(castlist{1}),1); 2*ones(length(castlist{2}),1); 3*ones(length(castlist{3}),1);...
    4*ones(length(castlist{4}),1); 5*ones(length(castlist{5}),1); 8*ones(length(castlist{8}),1); 9*ones(length(castlist{9}),1)];
casts.castnum = [castlist{1}'; castlist{2}'; castlist{3}'; castlist{4}'; castlist{5}'; castlist{8}'; castlist{9}'];
casts.lat = NaN*ones(length(casts.year),1);
casts.lon = NaN*ones(length(casts.year),1);
casts.time = NaN*ones(length(casts.year),1);
casts.numWinkl = NaN*ones(length(casts.year),1);
casts.HYPMdist = NaN*ones(length(casts.year),1);


%for each year, pull out each cast and plot lat/lon
figure(10); clf
M = 20;
Ftsz = 12;
colorlist = [nicecolor('rrry'); nicecolor('ryyym'); nicecolor('yyyyyym'); nicecolor('ggyc'); nicecolor('bcc'); nicecolor('c'); nicecolor('cbbbb'); nicecolor('mmmrc'); nicecolor('rrrrm')];

h5 = plot(-OOImoorings.SUMO5(2), OOImoorings.SUMO5(1),'sk','markersize',M,'markerfacecolor',nicecolor('ggk')); hold on;
h6 = plot(-OOImoorings.HYPM5(2), OOImoorings.HYPM5(1),'dk','markersize',M,'markerfacecolor',nicecolor('rb')); hold on;
h7 = plot(-OOImoorings.FLMA5(2), OOImoorings.FLMA5(1),'ok','markersize',M,'markerfacecolor',nicecolor('bw')); hold on;
h8 = plot(-OOImoorings.FLMB5(2), OOImoorings.FLMB5(1),'ok','markersize',M,'markerfacecolor',nicecolor('bbbkkw')); hold on;

for yr = [1:5,8,9]
    castsumyr = castsum{yr};
    castlistyr = castlist{yr};
    ind = find(casts.year == yr);
    for i = castlistyr
        lat(i) = nanmean(castsumyr{i}.lat);
        lon(i) = nanmean(castsumyr{i}.lon);
        numWinkl(i) = length(find(btlsum{yr}.Cast == i & isnan(btlsum{yr}.Winkler_umolkg) == 0));
        HYPMdist(i) = distlatlon(OOImoorings.HYPM5(1), lat(i), OOImoorings.HYPM5(2), lon(i)); %need to update with locations for each year
        time(i) = datenum(castsumyr{i}.StartTimeUTC(1));
    end
    casts.lat(ind) = lat(castlistyr); 
    casts.lon(ind) = lon(castlistyr);
    casts.time(ind) = time(castlistyr);
    casts.numWinkl(ind) = numWinkl(castlistyr);
    casts.HYPMdist(ind) = HYPMdist(castlistyr);
    leg{yr} = plot(-casts.lon(ind), casts.lat(ind), 'ko','markersize',M/2-yr/3,'markerfacecolor',colorlist(yr,:)); hold on;
end

xlim([39 40.1])
ylim([59.6 60.05])
xlabel('Longitude (^oW)','Fontsize', Ftsz)
ylabel('Latitude (^oN)','Fontsize', Ftsz)
box on
legend([h5 h6 h7 h8 leg{1} leg{2} leg{3} leg{4} leg{5} leg{8} leg{9}],...
    'Apex surface mooring', 'Apex sub-surface profiler mooring', 'Flanking mooring A', 'Flanking mooring B',...
    'Year 1 casts', 'Year 2 casts', 'Year 3 casts', 'Year 4 casts','Year 5 casts','Year 8 casts','Year 9 casts',...
    'location','northeastoutside','Fontsize', Ftsz);

%% Extract cruise SBE43 and Winkler data aligned with HYPM

tol_dist = 10; %Extract all casts within X km of HYPM
tol_time = 2; %Extract all HYPM profiles within X days of cast
ind = find(casts.HYPMdist < tol_dist & casts.numWinkl > 0);

figure(1); clf
for i = 1:length(ind)
    subplot(5,3,i)
    castyr = castsum{casts.year(ind(i))};
    %Interpolate the SBE43 profile onto the HYPM depth grid
    SBE_interp = interp1(castyr{casts.castnum(ind(i))}.prs, castyr{casts.castnum(ind(i))}.DOcorr_umolkg.*castyr{casts.castnum(ind(i))}.prho/1000, pres_grid); %umol/L
    %Extract the corresponding HYPM data
    indHYPM = find(abs(wggmerge.time - casts.time(ind(i))) < tol_time);
    for j = 1:length(indHYPM)
        deploy_yr = wggmerge.deploy_yr(indHYPM(j));
        HYPMprofile_O2umolL = wggmerge.doxy_lagcorr(:,indHYPM(j)).*wggmerge.pdens(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
        plot(HYPMprofile_O2umolL, pres_grid, '-', 'color', colorlist(deploy_yr,:)*0.6); hold on; 
        %Calculate gain correction
        gain = SBE_interp'./HYPMprofile_O2umolL;
        gainmeans(i,j) = nanmean(gain);
        gainstds(i,j) = nanstd(gain);
    end
    gain_cast(i) = mean(gainmeans(i,1:length(indHYPM)));
    %Replot with the corrected data
    for j = 1:length(indHYPM)
        deploy_yr = wggmerge.deploy_yr(indHYPM(j));
        HYPMprofile_O2umolL = wggmerge.doxy_lagcorr(:,indHYPM(j)).*wggmerge.pdens(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
        plot(HYPMprofile_O2umolL*gain_cast(i), pres_grid, '-', 'color', colorlist(deploy_yr,:)); hold on; 
    end
    %Plot the SBE43 profile and Winkler data
    castyr = castsum{casts.year(ind(i))};
    plot(castyr{casts.castnum(ind(i))}.DOcorr_umolkg.*castyr{casts.castnum(ind(i))}.prho/1000, castyr{casts.castnum(ind(i))}.prs, 'k-','linewidth',1.5); hold on;
    if casts.numWinkl(ind(i)) > 0
        btlsumyr = btlsum{casts.year(ind(i))};
        indbtl = find(btlsumyr.Cast == casts.castnum(ind(i)));
        plot(btlsumyr.Winkler_mLL(indbtl)*mltoumol, btlsumyr.prs(indbtl), 'bo'); hold on;
    end
    %Plot cleanup
    axis ij
    axis([240 320 0 2700])
    title(['Year ' num2str(casts.year(ind(i))) ', Cast ' num2str(casts.castnum(ind(i))) ', ' datestr(casts.time(ind(i)),1)])
    xlabel('Oxygen, \mumol/L','Fontsize',8)
    ylabel('Pressure (db)','Fontsize',8)
end

%% Approach of applying best initial gain correction for each year - plotted via pressure
tol_time = 2;
minnumpts = 300; %number of points in HYPM profile if usable for comparison

figure(100); clf
for yr = [1,2,3,5,8,9];
    
subplot(3,3,yr)

if yr == 1
    cast = 6;
elseif yr == 2
    cast = 12;
elseif yr == 3
    cast = 7;
elseif yr == 5
    cast = 6;
elseif yr == 8
    cast = 5;
elseif yr == 9
    cast = 13;
end

castyr = castsum{yr};
ind = find(casts.year == yr & casts.castnum == cast);

    %Interpolate the SBE43 profile onto the HYPM depth grid
    SBE_interp = interp1(castyr{casts.castnum(ind)}.prs, castyr{casts.castnum(ind)}.DOcorr_umolkg.*castyr{casts.castnum(ind)}.prho/1000, pres_grid);
    %Extract the corresponding HYPM data
    if yr < 6
        indHYPM = find(abs(wggmerge.time - casts.time(ind)) < tol_time & wggmerge.deploy_yr == yr);
    else
        indHYPM = find(abs(wggmerge.time - casts.time(ind)) < tol_time & wggmerge.deploy_yr == yr-1);
    end
    for j = 1:length(indHYPM)
        if sum(~isnan(wggmerge.doxy_lagcorr(:,indHYPM(j)))) > minnumpts
            HYPMprofile_O2umolL = wggmerge.doxy_lagcorr(:,indHYPM(j)).*wggmerge.pdens(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
            plot(HYPMprofile_O2umolL, pres_grid, '-', 'color', colorlist(yr,:)*0.6); hold on; 
            %Calculate gain correction
            gain = SBE_interp'./HYPMprofile_O2umolL;
            gainm(j) = nanmean(gain);
            gains(j) = nanstd(gain);
        else
            gainm(j) = NaN; gains(j) = NaN;
        end
    end
    gain_yr(yr,1) = nanmean(gainm);
    gain_yr(yr,2) = nanstd(gainm);
    clear gainm gains
    
    %Replot with the corrected data
    for j = 1:length(indHYPM)
        HYPMprofile_O2umolL = wggmerge.doxy_lagcorr(:,indHYPM(j)).*wggmerge.pdens(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
        plot(HYPMprofile_O2umolL*gain_yr(yr,1), pres_grid, '-', 'color', colorlist(yr,:)); hold on; 
    end
    %Plot the SBE43 profile and Winkler data
    plot(castyr{casts.castnum(ind)}.DOcorr_umolkg.*castyr{casts.castnum(ind)}.prho/1000, castyr{casts.castnum(ind)}.prs, 'k-','linewidth',1.5); hold on;
    if casts.numWinkl(ind) > 0
        btlsumyr = btlsum{casts.year(ind)};
        indbtl = find(btlsumyr.Cast == casts.castnum(ind));
        plot(btlsumyr.Winkler_mLL(indbtl)*mltoumol, btlsumyr.prs(indbtl), 'bo','markerfacecolor','b'); hold on;
    end
    %Plot cleanup
    axis ij
    axis([250 320 0 2700])
    title(['Year ' num2str(yr) ', Cast ' num2str(cast) ', ' datestr(casts.time(ind),1)])
    xlabel('Oxygen, \mumol/L','Fontsize',8)
    ylabel('Pressure (db)','Fontsize',8)

end

%Plot year 4 as a special case with no SBE data
yr = 4; cast = 8;
subplot(3,3,yr)

%Extract bottle IDs of Winkler samples
indbtl = find(btlsum{4}.Cast == cast);

%Extract the corresponding HYPM data
    indHYPM = find(abs(wggmerge.time - datenum(btlsum{4}.Date(indbtl(1))) < tol_time) & wggmerge.deploy_yr == yr);
    for j = 1:length(indHYPM)
        if sum(~isnan(wggmerge.doxy_lagcorr(:,indHYPM(j)))) > minnumpts
            HYPMprofile_O2umolL = wggmerge.doxy_lagcorr(:,indHYPM(j)).*wggmerge.pdens(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
            plot(HYPMprofile_O2umolL, pres_grid, '-', 'color', colorlist(yr,:)*0.6); hold on; 
            %Calculate gain correction
            HYPM_interp = interp1(pres_grid, HYPMprofile_O2umolL, btlsum{4}.prs(indbtl));
            gain = btlsum{4}.Winkler_mLL(indbtl)*mltoumol./HYPM_interp;
            gainm(j) = nanmean(gain);
            gains(j) = nanstd(gain);
        else
            gainm(j) = NaN; gains(j) = NaN;
        end
    end
    gain_yr(yr,1) = nanmean(gainm);
    gain_yr(yr,2) = nanstd(gainm);
    clear gainm gains
    
    %Replot with the corrected data
    for j = 1:length(indHYPM)
        HYPMprofile_O2umolL = wggmerge.doxy_lagcorr(:,indHYPM(j)).*wggmerge.pdens(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
        plot(HYPMprofile_O2umolL*gain_yr(yr,1), pres_grid, '-', 'color', colorlist(yr,:)); hold on; 
    end
    
    %Plot the Winkler data
    plot(btlsum{4}.Winkler_mLL(indbtl)*mltoumol, btlsum{4}.prs(indbtl), 'bo','markerfacecolor','b'); hold on;

    %Plot cleanup
    axis ij
    axis([250 320 0 2700])
    title(['Year ' num2str(yr) ', Cast ' num2str(cast) ', ' datestr(btlsum{4}.Date(indbtl(1)),1)])
    xlabel('Oxygen, \mumol/L','Fontsize',8)
    ylabel('Pressure (db)','Fontsize',8)
    
%% Approach of applying best initial gain correction for each year - plotted via temperature
tol_time = 2;
minnumpts = 300; %number of points in HYPM profile if usable for comparison

figure(101); clf
for yr = [1,2,3,5,8,9];
    
subplot(3,3,yr)

if yr == 1
    cast = 6;
elseif yr == 2
    cast = 12;
elseif yr == 3
    cast = 7;
elseif yr == 5
    cast = 6;
elseif yr == 8
    cast = 5;
elseif yr == 9
    cast = 13;
end

castyr = castsum{yr};
ind = find(casts.year == yr & casts.castnum == cast);

%Subset isotherms to fit gain specifically to deep depths
deep_ind = find(pt_grid >= 2 & pt_grid <= 3);

%Pull out relevant cast information to save for later use
castalign{yr}.castnum = cast;
idc = find(casts.year == yr & casts.castnum == cast);
castalign{yr}.time = casts.time(idc);

    %Interpolate the SBE43 profile onto the HYPM depth grid
    castalign{yr}.SBE_interp = interp1(castyr{casts.castnum(ind)}.pt, castyr{casts.castnum(ind)}.DOcorr_umolkg.*castyr{casts.castnum(ind)}.prho/1000, pt_grid);
    castalign{yr}.SP_interp = interp1(castyr{casts.castnum(ind)}.pt, castyr{casts.castnum(ind)}.SP, pt_grid);
    castalign{yr}.SA_interp = interp1(castyr{casts.castnum(ind)}.pt, castyr{casts.castnum(ind)}.SA, pt_grid);
    castalign{yr}.prho_interp = interp1(castyr{casts.castnum(ind)}.pt, castyr{casts.castnum(ind)}.prho, pt_grid);
    %Extract the corresponding HYPM data
    if yr < 6
        indHYPM = find(abs(wggmerge.time - casts.time(ind)) < tol_time & wggmerge.deploy_yr == yr);
    else
        indHYPM = find(abs(wggmerge.time - casts.time(ind)) < tol_time & wggmerge.deploy_yr == yr-1);
    end
    for j = 1:length(indHYPM)
        if sum(~isnan(wggmerge.doxy_lagcorr(:,indHYPM(j)))) > minnumpts
            HYPMprofile_O2umolL = wggmerge.doxy_lagcorr_pt(:,indHYPM(j)).*wggmerge.pdens_pt(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
            plot(HYPMprofile_O2umolL, pt_grid, '-', 'color', colorlist(yr,:)*0.6); hold on; 
            %Calculate gain correction
            gain = castalign{yr}.SBE_interp'./HYPMprofile_O2umolL;
            gain_deep = castalign{yr}.SBE_interp(deep_ind)'./HYPMprofile_O2umolL(deep_ind);
            gainm(j) = nanmean(gain);
            gains(j) = nanstd(gain);
            gainmd(j) = nanmean(gain_deep);
            gainsd(j) = nanstd(gain_deep);
        else
            gainm(j) = NaN; gains(j) = NaN;
        end
    end
    gain_yr_pt(yr,1) = nanmean(gainm);
    gain_yr_pt(yr,2) = nanstd(gainm);
    gain_yr_pt(yr,3) = nanmean(gainmd);
    gain_yr_pt(yr,4) = nanstd(gainmd);
    %clear gainm gains
    
    %Replot with the corrected data
    for j = 1:length(indHYPM)
        HYPMprofile_O2umolL = wggmerge.doxy_lagcorr_pt(:,indHYPM(j)).*wggmerge.pdens_pt(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
        plot(HYPMprofile_O2umolL*gain_yr_pt(yr,1), pt_grid, '-', 'color', colorlist(yr,:)); hold on; 
    end
    %Plot the SBE43 profile and Winkler data
    plot(castyr{casts.castnum(ind)}.DOcorr_umolkg.*castyr{casts.castnum(ind)}.prho/1000, castyr{casts.castnum(ind)}.pt, 'k-','linewidth',1.5); hold on;
    if casts.numWinkl(ind) > 0
        btlsumyr = btlsum{casts.year(ind)};
        indbtl = find(btlsumyr.Cast == casts.castnum(ind));
        plot(btlsumyr.Winkler_mLL(indbtl)*mltoumol, btlsumyr.pt(indbtl), 'bo','markerfacecolor','b'); hold on;
    end
    %Plot cleanup
    axis([250 320 1.5 5])
    title(['Year ' num2str(yr) ', Cast ' num2str(cast) ', ' datestr(casts.time(ind),1)])
    xlabel('Oxygen, \mumol/L','Fontsize',8)
    ylabel('Potential temp','Fontsize',8)

end

% %Plot year 4 as a special case with no SBE data
yr = 4; cast = 8;
subplot(3,3,yr)

%Extract bottle IDs of Winkler samples
indbtl = find(btlsum{4}.Cast == cast);

%Extract the corresponding HYPM data
    indHYPM = find(abs(wggmerge.time - datenum(btlsum{4}.Date(indbtl(1))) < tol_time) & wggmerge.deploy_yr == yr);
    for j = 1:length(indHYPM)
        if sum(~isnan(wggmerge.doxy_lagcorr(:,indHYPM(j)))) > minnumpts
            HYPMprofile_O2umolL = wggmerge.doxy_lagcorr_pt(:,indHYPM(j)).*wggmerge.pdens_pt(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
            plot(HYPMprofile_O2umolL, pt_grid, '-', 'color', colorlist(yr,:)*0.6); hold on; 
            %Calculate gain correction
            HYPM_interp = interp1(pt_grid, HYPMprofile_O2umolL, btlsum{4}.pt(indbtl));
            gain = btlsum{4}.Winkler_mLL(indbtl)*mltoumol./HYPM_interp;
            gainm(j) = nanmean(gain);
            gains(j) = nanstd(gain);
        else
            gainm(j) = NaN; gains(j) = NaN;
        end
    end
    gain_yr_pt(yr,1) = nanmean(gainm);
    gain_yr_pt(yr,2) = nanstd(gainm);
    clear gainm gains
    
    %Replot with the corrected data
    for j = 1:length(indHYPM)
        HYPMprofile_O2umolL = wggmerge.doxy_lagcorr_pt(:,indHYPM(j)).*wggmerge.pdens_pt(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
        plot(HYPMprofile_O2umolL*gain_yr_pt(yr,1), pt_grid, '-', 'color', colorlist(yr,:)); hold on; 
    end
    
    %Plot the Winkler data
    plot(btlsum{4}.Winkler_mLL(indbtl)*mltoumol, btlsum{4}.pt(indbtl), 'bo','markerfacecolor','b'); hold on;

    %Plot cleanup
    axis([250 320 1.5 5])
    title(['Year ' num2str(yr) ', Cast ' num2str(cast) ', ' datestr(btlsum{4}.Date(indbtl(1)),1)])
    xlabel('Oxygen, \mumol/L','Fontsize',8)
    ylabel('Potential temp','Fontsize',8)

