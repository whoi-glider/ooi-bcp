%% Read in corrected cruise oxygen data
% CTD SBE43 data are corrected by Kristen using Winklers

mltoumol = 44.661;

%% Read in processed data and reformat
addpath(genpath('G:\Shared drives\NSF_Irminger\OOI Cruises CTD Casts\CTD_Data\Alfresco'))
clear castlist

load Year1_Processed_KF.mat
castsum_yr1{2} = cast02; %cast02.units - prDM is pressure in db
castsum_yr1{3} = cast03;
castsum_yr1{4} = cast04;
castsum_yr1{5} = cast05;
castsum_yr1{6} = cast06;
castsum_yr1{7} = cast07;
castsum_yr1{8} = cast08;
castsum_yr1{9} = cast09;
castlist{1} = [2:9];

load Year2_Processed_KF.mat
castsum_yr2{1} = cast01;
castsum_yr2{2} = cast02;
castsum_yr2{3} = cast03;
castsum_yr2{4} = cast04;
castsum_yr2{5} = cast05;
castsum_yr2{6} = cast06;
castsum_yr2{7} = cast07;
castsum_yr2{8} = cast08;
castsum_yr2{9} = cast09;
castsum_yr2{10} = cast10;
castsum_yr2{11} = cast11;
castsum_yr2{12} = cast12;
castsum_yr2{13} = cast13;
castlist{2} = [1:13];

load Year3_Processed_KF.mat
castsum_yr3{1} = cast01;
castsum_yr3{2} = cast02;
castsum_yr3{3} = cast03;
castsum_yr3{4} = cast04;
castsum_yr3{5} = cast05;
castsum_yr3{6} = cast06;
castsum_yr3{7} = cast07;
castsum_yr3{8} = cast08;
castsum_yr3{9} = cast09;
castsum_yr3{10} = cast10;
castlist{3} = [1:10];

load Year4_Processed_KF.mat
castsum_yr4{9} = cast09;
castsum_yr4{10} = cast10;
castsum_yr4{11} = cast11;
castsum_yr4{12} = cast12;
castlist{4} = [9:12];

load Year5_Processed_KF.mat
castsum_yr5{1} = cast01;
castsum_yr5{2} = cast02;
castsum_yr5{3} = cast03;
castsum_yr5{4} = cast04;
castsum_yr5{5} = cast05;
castsum_yr5{6} = cast06;
castsum_yr5{7} = cast07;
castsum_yr5{8} = cast08;
castsum_yr5{9} = cast09;
castsum_yr5{10} = cast10;
castsum_yr5{11} = cast11;
castsum_yr5{12} = cast12;
castsum_yr5{13} = cast13;
castsum_yr5{14} = cast14;
castsum_yr5{15} = cast15;
castsum_yr5{16} = cast16;
castsum_yr5{17} = cast17;
castsum_yr5{18} = cast18;
castsum_yr5{19} = cast19;
castsum_yr5{20} = cast20;
castsum_yr5{21} = cast21;
castsum_yr5{22} = cast22;
castsum_yr5{23} = cast23;
castlist{5} = [1:23];

%Merge bottle summary tables
btlsum{1} = btlsum_yr1; %use prDM
btlsum{2} = btlsum_yr2;
btlsum{3} = btlsum_yr3;
btlsum{4} = btlsum_yr4;
btlsum{5} = btlsum_yr5;

%Merge cast structures
castsum{1} = castsum_yr1;
castsum{2} = castsum_yr2;
castsum{3} = castsum_yr3;
castsum{4} = castsum_yr4;
castsum{5} = castsum_yr5;

%% Loop over data from all years
for yr = 1:4
    castsumyr = castsum{yr};
    castlistyr = castlist{yr};
    
    %% Plot all casts and interpolate at Winkler depth
    figure(100*yr); clf
    for i = castlistyr
        if i < 17 %set for subplot max
        subplot(4,4,i)
            btlid = find(btlsum{yr}.Cast == i);
        plot(castsumyr{i}.DOcorr_umolkg, castsumyr{i}.pm, 'k-'); hold on;
        if length(btlid) > 0
            plot(btlsum{yr}.Winkler_umolkg(btlid), btlsum{yr}.PrDM(btlid), 'ro'); hold on;
            %btlsum{yr}.CTD_Oxygen_mLL_corr(btlid) = interp1(castsumyr{i}.pm, castsumyr{i}.sbeox0mL_L, btlsum{yr}.PrDM(btlid));  removed line because Kristen now calculates interpolated DOcorr_umolkg in bottle files
        end
        axis ij
        xlabel('SBE43 calibrated O_2, \mumol/kg')
        ylabel('Pressure (db)')
        end
    end

    %% Check goodness of calibration

    %Regression plot
    figure(100*yr + 1); clf
    scatter(btlsum{yr}.Winkler_umolkg, btlsum{yr}.DOcorr_umolkg, [], btlsum{yr}.PrDM, 'filled'); hold on;
    plot([140:350],[140:350],'k--')
    c = colorbar; ylabel(c, 'Depth (m)')
    axis([140 345 140 345])
    xlabel('Winkler oxygen, \mumol/L'); ylabel('Calibrated SBE43 oxygen, \mumol/kg')

    %Calculate uncertainty from residuals
    U = nanmean(abs(btlsum{yr}.Winkler_umolkg - btlsum{yr}.DOcorr_umolkg));
        deep_cutoff = 200; %subset deep values only
        deep = find(btlsum{yr}.PrDM > deep_cutoff);
    U_deep = nanmean(abs(btlsum{yr}.Winkler_umolkg(deep) - btlsum{yr}.DOcorr_umolkg(deep)));  

    %Plot residuals
    figure(100*yr + 2); clf
    plot(btlsum{yr}.PrDM, (btlsum{yr}.Winkler_umolkg - btlsum{yr}.DOcorr_umolkg), 'k.'); hold on;
    plot([0 3000], 0*[0 3000], 'k--')
    ylabel('Residual, Winkler - SBE O_2, \mumol/kg')
    xlabel('Pressure (db)')
    text(500, 4, ['Mean residual = ' num2str(U) ' \mumol/kg'])
    text(500, 3.4, ['Mean residual > ' num2str(deep_cutoff) 'm = ' num2str(U_deep) ' \mumol/kg'])

end

%% Create table to hold information about casts and plot locations

load OOImooringLocations.mat %lat-lon coordinates of each year's moorings
clear casts
casts = table;
casts.year = [1*ones(length(castlist{1}),1); 2*ones(length(castlist{2}),1); 3*ones(length(castlist{3}),1);...
    4*ones(length(castlist{4}),1); 5*ones(length(castlist{5}),1)];
casts.castnum = [castlist{1}'; castlist{2}'; castlist{3}'; castlist{4}'; castlist{5}';];
casts.lat = NaN*ones(length(casts.year),1);
casts.lon = NaN*ones(length(casts.year),1);
casts.time = NaN*ones(length(casts.year),1);
casts.numWinkl = NaN*ones(length(casts.year),1);
casts.HYPMdist = NaN*ones(length(casts.year),1);


%for each year, pull out each cast and plot lat/lon
figure(10); clf
M = 20;
Ftsz = 12;
colorlist = [nicecolor('rrry'); nicecolor('ryyym'); nicecolor('yyyyyym'); nicecolor('ggyc'); nicecolor('m')];

h5 = plot(-OOImoorings.SUMO5(2), OOImoorings.SUMO5(1),'sk','markersize',M,'markerfacecolor',nicecolor('ggk')); hold on;
h6 = plot(-OOImoorings.HYPM5(2), OOImoorings.HYPM5(1),'dk','markersize',M,'markerfacecolor',nicecolor('rb')); hold on;
h7 = plot(-OOImoorings.FLMA5(2), OOImoorings.FLMA5(1),'ok','markersize',M,'markerfacecolor',nicecolor('bw')); hold on;
h8 = plot(-OOImoorings.FLMB5(2), OOImoorings.FLMB5(1),'ok','markersize',M,'markerfacecolor',nicecolor('bbbkkw')); hold on;

for yr = 1:4
    castsumyr = castsum{yr};
    castlistyr = castlist{yr};
    ind = find(casts.year == yr);
    for i = castlistyr
        lat(i) = nanmean(castsumyr{i}.lat);
        lon(i) = nanmean(castsumyr{i}.lon);
        numWinkl(i) = length(find(btlsum{yr}.Cast == i & isnan(btlsum{yr}.Winkler_umolkg) == 0));
        HYPMdist(i) = distlatlon(HYPMlat(yr), lat(i), HYPMlon(yr), lon(i));
        time(i) = datenum(castsumyr{i}.instrumentheaders.SystemUTC);
    end
    casts.lat(ind) = lat(castlistyr); 
    casts.lon(ind) = lon(castlistyr);
    casts.time(ind) = time(castlistyr);
    casts.numWinkl(ind) = numWinkl(castlistyr);
    casts.HYPMdist(ind) = HYPMdist(castlistyr);
    leg{yr} = plot(-casts.lon(ind), casts.lat(ind), 'ko','markersize',M/2-yr,'markerfacecolor',colorlist(yr,:)); hold on;
end

xlim([39 40.1])
ylim([59.6 60.05])
xlabel('Longitude (^oW)','Fontsize', Ftsz)
ylabel('Latitude (^oN)','Fontsize', Ftsz)
box on
legend([h5 h6 h7 h8 leg{1} leg{2} leg{3} leg{4}],...
    'Apex surface mooring', 'Apex sub-surface profiler mooring', 'Flanking mooring A', 'Flanking mooring B',...
    'Year 1 casts', 'Year 2 casts', 'Year 3 casts', 'Year 4 casts',...
    'location','northeastoutside','Fontsize', Ftsz);

%% Extract cruise SBE43 and Winkler data aligned with HYPM

tol_dist = 10; %Extract all casts within X km of HYPM
tol_time = 2; %Extract all HYPM profiles within X days of cast
ind = find(casts.HYPMdist < tol_dist & casts.numWinkl > 0);

figure(1); clf
for i = 1:length(ind)
    subplot(4,3,i)
    castyr = castsum{casts.year(ind(i))};
    %Interpolate the SBE43 profile onto the HYPM depth grid
    SBE_interp = interp1(castyr{casts.castnum(ind(i))}.pm, castyr{casts.castnum(ind(i))}.DOcorr_umolkg.*castyr{casts.castnum(ind(i))}.prho/1000, pres_grid); %umol/L
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
    plot(castyr{casts.castnum(ind(i))}.DOcorr_umolkg.*castyr{casts.castnum(ind(i))}.prho/1000, castyr{casts.castnum(ind(i))}.pm, 'k-','linewidth',1.5); hold on;
    if casts.numWinkl(ind(i)) > 0
        btlsumyr = btlsum{casts.year(ind(i))};
        indbtl = find(btlsumyr.Cast == casts.castnum(ind(i)));
        plot(btlsumyr.Winkler_mLL(indbtl)*mltoumol, btlsumyr.PrDM(indbtl), 'bo'); hold on;
    end
    %Plot cleanup
    axis ij
    axis([240 320 0 2700])
    title(['Year ' num2str(casts.year(ind(i))) ', Cast ' num2str(casts.castnum(ind(i))) ', ' datestr(casts.time(ind(i)),1)])
    xlabel('Oxygen, \mumol/L','Fontsize',8)
    ylabel('Pressure (db)','Fontsize',8)
end

%% Approach of applying best initial gain correction for each year - plotted via pressure
tol_time = 1.5;
minnumpts = 500; %number of points in HYPM profile if usable for comparison

figure(100); clf
for yr = [1,2,3];
    
subplot(3,2,yr)

if yr == 1
    cast = 6;
elseif yr == 2
    cast = 12;
elseif yr == 3
    cast = 7;
elseif yr == 5
    cast = 6;
end

castyr = castsum{yr};
ind = find(casts.year == yr & casts.castnum == cast);

    %Interpolate the SBE43 profile onto the HYPM depth grid
    SBE_interp = interp1(castyr{casts.castnum(ind)}.pm, castyr{casts.castnum(ind)}.DOcorr_umolkg.*castyr{casts.castnum(ind)}.prho/1000, pres_grid);
    %Extract the corresponding HYPM data
    indHYPM = find(abs(wggmerge.time - casts.time(ind)) < tol_time & wggmerge.deploy_yr == yr);
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
    plot(castyr{casts.castnum(ind)}.DOcorr_umolkg.*castyr{casts.castnum(ind)}.prho/1000, castyr{casts.castnum(ind)}.pm, 'k-','linewidth',1.5); hold on;
    if casts.numWinkl(ind) > 0
        btlsumyr = btlsum{casts.year(ind)};
        indbtl = find(btlsumyr.Cast == casts.castnum(ind));
        plot(btlsumyr.Winkler_mLL(indbtl)*mltoumol, btlsumyr.PrDM(indbtl), 'bo','markerfacecolor','b'); hold on;
    end
    %Plot cleanup
    axis ij
    axis([250 320 0 2700])
    title(['Year ' num2str(yr) ', Cast ' num2str(cast) ', ' datestr(casts.time(ind),1)])
    xlabel('Oxygen, \mumol/L','Fontsize',8)
    ylabel('Pressure (db)','Fontsize',8)

end

% %Plot year 4 as a special case with no SBE data
% yr = 4; cast = 8;
% subplot(3,2,yr)
% 
% %Extract bottle IDs of Winkler samples
% indbtl = find(btlsum{4}.Cast == cast);
% 
% %Extract the corresponding HYPM data
%     indHYPM = find(abs(wggmerge.time - datenum(btlsum{4}.Date(indbtl(1))) < tol_time) & wggmerge.deploy_yr == yr);
%     for j = 1:length(indHYPM)
%         if sum(~isnan(wggmerge.doxy_lagcorr(:,indHYPM(j)))) > minnumpts
%             HYPMprofile_O2umolL = wggmerge.doxy_lagcorr(:,indHYPM(j)).*wggmerge.pdens(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
%             plot(HYPMprofile_O2umolL, pres_grid, '-', 'color', colorlist(yr,:)*0.6); hold on; 
%             %Calculate gain correction
%             HYPM_interp = interp1(pres_grid, HYPMprofile_O2umolL, btlsum{4}.PrDM(indbtl));
%             gain = btlsum{4}.Winkler_mLL(indbtl)*mltoumol./HYPM_interp;
%             gainm(j) = nanmean(gain);
%             gains(j) = nanstd(gain);
%         else
%             gainm(j) = NaN; gains(j) = NaN;
%         end
%     end
%     gain_yr(yr,1) = nanmean(gainm);
%     gain_yr(yr,2) = nanstd(gainm);
%     clear gainm gains
%     
%     %Replot with the corrected data
%     for j = 1:length(indHYPM)
%         HYPMprofile_O2umolL = wggmerge.doxy_lagcorr(:,indHYPM(j)).*wggmerge.pdens(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
%         plot(HYPMprofile_O2umolL*gain_yr(yr,1), pres_grid, '-', 'color', colorlist(yr,:)); hold on; 
%     end
%     
%     %Plot the Winkler data
%     plot(btlsum{4}.Winkler_mLL(indbtl)*mltoumol, btlsum{4}.PrDM(indbtl), 'bo','markerfacecolor','b'); hold on;
% 
%     %Plot cleanup
%     axis ij
%     axis([250 320 0 2700])
%     title(['Year ' num2str(yr) ', Cast ' num2str(cast) ', ' datestr(btlsum{4}.Date(indbtl(1)),1)])
%     xlabel('Oxygen, \mumol/L','Fontsize',8)
%     ylabel('Pressure (db)','Fontsize',8)
    
%% Approach of applying best initial gain correction for each year - plotted via temperature
tol_time = 1.5;
minnumpts = 500; %number of points in HYPM profile if usable for comparison

figure(101); clf
for yr = [1,2,3];
    
subplot(3,2,yr)

if yr == 1
    cast = 6;
elseif yr == 2
    cast = 12;
elseif yr == 3
    cast = 7;
elseif yr == 5
    cast = 6;
end

castyr = castsum{yr};
ind = find(casts.year == yr & casts.castnum == cast);

%Pull out relevant cast information to save for later use
castalign{yr}.castnum = cast;
idc = find(casts.year == yr & casts.castnum == cast);
castalign{yr}.time = casts.time(idc);

    %Interpolate the SBE43 profile onto the HYPM depth grid
        %convert from measured to potential temperature - note this is
        %calculating salinity straight from conductivity, no check with
        %bottle salts
    CTD_SP = gsw_SP_from_C(castyr{casts.castnum(ind)}.c0mScm, castyr{casts.castnum(ind)}.t090C, castyr{casts.castnum(ind)}.pm);
    [castalign{yr}.CTD_SA,~] = gsw_SA_from_SP(CTD_SP, castyr{casts.castnum(ind)}.pm, castyr{casts.castnum(ind)}.lon, castyr{casts.castnum(ind)}.lat);
    castalign{yr}.CTD_pt = gsw_pt0_from_t(castalign{yr}.CTD_SA, castyr{casts.castnum(ind)}.t090C, castyr{casts.castnum(ind)}.pm);
    castalign{yr}.SBE_interp = interp1(castalign{yr}.CTD_pt, castyr{casts.castnum(ind)}.DOcorr_umolkg.*castyr{casts.castnum(ind)}.prho/1000, pt_grid);
    %Extract the corresponding HYPM data
    indHYPM = find(abs(wggmerge.time - casts.time(ind)) < tol_time & wggmerge.deploy_yr == yr);
    for j = 1:length(indHYPM)
        if sum(~isnan(wggmerge.doxy_lagcorr(:,indHYPM(j)))) > minnumpts
            HYPMprofile_O2umolL = wggmerge.doxy_lagcorr_pt(:,indHYPM(j)).*wggmerge.pdens_pt(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
            plot(HYPMprofile_O2umolL, pt_grid, '-', 'color', colorlist(yr,:)*0.6); hold on; 
            %Calculate gain correction
            gain = castalign{yr}.SBE_interp'./HYPMprofile_O2umolL;
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
    %Plot the SBE43 profile and Winkler data
    plot(castyr{casts.castnum(ind)}.DOcorr_umolkg.*castyr{casts.castnum(ind)}.prho/1000, castalign{yr}.CTD_pt, 'k-','linewidth',1.5); hold on;
    if casts.numWinkl(ind) > 0
        btlsumyr = btlsum{casts.year(ind)};
        indbtl = find(btlsumyr.Cast == casts.castnum(ind));
        bottle_SP = gsw_SP_from_C(btlsumyr.C0mS_cm(indbtl), btlsumyr.T090C(indbtl), btlsumyr.PrDM(indbtl));
        [bottle_SA,~] = gsw_SA_from_SP(bottle_SP, btlsumyr.PrDM(indbtl), btlsumyr.Longitude(indbtl), btlsumyr.Latitude(indbtl));
        bottle_pt = gsw_pt0_from_t(bottle_SA, btlsumyr.T090C(indbtl), btlsumyr.PrDM(indbtl));
        plot(btlsumyr.Winkler_mLL(indbtl)*mltoumol, bottle_pt, 'bo','markerfacecolor','b'); hold on;
    end
    %Plot cleanup
    axis([250 320 1.5 5])
    title(['Year ' num2str(yr) ', Cast ' num2str(cast) ', ' datestr(casts.time(ind),1)])
    xlabel('Oxygen, \mumol/L','Fontsize',8)
    ylabel('Potential temp','Fontsize',8)

end

% %Plot year 4 as a special case with no SBE data
% yr = 4; cast = 8;
% subplot(3,2,yr)
% 
% %Extract bottle IDs of Winkler samples
% indbtl = find(btlsum{4}.Cast == cast);
% 
% %Calculate bottle potential temp
% bottle_pt = gsw_pt0_from_t(btlsum{4}.C0mS_cm(indbtl), btlsum{4}.T090C(indbtl), btlsum{4}.PrDM(indbtl));
% 
% %Extract the corresponding HYPM data
%     indHYPM = find(abs(wggmerge.time - datenum(btlsum{4}.Date(indbtl(1))) < tol_time) & wggmerge.deploy_yr == yr);
%     for j = 1:length(indHYPM)
%         if sum(~isnan(wggmerge.doxy_lagcorr(:,indHYPM(j)))) > minnumpts
%             HYPMprofile_O2umolL = wggmerge.doxy_lagcorr_pt(:,indHYPM(j)).*wggmerge.pdens_pt(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
%             plot(HYPMprofile_O2umolL, pt_grid, '-', 'color', colorlist(yr,:)*0.6); hold on; 
%             %Calculate gain correction
%             HYPM_interp = interp1(pt_grid, HYPMprofile_O2umolL, bottle_pt);
%             gain = btlsum{4}.Winkler_mLL(indbtl)*mltoumol./HYPM_interp;
%             gainm(j) = nanmean(gain);
%             gains(j) = nanstd(gain);
%         else
%             gainm(j) = NaN; gains(j) = NaN;
%         end
%     end
%     gain_yr_pt(yr,1) = nanmean(gainm);
%     gain_yr_pt(yr,2) = nanstd(gainm);
%     clear gainm gains
%     
%     %Replot with the corrected data
%     for j = 1:length(indHYPM)
%         HYPMprofile_O2umolL = wggmerge.doxy_lagcorr_pt(:,indHYPM(j)).*wggmerge.pdens_pt(:,indHYPM(j))/1000; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
%         plot(HYPMprofile_O2umolL*gain_yr_pt(yr,1), pt_grid, '-', 'color', colorlist(yr,:)); hold on; 
%     end
%     
%     %Plot the Winkler data
%     plot(btlsum{4}.Winkler_mLL(indbtl)*mltoumol, bottle_pt, 'bo','markerfacecolor','b'); hold on;
% 
%     %Plot cleanup
%     axis([250 320 1.5 5])
%     title(['Year ' num2str(yr) ', Cast ' num2str(cast) ', ' datestr(btlsum{4}.Date(indbtl(1)),1)])
%     xlabel('Oxygen, \mumol/L','Fontsize',8)
%     ylabel('Potential temp','Fontsize',8)

