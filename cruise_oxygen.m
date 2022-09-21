%% Read in corrected cruise oxygen data
% CTD SBE43 data are corrected by Kristen using Winklers

mltoumol = 44.661;

%% Read in processed data and reformat
addpath(genpath('G:\Shared drives\NSF_Irminger\OOI Cruises CTD Casts\CTD_Data\Alfresco'))
clear castlist

load Year1_Processed_KF.mat
castsum_yr1{2} = cast02;
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
btlsum{1} = btlsum_yr1;
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
for yr = 1:5
    castsumyr = castsum{yr};
    castlistyr = castlist{yr};
    
    %% Plot all casts and interpolate at Winkler depth
    figure(100*yr); clf
    for i = castlistyr
        if i < 17 %set for subplot max
        subplot(4,4,i)
            btlid = find(btlsum{yr}.Cast == i);
        plot(castsumyr{i}.sbeox0mL_L*mltoumol, castsumyr{i}.depSM, 'k-'); hold on;
        if length(btlid) > 0
            plot(btlsum{yr}.Winkler_mLL(btlid)*mltoumol, btlsum{yr}.CTD_Depth_m(btlid), 'ro'); hold on;
            btlsum{yr}.CTD_Oxygen_mLL_corr(btlid) = interp1(castsumyr{i}.depSM, castsumyr{i}.sbeox0mL_L, btlsum{yr}.CTD_Depth_m(btlid));
        end
        axis ij
        xlabel('SBE43 calibrated O_2, \mumol/L')
        ylabel('depSM - db? m?')
        end
    end

    %% Check goodness of calibration

    %Regression plot
    figure(100*yr + 1); clf
    scatter(btlsum{yr}.Winkler_mLL*mltoumol, btlsum{yr}.CTD_Oxygen_mLL_corr*mltoumol, [], btlsum{yr}.CTD_Depth_m, 'filled'); hold on;
    plot([140:350],[140:350],'k--')
    c = colorbar; ylabel(c, 'Depth (m)')
    axis([140 345 140 345])
    xlabel('Winkler oxygen, \mumol/L'); ylabel('Calibrated SBE43 oxygen, \mumol/L')

    %Calculate uncertainty from residuals
    U = nanmean(abs(btlsum{yr}.Winkler_mLL - btlsum{yr}.CTD_Oxygen_mLL_corr))*mltoumol;
        deep_cutoff = 200; %subset deep values only
        deep = find(btlsum{yr}.CTD_Depth_m > deep_cutoff);
    U_deep = nanmean(abs(btlsum{yr}.Winkler_mLL(deep) - btlsum{yr}.CTD_Oxygen_mLL_corr(deep)))*mltoumol;  

    %Plot residuals
    figure(100*yr + 2); clf
    plot(btlsum{yr}.CTD_Depth_m, (btlsum{yr}.Winkler_mLL - btlsum{yr}.CTD_Oxygen_mLL_corr)*mltoumol, 'k.'); hold on;
    plot([0 3000], 0*[0 3000], 'k--')
    ylabel('Residual, Winkler - SBE O_2, \mumol/L')
    xlabel('Depth (m)')
    text(500, 4, ['Mean residual = ' num2str(U) ' \mumol/L'])
    text(500, 3.4, ['Mean residual > ' num2str(deep_cutoff) 'm = ' num2str(U_deep) ' \mumol/L'])

end

%% Create table to hold information about casts and plot locations

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

for yr = 1:5
    castsumyr = castsum{yr};
    castlistyr = castlist{yr};
    ind = find(casts.year == yr);
    for i = castlistyr
        lat(i) = nanmean(castsumyr{i}.lat);
        lon(i) = nanmean(castsumyr{i}.lon);
        numWinkl(i) = length(find(btlsum{yr}.Cast == i & isnan(btlsum{yr}.Winkler_mLL) == 0));
        HYPMdist(i) = distlatlon(HYPMlat(yr), lat(i), HYPMlon(yr), lon(i));
        if numWinkl(i) > 0
            Cind = find(btlsum{yr}.Cast == i);
            time(i) = min(datenum(btlsum{yr}.Date(Cind)));
        else
            time(i) = NaN;
        end
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
legend([h5 h6 h7 h8 leg{1} leg{2} leg{3} leg{4} leg{5}],...
    'Apex surface mooring', 'Apex sub-surface profiler mooring', 'Flanking mooring A', 'Flanking mooring B',...
    'Year 1 casts', 'Year 2 casts', 'Year 3 casts', 'Year 4 casts', 'Year 5 casts',...
    'location','northeastoutside','Fontsize', Ftsz);

%% Extract cruise SBE43 and Winkler data aligned with HYPM

tol_dist = 10; %Extract all casts within X km of HYPM
tol_time = 1.5; %Extract all HYPM profiles within X days of cast
ind = find(casts.HYPMdist < tol_dist & casts.numWinkl > 0);

figure(1); clf
for i = 1:length(ind)
    subplot(4,3,i)
    %Plot the SBE43 profile and Winkler data
    castyr = castsum{casts.year(ind(i))};
    plot(castyr{casts.castnum(ind(i))}.sbeox0mL_L*mltoumol, castyr{casts.castnum(ind(i))}.depSM, 'k-'); hold on;
    if casts.numWinkl(ind(i)) > 0
        btlsumyr = btlsum{casts.year(ind(i))};
        indbtl = find(btlsumyr.Cast == casts.castnum(ind(i)));
        plot(btlsumyr.Winkler_mLL(indbtl)*mltoumol, btlsumyr.DepSM(indbtl), 'bo'); hold on;
    end
    %Extract the corresponding HYPM data
    indHYPM = find(abs(wggmerge.time - casts.time(ind(i))) < tol_time);
    for j = 1:length(indHYPM)
        deploy_yr = wggmerge.deploy_yr(indHYPM(j));
        plot(wggmerge.doxy_lagcorr(:,indHYPM(j)).*wggmerge.pdens(:,indHYPM(j))/1000,...
            depth_grid, '-', 'color', colorlist(deploy_yr,:)); hold on; %divide by potential density to convert to umol/L - will need to recalc b/c salinity issue
    end
    %Plot cleanup
    axis ij
    axis([240 320 0 2700])
    title(['Year ' num2str(casts.year(ind(i))) ', Cast ' num2str(casts.castnum(ind(i))) ', ' datestr(casts.time(ind(i)),1)])
    xlabel('Oxygen, \mumol/L','Fontsize',8)
    ylabel('DepSM','Fontsize',8)
end
