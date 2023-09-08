%% Read in corrected cruise discrete sample and calibrated cast data
% These data are also in BCO-DMO:
% Fogaren, K. E., Palevsky, H. I. (2023) Bottle-calibrated dissolved oxygen profiles from yearly turn-around cruises for the Ocean Observations Initiative (OOI) Irminger Sea Array 2014 – 2022. Biological and Chemical Oceanography Data Management Office (BCO-DMO). (Version 1) Version Date 2023-07-19. doi:10.26008/1912/bco-dmo.904721.1
% Palevsky, H. I., Fogaren, K. E., Nicholson, D. P., Yoder, M. (2023) Supplementary discrete sample measurements of dissolved oxygen, dissolved inorganic carbon, and total alkalinity from Ocean Observatories Initiative (OOI) cruises to the Irminger Sea Array 2018-2019. Biological and Chemical Oceanography Data Management Office (BCO-DMO). (Version 1) Version Date 2023-07-19. doi:10.26008/1912/bco-dmo.904722.1 

mltoumol = 44.661;

%% Read in processed data and reformat
addpath(genpath('G:\Shared drives\NSF_Irminger\OOI Cruises CTD Casts\CTD_Data\Processed'))
load('AllYears_Processed_KF.mat')

%Reassign variable names for casts
castlist = cast_num;
castsum = downcasts;

%% Create table to hold information about casts and plot locations
load OOImooringLocations.mat %lat-lon coordinates of each year's moorings
clear casts
casts = table;
casts.year = [1*ones(length(castlist{1}),1); 2*ones(length(castlist{2}),1); 3*ones(length(castlist{3}),1);...
    4*ones(length(castlist{4}),1); 5*ones(length(castlist{5}),1); 6*ones(length(castlist{6}),1); 7*ones(length(castlist{7}),1); ...
    8*ones(length(castlist{8}),1); 9*ones(length(castlist{9}),1)];
casts.castnum = [castlist{1}; castlist{2}; castlist{3}; castlist{4}; castlist{5}; castlist{6}; castlist{7}; castlist{8}; castlist{9}];
casts.lat = NaN*ones(length(casts.year),1);
casts.lon = NaN*ones(length(casts.year),1);
casts.time = NaN*ones(length(casts.year),1);
casts.numWinkl = NaN*ones(length(casts.year),1);
casts.HYPMdist = NaN*ones(length(casts.year),1);


%for each year, pull out each cast and plot lat/lon
figure; clf
M = 20;
Ftsz = 12;
colorlist = [nicecolor('rrry'); nicecolor('ryyym'); nicecolor('yyyyyym'); nicecolor('ggyc'); nicecolor('bcc'); nicecolor('c'); nicecolor('cbbbb'); nicecolor('mmmrc'); nicecolor('rrrrm')];

h5 = plot(-OOImoorings.SUMO5(2), OOImoorings.SUMO5(1),'sk','markersize',M,'markerfacecolor',nicecolor('ggk')); hold on;
h6 = plot(-OOImoorings.HYPM5(2), OOImoorings.HYPM5(1),'dk','markersize',M,'markerfacecolor',nicecolor('rb')); hold on;
h7 = plot(-OOImoorings.FLMA5(2), OOImoorings.FLMA5(1),'ok','markersize',M,'markerfacecolor',nicecolor('bw')); hold on;
h8 = plot(-OOImoorings.FLMB5(2), OOImoorings.FLMB5(1),'ok','markersize',M,'markerfacecolor',nicecolor('bbbkkw')); hold on;

for yr = [1:9]
    castsumyr = castsum{yr};
    castlistyr = castlist{yr};
    btlsumyr = btlsum{yr};
    ind = find(casts.year == yr);
    for i = 1:length(castlistyr)
        I = castlistyr(i);
        lat(I) = nanmean(castsumyr{I}.lat);
        lon(I) = nanmean(castsumyr{I}.lon);
        try
            numWinkl(I) = length(btlsumyr{I}.DOcorr_umolkg);
        catch
            numWinkl(I) = 0;
        end
        HYPMdist(I) = distlatlon(OOImoorings.HYPM5(1), lat(I), OOImoorings.HYPM5(2), lon(I)); %need to update with locations for each year
        time(I) = datenum(castsumyr{I}.StartTimeUTC(1));
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
legend([h5 h6 h7 h8 leg{1} leg{2} leg{3} leg{4} leg{5} leg{6} leg{7} leg{8} leg{9}],...
    'Apex surface mooring', 'Apex sub-surface profiler mooring', 'Flanking mooring A', 'Flanking mooring B',...
    'Year 1 casts', 'Year 2 casts', 'Year 3 casts', 'Year 4 casts','Year 5 casts','Year 6 casts','Year 7 casts','Year 8 casts','Year 9 casts',...
    'location','northeastoutside','Fontsize', Ftsz);

%% Pull all data on selected deep isotherm across all casts
%Loop over all casts in castsum to find values interpolated at pt == ISO
for i = 1:height(casts)
    casttbl = castsum{casts.year(i)}(casts.castnum(i));
    casts.ISOdepth(i) = interp1(casttbl{1}.pt, casttbl{1}.depth, ISO);
    casts.ISO_SP(i) = interp1(casttbl{1}.pt, casttbl{1}.SP, ISO);
    casts.ISO_SA(i) = interp1(casttbl{1}.pt, casttbl{1}.SA, ISO);
    casts.ISO_prho(i) = interp1(casttbl{1}.pt, casttbl{1}.prho, ISO);
    casts.ISO_DOcorr_umolkg(i) = interp1(casttbl{1}.pt, casttbl{1}.DOcorr_umolkg, ISO);
end

%% Filter deep isotherm data to remove outliers
%Identify all values within +/- 1 stdev of mean for depth, salinity, and density
ind_ISO_depth = find(casts.ISOdepth < nanmean(casts.ISOdepth) + nanstd(casts.ISOdepth) & casts.ISOdepth > nanmean(casts.ISOdepth) - nanstd(casts.ISOdepth));
ind_ISO_SP = find(casts.ISO_SP < nanmean(casts.ISO_SP) + nanstd(casts.ISO_SP) & casts.ISO_SP > nanmean(casts.ISO_SP) - nanstd(casts.ISO_SP));
ind_ISO_prho = find(casts.ISO_prho < nanmean(casts.ISO_prho) + nanstd(casts.ISO_prho) & casts.ISO_prho > nanmean(casts.ISO_prho) - nanstd(casts.ISO_prho));
ind_phys = intersect(intersect(ind_ISO_depth, ind_ISO_SP),ind_ISO_prho);

%Within the non-outlier values for physical properties, identify DO values within +/- 2 stdev of mean
ind_ISO_DO = find(casts.ISO_DOcorr_umolkg(ind_phys) < nanmean(casts.ISO_DOcorr_umolkg(ind_phys)) + 2*nanstd(casts.ISO_DOcorr_umolkg(ind_phys)) &...
    casts.ISO_DOcorr_umolkg(ind_phys) > nanmean(casts.ISO_DOcorr_umolkg(ind_phys)) - 2*nanstd(casts.ISO_DOcorr_umolkg(ind_phys)));

%Index with outliers removed
ind_ISO = ind_phys(ind_ISO_DO);

%Remove data from 2020 (no Winklers)
ind_no2020 = find(casts.time < datenum(2020,1,1) | casts.time > datenum(2021,1,1));
ind_ISOno2020 = intersect(ind_ISO, ind_no2020);

figure; clf
    subplot(4,1,1)
histogram(casts.ISOdepth(ind_ISOno2020))
    subplot(4,1,2)
histogram(casts.ISO_SP(ind_ISOno2020))
    subplot(4,1,3)
histogram(casts.ISO_prho(ind_ISOno2020))
    subplot(4,1,4)
histogram(casts.ISO_DOcorr_umolkg(ind_ISOno2020))

    M = 16;
figure; clf
    subplot(211)
scatter(casts.time(ind_ISOno2020), casts.ISO_DOcorr_umolkg(ind_ISOno2020),[],casts.HYPMdist(ind_ISOno2020),'filled');
%plot(casts.time(ind_ISOno2020), casts.ISO_DOcorr_umolkg(ind_ISOno2020),'k.','markersize',M);
text(datenum(2014,3,1), 278.5, ['mean = ' num2str(mean(casts.ISO_DOcorr_umolkg(ind_ISOno2020)),4) char(177) num2str(std(casts.ISO_DOcorr_umolkg(ind_ISOno2020)),2) ' \mumol/kg'])
datetick('x','keeplimits')
ylabel('\mumol/kg')
ylim([268 281])
title('Dissolved oxygen on 3.1 deg isotherm from calibrated oxygen profiles on cruise casts')
c = colorbar; caxis([0 50])
ylabel(c, 'Distance from HYPM')
    subplot(212)
scatter(casts.time(ind_ISOno2020), casts.ISO_SP(ind_ISOno2020),[],casts.HYPMdist(ind_ISOno2020),'filled');
text(datenum(2014,3,1), 34.94, ['mean = ' num2str(mean(casts.ISO_SP(ind_ISOno2020)),4) char(177) num2str(std(casts.ISO_SP(ind_ISOno2020)),2)])
datetick('x','keeplimits')
ylabel('SP')
ylim([34.91 34.95])
title('Practical salinity on 3.1 deg isotherm from profiles on cruise casts')
c = colorbar; caxis([0 50])
ylabel(c, 'Distance from HYPM')


%% Uncertainty quantification: how well do we know the deep isotherm?
iso_percenterr = std(casts.ISO_DOcorr_umolkg(ind_ISOno2020))./mean(casts.ISO_DOcorr_umolkg(ind_ISOno2020))*100;
iso_rmse = sqrt(mean((casts.ISO_DOcorr_umolkg(ind_ISOno2020) - mean(casts.ISO_DOcorr_umolkg(ind_ISOno2020))).^2));