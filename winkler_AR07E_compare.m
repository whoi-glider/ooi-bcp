
%% Pull out AR07E data using pre-existing functions (even though not using MLR part here)
[data_05, ~, ~, ~] = MLR_fromGOSHIPdata('C:/Users/palevsky/Dropbox/WHOI Postdoc/OUC seed proposal/AR07E-2014/64PE20050907.csv', 2); %note that 2nd input variable is fudge factor to deal with different arrangement of columns
[data_07, ~, ~, ~] = MLR_fromGOSHIPdata('C:/Users/palevsky/Dropbox/WHOI Postdoc/OUC seed proposal/AR07E-2014/64PE20070830.csv', 0);
[data_14, ~, ~, ~] = MLR_fromGOSHIPdata_netcdf('C:/Users/palevsky/Dropbox/WHOI Postdoc/OUC seed proposal/AR07E-2014/sam_jr302_all.nc');

close all %because don't need the graphs this creates

%% Create location and depth bounds to extract AR07E data from for comparison

latmin = 58;
latmax = 62;
lonmin = -43;
lonmax = -35;
presmin = 2000;

    ind_a = find(data_05.lon > lonmin & data_05.lon < lonmax);
    ind_b = find(data_05.lat > latmin & data_05.lat < latmax);
    ind_c = find(data_05.pres > presmin);
ind_05 = intersect(intersect(ind_a,ind_b), ind_c);

    ind_a = find(data_07.lon > lonmin & data_07.lon < lonmax);
    ind_b = find(data_07.lat > latmin & data_07.lat < latmax);
    ind_c = find(data_07.pres > presmin);
ind_07 = intersect(intersect(ind_a,ind_b), ind_c);

    ind_a = find(data_14.lon > lonmin & data_14.lon < lonmax);
    ind_b = find(data_14.lat > latmin & data_14.lat < latmax);
    ind_c = find(data_14.pres > presmin);
ind_14 = intersect(intersect(ind_a,ind_b), ind_c);

%% Read in Irminger-5 Winkler data
[Winkler5_casts] = xlsread('C:/Users/palevsky/Dropbox/Irminger5/Oxygen data/Irminger5_WinklerSamples.xlsx',2);
Winkler5.depth = Winkler5_casts(:,8);
Winkler5.T = Winkler5_casts(:,9);
Winkler5.S = Winkler5_casts(:,10);
Winkler5.O2_dave = Winkler5_casts(:,17);
Winkler5.O2_bcp = mean(Winkler5_casts(:,21:22),2);
Winkler5.O2_bcp_flag = Winkler5_casts(:,23);
%Calculate O2_equil and AOU
Winkler5.O2equil = gsw_O2sol_SP_pt(Winkler5.S, Winkler5.T);
Winkler5.AOU_dave = Winkler5.O2equil - Winkler5.O2_dave;
Winkler5.AOU_bcp = Winkler5.O2equil - Winkler5.O2_bcp;

%Indices for deep data to compare
ind_Winkler5 = find(Winkler5.depth > 2000);
ind_Winkler5_bcp = find(Winkler5.depth > 2000 & Winkler5.O2_bcp_flag == 1);

%Index for good values from both to compare
ind_compare5 = find(isnan(Winkler5.O2_dave) == 0 & Winkler5.O2_bcp_flag == 1);

%% Read in Irminger-6 Winkler data
[Winkler6_casts] = xlsread('C:/Users/palevsky/Dropbox/Irminger6/Oxygen data/Irminger6_Winkler.xlsx',2);
Winkler6.depth = Winkler6_casts(:,8);
Winkler6.T = Winkler6_casts(:,9);
Winkler6.S = Winkler6_casts(:,11);
Winkler6.O2_dave = Winkler6_casts(:,18);
Winkler6.O2_dave_flag = Winkler6_casts(:,19);
Winkler6.O2_bcp = Winkler6_casts(:,29);
Winkler6.O2_bcp_flag = Winkler6_casts(:,32);

%Calculate O2_equil and AOU
Winkler6.O2equil = gsw_O2sol_SP_pt(Winkler6.S, Winkler6.T);
Winkler6.AOU_dave = Winkler6.O2equil - Winkler6.O2_dave;
Winkler6.AOU_bcp = Winkler6.O2equil - Winkler6.O2_bcp;

%Indices for deep data to compare
ind_Winkler6 = find(Winkler6.depth > 2000);
ind_Winkler6_bcp = find(Winkler6.depth > 2000 & Winkler6.O2_bcp_flag > 0);

%Index for good values from both to compare
ind_compare6 = find(Winkler6.O2_dave_flag < 3 & Winkler6.O2_bcp_flag > 0);

%% Plot comparison between Irminger-5 Winklers from BCP team and Wellwood
figure(10); clf
M = 15;
    subplot(121)
plot(Winkler5.O2_dave(ind_compare5), Winkler5.O2_bcp(ind_compare5), 'k.','markersize',M); hold on;
plot([270:350],[270:350],'r--'); hold on;
xlabel('Wellwood O_2 (\mumol/kg)'); ylabel('BCP team O_2 (\mumol/kg)');
    offset = (Winkler5.O2_bcp(ind_compare5) - Winkler5.O2_dave(ind_compare5));
    ratio = (Winkler5.O2_bcp(ind_compare5)./Winkler5.O2_dave(ind_compare5));
text(275,340,{['offset = ' num2str(mean(offset)) ' +/- ' num2str(std(offset))]})
text(275,332,{['ratio = ' num2str(mean(ratio)) ' +/- ' num2str(std(ratio))]})
title('Irminger5 (June 2018) - Winkler lab intercomparison')
    subplot(122)
plot(Winkler6.O2_dave(ind_compare6), Winkler6.O2_bcp(ind_compare6), 'k.','markersize',M); hold on;
plot([270:310],[270:310],'r--'); hold on;
xlabel('Wellwood O_2 (\mumol/kg)'); ylabel('BCP team O_2 (\mumol/kg)');
    offset = (Winkler6.O2_bcp(ind_compare6) - Winkler6.O2_dave(ind_compare6));
    ratio = (Winkler6.O2_bcp(ind_compare6)./Winkler6.O2_dave(ind_compare6));
text(275,305,{['offset = ' num2str(mean(offset)) ' +/- ' num2str(std(offset))]})
text(275,301,{['ratio = ' num2str(mean(ratio)) ' +/- ' num2str(std(ratio))]})
title('Irminger6 (August 2019) - Winkler lab intercomparison')

%% Read in earlier Irminger Winkler data
%Currently just using years 2-3 (Yrs 1 & 4 lack bottle file for pot temp
%and S)

loadWinklerIrmingerYrs1to5
addpath('C:\Users\palevsky\Dropbox\Wellesley\OOI_Irminger_students\CruiseData_Yrs1to4')

%Calculate potential density
Yr2_disc.pdens = sw_dens0(Yr2_disc.S, Yr2_disc.potT);
Yr3_disc.pdens = sw_dens0(Yr3_disc.S, Yr3_disc.potT);

%Calculate O2_equil
Yr2_disc.O2equil = gsw_O2sol_SP_pt(Yr2_disc.S, Yr2_disc.potT);
Yr3_disc.O2equil = gsw_O2sol_SP_pt(Yr3_disc.S, Yr3_disc.potT);

%Calculate AOU - note that O2 values are in mmol/L and need to be converted to umol/kg
Yr2_disc.AOU = Yr2_disc.O2equil - Yr2_disc.oxy./(Yr2_disc.pdens/1000);
Yr3_disc.AOU = Yr3_disc.O2equil - Yr3_disc.oxy./(Yr3_disc.pdens/1000);

%Indices for deep data to compare
Yr2_disc.ind_deep = find(Yr2_disc.depth > 2000);
Yr3_disc.ind_deep = find(Yr3_disc.depth > 2000);

%% Visualize AOU relationship with T and S
%Note that I have read in the 

figure(1); clf
M1 = 5;
M2 = 18;
    subplot(211)
scatter(data_05.T(ind_05), data_05.S(ind_05), [], data_05.AOU(ind_05),'filled'); hold on;
scatter(data_07.T(ind_07), data_07.S(ind_07), [], data_07.AOU(ind_07),'filled'); hold on;
%scatter(data_14.T(ind_14), data_14.S(ind_14), [], data_14.AOU(ind_14),'filled'); hold on;
%scatter(Yr2_disc.potT(Yr2_disc.ind_deep), Yr2_disc.S(Yr2_disc.ind_deep), [], Yr2_disc.AOU(Yr2_disc.ind_deep),'filled'); hold on;
%scatter(Yr3_disc.potT(Yr3_disc.ind_deep), Yr3_disc.S(Yr3_disc.ind_deep), [], Yr3_disc.AOU(Yr3_disc.ind_deep),'filled'); hold on;
xlabel('T'); ylabel('S');
h = colorbar; ylabel(h,'AOU'); title('AR07E, 2005 and 2007, below 2000 m')
    subplot(212)
plot(sw_dens0(data_05.S(ind_05), data_05.T(ind_05)) - 1000, data_05.AOU(ind_05), 'r.','markersize',M1); hold on;
plot(sw_dens0(data_07.S(ind_07), data_07.T(ind_07)) - 1000, data_07.AOU(ind_07), 'b.','markersize',M1); hold on;
plot(sw_dens0(data_14.S(ind_14), data_14.T(ind_14)) - 1000, data_14.AOU(ind_14), 'k.','markersize',M1); hold on;
plot(sw_dens0(Winkler5.S(ind_Winkler5_bcp), Winkler5.T(ind_Winkler5_bcp)) - 1000, Winkler5.AOU_bcp(ind_Winkler5_bcp), '.','color',nicecolor('mmcw'),'markersize',M2); hold on;
plot(sw_dens0(Winkler5.S(ind_Winkler5), Winkler5.T(ind_Winkler5)) - 1000, Winkler5.AOU_dave(ind_Winkler5), 'c.','markersize',M2); hold on;
plot(sw_dens0(Winkler6.S(ind_Winkler6_bcp), Winkler6.T(ind_Winkler6_bcp)) - 1000, Winkler6.AOU_bcp(ind_Winkler6_bcp), '.','color',nicecolor('ryyw'),'markersize',M2); hold on;
plot(sw_dens0(Winkler6.S(ind_Winkler6), Winkler6.T(ind_Winkler6)) - 1000, Winkler6.AOU_dave(ind_Winkler6), 'r.','markersize',M2); hold on;
%plot(Yr2_disc.pdens(Yr2_disc.ind_deep) - 1000, Yr2_disc.AOU(Yr2_disc.ind_deep), 'g.','markersize',M2); hold on;
%plot(Yr3_disc.pdens(Yr3_disc.ind_deep) - 1000, Yr3_disc.AOU(Yr3_disc.ind_deep), 'm.','markersize',M2); hold on;
legend('2005','2007','2014','Irminger5, BCP','Irminger5, Wellwood','Irminger6, BCP','Irminger6, Wellwood','Irminger2','Irminger3','location','southwest')
xlabel('\sigma_0'); ylabel('AOU'); title('AR07E, below 2000 m')

