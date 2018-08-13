
%% Pull out AR07E data using pre-existing functions (even though not using MLR part here)
[data_05, ~, ~, ~] = MLR_fromGOSHIPdata('C:/Users/Hilary/Dropbox/WHOI Postdoc/OUC seed proposal/AR07E-2014/64PE20050907.csv', 2); %note that 2nd input variable is fudge factor to deal with different arrangement of columns
[data_07, ~, ~, ~] = MLR_fromGOSHIPdata('C:/Users/Hilary/Dropbox/WHOI Postdoc/OUC seed proposal/AR07E-2014/64PE20070830.csv', 0);
[data_14, ~, ~, ~] = MLR_fromGOSHIPdata_netcdf('C:/Users/Hilary/Dropbox/WHOI Postdoc/OUC seed proposal/AR07E-2014/sam_jr302_all.nc');

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
[Winkler_casts] = xlsread('C:/Users/Hilary/Google Drive/2018 Irminger/Data_Winkler/Irminger5_WinklerSamples.xlsx',2);
Winkler.depth = Winkler_casts(:,8);
Winkler.T = Winkler_casts(:,9);
Winkler.S = Winkler_casts(:,10);
Winkler.O2_dave = Winkler_casts(:,17);
Winkler.O2_bcp = mean(Winkler_casts(:,21:22),2);
Winkler.O2_bcp_flag = Winkler_casts(:,23);

%Calculate O2_equil and AOU
Winkler.O2equil = gsw_O2sol_SP_pt(Winkler.S, Winkler.T);
Winkler.AOU_dave = Winkler.O2equil - Winkler.O2_dave;
Winkler.AOU_bcp = Winkler.O2equil - Winkler.O2_bcp;

%Indices for deep data to compare
ind_Winkler = find(Winkler.depth > 2000);
ind_Winkler_bcp = find(Winkler.depth > 2000 & Winkler.O2_bcp_flag == 1);

%Index for good values from both to compare
ind_compare = find(isnan(Winkler.O2_dave) == 0 & Winkler.O2_bcp_flag == 1);

%% Plot comparison between Irminger-5 Winklers from BCP team and Wellwood
figure(10); clf
plot(Winkler.O2_dave(ind_compare), Winkler.O2_bcp(ind_compare), 'k.'); hold on;
plot([270:350],[270:350],'r--'); hold on;
xlabel('Wellwood O_2 (\mumol/kg)'); ylabel('BCP team O_2 (\mumol/kg)');
    offset = (Winkler.O2_bcp(ind_compare) - Winkler.O2_dave(ind_compare));
    ratio = (Winkler.O2_bcp(ind_compare)./Winkler.O2_dave(ind_compare));
text(275,340,{['offset = ' num2str(mean(offset)) ' +/- ' num2str(std(offset))]})
text(275,332,{['ratio = ' num2str(mean(ratio)) ' +/- ' num2str(std(ratio))]})

%% Read in earlier Irminger Winkler data
%Currently just using years 2-3 (Yr1 lacks bottle file for pot temp and S,
%Yr 4 spreadsheet not yet created)

loadWinklerIrminger

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
scatter(Yr2_disc.potT(Yr2_disc.ind_deep), Yr2_disc.S(Yr2_disc.ind_deep), [], Yr2_disc.AOU(Yr2_disc.ind_deep),'filled'); hold on;
scatter(Yr3_disc.potT(Yr3_disc.ind_deep), Yr3_disc.S(Yr3_disc.ind_deep), [], Yr3_disc.AOU(Yr3_disc.ind_deep),'filled'); hold on;
xlabel('T'); ylabel('S');
h = colorbar; ylabel(h,'AOU'); title('AR07E, 2005 and 2007, below 2000 m')
    subplot(212)
plot(sw_dens0(data_05.S(ind_05), data_05.T(ind_05)) - 1000, data_05.AOU(ind_05), 'r.','markersize',M1); hold on;
plot(sw_dens0(data_07.S(ind_07), data_07.T(ind_07)) - 1000, data_07.AOU(ind_07), 'b.','markersize',M1); hold on;
plot(sw_dens0(data_14.S(ind_14), data_14.T(ind_14)) - 1000, data_14.AOU(ind_14), 'k.','markersize',M1); hold on;
plot(sw_dens0(Winkler.S(ind_Winkler_bcp), Winkler.T(ind_Winkler_bcp)) - 1000, Winkler.AOU_bcp(ind_Winkler_bcp), '.','color',nicecolor('mmcw'),'markersize',M2); hold on;
plot(sw_dens0(Winkler.S(ind_Winkler), Winkler.T(ind_Winkler)) - 1000, Winkler.AOU_dave(ind_Winkler), 'c.','markersize',M2); hold on;
plot(Yr2_disc.pdens(Yr2_disc.ind_deep) - 1000, Yr2_disc.AOU(Yr2_disc.ind_deep), 'g.','markersize',M2); hold on;
plot(Yr3_disc.pdens(Yr3_disc.ind_deep) - 1000, Yr3_disc.AOU(Yr3_disc.ind_deep), 'm.','markersize',M2); hold on;
legend('2005','2007','2014','Irminger5, BCP','Irminger5, Wellwood','Irminger2','Irminger3','location','southwest')
xlabel('\sigma_0'); ylabel('AOU'); title('AR07E, below 2000 m')

