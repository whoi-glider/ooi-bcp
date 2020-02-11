%% Read in calibration cast data
[Winkler_casts, Winkler_text] = xlsread('C:\Users\palevsky\Dropbox\Irminger5\Oxygen data\Irminger5_WinklerSamples.xlsx',2);
    Winkler.depth = Winkler_casts(:,8);
    Winkler.T = Winkler_casts(:,9); %potential temp
    Winkler.S = Winkler_casts(:,10); %practical salinity
    Winkler.O2_Dave = Winkler_casts(:,17); %Dave Wellwood values in umol/kg
    Winkler.O2_Dave_flag = Winkler_casts(:,18); %Quality flag (2 for different Niskin from BCP group)
    Winkler.O2_BCP = Winkler_casts(:,21:22); %Replicate samples from BCP group, umol/kg
    Winkler.O2_BCP_flag = Winkler_casts(:,23); %Quality flag
    
    %Calculated properties
    Winkler.O2_equil = gsw_O2sol_SP_pt(Winkler.S,Winkler.T); %solubility of O2 in umol/kg
    Winkler.SA = gsw_SA_from_SP(Winkler.S,0,-40,60); %absolute salinity
    Winkler.dens0 = gsw_sigma0(Winkler.SA, gsw_CT_from_pt(Winkler.SA,Winkler.T)) + 1000;

%Identify the indices for the two glider calibration casts
indGliderCast1 = find(Winkler_casts(:,1) == 10);
indGliderCast2 = find(Winkler_casts(:,1) == 21);

%Identify the indices for BCP Winkler data with good quality flag
indGoodBCPdata = find(Winkler.O2_BCP_flag == 1);

%% Plot profiles aligned with calibration casts and make rough initial gain correction

    M = 12;
figure(6); clf
subplot(221)
endtime = datenum(2018,6,12,23,45,0); starttime = datenum(2018,6,12,19,0,0); %Times corresponding with first calibration cast
    hold on;
d_363 = find(G363.daten < endtime & G363.daten > starttime);
    plot(G363.O2_corr(d_363), G363.depth_interp(d_363), '.', 'color', nicecolor('bbbww'))
d_453 = find(G453.daten < endtime & G453.daten > starttime);
    plot(G453.O2_corr(d_453), G453.depth_interp(d_453), '.', 'color', nicecolor('ry'))
plot(mean(Winkler.O2_BCP(intersect(indGliderCast1, indGoodBCPdata)),2).*(Winkler.dens0(intersect(indGliderCast1, indGoodBCPdata))/1000),...
    Winkler.depth(intersect(indGliderCast1, indGoodBCPdata)), 'k.','markersize',M)
plot(Winkler.O2_Dave(indGliderCast1).*(Winkler.dens0(indGliderCast1)/1000),...
    Winkler.depth(indGliderCast1),'.','color', nicecolor('k'),'markersize',M)
set(gca,'YDir','reverse'); 
axis([260 380 1 1005])
legend('Glider 363 oxygen','Glider 453 oxygen','Winkler oxygen','location','southeast')
title('O_2 concentration data from first calibration cast (\mumol/L)')

subplot(223)
    hold on;
d = find(G363.daten < endtime & G363.daten > starttime);
    plot(G363.O2sat_corr(d), G363.depth_interp(d), '.', 'color', nicecolor('bbbww'))
d = find(G453.daten < endtime & G453.daten > starttime);
    plot(G453.O2sat_corr(d), G453.depth_interp(d), '.','color',nicecolor('ry'))
plot(mean(Winkler.O2_BCP(intersect(indGliderCast1, indGoodBCPdata)),2)./(Winkler.O2_equil(intersect(indGliderCast1, indGoodBCPdata)))*100,...
    Winkler.depth(intersect(indGliderCast1, indGoodBCPdata)), 'k.','markersize',M)
plot(Winkler.O2_Dave(indGliderCast1)./(Winkler.O2_equil(indGliderCast1))*100,...
    Winkler.depth(indGliderCast1),'.','color', nicecolor('k'),'markersize',M)
set(gca,'YDir','reverse'); 
axis([80 120 1 1005])
legend('Glider 363 oxygen','Glider 453 oxygen','Winkler oxygen','location','southeast')
title('O_2 saturation data from first calibration cast')

subplot(222)
starttime = datenum(2018,6,19,16,35,0); endtime = datenum(2018,6,19,21,0,0); %Times corresponding with second glider calibration cast
    hold on;
d = find(G363.daten < endtime & G363.daten > starttime);
    plot(G363.O2_corr(d), G363.depth_interp(d), '.', 'color', nicecolor('bbbww'))
d = find(G453.daten < endtime & G453.daten > starttime);
    plot(G453.O2_corr(d), G453.depth_interp(d), '.', 'color', nicecolor('ry'))
plot(mean(Winkler.O2_BCP(intersect(indGliderCast2, indGoodBCPdata)),2).*(Winkler.dens0(intersect(indGliderCast2, indGoodBCPdata))/1000),...
    Winkler.depth(intersect(indGliderCast2, indGoodBCPdata)), 'k.','markersize',M)
set(gca,'YDir','reverse'); 
axis([260 380 1 1005])
legend('Glider 363 oxygen','Glider 453 oxygen','Winkler oxygen','location','southeast')
title('O_2 concentration data from second calibration cast (\mumol/L)')

subplot(224)
    hold on;
d = find(G363.daten < endtime & G363.daten > starttime);
    plot(G363.O2sat_corr(d), G363.depth_interp(d), '.', 'color', nicecolor('bbbww'))
d = find(G453.daten < endtime & G453.daten > starttime);
    plot(G453.O2sat_corr(d), G453.depth_interp(d), '.','color', nicecolor('ry'))
plot(mean(Winkler.O2_BCP(intersect(indGliderCast2, indGoodBCPdata)),2)./Winkler.O2_equil(intersect(indGliderCast2, indGoodBCPdata))*100,...
    Winkler.depth(intersect(indGliderCast2, indGoodBCPdata)), 'k.','markersize',M)
set(gca,'YDir','reverse'); 
axis([80 120 1 1005])
legend('Glider 363 oxygen','Glider 453 oxygen','Winkler oxygen','location','southeast')
title('O_2 saturation data from second calibration cast')

%% Rough calculation of gain correction from Winkler data
    tol = 2; %pull glider data from 1 m above and below

    %%% Calculate for first calibration cast, just with Wellwood Winkler data
ind_good = intersect(find(Winkler.O2_Dave_flag > 0), indGliderCast1);
    gain = NaN*ones(length(ind_good),2);
    gain_sat = NaN*ones(length(ind_good),2);
for i = 2:length(ind_good)
    ind_363 = find(G363.depth_interp(d_363) < Winkler.depth(ind_good(i)) + tol & G363.depth_interp(d_363) > Winkler.depth(ind_good(i)) - tol);
    ind_453 = find(G453.depth_interp(d_453) < Winkler.depth(ind_good(i)) + tol & G453.depth_interp(d_453) > Winkler.depth(ind_good(i)) - tol);
    if length(ind_453) > 0
        gain(i,1) = (Winkler.O2_Dave(ind_good(i)).*(Winkler.dens0(ind_good(i))/1000))./nanmean(G453.O2_corr(d_453(ind_453)));
        gain_sat(i,1) = (Winkler.O2_Dave(ind_good(i))./Winkler.O2_equil(ind_good(i))*100)./nanmean(G453.O2sat_corr(d_453(ind_453)));
    end
    if length(ind_363 > 0)
        gain(i,2) = (Winkler.O2_Dave(ind_good(i)).*(Winkler.dens0(ind_good(i))/1000))./nanmean(G363.O2_corr(d_363(ind_363)));
        gain_sat(i,2) = (Winkler.O2_Dave(ind_good(i))./Winkler.O2_equil(ind_good(i))*100)./nanmean(G363.O2sat_corr(d_363(ind_363)));
    end
end

    %%% Calculate for first calibration cast, just with BCP Winkler data
ind_good = intersect(indGoodBCPdata, indGliderCast1);
    gain_bcp = NaN*ones(length(ind_good),2);
for i = 2:length(ind_good)
    ind_363 = find(G363.depth_interp(d_363) < Winkler.depth(ind_good(i)) + tol & G363.depth_interp(d_363) > Winkler.depth(ind_good(i)) - tol);
    ind_453 = find(G453.depth_interp(d_453) < Winkler.depth(ind_good(i)) + tol & G453.depth_interp(d_453) > Winkler.depth(ind_good(i)) - tol);
    if length(ind_453) > 0
        gain_bcp(i,1) = mean(Winkler.O2_BCP(ind_good(i),:).*(Winkler.dens0(ind_good(i))/1000))./nanmean(G453.O2_corr(d_453(ind_453)));
    end
    if length(ind_363 > 0)
        gain_bcp(i,2) = mean(Winkler.O2_BCP(ind_good(i),:).*(Winkler.dens0(ind_good(i))/1000))./nanmean(G363.O2_corr(d_363(ind_363)));
    end
end

    %%% Calculate for second calibration cast (only BCP Winkler data)
ind_good = intersect(indGoodBCPdata, indGliderCast2);
    gain_bcp2 = NaN*ones(length(ind_good),2);
for i = 2:length(ind_good)
    ind_363 = find(G363.depth_interp(d_363) < Winkler.depth(ind_good(i)) + tol & G363.depth_interp(d_363) > Winkler.depth(ind_good(i)) - tol);
    ind_453 = find(G453.depth_interp(d_453) < Winkler.depth(ind_good(i)) + tol & G453.depth_interp(d_453) > Winkler.depth(ind_good(i)) - tol);
    if length(ind_453) > 0
        gain_bcp2(i,1) = mean(Winkler.O2_BCP(ind_good(i),:).*(Winkler.dens0(ind_good(i))/1000))./nanmean(G453.O2_corr(d_453(ind_453)));
    end
    if length(ind_363 > 0)
        gain_bcp2(i,2) = mean(Winkler.O2_BCP(ind_good(i),:).*(Winkler.dens0(ind_good(i))/1000))./nanmean(G363.O2_corr(d_363(ind_363)));
    end
end

%% Calculate statistics for separate and combined set of Winklers
gain_stats_dave = [nanmean(gain); ...
    nanstd(gain); ...
    2*nanstd(gain)./sqrt(sum(~isnan(gain)))];

gain_stats_bcp = [nanmean(gain_bcp); ...
    nanstd(gain_bcp); ...
    2*nanstd(gain_bcp)./sqrt(sum(~isnan(gain_bcp)))];

gain_stats_all = [nanmean([gain; gain_bcp]); ...
    nanstd([gain; gain_bcp]); ...
    2*nanstd([gain; gain_bcp])./sqrt(sum(~isnan([gain; gain_bcp])))];

gain_stats_C2 = [nanmean(gain_bcp2); ...
    nanstd(gain_bcp2); ...
    2*nanstd(gain_bcp2)./sqrt(sum(~isnan(gain_bcp2)))];

%Visualize mean & errorbar range
    %All three options for C1 overlap so use the combination (stats_all),
    %but yields a different gain correction than the 2nd calibration cast
figure(7); clf
er = errorbar([1:9],[gain_stats_dave(1,1) gain_stats_bcp(1,1) gain_stats_all(1,1) gain_stats_C2(1,1) NaN...
    gain_stats_dave(1,2) gain_stats_bcp(1,2) gain_stats_all(1,2) gain_stats_C2(1,2)],...
    [gain_stats_dave(3,1) gain_stats_bcp(3,1) gain_stats_all(3,1) gain_stats_C2(3,1) NaN...
    gain_stats_dave(3,2) gain_stats_bcp(3,2) gain_stats_all(3,2) gain_stats_C2(3,2)]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xlim([0.5 9.5])
