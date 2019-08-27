%% Code to process CTD data from Irminger 6
% Due to failure of rosette, we only have sensor data, no bottle data
% Added an OOI Aanderaa optode (a 4831, serial #502) onto the CTD
% to collect data and be able to characterize in post-cruise measurements

% Note that there is also an issue with the primary CTD data, so need to
% make sure that this is using secondary T and S data (T2, not T1 - last
% year used secondary salinity already as Sal)

%% INPUT
% volts: raw optode output from CTD in volts
% tempc: temperature (deg C) from CTD
% salin: salinity from CTD
% press: pressure from CTD

%% Coefficients for processing data from Aanderaa optode SN502
%Coefficients for converting from voltage to calphase
% % %  0-5V  Output 1: CalPhase(Deg)             1.664V, use scaling coef.
% % % A:=1.00000                                    0E+01 B:=1.200000E+01
% % % 4-20mA Output 1: CalPhase(Deg)             9.324mA, use scaling coef.
% % % A:=-5.0000                                    00E+00 B:=3.750000E+00
A = 10; B = 12;
%calphase = A + B.*volts; %assuming linear equation to apply coefficients

foilcoeff_502 = [2.82567E-3 1.20716E-4 2.4593E-6 2.30757E2 -3.09502E-1 -5.60627E1 4.5615E0];
conccoeff = [-1.28596 1.039998];

%% Check of calphase conversion with factory calibration data on calsheet
%Values from factory 2-point calibration, 2018-August-17
phasereading_tests = [33.01 60.96];
T_tests = [9.96 22.36];

%Check values from 2-point calibration
[optode_uM, optode_umolkg] = aaoptode_sternvolmer(foilcoeff_502, phasereading_tests, T_tests, [0 0], [0 0]);

%% Read data in from Excel compilation of CTD casts
path = ['C:/Users/palevsky/Dropbox/Irminger6/'];
file = ['Irminger6_ASC_Continuous_CTD_Data.xlsx'];
filename= [path file]
castnums = [3:6,8:13];
castnames = {'C003: HYPM, no bottles', 'C004: Glider deployment, no bottles', 'C005: FLMA, no bottles', 'C006: FLMB, no bottles', 'C008: Glider cast, anchor release bottles', 'C009: SUMO/HYPM, anchor release bottles', 'C010: FLMB, anchor release bottles',...
    'C011: HYPM/SUMO/Glider', 'C012: FLMB', 'C013: FLMA'};
%initialize cast
cast{1} = [];
for i = 1:length(castnums)
    ASCContinuousCTDData = xlsread(filename,i);
    cast{castnums(i)}.Pres= ASCContinuousCTDData (:,1);
    cast{castnums(i)}.T1= ASCContinuousCTDData (:,2); %T090C Temperature (ITS-90, deg C)
    cast{castnums(i)}.T2= ASCContinuousCTDData (:,3); %T190C Temperature (ITS-90, deg C)
    cast{castnums(i)}.C1= ASCContinuousCTDData (:,4); %COS/m Conductivity (S/m)
    cast{castnums(i)}.C2= ASCContinuousCTDData (:,5); %C1S/m Conductivity (S/m)
    cast{castnums(i)}.OR= ASCContinuousCTDData (:,6); %sbeox0V: Oxygen raw, SBE 43 [V]
    cast{castnums(i)}.F = ASCContinuousCTDData (:,7); %flECO-AFL: Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]
    cast{castnums(i)}.Turb = ASCContinuousCTDData (:,8); %turbWETntu0: Turbidity, WET Labs ECO [NTU]
    cast{castnums(i)}.O2_volt = ASCContinuousCTDData(:,11); %Voltage output from Aanderaa optode
    cast{castnums(i)}.D= ASCContinuousCTDData (:,12); %Depth, meters
        [~, cast{castnums(i)}.maxindex] = max (cast{castnums(i)}.D);
    cast{castnums(i)}.Sal= ASCContinuousCTDData (:,14); %sal11: Salinity, Practical [PSU]
    cast{castnums(i)}.O2 = ASCContinuousCTDData (:,15)/(.022391); %Sbeox0ML/L converted to micromoles/L
    cast{castnums(i)}.SvCM = ASCContinuousCTDData (:,16); 
    cast{castnums(i)}.SA = gsw_SA_from_SP(cast{castnums(i)}.Sal,cast{castnums(i)}.Pres,60,19);
    cast{castnums(i)}.CT = gsw_CT_from_t(cast{castnums(i)}.SA,cast{castnums(i)}.T2,cast{castnums(i)}.Pres);
    cast{castnums(i)}.rho0 = gsw_rho(cast{castnums(i)}.SA,cast{castnums(i)}.CT,0);
    cast{castnums(i)}.O2sol = gsw_O2sol(cast{castnums(i)}.SA,cast{castnums(i)}.CT,0,60,19);
    cast{castnums(i)}.O2sol2 = gsw_O2sol_SP_pt(cast{castnums(i)}.Sal,cast{castnums(i)}.T2);
    cast{castnums(i)}.aou = cast{castnums(i)}.O2sol - cast{castnums(i)}.O2 ;
    %Aanderaa processing calculations
    cast{castnums(i)}.calphase = cast{castnums(i)}.O2_volt*B + A;
    [optode_uM, ~] = aaoptode_sternvolmer(foilcoeff_502, cast{castnums(i)}.calphase, cast{castnums(i)}.T2, cast{castnums(i)}.Sal, cast{castnums(i)}.Pres);
    optode_uM = conccoeff(1) + conccoeff(2).*optode_uM;
    cast{castnums(i)}.O2corr = aaoptode_salpresscorr(optode_uM, cast{castnums(i)}.T2, cast{castnums(i)}.Sal, cast{castnums(i)}.Pres, 0);
end

%% Read in Irminger-6 Winkler data
[Winkler6_casts] = xlsread('C:/Users/palevsky/Dropbox/Irminger6/Oxygen data/Irminger6_Winkler.xlsx',2);
Winkler6.cast = Winkler6_casts(:,1);
Winkler6.depth = Winkler6_casts(:,8);
Winkler6.pres = Winkler6_casts(:,12);
Winkler6.T = Winkler6_casts(:,9);
Winkler6.S = Winkler6_casts(:,11);
Winkler6.O2_dave = Winkler6_casts(:,18);
Winkler6.O2_dave_flag = Winkler6_casts(:,19);
Winkler6.O2_bcp = Winkler6_casts(:,29);
Winkler6.O2_bcp_flag = Winkler6_casts(:,32);
    ind_bcp = find(Winkler6.O2_bcp_flag > 0);
    ind_dave = find(Winkler6.O2_dave_flag > 0 & Winkler6.O2_dave_flag < 3);

%% Plot O2 data from casts
tol = 1.5; %tolerance for extracting cast sensor data aligned with bottle samples

Winkler6.optode_align_down = NaN*ones(size(Winkler6.depth));
Winkler6.optode_align_down_std = NaN*ones(size(Winkler6.depth));
Winkler6.optode_align_up = NaN*ones(size(Winkler6.depth));
Winkler6.optode_align_up_std = NaN*ones(size(Winkler6.depth));

for i = 1:length(castnums)
    figure(castnums(i)); clf
%subplot(121)
plot(cast{castnums(i)}.O2corr, cast{castnums(i)}.D, 'k.'); hold on;
plot(cast{castnums(i)}.O2corr(cast{castnums(i)}.maxindex:end), cast{castnums(i)}.D(cast{castnums(i)}.maxindex:end), 'r.'); hold on;
    indWink = find(Winkler6.cast == castnums(i));
    for j = 1:length(indWink)
        indCastAlignDown = find(cast{castnums(i)}.D(1:cast{castnums(i)}.maxindex) > Winkler6.depth(indWink(j)) - tol...
            & cast{castnums(i)}.D(1:cast{castnums(i)}.maxindex) < Winkler6.depth(indWink(j)) + tol);
        indCastAlignUp = find(cast{castnums(i)}.D(cast{castnums(i)}.maxindex:end) > Winkler6.depth(indWink(j)) - tol...
            & cast{castnums(i)}.D(cast{castnums(i)}.maxindex:end) < Winkler6.depth(indWink(j)) + tol);
        Winkler6.optode_align_up(indWink(j)) = mean(cast{castnums(i)}.O2corr(indCastAlignUp + (cast{castnums(i)}.maxindex -1)));
        Winkler6.optode_align_up_std(indWink(j)) = std(cast{castnums(i)}.O2corr(indCastAlignUp + (cast{castnums(i)}.maxindex -1)));
        Winkler6.optode_align_down(indWink(j)) = mean(cast{castnums(i)}.O2corr(indCastAlignDown));
        Winkler6.optode_align_down_std(indWink(j)) = std(cast{castnums(i)}.O2corr(indCastAlignDown));
    end
plot(Winkler6.optode_align_up(intersect(indWink, ind_bcp)), Winkler6.depth(intersect(indWink, ind_bcp)), 'k.','markersize',15); hold on;
plot(Winkler6.O2_bcp(intersect(indWink, ind_bcp)), Winkler6.depth(intersect(indWink, ind_bcp)), 'c.','markersize',15); hold on;
plot(Winkler6.O2_dave(intersect(indWink, ind_dave)), Winkler6.depth(intersect(indWink, ind_dave)), 'm.','markersize',15); hold on;
axis ij
xlabel('Oxygen (\muM) from Aanderaa optode')
ylabel('Depth (m)')
legend('Downcast','Upcast','Depth-aligned upcast point','BCP Winkler','Wellwood Winkler')
title(castnames(i))

% subplot(122)
% %plot(cast{castnums(i)}.Sal, cast{castnums(i)}.D, 'r.'); hold on; 
% plot(gsw_pt_from_CT(cast{castnums(i)}.SA,cast{castnums(i)}.CT), cast{castnums(i)}.D, 'k.'); hold on; 
% axis ij
% xlabel('Potential temperature (deg C)')
% ylabel('Depth, m')
end

%% Plot data from both HYPM casts before pylon

figure(1); clf
    i = 1; %First HYPM cast
plot(cast{castnums(i)}.O2corr(cast{castnums(i)}.maxindex:end), cast{castnums(i)}.D(cast{castnums(i)}.maxindex:end), 'r-'); hold on;
    i = 6; %HYPM/SUMO cast
plot(cast{castnums(i)}.O2corr(cast{castnums(i)}.maxindex:end), cast{castnums(i)}.D(cast{castnums(i)}.maxindex:end), 'b-'); hold on;
axis ij
xlabel('Oxygen, \muM')
ylabel('Depth, m')
legend('C003: First HYPM cast','C009: HYPM/SUMO cast','location','northwest')

%% Plot all Winkler data against up-cast depth-aligned optode data
M = 8;    
ind_align_good = find(Winkler6.optode_align_up_std < 0.75); % & Winkler6.depth > 100 & Winkler6.depth < 2400);
    figure(2); clf
subplot(121)
scatter(Winkler6.optode_align_up(ind_bcp), Winkler6.O2_bcp(ind_bcp), [], Winkler6.cast(ind_bcp),'filled'); hold on;
plot(Winkler6.optode_align_up(intersect(ind_align_good,ind_bcp)), Winkler6.O2_bcp(intersect(ind_align_good,ind_bcp)), 'r.','markersize',M); hold on;
cbar = colorbar; colormap(cmocean('haline',6));
ylabel(cbar,'Cast #');
%caxis([0 1])
plot([270:310],[270:310],'r--'); hold on;
ylabel('BCP Winkler O_2 (\mumol/kg)'); xlabel('Optode upcast O_2 (\mumol/kg)');
    offset = (Winkler6.O2_bcp(intersect(ind_align_good,ind_bcp)) - Winkler6.optode_align_up(intersect(ind_align_good,ind_bcp)));
    ratio = (Winkler6.O2_bcp(intersect(ind_align_good,ind_bcp))./Winkler6.optode_align_up(intersect(ind_align_good,ind_bcp)));
text(275,305,{['offset = ' num2str(mean(offset)) ' +/- ' num2str(std(offset))]})
text(275,300,{['ratio = ' num2str(mean(ratio)) ' +/- ' num2str(std(ratio)*2/sqrt(length(intersect(ind_align_good,ind_bcp))))]})
title('Irminger6 - calibration of CTD package Aanderaa optode')

subplot(122)
scatter(Winkler6.optode_align_up(ind_dave), Winkler6.O2_dave(ind_dave), [], Winkler6.cast(ind_dave),'filled'); hold on;
plot(Winkler6.optode_align_up(intersect(ind_align_good,ind_dave)), Winkler6.O2_dave(intersect(ind_align_good,ind_dave)), 'r.','markersize',M); hold on;
cbar = colorbar; colormap(cmocean('haline',6));
ylabel(cbar,'Cast #');
%caxis([0 1])
plot([270:310],[270:310],'r--'); hold on;
ylabel('Wellwood Winkler O_2 (\mumol/kg)'); xlabel('Optode upcast O_2 (\mumol/kg)');
    offset = (Winkler6.O2_dave(intersect(ind_align_good,ind_dave)) - Winkler6.optode_align_up(intersect(ind_align_good,ind_dave)));
    ratio = (Winkler6.O2_dave(intersect(ind_align_good,ind_dave))./Winkler6.optode_align_up(intersect(ind_align_good,ind_dave)));
text(275,305,{['offset = ' num2str(mean(offset)) ' +/- ' num2str(std(offset))]})
text(275,300,{['ratio = ' num2str(mean(ratio)) ' +/- ' num2str(std(ratio)*2/sqrt(length(intersect(ind_align_good,ind_dave))))]})
title('Irminger6 - calibration of CTD package Aanderaa optode')
