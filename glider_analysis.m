%% Main script for analyzing OOI Irminger glider data

%% Add paths to data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\METBKA') %SUMO met buoy data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\From_Roo') %Processed glider data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\From_Roo\glider_reprocessed') %Processed glider data

%% Add path to optode lag correction toolbox
addpath(genpath('C:\Users\palevsky\Documents\GitHub\optode-response-time'))

%% Provide boundary layer thickness for glider lag correction
%Value determined from glider_lagAssess (note that to run this, needs data read in from gliders, included later in this script)
tau_in = 66.8878;

%% Identify METBKA files to use for air calibration
% See metbka_assess for analysis of METBKA data to determine coverage of recovered vs telemetered files

Yr5.met = 'deployment0005_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20180608T172109.234000-20190809T080351.721000.nc';
Yr6.met = 'deployment0006_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20190805T152914.351000-20200826T104006.497000.nc';
Yr7.met = 'deployment0007_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20200817T173428.501000-20210819T064814.736000.nc';
Yr8.met = 'deployment0008_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20210812T170400.699000-20211103T161130.307000.nc';

%% Settings used for all air calibration
rhcorr = 1; 
smoothval = 90;
outtolsd = 2;

%% Year 5 (2018 deployment) - corrections for S & pressure, lag correction, and air calibration
load GI05MOAS-GL453.mat %load G453.mat
G453 = T; clear T
load GI05MOAS-GL363.mat % load G363.mat; G363 = G363R;
G363 = T; clear T 

%Interpolation and S & P corrections - Note that internal salinity setting is 35 for 363 and 0 for 453
Yr5.G363 = glider_interpCorrFun(G363, 35);
Yr5.G453 = glider_interpCorrFun(G453, 0);

%For each glider dataset, reshape data into format for Gordon functions, then apply glider_lagCorrectFun
[Yr5.glg363] = glider_reshape(Yr5.G363);
[Yr5.glg363.doxy_lagcorr] = glider_lagCorrectFun(Yr5.glg363, tau_in);

[Yr5.glg453] = glider_reshape(Yr5.G453);
[Yr5.glg453.doxy_lagcorr] = glider_lagCorrectFun(Yr5.glg453, tau_in);

%Air calibration
[Yr5.T_363, Yr5.met_363] = aircalfun(Yr5.G363, Yr5.glg363, -1, Yr5.met, rhcorr, smoothval, outtolsd);
[Yr5.T_453, Yr5.met_453] = aircalfun(Yr5.G453, Yr5.glg453, 1, Yr5.met, rhcorr, smoothval, outtolsd);

clear G363 G363R G453

%% Year 6 (2019 deployment) - corrections for S & pressure, lag correction, and air calibration
% load G525.mat %data until 3 April 2020
load GI05MOAS-GL525.mat
G525 = T; clear T
load GI05MOAS-GL560.mat % load G560.mat %data until Oct. 2019 (O2 sensor failed - other data ok through)
G560 = T; clear T

%Interpolation and S & P corrections - Internal salinity setting is 0
Yr6.G525 = glider_interpCorrFun(G525, 0);
Yr6.G560 = glider_interpCorrFun(G560, 0);

%For each glider dataset, reshape data into format for Gordon functions, then apply glider_lagCorrectFun
[Yr6.glg525] = glider_reshape(Yr6.G525);
[Yr6.glg525.doxy_lagcorr] = glider_lagCorrectFun(Yr6.glg525, tau_in);

[Yr6.glg560] = glider_reshape(Yr6.G560);
[Yr6.glg560.doxy_lagcorr] = glider_lagCorrectFun(Yr6.glg560, tau_in);

%Air calibration
[Yr6.T_560, Yr6.met_560] = aircalfun(Yr6.G560, Yr6.glg560, 1, Yr6.met, rhcorr, smoothval, outtolsd);
[Yr6.T_525, Yr6.met_525] = aircalfun(Yr6.G525, Yr6.glg525, -1, Yr6.met, rhcorr, smoothval, outtolsd);

clear G525 G560

%% Year 7 (2020 deployment) - corrections for S & pressure, lag correction, and air calibration
% load G515.mat %data through 15 June 2021
load GI05MOAS-PG515.mat %slightly different results than earlier processing
G515 = T; clear T
load GI05MOAS-GL365.mat % load G365.mat %data through 21 Nov 2020
G365 = T; clear T

%Interpolation and S & P corrections - Internal salinity setting should be 0
Yr7.G515 = glider_interpCorrFun(G515, 0);
Yr7.G365 = glider_interpCorrFun(G365, 0);

%For each glider dataset, reshape data into format for Gordon functions, then apply glider_lagCorrectFun
[Yr7.glg515] = glider_reshape(Yr7.G515);
[Yr7.glg515.doxy_lagcorr] = glider_lagCorrectFun(Yr7.glg515, tau_in);

[Yr7.glg365] = glider_reshape(Yr7.G365);
[Yr7.glg365.doxy_lagcorr] = glider_lagCorrectFun(Yr7.glg365, tau_in);

%Air calibration
[Yr7.T_515, Yr7.met_515] = aircalfun(Yr7.G515, Yr7.glg515, 1, Yr7.met, rhcorr, smoothval, outtolsd);
[Yr7.T_365, Yr7.met_365] = aircalfun(Yr7.G365, Yr7.glg365, 1, Yr7.met, rhcorr, smoothval, outtolsd);

clear G515 G365

%% Year 8 (2021 deployment) - corrections for S & pressure, lag correction, and air calibration
load GI05MOAS-GL469.mat % load GL469.mat; G469 = T; %11 Aug to 18 Oct 2021
G469 = T; clear T
load GI05MOAS-GL537.mat % load GL537.mat; G537 = T; %11 Aug to 5 Sept 2021
G537 = T; clear T
load GI05MOAS-PG565.mat % load PG565.mat; G565 = T; %30 July 2021 to 10 Jan 2022 - no air cal
G565 = T; clear T

%Interpolation and S & P corrections - Internal salinity setting should be 0
Yr8.G469 = glider_interpCorrFun(G469, 0);
Yr8.G537 = glider_interpCorrFun(G537, 0);
Yr8.G565 = glider_interpCorrFun(G565, 0);

%For each glider dataset, reshape data into format for Gordon functions, then apply glider_lagCorrectFun
[Yr8.glg469] = glider_reshape(Yr8.G469);
[Yr8.glg469.doxy_lagcorr] = glider_lagCorrectFun(Yr8.glg469, tau_in);

[Yr8.glg537] = glider_reshape(Yr8.G537);
[Yr8.glg537.doxy_lagcorr] = glider_lagCorrectFun(Yr8.glg537, tau_in);

[Yr8.glg565] = glider_reshape(Yr8.G565);
[Yr8.glg565.doxy_lagcorr] = glider_lagCorrectFun(Yr8.glg565, tau_in);

%Air calibration
[Yr8.T_469, Yr8.met_469] = aircalfun(Yr8.G469, Yr8.glg469, 1, Yr8.met, rhcorr, smoothval, outtolsd);
[Yr8.T_537, Yr8.met_537] = aircalfun(Yr8.G537, Yr8.glg537, 1, Yr8.met, rhcorr, smoothval, outtolsd);

clear G469 G537 G565

%% Year 2-4 data

%Load data and then apply interpolation and S & P corrections - assuming internal setting is 0 for all
load GI05MOAS-GL493-D00002.mat %Year 4 - 8/9/2017-10/14/2017
    Yr4.G493 = glider_interpCorrFun(T, 0); clear T
load GI05MOAS-GL559-R00001.mat %Year 3 - 7/9/2016-1/8/2017
    Yr3.G559 = glider_interpCorrFun(T, 0); clear T
load GI05MOAS-GL484-D00002.mat %Year 2 - 8/17/2015-10/18/2016
    Yr2.G484 = glider_interpCorrFun(T, 0); clear T
load GI05MOAS-PG528-R00001.mat %Year 2 - 8/17/2015-5/30/2016
    Yr2.G528 = glider_interpCorrFun(T, 0); clear T
load GI05MOAS-GL485-D00002.mat %Year 2 - 8/17/2015-11/22/2015
    Yr2.G485 = glider_interpCorrFun(T, 0); clear T
load GI05MOAS-GL495-D00002.mat %Year 2 - 8/17/2015-5/11/2015
    Yr2.G495 = glider_interpCorrFun(T, 0); clear T
    
%For each glider dataset, reshape data into format for Gordon functions, then apply glider_lagCorrectFun
[Yr4.glg493] = glider_reshape(Yr4.G493);
[Yr4.glg493.doxy_lagcorr] = glider_lagCorrectFun(Yr4.glg493, tau_in);

[Yr3.glg559] = glider_reshape(Yr3.G559);
[Yr3.glg559.doxy_lagcorr] = glider_lagCorrectFun(Yr3.glg559, tau_in);

[Yr2.glg484] = glider_reshape(Yr2.G484);
[Yr2.glg484.doxy_lagcorr] = glider_lagCorrectFun(Yr2.glg484, tau_in);

[Yr2.glg528] = glider_reshape(Yr2.G528);
[Yr2.glg528.doxy_lagcorr] = glider_lagCorrectFun(Yr2.glg528, tau_in);

[Yr2.glg485] = glider_reshape(Yr2.G485);
[Yr2.glg485.doxy_lagcorr] = glider_lagCorrectFun(Yr2.glg485, tau_in);

[Yr2.glg495] = glider_reshape(Yr2.G495);
[Yr2.glg495.doxy_lagcorr] = glider_lagCorrectFun(Yr2.glg495, tau_in);

%% Identify and flag outliers and spikes in data, and regrid on pressure surfaces (following same approach as wfp_analysis)
%tolerances for spikes (2 x those chosen for WFP for oxy and 20x for
%temp and sal, given higher variability near surface)
oxyspike = 6;
tempspike = 1;
salspike = 0.1;
oxymin = 240;
oxymin2 = 180;
oxymax = 360;
tempmin = 3;
tempmax = 12;
salmin = 32;
salmax = 35.1;

%Select depth resolution and smoothing - current setting is 1 m resolution
%w/ 5-m smoothing
pres_grid = [1:1:1000];
smval = 5; %points to smooth over

[Yr5.glg363,Yr5.glg363.fract_flag] = glider_flag_regrid(Yr5.glg363, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr5.glg453,Yr5.glg453.fract_flag] = glider_flag_regrid(Yr5.glg453, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr6.glg525,Yr6.glg525.fract_flag] = glider_flag_regrid(Yr6.glg525, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr6.glg560,Yr6.glg560.fract_flag] = glider_flag_regrid(Yr6.glg560, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr7.glg515,Yr7.glg515.fract_flag] = glider_flag_regrid(Yr7.glg515, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr7.glg365,Yr7.glg365.fract_flag] = glider_flag_regrid(Yr7.glg365, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr8.glg469,Yr8.glg469.fract_flag] = glider_flag_regrid(Yr8.glg469, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr8.glg537,Yr8.glg537.fract_flag] = glider_flag_regrid(Yr8.glg537, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr8.glg565,Yr8.glg565.fract_flag] = glider_flag_regrid(Yr8.glg565, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);

[Yr4.glg493,Yr4.glg493.fract_flag] = glider_flag_regrid(Yr4.glg493, oxymin2, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr3.glg559,Yr3.glg559.fract_flag] = glider_flag_regrid(Yr3.glg559, oxymin2, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr2.glg484,Yr2.glg484.fract_flag] = glider_flag_regrid(Yr2.glg484, oxymin2, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr2.glg528,Yr2.glg528.fract_flag] = glider_flag_regrid(Yr2.glg528, oxymin2, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr2.glg485,Yr2.glg485.fract_flag] = glider_flag_regrid(Yr2.glg485, oxymin2, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);
[Yr2.glg495,Yr2.glg495.fract_flag] = glider_flag_regrid(Yr2.glg495, oxymin2, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval);

%% Combine glider datasets
%Note - Year 8 glider 537 removed from further analysis because all salinity data flagged as outlier, and only lasted a few weeks anyway
glgmerge{1} = Yr5.glg363;
glgmerge{2} = Yr5.glg453;
glgmerge{3} = Yr6.glg525;
glgmerge{4} = Yr6.glg560;
glgmerge{5} = Yr7.glg515;
glgmerge{6} = Yr7.glg365;
glgmerge{7} = Yr8.glg469;
glgmerge{8} = Yr8.glg565;

glgmerge{9} = Yr2.glg495;
glgmerge{10} = Yr2.glg485;
glgmerge{11} = Yr2.glg528;
glgmerge{12} = Yr2.glg484; %no usable data
glgmerge{13} = Yr3.glg559;
glgmerge{14} = Yr4.glg493;

glidertitles = [{'Glider 363, Year 5','Glider 453, Year 5','Glider 525, Year 6','Glider 560, Year 6',...
    'Glider 515, Year 7','Glider 365, Year 7','Glider 469, Year 8','Glider 565, Year 8'...
    'Glider 495, Year 2', 'Glider 485, Year 2', 'Glider 528, Year 2', 'Glider 484, Year 2', 'Glider 559, Year 3', 'Glider 493, Year 4'}];

%% Add air cal output to combined glider datasets
glgmerge{1}.Taircal = Yr5.T_363;
glgmerge{2}.Taircal = Yr5.T_453;
glgmerge{3}.Taircal = Yr6.T_525;
glgmerge{4}.Taircal = Yr6.T_560;
glgmerge{5}.Taircal = Yr7.T_515;
glgmerge{6}.Taircal = Yr7.T_365;
glgmerge{7}.Taircal = Yr8.T_469;

%% Add met output to combined glider datasets
glgmerge{1}.met = Yr5.met_363;
glgmerge{2}.met = Yr5.met_453;
glgmerge{3}.met = Yr6.met_525;
glgmerge{4}.met = Yr6.met_560;
glgmerge{5}.met = Yr7.met_515;
glgmerge{6}.met = Yr7.met_365;
glgmerge{7}.met = Yr8.met_469;

%% Assessment histograms
figure(1); clf
subplot(311)
for i = 1:length(glgmerge)
    histogram(glgmerge{i}.doxy_lagcorr_grid(:)); hold on;
end
legend(glidertitles(6:14),'location','NW')
xlabel('Dissolved oxygen, uncalibrated (\mumol/kg)');
xlim([oxymin2 oxymax])

subplot(312)
for i = 1:length(glgmerge)
    histogram(glgmerge{i}.temp_grid(:)); hold on;
end
legend(glidertitles,'location','NW')
xlabel('Temperature (^oC)');
xlim([tempmin tempmax])

subplot(313)
for i = 1:length(glgmerge)
    histogram(glgmerge{i}.sal_grid(:)); hold on;
end
legend(glidertitles,'location','NW')
xlabel('Salinity (PSU)');
xlim([salmin + 2.5 salmax])

%% Assess sensitivity to selection of air measurement value from distribution within each surface interval
figure(2); clf
n = 7; %number of gliders in glgmerge with aircal data
clear A B
for ii = 1:n
    subplot(2,n,ii)
    for i = 1:height(glgmerge{ii}.Taircal)
       plot([0.2:0.05:0.8], glgmerge{ii}.Taircal.air_meas_dist(i,4:16) - glgmerge{ii}.Taircal.air_meas_dist(i,10)); hold on;
       A(i,:) = glgmerge{ii}.Taircal.air_meas_dist(i,4:16) - glgmerge{ii}.Taircal.air_meas_dist(i,10);
    end
    plot([0.2:0.05:0.8], nanmean(A),'k-','linewidth',5); hold on;
    title([glidertitles{ii}])
    xlabel('Percentile'); ylabel('\DeltaO_{2,a}^{meas} - median');
    ylim([-3 3])

    subplot(2,n,ii+n)
    for i = 1:height(glgmerge{ii}.Taircal)
       plot([0.225:0.05:0.8], diff(glgmerge{ii}.Taircal.air_meas_dist(i,4:16))); hold on;
       B(i,:) = diff(glgmerge{ii}.Taircal.air_meas_dist(i,4:16));
    end
    plot([0.225:0.05:0.8], nanmean(B),'k-','linewidth',5); hold on;
    title([glidertitles{ii}])
    xlabel('Percentile'); ylabel('d\DeltaO_{2,a}^{meas}');
    ylim([0 1])
    derivmean(ii,:) = nanmean(B);
    [~,id(ii)] = min(nanmean(B));
end

%% Recalculate glider air corrections to test sensitivity to selection of air measurement value from distribution & lag corrections to splashed water
percentiles = [0.05:0.05:0.95];
for ii = 1:n
    %Sensitivity to selection of air measurement value from distribution
    for i = 1:length(percentiles)
        LM = fitlm(glgmerge{ii}.Taircal.ml_o2sat,glgmerge{ii}.Taircal.air_meas_dist(:,i));
        slope(ii,i) = table2array(LM.Coefficients(2,1));
        slope_SE(ii,i) = table2array(LM.Coefficients(2,2));
        intercept(ii,i) = table2array(LM.Coefficients(1,1));
    end
    %Compare with median val using uncorrected splash water
    LM2 = fitlm(glgmerge{ii}.Taircal.ml_o2conc_nocorr./gsw_O2sol_SP_pt(glgmerge{ii}.Taircal.ml_sal, glgmerge{ii}.Taircal.ml_tem)*100,...
        glgmerge{ii}.Taircal.air_meas);
    slope_nocorr(ii) = table2array(LM2.Coefficients(2,1));
    slope_nocorr_SE(ii) = table2array(LM2.Coefficients(2,2));
    intercept_nocorr(ii) = table2array(LM2.Coefficients(1,1));
    %Assess quantitative difference between corrected and uncorrected splash water
    watercorr = glgmerge{ii}.Taircal.ml_o2conc - glgmerge{ii}.Taircal.ml_o2conc_nocorr;
    splash_o2conccorr(ii,1) = nanmean(abs(watercorr));
    splash_o2conccorr(ii,2) = nanstd(abs(watercorr));
end

%For each glider, mean across slopes from the 30th-70th percentile
slope_mean = mean(slope(:,6:14),2);
%Range of slopes using the 30th-70th percentiles
R = max(slope(:,6:14),[],2) - min(slope(:,6:14),[],2);
%Maximum of the standard errors of the slope fits across the 30th-70th percentiles
SE = max(slope_SE(:,6:14),[],2);

slope_err = max(R,SE);

%% Test different outcomes for range of slopes
slope_vals = [0.2:0.05:0.6];
percent_list = [12,8,10];
C = parula(length(slope_vals));
slope_pick = mean(slope_mean);

figure(3); clf
for ii = [1,3,5,7];
    for i = 1:length(slope_vals)
        for j = 1:length(percent_list)
            A(:,i,j) = (glgmerge{ii}.Taircal.air_meas_dist(:,percent_list(j))-slope_vals(i).*glgmerge{ii}.Taircal.ml_o2sat)./(1-slope_vals(i)); %eqn 5 in Nicholson and Feen 2017
            plot(glgmerge{ii}.Taircal.ml_daten, movmean(glgmerge{ii}.Taircal.met_o2sat./A(:,i,j), 60),'-','linewidth',j,'color',C(i,:)); hold on; %account for surface water splashing
        end
    end
    corr_to_plot = (glgmerge{ii}.Taircal.air_meas_dist(:,10)-slope(ii,10).*glgmerge{ii}.Taircal.ml_o2sat)./(1-slope(ii,10));
    corr_to_plot_L = (glgmerge{ii}.Taircal.air_meas_dist(:,10)-(slope(ii,10)-slope_err(ii)).*glgmerge{ii}.Taircal.ml_o2sat)./(1-(slope(ii,10)-slope_err(ii)));
    corr_to_plot_H = (glgmerge{ii}.Taircal.air_meas_dist(:,10)-(slope(ii,10)+slope_err(ii)).*glgmerge{ii}.Taircal.ml_o2sat)./(1-(slope(ii,10)+slope_err(ii)));
    corr_to_plot2 = (glgmerge{ii}.Taircal.air_meas_dist(:,10)-slope_pick.*glgmerge{ii}.Taircal.ml_o2sat)./(1-slope_pick);
    plot(glgmerge{ii}.Taircal.ml_daten, movmean(glgmerge{ii}.Taircal.met_o2sat./corr_to_plot, 60),'k-','linewidth',4); hold on;
    plot(glgmerge{ii}.Taircal.ml_daten, movmean(glgmerge{ii}.Taircal.met_o2sat./corr_to_plot_L, 60),'k-','linewidth',1); hold on;
    plot(glgmerge{ii}.Taircal.ml_daten, movmean(glgmerge{ii}.Taircal.met_o2sat./corr_to_plot_H, 60),'k-','linewidth',1); hold on;
    %plot(glgmerge{ii}.Taircal.ml_daten, movmean(glgmerge{ii}.Taircal.met_o2sat./corr_to_plot2, 60),'m-','linewidth',4); hold on;
    clear A
end
datetick('x',2,'keeplimits')
title('Glider gain corrections: Sensitivity to air-water slopes');

%% Synthesize all air vs water empirical slope data

p_ind = 10; %index for 50th percentile
ftsz = 10;
lnw = 1.5;
mrkr = 10;
spc = 2;

figure(4); clf
for ii = 1:7
subplot(4,2,ii)
hold all; box on;
xrng = [floor(min(glgmerge{ii}.Taircal.ml_o2sat)) - spc, ceil(max(glgmerge{ii}.Taircal.ml_o2sat)) + spc];
yrng = [floor(min(glgmerge{ii}.Taircal.air_meas)) - spc, ceil(max(glgmerge{ii}.Taircal.air_meas)) + spc];
    plot(glgmerge{ii}.Taircal.ml_o2sat,glgmerge{ii}.Taircal.air_meas,'.','MarkerSize',mrkr);
    plot([xrng],slope(ii,p_ind).*xrng+intercept(ii,p_ind),'-k','LineWidth',lnw);
    xlabel('\DeltaO_{2,w}^{meas}');
    ylabel('\DeltaO_{2,a}^{meas}');
    title([glidertitles(ii) ' air vs. surface water'],'Fontsize',ftsz);
    xlim(xrng)
    ylim(yrng)
    text(xrng(1) + [(xrng(2)-xrng(1))/20], yrng(2) - [(yrng(2)-yrng(1))/10], ['Slope = ' num2str(slope(ii,p_ind),2) ' ' char(177) ' ' num2str(slope_err(ii),1)])
end

%% Synthesize all time series gain correction
tgap = 5;
lnw1 = 2;
lnw2 = 2.5;

figure(5); clf
for ii = 1:7
subplot(4,2,ii)
ax2 = gca; cols = ax2.ColorOrder;
hold all; box on;
    air_corr = (glgmerge{ii}.Taircal.air_meas-slope(ii,p_ind).*glgmerge{ii}.Taircal.ml_o2sat)./(1-slope(ii,p_ind)); %eqn 5 in Nicholson and Feen 2017    
    h1 = plot(glgmerge{ii}.Taircal.ml_daten,glgmerge{ii}.Taircal.ml_o2sat,'-','LineWidth',lnw1,'Color',cols(1,:));
    h2 = plot(glgmerge{ii}.Taircal.air_daten,glgmerge{ii}.Taircal.air_meas,'-','LineWidth',lnw1,'Color',cols(3,:));
    h3 = plot(glgmerge{ii}.Taircal.air_daten,air_corr,'-','MarkerSize',mrkr,'LineWidth',lnw2,'Color',cols(4,:));
    h4 = plot(glgmerge{ii}.met.daten,glgmerge{ii}.met.O2satcorr,'-','LineWidth',lnw2,'Color','k');

title(glidertitles(ii),'Fontsize',ftsz+2);
ylabel('Oxygen saturation (%)')
xlim([min(glgmerge{ii}.Taircal.air_daten) - tgap max(glgmerge{ii}.Taircal.air_daten) + tgap])
datetick('x',2,'keeplimits');
if ii == 1
    legend([h1 h2 h3 h4],'\DeltaO_{2,w}^{meas}','\DeltaO_{2,a}^{meas}','\DeltaO_{2,a}^{splash corr}','\DeltaO_{2}^{met}',...
        'location','northeast','orientation','horizontal');
end
end

figure(6); clf
ax = gca; cols = ax.ColorOrder;
for ii = 1:7
subplot(4,2,ii)
hold all; box on;
    air_corr = (glgmerge{ii}.Taircal.air_meas-slope(ii,p_ind).*glgmerge{ii}.Taircal.ml_o2sat)./(1-slope(ii,p_ind)); %eqn 5 in Nicholson and Feen 2017    
    plot(glgmerge{ii}.Taircal.met_o2sat,air_corr,'.','MarkerSize',mrkr,'Color',cols(ii,:)); hold on;
    title([glidertitles(ii) 'corrected air vs. MET data'],'Fontsize',ftsz);
    xlabel('\DeltaO_{2}^{met}');
    ylabel('\DeltaO_{2,a}^{splash corr}');
end
subplot(4,2,8)
hold all; box on;
for ii = 1:7
    air_corr = (glgmerge{ii}.Taircal.air_meas-slope(ii,p_ind).*glgmerge{ii}.Taircal.ml_o2sat)./(1-slope(ii,p_ind)); %eqn 5 in Nicholson and Feen 2017    
    plot(glgmerge{ii}.Taircal.met_o2sat,air_corr,'.','MarkerSize',mrkr,'Color',cols(ii,:)); hold on;
end
title({'All Gliders' 'corrected air vs. MET data'},'Fontsize',ftsz);
xlabel('\DeltaO_{2}^{met}');
ylabel('\DeltaO_{2,a}^{splash corr}');

%% Visualize lag corrected, gridded glider oxygen data from near surface
figure(7); clf
for ii = 1:8
subplot(4,2,ii)
    indplot = find(sum(~isnan(glgmerge{ii}.doxy_lagcorr_grid(1:50,:))) > 0);
    imagesc(glgmerge{ii}.doxy_lagcorr_grid(1:50, indplot))
    colorbar
    title(glidertitles(ii))
    ylabel('Depth (m)','Fontsize',8)
end

%%
figure(7); clf
for ii = 1:14
subplot(4,4,ii)
    indplot = find(sum(~isnan(glgmerge{ii}.doxy_lagcorr_grid(1:100,:))) > 0);
    imagesc((glgmerge{ii}.doxy_lagcorr_grid(:, indplot)))
    colorbar
    title(glidertitles(ii))
    ylabel('Depth (m)','Fontsize',8)
end

%% Save output
save('glider_output.mat', 'glgmerge', '-v7.3')