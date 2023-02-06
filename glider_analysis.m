%% Main script for analyzing OOI Irminger glider data
%Mirrors analysis in aircal_all, with additions for further processing

%% Add paths to data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\METBKA') %SUMO met buoy data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\From_Roo') %Processed glider data

%% Identify METBKA files to use for air calibration
% See metbka_assess for analysis of METBKA data to determine coverage of recovered vs telemetered files

Yr5.met = 'deployment0005_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20180608T172109.234000-20190809T080351.721000.nc';
Yr6.met = 'deployment0006_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20190805T152914.351000-20200826T104006.497000.nc';
Yr7.met = 'deployment0007_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20200817T173428.501000-20210819T064814.736000.nc';
Yr8.met = 'deployment0008_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20210812T170400.699000-20211103T161130.307000.nc';

%% Year 5 (2018 deployment) - initial air calibration & corrections for S & pressure
load G453.mat
load G363.mat
G363 = G363R;

%Air calibration
rhcorr = 1; mindateplot = datenum(2018,6,1);
[Yr5.T_363, Yr5.med_gain_363] = aircalfun(G363, 'Glider 363', -1, Yr5.met, mindateplot, rhcorr, mindateplot, datenum(2019,6,1));
[Yr5.T_453, Yr5.med_gain_453] = aircalfun(G453, 'Glider 453', 1, Yr5.met, mindateplot, rhcorr, mindateplot, datenum(2019,2,1));

%Interpolation and S & P corrections - Note that internal salinity setting is 35 for 363 and 0 for 453
Yr5.G363 = glider_interpCorrFun(G363, 35);
Yr5.G453 = glider_interpCorrFun(G453, 0);

clear G363 G363R G453

%% Year 6 (2019 deployment) - initial air calibration & corrections for S & pressure
load G525.mat %data until 3 April 2020
load G560.mat %data until Oct. 2019 (O2 sensor failed - other data ok through)

rhcorr = 1; mindateplot = datenum(2019,8,1);

[Yr6.T_560, Yr6.med_gain_560] = aircalfun(G560, 'Glider 560', 1, Yr6.met, mindateplot, rhcorr, mindateplot, datenum(2019,10,15));
[Yr6.T_525, Yr6.med_gain_525] = aircalfun(G525, 'Glider 525', -1, Yr6.met, mindateplot, rhcorr, mindateplot, datenum(2020,4,15));

%Interpolation and S & P corrections - Internal salinity setting is 0
Yr6.G525 = glider_interpCorrFun(G525, 0);
Yr6.G560 = glider_interpCorrFun(G560, 0);

clear G525 G560

%% Year 7 (2020 deployment) - initial air calibration & corrections for S & pressure
load G515.mat %data through 15 June 2021
load G365.mat %data through 21 Nov 2020

rhcorr = 1; mindateplot = datenum(2020,8,1);

[Yr7.T_515, Yr7.med_gain_515] = aircalfun(G515, 'Glider 515', 1, Yr7.met, mindateplot, rhcorr, mindateplot, datenum(2021,6,30));
[Yr7.T_365, Yr7.med_gain_365] = aircalfun(G365, 'Glider 365', 1, Yr7.met, mindateplot, rhcorr, mindateplot, datenum(2020,11,30));

%Interpolation and S & P corrections - Internal salinity setting should be 0
Yr7.G515 = glider_interpCorrFun(G515, 0);
Yr7.G365 = glider_interpCorrFun(G365, 0);

clear G515 G365

%% Year 8 (2021 deployment) - initial air calibration & corrections for S & pressure
load GL469.mat; G469 = T; %11 Aug to 18 Oct 2021
load GL537.mat; G537 = T; %11 Aug to 5 Sept 2021
load PG565.mat; G565 = T; %30 July 2021 to 10 Jan 2022 - no air cal

rhcorr = 1; mindateplot = datenum(2021,7,30);
[Yr8.T_469, Yr7.med_gain_515] = aircalfun(G469, 'Glider 469', 1, Yr8.met, mindateplot, rhcorr, mindateplot, datenum(2021,10,20));
[Yr8.T_537, Yr7.med_gain_537] = aircalfun(G537, 'Glider 537', 1, Yr8.met, mindateplot, rhcorr, datenum(2021,8,11), datenum(2021,8,24));
%[Yr8.T_565, Yr7.med_gain_565] = aircalfun(G565, 'Glider 565', 1, Yr8.met, mindateplot, rhcorr, mindateplot, datenum(2022,1,12)); %doesn't actuallyhave air cal, so this just shows that air & surf water match

%Interpolation and S & P corrections - Internal salinity setting should be 0
Yr8.G469 = glider_interpCorrFun(G469, 0);
Yr8.G537 = glider_interpCorrFun(G537, 0);
Yr8.G565 = glider_interpCorrFun(G565, 0);

clear G469 G537 G565 T

%% Lag correction - only for years 5 and 6 (no paired profiles for lag assessment in years 7 and 8)
% Correct for historesis effects between up and down profiles using Gordon et al. 2020 approach, as in wfp_lag
addpath(genpath('C:\Users\Palevsky\Documents\GitHub\optode-response-time'))

% only one set of paired dive and climb, aligns with OSM 2020 analysis and 363 set b
[Yr5.Lag453] = glider_lagAssessFun(Yr5.G453, 26000, 39400, 300, 'Glider 453, Year 5', [43:51]);

% two sets of paired dive and climb, second aligns with cal cast analyzed for OSM 2020 lag analysis
[Yr5.Lag363a] = glider_lagAssessFun(Yr5.G363, 10000, 210000, 300, 'Glider 363, Year 5', [13:20]); %just first set
[Yr5.Lag363b] = glider_lagAssessFun(Yr5.G363, 10000, 210000, 300, 'Glider 363, Year 5', [42:46]); %just second set

%only one set of paired dive and climb
[Yr6.Lag525] = glider_lagAssessFun(Yr6.G525, 800, 68200, 300, 'Glider 525, Year 6', [7:17]);

%long period of paired dive & climb (I think was a mistake...) after Glider 560 deployment
[Yr6.Lag560a] = glider_lagAssessFun(Yr6.G560, 2000, 109000, 300, 'Glider 560, Year 6', [10:42]); %a = lag just for first set, before dives reduced to only 500 m
[Yr6.Lag560b] = glider_lagAssessFun(Yr6.G560, 2000, 109000, 300, 'Glider 560, Year 6', [56:108]); %b = lag just for second set, after dives reduced to only 500 m

%% Plot lag correction histograms & summary stats
    figure(10); clf
subplot(221)
histogram(Yr5.Lag453.thickness,[0:25:200])
title(['Glider 453, Year 5; ' num2str(nanmean(Yr5.Lag453.thickness),3) ' +/- ' num2str(nanstd(Yr5.Lag453.thickness),2)])
xlabel('Boundary layer thickness, \mum')

subplot(222)
histogram([Yr5.Lag363a.thickness Yr5.Lag363b.thickness],[0:25:200])
title(['Glider 363, Year 5; ' num2str(nanmean([Yr5.Lag363a.thickness Yr5.Lag363b.thickness]),3) ' +/- ' num2str(nanstd([Yr5.Lag363a.thickness Yr5.Lag363b.thickness]),2)])
xlabel('Boundary layer thickness, \mum')

subplot(223)
histogram(Yr6.Lag525.thickness,[0:25:200])
title(['Glider 525, Year 6; ' num2str(nanmean(Yr6.Lag525.thickness),3) ' +/- ' num2str(nanstd(Yr6.Lag525.thickness),2)])
xlabel('Boundary layer thickness, \mum')

subplot(224)
histogram([Yr6.Lag560a.thickness Yr6.Lag560b.thickness],[0:25:200])
title('Glider 560, Year 6')
title(['Glider 560, Year 6; ' num2str(nanmean([Yr6.Lag560a.thickness Yr6.Lag560b.thickness]),3) ' +/- ' num2str(nanstd([Yr6.Lag560a.thickness Yr6.Lag560b.thickness]),2)])
xlabel('Boundary layer thickness, \mum')
%%
    figure(11); clf
secinday = 60*60*24;
subplot(221)
    Yr5.Lag453.vert_velocity = (diff(Yr5.Lag453.pres')')./(diff(Yr5.Lag453.mtime')'*secinday);
histogram(abs(Yr5.Lag453.vert_velocity(:)))
title(['Glider 453, Year 5; ' num2str(nanmean(abs(Yr5.Lag453.vert_velocity(:))),3) ' +/- ' num2str(nanstd(abs(Yr5.Lag453.vert_velocity(:))),2) 'dbar s^{-1}'])
xlabel('Vertical velocity, dbar s^{-1}')

subplot(222)
    Yr5.Lag363a.vert_velocity = (diff(Yr5.Lag363a.pres')')./(diff(Yr5.Lag363a.mtime')'*secinday);
    Yr5.Lag363b.vert_velocity = (diff(Yr5.Lag363b.pres')')./(diff(Yr5.Lag363b.mtime')'*secinday);
histogram([abs(Yr5.Lag363a.vert_velocity(:)); abs(Yr5.Lag363b.vert_velocity(:))])
title(['Glider 363, Year 5; ' num2str(nanmean([abs(Yr5.Lag363a.vert_velocity(:)); abs(Yr5.Lag363b.vert_velocity(:))]),3) ' +/- ' num2str(nanstd([abs(Yr5.Lag363a.vert_velocity(:)); abs(Yr5.Lag363b.vert_velocity(:))]),2) 'dbar s^{-1}'])
xlabel('Vertical velocity, dbar s^{-1}')

subplot(223)
    Yr6.Lag525.vert_velocity = (diff(Yr6.Lag525.pres')')./(diff(Yr6.Lag525.mtime')'*secinday);
histogram(abs(Yr6.Lag525.vert_velocity(:)))
title(['Glider 525, Year 6; ' num2str(nanmean(abs(Yr6.Lag525.vert_velocity(:))),3) ' +/- ' num2str(nanstd(abs(Yr6.Lag525.vert_velocity(:))),2) 'dbar s^{-1}'])
xlabel('Vertical velocity, dbar s^{-1}')

subplot(224)
    Yr6.Lag560a.vert_velocity = (diff(Yr6.Lag560a.pres')')./(diff(Yr6.Lag560a.mtime')'*secinday);
    Yr6.Lag560b.vert_velocity = (diff(Yr6.Lag560b.pres')')./(diff(Yr6.Lag560b.mtime')'*secinday);
histogram([abs(Yr6.Lag560a.vert_velocity(:)); abs(Yr6.Lag560b.vert_velocity(:))])
title(['Glider 560, Year 6; ' num2str(nanmean([abs(Yr6.Lag560a.vert_velocity(:)); abs(Yr6.Lag560b.vert_velocity(:))]),3) ' +/- ' num2str(nanstd([abs(Yr6.Lag560a.vert_velocity(:)); abs(Yr6.Lag560b.vert_velocity(:))]),2) 'dbar s^{-1}'])
xlabel('Vertical velocity, dbar s^{-1}')

%%
figure(12); clf
thickness_all = [Yr5.Lag453.thickness Yr5.Lag363a.thickness Yr5.Lag363b.thickness Yr6.Lag525.thickness Yr6.Lag560a.thickness Yr6.Lag560b.thickness];
histogram(thickness_all)
xlabel('Boundary layer thickness, \mum')
title(['All paired up & down profiles, Yr 6 & 7 gliders: mean = ' num2str(nanmean(thickness_all),3) ' +/- ' num2str(nanstd(thickness_all),2)])


%% Apply tau correction to glider data that has been reprocessed for lag correction

tau_in = nanmean(thickness_all);

[Yr5.Lag453.doxy_lagcorr] = glider_lagCorrectFun(Yr5.Lag453, tau_in);
[Yr5.Lag363a.doxy_lagcorr] = glider_lagCorrectFun(Yr5.Lag363a, tau_in);
[Yr5.Lag363b.doxy_lagcorr] = glider_lagCorrectFun(Yr5.Lag363b, tau_in);
[Yr6.Lag525.doxy_lagcorr] = glider_lagCorrectFun(Yr6.Lag525, tau_in);
[Yr6.Lag560a.doxy_lagcorr] = glider_lagCorrectFun(Yr6.Lag560a, tau_in);
[Yr6.Lag560b.doxy_lagcorr] = glider_lagCorrectFun(Yr6.Lag560b, tau_in);


%% Grid lag corrected (and uncorrected) data for plotting

    pmin = 0; pmax = 1000; pinterval = 5;
[Yr5.Lag453] = gliderGrid(Yr5.Lag453, pmin, pmax, pinterval);
[Yr5.Lag363a] = gliderGrid(Yr5.Lag363a, pmin, pmax, pinterval);
[Yr5.Lag363b] = gliderGrid(Yr5.Lag363b, pmin, pmax, pinterval);
[Yr6.Lag525] = gliderGrid(Yr6.Lag525, pmin, pmax, pinterval);
[Yr6.Lag560a] = gliderGrid(Yr6.Lag560a, pmin, pmax, pinterval);
[Yr6.Lag560b] = gliderGrid(Yr6.Lag560b, pmin, pmax, pinterval);

%% Plot output data

downC = nicecolor('bbbc');
down = nicecolor('bbcww');
upC = nicecolor('gby');
up = nicecolor('gbyww');
L1 = 0.5;
L2 = 2;
L3 = 3;
pplotmax = 1000;
titlestr = {'Glider 453, Year 5, June 12-13, 2018', 'Glider 363, Year 5, June 9-10, 2018', 'Glider 363, Year 5, June 12-13, 2018',...
    'Glider 525, Year 6, Aug. 6-7, 2019', 'Glider 560, Year 6, Aug. 7-10, 2019', 'Glider 560, Year 6, Aug. 11-16, 2019'};

figure(1); clf

for i = 1:6
    if i == 1
        Gin = Yr5.Lag453; d1 = 1;
    elseif i == 2
        Gin = Yr5.Lag363a; d1 = 1;
    elseif i == 3
        Gin = Yr5.Lag363b; d1 = 2;
    elseif i == 4
        Gin = Yr6.Lag525; d1 = 1;
    elseif i == 5
        Gin = Yr6.Lag560a; d1 = 2;
    elseif i == 6
        Gin = Yr6.Lag560b; d1 = 2;
    end

subplot(3,4,1+2*(i-1))
plot(Gin.doxy_gridmean(d1:2:end,:),Gin.pgrid,'-','color',down,'linewidth',L1); hold on;
plot(Gin.doxy_lagcorr_gridmean(d1:2:end,:),Gin.pgrid,'-','color',downC,'linewidth',L1); hold on;
plot(Gin.doxy_gridmean(d1+1:2:end,:),Gin.pgrid,'-','color',up,'linewidth',L1); hold on;
plot(Gin.doxy_lagcorr_gridmean(d1+1:2:end,:),Gin.pgrid,'-','color',upC,'linewidth',L1); hold on;

h1 = plot(nanmean(Gin.doxy_gridmean(d1:2:end,:)),Gin.pgrid,'-','color',down,'linewidth',L2); hold on;
h2 = plot(nanmean(Gin.doxy_lagcorr_gridmean(d1:2:end,:)),Gin.pgrid,'-','color',downC,'linewidth',L3); hold on;
h3 = plot(nanmean(Gin.doxy_gridmean(d1+1:2:end,:)),Gin.pgrid,'-','color',up,'linewidth',L2); hold on;
h4 = plot(nanmean(Gin.doxy_lagcorr_gridmean(d1+1:2:end,:)),Gin.pgrid,'-','color',upC,'linewidth',L3); hold on;

axis ij; ylim([-10 pplotmax])
xlabel('Oxygen % (P-corr only)')
ylabel('Pressure (db)')
title(titlestr(i))
legend([h1 h2 h3 h4], 'Dive','Dive, lag corr','Climb','Climb, lag corr','location','SE')

subplot(3,4,2+2*(i-1))
plot([0 0],[-10 pplotmax],'k--'); hold on;
h1 = plot(nanmean(Gin.doxy_gridmean(d1:2:end,:)) - nanmean(Gin.doxy_gridmean(d1+1:2:end,:)),Gin.pgrid,'-','color',nicecolor('kw'),'linewidth',L2); hold on;
h2 = plot(nanmean(Gin.doxy_lagcorr_gridmean(d1:2:end,:)) - nanmean(Gin.doxy_lagcorr_gridmean(d1+1:2:end,:)),Gin.pgrid,'-','color',nicecolor('kkkw'),'linewidth',L3); hold on;

axis ij; ylim([-10 pplotmax])
xlabel('Dive - Climb, Oxygen %')
ylabel('Pressure (db)')
xlim([-8 8])
legend([h1 h2],'No lag corr','After lag corr','location','SE')

end



%% Next steps

% 2) Identify cross-calibration opportunities
    % a) WFP aligned profiles
    % b) Turn-around cruise casts (Winkler + SBE43 on CTD)
    % c) Intercalibration with FLMA and FLMB

% 
% %% Compare with wire-following profiler
% wfp_analysis %Load wfp data and do initial processings
%     dist_compare = 4; %Select radius around mooring to use for comparison (kilometers)
%     time_compare = 1; %Select time difference between mooring and glider profiles to keep for comparison (days)
%     % Need to select profile direction (-1 == up 1 == down); -1 for 363 and 1 for 453
% [G363_profile_summary, wfp_profile_summary_363] = glider_wfp_compare(G363, Yr5_wfp, dist_compare, time_compare, -1);
% [G453_profile_summary, wfp_profile_summary_453] = glider_wfp_compare(G453, Yr5_wfp, dist_compare, time_compare, 1);

