%% Lag correction - only for years 5 and 6 (no paired profiles for lag assessment in years 7 and 8)
% Correct for historesis effects between up and down profiles using Gordon et al. 2020 approach, as in wfp_lag
addpath(genpath('C:\Users\Palevsky\Documents\GitHub\optode-response-time'))

% only one set of paired dive and climb, aligns with OSM 2020 analysis and 363 set b - NO TEMP DATA FOR DIVE!
[Yr5.Lag453] = glider_lagAssessFun(Yr5.G453, 26000, 39400, 300, 'Glider 453, Year 5', [43:51]);

% two sets of paired dive and climb, second aligns with cal cast analyzed for OSM 2020 lag analysis
[Yr5.Lag363a] = glider_lagAssessFun(Yr5.G363, 10000, 210000, 300, 'Glider 363, Year 5', [13:20]); %just first set
[Yr5.Lag363b] = glider_lagAssessFun(Yr5.G363, 10000, 210000, 300, 'Glider 363, Year 5', [42:46]); %just second set - NO TEMP DATA FOR DIVE!

%only one set of paired dive and climb
[Yr6.Lag525] = glider_lagAssessFun(Yr6.G525, 800, 68200, 300, 'Glider 525, Year 6', [9:16]);

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
histogram([Yr5.Lag363a.thickness],[0:25:200])
title(['Glider 363, Year 5; ' num2str(nanmean([Yr5.Lag363a.thickness]),3) ' +/- ' num2str(nanstd([Yr5.Lag363a.thickness]),2)])
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
thickness_all = [Yr5.Lag363a.thickness Yr6.Lag525.thickness Yr6.Lag560a.thickness Yr6.Lag560b.thickness]; %Yr5.Lag453.thickness Yr5.Lag363b.thickness 
histogram(thickness_all,[0:10:150])
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
titlestr = {'Glider 363, Year 5, June 9-10, 2018', 'Glider 525, Year 6, Aug. 6-7, 2019',...
    'Glider 560, Year 6, Aug. 7-10, 2019', 'Glider 560, Year 6, Aug. 11-16, 2019'};

figure(1); clf

for i = 1:4
    if i == 1
        %Gin = Yr5.Lag453; d1 = 1;
        %Gin = Yr5.Lag363b; d1 = 2;
        Gin = Yr5.Lag363a; d1 = 1;
    elseif i == 2
        Gin = Yr6.Lag525; d1 = 1;
    elseif i == 3
        Gin = Yr6.Lag560a; d1 = 2;
    elseif i == 4
        Gin = Yr6.Lag560b; d1 = 2;
    end

subplot(2,4,1+2*(i-1))
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

subplot(2,4,2+2*(i-1))
plot([0 0],[-10 pplotmax],'k--'); hold on;
h1 = plot(nanmean(Gin.doxy_gridmean(d1:2:end,:)) - nanmean(Gin.doxy_gridmean(d1+1:2:end,:)),Gin.pgrid,'-','color',nicecolor('kw'),'linewidth',L2); hold on;
h2 = plot(nanmean(Gin.doxy_lagcorr_gridmean(d1:2:end,:)) - nanmean(Gin.doxy_lagcorr_gridmean(d1+1:2:end,:)),Gin.pgrid,'-','color',nicecolor('kkkw'),'linewidth',L3); hold on;

axis ij; ylim([-10 pplotmax])
xlabel('Dive - Climb, Oxygen %')
ylabel('Pressure (db)')
xlim([-6 6])
legend([h1 h2],'No lag corr','After lag corr','location','SE')

end

%% Diagnostic plot using temperature, showing match for used data and missing data for unused data from June 12-13, 2018
figure(2); clf
% titlestr = {'Glider 453, Year 5, June 12-13, 2018', 'Glider 363, Year 5, June 9-10, 2018', 'Glider 363, Year 5, June 12-13, 2018',...
%     'Glider 525, Year 6, Aug. 6-7, 2019', 'Glider 560, Year 6, Aug. 7-10, 2019', 'Glider 560, Year 6, Aug. 11-16, 2019'};

for i = 1:4
    if i == 1
        %Gin = Yr5.Lag453; d1 = 1;
        %Gin = Yr5.Lag363b; d1 = 2;
        Gin = Yr5.Lag363a; d1 = 1;
    elseif i == 2
        Gin = Yr6.Lag525; d1 = 1;
    elseif i == 3
        Gin = Yr6.Lag560a; d1 = 2;
    elseif i == 4
        Gin = Yr6.Lag560b; d1 = 2;
    end

subplot(2,4,1+2*(i-1))
plot(Gin.temp_gridmean(d1:2:end,:),Gin.pgrid,'-','color',down,'linewidth',L1); hold on;
plot(Gin.temp_gridmean(d1+1:2:end,:),Gin.pgrid,'-','color',up,'linewidth',L1); hold on;

h1 = plot(nanmean(Gin.temp_gridmean(d1:2:end,:)),Gin.pgrid,'-','color',down,'linewidth',L2); hold on;
h2 = plot(nanmean(Gin.temp_gridmean(d1+1:2:end,:)),Gin.pgrid,'-','color',up,'linewidth',L2); hold on;

axis ij; ylim([-10 pplotmax])
xlabel('Temp (^oC)')
ylabel('Pressure (db)')
title(titlestr(i))
legend([h1 h2], 'Dive','Climb','location','SE')

subplot(2,4,2+2*(i-1))
plot([0 0],[-10 pplotmax],'k--'); hold on;
h1 = plot(nanmean(Gin.temp_gridmean(d1:2:end,:)) - nanmean(Gin.temp_gridmean(d1+1:2:end,:)),Gin.pgrid,'-','color',nicecolor('kw'),'linewidth',L2); hold on;

axis ij; ylim([-10 pplotmax])
xlabel('Dive - Climb, Temp (^oC)')
ylabel('Pressure (db)')
xlim([-0.5 0.5])

end
