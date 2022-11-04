
%% Analysis of glider air calibration test data from the Pioneer array
addpath('G:\Shared drives\NSF_Irminger\Data_Files\From_Roo\Pioneer') %Processed glider data
load G559.mat %data from Pioneer test of new air cal design
G559 = T;

%Need met data to do air calibration - downloaded data from offshore surface mooring - CP04OSSM (central surface mooring seems to have failed on 9/13)
%Downloaded via OOINet on 10/31/22
addpath('G:\Shared drives\NSF_Irminger\Data_Files\METBKA\Pioneer') %MET data
P15.met = 'deployment0015_CP04OSSM-SBD11-06-METBKA000-telemetered-metbk_a_dcl_instrument_20220801T080515.187000-20221031T140452.149000.nc';

%% Interpolate GPS and CTD data onto all time points
G559.lon_interp = naninterp1(G559.time, G559.longitude, G559.time);
G559.lat_interp = naninterp1(G559.time, G559.latitude, G559.time);
G559.salinity_interp = naninterp1(G559.time, G559.salinity, G559.time);
G559.pressure_interp = naninterp1(G559.time, G559.pressure, G559.time);
G559.temperature_interp = naninterp1(G559.time, G559.temperature, G559.time);

%% Note - haven't applied a lag correction
%Correcting assuming that internal salinity setting is 0

G559.O2_corr = aaoptode_salpresscorr(G559.oxygen_concentration, G559.temperature_interp, G559.salinity_interp, G559.pressure_interp, 0);
G559.O2sat_corr = G559.oxygen_saturation.*(1+G559.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

%% Plot full depth data over entire deployments
plotstart = datenum(2022,8,21,0,0,0);
plotend = datenum(2022,10,4,0,0,0);

figure(1); clf
    subplot(211)
plot(G559.daten, G559.depth_interp, 'k.'); hold on;
scatter(G559.daten, G559.depth_interp, [], G559.O2sat_corr,'filled'); colorbar; caxis([48 122])
set(gca,'YDir','reverse'); 
xlim([plotstart plotend])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 559, oxygen saturation')
    subplot(212)
plot(G559.daten, G559.depth_interp, 'k.'); hold on;
scatter(G559.daten, G559.depth_interp, [], G559.temperature,'filled'); colorbar; caxis([4 28])
set(gca,'YDir','reverse'); 
xlim([plotstart plotend])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Glider 559, temperature')

%% Calculate statistics for each air interval and keep only those with low stdev
    spacer = 10; %number of air measurements at beginning of each surfacing to cut
    tol = 1; %only keep air data if standard deviation of oxygen saturation measurements after cutting spacer is less than this value
[stats_559] = glideraircal_stats(G559, spacer, tol);

figure(2); clf
    d = rem(G559.profile_index,1) == 0.5;
histogram(G559.oxygen_saturation(d),[99:0.2:110]); hold on;
title('Histogram of all raw G559 oxygen saturation data during air intervals')

figure(3); clf
    d = rem(G559.profile_index,1) == 0.5;
plot(G559.daten(d), G559.oxygen_saturation(d),'.', 'color', 'b'); hold on;
    ind = find(stats_559(:,7) == 1);
plot(stats_559(ind,6), stats_559(ind,3), 'ro','markerfacecolor','r'); hold on;
ylim([97 115])
xlim([plotstart plotend])
datetick('x',2,'keeplimits')
legend('Glider 559, all','Glider 559, "good" interval mean','location','northeast')
ylabel('Raw oxygen saturation (%)')
ylim([97 106])
title('Oxygen measurements during glider surface intervals measuring air, Pioneer Array test of new mount')

%% Draw from aircalfun but do in script to enable modifications

prof_dir = 1; %sampling on dive (downwards)
mindateplot = datenum(2022,8,20);
maxdateplot = datenum(2022,10,5);
rhcorr = 1; 
filename = P15.met;
G = G559;
tref = mindateplot;

%% Set values used throughout
%Constants
mbar2atm = 1013.25;
sec2day = 60*60*24;

% boolean flag for profile direction
isup = prof_dir == -1;

%% Pull out desired meteorological variables from input SUMO file
met.daten = datenum(1900,1,1,00,0,ncread(filename,'time')); % time is in sec since 190000
met.barometric_pressure = ncread(filename,'barometric_pressure'); %mbar
met.patm = met.barometric_pressure / mbar2atm;
met.rh = ncread(filename,'relative_humidity');
met.water_temp = ncread(filename,'sea_surface_temperature');
met.salinity = ncread(filename,'met_salsurf');

%% Calculate O2 measurement expected in air based on SUMO data

met.SVP_S = vpress(met.salinity,met.water_temp); % saturated water vapor pressure
if rhcorr == 1
    %use observed water vapor pressure
    met.ph2o = met.rh.*met.SVP_S./100;
    met.O2satcorr = 100.*(met.patm - met.ph2o)./(1 - met.SVP_S);
else
    %assume saturated water vapor pressure
    met.O2satcorr = 100.*(met.patm - met.SVP_S)./(1 - met.SVP_S);
end

%% Filter for surface measurements in upper 10m

%%%% Constants used in calculations below
    %Percentile of air measurement distribution to use
qntl = 0.32; %(value from Nicholson and Feen, 2017, p. 499)
    %Depths used to define near-surface oxygen measurements
surf_mindepth = 0.5;
surf_maxdepth = 10;
O2satmin = 20;
    %Time window for surface data: twin_beg sec < obs < twin_end sec
twin_beg = 90;
twin_end = 800;

%Create a table to hold output
np = max(floor(G.profile_index)); %number of profiles
vars = {'air_daten','air_meas','air_corr','met_o2sat','ml_daten','ml_o2sat','ml_tem','ml_sal','nsurf'};
T = array2table(nan(np,length(vars)));
T.Properties.VariableNames = vars;

%Loop over all profiles
for ii = 1:np   
    % select near-surface measurements from glider profile
    u = G.profile_index == ii & G.profile_direction == prof_dir & G.depth_interp < surf_maxdepth & G.depth_interp > surf_mindepth & G.oxygen_saturation > O2satmin;
    uts = G.profile_index == ii & G.profile_direction == prof_dir & G.depth_interp < surf_maxdepth & G.depth_interp > surf_mindepth;
    
    % select surface interval
    s = G.profile_index == ii-0.5+isup & ~isnan(G.oxygen_saturation) & G.depth_interp < surf_mindepth & G.oxygen_saturation > 80;

    % extract air oxygen measurements from surface interval
    t0 = find(s > 0,1);
    if ~isempty(t0)
        % time window for surface data: twin_beg sec < obs < twin_end sec
        s2 = s & G.daten - G.daten(t0) > twin_beg/sec2day & G.daten - G.daten(t0) < twin_end/sec2day;
        O2air = G.oxygen_saturation(s2);
        T.nsurf(ii) = sum(s2);
    else
        O2air = nan;
    end

    %Save output value for each profile
    T.ml_tem(ii) = nanmean(G.temperature(uts));
    T.ml_sal(ii) = nanmean(G.salinity(uts));
    T.ml_o2sat(ii) = nanmedian(G.oxygen_saturation(u));
    T.ml_daten(ii) = nanmean(G.daten(u));
    T.air_meas(ii) = quantile(O2air,qntl);
    T.air_daten(ii,1) = nanmean(G.daten(s));
end

%% Determine remaining output values
T.met_o2sat = naninterp1(met.daten,met.O2satcorr,T.air_daten);

%Remove lines with missing data
d = ~isnan(T.air_meas+T.ml_o2sat) & T.air_meas > 0;
T(~d,:) = [];

%Correct air measurements for surface water splashing
p = polyfit(T.ml_o2sat,T.air_meas,1);
T.air_corr = (T.air_meas-p(1).*T.ml_o2sat)./(1-p(1));

%% Calculate gain corrections

%This is a version not accounting for surface water splashing
px = polyfit(T.ml_daten-tref,T.met_o2sat./T.air_meas,2);
T.air_meas_corr = T.air_meas.*polyval(px,T.ml_daten-tref);
T.ml_o2sat_corr = T.ml_o2sat.*polyval(px,T.ml_daten-tref);

%Account for surface water splashing
px2 = polyfit(T.ml_daten-tref,T.met_o2sat./T.air_corr,1);
med_gain = median(T.met_o2sat(~isnan(T.met_o2sat))./T.air_corr(~isnan(T.met_o2sat)));

%% Figures
    %set figure options
ftsz = 14;
lnw = 1.5;
mrkr = 10;
gliderstring = 'Pioneer G559';

%Plot time series of gain data with and without splash corrections
figure(4); clf
    ax1 = gca;
    hold all;
    ax1.FontSize = ftsz;
    dateplot = [floor(min(T.ml_daten)):ceil(max(T.ml_daten))];
plot(T.ml_daten, T.met_o2sat./T.air_meas, 'k.'); %not accounting for surface water splashing
    plot(dateplot, polyval(px,dateplot-tref),'k-','linewidth',lnw);
plot(T.ml_daten, T.met_o2sat./T.air_corr, 'b.'); %account for surface water splashing
    plot(dateplot, polyval(px2,dateplot-tref),'b-','linewidth',lnw);
    plot(dateplot, med_gain*ones(size(dateplot)),'b--','linewidth',lnw);
datetick('x',2,'keeplimits')
ylabel('Gain correction')
legend('Gain data w/ no splash correction','Gain quadratic fit w/ no splash correction',...
    'Gain data w/ splash correction','Gain linear fit w/ splash correction','Median gain w/ splash correction')
title([gliderstring ' gain corrections'])

%% Plot time series of glider data for gain calculations and corresponding MET data
figure(5); clf
    ax2 = gca;
    hold all;
    box on;
    cols = ax2.ColorOrder;
    ax2.FontSize = ftsz;
plot(T.ml_daten,T.ml_o2sat,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(1,:));
plot(T.air_daten,T.air_meas,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(3,:));
%plot(T.air_daten,T.air_corr,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(4,:));
plot(T.air_daten,med_gain.*T.air_corr,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(5,:));
plot(met.daten,met.O2satcorr,'.','LineWidth',lnw,'Color','k');
xlim([mindateplot maxdateplot])
title([gliderstring ' air calibration: gain = ' num2str(med_gain)]);
ylabel('Oxygen saturation (%)')
ylim([98 106])
datetick('x',2,'keeplimits');
legend('\DeltaO_{2,w}^{meas}','\DeltaO_{2,a}^{meas}','\DeltaO_{2,a}^{splash & gain corr}','\DeltaO_{2}^{met}');

%% Plot relationship between glider surface water and air measurements
%(evidence of splashing)
figure(6);
    ax3 = gca;
    ax3.FontSize = ftsz;
    hold all;
    box on;
plot(T.ml_o2sat,T.air_meas,'.','MarkerSize',mrkr);
plot(ax3.XLim,p(1).*ax3.XLim+p(2),'-k','LineWidth',lnw);
xlabel('\DeltaO_{2,w}^{meas}');
ylabel('\DeltaO_{2,a}^{meas}');
    LM = fitlm(array2table([T.ml_o2sat,T.air_meas]),'linear');
title([gliderstring ' air vs. surface water, R^2 = ' num2str(LM.Rsquared.Adjusted,3)]);

%% Relationship between MET and corrected glider air data
figure(7);
    ax4 = gca;
    ax4.FontSize = ftsz;
    hold all;
    box on;
plot(T.met_o2sat,T.air_corr,'.','MarkerSize',mrkr); hold on;
    ind = find(isnan(T.met_o2sat + T.air_corr) == 0);
    p2 = polyfit(T.met_o2sat(ind),T.air_corr(ind),1);
plot(ax4.XLim,p2(1).*ax4.XLim+p2(2),'-k','LineWidth',lnw);
xlabel('\DeltaO_{2}^{met}');
ylabel('\DeltaO_{2,a}^{splash & gain corr}');
    LM = fitlm(array2table([T.met_o2sat,T.air_corr]),'linear');
title([gliderstring ' corrected air vs. MET data, R^2 = ' num2str(LM.Rsquared.Adjusted,3)]);