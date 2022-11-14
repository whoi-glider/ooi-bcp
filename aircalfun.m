function [T, med_gain] = aircalfun(G, gliderstring, prof_dir, filename, tref, rhcorr, mindateplot, maxdateplot)

%===================================================================
%
% DESCRIPTION:
%   Use this function to analyze glider oxygen air calibration measurements
%   in combination with meteorological data to calculate gain corrections
%   and create diagnostic plots
%
% INPUT:
%    G: Table of input glider data
%    gliderstring: string with glider name (used for plot titles)
%    prof_dir: profile direction (-1 == up 1 == down)
%    filename: string with name of netcdf file with meterological data
%    tref: time reference in datenum format (generally choose beginning of deployment)
%    rhcorr: rhcorr = 1 to use observed relative humidy rhcorr = 0 to assume saturated water vapor
%    mindateplot: time in datenum format of minimum date for plotting
%    maxdateplot: time in datenum format of maximum date for plotting
%
% OUTPUT:
%    T: Table with data from all air calibration intervals over the course of the deployment
%    med_gain: median gain correction over the course of the deployment
%
% AUTHOR:   Hilary Palevsky, 3 February 2020
%
% REFERENCE:
%    Based on code by D. P. Nicholson and analysis published in:
%    Nicholson, D. P., & Feen, M. L. (2017). Air calibration of an oxygen optode
%    on an underwater glider. Limnol. Oceanogr.: Methods, 15, 495–502. https://doi.org/10.1002/lom3.10177

%==================================================================

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

%Plot time series of gain data with and without splash corrections
figure; clf
subplot(3,2,1:2)
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

%Plot time series of glider data for gain calculations and corresponding MET data
subplot(3,2,3:4)
    ax2 = gca;
    hold all;
    box on;
    cols = ax2.ColorOrder;
    ax2.FontSize = ftsz;
plot(T.ml_daten,T.ml_o2sat,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(1,:));
plot(T.air_daten,T.air_meas,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(3,:));
plot(T.air_daten,T.air_corr,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(4,:));
plot(T.air_daten,med_gain.*T.air_corr,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(5,:));
plot(met.daten,met.O2satcorr,'-','LineWidth',lnw,'Color','k');
xlim([mindateplot maxdateplot])
title([gliderstring ' air calibration: median gain = ' num2str(med_gain)]);
ylabel('Oxygen saturation (%)')
datetick('x',2,'keeplimits');
legend('\DeltaO_{2,w}^{meas}','\DeltaO_{2,a}^{meas}','\DeltaO_{2,a}^{splash corr}','\DeltaO_{2,a}^{splash & gain corr}','\DeltaO_{2}^{met}');

%Plot relationship between glider surface water and air measurements
%(evidence of splashing)
subplot(3,2,5)
    ax3 = gca;
    ax3.FontSize = ftsz;
    hold all;
    box on;
plot(T.ml_o2sat,T.air_meas,'.','MarkerSize',mrkr);
plot(ax3.XLim,p(1).*ax3.XLim+p(2),'-k','LineWidth',lnw);
xlabel('\DeltaO_{2,w}^{meas}');
ylabel('\DeltaO_{2,a}^{meas}');
    LM = fitlm(array2table([T.ml_o2sat,T.air_meas]),'linear');
title([gliderstring ' air vs. surface water, R^2 = ' num2str(LM.Rsquared.Adjusted,3) ', slope = ' num2str(p(1),3)]);

%Relationship between MET and corrected glider air data
subplot(3,2,6)
    ax4 = gca;
    ax4.FontSize = ftsz;
    hold all;
    box on;
plot(T.met_o2sat,T.air_corr,'.','MarkerSize',mrkr); hold on;
    ind = find(isnan(T.met_o2sat + T.air_corr) == 0);
    p2 = polyfit(T.met_o2sat(ind),T.air_corr(ind),1);
plot(ax4.XLim,p2(1).*ax4.XLim+p2(2),'-k','LineWidth',lnw);
xlabel('\DeltaO_{2}^{met}');
ylabel('\DeltaO_{2,a}^{splash corr}');
    LM = fitlm(array2table([T.met_o2sat,T.air_corr]),'linear');
title([gliderstring ' corrected air vs. MET data, R^2 = ' num2str(LM.Rsquared.Adjusted,3)]);

end