function [T, met] = aircalfun(G, glg, prof_dir, filename, rhcorr, smoothval, outtolsd)

%===================================================================
%
% DESCRIPTION:
%   Use this function to analyze glider oxygen air calibration measurements
%   in combination with meteorological data to calculate gain corrections
%   and create diagnostic plots
%
% INPUT:
%    G: Table of input glider data before applying any corrections
%    glg: Table of input glider data, after having applied S, P, and lag
%       corrections (using glider_interpCorrFun, glider_reshape, and glider_lagCorrectFun)
%    prof_dir: profile direction (-1 == up 1 == down)
%    filename: string with name of netcdf file with meterological data
%    rhcorr: rhcorr = 1 to use observed relative humidy rhcorr = 0 to assume saturated water vapor
%    smoothval: number of points over which to apply moving mean before removing outliers
%    tolsd: number of standard deviations outside moving mean to use for outlier flagging
%
% OUTPUT:
%    T: Table with data from all air calibration intervals over the course of the deployment
%    met: Table with data from associated meteorological measurements
%
% AUTHOR:   Hilary Palevsky, updated 24 April 2023
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

%% Remove bad barometric pressure data
tol_low = 920; %identified based on metbka_assess plot
idout = find(met.barometric_pressure < tol_low);
met.barometric_pressure(idout) = NaN;
met.patm(idout) = NaN;

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
qntl = 0.5; %(value from Nicholson and Feen, 2017, p. 499 is 0.32 but using median here)
percentiles = [0.05:0.05:0.95];
    %Depths used to define near-surface oxygen measurements
surf_mindepth = 0.5;
surf_maxdepth = 10;
O2min = 60;
O2satmin_air = 80;
    %Time window for surface data: twin_beg sec < obs < twin_end sec
twin_beg = 90;
twin_end = 800;

%Create a table to hold output
np = max(floor(G.profile_index)); %number of profiles
vars = {'air_daten','air_meas','air_corr','met_o2sat','ml_daten','ml_o2conc','ml_o2conc_nocorr','ml_o2sat','ml_tem','ml_sal','nsurf'};
T = array2table(nan(np,length(vars)));
T.Properties.VariableNames = vars;

%Loop over all profiles
for ii = 1:np   
    % select near-surface measurements from glider profile (from corrected data)
    glgind = find((glg.profilelist == ii) == 1);
    if glg.profile_direction(glgind) == prof_dir
        u = glg.depth(glgind,:) < surf_maxdepth & glg.depth(glgind,:) > surf_mindepth & glg.doxy_lagcorr(glgind,:) > O2min;
        uts = glg.depth(glgind,:) < surf_maxdepth & glg.depth(glgind,:) > surf_mindepth;
    else
        u = NaN;
        uts = NaN;
    end
    
    % select surface interval (from raw data table)
    s = G.profile_index == ii-0.5+isup & ~isnan(G.oxygen_saturation) & G.depth_interp < surf_mindepth & G.oxygen_saturation > O2satmin_air;

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
    if sum(~isnan(uts)) > 0
        T.ml_tem(ii) = nanmean(glg.temp(glgind,uts)); %from corrected data
        T.ml_sal(ii) = nanmean(glg.sal(glgind,uts)); %from corrected data
        T.ml_o2conc(ii) = nanmedian(glg.doxy_lagcorr(glgind,u)); %from corrected data
        T.ml_o2conc_nocorr(ii) = nanmedian(glg.doxy(glgind,u)); %from corrected data
        T.ml_daten(ii) = nanmean(glg.mtime(glgind,u)); %from corrected data
    end
    T.air_meas(ii) = quantile(O2air,qntl); %from raw data table
    T.air_daten(ii,1) = nanmean(G.daten(s)); %from raw data table
    T.air_meas_dist(ii,:) = quantile(O2air,percentiles); %from raw data table
end

%% Remove lines with missing data
d = ~isnan(T.air_meas+T.ml_o2conc) & T.air_meas > 0;
T(~d,:) = [];

%% Flag and remove outliers
ind_mlout = find(abs(T.ml_o2conc - movmean(T.ml_o2conc,smoothval,'omitnan','endpoints','fill')) > outtolsd*nanstd((T.ml_o2conc)));
    T.ml_o2conc(ind_mlout) = NaN;
    T.ml_o2conc_nocorr(ind_mlout) = NaN;

T.ml_tem(T.ml_tem < 0.1) = NaN; %removing data approaching zero
ind_temout = find(abs(T.ml_tem - movmean(T.ml_tem ,smoothval,'omitnan','endpoints','fill')) > outtolsd*nanstd((T.ml_tem)));
    T.ml_tem(ind_temout) = NaN;

T.ml_sal(T.ml_sal < 0.1) = NaN; %removing data approaching zero
ind_salout = find(abs(T.ml_sal - movmean(T.ml_sal ,smoothval,'omitnan','endpoints','fill')) > outtolsd*nanstd((T.ml_sal)));
    T.ml_sal(ind_salout) = NaN;
    
ind_airout = find(abs(T.air_meas - movmean(T.air_meas,smoothval,'omitnan','endpoints','fill')) > outtolsd*nanstd((T.air_meas)));
    T.air_meas(ind_airout) = NaN;

ind_metout = find(abs(met.O2satcorr - movmean(met.O2satcorr,smoothval,'omitnan','endpoints','fill')) > outtolsd*nanstd((met.O2satcorr)));
    met.O2satcorr(ind_metout) = NaN;

%% Determine remaining output values
T.ml_o2sat = T.ml_o2conc./gsw_O2sol_SP_pt(T.ml_sal, T.ml_tem)*100;
T.met_o2sat = naninterp1(met.daten,met.O2satcorr,T.air_daten);

%% Again remove lines with missing data
d = ~isnan(T.air_meas+T.ml_o2sat) & T.air_meas > 0;
T(~d,:) = [];

%% Correct air measurements for surface water splashing
p = polyfit(T.ml_o2sat,T.air_meas,1);
T.air_corr = (T.air_meas-p(1).*T.ml_o2sat)./(1-p(1)); %eqn 5 in Nicholson and Feen 2017

end