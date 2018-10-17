
filename = 'deployment0005_GI01SUMO-SBD11-06-METBKA000-telemetered-metbk_a_dcl_instrument_20180608T172109.234000-20180920T090824.973000.nc';

% constants
mbar2atm = 1013.25;
sec2day = 60*60*24;

% options
load('latest');
% rhcorr = 1 to use observed relative humidy rhcorr = 0 to assume saturated
% water vapor
rhcorr = 1;
% profile direction (-1 == up 1 == down)
% choose which glider to process by commenting/uncommenting here
% G = G363;
% prof_dir = -1;

G = G453;
prof_dir = 1;

% boolean flag for profile direction
isup = prof_dir == -1;
%%

%Pull out desired met variables
% time is in sec since 190000
met.daten = datenum(1900,1,1,00,0,ncread(filename,'time'));
met.barometric_pressure = ncread(filename,'barometric_pressure'); %mbar
met.patm = met.barometric_pressure / mbar2atm;
met.rh = ncread(filename,'relative_humidity');
met.water_temp = ncread(filename,'sea_surface_temperature');
met.salinity = 34.8;%ncread(filename,'met_salsurf');

%% filter for surface measurements in upper 10m
% saturated water vapor pressure
met.SVP_S = vpress(met.salinity,met.water_temp);
if rhcorr == 1
    met.ph2o = met.rh.*met.SVP_S./100;
    met.O2satcorr = 100.*(met.patm - met.ph2o)./(1 - met.SVP_S);

else
    met.O2satcorr = 100.*(met.patm - met.SVP_S)./(1 - met.SVP_S);
end


%%
qntl = 0.32;
%q
GAIN = 1.00;
G.oxygen_concentration = GAIN.*(G.oxygen_concentration);
G.oxygen_saturation = GAIN.*(G.oxygen_saturation);

np = max(floor(G.profile_index));
vars = {'air_daten','air_meas','air_corr','met_o2sat','ml_daten','ml_o2sat','ml_tem','ml_sal','nsurf'};
T = array2table(nan(np,length(vars)));
T.Properties.VariableNames = vars;
% buoy correction

tsurf = [];
tsurf2 = [];
airmeas =[];
profid = [];

for ii = 1:np

    %near-surface measurements from glider profile
    u = G.profile_index == ii & G.profile_direction == prof_dir & G.depth_interp < 10 & G.depth_interp > 0.5 & G.oxygen_saturation > 20;
    uts = G.profile_index == ii & G.profile_direction == prof_dir & G.depth_interp < 10 & G.depth_interp > 0.5;
    % select surface interval
    s = G.profile_index == ii-0.5+isup & ~isnan(G.oxygen_saturation) & G.depth_interp < 0.5;


    t0 = find(s > 0,1);
    if ~isempty(t0)
        % time window for surface data: 90 sec < obs < 800 sec
        s2 = s & G.daten - G.daten(t0) > 90/(24*60*60) & G.daten - G.daten(t0) < 800/(24*60*60);


        % accumulates all individual air measurments
        tsurf = [tsurf; 24*60*60*(G.daten(s)-G.daten(t0))];
        airmeas = [airmeas; G.oxygen_saturation(s)-naninterp1(met.daten,100.*(met.patm),G.daten(s))];
        profid = [profid; 0*G.oxygen_saturation(s)+ii-0.5+isup];
        O2air = G.oxygen_saturation(s2);
        T.nsurf(ii) = sum(s2);
    else
        O2air = nan;
    end


    T.ml_tem(ii) = nanmean(G.temperature(uts));
    T.ml_sal(ii) = nanmean(G.salinity(uts));
    T.ml_o2sat(ii) = nanmedian(G.oxygen_saturation(u));
    T.ml_daten(ii) = nanmean(G.daten(u));
    T.air_meas(ii) = quantile(O2air,qntl);
    T.air_daten(ii,1) = nanmean(G.daten(s));

end
T.met_o2sat = naninterp1(met.daten,met.O2satcorr,T.air_daten);
d = ~isnan(T.air_meas+T.ml_o2sat) & T.air_meas > 0;
T(~d,:) = [];
p = polyfit(T.ml_o2sat,T.air_meas,1);
T.air_corr = (T.air_meas-p(1).*T.ml_o2sat)./(1-p(1));

%% Figures

ftsz = 14;
lnw = 1;
med_gain = median(T.met_o2sat(~isnan(T.met_o2sat))./T.air_corr(~isnan(T.met_o2sat)));
%med_gain = 1.8;
figure;
ax1 = gca;
hold all;
box on;
cols = ax1.ColorOrder;
ax1.LineWidth = 1;
ax1.FontSize = ftsz;
title(['median gain: ', num2str(med_gain)]);
plot(T.ml_daten,T.ml_o2sat,'.-','MarkerSize',12,'LineWidth',2,'Color',cols(1,:));
hold on;
plot(T.air_daten,T.air_meas,'.-','MarkerSize',12,'LineWidth',2,'Color',cols(3,:));
plot(met.daten,met.O2satcorr,'-','LineWidth',4,'Color',cols(2,:));

plot(T.air_daten,T.air_corr,'.-','MarkerSize',12,'LineWidth',2,'Color',cols(4,:));
plot(T.air_daten,med_gain.*T.air_corr,'.-','MarkerSize',12,'LineWidth',2,'Color',cols(5,:));
datetick;
legend('\DeltaO_{2,w}^{meas}','\DeltaO_{2,a}^{meas}','\DeltaO_{2,w}^{met}','\DeltaO_{2,w}^{mcorr}','\DeltaO_{2,w}^{gaincorr}');
%%
figure;
histogram(T.air_corr-T.met_o2sat,-5:.2:10);

figure;
plot(T.ml_o2sat,T.air_meas,'.','MarkerSize',12);
hold all;
box on;
ax3 = gca;
plot(ax3.XLim,p(1).*ax3.XLim+p(2),'-k','LineWidth',3);
ax3.LineWidth = 1;
ax3.FontSize = ftsz;
xlabel('\DeltaO_{2,w}^{meas}');
ylabel('\DeltaO_{2,a}^{meas}');

figure;
    %new_points = find(T.air_daten > datenum(2018,9,4));
plot(T.met_o2sat,T.air_corr,'.','MarkerSize',12); hold on;
%plot(T.met_o2sat(new_points),T.air_corr(new_points),'m.','MarkerSize',18); hold on;
    ind = find(isnan(T.met_o2sat + T.air_corr) == 0);
    p2 = polyfit(T.met_o2sat(ind),T.air_corr(ind),1);
hold all;
box on;
ax3 = gca;
plot(ax3.XLim,p2(1).*ax3.XLim+p2(2),'-k','LineWidth',3);
ax3.LineWidth = 1;
ax3.FontSize = ftsz;
xlabel('\DeltaO_{2,w}^{met}');
ylabel('\DeltaO_{2,w}^{mcorr}');

[rho,df,rho_sig95] = correlate(T.met_o2sat(ind),T.air_corr(ind))
