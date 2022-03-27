%% Look at the MET data
% Note that there are at least slight differences between telemetered data
% downloaded 3/14/2022 vs earlier download of same data - use new data

%Next step is to check telemetered vs recovered
%In Year 5, telemetered and recovered data are the same

Yr5_met_t = 'deployment0005_GI01SUMO-SBD11-06-METBKA000-telemetered-metbk_a_dcl_instrument_20180608T172109.234000-20190809T080351.721000.nc';
Yr5_met_r = 'deployment0005_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20180608T172109.234000-20190809T080351.721000.nc';

met5t.daten = datenum(1900,1,1,00,0,ncread(Yr5_met_t,'time')); % time is in sec since 190000
met5t.barometric_pressure = ncread(Yr5_met_t,'barometric_pressure'); %mbar
met5r.daten = datenum(1900,1,1,00,0,ncread(Yr5_met_r,'time')); % time is in sec since 190000
met5r.barometric_pressure = ncread(Yr5_met_r,'barometric_pressure'); %mbar

%In Year 6, need to merge the two for full coverage: period from late April
%through most of May only in telemetered, and long period from mid-Oct
%through March only in recovered; however no gliders past early April, so
%can just use recovered
Yr6_met_t = 'deployment0006_GI01SUMO-SBD11-06-METBKA000-telemetered-metbk_a_dcl_instrument_20190805T152914.351000-20200826T101155.752000.nc';
Yr6_met_r = 'deployment0006_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20190805T152914.351000-20200826T104006.497000.nc';

met6t.daten = datenum(1900,1,1,00,0,ncread(Yr6_met_t,'time')); % time is in sec since 190000
met6t.barometric_pressure = ncread(Yr6_met_t,'barometric_pressure'); %mbar
met6r.daten = datenum(1900,1,1,00,0,ncread(Yr6_met_r,'time')); % time is in sec since 190000
met6r.barometric_pressure = ncread(Yr6_met_r,'barometric_pressure'); %mbar

%Year 7 - a bit extra in the recovered, but nothing telemetered not in
%recovered
Yr7_met_t = 'deployment0007_GI01SUMO-SBD11-06-METBKA000-telemetered-metbk_a_dcl_instrument_20200817T173428.501000-20210819T064814.736000.nc';
Yr7_met_r = 'deployment0007_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20200817T173428.501000-20210819T064814.736000.nc';

met7t.daten = datenum(1900,1,1,00,0,ncread(Yr7_met_t,'time')); % time is in sec since 190000
met7t.barometric_pressure = ncread(Yr7_met_t,'barometric_pressure'); %mbar
met7r.daten = datenum(1900,1,1,00,0,ncread(Yr7_met_r,'time')); % time is in sec since 190000
met7r.barometric_pressure = ncread(Yr7_met_r,'barometric_pressure'); %mbar

figure(1); clf

subplot(3,1,1)
plot(met5r.daten, met5r.barometric_pressure, 'b.','markersize',10); hold on;
plot(met5t.daten, met5t.barometric_pressure, 'k.','markersize',5)
datetick('x')
legend('recovered','telemetered','location','southeast')
title('OOI Irminger METBKA Barometric Pressure, Year 5: 2018-2019')

subplot(3,1,2)
plot(met6r.daten, met6r.barometric_pressure, 'b.','markersize',10); hold on;
plot(met6t.daten, met6t.barometric_pressure, 'k.','markersize',5)
datetick('x')
ylabel('mbar')
ylabel('mbar')
legend('recovered','telemetered','location','southeast')
title('OOI Irminger METBKA Barometric Pressure, Year 6: 2019-2020')

subplot(3,1,3)
plot(met7r.daten, met7r.barometric_pressure, 'b.','markersize',10); hold on;
plot(met7t.daten, met7t.barometric_pressure, 'k.','markersize',5)
datetick('x')
ylabel('mbar')
legend('recovered','telemetered','location','southeast')
title('OOI Irminger METBKA Barometric Pressure, Year 7: 2020-2021')