%% Values for separating out only nighttime values
    lat = 60;
    lon = -40;
    tol = 2;
%%    
filename = 'deployment0001_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20140910T185041.465000-20150216T054857.646000.nc';
sumo{1}.daten = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{1}.SST = ncread(filename,'sea_surface_temperature');

%No NSIF DOSTA data

filename = 'deployment0001_GI01SUMO-SBD11-04-DOSTAD000-telemetered-dosta_abcdjm_dcl_instrument_20140910T190020.500000-20140918T150300.020000.nc';
sumo{1}.datenO2buoy = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{1}.O2buoy = ncread(filename,'dissolved_oxygen');
[~,sumo{1}.nightindO2buoy] = indexDayNight(lat,lon,0,sumo{1}.datenO2buoy,tol);
%%
filename = 'deployment0002_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20150815T192214.032000-20160127T080842.876000.nc';
sumo{2}.daten = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{2}.SST = ncread(filename,'sea_surface_temperature');

filename = 'deployment0002_GI01SUMO-RID16-06-DOSTAD000-telemetered-dosta_abcdjm_dcl_instrument_20150815T193009.886000-20160718T234802.768000.nc';
sumo{2}.datenO2nsif = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{2}.O2nsif = ncread(filename,'dissolved_oxygen');
[~,sumo{2}.nightindO2nsif] = indexDayNight(lat,lon,0,sumo{2}.datenO2nsif,tol);

filename = 'deployment0002_GI01SUMO-SBD11-04-DOSTAD000-telemetered-dosta_abcdjm_dcl_instrument_20150815T193013.458000-20160130T083306.921000.nc';
sumo{2}.datenO2buoy = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{2}.O2buoy = ncread(filename,'dissolved_oxygen');
[~,sumo{2}.nightindO2buoy] = indexDayNight(lat,lon,0,sumo{2}.datenO2buoy,tol);
%%
filename = 'deployment0003_GI01SUMO-SBD11-06-METBKA000-recovered_host-metbk_a_dcl_instrument_recovered_20160710T174403.350000-20170812T200352.175000.nc';
sumo{3}.daten = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{3}.SST = ncread(filename,'sea_surface_temperature');

filename = 'deployment0003_GI01SUMO-RID16-06-DOSTAD000-telemetered-dosta_abcdjm_dcl_instrument_20160710T174508.675000-20170813T140300.440000.nc';
sumo{3}.datenO2nsif = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{3}.O2nsif = ncread(filename,'dissolved_oxygen');
[~,sumo{3}.nightindO2nsif] = indexDayNight(lat,lon,0,sumo{3}.datenO2nsif,tol);

filename = 'deployment0003_GI01SUMO-SBD11-04-DOSTAD000-telemetered-dosta_abcdjm_dcl_instrument_20160710T174509.583000-20170814T044804.454000.nc';
sumo{3}.datenO2buoy = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{3}.O2buoy = ncread(filename,'dissolved_oxygen');
[~,sumo{3}.nightindO2buoy] = indexDayNight(lat,lon,0,sumo{3}.datenO2buoy,tol);
%%
filename = 'deployment0004_GI01SUMO-SBD11-06-METBKA000-telemetered-metbk_a_dcl_instrument_20170805T181721.211000-20171012T090753.398000.nc';
sumo{4}.daten = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{4}.SST = ncread(filename,'sea_surface_temperature');

filename = 'deployment0004_GI01SUMO-RID16-06-DOSTAD000-telemetered-dosta_abcdjm_dcl_instrument_20170805T181700.505000-20171012T090301.142000.nc';
sumo{4}.datenO2nsif = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{4}.O2nsif = ncread(filename,'dissolved_oxygen');
[~,sumo{4}.nightindO2nsif] = indexDayNight(lat,lon,0,sumo{4}.datenO2nsif,tol);

filename = 'deployment0004_GI01SUMO-SBD11-04-DOSTAD000-telemetered-dosta_abcdjm_dcl_instrument_20170805T181701.547000-20171012T090304.674000.nc';
sumo{4}.datenO2buoy = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{4}.O2buoy = ncread(filename,'dissolved_oxygen');
[~,sumo{4}.nightindO2buoy] = indexDayNight(lat,lon,0,sumo{4}.datenO2buoy,tol);
%%
filename = 'deployment0005_GI01SUMO-SBD12-06-METBKA000-telemetered-metbk_a_dcl_instrument_20180608T172154.969000-20190809T080328.237000.nc';
sumo{5}.daten = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{5}.SST = ncread(filename,'sea_surface_temperature');

filename = 'deployment0005_GI01SUMO-RID16-06-DOSTAD000-telemetered-dosta_abcdjm_dcl_instrument_20180608T173009.294000-20190809T080301.690000.nc';
sumo{5}.datenO2nsif = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{5}.O2nsif = ncread(filename,'dissolved_oxygen');
[~,sumo{5}.nightindO2nsif] = indexDayNight(lat,lon,0,sumo{5}.datenO2nsif,tol);

filename = 'deployment0005_GI01SUMO-SBD11-04-DOSTAD000-telemetered-dosta_abcdjm_dcl_instrument_20180608T173012.719000-20190809T080318.357000.nc';
sumo{5}.datenO2buoy = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{5}.O2buoy = ncread(filename,'dissolved_oxygen');
[~,sumo{5}.nightindO2buoy] = indexDayNight(lat,lon,0,sumo{5}.datenO2buoy,tol);
%%
filename = 'deployment0006_GI01SUMO-SBD11-06-METBKA000-telemetered-metbk_hourly_20190805T155704.097000-20190827T163006.467000.nc';
sumo{6}.daten = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{6}.SST = ncread(filename,'sea_surface_temperature');
%sumo{6}.wind10m = ncread(filename,'met_wind10m');

filename = 'deployment0006_GI01SUMO-RID16-06-DOSTAD000-telemetered-dosta_abcdjm_dcl_instrument_20190805T153008.599000-20190914T150302.580000.nc';
sumo{6}.datenO2nsif = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{6}.O2nsif = ncread(filename,'dissolved_oxygen');
[~,sumo{6}.nightindO2nsif] = indexDayNight(lat,lon,0,sumo{6}.datenO2nsif,tol);

filename = 'deployment0006_GI01SUMO-SBD11-04-DOSTAD000-telemetered-dosta_abcdjm_dcl_instrument_20190805T153015.838000-20190916T180319.640000.nc';
sumo{6}.datenO2buoy = datenum(1900,1,1,00,0,ncread(filename,'time'));
sumo{6}.O2buoy = ncread(filename,'dissolved_oxygen');
[~,sumo{6}.nightindO2buoy] = indexDayNight(lat,lon,0,sumo{6}.datenO2buoy,tol);

%%
datemin = datenum(2014,9,1);
datemax = datenum(2019,9,16);
set(0,'defaultAxesFontSize',14)

figure(1); clf
subplot(211)
for i = 1:6
    plot(sumo{i}.daten, sumo{i}.SST, 'k.'); hold on;
end
%plot(G453grid.time_start(1:end-100), squeeze(G453grid.scivars(2,1,1:end-100)),'g.'); hold on;
h3 = plot(G363grid.time_start(1:end-100), squeeze(G363grid.scivars(2,1,1:end-100)),'r.')
xlim([datemin datemax])
ylabel('SST (^oC)')
title('OOI Irminger Sea sea surface temperature')
datetick('x',2,'keeplimits')
ylim([2 12])

subplot(212)
for i = 2:6
    %plot(sumo{i}.datenO2nsif, sumo{i}.O2nsif, 'k.'); hold on;
    h1 = plot(sumo{i}.datenO2nsif(sumo{i}.nightindO2nsif), sumo{i}.O2nsif(sumo{i}.nightindO2nsif), 'b.'); hold on;
end
for i = 1:6
    %plot(sumo{i}.datenO2buoy, sumo{i}.O2buoy, 'r.'); hold on;
    h2 = plot(sumo{i}.datenO2buoy(sumo{i}.nightindO2buoy), sumo{i}.O2buoy(sumo{i}.nightindO2buoy), 'k.'); hold on;
end
h3 = plot(G363grid.time_start(1:end-100), squeeze(G363grid.scivars(2,3,1:end-100))*med_gain_363,'r.')
%plot(G453grid.time_start(1:end-100), squeeze(G453grid.scivars(4,3,1:end-100))*med_gain_453,'g.')
xlim([datemin datemax])
ylabel('\mumol/kg')
legend([h1 h2 h3], 'Mooring, 12 m', 'Mooring, 1 m', 'Glider, 10 m','location', 'northwest')
title('OOI Irminger Sea surface dissolved oxygen')
datetick('x',2,'keeplimits')

