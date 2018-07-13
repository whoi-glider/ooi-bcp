%% Extract meteorological data from Apex surface mooring
  
filename = ['deployment0005_GI01SUMO-SBD11-06-METBKA000-telemetered-metbk_a_dcl_instrument_20180608T172109.234000-20180713T150853.321000.nc'];
%ncdisp(filename) %To show all variablies in the netcdf file

%Pull out desired variables
    SUMO_met.time = ncread(filename,'time');
        SUMO_met.time = convertTime(SUMO_met.time); %convert to matlab time
    SUMO_met.barometric_pressure = ncread(filename,'barometric_pressure'); %mbar
    SUMO_met.relative_humidity = ncread(filename,'relative_humidity'); % in percent
    
%Plot data to check
figure(1); clf
    subplot(211)
plot(SUMO_met.time, SUMO_met.barometric_pressure, 'k.');
title('SUMO barometric pressure, mbar'); datetick('x');
    subplot(212)
plot(SUMO_met.time, SUMO_met.relative_humidity, 'k.');
title('SUMO relative humidity, %'); datetick('x');
