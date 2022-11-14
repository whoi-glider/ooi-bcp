function G_out = glider_interpCorrFun(G_in, sal_set)

%===================================================================
%
% DESCRIPTION:
%   Use this function to interpolate OOI pre-processed glider data to
%   1) interpolate needed variables, and 2) perform salinity and pressure
%   corrections to O2conc & pressure correction to O2sat
%
% INPUT:
%    G_in: Table of input glider data
%    sal_set: internal salinity setting on glider
%
% OUTPUT:
%    G_out: Original glider input table with added columns:
%           lon_interp, lat_interp, salinity_interp, pressure_interp,
%           temperature_interp - interpolated variables
%           O2_corr, O2sat_corr - salinity and pressure corrected
%
% DEPENDENCY:
%   AAOPTODE_SALPRESSCORR
%
% AUTHOR:   Hilary Palevsky, 14 November 2022


G_in.lon_interp = naninterp1(G_in.time, G_in.longitude, G_in.time);
G_in.lat_interp = naninterp1(G_in.time, G_in.latitude, G_in.time);
G_in.salinity_interp = naninterp1(G_in.time, G_in.salinity, G_in.time);
G_in.pressure_interp = naninterp1(G_in.time, G_in.pressure, G_in.time);
G_in.temperature_interp = naninterp1(G_in.time, G_in.temperature, G_in.time);
G_in.O2_corr = aaoptode_salpresscorr(G_in.oxygen_concentration, G_in.temperature_interp, G_in.salinity_interp, G_in.pressure_interp, sal_set);
G_in.O2sat_corr = G_in.oxygen_saturation.*(1+G_in.pressure_interp.*0.032./1000); %applies pressure but not salinity correction for saturation

G_out = G_in;

end
