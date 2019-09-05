function [optode_uM, optode_umolkg] = aaoptode_sternvolmer(foil_coef, optode_phase, tempc, salin, press)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AAOPTODE_STERNVOLMER: For Aanderaa optode 4XXX series, to interpret
% calphase measurements based on the Stern Volmer equation. To be used in
% combination with aaoptode_salpresscorr, which handles salinity and
% pressure corrections.

%% INPUTS
% foil_coef:    Struct of calibration coefficients
% optode_phase: vector of optode phase measurements (calphase or dphase)
% tempc:        temperature in deg C
% salin:        Salinity
% press:        pressure (db) for pressure correction

%% OUTPUTS
% optode_uM: Oxygen concentration in uM
% optode_umolkg: Oxygen concentration in umol/kg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H. Palevsky, 8/6/2019, based on sg_aaoptode.m


% Absolute salinity scale
SA = 35.16504./35.*salin;
% scaled temperature
temps = log((298.15-tempc)./(273.15+tempc));

% Pressure compensation coefficient should be 0.032 for all 4XXX sensors
    D = 0.032;
    % Uchida et al. 2008 Stern-Volmer based calbration mode
        C = foil_coef;
        Ksv = C(1) + C(2).*tempc+C(3).*tempc.^2;
        P0 = C(4) + C(5).*tempc;
        PC = C(6) + C(7).*optode_phase;
        optode_uM = ((P0./PC)-1)./Ksv;
        % note - this uses SR, reference salinity 
        SA = 35.16504./35.*salin;
        CT = gsw_CT_from_t(SA,tempc,press);
        optode_umolkg = 1000.*optode_uM./(1000+ gsw_sigma0(SA,CT));     
end