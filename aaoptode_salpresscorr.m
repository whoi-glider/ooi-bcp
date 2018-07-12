function O2corr = aaoptode_salpresscorr(O2raw,temp,sal,press,S0)

% AAOPTODE_SALPRESSCORR:  Calculates oxygen concentration in uM from raw O2
% concentration, correcting for temperature and salinity. Based on Salinity
% and Depth compensation equations in the Aanderaa operating manual (also
% see function sg_aaoptode.m)
% 
% INPUTS:
% O2raw:       O2 concentration output from optode (in uM)
% temp:        temperature in deg C
% sal:         salinity (measured)
% press:       pressure (db) for pressure correction
% S0:          internal salinity setting for optode (default is zero)
%
% OUTPUT:
% O2corr:       Oxygen concentration in uM
% O2sat:        Oxygen percent saturation (%) (100 = saturated)


temps = log((298.15-temp)./(273.15+temp)); %scaled temperature
SB = [-6.24097E-3; %salinity correction coefficients (B0, B1, B2, and B3)
        -6.93498E-3;
        -6.90358E-3;
        -4.29155E-3];
SC = -3.11680E-7; %salinity correction coefficient C0    
D = 0.032; %depth correction coefficient


O2corr = O2raw.*exp((sal - S0).*(SB(1)+SB(2)*temps+SB(3)*temps.^2+SB(4)*temps.^3)...
    + SC.*(sal.^2 - S0.^2))...
    .*(1+press.*D./1000);

end