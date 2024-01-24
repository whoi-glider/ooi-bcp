function [O2corr, D] = aaoptode_salpresscorr_empiricalD(O2raw,temp,sal,press,S0,calcast_press,calcast_O2uM,numcalcastalign)

% AAOPTODE_SALPRESSCORR:  Calculates oxygen concentration in uM from raw O2
% concentration, correcting for temperature and salinity. Uses Salinity
% compensation equations in the Aanderaa operating manual, but rather than
% using a fixed value of D (depth correction coefficient), finds an
% empirical fit based on alignment with a calibration cast.
% 
% INPUTS:
% O2raw:                O2 concentration output from optode (in uM)
% temp:                 temperature in deg C
% sal:                  salinity (measured)
% press:                pressure (db) for pressure correction
% S0:                   internal salinity setting for optode (default is zero)
% calcast_press         pressure (db) from associated SBE43 calibration cast
% calcast_O2uM          corrected oxygen from associated SBE43 calibration cast
% press_calcastalign    pressure profile from optode profile to compare with calcast
% O2raw_calcastalign    O2 concentration profile to compare with calcast
%
% OUTPUT:
% O2corr:       Oxygen concentration in uM
% D:            Empirically determined pressure coefficient


temps = log((298.15-temp)./(273.15+temp)); %scaled temperature
SB = [-6.24097E-3; %salinity correction coefficients (B0, B1, B2, and B3)
        -6.93498E-3;
        -6.90358E-3;
        -4.29155E-3];
SC = -3.11680E-7; %salinity correction coefficient C0    

O2salcorr = O2raw.*exp((sal - S0).*(SB(1)+SB(2)*temps+SB(3)*temps.^2+SB(4)*temps.^3)...
    + SC.*(sal.^2 - S0.^2)); %remove pressure factor

%Interpolate data from 1st cast in dataset onto calcast pressure grid
    indnonan = ~isnan(O2salcorr(numcalcastalign,:) + press(numcalcastalign,:));
O2salcorr_init_regrid = interp1(press(numcalcastalign,indnonan), O2salcorr(numcalcastalign,indnonan), calcast_press);

%Test range of possible pressure coefficients & calculate standard deviation between calcast and P-corrected profiles
Pc_grid = [0.001:0.001:0.1];
inddeep = find(calcast_press > 300);
M = NaN*ones(length(calcast_press),length(Pc_grid)); %initialize array
err = NaN*ones(1,length(Pc_grid));
for i = 1:length(Pc_grid)
    M(:,i) = O2salcorr_init_regrid.*(1+calcast_press.*Pc_grid(i)./1000);
    err(i) = nanstdev(M(inddeep,i) - calcast_O2uM(inddeep));
end
[~,indPc_opt] = min(err); %calculate minimum
D = Pc_grid(indPc_opt);

%Plot relationship between stdev and and pressure coefficients
figure; clf
plot(Pc_grid, err, 'k.'); hold on;
xline(D)
xlim([min(Pc_grid) max(Pc_grid)])
xlabel('Pressure coefficients')
ylabel('Stdev between optode and CTD oxygen')
title(['Optimal Pc for this optode is ' num2str(D)])

%Plot calcast along with 1st cast in dataset
figure; clf;
plot(calcast_O2uM, calcast_press, 'k.'); hold on;
plot(O2salcorr_init_regrid, calcast_press,'b.'); hold on;
plot(M(:,indPc_opt), calcast_press, 'r.');
axis ij
legend('SBE43 calcast','Optode profile, no P corr','Optode profile, P-compensated')
ylabel('Pressure (db)')
xlabel('Oxygen (\muM)')

%Calculate pressure-compensated oxygen over full dataset
O2corr = O2salcorr.*(1+press.*D./1000);

end