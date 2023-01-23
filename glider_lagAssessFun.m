function [glg] = glider_lagAssessFun(GLin, beg, stop, ymax, name, profilelist)

%===================================================================
%
% DESCRIPTION:
%   Use this function to interpolate OOI pre-processed glider data to
%   calculate the lag correction using the temperature-dependent function
%   from the optode-response-time library accompanying Gordon et al. 2020
%
% INPUT:
%    GLin: Table of input glider data
%    beg: index of first data point for plotting
%    stop: index final data point for plotting
%    ymax: depth maximum of glider data for analysis & plotting
%    name: string with name of glider, for plot title
%    profilelist: list of profile numbers to use for paired lag analysis
%
% OUTPUT:
%    glg: Structure containing
%           Reshaped versions of original input data for each profile in profilelist:
%               Variables are mtime, pres, doxy, temp
%           Output from lag analysis with calculate_tau_w_Temp:
%               thickness - boundary layer output in um for each profile pair
%               thickness_constants - input thicknesses analyzed
%               rmsdt - root mean squared error values for each thickness analyzed for each profile pair
%               tau_Tref - tau (sec) at chosen ref temp, set to 4 deg C
%
% DEPENDENCY:
%   optode-response-time library with updates on Kristen Fogaren-owned branch to
%   calculate_tau_w_Temp
%
% AUTHOR:   Hilary Palevsky, 23 January 2023
% -----------------------------------------------------------------------------

%Selection of depth range over which to compare profiles to calculate lag
zlim_depth = [20 ymax*.95]; %use 95% of maximum depth, exclude top 20 m

%Choices kept constant for this analysis but could be varied
tol = 50; %only use profiles with at least 50 points
tref = 4; %reference temperature for tau in T-dependent lag function (but will use thickness, not time)
jump = 2; %plotting interval

%Reshape data into format that goes into Gordon function
len = 10000;
    glg.mtime = NaN(length(profilelist),len);
    glg.pres = NaN(length(profilelist),len);
    glg.doxy = NaN(length(profilelist),len);
    glg.temp = NaN(length(profilelist),len);
for i = 1:length(profilelist)
    indp = find(GLin.profile_index == profilelist(i));
    if length(indp) > tol
        glg.mtime(i,1:length(indp)) = GLin.daten(indp);
        glg.pres(i,1:length(indp)) = GLin.pressure_interp(indp);
        glg.doxy(i,1:length(indp)) = GLin.O2_corr(indp);
        glg.temp(i,1:length(indp)) = GLin.temperature_interp(indp);
    end
end

%Calculate tau usinv version with temperature term, in depth space
    %Omitting version in density space because yields more NaN results
[glg.thickness, glg.tau_Tref, glg.thickness_constants, glg.rmsdt] = calculate_tau_wTemp(glg.mtime, glg.pres, glg.doxy, glg.temp,...
    'zlim',zlim_depth,'tref', tref);



figure; clf
    subplot(2,2,[1:2])
plot(GLin.daten(beg:jump:stop), GLin.depth_interp(beg:jump:stop), 'k.'); hold on;
scatter(GLin.daten(beg:jump:stop), GLin.depth_interp(beg:jump:stop), [], GLin.oxygen_saturation(beg:jump:stop),'filled'); colorbar; caxis([85 100])
set(gca,'YDir','reverse'); 
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
ylim([-10 ymax])
title({[name ', Initial paired up and down profiles (oxygen % saturation)']})

    subplot(2,2,3)
plot(glg.thickness_constants, glg.rmsdt,'.')
xlabel('thickness (\mum)')
ylabel('RMSD (%)')
title('w/ T term, depth aligned')

    subplot(2,2,4)
plot(glg.doxy(1:2:end), glg.pres(1:2:end),'k.'); hold on;
plot(glg.doxy(2:2:end), glg.pres(2:2:end),'b.'); hold on;
axis ij
xlabel('Oxygen (uncorr)')
ylabel('Pressure (db)')
title('Individual up & down profiles')


end