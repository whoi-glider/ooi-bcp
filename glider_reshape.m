function [glg] = glider_reshape(GLin)

%===================================================================
%
% DESCRIPTION:
%   Use this function to reshape OOI pre-processed glider data to
%   the shape for implementing a lag correction using the
%   optode-response-time library accompanying Gordon et al. 2020
%
% INPUT:
%    GLin: Table of input glider data
%
% OUTPUT:
%    glg: Structure containing profilelist, profile_direction
%           Reshaped versions of original input data for each profile in profilelist:
%               Variables are mtime, pres, doxy, temp, sal, lat, lon,
%               depth, chl, backscatter
%
%
% AUTHOR:   Hilary Palevsky, 15 March 2023
% -----------------------------------------------------------------------------

%Identify profiles (only during dive/climb, not the intermediate data while
%turning at bottom or at surface)
profilelist = unique(floor(unique(GLin.profile_index)));
profilelist = profilelist(find(profilelist > 0)); %exclude value of 0
%Include profilelist in output
glg.profilelist = profilelist;

%Choices kept constant for this analysis but could be varied
tol = 10; %only use profiles with at least tol points

%Add interpolation for oxygen values that aren't already interpolated
%     ind_nonan = find(isnan(GLin.O2_corr) == 0);
% GLin.O2_corr_interp = interp1(GLin.daten(ind_nonan), GLin.O2_corr(ind_nonan), GLin.daten);
% GLin.O2sat_corr_interp = interp1(GLin.daten(ind_nonan), GLin.O2sat_corr(ind_nonan), GLin.daten);

%Reshape data into format that goes into Gordon function
len = 10000;
    glg.mtime = NaN(length(profilelist),len);
    glg.pres = NaN(length(profilelist),len);
    glg.doxy = NaN(length(profilelist),len);
    glg.O2 = NaN(length(profilelist),len);
    glg.O2sat = NaN(length(profilelist),len);
    glg.temp = NaN(length(profilelist),len);
    glg.sal = NaN(length(profilelist),len);
    glg.lat = NaN(length(profilelist),len);
    glg.lon = NaN(length(profilelist),len);
    glg.depth = NaN(length(profilelist),len);
    glg.chl = NaN(length(profilelist),len);
    glg.backscatter = NaN(length(profilelist),len);
    glg.profile_direction = NaN(length(profilelist),1);
for i = 1:length(profilelist)
    indp = find(GLin.profile_index == profilelist(i));
    if length(indp) > tol
        glg.profile_direction(i) = GLin.profile_direction(indp(1));
        glg.mtime(i,1:length(indp)) = GLin.daten(indp);
        glg.pres(i,1:length(indp)) = GLin.pressure_interp(indp);
        glg.doxy(i,1:length(indp)) = GLin.O2_corr(indp);
        glg.O2(i,1:length(indp)) = GLin.O2_corr(indp);
        glg.O2sat(i,1:length(indp)) = GLin.O2sat_corr(indp);
        glg.temp(i,1:length(indp)) = GLin.temperature_interp(indp);
        glg.sal(i,1:length(indp)) = GLin.salinity_interp(indp);
        glg.lat(i,1:length(indp)) = GLin.lat_interp(indp);
        glg.lon(i,1:length(indp)) = GLin.lon_interp(indp);
        glg.depth(i,1:length(indp)) = GLin.depth_interp(indp);
        try
            glg.chl(i,1:length(indp)) = GLin.chlorophyll(indp);
            glg.backscatter(i,1:length(indp)) = GLin.backscatter(indp);
        end
    end
end