function [doxy_lagcorr] = glider_lagCorrectFun(glg, tau_in)

%===================================================================
%
% DESCRIPTION:
%   Use this function to apply a lag correction to OOI pre-processed glider
%   data using the temperature-dependent function from the
%   optode-response-time library accompanying Gordon et al. 2020 after
%   having calculated tau using glider_lagAssessFun
%
% INPUT:
%    glg: Structure containing output from glider_lagAssessFun
%    tau_in: boundary layer thickness to use for correction
%
% OUTPUT:
%    doxy_lagcorr: Lag-corrected doxy, calculated from input doxy, mtime,
%    and temp in the glg structure
%
% DEPENDENCY:
%   optode-response-time library with updates on Kristen Fogaren-owned branch to
%   calculate_tau_w_Temp
%
% AUTHOR:   Hilary Palevsky, 6 February 2023
% -----------------------------------------------------------------------------

[m,n] = size(glg.doxy);

%Initialize array to hold corrected oxygen data
doxy_lagcorr = NaN(size(glg.doxy));

for i = 1:m
    ind = ~(isnan(glg.mtime(i,:)) | isnan(glg.doxy(i,:)) | isnan(glg.temp(i,:)));
    if sum(ind) > 10
        doxy_lagcorr(i,ind) = correct_oxygen_profile_wTemp(glg.mtime(i,ind), glg.doxy(i,ind), glg.temp(i,ind), tau_in);
    end
end


end