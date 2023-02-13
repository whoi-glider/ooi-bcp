function [G_out] = gliderGrid(G_in, pmin, pmax, pinterval)

%===================================================================
%
% DESCRIPTION:
%   Use this function to grid glider output data in the format output by
%   glider_lagAssessFun and glider_lagCorrectFun (based on the format for
%   optode-response-time library) into discrete bins, taking the mean over
%   each bin defined by a specific pressure interval
%
% INPUT:
%    G_in: Structure output by glider_lagAssessFun, with doxy_lagcorr added
%    from glider_lagCorrectFun
%    pmin: minimum pressure for gridding
%    pmax: maximum pressure for gridding
%    pinterval: pressure interval for gridding
%
% OUTPUT:
%    G_out: G_in with added gridded temp, doxy, and doxy_lagcorr variables (mean
%    and standard deviation) along pgrid
%
%
% AUTHOR:   Hilary Palevsky, 6 February 2023
% -----------------------------------------------------------------------------

pgrid = [pmin + pinterval/2: pinterval: pmax];

[m,n] = size(G_in.doxy);
    G_in.doxy_gridmean = NaN*ones(m,length(pgrid));
    G_in.doxy_gridstd = NaN*ones(m,length(pgrid));
    G_in.doxy_lagcorr_gridmean = NaN*ones(m,length(pgrid));
    G_in.doxy_lagcorr_gridstd = NaN*ones(m,length(pgrid));
    G_in.temp_gridmean = NaN*ones(m,length(pgrid));
    G_in.temp_gridstd = NaN*ones(m,length(pgrid));
for i = 1:m
    for j = 1:length(pgrid)
        ind = find(G_in.pres(i,:) > pgrid(j) - 2.5 & G_in.pres(i,:) <= pgrid(j) + 2.5);
        if length(ind) > 1
            G_in.doxy_gridmean(i,j) = nanmean(G_in.doxy(i,ind));
            G_in.doxy_gridstd(i,j) = nanstd(G_in.doxy(i,ind));
            G_in.doxy_lagcorr_gridmean(i,j) = nanmean(G_in.doxy_lagcorr(i,ind));
            G_in.doxy_lagcorr_gridstd(i,j) = nanstd(G_in.doxy_lagcorr(i,ind));
            G_in.temp_gridmean(i,j) = nanmean(G_in.temp(i,ind));
            G_in.temp_gridstd(i,j) = nanstd(G_in.temp(i,ind));
        end
    end
end

G_out = G_in;
G_out.pgrid = pgrid;

end