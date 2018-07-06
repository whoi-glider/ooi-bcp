function [stats] = glideraircal_stats(glider, spacer, tol)
% This function identifies data from air calibration intervals, removes
% data (length spacer) from the beginning and end of the interval, and calculates
% mean and stdev for remaining data, and flags as usable all intervals with
% 3 or more data points and stdev of O2 saturation below a set tolerance
%-----------------------------------------------------------------------
%%%% INPUTS
% glider = glider data in table form, including:
    %oxygen_concentration
    %oxygen_saturation
    %profile_index
    %daten
% spacer = number of points to remove from both beginning and end of air
    % interval as contaminated data (drying off/losing buoyance prior to dive)
% tol = tolerance for maximum stdev of O2 conc for data to be flagged as
% good

%%%% OUTPUTS
% stats = matrix with row for each air calibration interval, after removing
% spacer points. Columns are:
    %1 = mean corrected oxygen (O2_corr)
    %2 = stdev orrected oxygen (O2_corr)
    %3 = mean oxygen saturation
    %4 = stdev oxygen saturation
    %5 = number of good points
    %6 = mean date/time of air interval
    %7 = flag (1 = good, NaN otherwise)

%-------------------------------------------------------------------------


d = rem(glider.profile_index,1) == 0.5;
intervals = unique(glider.profile_index(d));
stats = NaN*ones(length(intervals), 7);
for i = 1:length(intervals)
    ind = find(glider.profile_index == intervals(i));
    stats(i,1) = nanmean(glider.O2_corr(ind(spacer + 1:end - spacer)));
    stats(i,2) = nanstd(glider.O2_corr(ind(spacer + 1:end - spacer)));
    stats(i,3) = nanmean(glider.oxygen_saturation(ind(spacer + 1:end - spacer)));
    stats(i,4) = nanstd(glider.oxygen_saturation(ind(spacer + 1:end - spacer)));
    stats(i,5) = sum(~isnan(glider.O2_corr(ind(spacer + 1:end - spacer))));
    stats(i,6) = nanmean(glider.daten(ind(spacer + 1:end - spacer)));
    if stats(i,4) < tol & stats(i,5) > 2 % count data as good if 3 or more measurements and std of O2 conc < tolerance
        stats(i,7) = 1; %put a "good" flag on usable data
    end
end