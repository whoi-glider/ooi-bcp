%% Comparison of mixed layer oxygen across multiple approaches

%% Load and plot the fixed asset oxygen product for yrs 5-8 produced by Kristen
addpath('G:\Shared drives\NSF_Irminger\OOI_DO_fixed_depth\Data\mixed_layer')
load mixed_layer_calibrated_oxygen.mat

figure(100); clf
h1 = plot(ML_DO.DOdn, ML_DO.DO_umolkg_final, 'k.'); hold on;

%% Glider mixed layer oxygen, deep isotherm-corrected
mlmax = 5; %depth max to use for extracting glider data
slope_pick = 0.2; %slope for glider splash correction
smoothval = 30;

for i = 1:length(glgmerge)
    glg = glgmerge{i};
    indt = find(isnan(glg.time_start) == 0);
    glgmerge{i}.doxy_lagcorr_ml = nanmedian(glg.doxy_lagcorr_grid(1:mlmax, indt)); %pulls out surface data, don't have gain yet
    
    %Deep Isotherm Gain: interpolate onto glider times
    try
         num_isotherm_gains(i) = length(glgmerge{i}.deepisotherm_gains);
    end
    if num_isotherm_gains(i) == 0
        glgmerge{i}.deepisotherm_gains_interptime = NaN;
    elseif num_isotherm_gains(i) < 10 %too short for a full time series
        glgmerge{i}.deepisotherm_gains_interptime = median(glgmerge{i}.deepisotherm_gains)*ones(length(indt),1);
    else %long enough for full time series
        deepisotherm_movmedian = movmedian(glgmerge{i}.deepisotherm_gains, 10);
        [uniquetimes, uniquetimeind, ~] = unique(glgmerge{i}.deepisotherm_times);
        glgmerge{i}.deepisotherm_gains_interptime = interp1(uniquetimes, deepisotherm_movmedian(uniquetimeind), glgmerge{i}.time_start(indt), 'nearest', 'extrap');
    end
    
    %Aircal Gain: -- need to interpolate times
    try
        air_corr_slopeset = (glgmerge{i}.Taircal.air_meas_dist(:,10)-slope_pick.*glgmerge{i}.Taircal.ml_o2sat)./(1-slope_pick);
        air_corr_gain = movmean(glgmerge{i}.Taircal.met_o2sat./air_corr_slopeset, 60);
        glgmerge{i}.aircal_gains_interptime = interp1(glgmerge{i}.Taircal.air_daten, air_corr_gain, glgmerge{i}.time_start(indt), 'nearest', 'extrap');
    catch
        glgmerge{i}.aircal_gains_interptime = NaN;
    end
        
    h2 = plot(glgmerge{i}.time_start(indt), movmean(glgmerge{i}.doxy_lagcorr_ml,smoothval,'omitnan'), 'b.'); hold on;
    h3 = plot(glgmerge{i}.time_start(indt), movmean(glgmerge{i}.doxy_lagcorr_ml.*glgmerge{i}.deepisotherm_gains_interptime',smoothval,'omitnan'), 'm.'); hold on;
    h4 = plot(glgmerge{i}.time_start(indt), movmean(glgmerge{i}.doxy_lagcorr_ml.*glgmerge{i}.aircal_gains_interptime',smoothval,'omitnan'), 'g.'); hold on;
end

%% Plot mixed layer winklers
tol = 12;
for yr = 1:length(btlsum)
    btlsumyr = btlsum{yr};
    if length(btlsumyr) > 0
        for i = 1:length(btlsumyr)
            btlsumcast = btlsumyr(i);
            try
            indsurf = find(btlsumcast{1}.depth < tol);
            %Plot Winkler values
                  try
                       plot(datenum(btlsumcast{1}.Date(indsurf)), btlsumcast{1}.Winkler1_umolkg(indsurf), 'ko','markerfacecolor','r'); hold on;
                  end
                  try
                       h5 = plot(datenum(btlsumcast{1}.Date(indsurf)), btlsumcast{1}.Winkler2_umolkg(indsurf), 'ko','markerfacecolor','r'); hold on;
                  end
                  try
                      plot(datenum(btlsumcast{1}.Date(indsurf)), btlsumcast{1}.Winkler_umolkg(indsurf), 'ko','markerfacecolor','r'); hold on;
                  end
                  try
                      plot(datenum(btlsumcast{1}.Date(indsurf)), btlsumcast{1}.Winkler1_HIP_umolkg(indsurf), 'ko','markerfacecolor','r'); hold on;
                  end
            end
        end
    end
end

%% Add labels to plot
title('OOI Irminger mixed layer dissolved oxygen')
legend([h1 h2 h3 h4 h5], 'Calibrated buoy/NSIF moored sensors','Uncorrected glider','Deep isotherm-corrected glider','Air cal-corrected glider','Winklers','location','northwest')
datetick('x',2)
ylabel('\mumol/kg')

