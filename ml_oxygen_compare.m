%% Comparison of mixed layer oxygen across multiple approaches

%% Load and plot the fixed asset oxygen product for yrs 5-8 produced by Kristen
addpath('G:\Shared drives\NSF_Irminger\OOI_DO_fixed_depth\Data\mixed_layer')
load mixed_layer_calibrated_oxygen.mat

figure(100); clf
h1 = plot(ML_DO.DOdn, ML_DO.DO_umolkg_final, 'k.'); hold on;

%% Glider mixed layer oxygen, deep isotherm-corrected
mlmax = 5; %depth max to use for extracting glider data
slope_pick = 0.2; %slope for glider splash correction
smoothval = 10;

for i = 1:length(glgmerge)
    glg = glgmerge{i};
    indt = find(isnan(glg.time_start) == 0);
    glgmerge{i}.doxy_lagcorr_ml = glg.doxy_lagcorr_grid(5, indt); %pulls out surface data before gain is applied
        
    %Aircal Gain: -- need to interpolate times
    try
        air_corr_slopeset = (glgmerge{i}.Taircal.air_meas_dist(:,10)-slope_pick.*glgmerge{i}.Taircal.ml_o2sat)./(1-slope_pick);
        air_corr_gain = movmean(glgmerge{i}.Taircal.met_o2sat./air_corr_slopeset, 60);
        glgmerge{i}.aircal_gains_interptime = interp1(glgmerge{i}.Taircal.air_daten, air_corr_gain, glgmerge{i}.time_start(indt), 'nearest', 'extrap');
    catch
        glgmerge{i}.aircal_gains_interptime = NaN;
    end        
    h2 = plot(glgmerge{i}.time_start(indt), movmean(glgmerge{i}.doxy_lagcorr_ml,smoothval,'omitnan'), 'b.'); hold on;
    h4 = plot(glgmerge{i}.time_start(indt), movmean(glgmerge{i}.doxy_lagcorr_ml.*glgmerge{i}.aircal_gains_interptime',smoothval,'omitnan'), '.','color',nicecolor('cccbw')); hold on;
end

h3 = plot(glidermerge.time, movmean(glidermerge.doxy(5,:),smoothval,'omitnan'), '.','color',nicecolor('bbbbbrmmww')); hold on;

%% Plot mixed layer winklers
tol = 20;
winkmin = 270;
for yr = 1:length(btlsum)
    btlsumyr = btlsum{yr};
    if length(btlsumyr) > 0
        for i = 1:length(btlsumyr)
            btlsumcast = btlsumyr(i);
            try
            indsurf = find(btlsumcast{1}.depth < tol);
            %Plot Winkler values
                  try
                      %indgood = find(btlsumcast{1}.NLMR_Outlier1 == 2);
                      %indplot = intersect(indsurf,indgood); 
                      plot(datenum(btlsumcast{1}.Date(indsurf)), btlsumcast{1}.Winkler1_umolkg(indsurf), 'ko','markerfacecolor','r'); hold on;
                  end
                  try
                      %indgood = find(btlsumcast{1}.NLMR_Outlier2 == 2);
                      %indplot = intersect(indsurf,indgood);  
                      h5 = plot(datenum(btlsumcast{1}.Date(indsurf)), btlsumcast{1}.Winkler2_umolkg(indsurf), 'ko','markerfacecolor','r'); hold on;
                  end
                  try
                      indgood = find(btlsumcast{1}.NLMR_Outlier < 3 & btlsumcast{1}.Winkler_umolkg > winkmin);
                      indplot = intersect(indsurf,indgood);
                      plot(datenum(btlsumcast{1}.Date(indplot)), btlsumcast{1}.Winkler_umolkg(indplot), 'ko','markerfacecolor','r'); hold on;
                  end
                  try
                      indgood = find(btlsumcast{1}.NLMR_HIP1_Outlier == 2);
                      indplot = intersect(indsurf,indgood);
                      plot(datenum(btlsumcast{1}.Date(indplot)), btlsumcast{1}.Winkler1_HIP_umolkg(indplot), 'ko','markerfacecolor','r'); hold on;
                  end
            end
        end
    end
end

%% Add labels to plot
title('OOI Irminger mixed layer dissolved oxygen')
legend([h1 h2 h3 h4 h5], 'Calibrated buoy/NSIF moored sensors','Uncorrected glider','Deep isotherm-corrected glider, 5m','Aircal-corrected glider, 5m','Winklers','location','northwest')
xlim([datenum(2014,6,1) datenum(2022, 9, 1)])
datetick('x',2,'keeplimits')
ylabel('\mumol/kg')

