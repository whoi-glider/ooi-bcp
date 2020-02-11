%% Modified from Lucy's thesis

%Make initial gain correction    
Yr5_wfp.oxygen_corr = Yr5_wfp.oxygen * gain_hypm; 
Yr5_wfpgrid.oxygen_corr = Yr5_wfpgrid.O2conc * gain_hypm; 
Yr5_wfpgrid_therm.oxygen_corr = Yr5_wfpgrid_therm.O2conc * gain_hypm; 
Yr5_wfpgrid.O2satcorr = (Yr5_wfpgrid.oxygen_corr./Yr5_wfpgrid.O2equil - 1)*100;

%% Pin everything to stable deep ocean measurement in deep WFP data - plotting on isotherms
% Use rate of change on deep isotherms to correct for drift in WFP O2 and salinity
figure(1); clf
thermplot = Yr5_wfpgrid_therm.therm_grid(25:8:43); %select stable deep isotherms
thermstr = {'2.3 deg C', '2.7 deg C', '3.1 deg C'};
clear C h; bot = nicecolor('rrywwwwww'); top = nicecolor('rrykkkkkk');
C = [linspace(top(1),bot(1),length(thermplot))' linspace(top(2),bot(2),length(thermplot))' linspace(top(3),bot(3),length(thermplot))'];
yearspan = {'2018-2019'};
M = 15;

plotting = Yr5_wfpgrid_therm;
plotting.time_start = Yr5_wfpgrid_therm.time_start(Yr5_wfpgrid_therm.ind_pair);
O2gaincorr = plotting.oxygen_corr; %plot gain corr O2 -- you want what was originally Yr5_wfpgrid_therm.oxygen_corr, which is now also plotting.oxygen_corr
%Initialize arrays to hold slopes
    O2slope = NaN*ones(length(plotting.therm_grid),2);
    O2slope_err = NaN*ones(length(plotting.therm_grid),2);
    Tslope = NaN*ones(length(plotting.therm_grid),2);
    Tslope_err = NaN*ones(length(plotting.therm_grid),2);
    Sslope = NaN*ones(length(plotting.therm_grid),2);
    Sslope_err = NaN*ones(length(plotting.therm_grid),2);

%Set up to calculate varying drift rate over time
tstep = 90; %number of days in each time period for derivative plotting
tmin = min(plotting.time_start); tmax = max(plotting.time_start);
t = [floor(tmin)+1:ceil(tmax)];
tgrid = [floor(tmin)+1:tstep:tmax]-floor(tmin); tgrid = [tgrid max(t)-min(t)+1];
tgrid2 = [1 60 max(t)-min(t)+1];

%%% Look for trends over time in individual isotherms
for k = 1:length(plotting.therm_grid)
    %ind = find(~isnan(O2gaincorr(k,:)));
    ind = find(O2gaincorr(k,:) > nanmean(O2gaincorr(k,:)) - 2*nanstd(O2gaincorr(k,:)) & O2gaincorr(k,:) < nanmean(O2gaincorr(k,:)) + 2*nanstd(O2gaincorr(k,:)));
    [P,S] = polyfit(plotting.time_start(ind), O2gaincorr(k,ind)',1);
    O2slope(k,i) = P(1); O2int(k) = P(2);
    for j = 1:length(tgrid2)-1
        tind = find(plotting.time_start >= tgrid2(j) + floor(tmin) & plotting.time_start < tgrid2(j+1) + floor(tmin));
        intind = intersect(tind,ind);
        numpts_pieceslope(k,j,i) = length(intind);
        P = polyfit(plotting.time_start(intind), O2gaincorr(k,intind)',1);
        O2pieceslope(k,j,i) = P(1);
    end
    [P,S] = polyfit(plotting.time_start(ind), plotting.T(k,ind)',1);
    Tslope(k,i) = P(1); Tint(k) = P(2);
    if sum(size(S.R)) == 4
        err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
        Tslope_err(k,i) = err(1); Tint_err(k) = err(2);
    end
    [P,S] = polyfit(plotting.time_start(ind), plotting.S(k,ind)',1);
    Sslope(k,i) = P(1); Sint(k) = P(2);
    if sum(size(S.R)) == 4
        err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
        Sslope_err(k,i) = err(1); Sint_err(k) = err(2);
    end
end
%%% Plot trends at specific isotherms
%subplot(1,3,1)
    plottime = [1:ceil(length(plotting.time_start)*1)];
for j = 1:length(thermplot)
    idtherm = find(plotting.therm_grid == thermplot(j));
    if length(idtherm) == 1
        h(j) = plot(plotting.time_start, O2gaincorr(idtherm,:),'.','color',C(j,:),'markersize',M); hold on;
        ylabel ('O_2 (\mumol/kg)', 'Fontsize', 15)
        %set(get(gca,'ylabel'),'rotation',0)
    end
end
datetick('x',3); ylim([268 295]); %legend(h, thermstr)
title('Raw oxygen concentration on WFP isotherms below winter ventilation')
legend(h,thermstr,'location','northeast')

% subplot(1,3,2)
% for j = 1:length(thermplot)
%     idtherm = find(plotting.therm_grid == thermplot(j));
%     if length(idtherm) == 1
%         h(j) = plot(plotting.time_start, plotting.S(idtherm,:),'.','color',C(j,:),'markersize',M); hold on;
%         plot(plotting.time_start, Sint(idtherm) + plotting.time_start*Sslope(idtherm,i),'color',C(j,:),'linewidth',1); hold on;
%     end
% end
% datetick('x',3); ylim([34.86 34.96]); %legend(h,thermstr,'location','southwest')
% title(['Salinity from WFP isotherms below winter ventilation, ' yearspan])
% 
% subplot(1,3,3)
% for j = 1:length(thermplot)
%     idtherm = find(plotting.therm_grid == thermplot(j));
%     if length(idtherm) == 1
%         h(j) = plot(plotting.time_start, plotting.depth(idtherm,:),'.','color',C(j,:),'markersize',M); hold on;
%     end
% end
% datetick('x',3); set(gca,'ydir','reverse'); ylim([1400 2700]); legend(h,thermstr,'location','southeast')
% title(['Depth of WFP isotherms below winter ventilation, ' yearspan{i}])

%% Calculate mean drift rate based on all non-outlier stable deep isotherms
stable = find(plotting.therm_grid >= min(thermplot) & plotting.therm_grid <= max(thermplot)); %indices of stable deep isotherms
for i = 1
    wfp_O2drift(i) = nanmean(O2slope(stable,i));
    wfp_O2drift_std(i) = nanstd(O2slope(stable,i));
    wfp_O2drift_num(i) = sum(~isnan(O2slope(stable,i)));
end


