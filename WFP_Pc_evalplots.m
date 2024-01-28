%% Look at oxygen profile on ~Feb 1. of each year

cmap = jet(6);
mldyr = str2num(datestr(chldt_mat,10));
pad = -1;

figure(1); clf
for i = 1:6
    yr = 2014 + i;
    yrindmld = find(mldyr == yr);
    [mldmaxyr, mldmaxind] = max(chlmld(yrindmld));
    datestr(chldt_mat(yrindmld(mldmaxind)));
    [d, t_ind] = min(abs(wggmerge.time - chldt_mat(yrindmld(mldmaxind)-pad)));
    h(i) = plot(movmean(wggmerge.doxy(:,t_ind),20), pres_grid_hypm,'-','linewidth',2,'color',cmap(i,:)); hold on;
    if i == 5
        oxy_mldmaxstd(1,i) = nanstdev(movmean(wggmerge.doxy(101:201,t_ind),20));
    else
        oxy_mldmaxstd(1,i) = nanstdev(movmean(wggmerge.doxy(101:341,t_ind),20));
    end
end
axis ij
ylim([250 1300])
grid on
legend('2015', '2016', '2017', '2018', '2019', '2020','location','SW')
title('WFP oxygen profile at MLDmax')
xlabel('Oxygen, \mumol/kg')
ylabel('Pressure, db')

%% Plots used for comparing selected Pc approach with others
% figure(100); clf
% subplot(211)
% % plot(wggmerge_empD.time, wggmerge_empD.doxy(1501,:),'k.'); hold on;
% % plot(wggmerge_empD.time, movmean(wggmerge_empD.doxy(1501,:),45,'omitnan'),'k-','linewidth',2); hold on;
% plot(wggmerge_fixedD.time, wggmerge_fixedD.doxy(1501,:),'r.'); hold on;
% plot(wggmerge_fixedD.time, movmean(wggmerge_fixedD.doxy(1501,:),45,'omitnan'),'r-','linewidth',2); hold on;
% plot(wggmerge.time, wggmerge.doxy(1501,:),'b.'); hold on;
% plot(wggmerge.time, movmean(wggmerge.doxy(1501,:),45,'omitnan'),'b-','linewidth',2); hold on;
% datetick('x',2)
% ylim([267 282])
% title('Fixed empirical times series mean (red) vs empirical by year (blue) pressure coefficient, 1650 m')
% 
% subplot(212)
%     ind = 101;
% % plot(wggmerge_empD.time, wggmerge_empD.doxy(ind,:),'k.'); hold on;
% % plot(wggmerge_empD.time, movmean(wggmerge_empD.doxy(ind,:),45,'omitnan'),'k-','linewidth',2); hold on;
% plot(wggmerge_fixedD.time, wggmerge_fixedD.doxy(ind,:),'r.'); hold on;
% plot(wggmerge_fixedD.time, movmean(wggmerge_fixedD.doxy(ind,:),45,'omitnan'),'r-','linewidth',2); hold on;
% plot(wggmerge.time, wggmerge.doxy(ind,:),'b.'); hold on;
% plot(wggmerge.time, movmean(wggmerge.doxy(ind,:),45,'omitnan'),'b-','linewidth',2); hold on;
% datetick('x',2)
% %ylim([265 285])
% title('Fixed empirical times series mean (red) vs empirical by year (blue) pressure coefficient, 250 m')

% subplot(1,2,2)
% for i = 1:6
%     yr = 2014 + i;
%     yrindmld = find(mldyr == yr);
%     [mldmaxyr, mldmaxind] = max(chlmld(yrindmld));
%     datestr(chldt_mat(yrindmld(mldmaxind)));
%     [d, t_ind] = min(abs(wggmerge.time - chldt_mat(yrindmld(mldmaxind)-pad)));
%     h(i) = plot(movmean(wggmerge_fixedD.doxy(:,t_ind),20), pres_grid_hypm,'-','linewidth',2,'color',cmap(i,:)); hold on;
%     if i == 5
%         oxy_mldmaxstd(2,i) = nanstdev(movmean(wggmerge_fixedD.doxy(101:201,t_ind),20));
%     else
%         oxy_mldmaxstd(2,i) = nanstdev(movmean(wggmerge_fixedD.doxy(101:341,t_ind),20));
%     end
% end
% axis ij
% ylim([250 1300])
% grid on
% legend('2015', '2016', '2017', '2018', '2019', '2020','location','SW')
% title('WFP oxygen profile at MLDmax - Pc = 0.040 (overall mean)')
% xlabel('Oxygen, \mumol/kg')
% ylabel('Pressure, db')