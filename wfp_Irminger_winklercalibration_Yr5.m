
%% Plot locations of Year 5 casts that align with SUMO/HYPM
load OOImooringLocations.mat

ind = find(Winkler_casts(:,1) >=5 & Winkler_casts(:,1) <=6);
M = 10;
    figure(10); clf
plot(-OOImoorings.SUMO4(2), OOImoorings.SUMO4(1),'^k','markersize',M); hold on;
plot(-OOImoorings.HYPM4(2), OOImoorings.HYPM4(1),'^k','markersize',M); hold on;
plot(-OOImoorings.FLMA4(2), OOImoorings.FLMA4(1),'^k','markersize',M); hold on;
plot(-OOImoorings.FLMB4(2), OOImoorings.FLMB4(1),'^k','markersize',M); hold on;
plot(Winkler_casts(ind,4), Winkler_casts(ind,3),'r.','markersize',M*2); hold on;

%% Extract discrete Winkler values to use for initial gain correction
Yr5_disc.cast = Winkler_casts(:,1);
Yr5_disc.day = datenum(Winkler_text(2:end,4));
Yr5_disc.time = Winkler_casts(:,5);
Yr5_disc.lat = Winkler_casts(:,3);
Yr5_disc.lon360 = -1*Winkler_casts(:,4) + 360; %degrees W - make negative to match convention
Yr5_disc.depth = Winkler.depth; 
Yr5_disc.oxy = Winkler.O2_BCP;

%% Plot Year 5 discrete lined up with corresponding hypm profiles
C(1,:) = nicecolor('rrry'); C(2,:) = nicecolor('gggkb'); C(3,:) = nicecolor('bbc'); C(4,:) = nicecolor('rbm'); %color
time_window = 1; %all profiles within X days before or after the cast
depth_window = 20; %compare points in profile only if within X meters of Winkler sampling depth
tol_O2std = 5; %tolerance - only calculate gain if std of profile O2 measurements at a given depth are within X
tol_dist = 30; %limit of distance between cast and glider in order to use for calibration

k = 5;
    cast_list = [5,6]; %casts at SUMO and HYPM
    wfp_plot = Yr5_wfpgrid;
    disc_plot = Yr5_disc;
figure(k); clf
set(gcf,'color','w')
    x0=1;
    y0=1;
    width=28;
    height=20;
    set(gcf,'units','centimeters','position',[x0,y0,width,height]) 
G = NaN*ones(length(disc_plot.cast),1);
    for i = 1:length(cast_list)
        subplot(2,length(cast_list)/2,i)
        ind = find(disc_plot.cast == cast_list(i));
            A = disc_plot.time(ind(1)) + disc_plot.day(ind(1));
            ind_time = find(wfp_plot.time_start < A + time_window & wfp_plot.time_start > A - time_window);
            dist = distlatlon(wfp_plot.lat(ind_time(1)),disc_plot.lat(ind(end)),wfp_plot.lon(ind_time(1))+360,disc_plot.lon360(ind(end)));
        %plot(disc_plot.oxy(ind),disc_plot.depth(ind),'k*','color',C(1,:),'linewidth',4); hold on;
        h2 = plot(wfp_plot.O2conc(:,ind_time),wfp_plot.depth_grid,'k.'); hold on;
        set(gca,'YDir','reverse'); axis([240 315 0 2500]); ylabel('Depth (m)'); xlabel('O_2 concentration')
        title({['Cast ' num2str(cast_list(i)) ', Distance to mooring = ' num2str(dist,3) ' km']});
        for j = 1:length(ind) %loop over all discrete samples in cast
            [depth_dif,depth_ind] = min(abs(wfp_plot.depth_grid - disc_plot.depth(ind(j))));
            if depth_dif < depth_window %only use points within X meters of sampling depth
                O2raw = nanmean(wfp_plot.O2conc(depth_ind,ind_time));
                O2raw_std = nanstd(wfp_plot.O2conc(depth_ind,ind_time));
                if O2raw_std < tol_O2std & dist < tol_dist %don't calculate gain if the profiles at this depth/time are too variable
                    G(ind(j)) = disc_plot.oxy(ind(j))./O2raw;
                    h1 = plot(disc_plot.oxy(ind(j)),disc_plot.depth(ind(j)),'k*','color',C(1,:),'linewidth',4); hold on;
                end
            end
        end
    end   

%% Calculate gain correction
gain_hypm = nanmean(G);
gainstd_hypm = nanstd(G);
gainnum_hypm = sum(~isnan(G));

% Replot WFP data, now making gain correction
    for i = 1:length(cast_list)
        subplot(2,length(cast_list)/2,i)
        ind = find(disc_plot.cast == cast_list(i));
            A = disc_plot.time(ind(1)) + disc_plot.day(ind(1));
            ind_time = find(wfp_plot.time_start < A + time_window & wfp_plot.time_start > A - time_window);
            dist = distlatlon(wfp_plot.lat(ind_time(1)),disc_plot.lat(ind(end)),wfp_plot.lon(ind_time(1))+360,disc_plot.lon360(ind(end)));
        h3 = plot(wfp_plot.O2conc(:,ind_time)*nanmean(G),wfp_plot.depth_grid,'.','color',nicecolor('kw')); hold on; %gain correction
        %plot(wfp_plot.O2conc(:,ind_time)*(nanmean(G)-nanstd(G)),wfp_plot.depth_grid,'.','color',nicecolor('kww')); hold on; %std of gain lower
        %plot(wfp_plot.O2conc(:,ind_time)*(nanmean(G)+nanstd(G)),wfp_plot.depth_grid,'.','color',nicecolor('kww')); hold on; %std of gain higher
        %plot(disc_plot.oxy(ind),disc_plot.depth(ind),'k*','color',C(1,:),'linewidth',4); hold on; %original
        %legend([h1, h2, h3], 'Winkler measurements','Raw optode profile','Gain-corrected optode profile')
    end  
    