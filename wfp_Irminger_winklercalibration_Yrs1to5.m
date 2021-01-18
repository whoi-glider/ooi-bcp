
%% Plot discrete data lined up with corresponding hypm profiles
C(1,:) = nicecolor('rrry'); C(2,:) = nicecolor('gggkb'); C(3,:) = nicecolor('bbc'); C(4,:) = nicecolor('rbm'); %color
time_window = 3; %all profiles within X days before or after the cast
depth_window = 20; %compare points in profile only if within X meters of Winkler sampling depth
tol_O2std = 5; %tolerance - only calculate gain if std of profile O2 measurements at a given depth are within X
tol_dist = 30; %limit of distance between cast and HYPM in order to use for calibration

for k = 1:5
    if k == 1
        cast_list = [5,6,7,9]; %casts in OOI site region (no test casts)
    elseif k == 2
        cast_list = [6,9,10,11,12,13]; %casts in OOI site region (no test casts, and not 4 and 5 b/c before gliders in water) %casts 11,12
    elseif k == 3
        cast_list = [7,8,9,10]; %casts in OOI site region - note that casts before #7 are before the Yr 3 profiler data has started
    elseif k == 4
        cast_list = [8,10,11,12]; %casts in OOI site region - note that casts before #8 are before the Yr 4 profiler data has started
    elseif k == 5
        cast_list = [5,6]; %casts at SUMO and HYPM
    end
    wfp_plot = wfpgrid{k};
    disc_plot = disc{k};

figure(k + 10); clf
set(gcf,'color','w')
    x0=1;
    y0=1;
    width=28;
    height=18;
    set(gcf,'units','centimeters','position',[x0,y0,width,height]) 
G = NaN*ones(length(disc_plot.cast),1);
    for i = 1:length(cast_list)
        subplot(2,length(cast_list)/2,i)
        ind = find(disc_plot.cast == cast_list(i));
            A = disc_plot.time(ind(1)) + disc_plot.day(ind(1));
            ind_time = find(wfp_plot.time_start < A + time_window & wfp_plot.time_start > A - time_window);
        if length(ind_time) >= 1
            dist = distlatlon(wfp_plot.lat(ind_time(1)),disc_plot.lat(ind(end)),wfp_plot.lon(ind_time(1))+360,disc_plot.lon360(ind(end)));
        h2 = plot(wfp_plot.O2conc(:,ind_time),wfp_plot.depth_grid,'k.'); hold on;
        %plot(disc_plot.oxy(ind),disc_plot.depth(ind),'k*','color',C(1,:),'linewidth',4); hold on;
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
                    %plot(disc_plot.oxy(ind(j)),disc_plot.depth(ind(j)),'bo','markersize',10,'linewidth',2); hold on;
                end
            end
        end
        end
    end   

%% Calculate gain correction
gain_hypm(k) = nanmean(G);
gainstd_hypm(k) = nanstd(G);
gainnum_hypm(k) = sum(~isnan(G));

% Replot WFP data, now making gain correction
    for i = 1:length(cast_list)
        subplot(2,length(cast_list)/2,i)
        ind = find(disc_plot.cast == cast_list(i));
            A = disc_plot.time(ind(1)) + disc_plot.day(ind(1));
            ind_time = find(wfp_plot.time_start < A + time_window & wfp_plot.time_start > A - time_window);
        if length(ind_time) >= 1
            dist = distlatlon(wfp_plot.lat(ind_time(1)),disc_plot.lat(ind(end)),wfp_plot.lon(ind_time(1))+360,disc_plot.lon360(ind(end)));
        h3 = plot(wfp_plot.O2conc(:,ind_time)*nanmean(G),wfp_plot.depth_grid,'.','color',nicecolor('kw')); hold on; %gain correction
        %plot(wfp_plot.O2conc(:,ind_time)*(nanmean(G)-nanstd(G)),wfp_plot.depth_grid,'.','color',nicecolor('kww')); hold on; %std of gain lower
        %plot(wfp_plot.O2conc(:,ind_time)*(nanmean(G)+nanstd(G)),wfp_plot.depth_grid,'.','color',nicecolor('kww')); hold on; %std of gain higher
        %plot(disc_plot.oxy(ind),disc_plot.depth(ind),'k*','color',C(1,:),'linewidth',4); hold on; %original
        %legend([h1, h2, h3], 'Winkler measurements','Raw optode profile','Gain-corrected optode profile')
        end
    end  
    
end
    