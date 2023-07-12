% Creating to prep for URI talk, March 2023, using results so far

%Load Kristen's MLD calculated from WFP chl data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\From_Kristen')
load MLD_CHL_WFP.mat
load MLD_CHL_WFP7.mat
%chldt_mat = datenum(floor(chldt),0,0) + 365*(chldt - floor(chldt));
chldt_mat = [chldt; chldt7];
chlmld = [chlmld; chlmld7];

%% Calculate O2 saturation for plotting

%Calculate for glider
for i = 1:length(glgmerge)
    glgmerge{i}.O2sat = gsw_O2sol(glgmerge{i}.SA_grid, glgmerge{i}.CT_grid,...
        repmat(pres_grid_glider', 1, length(glgmerge{i}.lat_profile)), repmat(glgmerge{i}.lon_profile', length(pres_grid_glider), 1),...
        repmat(glgmerge{i}.lat_profile', length(pres_grid_glider), 1));
end

wggmerge.O2sat = gsw_O2sol(wggmerge.SA, wggmerge.CT, repmat(pres_grid_hypm', 1, length(wggmerge.lat)),...
    repmat(wggmerge.lon', length(pres_grid_hypm), 1), repmat(wggmerge.lat', length(pres_grid_hypm), 1));

%% Scatter plot with WFP and glider data
profilerng = [1:5:3371];
sz = 1;
ymax = 1800;

figure(1); clf %Oxygen concentration
C = cmocean('Dense'); %set colormap
for i = 1:length(glgmerge)
    glg = glgmerge{i};
    [X,Y] = meshgrid(glg.time_start(1:5:end), pres_grid_glider);
    GL_doxy_scat = glg.doxy_lagcorr_grid(:,1:5:end)./nanmean(glgmerge{i}.HYPMalign_stats.O2_presA_deepcor_mean);
    scatter(X(:),Y(:),5,GL_doxy_scat(:),'filled'); hold on;
end
doxy_scat = wggmerge.doxy_lagcorr(:,profilerng).*wfp_profilegain_interp(profilerng)';
[X,Y] = meshgrid(wggmerge.time(profilerng), pres_grid_hypm);
scatter(X(:),Y(:),5,doxy_scat(:),'filled'); hold on;

plot(chldt_mat, chlmld, 'k.','markersize',5); hold on; %Put in Kristen's new MLD data

axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]);
colormap(C); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Merged OOI Irminger glider and WFP oxygen concentration (\mumol/kg)', 'Fontsize', 12)
%%
figure(2); clf %Temperature
C = cmocean('Thermal'); %set colormap
for i = 1:length(glgmerge)
    glg = glgmerge{i};
    [X,Y] = meshgrid(glg.time_start(1:5:end), pres_grid_glider);
    GL_temp_scat = glg.temp_grid(:,1:5:end);
    scatter(X(:),Y(:),5,GL_temp_scat(:),'filled'); hold on;
end
temp_scat = wggmerge.temp(:,profilerng);
[X,Y] = meshgrid(wggmerge.time(profilerng), pres_grid_hypm);
scatter(X(:),Y(:),5,temp_scat(:),'filled'); hold on;

plot(chldt_mat, chlmld, 'k.','markersize',5); hold on; %Put in Kristen's new MLD data

axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]);
colormap(C); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Merged OOI Irminger glider and WFP potential temperature (^oC)', 'Fontsize', 12)

%% 
figure(3); clf %Oxygen saturation state
for i = 1:length(glgmerge)
    glg = glgmerge{i};
    [X,Y] = meshgrid(glg.time_start(1:5:end), pres_grid_glider);
    GL_doxy_scat = ((glg.doxy_lagcorr_grid(:,1:5:end)./nanmean(glgmerge{i}.HYPMalign_stats.O2_presA_deepcor_mean))./glg.O2sat(:,1:5:end) - 1)*100;
    scatter(X(:),Y(:),5,GL_doxy_scat(:),'filled'); hold on;
end
doxy_scat = ((wggmerge.doxy_lagcorr(:,profilerng).*wfp_profilegain_interp(profilerng)')./wggmerge.O2sat(:,profilerng) -1)*100;
[X,Y] = meshgrid(wggmerge.time(profilerng), pres_grid_hypm);
scatter(X(:),Y(:),5,doxy_scat(:),'filled'); hold on;

plot(chldt_mat, chlmld, 'k.','markersize',5); hold on; %Put in Kristen's new MLD data

axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]);
colormap(cmocean('Balance','pivot',0)); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Merged OOI Irminger glider and WFP oxygen saturation (%)', 'Fontsize', 12)

%% 
figure(4); clf %Chlorophyll - note that background at depth should be subtracted, and differs by deployment
for i = 1:length(glgmerge)
    glg = glgmerge{i};
    [X,Y] = meshgrid(glg.time_start, pres_grid_glider);
    GL_chl_scat = log10(glg.chl_grid);
    scatter(X(:),Y(:),5,GL_chl_scat(:),'filled'); hold on;
end
chl_scat = log10(wggmerge_fl.chla(:,1:3:end));
[X,Y] = meshgrid(wggmerge_fl.time(1:3:end), pres_grid_hypm);
scatter(X(:),Y(:),5,chl_scat(:),'filled'); hold on;

plot(chldt_mat, chlmld, 'k.','markersize',5); hold on; %Put in Kristen's new MLD data

axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]); caxis([-3 0]);
colormap(cmocean('Algae')); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('OOI Irminger WFP log10(chlorophyll a)', 'Fontsize', 12)

%% 
figure(5); clf %Backscatter spikes
spike_scat = wggmerge_fl.spikes(:,1:5:end);
[X,Y] = meshgrid(wggmerge_fl.time(1:5:end), pres_grid_hypm);
scatter(X(:),Y(:),5,spike_scat(:),'filled'); hold on;

plot(chldt_mat, chlmld, 'k.','markersize',8); hold on; %Put in Kristen's new MLD data

axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]); caxis([0 7E-5]);
colormap(cmocean('Algae')); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('OOI Irminger WFP backscatter spikes (large particles)', 'Fontsize', 12)

%% Plot Meg's DIC figure
addpath 'C:\Users\palevsky\Documents\GitHub\CANYON-B-MFY\For Hilary review'
load workspace20230321.mat %Meg's workspace from analysis
addpath(genpath('C:\Users\palevsky\Documents\GitHub\boundedline-pkg'))

colorlist = [nicecolor('rrry'); nicecolor('ryyym'); nicecolor('yyyyyym');  nicecolor('ggyc'); nicecolor('ggbb'); nicecolor('bccc'); nicecolor('ccbb'); nicecolor('bbbbb'); nicecolor('m'); nicecolor('rrmmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb')];
  clear leg
red = nicecolor('rrryb');
cruise1 = nicecolor('ggycc'); 
cruise2 = nicecolor('gbbb'); 
  figure(10); clf 
    
      hold on 
for yr = 2:8
    %xline(deployment_dates{yr}(1))
    castsumyr = castsum{yr};
    castlistyr = castlist{yr};
    
            v = b_shift{yr};

            bl_er = 7.5*ones(length(deep{yr}.longdic), 1);
%             if yr == 2 | yr == 5 | yr == 6 | yr == 7
%                 bl_er(deep{yr}.longdic == NaN) = NaN;
%                 bl2 = boundedline(deep{yr}.longdate, deep{yr}.longdic,bl_er, 'nan', 'gap', 'alpha'); 
%             end 
            bl = boundedline(deployment_dates{yr}, dic_mean{yr}.dic+v,dic_mean{yr}.dic_u+nanmean(bl_er), 'nan', 'gap', 'alpha'); 
            h1 = plot(deployment_dates{yr}, dic_mean{yr}.dic+v,'.','markersize',6,'color',nicecolor('bbckkkw')); 
            
            %h2 =  plot(deep{yr}.longdate, deep{yr}.longdic,'.','markersize',10,'color','k');

                %legend([h1 h2],...
                   %'Average 0-50m, adjusted', 'Average 80-130m', 'orientation','horizontal','location','south')
  
        %formatting 
        ylabel('DIC (\mumol/kg)')
        xlim([datenum(2015,6,10) datenum(2022,8,1)]); 
        datetick('x',2,'keeplimits');
        title('OOI Irminger mixed layer DIC, Deployments 2-8','Fontsize',12)
         ylim([2020 2200]) 
                
                
    end  

%% Next steps
% 1. Add temperature plot with glider & WFP - DONE
% 2. Calculate O2 saturation with glider & WFP - DONE
% 3. Ask Kristen for updated MLDs with chlorophyll - DONE
% 4. Make chl and backscatter plots, building on Jose's pipeline and my
% version of it - DONE
% 5. Make the version of Meg's plot
% Not for presentation, but look back at Yr 1-4 gliders, think I might be
% able to use