%% Plot synthesis of glider, WFP, and fixed depth assets

%Load Kristen's MLD calculated from WFP chl data
addpath('G:\Shared drives\NSF_Irminger\Data_Files\From_Kristen')
load MLD_CHL_WFP.mat
load MLD_CHL_WFP7.mat
chldt_mat = [chldt; chldt7];
chlmld = [chlmld; chlmld7];

% Load the fixed asset oxygen product for yrs 5-8 produced by Kristen
addpath('G:\Shared drives\NSF_Irminger\OOI_DO_fixed_depth\Data\mixed_layer')
load mixed_layer_calibrated_oxygen.mat

%% Calculate O2 saturation for plotting

%Calculate for merged glider output
glidermerge.O2sat = gsw_O2sol(glidermerge.SA, glidermerge.CT,...
    repmat(pres_grid_glider', 1, length(glidermerge.lat)), repmat(glidermerge.lon', length(pres_grid_glider), 1),...
    repmat(glidermerge.lat', length(pres_grid_glider), 1));

%Calculate for wfp output
wggmerge.O2sat = gsw_O2sol(wggmerge.SA, wggmerge.CT, repmat(pres_grid_hypm', 1, length(wggmerge.lat)),...
    repmat(wggmerge.lon', length(pres_grid_hypm), 1), repmat(wggmerge.lat', length(pres_grid_hypm), 1));

%% Scatter plot with WFP and glider oxygen data
sz = 1;
ymax = 2000;
cmin = 265; cmax = 320;
skipval = 5;

figure(1); clf %Oxygen concentration
C = cmocean('Dense'); %set colormap
doxy_scat_wfp = wggmerge.doxy(:,1:skipval:end);
[Xw,Yw] = meshgrid(wggmerge.time(1:skipval:end), pres_grid_hypm);
doxy_scat_gl = glidermerge.doxy(:,1:skipval:end);
[Xg,Yg] = meshgrid(glidermerge.time(1:skipval:end), pres_grid_glider);

    subplot(311)
scatter(Xg(:),Yg(:),5,doxy_scat_gl(:),'filled'); hold on;
axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]);
colormap(C); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside'); caxis([cmin cmax]);
datetick('x',2,'keeplimits');
title('OOI Irminger glider oxygen concentration (\mumol/kg)', 'Fontsize', 12)

    subplot(312)
scatter(Xw(:),Yw(:),5,doxy_scat_wfp(:),'filled'); hold on;
axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]);
colormap(C); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside'); caxis([cmin cmax]);
datetick('x',2,'keeplimits');
title('OOI Irminger WFP oxygen concentration (\mumol/kg)', 'Fontsize', 12)

    subplot(313)
scatter(Xg(:),Yg(:),5,doxy_scat_gl(:),'filled'); hold on;
scatter(Xw(:),Yw(:),5,doxy_scat_wfp(:),'filled'); hold on;
scatter(ML_DO.DOdn, 5*ones(length(ML_DO.DOdn),1), 5, ML_DO.DO_umolkg_final,'filled');
plot(chldt_mat, chlmld, 'k.','markersize',5); hold on; %Kristen's MLD calcs from spring 2023
axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]);
colormap(C); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside'); caxis([cmin cmax]);
datetick('x',2,'keeplimits');
title('Merged OOI Irminger glider, WFP, and fixed depth oxygen concentration (\mumol/kg)', 'Fontsize', 12)

%% Merged glider and WFP records for temperature, salinity, O2sat, and chla
figure(2); clf %Temperature
C = cmocean('Thermal'); %set colormap
temp_scat_wfp = wggmerge.temp(:,1:skipval:end);
temp_scat_gl = glidermerge.temp(:,1:skipval:end);

scatter(Xg(:),Yg(:),5,temp_scat_gl(:),'filled'); hold on;
scatter(Xw(:),Yw(:),5,temp_scat_wfp(:),'filled'); hold on;
scatter(ML_DO.DOdn, 5*ones(length(ML_DO.DOdn),1), 5, ML_DO.sea_water_temperature,'filled');
plot(chldt_mat, chlmld, 'k.','markersize',5); hold on; %Kristen's MLD calcs from spring 2023

axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]);
colormap(C); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside'); caxis([1.5 10.5]);
datetick('x',2,'keeplimits');
title('Merged OOI Irminger glider, WFP, and fixed depth potential temperature (^oC)', 'Fontsize', 12)
%%
figure(3); clf %Salinity
C = cmocean('Haline'); %set colormap
sal_scat_wfp = wggmerge.pracsal(:,1:skipval:end);
sal_scat_gl = glidermerge.pracsal(:,1:skipval:end);

scatter(Xg(:),Yg(:),5,sal_scat_gl(:),'filled'); hold on;
scatter(Xw(:),Yw(:),5,sal_scat_wfp(:),'filled'); hold on;
scatter(ML_DO.DOdn, 5*ones(length(ML_DO.DOdn),1), 5, ML_DO.sea_water_practical_salinity,'filled');
plot(chldt_mat, chlmld, 'k.','markersize',5); hold on; %Kristen's MLD calcs from spring 2023

axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]);
colormap(C); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside'); caxis([34.6 35.05])
datetick('x',2,'keeplimits');
title('Merged OOI Irminger glider, WFP, and fixed depth practical salinity', 'Fontsize', 12)

%% 
figure(4); clf %Chlorophyll
C = cmocean('Algae'); %set colormap
chl_scat_wfp = log10(wggmerge_fl.chla(:,1:skipval:end));
[Xw_fl,Yw_fl] = meshgrid(wggmerge_fl.time(1:skipval:end), pres_grid_hypm);
chl_scat_gl = log10(glidermerge.chla(:,1:skipval:end));

scatter(Xg(:),Yg(:),5,chl_scat_gl(:),'filled'); hold on;
scatter(Xw_fl(:),Yw_fl(:),5,chl_scat_wfp(:),'filled'); hold on;
plot(chldt_mat, chlmld, 'k.','markersize',5); hold on; %Kristen's MLD calcs from spring 2023

axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]);
colormap(C); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside'); caxis([-3 0.5])
datetick('x',2,'keeplimits');
title('Merged OOI Irminger glider and WFP log10 chlorophyll', 'Fontsize', 12)

%%
figure(5); clf %Oxygen saturation
O2sat_scat_wfp = (doxy_scat_wfp./wggmerge.O2sat(:,1:skipval:end)-1)*100;
O2sat_scat_gl = (doxy_scat_gl./glidermerge.O2sat(:,1:skipval:end)-1)*100;

scatter(Xg(:),Yg(:),5,O2sat_scat_gl(:),'filled'); hold on;
scatter(Xw(:),Yw(:),5,O2sat_scat_wfp(:),'filled'); hold on;
scatter(ML_DO.DOdn, 5*ones(length(ML_DO.DOdn),1), 5, (ML_DO.DO_umolkg_final./ML_DO.O2sol_umolkg - 1)*100,'filled');
plot(chldt_mat, chlmld, 'k.','markersize',5); hold on; %Kristen's MLD calcs from spring 2023

axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]);
colormap(cmocean('Balance','pivot',0)); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Merged OOI Irminger glider, WFP, and fixed depth oxygen saturation (%)', 'Fontsize', 12)

%% 
figure(6); clf %Backscatter spikes - just wfp for now
spike_scat = wggmerge_fl.spikes(:,1:skipval:end);
scatter(Xw_fl(:),Yw_fl(:),5,spike_scat(:),'filled'); hold on;
plot(chldt_mat, chlmld, 'k.','markersize',8); hold on; %Kristen's MLD calcs from spring 2023

axis ij; axis tight; xlim([datenum(2014,9,10) datenum(2022,1,1)]); ylim([0 ymax]); caxis([0 7E-5]);
colormap(cmocean('Algae')); ylabel('Pressure (db)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('OOI Irminger WFP backscatter spikes (large particles)', 'Fontsize', 12)

%% Plot Meg's DIC figure
addpath 'C:\Users\palevsky\Documents\GitHub\CANYON-B-MFY\For Hilary review'
load workspace20230321.mat %Meg's workspace from analysis
addpath(genpath('C:\Users\palevsky\Documents\GitHub\boundedline-pkg'))

colorlist = [nicecolor('rrry'); nicecolor('ryyym'); nicecolor('yyyyyym');  nicecolor('ggyc'); nicecolor('ggbb'); nicecolor('bccc'); nicecolor('ccbb'); nicecolor('bbbbb'); nicecolor('m'); nicecolor('rrmmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb'); nicecolor('mmbb')];
red = nicecolor('rrryb'); cruise1 = nicecolor('ggycc'); cruise2 = nicecolor('gbbb'); 
  figure(7); clf; hold on 
for yr = 2:8
    castsumyr = castsum{yr};
    castlistyr = castlist{yr};  
            v = b_shift{yr};
            bl_er = 7.5*ones(length(deep{yr}.longdic), 1);
            bl = boundedline(deployment_dates{yr}, dic_mean{yr}.dic+v,dic_mean{yr}.dic_u+nanmean(bl_er), 'nan', 'gap', 'alpha'); 
            h1 = plot(deployment_dates{yr}, dic_mean{yr}.dic+v,'.','markersize',6,'color',nicecolor('bbckkkw'));         
end  
%formatting 
ylabel('DIC (\mumol/kg)')
xlim([datenum(2015,6,10) datenum(2022,8,1)]); 
datetick('x',2,'keeplimits');
title('OOI Irminger mixed layer DIC, Deployments 2-8','Fontsize',12)
ylim([2020 2200]) 
