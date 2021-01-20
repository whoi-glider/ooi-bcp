%Make test plots of wfp oxygen gain-corrected data

%Create a merged data set
wfpmerge.time = wfpgrid{1}.time_start(wfpgrid{1}.ind_pair);
wfpmerge.pdens = wfpgrid{1}.pdens;
wfpmerge.T = wfpgrid{1}.T;
wfpmerge.S = wfpgrid{1}.S;
wfpmerge.oxygen_gaincorr = wfpgrid{1}.oxygen_gaincorr;
for i = 2:5
    wfpmerge.time = [wfpmerge.time; wfpgrid{i}.time_start(wfpgrid{i}.ind_pair)];
    wfpmerge.pdens = [wfpmerge.pdens wfpgrid{i}.pdens];
    wfpmerge.T = [wfpmerge.T wfpgrid{i}.T];
    wfpmerge.S = [wfpmerge.S wfpgrid{i}.S];
    wfpmerge.oxygen_gaincorr = [wfpmerge.oxygen_gaincorr wfpgrid{i}.oxygen_gaincorr];
end


%Adjustable parameters for plotting
    mindepth = 150; maxdepth = 2600;
    cints = 60; %number of contour intervals
    C = cmocean('Dense'); %set colormap
    C2 = cmocean('Algae'); 

%Make plotting grid
[X,Y] = meshgrid(wfpmerge.time, depth_grid);

figure(1); clf;
    subplot(311) %Density
cmin = 27.5; cmax = 27.8; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.pdens - 1000,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('\sigma_\theta', 'Fontsize', 12)

    subplot(312) %Temperature
cmin = 2; cmax = 6; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.T,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Temperature (deg C)', 'Fontsize', 12)

    subplot(313) %Oxygen_corr concentration
cmin = 260; cmax = 320; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.oxygen_gaincorr,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Oxygen concentration (\mumol/kg)', 'Fontsize', 12)