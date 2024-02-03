%% Grid all glider data (not just the single pick per year)

yrs = [2 2 2 3 4 5 5 6 6 7 7 8 8];
glidernum = [495 485 528 559 493 363 453 525 560 515 365 469 565];
i = 0;
for glgid = [9:11,13:14, 1:8]
    i = i + 1; %increment the index
    ind = find(isnan(glgmerge{glgid}.time_start) == 0);
    glidergrid{i}.glidernum = glidernum(i);
    glidergrid{i}.time = glgmerge{glgid}.time_start(ind);
    glidergrid{i}.duration = glgmerge{glgid}.duration(ind);
    glidergrid{i}.lat = glgmerge{glgid}.lat_profile(ind);
    glidergrid{i}.lon = glgmerge{glgid}.lon_profile(ind);
    glidergrid{i}.deploy_yr = yrs(i)*ones(length(ind),1);
    glidergrid{i}.temp = glgmerge{glgid}.temp_grid(:,ind);
    glidergrid{i}.pracsal = glgmerge{glgid}.sal_grid(:,ind);
    glidergrid{i}.SA = glgmerge{glgid}.SA_grid(:,ind);
    glidergrid{i}.CT = glgmerge{glgid}.CT_grid(:,ind);
    glidergrid{i}.pdens = glgmerge{glgid}.pdens_grid(:,ind);
    glidergrid{i}.doxy = glgmerge{glgid}.doxy_lagcorr_grid(:,ind).*glgmerge{glgid}.oxygain_deepisotherm_linear(ind)';
    glidergrid{i}.chla = glgmerge{glgid}.chl_grid(:,ind);
    glidergrid{i}.backscatter = glgmerge{glgid}.backscatter_grid(:,ind);
end
