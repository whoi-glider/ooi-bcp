%% Create a merged glider data product for analysis
%All glider data is in glgmerge - extract gridded data from longest record in each year to merge for analysis

%Overarching data
glgstart = 11; %selected glider for year 2 - PG528
yr = 2;
    ind = find(isnan(glgmerge{glgstart}.time_start) == 0);
glidermerge.time = glgmerge{glgstart}.time_start(ind);
glidermerge.duration = glgmerge{glgstart}.duration(ind);
glidermerge.lat = glgmerge{glgstart}.lat_profile(ind);
glidermerge.lon = glgmerge{glgstart}.lon_profile(ind);
glidermerge.deploy_yr = yr*ones(length(ind),1);

%Data gridded on pressure surfaces
glidermerge.temp = glgmerge{glgstart}.temp_grid(:,ind);
glidermerge.pracsal = glgmerge{glgstart}.sal_grid(:,ind).*G_sal(glgstart,1);
glidermerge.SA = glgmerge{glgstart}.SA_grid(:,ind).*G_sal(glgstart,1);
glidermerge.CT = glgmerge{glgstart}.CT_grid(:,ind);
glidermerge.pdens = glgmerge{glgstart}.pdens_grid(:,ind);
glidermerge.doxy = glgmerge{glgstart}.doxy_lagcorr_grid(:,ind).*glgmerge{glgstart}.oxygain_deepisotherm_linear(ind)';
glidermerge.chla = glgmerge{glgstart}.chl_grid(:,ind);
glidermerge.backscatter = glgmerge{glgstart}.backscatter_grid(:,ind);

for glgid = [13 14 1 3 5 8] %set to pull longest record each year - in years 2, 7, and 8, the longest record is from a profiling glider (only to 200 db)
    yr = yr + 1; %increment the year
    ind = find(isnan(glgmerge{glgid}.time_start) == 0);
    glidermerge.time = [glidermerge.time; glgmerge{glgid}.time_start(ind)];
    glidermerge.duration = [glidermerge.duration; glgmerge{glgid}.duration(ind)];
    glidermerge.lat = [glidermerge.lat; glgmerge{glgid}.lat_profile(ind)];
    glidermerge.lon = [glidermerge.lon; glgmerge{glgid}.lon_profile(ind)];
    glidermerge.deploy_yr = [glidermerge.deploy_yr; yr*ones(length(ind),1)];
    glidermerge.temp = [glidermerge.temp glgmerge{glgid}.temp_grid(:,ind)];
    glidermerge.pracsal = [glidermerge.pracsal glgmerge{glgid}.sal_grid(:,ind).*G_sal(glgid,1)];
    glidermerge.SA = [glidermerge.SA glgmerge{glgid}.SA_grid(:,ind).*G_sal(glgid,1)];
    glidermerge.CT = [glidermerge.CT glgmerge{glgid}.CT_grid(:,ind)];
    glidermerge.pdens = [glidermerge.pdens glgmerge{glgid}.pdens_grid(:,ind)];
    glidermerge.doxy = [glidermerge.doxy glgmerge{glgid}.doxy_lagcorr_grid(:,ind).*glgmerge{glgid}.oxygain_deepisotherm_linear(ind)'];
    glidermerge.chla = [glidermerge.chla glgmerge{glgid}.chl_grid(:,ind)];
    glidermerge.backscatter = [glidermerge.backscatter glgmerge{glgid}.backscatter_grid(:,ind)];
end
