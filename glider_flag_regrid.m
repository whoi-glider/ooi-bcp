function [GLout,fract_flag] = glider_flag_regrid(glg, oxymin, oxymax, oxyspike, tempmin, tempmax, tempspike, salmin, salmax, salspike, pres_grid, smval)

%% Flag spikes and outliers, and regrid data excluding those points
c = 0; %counter to keep track of # datapoints flagged

%Number of profiles
num_profiles = length(glg.profilelist);

%Create array to hold flag
glg.flag = 0*ones(size(glg.mtime));

for i = 1:num_profiles
    %Identify spikes
    O = diff(glg.doxy_lagcorr(i,:));
    Oout = find(abs(O) > oxyspike);
    T = diff(glg.temp(i,:));
    Tout = find(abs(T) > tempspike);
    S = diff(glg.sal(i,:));
    Sout = find(abs(S) > salspike);
    %Identify outliers
    O_outrange = find(glg.doxy_lagcorr(i,:) < oxymin | glg.doxy_lagcorr(i,:) > oxymax);
    T_outrange = find(glg.temp(i,:) < tempmin | glg.temp(i,:) > tempmax);
    S_outrange = find(glg.sal(i,:) < salmin | glg.sal(i,:) > salmax);
    gps_outrange = find(glg.lon(i,:) > 360 | glg.lat(i,:) > 90);
    %Flag spikes and outliers
    glg.flag(i,Oout) = glg.flag(i,Oout) + 1;
    glg.flag(i,Tout) = glg.flag(i,Tout) + 10;
    glg.flag(i,Sout) = glg.flag(i,Sout) + 100;
    glg.flag(i,O_outrange) = glg.flag(i,O_outrange) + 1000;
    glg.flag(i,T_outrange) = glg.flag(i,T_outrange) + 10000;
    glg.flag(i,S_outrange) = glg.flag(i,S_outrange) + 10000;
    glg.lon(i,gps_outrange) = NaN; glg.lat(i,gps_outrange) = NaN;
    %Keep track out # datapoints flagged
    c = c + length(find(glg.flag(i,:) > 0));
end

%Calculate the fraction of total oxygen data flagged as spikes out outliers
fract_flag = c./sum(~isnan(glg.doxy_lagcorr(:)));

%% Calculate density in raw profiles prior to gridding
[glg.SA, ~] = gsw_SA_from_SP(glg.sal, glg.pres, glg.lon, glg.lat); %absolute salinity from practical salinity - [SA, ~] = gsw_SA_from_SP(SP,p,long,lat)
glg.CT = gsw_CT_from_t(glg.SA, glg.temp, glg.pres); %Conservative Temperature from in-situ temperature - CT = gsw_CT_from_t(SA,t,p)
glg.pdens = gsw_rho(glg.SA, glg.CT, 0); %calculate potential density at reference pressure of 0 (surface)

%% Create arrays to hold gridded output
glg.time_start = NaN*ones(num_profiles,1);
glg.duration = NaN*ones(num_profiles,1);
glg.lat_profile = NaN*ones(num_profiles,1);
glg.lon_profile = NaN*ones(num_profiles,1);
glg.doxy_lagcorr_grid = NaN*ones(length(pres_grid),num_profiles);
glg.SA_grid = NaN*ones(length(pres_grid),num_profiles);
glg.CT_grid = NaN*ones(length(pres_grid),num_profiles);
glg.sal_grid = NaN*ones(length(pres_grid),num_profiles);
glg.temp_grid = NaN*ones(length(pres_grid),num_profiles);
glg.pdens_grid = NaN*ones(length(pres_grid),num_profiles);
glg.chl_grid = NaN*ones(length(pres_grid),num_profiles);
glg.backscatter_grid = NaN*ones(length(pres_grid),num_profiles);

    for i = 1:num_profiles
        ind = find(~isnan(glg.pres(i,:)) & ~isnan(glg.doxy_lagcorr(i,:)) & glg.flag(i,:) == 0); %no nan values for depth or oxygen and no range or spike flags
        if length(ind) > 10
            [~,ind_u,~] = unique(glg.pres(i,ind),'stable');
            ind = ind(ind_u);
            if length(ind) > 10
                glg.doxy_lagcorr_grid(:,i) = movmean(interp1(glg.pres(i,ind), glg.doxy_lagcorr(i,ind), pres_grid),smval);
                glg.SA_grid(:,i) = movmean(interp1(glg.pres(i,ind), glg.SA(i,ind), pres_grid),smval);
                glg.CT_grid(:,i) = movmean(interp1(glg.pres(i,ind), glg.CT(i,ind), pres_grid),smval);
                glg.sal_grid(:,i) = movmean(interp1(glg.pres(i,ind), glg.sal(i,ind), pres_grid),smval);
                glg.temp_grid(:,i) = movmean(interp1(glg.pres(i,ind), glg.temp(i,ind), pres_grid),smval);
                glg.pdens_grid(:,i) = movmean(interp1(glg.pres(i,ind), glg.pdens(i,ind), pres_grid),smval);
                glg.chl_grid(:,i) = movmean(interp1(glg.pres(i,ind), glg.chl(i,ind), pres_grid),smval);
                glg.backscatter_grid(:,i) = movmean(interp1(glg.pres(i,ind), glg.backscatter(i,ind), pres_grid),smval);
                glg.time_start(i) = nanmin(glg.mtime(i,:));
                glg.duration(i) = nanmax(glg.mtime(i,:)) - nanmin(glg.mtime(i,:));
                glg.lat_profile(i) = nanmean(glg.lat(i,:));
                glg.lon_profile(i) = nanmean(glg.lon(i,:)); 
            end
        end
    end

GLout = glg;
