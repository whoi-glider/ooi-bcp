%%
%pfit = [-6.42714604659502e-07,0.000381938037536471,0.988215184315801];
pfit = [0.000243430208186408,0.993179458863815];
load('latest');
TG.oxygen_saturation_corr = polyval(pfit3,TG.daten-datenum(2018,6,1)).*TG.oxygen_saturation;
TG.oxygen_concentration_corr = polyval(pfit3,TG.daten-datenum(2018,6,1)).*TG.oxygen_concentration;

%%
ng = length(TG.zgrid);
TG.prof_daten = nanmean(TG.daten);
d = TG.prof_daten > datenum(2018,06,10) & TG.prof_daten < datenum(2018,10,10);
meanz = nanmean(TG.depth_interp(:,d),2);

rc = nan(ng,1);
rp = rc;
for ii = 1:ng
    dd = TG.prof_daten > datenum(2018,06,10) & TG.prof_daten < datenum(2018,10,10) & ~isnan(TG.oxygen_concentration(ii,:));
    pc = polyfit(TG.prof_daten(dd),TG.oxygen_concentration_corr(ii,dd),1);
    pp = polyfit(TG.prof_daten(dd),TG.oxygen_saturation_corr(ii,dd),1);
    rc(ii) = pc(1);
    rp(ii) = pp(1);
end

%% convert to mol m-3 y-1
r_season = 1000.*(datenum(2018,6,10)-datenum(2018,10,10)).*rc./1027;