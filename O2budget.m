load('latest.mat','G363');
% G363
%pfit = [0.000321169430436814,-235.809687886004];
pfit3 = [0.000320654357886349,0.960190376123461];
%pfit = [-1.19699085661010e-06,0.000620540002370935,0.948746118346959];
p3 = [8.56281474750741e-07,-0.000471952156923902,1.03836296648514];
G363.oxygen_saturation = polyval(pfit3,G363.daten-datenum(2018,6,1)).*G363.oxygen_saturation;zgrid = 1:2:1000;

d = G363.oxygen_saturation == 0 | G363.oxygen_concentration == 0;
G363.oxygen_saturation(d) = NaN;
G363.oxygen_concentration(d) = NaN;
[TG] = ws_grid(G363,zgrid,'UoD',-1);

%%
TG.salinity(TG.salinity < 30) = NaN;
O2eq = gasmoleq(TG.salinity,TG.temperature,'O2');
O2 = O2eq.*TG.oxygen_saturation./100;

dens = TG.density;
dens = fillsurf(dens);
mld = calcmld(dens,zgrid,0.05);
O2 = fillsurf(O2);
O2eq = fillsurf(O2eq);
S = fillsurf(TG.salinity);
T = fillsurf(TG.temperature);

O2surf = nanmean(O2(4:5,:));
O2eq_surf = nanmean(O2eq(4:5,:));

% mol O2 m-2
O2_int100 = 100.*mean(O2(1:50,:));
O2_int200 = 200.*mean(O2(1:100,:));
O2_int1000 = 1000.*mean(O2);


mbar2atm = 1013.25;
sec2day = 60*60*24;
%% get wind data
ws = ncread('dp0219.nc','met_relwind_speed');
slp = ncread('dp0219.nc','barometric_pressure')./mbar2atm;
dn = datenum(1900,1,1,0,0,ncread('dp0219.nc','time'));
d = ~isnan(O2surf);
TG.prof_daten = nanmean(TG.daten);
o2ge = interp1(TG.prof_daten(d),O2surf(d),dn);
%%
d = ~isnan(S(1,:));
Ssurf = interp1(TG.prof_daten(d),mean(S(1:4,d)),dn);
Tsurf = interp1(TG.prof_daten(d),mean(T(1:4,d)),dn);
[Fd, Fc, Fp, Deq, k] = fas(o2ge,ws,Ssurf,Tsurf,slp,'O2','L13');

[~,first] = find(~isnan(o2ge),1);


Ftot = sec2day.*(Fd+Fc+Fp);
d = ~isnan(Ftot);
cumflux = cumtrapz(dn(d),Ftot(d));

cumas = -interp1(dn(d),cumflux,TG.prof_daten);


function V = fillsurf(V)
    [nz,nt] = size(V);
    for ii = 1:nt
        for zz = 10:-1:1
            if isnan(V(zz,ii))
                V(zz,ii) = V(zz+1,ii);
            end
        end
    end
end
