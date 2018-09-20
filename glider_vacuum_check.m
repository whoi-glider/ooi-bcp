%engineering data for glider 363 from Collin
A = xlsread('GI_363 Irminger-5-m_vacuum.csv');
G363eng.daten = datenum(1970,1,1,00,0,A(:,1));
G363eng.m_vac = A(:,2);

%find when air bladder was fully inflated
    threshold = 9;
ind_infl = find(G363eng.m_vac > threshold);

%identify individual fully-inflated intervals
    mintime = 2/24; %find indices for full inflation after a gap of >2 hrs (i.e. not part of last surfacing)
timegaps = diff(G363eng.daten(ind_infl));
ind = find(timegaps > (mintime));

%plot m_vac data data over time, highlighting fully inflated times
figure; clf
plot(G363eng.daten, G363eng.m_vac, 'k.'); hold on;
plot(G363eng.daten(ind_infl), G363eng.m_vac(ind_infl), 'b.'); hold on;
plot(G363eng.daten(ind_infl(ind + 1)), G363eng.m_vac(ind_infl(ind + 1)), 'r.'); hold on;
datetick('x')

%% combine with aircal results for GL363
    %run aircal with options for 363 selected to output table T
    tol = 0.5/24; %needs to have a fully inflated air value within half hour of the mean time of the surface interval
for i = 1:height(T)
    ind = find(abs(G363eng.daten(ind_infl) - T.air_daten(i)) < tol);
    if size(ind) > 0
        check_infl(i) = length(ind);
    else
        check_infl(i) = NaN;
    end
end

figure;
plot(T.met_o2sat,T.air_corr,'.','MarkerSize',8); hold on;
    ind_infl = find(isnan(check_infl) == 0);
    new_points = find(T.air_daten > datenum(2018,9,4));
    new_plot = intersect(ind_infl, new_points);
plot(T.met_o2sat(ind_infl),T.air_corr(ind_infl),'r.','MarkerSize',15); hold on;
plot(T.met_o2sat(new_plot),T.air_corr(new_plot),'m.','MarkerSize',18); hold on;
    ind = intersect(find(isnan(T.met_o2sat + T.air_corr) == 0), ind_infl);
    p2 = polyfit(T.met_o2sat(ind),T.air_corr(ind),1);
hold all;
box on;
ax3 = gca;
plot(ax3.XLim,p2(1).*ax3.XLim+p2(2),'-k','LineWidth',3);
ax3.LineWidth = 1;
ax3.FontSize = ftsz;
xlabel('\DeltaO_{2,w}^{met}');
ylabel('\DeltaO_{2,w}^{mcorr}');

[rho,df,rho_sig95] = correlate(T.met_o2sat(ind),T.air_corr(ind))
