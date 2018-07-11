%% Read in underway Winkler samples value times and corresponding optode values

[Winkler_underway, Winkler_underwayTxt] = xlsread('C:/Users/Hilary/Dropbox/Irminger5/Irminger5_WinklerSamples.xlsx',3);
Winkler_underwayTxt = Winkler_underwayTxt(2:53,:);
Winkler_underway = Winkler_underway(1:52,:);
Winkler_timeUnderway(:,1) = datenum(Winkler_underwayTxt(:,3)) + Winkler_underway(:,8);
Winkler_timeUnderway(:,2) = datenum(Winkler_underwayTxt(:,3)) + Winkler_underway(:,9);

optode_UnderwayPts(:,1) = Winkler_underway(:,14);
optode_UnderwayPts(:,2) = Winkler_underway(:,15);

% Used once to extract underway optode values at times corresponding to
% sample collection --> now extract back from Excel file
% optode_UnderwayPts = NaN*ones(length(Winkler_timeUnderway),2);
% timepad = 1/60/24;
% for i = 1:length(Winkler_timeUnderway)
%     %find all optode times from timepad minutes before and after sampling
%     ind_optode = find(optode.time > (Winkler_timeUnderway(i,1) - timepad) & optode.time < (Winkler_timeUnderway(i,2) + timepad));
%     if length(ind_optode > 1)
%         optode_UnderwayPts(i,1) = nanmean(optode.O2_nospike_salcorr(ind_optode));
%         optode_UnderwayPts(i,2) = nanstd(optode.O2_nospike_salcorr(ind_optode));
%     end
% end

%% Plot gain correction over full cruise

%Pull out good Winkler values
Winkler_means = mean(Winkler_underway(:,10:11),2);
ind_goodWinkler = find(Winkler_underway(:,12) < 3 & isnan(Winkler_means) == 0);

%Approximate times of possible biofouling in underway system to exclude
%%% Between time of major chl spike and time we noticed biofouling and
%%% cleaned system
ind_questionable = find(Winkler_timeUnderway(:,1) > datenum(2018,6,17,12,0,0) & Winkler_timeUnderway(:,1) < datenum(2018,6,20));

%Find late cruise values when standards ran higher
%%% Fact that these do not look like there is any offset suggests that
%%% Winkler analysis batches with offset standards do not require different
%%% thiosulfate value
ind_endcruise = find(Winkler_timeUnderway(:,1) > datenum(2018,6,20));

%Points to keep
ind_fit = setdiff(ind_goodWinkler, ind_questionable);

%Fit a line to the good data
P = polyfit(Winkler_means(ind_fit), optode_UnderwayPts(ind_fit,1), 1);
[rho,df,rho_sig95] = correlate(Winkler_means(ind_fit), optode_UnderwayPts(ind_fit,1));

figure(1); clf
plot([300:350], P(1)*[300:350] + P(2), 'k--'); hold on;
plot(Winkler_means(ind_goodWinkler), optode_UnderwayPts(ind_goodWinkler,1),'r.','markersize', 15); hold on;
plot(Winkler_means(intersect(ind_goodWinkler,ind_endcruise)), optode_UnderwayPts(intersect(ind_goodWinkler,ind_endcruise),1),'c.','markersize', 12); hold on;
plot(Winkler_means(intersect(ind_questionable,ind_goodWinkler)), optode_UnderwayPts(intersect(ind_questionable,ind_goodWinkler),1),'b.','markersize',12); hold on;
xlabel('Winkler O_2 (\mumol/kg)'); ylabel('Optode O_2 (\mumol/L)'); title('Calibration of underway optode')
