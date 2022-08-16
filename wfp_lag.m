%% Attempt to calculate and correct for lag based on determination of tau
% Rather than pairing up-down profiles as I have done previously
% Runs after 1st code block in wfp_analysis.m, which reads in all years of
% WFP O2 data and passes through load_HYPM_DOSTA_fun

%% Loop over all years - note that this is very time intensive to output saved
% Can comment out this full loop section and instead load data:
load lagyr1to7.mat %has everything in current code version - note that years 6 and 7 only have beginning of year

for yr = 1:7

%% Gordon et al. 2020 version of correction
% Format WFP output into matrices that match the required argument formats
% listed in calculate_tau.m and calculate_tau_wTemp.m

tol = 50; %only use profiles with at least 50 points

numprofiles = max(wfp{yr}.profile_index);
    wgg{yr}.mtime = NaN(numprofiles,2000);
    wgg{yr}.pres = NaN(numprofiles,2000);
    wgg{yr}.doxy = NaN(numprofiles,2000);
    wgg{yr}.temp = NaN(numprofiles,2000);
    wgg{yr}.pdens = NaN(numprofiles,2000);
    wgg{yr}.updown = NaN(numprofiles,1);
for i = 1:numprofiles
    indp = find(wfp{yr}.profile_index == i);
    if length(indp) > tol
        wgg{yr}.mtime(i,1:length(indp)) = wfp{yr}.time_dosta_mat(indp);
        wgg{yr}.pres(i,1:length(indp)) = wfp{yr}.pressure_dosta(indp);
        wgg{yr}.doxy(i,1:length(indp)) = wfp{yr}.oxygen(indp);
        wgg{yr}.temp(i,1:length(indp)) = wfp{yr}.temperature_dosta(indp);
        wgg{yr}.pdens(i,1:length(indp)) = wfp{yr}.pdens(indp);
        wgg{yr}.updown(i) = wfp{yr}.updown_index(indp(1+tol));
    end
    indzero = find(wgg{yr}.doxy(i,:) == 0);
    wgg{yr}.doxy(i,indzero) = NaN;
end

%% Cut out rows and columns of all NaNs
ind_good = find(~all(isnan(wgg{yr}.doxy')) & ~isnan(wgg{yr}.updown)' == 1); %includes > tol pts and has up/down index
    wgg{yr}.mtime = wgg{yr}.mtime(ind_good,~all(isnan(wgg{yr}.doxy)));
    wgg{yr}.pres = wgg{yr}.pres(ind_good,~all(isnan(wgg{yr}.doxy)));
    wgg{yr}.doxy = wgg{yr}.doxy(ind_good,~all(isnan(wgg{yr}.doxy)));
    wgg{yr}.temp = wgg{yr}.temp(ind_good,~all(isnan(wgg{yr}.doxy)));
    wgg{yr}.pdens = wgg{yr}.pdens(ind_good,~all(isnan(wgg{yr}.doxy)));
    wgg{yr}.updown = wgg{yr}.updown(ind_good);
    
%% Identify gaps in continuous paired profile sampling
    ttol = 2; %tolerance for gap between profiles, in days
wgg{yr}.tgap = find(diff(wgg{yr}.mtime(:,1)) > ttol); %Calculate time between profiles
wgg{yr}.pgap = find(abs(diff(wgg{yr}.updown)) < 2); %== 2 when switching from -1 to 1
wgg{yr}.gaps = union(wgg{yr}.tgap,wgg{yr}.pgap);
    
%% Calculate tau using Gordon et al. 2020 method
addpath(genpath('C:\Users\palevsky\Documents\GitHub\optode-response-time'))

if yr == 1
    wgg{yr}.rng = [1:69, 71:112, 120:183, 201:241, 243:273, 276:351, 353:393];
elseif yr == 2
    wgg{yr}.rng = [1:89,91:382];
elseif yr == 3
    wgg{yr}.rng = [1:238, 242:388, 392:length(wgg{yr}.updown)];
elseif yr == 4
    wgg{yr}.rng = [1:185, 188:224, 256:268, 293:362, 375:length(wgg{yr}.updown)];
elseif yr == 6
    %rng = [1:length(wgg{yr}.updown)]; rng(wgg{yr}.gaps) = NaN; rng(240:279) = NaN; rng = rng(~isnan(rng));
    %rng = [1:57, 59:93, 102:188, 209:222, 224:272, 286:327];
    wgg{yr}.rng = [1:57];
elseif yr == 7
    wgg{yr}.rng = [1:79];
else
    wgg{yr}.rng = [1:length(wgg{yr}.updown)]; wgg{yr}.rng(wgg{yr}.gaps) = NaN; wgg{yr}.rng = wgg{yr}.rng(~isnan(wgg{yr}.rng));
end
tref = 4; %reference temperature for tau
zlim_depth = [300 2500];
zlim_dens = [1027.62 1027.92]; %density bounds over which to calculate tau (~+/- 2.5 std from mean in yr1)
zres = 0.0005;

% %Standard version, no temperature term
% [wgg.tau, wgg.time_constants, wgg.rmsd] = calculate_tau( wgg.mtime(rng,:), wgg.pres(rng,:), wgg.doxy(rng,:), 'zlim',zlim_depth);
% 
% %Standard version but in density space, no temperature term
% [wgg.taud, wgg.time_constants, wgg.rmsdd] = calculate_tau( wgg.mtime(rng,:), wgg.pdens(rng,:), wgg.doxy(rng,:),...
%     'zlim',zlim_dens, 'zres', zres);
% 
% %Version with temperature term, pressure space
% [wgg.thickness, wgg.tau_Tref, wgg.thickness_constants, wgg.rmsdt] = calculate_tau_wTemp( wgg.mtime(rng,:), wgg.pres(rng,:), wgg.doxy(rng,:), wgg.temp(rng,:),...
%     'zlim',zlim_depth, 'tref', tref);

%Version with temperature term, also in density space
tic
[wgg{yr}.thicknessd, wgg{yr}.tau_Trefd, wgg{yr}.thickness_constants, wgg{yr}.rmsdtd] = calculate_tau_wTemp( wgg{yr}.mtime(rng,:), wgg{yr}.pdens(rng,:), wgg{yr}.doxy(rng,:), wgg{yr}.temp(rng,:),...
    'zlim',zlim_dens, 'zres', zres, 'tref', tref);
toc
length(rng)

% Calculate annual stats on datapoints excluding the gaps
wgg{yr}.stats = [nanmean(wgg{yr}.thicknessd) nanstd(wgg{yr}.thicknessd) length(wgg{yr}.thicknessd);...
    nanmean(wgg{yr}.tau_Trefd) nanstd(wgg{yr}.tau_Trefd) length(wgg{yr}.tau_Trefd)];
wgg{yr}.stats(:,4) = wgg{yr}.stats(:,2)./sqrt(wgg{yr}.stats(:,3));


%% Apply tau correction using version w/ temperature term & density space
[m,n] = size(wgg{yr}.doxy);

%Initialize array to hold corrected oxygen data
wgg{yr}.doxy_lagcorr = NaN(size(wgg{yr}.doxy)); %with temperature correction (Gordon et al. 2020 version w/ Bittig add on)
wgg{yr}.doxy_lagcorr2 = NaN(size(wgg{yr}.doxy)); %Roo's implementation, no temperature correction

for i = 1:m
    ind = ~(isnan(wgg{yr}.mtime(i,:)) | isnan(wgg{yr}.doxy(i,:)) | isnan(wgg{yr}.temp(i,:)));
    wgg{yr}.doxy_lagcorr(i,ind) = correct_oxygen_profile_wTemp(wgg{yr}.mtime(i,ind), wgg{yr}.doxy(i,ind), wgg{yr}.temp(i,ind), nanmean(wgg{yr}.thicknessd));
    wgg{yr}.doxy_lagcorr2(i,ind) = aa_deconv(wgg{yr}.mtime(i,ind), wgg{yr}.doxy(i,ind), nanmean(wgg{yr}.tau_Trefd));
end

end

%% Assess temporal variability of thickness/tau

% indgap = find(diff(rng) > 1);
% indnotgap = find(diff(rng) == 1);
%Don't fill gaps because it prevents movmean plotting

figure(100); clf
    smthval = 60;
subplot(211)
for yr = 1:7
plot(wgg{yr}.mtime(wgg{yr}.rng(1:end-1),1), wgg{yr}.thicknessd,'k.'); hold on;
plot(wgg{yr}.mtime(wgg{yr}.rng(1:end-1),1), movmean(wgg{yr}.thicknessd,smthval),'r.'); hold on;
end
xlim([datenum(2014,8,15) datenum(2021,1,1)])
ylim([0 120])
datetick('x','keeplimits')
ylabel('Thickness (\mum)')
legend('Individual points','60-pt moving mean')
title(['OOI Irminger WFP lag, Gordon et al. 2020 T-dependent method calc. in density space'])

subplot(212)
for yr = 1:7
plot(wgg{yr}.mtime(wgg{yr}.rng(1:end-1),1), wgg{yr}.tau_Trefd,'k.'); hold on;
plot(wgg{yr}.mtime(wgg{yr}.rng(1:end-1),1), movmean(wgg{yr}.tau_Trefd,smthval),'r.'); hold on;
end
xlim([datenum(2014,8,15) datenum(2021,1,1)])
ylim([0 80])
datetick('x','keeplimits')
ylabel('\tau (s) at 4^oC')
legend('Individual points','60-pt moving mean')
title(['OOI Irminger WFP lag, Gordon et al. 2020 T-dependent method calc. in density space'])

%% Check lag-corr output

%Check for outliers/range test
figure(101); clf
for yr = 1:7
    A = wgg{yr}.doxy_lagcorr(:);
    histogram(A); hold on;
    r(yr,1) = nanmean(A) - 4*nanstd(A);
    r(yr,2) = nanmean(A) + 4*nanstd(A);
    r(yr,3) = length(find(A < r(yr,1)));
    r(yr,4) = length(find(A > r(yr,2)));
end
xlim([220 320])
legend('Year 1','Year 2','Year 3', 'Year 4', 'Year 5', 'Year 6', 'Year 7')
title('Histogram all OOI Irminger WFP L2-oxygen, lag corrected only')

%% Plot examples of paired up & down profiles
smthval = 1;
M = 4;
int = 80;
yr = 6;
figure(1); clf
for i = 1:4
    subplot(2,2,i)
    plot(wgg{yr}.doxy(i*int,:), wgg{yr}.pres(i*int,:),'k.'); hold on;
    plot(wgg{yr}.doxy(i*int + 1,:), wgg{yr}.pres(i*int + 1,:),'.'); hold on;
    plot(movmean(wgg{yr}.doxy_lagcorr(i*int,:),smthval), wgg{yr}.pres(i*int,:),'.','color',nicecolor('kww'),'markersize',M); hold on;
    plot(movmean(wgg{yr}.doxy_lagcorr(i*int + 1,:),smthval), wgg{yr}.pres(i*int + 1,:),'.','color',nicecolor('rrykkwwww'),'markersize',M); hold on;
    xlabel('Oxygen-L2, (\mumol/kg)')
    ylabel('dbar')
    title(datestr(wgg{yr}.mtime(i*int + 1,1),1))
    legend('Up','Down','location','northwest')
    axis ij
end

%% Plot RMSD for each of the 4 approaches to calculating tau
% Tested for 1st 60 profiles in year 1
figure(2); clf
%     subplot(221)
% plot(wgg.time_constants,wgg.rmsd,'.')
% xlabel('\tau (sec)')
% ylabel('RMSD')
% title('No T term, depth aligned')
%     subplot(222)
% plot(wgg.time_constants,wgg.rmsd,'.')
% xlabel('\tau (sec)')
% ylabel('RMSD')
% title('No T term, density aligned')
%     subplot(223)
% plot(wgg.thickness_constants,wgg.rmsdt,'.')
% xlabel('thickness (\mum)')
% ylabel('RMSD')
% title('w/ T term, depth aligned')
%     subplot(224)
plot(wgg{6}.thickness_constants,wgg{6}.rmsdtd,'.')
xlabel('thickness (\mum)')
ylabel('RMSD')
title('w/ T term, density aligned')

%% Calculate WFP velocity in dbar s-1 to compare with Bittig & Kortzinger 2017
% v = diff(wfp{testyr}.pressure_dosta)./(diff(wfp{testyr}.time_dosta_mat)*(24*60*60));
% figure(3); clf
% histogram(abs(v))
% title({'WFP Yr 1 velocity histogram, median = ' num2str(nanmedian(abs(v)),3)})
% xlabel('|dbar s^{-1}|')

