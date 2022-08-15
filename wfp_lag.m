%% Attempt to calculate and correct for lag based on determination of tau
% Rather than pairing up-down profiles as I have done previously
% Runs after 1st code block in wfp_analysis.m, which reads in all years of
% WFP O2 data and passes through load_HYPM_DOSTA_fun

%Make first attempt to do this just for Year 1 wfp data
testyr = 1;

%% Roo's version of correction, example in ws_O2tau.m for glider
% load('Tau4330_lookup.mat','T4330');
% [tau] = optode_tau(wfp{testyr}.time_dosta_mat, wfp{testyr}.pressure_dosta, wfp{testyr}.temperature_dosta, T4330);
% Note that there is an error that leads to all NaNs in this calc - haven't
% gone further down the path of troubleshooting, since this approach
% doesn't empirically determine the boundary layer thickness (uses values
% from Bittig & Kortzinger, 2017, section 3 - Figure 1 and eqns 1 & 2),
% which may not be applicable to WFP since determined for floats & gliders

%% Gordon et al. 2020 version of correction
% Format WFP output into matrices that match the required argument formats
% listed in calculate_tau.m and calculate_tau_wTemp.m

tol = 50; %only use profiles with at least 50 points

numprofiles = max(wfp{testyr}.profile_index);
    wgg.mtime = NaN(numprofiles,2000);
    wgg.pres = NaN(numprofiles,2000);
    wgg.doxy = NaN(numprofiles,2000);
    wgg.temp = NaN(numprofiles,2000);
    wgg.pdens = NaN(numprofiles,2000);
    wgg.updown = NaN(numprofiles,1);
    
for i = 1:numprofiles
    indp = find(wfp{testyr}.profile_index == i);
    if length(indp) > tol
        wgg.mtime(i,1:length(indp)) = wfp{testyr}.time_dosta_mat(indp);
        wgg.pres(i,1:length(indp)) = wfp{testyr}.pressure_dosta(indp);
        wgg.doxy(i,1:length(indp)) = wfp{testyr}.oxygen(indp);
        wgg.temp(i,1:length(indp)) = wfp{testyr}.temperature_dosta(indp);
        wgg.pdens(i,1:length(indp)) = wfp{testyr}.pdens(indp);
        wgg.updown(i) = wfp{testyr}.updown_index(indp(1+tol));
    end
    indzero = find(wgg.doxy(i,:) == 0);
    wgg.doxy(i,indzero) = NaN;
end

%% Cut out rows and columns of all NaNs
ind_good = find(~all(isnan(wgg.doxy')) & ~isnan(wgg.updown)' == 1); %includes > tol pts and has up/down index
    wgg.mtime = wgg.mtime(ind_good,~all(isnan(wgg.doxy)));
    wgg.pres = wgg.pres(ind_good,~all(isnan(wgg.doxy)));
    wgg.doxy = wgg.doxy(ind_good,~all(isnan(wgg.doxy)));
    wgg.temp = wgg.temp(ind_good,~all(isnan(wgg.doxy)));
    wgg.pdens = wgg.pdens(ind_good,~all(isnan(wgg.doxy)));
    wgg.updown = wgg.updown(ind_good);
    
%% Identify gaps in continuous paired profile sampling
    ttol = 1; %tolerance for gap between profiles, in days
tgap = find(diff(wgg.mtime(:,1)) > ttol); %Calculate time between profiles
pgap = find(abs(diff(wgg.updown)) < 2); %== 2 when switching from -1 to 1
gaps = union(tgap,pgap);
    
%% Calculate tau using Gordon et al. 2020 method
addpath(genpath('C:\Users\palevsky\Documents\GitHub\optode-response-time'))

rng = [1:69, 71:112, 120:183, 201:241, 243:273, 276:351, 353:393]; %use full range now that have identified which version to go forward with - hard coded based on gaps output
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
[wgg.thicknessd, wgg.tau_Trefd, thickness_constants, wgg.rmsdtd] = calculate_tau_wTemp( wgg.mtime(rng,:), wgg.pdens(rng,:), wgg.doxy(rng,:), wgg.temp(rng,:),...
    'zlim',zlim_dens, 'zres', zres, 'tref', tref);
toc
length(rng)

%% Assess temporal variability of thickness/tau

indgap = find(diff(rng) > 1);
indnotgap = find(diff(rng) == 1);
%Don't fill gaps because it prevents movmean plotting

figure(100); clf
    smthval = 30;
subplot(211)
plot(wgg.mtime(rng(1:end-1),1), wgg.thicknessd,'k.'); hold on;
plot(wgg.mtime(rng(indgap),1), wgg.thicknessd(indgap),'b.','markersize',10); hold on;
plot(wgg.mtime(rng(1:end-1),1), movmean(wgg.thicknessd,smthval),'r.'); hold on;
plot(wgg.mtime(rng(1:end-1),1), movmean(wgg.thicknessd,smthval) + movstd(wgg.thicknessd,smthval),'r-'); hold on;
plot(wgg.mtime(rng(1:end-1),1), movmean(wgg.thicknessd,smthval) - movstd(wgg.thicknessd,smthval),'r-'); hold on;
datetick('x')
ylabel('Thickness (\mum)')
legend('Individual points','Gap data','30-pt moving mean','+/- 30-pt moving stdev')
title('OOI Irminger WFP Year 1 lag, Gordon et al. 2020 T-dependent method calc. in density space')


subplot(212)
plot(wgg.mtime(rng(1:end-1),1), wgg.tau_Trefd,'k.'); hold on;
plot(wgg.mtime(rng(indgap),1), wgg.tau_Trefd(indgap),'b.','markersize',10); hold on;
plot(wgg.mtime(rng(1:end-1),1), movmean(wgg.tau_Trefd,smthval),'r.'); hold on;
plot(wgg.mtime(rng(1:end-1),1), movmean(wgg.tau_Trefd,smthval) + movstd(wgg.tau_Trefd,smthval),'r-'); hold on;
plot(wgg.mtime(rng(1:end-1),1), movmean(wgg.tau_Trefd,smthval) - movstd(wgg.tau_Trefd,smthval),'r-'); hold on;
datetick('x')
ylabel('\tau (s) at 4^oC')
legend('Individual points','Gap data','30-pt moving mean','+/- 30-pt moving stdev')

%Calculate annual stats on datapoints excluding the gaps
wgg.stats = [nanmean(wgg.thicknessd(indnotgap)) nanstd(wgg.thicknessd(indnotgap)) length(wgg.thicknessd(indnotgap));...
    nanmean(wgg.tau_Trefd(indnotgap)) nanstd(wgg.tau_Trefd(indnotgap)) length(wgg.tau_Trefd(indnotgap))];
wgg.stats(:,4) = wgg.stats(:,2)./sqrt(wgg.stats(:,3));


%% Apply tau correction using version w/ temperature term & density space

%Initialize array to hold corrected oxygen data
wgg.doxy_lagcorr = NaN(size(wgg.doxy)); %with temperature correction (Gordon et al. 2020 version w/ Bittig add on)
wgg.doxy_lagcorr2 = NaN(size(wgg.doxy)); %Roo's implementation, no temperature correction

for i = 1:404
    ind = ~(isnan(wgg.mtime(i,:)) | isnan(wgg.doxy(i,:)) | isnan(wgg.temp(i,:)));
    wgg.doxy_lagcorr(i,ind) = correct_oxygen_profile_wTemp(wgg.mtime(i,ind), wgg.doxy(i,ind), wgg.temp(i,ind), nanmean(wgg.thicknessd));
    wgg.doxy_lagcorr2(i,ind) = aa_deconv(wgg.mtime(i,ind), wgg.doxy(i,ind), nanmean(wgg.tau_Tref));
end

%% Plot examples of paired up & down profiles
smthval = 10;
M = 4;
figure(1); clf
for i = 1:4
    subplot(2,2,i)
    plot(wgg.doxy(i*10,:), wgg.pres(i*10,:),'k.'); hold on;
    plot(wgg.doxy(i*10 + 1,:), wgg.pres(i*10 + 1,:),'.'); hold on;
    plot(movmean(wgg.doxy_lagcorr(i*10,:),smthval), wgg.pres(i*10,:),'.','color',nicecolor('kww'),'markersize',M); hold on;
    plot(movmean(wgg.doxy_lagcorr(i*10 + 1,:),smthval), wgg.pres(i*10 + 1,:),'.','color',nicecolor('rrykkwwww'),'markersize',M); hold on;
    xlabel('Oxygen-L2, (\mumol/kg)')
    ylabel('dbar')
    title(datestr(wgg.mtime(i*10 + 1,1),1))
    legend('Up','Down','location','northwest')
    axis ij
end

%% Plot RMSD for each of the 4 approaches to calculating tau
% Tested for 1st 60 profiles in year 1
figure(2); clf
    subplot(221)
plot(wgg.time_constants,wgg.rmsd,'.')
xlabel('\tau (sec)')
ylabel('RMSD')
title('No T term, depth aligned')
    subplot(222)
plot(wgg.time_constants,wgg.rmsd,'.')
xlabel('\tau (sec)')
ylabel('RMSD')
title('No T term, density aligned')
    subplot(223)
plot(wgg.thickness_constants,wgg.rmsdt,'.')
xlabel('thickness (\mum)')
ylabel('RMSD')
title('w/ T term, depth aligned')
    subplot(224)
plot(wgg.thickness_constants,wgg.rmsdtd,'.')
xlabel('thickness (\mum)')
ylabel('RMSD')
title('w/ T term, density aligned')

%% Calculate WFP velocity in dbar s-1 to compare with Bittig & Kortzinger 2017
v = diff(wfp{testyr}.pressure_dosta)./(diff(wfp{testyr}.time_dosta_mat)*(24*60*60));
figure(3); clf
histogram(abs(v))
title({'WFP Yr 1 velocity histogram, median = ' num2str(nanmedian(abs(v)),3)})
xlabel('|dbar s^{-1}|')

%% NEXT STEPS
% Try actually applying the correction (correct_oxygen_profile.m)
% Expand to include the full year - probably just do for v. w/ dens & T