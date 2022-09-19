%% Read in corrected cruise oxygen data
% CTD SBE43 data are corrected by Kristen using Winklers

mltoumol = 44.661;

%% Read in processed data and reformat
addpath(genpath('G:\Shared drives\NSF_Irminger\OOI Cruises CTD Casts\CTD_Data\Alfresco'))
clear castlist

load Year1_Processed_KF.mat
castsum_yr1{2} = cast02;
castsum_yr1{3} = cast03;
castsum_yr1{4} = cast04;
castsum_yr1{5} = cast05;
castsum_yr1{6} = cast06;
castsum_yr1{7} = cast07;
castsum_yr1{8} = cast08;
castsum_yr1{9} = cast09;
castlist{1} = [2:9];

load Year2_Processed_KF.mat
castsum_yr2{1} = cast01;
castsum_yr2{2} = cast02;
castsum_yr2{3} = cast03;
castsum_yr2{4} = cast04;
castsum_yr2{5} = cast05;
castsum_yr2{6} = cast06;
castsum_yr2{7} = cast07;
castsum_yr2{8} = cast08;
castsum_yr2{9} = cast09;
castsum_yr2{10} = cast10;
castsum_yr2{11} = cast11;
castsum_yr2{12} = cast12;
castsum_yr2{13} = cast13;
castlist{2} = [1:13];

load Year3_Processed_KF.mat
castsum_yr3{1} = cast01;
castsum_yr3{2} = cast02;
castsum_yr3{3} = cast03;
castsum_yr3{4} = cast04;
castsum_yr3{5} = cast05;
castsum_yr3{6} = cast06;
castsum_yr3{7} = cast07;
castsum_yr3{8} = cast08;
castsum_yr3{9} = cast09;
castsum_yr3{10} = cast10;
castlist{3} = [1:10];

load Year4_Processed_KF.mat
castsum_yr4{9} = cast09;
castsum_yr4{10} = cast10;
castsum_yr4{11} = cast11;
castsum_yr4{12} = cast12;
castlist{4} = [9:12];

load Year5_Processed_KF.mat
castsum_yr5{1} = cast01;
castsum_yr5{2} = cast02;
castsum_yr5{3} = cast03;
castsum_yr5{4} = cast04;
castsum_yr5{5} = cast05;
castsum_yr5{6} = cast06;
castsum_yr5{7} = cast07;
castsum_yr5{8} = cast08;
castsum_yr5{9} = cast09;
castsum_yr5{10} = cast10;
castsum_yr5{11} = cast11;
castsum_yr5{12} = cast12;
castsum_yr5{13} = cast13;
castsum_yr5{14} = cast14;
castsum_yr5{15} = cast15;
castsum_yr5{16} = cast16;
castsum_yr5{17} = cast17;
castsum_yr5{18} = cast18;
castsum_yr5{19} = cast19;
castsum_yr5{20} = cast20;
castsum_yr5{21} = cast21;
castsum_yr5{22} = cast22;
castsum_yr5{23} = cast23;
castlist{5} = [1:16];

%Merge bottle summary tables
btlsum{1} = btlsum_yr1;
btlsum{2} = btlsum_yr2;
btlsum{3} = btlsum_yr3;
btlsum{4} = btlsum_yr4;
btlsum{5} = btlsum_yr5;

%Merge cast structures
castsum{1} = castsum_yr1;
castsum{2} = castsum_yr2;
castsum{3} = castsum_yr3;
castsum{4} = castsum_yr4;
castsum{5} = castsum_yr5;

%% Loop over data from all years
for yr = 1:5
    castsumyr = castsum{yr};
    castlistyr = castlist{yr};
    
%% Plot all casts and interpolate at Winkler depth
figure(100*yr); clf
for i = castlistyr
subplot(4,4,i)
    btlid = find(btlsum{yr}.Cast == i);
plot(castsumyr{i}.sbeox0mL_L*mltoumol, castsumyr{i}.depSM, 'k-'); hold on;
if length(btlid) > 0
    plot(btlsum{yr}.Winkler_mLL(btlid)*mltoumol, btlsum{yr}.CTD_Depth_m(btlid), 'ro'); hold on;
    btlsum{yr}.CTD_Oxygen_mLL_corr(btlid) = interp1(castsumyr{i}.depSM, castsumyr{i}.sbeox0mL_L, btlsum{yr}.CTD_Depth_m(btlid));
end
axis ij
xlabel('SBE43 calibrated O_2, \mumol/L')
ylabel('depSM - db? m?')
end

%% Check goodness of calibration

%Regression plot
figure(100*yr + 1); clf
scatter(btlsum{yr}.Winkler_mLL*mltoumol, btlsum{yr}.CTD_Oxygen_mLL_corr*mltoumol, [], btlsum{yr}.CTD_Depth_m, 'filled'); hold on;
plot([140:350],[140:350],'k--')
c = colorbar; ylabel(c, 'Depth (m)')
axis([140 345 140 345])
xlabel('Winkler oxygen, \mumol/L'); ylabel('Calibrated SBE43 oxygen, \mumol/L')

%Calculate uncertainty from residuals
U = nanmean(abs(btlsum{yr}.Winkler_mLL - btlsum{yr}.CTD_Oxygen_mLL_corr))*mltoumol;
    deep_cutoff = 200; %subset deep values only
    deep = find(btlsum{yr}.CTD_Depth_m > deep_cutoff);
U_deep = nanmean(abs(btlsum{yr}.Winkler_mLL(deep) - btlsum{yr}.CTD_Oxygen_mLL_corr(deep)))*mltoumol;  

%Plot residuals
figure(100*yr + 2); clf
plot(btlsum{yr}.CTD_Depth_m, (btlsum{yr}.Winkler_mLL - btlsum{yr}.CTD_Oxygen_mLL_corr)*mltoumol, 'k.'); hold on;
plot([0 3000], 0*[0 3000], 'k--')
ylabel('Residual, Winkler - SBE O_2, \mumol/L')
xlabel('Depth (m)')
text(500, 4, ['Mean residual = ' num2str(U) ' \mumol/L'])
text(500, 3.4, ['Mean residual > ' num2str(deep_cutoff) 'm = ' num2str(U_deep) ' \mumol/L'])

end
  
