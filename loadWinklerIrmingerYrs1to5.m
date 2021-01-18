%% Load discrete oxygen data from Irminger deployment cruises
%Year 1
    %For calculating precision, note that first 12 samples are 6 pairs of
    %duplicates (test samples)
[num,txt,~]=xlsread('IrmingerYr1_KN221_CTDWaterSamplingData2.xlsx');

disc{1}.cast = num(:,1);
disc{1}.day = num(:,2);
disc{1}.time = num(:,3);
disc{1}.lat = num(:,4);
disc{1}.lon360 = -1*num(:,5) + 360; %degrees W - make negative to match convention
disc{1}.depth = num(:,8); 
disc{1}.oxy = num(:,10)/(0.0223916); %convert from mL/L to mmol/L; O2 = 22391.6 ml/mol (Garcia and Gordon, 1992)

%Calculate uncertainty in Year 1 samples (take stdev of all 6 duplicates,
%and use mean stdev value)
for i = 1:6
    dups = disc{1}.oxy(2*i-1:2*i);
    duperr(i) = std(dups);
end
disc{1}.oxy_err = mean(duperr);
clear duperr

%Year 2
[num,txt,~]=xlsread('IrmingerYr2_AT30_CTDWaterSamplingData2.xlsx');

disc{2}.cast = num(:,1);
disc{2}.day = num(:,2);
disc{2}.time = num(:,3);
disc{2}.lat = num(:,4);
disc{2}.lon360 = -1*num(:,5) + 360; %degrees W - make negative to match convention
disc{2}.depth = num(:,8); 
disc{2}.oxy = num(:,10)/(0.0223916); %convert from mL/L to mmol/L; O2 = 22391.6 ml/mol (Garcia and Gordon, 1992)
disc{2}.nitrate = num(:,12); %uM
disc{2}.potT = num(:,13); %potential temp from CTD
disc{2}.S = num(:,14); %salinity from CTD

%Calculate uncertainty in Year 2 samples (take stdev of all duplicates,
%and use mean stdev value)
for i = 1:6
    dups = disc{2}.oxy(2*i-1:2*i);
    duperr(i) = std(dups);
end
for i = 18:22
    dups = disc{2}.oxy(2*i:2*i+1);
    duperr(i-11) = std(dups);
end
disc{2}.oxy_err = mean(duperr);
clear duperr

%Year 3
[num,txt,~]=xlsread('IrmingerYr3_AR07_CTDWaterSamplingData2.xlsx');

disc{3}.cast = num(:,1);
disc{3}.day = num(:,2);
disc{3}.time = num(:,3); %note that this is currently not formatted properly
disc{3}.lat = num(:,4);
disc{3}.lon360 = -1*num(:,5) + 360; %degrees W - make negative to match convention
disc{3}.depth = num(:,8); 
disc{3}.oxy = num(:,10)/(0.0223916); %convert from mL/L to mmol/L; O2 = 22391.6 ml/mol (Garcia and Gordon, 1992)
disc{3}.nitrate = num(:,14); %uM
disc{3}.potT = num(:,11); %potential temp from CTD
disc{3}.S = num(:,12); %salinity from CTD

%Calculate uncertainty in Year 3 samples (take stdev of all duplicates,
%and use mean stdev value)
for i = 1:10
    dups = disc{3}.oxy(2*i-1:2*i);
    duperr(i) = std(dups);
end
disc{3}.oxy_err = mean(duperr);
clear duperr

%Year 4
[num,txt,~]=xlsread('IrmingerYr4_AR21_CTDWaterSamplingData2.xlsx');

disc{4}.cast = num(:,1);
disc{4}.day = num(:,2);
disc{4}.time = num(:,3); %note that this is currently not formatted properly
disc{4}.lat = num(:,4);
disc{4}.lon360 = -1*num(:,5) + 360; %degrees W - make negative to match convention
disc{4}.depth = num(:,8); 
disc{4}.oxy = num(:,10)/(0.0223916); %convert from mL/L to mmol/L; O2 = 22391.6 ml/mol (Garcia and Gordon, 1992)

%Calculate uncertainty in Year 4 samples (take stdev of all duplicates,
%and use mean stdev value)
for i = 23:28 
    dups = disc{4}.oxy(2*i:2*i+1);
    duperr(i) = std(dups);
end
disc{4}.oxy_err = mean(duperr(23:28));
clear duperr

% Year 5
[Winkler_casts, Winkler_text] = xlsread('C:\Users\palevsky\Dropbox\Irminger5\Oxygen data\Irminger5_WinklerSamples.xlsx',2);
    Winkler.depth = Winkler_casts(:,8);
    Winkler.T = Winkler_casts(:,9); %potential temp
    Winkler.S = Winkler_casts(:,10); %practical salinity
    Winkler.O2_Dave = Winkler_casts(:,17); %Dave Wellwood values in umol/kg
    Winkler.O2_Dave_flag = Winkler_casts(:,18); %Quality flag (2 for different Niskin from BCP group)
    Winkler.O2_BCP = Winkler_casts(:,21:22); %Replicate samples from BCP group, umol/kg
    Winkler.O2_BCP_flag = Winkler_casts(:,23); %Quality flag

disc{5}.cast = Winkler_casts(:,1);
disc{5}.day = datenum(Winkler_text(2:end,4));
disc{5}.time = Winkler_casts(:,5);
disc{5}.lat = Winkler_casts(:,3);
disc{5}.lon360 = -1*Winkler_casts(:,4) + 360; %degrees W - make negative to match convention
disc{5}.depth = Winkler_casts(:,8);
disc{5}.oxy = Winkler_casts(:,21:22); %Replicate samples from BCP group, umol/kg
disc{5}.oxy_flag = Winkler_casts(:,23); %Quality flag (1 and 2 are good)