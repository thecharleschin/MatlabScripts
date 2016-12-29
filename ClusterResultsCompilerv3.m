%Script which compiles all results which are outputted from the cluster

%Number of runs
Runs = 100;
%Maximum length of a storage vector
MaxStore = 100000;
%Matix to contain all mRNA time data
AllmRNATrack = zeros(12001,Runs);
%Matrix to contain all protien time data
AllProteinTrack = zeros(12001,Runs);
%vector to track all first passage times of the two diffusing particles
FirstPassageTime = zeros(Runs,1);
%vector to track the number of times the two particles encounter each other
AllPassageCounts = zeros(Runs,1);
%matrix to track the duration of each encounter
AllPassageDurations = zeros(MaxStore,Runs);
%matrix to track the times at which each encounter occured
AllPassageTimesTrack = zeros(MaxStore,Runs);
%matrix to track the time between each encounter
AllPassageTimesBetween = zeros(MaxStore,Runs);
%matrix to track the number of mRNA produced in a given encounter
AllmRNAperBurst = zeros(MaxStore,Runs);
%matrix to track the times each mRNA is produced
AllmRNATimes = zeros(MaxStore,Runs);

%iterate through all runs
for ii = 1:Runs
    disp(ii)
    %load the next file
    %(must be changed to match files in folder)
    DynamicFilename = sprintf('16x16x16Crowding50per0006DResults%g',ii);
    load(DynamicFilename)
    
    
    %Store data into corresponding matrix
    AllPassageCounts(ii) = PassageCounterMatrix-1;
    PassageDurationTrack = PassageDurationTrack(:);
    AllPassageDurations(1:AllPassageCounts(ii),ii) = PassageDurationTrack(1:AllPassageCounts(ii));
    
    PassageTimesTrack = PassageTimesTrack(:);
    AllPassageTimesTrack(1:AllPassageCounts(ii),ii) = PassageTimesTrack(1:AllPassageCounts(ii));
    temp = diff(PassageTimesTrack(1:AllPassageCounts(ii)));
    AllPassageTimesBetween(1:length(temp),ii) = temp;
    AllmRNATrack(:,ii) = mRNATrack(:);
    AllProteinTrack(:,ii) = ProteinTrack(:);
    FirstPassageTime(ii) = FirstPassageMatrix;
    AllmRNATimes(:,ii) = mRNATimeTrack(1:MaxStore);
    AllmRNAperBurst(:,ii) = mRNAperBurstTrack(1:MaxStore);
    AllmRNAperBurst(AllPassageCounts(ii):end,ii) = NaN;
end

%remove 0 values from calculating means
AllmRNATimes(AllmRNATimes == 0) = NaN;
AllPassageDurations(AllPassageDurations == 0) = NaN;
AllPassageTimesTrack(AllPassageTimesTrack == 0) = NaN;
AllPassageTimesBetween(AllPassageTimesBetween == 0) = NaN;

AllmeanDurations = nanmean(AllPassageDurations,1);
AllmeanTimesBetween = nanmean(AllPassageTimesBetween,1);

BurstSize = AllmeanDurations;
BurstFreq = 1./AllmeanTimesBetween;

save CompiledResults16x16x16Crowding50per0006D