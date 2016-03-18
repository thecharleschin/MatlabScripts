%Processing from Traces to noise data
%Outputs - noise data, autocorrelations, steady state values for all
%days, buckets, and ROIs.

%%%%%%
%NOTE%
%Additional Toolboxes may be required to run all functions. Known required
%toolboxes include:
%Statistics Toolbox


%Reset workspace state
clear all
close all
rand('state',sum(100*clock)); %#ok<RAND>


%%
%%%%%%%%
%INPUTS%
%%%%%%%%
%Pictures are in same folder as script, named as 'size'um'date'.tif
%example - 10um20140312.tif

%Centers are in same folder as script, named as centers'size'um'date.mat
%example - centers10um20140423.mat

%Dates of files being considered of the form YYMMDD. 20 is appended later
dates = [140317,140319,140402,140411,140423];
%dates = 140416;%[140317,140319,140326,140402,140411,140423];
%size of bucket being analyzed, in um.
bucketsize = 2;

%Range of ROIs to analyze. Radius of ROI, such that diameter = Radius*2 + 1
%example, ROI 10 has radius 10 and diameter 21. start and end range inclusive.
ROIRange = 5;%1:1:70;

%Threshold where white noise dominated traces are thrown away. Higher value
%of threshold corresponds with less traces thrown away.
%example: 0.5 means a trace that drops below 0.5 of 0 lag value is thrown
%away.
WhiteNoiseThreshold = -inf;

%%%%%%%%%%%
%MAIN LOOP%
%%%%%%%%%%%
%For each day being analyzed
for ii = 1:length(dates)
    
    %Grab the date being used
    date = sprintf('20%g',dates(ii)); %04022014, 04162014, 04232014
    disp(date)
    
    %Certain bucket sizes of a certain day only have 22 time points
    if strcmp(date,'20140312') == 1 && bucketsize ~= 5
        tMax = 22;
    else
        tMax = 23;
    end

    %For each ROI value in ROI range
    for n = 1:length(ROIRange)
        disp('ROI')
        disp(n);
        %Load storage file
        ROI = ROIRange(n);
        ROIArea(n) = (ROI+.5)^2*pi();
        DynamicFileName = sprintf('%gum%sROI%gTraces.mat',bucketsize,date,ROI);
        load(DynamicFileName);
        
        %Preallocate arrays
        buckets = length(BucketIdx);
        Avgcv2Tot = zeros(buckets,1);
        t50ArrayTot = Avgcv2Tot;
        GillespieDatas = zeros(tMax,buckets);
        AutosAll = zeros(2*tMax-1,buckets);
        AllGillespieDataA = zeros(tMax,buckets);
        tspan = 0:3:67;

        %remove first point for each trace
        for j = 1:buckets
            GillespieData(:,j) = GillespieData(:,j) - GillespieData(1,j);
        end


        len = length(GillespieData(:,1));
        %calculate mean
        GillespieDataAvg = mean(GillespieData,2);

        %Preallocated arrays
        GillespieDataA = zeros(len,buckets);
        SSlevelsArray = zeros(buckets,1);
        SSnormArray = zeros(buckets,1);
        Avgcv2 = zeros(buckets,1);
        Avgcv2Norm = zeros(buckets,1);
        AutoArray = zeros(len*2-1,buckets);
        NoiseStr = Avgcv2;
        t50Array = Avgcv2;
        GainStore = Avgcv2;
        VarianceArray = Avgcv2;

        %%%%%%%%%%%%%%%
        %GAIN ANALYSIS%
        %%%%%%%%%%%%%%%
        %Removing the general trend with a gain from each bucket
        

        for k = 1:buckets
        %         disp(k)

            SSlevelsArray(k) = GillespieData(end,k);

            %find gain that minimizes trace
            %coarse grain
            GainRange = 0:.001:2;
            GainsTemp = zeros(1,length(GainRange));

            for h = 1:length(GainRange)
                DataTemp = GillespieData(:,k) - GainRange(h).*GillespieDataAvg;
                DataTemp2 = xcorr(DataTemp,GillespieDataAvg,'biased');

                GainsTemp(h) = abs(DataTemp2(len));
            end

            [~, GainMinIdx] = min(GainsTemp);

            %Finer grain 
            GainRange = GainRange(GainMinIdx) - .002:.00001: GainRange(GainMinIdx) + .002;
            GainsTemp = zeros(1,length(GainRange));

            for h = 1:length(GainRange)
                DataTemp = GillespieData(:,k) - GainRange(h).*GillespieDataAvg;
                DataTemp2 = xcorr(DataTemp,GillespieDataAvg,'biased');

                GainsTemp(h) = abs(DataTemp2(len));
            end

            [GainMin, GainMinIdx] = min(GainsTemp);

            %Store the gain values calculated
            GainStore(k) = GainRange(GainMinIdx);
            %Generate the noise traces
            GillespieDataA(:,k) = GillespieData(:,k) - GainRange(GainMinIdx).*GillespieDataAvg;
            %generate Autocorrelations
            AutoArrayTemp = xcorr(GillespieDataA(:,k),'biased');%,length(tspan2)-1);
            AutoArray(:,k) = AutoArrayTemp;
            VarianceArray(k) = AutoArrayTemp(len);
            Avgcv2(k) = AutoArrayTemp(len)/(SSlevelsArray(k))^2;

            %calculate t50
            halfmax = AutoArray(len,k)/2;
            t50Idx = len;
            for j = len:length(AutoArray(:,k))
                if AutoArray(j,k) < halfmax
                    t50Idx = j;
                    break
                end
            end

            t50Array(k) = t50Idx-len;

        end
        
        %Remove traces with large white noise component until all traces
        %have sufficient real noise
        WhiteNoiseRemoval = 1;
        BucketIdxStart = BucketIdx;
        BucketIdxTemp = BucketIdx;
        tempAutoArray = AutoArray;
        tempAvgcv2 = Avgcv2;
        tempSSArray = SSlevelsArray;
        tempVar = VarianceArray;
        tempGillespieDataA = GillespieDataA;
        tempGillespieData = GillespieData;
        counter = 0;
        while WhiteNoiseRemoval == 1
            counter = counter + 1;
            %Remove traces with large white noise component
            Newcv2 = 0;
            NewSS = 0;
            NewVar = 0;
            NewGillespieDataA = zeros(tMax,1);
            NewGillespieData = zeros(tMax,1);
            NewAutoArray = zeros(tMax*2-1,1);
            count = 1;

            BucketIdxFinal = 0;
            for k = 1:length(BucketIdxTemp)
                Threshold = tempAutoArray(len,k)*WhiteNoiseThreshold;
                if tempAutoArray(len+1,k) < Threshold
                else
                    NewGillespieData(:,count) = tempGillespieData(:,k);
                    NewGillespieDataA(:,count) = tempGillespieDataA(:,k);
                    NewAutoArray(:,count) = tempAutoArray(:,k);
                    Newcv2(count) = tempAvgcv2(k);
                    NewSS(count) = tempSSArray(k);
                    NewVar(count) = tempVar(k);
                    BucketIdxFinal(count) = BucketIdxTemp(k);
                    count = count + 1;
                end  
            end
            
            if length(BucketIdxFinal) == length(BucketIdxTemp)
                WhiteNoiseRemoval = 0;
                Trig = 1;
                break
            end
            
            if length(BucketIdxFinal) <= 1
                WhiteNoiseRemoval = 0;
                if counter == 1
                    
                end
                Trig = 2;
                break
            end
            
            %Redo gain calculation with new traces
            BucketIdxTemp = BucketIdxFinal;
            tempGillespieData = zeros(tMax,length(BucketIdxTemp));
            tempGillespieDataA = zeros(tMax,length(BucketIdxTemp));
             
            for k = 1:length(BucketIdxTemp)
                [temp, tempIdx] = find(BucketIdx == BucketIdxTemp(k));
                tempGillespieData(:,k) = GillespieData(:,tempIdx);
            end
            
            tempGillespieAvg = mean(tempGillespieData,2);
            tempAutoArray = zeros(tMax*2-1,length(BucketIdxTemp));
            tempAvgcv2 = zeros(length(BucketIdxTemp),1);
            tempSSArray = zeros(length(BucketIdxTemp),1);
            tempVar = zeros(length(BucketIdxTemp),1);
            
            for k = 1:length(BucketIdxTemp)
            %         disp(k)

                tempSSArray(k) = tempGillespieData(end,k);

                %find gain that minimizes trace
                %coarse grain
                GainRange = 0:.001:2;
                GainsTemp = zeros(1,length(GainRange));

                for h = 1:length(GainRange)
                    DataTemp = tempGillespieData(:,k) - GainRange(h).*tempGillespieAvg;
                    DataTemp2 = xcorr(DataTemp,tempGillespieAvg,'biased');

                    GainsTemp(h) = abs(DataTemp2(len));
                end

                [~, GainMinIdx] = min(GainsTemp);

                %Finer grain 
                GainRange = GainRange(GainMinIdx) - .002:.00001: GainRange(GainMinIdx) + .002;
                GainsTemp = zeros(1,length(GainRange));

                for h = 1:length(GainRange)
                    DataTemp = tempGillespieData(:,k) - GainRange(h).*tempGillespieAvg;
                    DataTemp2 = xcorr(DataTemp,tempGillespieAvg,'biased');

                    GainsTemp(h) = abs(DataTemp2(len));
                end

                [GainMin, GainMinIdx] = min(GainsTemp);

                %Generate the noise traces
                tempGillespieDataA(:,k) = tempGillespieData(:,k) - GainRange(GainMinIdx).*tempGillespieAvg;
                %generate Autocorrelations
                AutoArrayTemp = xcorr(tempGillespieDataA(:,k),'biased');%,length(tspan2)-1);
                tempAutoArray(:,k) = AutoArrayTemp;
                tempVar(k) = AutoArrayTemp(len);
                tempAvgcv2(k) = AutoArrayTemp(len)/(tempSSArray(k))^2;

            end

            
        end
                

        %Restore new data
        GoodVarArray = NewVar;
        Goodcv2 = Newcv2;
        GoodSSArray = NewSS;
        GoodAutoArray = NewAutoArray;
        GoodGillespieDataA = NewGillespieDataA;
        GoodGillespieData = NewGillespieData;
        

        %Store Data

        CentroidVar = mean(VarianceArray);
        CentroidVarStd = std(GoodVarArray);
        CentroidSS = mean(SSlevelsArray);
        CentroidSSArea = mean(GoodSSArray)/ROIArea(n);

        CentroidSSStd = std(GoodSSArray);
        Centroidcv2 = CentroidVar/CentroidSS^2;
        Centroidcv2Std = std(GoodVarArray./GoodSSArray.^2);
        %pause
        DynamicFileName = sprintf('%gum%sROI%gResults.mat',bucketsize,date,ROI);
        save(DynamicFileName,'CentroidSS','Centroidcv2','Avgcv2','SSlevelsArray','GillespieData'...
            ,'GillespieDataA','AutoArray','BucketIdxStart','BucketGroup','WhiteNoiseThreshold',...
            'BucketIdxFinal','GoodVarArray','Goodcv2','GoodSSArray','GoodAutoArray',...
            'GoodGillespieDataA','GoodGillespieData','CentroidVar','VarianceArray');
    end
    
end

