%Processing Noise Data to plots
%Outputs - Plots and visualizations by days, buckets, and ROIs.

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

%dates = 140317;%[140402,140411,140423];
%size of bucket being analyzed, in um.
bucketsize = 2;
maxbuckets = 9;

%Range of ROIs to analyze. Radius of ROI, such that diameter = Radius*2 + 1
%example, ROI 10 has radius 10 and diameter 21. start and end range inclusive.
if bucketsize == 10
    ROIRange = 23;%1:1:70;
elseif bucketsize == 5
    ROIRange = 10;
else
    ROIRange = 5;
end
  
PlotIndividual = 0;


count = 0;
count2 = 0;
%%%%%%%%%%%%%%%%%
%Data Agregation%
%%%%%%%%%%%%%%%%%
%For each day being analyzed
TotalBuckets = zeros(length(ROIRange),length(dates));
TotalGoodBuckets = TotalBuckets;
for ii = 1:length(dates)
    %Grab the date being used
    date = sprintf('20%g',dates(ii)); %04022014, 04162014, 04232014
    disp(date)
    
    %Certain bucket sizes of a certain day only have 22 time points
    if strcmp(date,'20140312') == 1 && bucketsize ~= 5
        tMax = 22;
        tspan = 0:3:64;
    else
        tMax = 23;
        tspan = 0:3:67;
    end

    %For each ROI value in ROI range
    for n = 1:length(ROIRange)
        
        
        ROI = ROIRange(n);
        DynamicFileName = sprintf('%gum%sROI%gResults.mat',bucketsize,date,ROI);
        load(DynamicFileName)
        
        FinalBucketNumber(n,ii) = length(BucketIdxFinal)/length(BucketIdxStart);
        Allcv2(1:length(Avgcv2),n,ii) = Avgcv2;
        for jj = 1:length(SSlevelsArray)
            count = count + 1;
            AllSS(count) = SSlevelsArray(jj);
        end
        AllVar(1:length(VarianceArray),n,ii) = VarianceArray;
        AllCentroidcv2(n,ii) = Centroidcv2;
        AllCentroidSS(n,ii) = CentroidSS;
        AllCentroidvar(n,ii) = CentroidVar;

        AllGillespieData(1:length(GillespieData(:,1)),1:length(GillespieData(1,:)),ii,n) = GillespieData;
        for jj = 1:length(GillespieDataA(1,:))
            count2 = count2 + 1;
            AllNoiseTraces(1:length(GillespieDataA(:,1)),count2) = GillespieDataA(:,jj);
        end
        AllAutocorrs(1:length(AutoArray(:,1)),1:length(AutoArray(1,:)),ii,n) = AutoArray;
        AllBucketIdx(1:length(BucketIdxStart),n,ii) = BucketIdxStart;
        AllBucketNum(n,ii) = length(BucketIdxStart);
        TotalBuckets(n,ii) = TotalBuckets(n,ii) + length(BucketIdxStart);
        TotalGoodBuckets(n,ii) = TotalGoodBuckets(n,ii) + length(BucketIdxFinal);
        
        %Plot idividual figures
        if PlotIndividual == 1
            %Plot Abundance Curves
            hold off
            LineStore = zeros(length(BucketIdxStart),1);
            c = colormap(hsv(maxbuckets));
            for j = 1:length(BucketIdxStart)
                DynamicName = sprintf('bucket%g',BucketIdxStart(j));
                LineStore(j) = plot(tspan,GillespieData(:,j),'DisplayName',DynamicName,'color',c(BucketIdxStart(j),:));
                hold on
            end
            legend(LineStore,'Location', 'northwest')

            %plot "good" curves thicker
            if BucketIdxFinal == 0
            else
                for j = 1:length(BucketIdxFinal)
                    [temp, tempIdx] = find(BucketIdxStart == BucketIdxFinal(j));
                    plot(tspan,GillespieData(:,tempIdx),'linewidth',3,'color',c(BucketIdxFinal(j),:));
                end
            end

            set(gca,'fontsize',10);
            xlabel('Time','FontSize',10)
            ylabel('Sum Flourescence','FontSize',10)
            %save image
            DynamicFileName = sprintf('Sum F Traces %g ROI %s',ROIRange(n),date);
            title(DynamicFileName)  
            DynamicFileName = sprintf('FTracesfor%gum%sROI%g.jpg',bucketsize,date,ROIRange(n));
            saveas(gcf,DynamicFileName)

            %Plot Noise Curves
            hold off
            LineStore = zeros(length(BucketIdxStart),1);
            for j = 1:length(BucketIdxStart)
                DynamicName = sprintf('bucket%g',BucketIdxStart(j));
                LineStore(j) = plot(tspan,GillespieDataA(:,j),'DisplayName',DynamicName,'color',c(BucketIdxStart(j),:));
                hold on
            end
            %legend(LineStore,'Location', 'northeast')
            %plot "good" curves thicker
            if BucketIdxFinal == 0
            else
            for j = 1:length(BucketIdxFinal)
                plot(tspan,GoodGillespieDataA(:,j),'linewidth',3,'color',c(BucketIdxFinal(j),:));
            end
            end
            set(gca,'fontsize',10);
            xlabel('Time','FontSize',10)
            ylabel('Noise','FontSize',10)
            %save image
            DynamicFileName = sprintf('Noise Traces %g ROI %s',ROIRange(n),date);
            title(DynamicFileName)  
            DynamicFileName = sprintf('NoiseTracesfor%gum%sROI%g.jpg',bucketsize,date,ROIRange(n));
            saveas(gcf,DynamicFileName)

            %Plot Autocorr Curves
            hold off
            LineStore = zeros(length(BucketIdxStart),1);
            lags = 0:tMax-1;
            for j = 1:length(BucketIdxStart)
                DynamicName = sprintf('bucket%g',BucketIdxStart(j));
                LineStore(j) = plot(lags,AutoArray(tMax:end,j),'DisplayName',DynamicName,'color',c(BucketIdxStart(j),:));
                hold on
            end
            legend(LineStore,'Location', 'northeast')
            %plot "good" curves thicker
            if BucketIdxFinal == 0
            else
            for j = 1:length(BucketIdxFinal)
                plot(lags,GoodAutoArray(tMax:end,j),'linewidth',3,'color',c(BucketIdxFinal(j),:));
            end
            end
            set(gca,'fontsize',10);
            xlabel('lags','FontSize',10)
            ylabel('Autocorr','FontSize',10)
            %save image
            DynamicFileName = sprintf('AutoCorr Traces %g ROI %s',ROIRange(n),date);
            title(DynamicFileName)  
            DynamicFileName = sprintf('AutoTracesfor%gum%sROI%g.jpg',bucketsize,date,ROIRange(n));
            saveas(gcf,DynamicFileName)
        end
        
        
        
    end
    
end
    DynamicFileName = sprintf('Compiled%gumData',bucketsize);
    save(DynamicFileName)
    DynamicFileName = sprintf('Compiled%gumSSNoiseData',bucketsize);
    save(DynamicFileName,'AllSS','AllNoiseTraces')
    
%%
%%%%%%%%%%%%%%%
%Data Analysis%
%%%%%%%%%%%%%%%
len = tMax;
temparray = AllAutocorrs(len:end,1,:,1);
i = 0;
Allcv2(Allcv2 == 0) = NaN;
for i = 1:length(dates)
    for j = 1:AllBucketNum(i,1)
        
        
        
        [mincv2(j,i), mincv2Idx(j,i)] = min(Allcv2(j,:,i));
        minROIbyBucket(j,i) = ROIRange(mincv2Idx(j,i))*2+1;
    end
    [MinByDay(i),minROIbyDayIdx(i)] = min(AllCentroidcv2(:,i));
    minROIbyDay(i) = ROIRange(minROIbyDayIdx(i))*2+1;
end

[~,minROITotal] = min(mean(AllCentroidcv2,2));
minROITotal = ROIRange(minROITotal)*2+1;

%Individual days cv2 vs ROI
% figure
c = colormap(hsv(maxbuckets));
for i = 1:length(dates)
    figure
    hold on

    LineStore = 0;
    for j = 1:length(BucketIdxStart)
        SStemp = zeros(length(ROIRange),1);
        cv2temp = SStemp;
        %SStemp(:) = AllSS(j,:,i);
        cv2temp(:) = Allcv2(j,:,i);
        DynamicName = sprintf('bucket%g',j);
        LineStore(j) = plot(ROIRange(:),cv2temp,'linestyle','none','marker','.',...
            'markersize',10,'color',c(BucketIdxStart(j),:),'DisplayName',DynamicName);
        hold on
    end
         plot(ROIRange(:),AllCentroidcv2(:,i),'linestyle','none',...
        'marker','.','markersize',12,'color','k')
    legend(LineStore,'Location','northeast');
    set(gca,'fontsize',10);
   % set(gca,'XScale','log');
    set(gca,'YScale','log');
    xlabel('ROI','FontSize',15)
    ax = gca;
    ax.XTick = 2:2:70;
    grid on
    ylabel('cv2','FontSize',15)
    tempdate = sprintf('20%g',dates(i));
    DynamicFileName = sprintf('Abundance vs ROI %s',tempdate);
    title(DynamicFileName)  
    DynamicFileName = sprintf('AbundanceROIfor%gum%s.jpg',bucketsize,tempdate);
    saveas(gcf,DynamicFileName)
end