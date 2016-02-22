clear all
close all
rand('state',sum(100*clock)); %#ok<RAND>

FullTrials = 10;

GeneArray = [5,10,15,20,25,30,35,40,45,50];
RunsArray = ones(50,1);%[20,20,20,20,20,20,20,20,20,20];%[120,60,40,30,24,20,17,14,13,12];
RunsArray = RunsArray .*50;

AllRibosomesPerGene = zeros(FullTrials,max(RunsArray),max(GeneArray));
AllmRNAperGene = zeros(FullTrials,max(RunsArray),max(GeneArray));
FirstBurstSizeTrack = zeros(FullTrials,max(RunsArray));
FirstBurstTimeTrack = zeros(FullTrials,max(RunsArray));
AllBurstTrack = zeros(FullTrials,max(RunsArray),max(GeneArray));
for k = 1:FullTrials
    Case = k;
    disp(Case)

    

    NumGenes = GeneArray(Case);
    RibosomesInitial = 300*NumGenes;
    Runs = RunsArray(k);
    
    kB = .001;
    kUb = .0001;
    kON = .001;
    kOFF = 2;
    alpha = 5; %20 seconds a transcript, 3 per min
    
    tMax = 100;
    dt = 1;
    tspan = 0:dt:tMax;
    tspan2 = length(tspan);



    a = zeros(NumGenes*2,1);
    TraceVar = zeros(Runs,NumGenes);
    TraceAuto = zeros(tspan2,Runs*NumGenes);
    TraceAutoNorm = TraceAuto;

    Genes = zeros(1,NumGenes*2);
    Genes(NumGenes+1:end) = 1;
    BoundRibosomeArray = zeros(1000,1);
    mRNATotals = zeros(1,NumGenes);
    aGenes = zeros(NumGenes*2,1);
    
    GeneRuns = zeros(Runs,tspan2,NumGenes);
    BurstTrackRuns = zeros(Runs,NumGenes);
    BurstSizeRuns = zeros(Runs,1000,NumGenes);
    BurstTimesRuns = zeros(Runs,1000,NumGenes);
    mRNATotalsRuns = zeros(Runs,NumGenes);
    RibosomesPerGene = zeros(Runs,NumGenes);
    
    BurstDurationTrack = zeros(100000,Runs);
    
    %Build Rxn Matrix and Propensity for genes
    RxnMatrixGene = zeros(NumGenes*2,NumGenes*2);
    for i = 1:NumGenes
        RxnMatrixGene(i,i) = 1;
        RxnMatrixGene(i,NumGenes+i) = -1;
        aGenes(i) = kON;
    end
    for i = 1:NumGenes
        RxnMatrixGene(i+NumGenes,i) = -1;
        RxnMatrixGene(i+NumGenes,NumGenes+i) = 1;
        aGenes(i+NumGenes) = kOFF;
    end
    
    %
    
    for i = 1:Runs

        MaxOutput = 100000; %Maximum size expected of the output file
        %T = zeros(MaxOutput, 1); %Time tracking
        %X = zeros(1, NumSpecies); %Species number tracking
        %T(1)     = 0;
        %X(1,:)   = x0;
        RxnCount = 1;
        T = 0; %Time
        mRNAArray = zeros(1000,3);
        Ttrack = zeros(length(tspan), 1); %Time tracking
        XGenes = zeros(length(tspan), NumGenes*2); %Species number tracking
        XmRNA = zeros(length(tspan),length(mRNAArray(:,1)));
        XRibosomes = XmRNA;
        Ribosomes = RibosomesInitial;
        XGenes(1,:)   = Genes;
        xCurrentGenes = Genes;
        xCurrentmRNA = zeros(1000,1);
        xCurrentRibosomes = zeros(1000,1);
        mRNATotals = zeros(1,NumGenes);
        BurstTrack = zeros(1,NumGenes); %track burst index
        BurstSizeTrack = zeros(1000,NumGenes);%brst size track
        BurstTimes = zeros(1000,NumGenes);%burst time track
        OnTrack = zeros(1,NumGenes); %Burst state track
    
        RecordTime = dt; %Recording time
        RecordCount = 2;
        mRNAnum = 0;
        count = 1;
        OnDuration = zeros(MaxOutput,1);
        OnTimes = zeros(MaxOutput,1);
        %%%%%%%%%%%%%%%%%%%%%%
        %Gillespie Simulation%
        %%%%%%%%%%%%%%%%%%%%%%

        while T <= tMax

            % Calculate reaction propensities
            x = xCurrentGenes;

            for j = 1:NumGenes
                aGenes(j) = kON*xCurrentGenes(NumGenes+j);
                aGenes(j+NumGenes) = kOFF*xCurrentGenes(j);
                aAlpha(j) = xCurrentGenes(j)*alpha;
            end
            
            %binding propensities
            aRibosomes = xCurrentmRNA(:,1) .* Ribosomes .* kB;
            %unbinding propensities
            aRibosomesU = xCurrentRibosomes .* kUb;
            
            
            % Compute tau and mu using random variables
            a0 = sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + sum(aRibosomesU);
            r = rand(1,2);
            %tau = -log(r(1))/a0;
            tau = (1/a0)*log(1/r(1));

            %Store information if at time
            tau0 = tau;
            if tau0 > dt
                while tau0 > dt
                    XGenes(RecordCount,:) = xCurrentGenes;
                    XmRNA(RecordCount,:) = xCurrentmRNA(:,1);
                    XRibosomes(RecordCount,:) = xCurrentRibosomes;
                    RecordCount = RecordCount + 1;
                    RecordTime = RecordTime + dt;
                    tau0 = tau0 - dt;
                end
            end

            if T + tau > RecordTime
                XGenes(RecordCount,:) = xCurrentGenes;
                XmRNA(RecordCount,:) = xCurrentmRNA(:,1);
                XRibosomes(RecordCount,:) = xCurrentRibosomes;
                RecordCount = RecordCount + 1;
                RecordTime = RecordTime + dt;

            end
            %[~, mu] = histc(r(2)*a0, [0;cumsum(a(:))]);
            if r(2)*a0 < sum(aGenes)
                
                mu  = find((cumsum(aGenes) >= r(2)*a0),1,'first');
                T   = T  + tau;
                xCurrentGenes = xCurrentGenes + RxnMatrixGene(mu,:);
                if mu <= NumGenes
                    OnTrack(mu) = 1;
                    BurstTrack(mu) = BurstTrack(mu) + 1;
                    
                    
                else
                    OnTrack(mu-NumGenes) = 0;
                end
                
            elseif r(2)*a0 < sum(aGenes) + sum(aAlpha)
                mu  = find((sum(aGenes) + cumsum(aAlpha) >= r(2)*a0),1,'first');
                T   = T  + tau;
                mRNATotals(mu) = mRNATotals(mu) + 1;
                mRNAnum = mRNAnum + 1;
                xCurrentmRNA(mRNAnum,1) = 1;
                mRNAArray(mRNAnum,1) = 1;
                mRNAArray(mRNAnum,2) = mu;
                mRNAArray(mRNAnum,3) = T;
                if BurstSizeTrack(BurstTrack(mu),mu) == 0
                    BurstTimes(BurstTrack(mu),mu) = T;
                end
                BurstSizeTrack(BurstTrack(mu),mu) = BurstSizeTrack(BurstTrack(mu),mu) + 1;
                
            elseif r(2)*a0 < sum(aGenes) + sum(aAlpha) + sum(aRibosomes)
                mu  = find((sum(aGenes) + sum(aAlpha) + cumsum(aRibosomes) >= r(2)*a0),1,'first');
                %find the next change in state before tau
                T   = T  + tau;
                %xCurrentGenes = xCurrentGenes + RxnMatrixGene(mu,:);
                xCurrentRibosomes(mu) = xCurrentRibosomes(mu) + 1;
                Ribosomes = Ribosomes - 1;
            else
                mu  = find((sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + cumsum(aRibosomesU) >= r(2)*a0),1,'first');
                %find the next change in state before tau
                T   = T  + tau;
                %xCurrentGenes = xCurrentGenes + RxnMatrixGene(mu,:);
                xCurrentRibosomes(mu) = xCurrentRibosomes(mu) - 1;
                Ribosomes = Ribosomes + 1;
            end
            
            

        end


        % Record output
         XGenes(RecordCount,:) = xCurrentGenes;
         XmRNA(RecordCount,:) = xCurrentmRNA(:,1);
         XRibosomes(RecordCount,:) = xCurrentRibosomes;
         
         mRNATotnum(i) = sum(mRNATotals);
         for j = 1:mRNATotnum(i)
            RibosomesPerGene(i,mRNAArray(j,2)) = RibosomesPerGene(i,mRNAArray(j,2)) + xCurrentRibosomes(j);
         end
         
        % Record output
        count = count - 1;
        mRNATotalsRuns(i,:) = mRNATotals;
        GeneRuns(i,:,:) = XGenes(1:tspan2,1:NumGenes);
        mRNARuns(i,:,:) = XmRNA(1:tspan2,:);
        RibosomeRuns(i,:,:) = XRibosomes(1:tspan2,:);
        mRNAIdx(i,:) = mRNAArray(:,2); 
        mRNAArrayRuns(i,:,:) = mRNAArray;
        BurstTrackRuns(i,:) = BurstTrack;
        BurstSizeRuns(i,:,:) = BurstSizeTrack;
        BurstTimesRuns(i,:,:) = BurstTimes;
        FreeRibosomes(i) = Ribosomes;
        
%         [temp1,Idx] = find(BurstSizeTrack);
%         
%         if isempty(Idx)
%         else
%         FirstBurstSizeTrack(k,i) = BurstSizeTrack(1,Idx);
%         FirstBurstTimeTrack(k,i) = BurstTimes(1,Idx);
%         end


    end
   
   AllRibosomeRuns(k,:,:,:) = RibosomeRuns;
   AllmRNARuns(k,:,:,:) = mRNARuns;
   AllmRNAArray(k,:,:,:) = mRNAArrayRuns;
   AllRibosomesPerGene(k,1:Runs,1:NumGenes) = RibosomesPerGene;
   AllmRNATot(k,:) = mRNATotnum;
   AllmRNAperGene(k,1:Runs,1:NumGenes) = mRNATotalsRuns;
   AllmRNAIdx(k,:,:) = mRNAIdx;
   AllBurstTrack(k,1:Runs,1:NumGenes) = BurstTrackRuns;
   AllBurstSize(k,1:Runs,:,1:NumGenes) = BurstSizeRuns;
   AllBurstTimes(k,1:Runs,:,1:NumGenes) = BurstTimesRuns;
   AllFreeRibosomes(k,:) = FreeRibosomes;
   
end
save AllData
%%
%Calculations
for k = 1:FullTrials
    for j = 1:Runs
        temp = AllRibosomeRuns(k,j,:,:);
        temp = reshape(temp,[tspan2,1000]);
        tempsum = sum(temp,1);
        [RibosomeIdx, misc] = find(tempsum);
        for i = 1:length(RibosomeIdx)
            EndRibosomeCount(k,j,i) = temp(end,i);
            
        end
        tempsum = sum(temp,2);
        AllTotalRibosomes(k,j,:) = tempsum;
        EndmRNAwRibosomes(k,j) = length(RibosomeIdx); 
    end
    
    for j = 1:Runs
        temp = AllmRNARuns(k,j,:,:);
        temp = reshape(temp,[tspan2,1000]);
        tempsum = sum(temp,1);
    end
    
    temp = mean(AllTotalRibosomes(k,:,:),2);
    TotalRibosomesGenTrend(k,:) = temp(:);
    
   
    for j = 1:Runs
        TotalRibosomeSS(k,j) = AllTotalRibosomes(k,j,end);
        temp = AllTotalRibosomes(k,j,:);
        temp2 = TotalRibosomesGenTrend(k,:);
        TotalRibosomeA(k,j,:) = temp(:) - temp2(:);
        TotalRibosomeVar(k,j) = var(TotalRibosomeA(k,j,:));
        TotalRibosomecv2(k,j) =  TotalRibosomeVar(k,j)/TotalRibosomeSS(k,j)^2;
        
    end
    

end

%%
%find all bursts with 0s included
hold on
c = colormap(jet(FullTrials));
edges = 1:1:20;
AllBursts = 0;
for k = 1:FullTrials
    count = 1;
    
    for j = 1:Runs
        %find the indecies and values of all bursts
         temp = AllBurstTrack(k,j,1:GeneArray(k));
         temp = temp(:);
         [Idx,~,Bursts] = find(temp);
         %find the number of mRNA made in a burst
         temp2 = AllBurstSize(k,j,:,1:GeneArray(k));
         temp2 = reshape(temp2,[1000,GeneArray(k)]);
         %find total mRNA made per Gene
         temp1 = AllmRNAperGene(k,j,:);
         temp1 = temp1(:);
         %find end ribosome values per gene
         temp3 = AllRibosomesPerGene(k,j,:);
         temp3 = temp3(:);
         for i = 1:GeneArray(k)
             if temp1(i) == 0;
             else 
             AllBursts(i,j,k) = temp1(i)/temp(i);
             AllRibosomespermRNA(i,j,k) = temp3(i)/temp1(i);
             end
         end
    end
%     h = histogram(AllBursts(:,:,k),'Normalization','probability');
%     h.BinWidth = 1;
%     h.FaceColor = c(k,:);

end

% xlabel('Burst Size (mRNA per burst)','FontSize',15)
% ylabel('Probability','FontSize',15)
% title('Transcriptional BS with 0s')
% saveas(gcf,'mRNABSwith0sHistogram.jpg')

%%
figure
hold on
AllBurstsNo0 = AllBursts;
AllBurstsNo0(AllBurstsNo0 == 0) = NaN;

for k = 1:FullTrials
     h = histogram(AllBurstsNo0(:,:,k),'Normalization','probability');
    h.BinWidth = 1;
    h.FaceColor = c(k,:);
end

xlabel('Burst Size (mRNA per burst)','FontSize',15)
ylabel('Probability','FontSize',15)
title('Transcriptional BS without 0s')
saveas(gcf,'mRNABSno0sHistogram.jpg')

%%
figure
hold on
for k = 1:FullTrials
    tempx = AllBursts(:,:,k);
    tempy = AllRibosomespermRNA(:,:,k);
    tempx(tempx==0) = NaN;
    tempy(tempy==0) = NaN;
    tempx = tempx(:);
    tempy = tempy(:);
    meanX(k) = nanmean(tempx);
    meanY(k) = nanmean(tempy);
     plot(tempx,tempy,'linestyle','none','marker','.','markersize',10,...
         'color',c(k,:));
end

for k = 1:FullTrials
     plot(meanX(k),meanY(k),'linestyle','none','marker','o','markersize',10,...
         'markerfacecolor',c(k,:),'markeredgecolor','k');
end
set(gca,'YScale','log');
xlabel('Burst Size (mRNA per burst)','FontSize',15)
ylabel('Ribosomes per mRNA','FontSize',15)
title('TBS without 0s vs Ribosomes per mRNA')
saveas(gcf,'mRNABSvsRibosomespermRNA.jpg')
%%
%plot ribosomes per gene vs time of first mRNA
hold on
timeON = zeros(FullTrials, Runs, max(GeneArray));
RibosomeEnd = zeros(FullTrials, Runs, max(GeneArray));
c = colormap(jet(10));
for k = FullTrials:-1:1
    for j = 1:Runs
        temp = AllRibosomeRuns(k,j,:,1:GeneArray(k));
        tempplot = reshape(temp,[tspan2,GeneArray(k)]);
        for i = 1:GeneArray(k)
            if temp(1,1,end,i) > 0
                timeON(k,j,i) = find(temp(1,1,:,i) > 0,1,'first');
                RibosomeEnd(k,j,i) = temp(1,1,end,i);
            else
            end
        end
        %plot(tempplot,'color',c(k,:))
    end
    %timeON(k,:,GeneArray(k)+1:end)
end

for k = 1:FullTrials
    for j = 1:Runs
        RibosomeEnd(k,j,GeneArray(k)+1:end) = NaN;
    end
end

figure
hold on
timeOn(timeON == 0) = NaN;
for k = 1:FullTrials
    name = sprintf('%g Genes',GeneArray(k));
    temp1 = timeON(k,:,:);
    temp2 = RibosomeEnd(k,:,:);
    tempm1 = temp1;
    tempm1(tempm1 == 0) = NaN;
    tempm2 = temp2;
    tempm2(tempm2 == 0) = [];
    
    meanTime(k) = nanmean(tempm1(:));
    meanRib(k) = nanmean(tempm2(:));
    medianRib(k) = nanmedian(tempm2(:));
    modeRib(k) = mode(tempm2(:));
    linestore2(k) = plot(temp1(:),temp2(:),'linestyle','none','marker','.','markersize',12,...
        'color',c(k,:),'displayname',name);
end

for k = 1:FullTrials
    plot(meanTime(k),meanRib(k),'linestyle','none','marker','o','markersize',10,...
        'markerfacecolor',c(k,:),'markeredgecolor','k');
end
legend(linestore2,'location','northeast')
xlabel('Time of First Burst','FontSize',15)
ylabel('Total Ribosomes Per Gene','FontSize',15)
title('Ribosomes per Gene vs First Burst')
saveas(gcf,'RibsomesTotalvsGeneTime.jpg')
%% plot ribosome per gene vs mRNA per gene
figure
hold on
for k = 1:FullTrials
    name = sprintf('%g Genes',GeneArray(k));
    temp1 = AllmRNAperGene(k,:,:);
    temp2 = RibosomeEnd(k,:,:);
    tempm1 = temp1;
    tempm1(tempm1 == 0) = NaN;
    tempm2 = temp2;
    tempm2(tempm2 == 0) = NaN;
    
    meanmRNA(k) = nanmean(tempm1(:));
    meanRib(k) = nanmean(tempm2(:));
    linestore2(k) = plot(tempm1(:),tempm2(:),'linestyle','none','marker','.','markersize',12,...
        'color',c(k,:),'displayname',name);
end

for k = 1:FullTrials
    plot(meanmRNA(k),meanRib(k),'linestyle','none','marker','o','markersize',10,...
        'markerfacecolor',c(k,:),'markeredgecolor','k');
end
legend(linestore2,'location','northeast')
xlabel('Total mRNA per Gene','FontSize',15)
ylabel('Total Ribosomes Per Gene','FontSize',15)
title('Ribosomes per Gene vs mRNA per Gene')
saveas(gcf,'RibsomesperGenevsmRNAperGene.jpg')
%%
% %plot Total Ribosomes over time
% figure
% c = colormap(jet(FullTrials));
% hold on
% for k = 1:FullTrials
%     for j = 1:Runs
%         temp = AllTotalRibosomes(k,j,:);
%         plot(temp(:),'color',c(k,:))
%     end
% end
% xlabel('Time','FontSize',15)
% ylabel('Total Ribosomes','FontSize',15)
% title('Total Ribosome over time')
% saveas(gcf,'RibsomesTotalOverTime.jpg')
%%
%plot cv2 vs abundance
figure
hold on
c = colormap(jet(FullTrials));
for k = 1:FullTrials
    
    plot(TotalRibosomeSS(k,:),TotalRibosomecv2(k,:),'color',c(k,:),...
        'linestyle','none','marker','.','markersize',10)

end

for k = 1:FullTrials
    tempSS = TotalRibosomeSS(k,:);
    tempSS(tempSS == 0) = NaN;
    s = isnan(tempSS);
    SSmean(k) = geomean(tempSS(~s));
    tempcv2 = TotalRibosomecv2(k,:);
    tempcv2(tempcv2 == Inf) = NaN;
    
    s = isnan(tempcv2);
    cv2mean(k) = geomean(tempcv2(~s));
    plot(SSmean(k),cv2mean(k),'marker','o','markerfacecolor',c(k,:),...
        'linestyle','none','markersize',8,'markeredgecolor','k')

end

offset = [10^1,10^5,10^6,10^7,10^8];
for k = 1
    xrange = 100:10000;
    yrange = offset(k)./(xrange)  ;
    plot(xrange,yrange,'color','k')
end
set(gca,'XScale','log');
set(gca,'YScale','log');
axis([10 100000 10^-4 100]);
xlabel('Ribosome Abundance','FontSize',15)
ylabel('Ribosome cv^2','FontSize',15)
title('Total Ribosome cv^2 vs Abundance')
saveas(gcf,'Ribsomescv2abundance.jpg')
%%
figure
hold on
for k = 1:FullTrials
    AllPercentFreeRibosomes(k,:) = AllFreeRibosomes(k,:)./(GeneArray(k)*300);
    plot(AllmRNATot(k,:),AllFreeRibosomes(k,:)./(GeneArray(k)*300),'linestyle','none','marker',...
        '.','color',c(k,:),'markersize',10)
end

%set(gca,'YScale','log');
xlabel('Total mRNA','FontSize',15)
ylabel('Total Percentage Free Ribosomes','FontSize',15)
title('Free Ribosomes vs mRNA Abundance')
saveas(gcf,'FreeRibsomesvsmRNANorm.jpg')
%%
figure
hold on
for k = 1:FullTrials
    plot(AllmRNATot(k,:),AllFreeRibosomes(k,:),'linestyle','none','marker',...
        '.','color',c(k,:),'markersize',10)
end
%set(gca,'YScale','log');
xlabel('Total mRNA','FontSize',15)
ylabel('Total Free Ribosomes','FontSize',15)
title('Free Ribosomes vs mRNA Abundance')
saveas(gcf,'FreeRibsomesvsmRNA.jpg')
% %%
% %%
% figure
% hold on
% c = colormap(jet(FullTrials));
% for k = 1:FullTrials
%     for j = 1:Runs
%         temp = AllRibosomeRuns(k,j,:,:);
%         temp = reshape(temp,[tspan2,1000]);
%         tempsum = sum(temp,1);
%         [RibosomeIdx, misc] = find(tempsum);
%         for i = 1:length(RibosomeIdx)
%             EndRibosomeCount(k,j,i) = temp(end,i);
%             
%         end
%         tempsum = sum(temp,2);
%         AllTotalRibosomes(k,j,:) = tempsum;
%         %[f,x] = ecdf(temp(end,1:length(RibosomeIdx)));
%         %plot(x,f,'color',c(k,:))
%         if isempty(RibosomeIdx)
%         else
%         plot(GeneArray(k),temp(end,1:length(RibosomeIdx)),'color',c(k,:),...
%             'linestyle','none','marker','.','markersize',10)
%         EndmRNAwRibosomes(k,j) = length(RibosomeIdx); 
%         end
%     end
% end
% xlabel('Concurrent Genes','FontSize',15)
% ylabel('Ribosomes per Transcript','FontSize',15)
% title('Number of Ribsomes per transcript')
% saveas(gcf,'RibsomespermRNA.jpg')
%%
% figure
% hold on
% for k = 1:FullTrials
% 
%     for j = 1:Runs
%         %AllRibosomesPerGene(k,j,GeneArray(k)+1:end) = NaN;
%         RibosomeEnd(k,j,GeneArray(k)+1:end) = NaN;
%         timeON(k,j,GeneArray(k)+1:end) = NaN;
%         AllmRNAperGene(k,j,GeneArray(k)+1:end) = NaN;
% %         temp = AllRibosomesPerGene(k,j,1:GeneArray(k));
% %         temp = temp(:);
% %         plot(GeneArray(k),temp,'linestyle','none',...
% %             'marker','.','color',c(k,:),'markersize',10)
%     end
% end
% %AllRibosomesPerGeneBox = reshape(AllRibosomesPerGene,[10,Runs*50]);
% AllRibosomesPerGeneBox = reshape(RibosomeEnd,[10,Runs*50]);
% TimeONPerGeneBox = reshape(timeON,[10,Runs*50]);
% AllRibosomesPerGeneBox = AllRibosomesPerGeneBox';
% TimeONPerGeneBox = TimeONPerGeneBox';
% 
% for k = 1:FullTrials
%     AllRibosomesPerGeneNorm(:,k) = AllRibosomesPerGeneBox(:,k) ./ max(AllRibosomesPerGeneBox(:,k));
% end
% 
% AllRibosomesPerGeneNorm(AllRibosomesPerGeneNorm == 0) = NaN;
% AllRibosomesPerGeneBox(AllRibosomesPerGeneBox == 0) = NaN;
% 
% for k = 1:FullTrials
%     [num, idx] = find(AllRibosomesPerGeneNorm(:,k)> 0);
%     OnRibosomes(k) = length(idx);
% end
% BoundFraction = OnRibosomes(k) ./ (GeneArray .* Runs);
% 
% boxplot(AllRibosomesPerGeneBox,GeneArray)
% xlabel('Concurrent Genes','FontSize',15)
% ylabel('Ribosomes per Gene','FontSize',15)
% title('Ribosomes per Gene Box plot')
% saveas(gcf,'RibsomesperGeneBox.jpg')
% figure
% 
% boxplot(AllRibosomesPerGeneNorm,GeneArray)
% xlabel('Concurrent Genes','FontSize',15)
% ylabel('Normalized Ribosomes per Gene','FontSize',15)
% title('Ribosomes per gene Box plot')
% saveas(gcf,'RibsomesperGeneBoxNorm.jpg')
% 
% %%
% figure
% hold on
% edges = 0:30:1500;
% for k = 1:FullTrials
%     h = histogram(AllRibosomesPerGeneBox(1:GeneArray(k)*Runs,k),edges,'Normalization','probability');
%     h.FaceColor = c(k,:);
% 
% end
% xlabel('Bound Ribosomes per Gene','FontSize',15)
% ylabel('Probability','FontSize',15)
% title('Ribosomes per gene Histogram')
% saveas(gcf,'RibsomesperGeneHistR.jpg')
% %%
% figure
% hold on
% c = colormap(jet(FullTrials));
% for k = 1:FullTrials
%     for j = 1:Runs
%         temp = AllRibosomesPerGene(k,j,1:GeneArray(k));
%         temp = temp(:);
%         temp(temp == 0) = [];
%         meanBS(k,j) = mean(temp);
%         meanBF(k,j) = length(temp);
%         plot(meanBF(k,:), meanBS(k,:),'linestyle','none',...
%             'marker','.','color',c(k,:),'markersize',12)
%     end
% end
% for k = 1:FullTrials
%     name = sprintf('%g Genes',GeneArray(k));
%     temp = meanBS(k,:);
%     temp(temp == 0) = NaN;
%     tempBS = nanmean(temp);
%     temp = meanBF(k,:);
%     temp(temp == 0) = NaN;
%     tempBF = nanmean(temp);
%     linestore(k) = plot(tempBF, tempBS,'linestyle','none',...
%         'marker','o','markerfacecolor',c(k,:),'markersize',8,...
%         'markeredgecolor','k','displayname',name);
% 
% end
% legend(linestore,'location','northeast')
% % set(gca,'XScale','log');
% % set(gca,'YScale','log');
% ylabel('Mean Intensity ON Gene(BS)','FontSize',15)
% xlabel('Mean Number ON (BF)','FontSize',15)
% title('Burst Frequency Total Ribosome Per Gene ')
% saveas(gcf,'BFvsBS.jpg')
%figure
%%
% figure
% hold on
% c = colormap(hsv(Runs));
% for k = 10
%     for j = 1:Runs
%         temp = AllRibosomesPerGene(k,j,1:GeneArray(k));
%         temp = temp(:);
%         temp(temp == 0) = [];
%         meanBS = mean(temp);
%         meanBF = length(temp);
%         plot(meanBF, meanBS,'linestyle','none',...
%             'marker','.','color',c(j,:),'markersize',10)
%     end
% end
% % set(gca,'XScale','log');
% % set(gca,'YScale','log');
% ylabel('Mean Intensity ON Gene(BS)','FontSize',15)
% xlabel('Mean Number ON (BF)','FontSize',15)
% title('Burst Frequency Total Ribosome Per Gene ')
% %for a particular number of genes and run
% c = colormap(hsv(Runs));
% figure
% hold on
% for k = 10
%     for j = 1:Runs
%         %grab the total ribosomes per gene
%         temp = AllRibosomesPerGene(k,j,1:GeneArray(k));
%         %reshape to make pretty
%         temp = temp(:);
%         %find the corresponding mRNAs to this gene index
%         temp2 = find(AllmRNAIdx(k,j,:)==14);
%         temp3 = zeros(101,1);
%         %for each mRNA
%         for i = 1:length(temp2)
%             %Find the time trace of the ribosomes to that mRNA
%             temp4 = AllRibosomeRuns(k,j,:,temp2(i));
%             %reshape
%             temp4 = temp4(:);
%             %sum them together to get a ribosome per gene trace
%             temp3 = temp3 + temp4;
%         end
%         plot(temp3,'color',c(j,:))
%     end
% end


% hold on
% for i = 1:Runs
%     tempplot = AllRibosomeRuns(10,1,:,i);
%     plot(tempplot(:))
% end
% edges = 0:20:400;
% for k = 1:FullTrials
%     figure
%     temp = AllRibosomesPerGene(k,:,1:GeneArray(k));
%     temp = temp(:);
%     histogram(temp(:),edges,'normalization','probability')
%     xlabel('Ribosomes Per Gene','FontSize',15)
%     ylabel('Count','FontSize',15)
%     name = sprintf('Total Ribosomes for %g Genes',GeneArray(k));
%     title(name)
%     name = sprintf('HistogramTotRibosome%gGenes.jpg',GeneArray(k));
%     saveas(gcf,name)
% end
