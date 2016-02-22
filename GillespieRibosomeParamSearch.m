clear all
close all
rand('state',sum(100*clock)); %#ok<RAND>

%testParam = 'RibosomeNum';
testParam = 'alpha';

%ParamArray = [200,250,300,350,400]; %Number of Ribosomes per Gene
%ParamArray = [.0005,.001,.002,.005,.01]; %kON Value
ParamArray = [.1,1,5,10,15]; %alpha Value
%ParamArray = [.0001,.0005,.001,.005,.01]; %kUb Value
%ParamArray = [.1,.5,1,5,10]; %kOFF Value

GeneArray = [5,10,15,20,25,30,35,40,45,50];
%GeneArray = [10,20,30,40,50,60,70,80,90,100];
RunsArray = ones(50,1);%[20,20,20,20,20,20,20,20,20,20];%[120,60,40,30,24,20,17,14,13,12];
RunsArray = RunsArray .*30;

for n = 1:length(ParamArray)
    
    

    kB =  .001;
    kUb = 0.0001;
    kON = 0.01;
    kOFF =  10;
    alpha = ParamArray(n);

    FullTrials = 10;

    

    AllRibosomesPerGene = zeros(3,max(RunsArray),max(GeneArray));
    for k = 1:FullTrials
        Case = k;
        disp(Case)
        NumGenes = GeneArray(Case);
        RibosomesInitial = 300*NumGenes;
        Runs = RunsArray(k);

        
        tMax = 100;
        dt = 1;
        tspan = 0:dt:tMax;
        tspan2 = length(tspan);

        MaxOutput = 2000; %Maximum size expected of the output file
        BurstDurationTrack = zeros(100000,Runs);
        BurstTimesTrack = zeros(100000,Runs);
        a = zeros(NumGenes*2,1);
        TraceVar = zeros(Runs,NumGenes);
        TraceAuto = zeros(tspan2,Runs*NumGenes);
        TraceAutoNorm = TraceAuto;

        Genes = zeros(1,NumGenes*2);
        Genes(NumGenes+1:end) = 1;
        BoundRibosomeArray = zeros(MaxOutput,1);
        mRNATotals = zeros(1,NumGenes);
        aGenes = zeros(NumGenes*2,1);

        GeneRuns = zeros(Runs,tspan2,NumGenes);
        RibosomesPerGene = zeros(Runs,NumGenes);

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
        aAlpha = zeros(1,NumGenes);
        %

        for i = 1:Runs

           
            %T = zeros(MaxOutput, 1); %Time tracking
            %X = zeros(1, NumSpecies); %Species number tracking
            %T(1)     = 0;
            %X(1,:)   = x0;
            RxnCount = 1;
            T = 0; %Time
            mRNAArray = zeros(MaxOutput,2);
            Ttrack = zeros(length(tspan), 1); %Time tracking
            XGenes = zeros(length(tspan), NumGenes*2); %Species number tracking
            XmRNA = zeros(length(tspan),length(mRNAArray(:,1)));
            XRibosomes = XmRNA;
            Ribosomes = RibosomesInitial;
            XGenes(1,:)   = Genes;
            xCurrentGenes = Genes;
            xCurrentmRNA = zeros(MaxOutput,1);
            xCurrentRibosomes = zeros(MaxOutput,1);
            mRNATotals = zeros(1,NumGenes);

            RecordTime = dt; %Recording time
            RecordCount = 2;
            OnTrack = 0;
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

                elseif r(2)*a0 < sum(aGenes) + sum(aAlpha)
                    mu  = find((sum(aGenes) + cumsum(aAlpha) >= r(2)*a0),1,'first');
                    T   = T  + tau;
                    mRNATotals(mu) = mRNATotals(mu) + 1;
                    mRNAnum = mRNAnum + 1;
                    xCurrentmRNA(mRNAnum,1) = 1;
                    mRNAArray(mRNAnum,1) = 1;
                    mRNAArray(mRNAnum,2) = mu;

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

             mRNAIdx = sum(mRNATotals);
             for j = 1:mRNAIdx
                RibosomesPerGene(i,mRNAArray(j,2)) = RibosomesPerGene(i,mRNAArray(j,2)) + xCurrentRibosomes(j);
             end

            % Record output
            count = count - 1;

            GeneRuns(i,:,:) = XGenes(1:tspan2,1:NumGenes);
            mRNARuns(i,:,:) = XmRNA(1:tspan2,:);
            RibosomeRuns(i,:,:) = XRibosomes(1:tspan2,:);



        end

       AllRibosomeRuns(k,:,:,:) = RibosomeRuns;
       AllRibosomesPerGene(k,1:Runs,1:NumGenes) = RibosomesPerGene;

    end
    %save AllData
    %%
    %Calculations
    for k = 1:FullTrials
        for j = 1:Runs
            temp = AllRibosomeRuns(k,j,:,:);
            temp = reshape(temp,[tspan2,MaxOutput]);
            tempsum = sum(temp,1);
            [RibosomeIdx, misc] = find(tempsum);
            for i = 1:length(RibosomeIdx)
                EndRibosomeCount(k,j,i) = temp(end,i);

            end
            tempsum = sum(temp,2);
            AllTotalRibosomes(k,j,:) = tempsum;
            EndmRNAwRibosomes(k,j) = length(RibosomeIdx); 
        end


        temp = mean(AllTotalRibosomes(k,:,:),2);
        TotalRibosomesGenTrend(k,:) = temp(:);


        for j = 1:Runs
            TotalRibosomeSS(n,k,j) = AllTotalRibosomes(k,j,end);
            temp = AllTotalRibosomes(k,j,:);
            temp2 = TotalRibosomesGenTrend(k,:);
            TotalRibosomeA(k,j,:) = temp(:) - temp2(:);
            TotalRibosomeVar(n,k,j) = var(TotalRibosomeA(k,j,:));
            TotalRibosomecv2(n,k,j) =  TotalRibosomeVar(n,k,j)/TotalRibosomeSS(n,k,j)^2;
            if isinf(TotalRibosomecv2(n,k,j))
                TotalRibosomecv2(n,k,j) = NaN;
            end

        end
        
        for j = 1:Runs
            temp = AllRibosomesPerGene(k,j,1:GeneArray(k));
            temp = temp(:);
            temp(temp == 0) = [];
            TotalRibosomeBS(n,k,j) = mean(temp);
            TotalRibosomeBF(n,k,j) = length(temp);
            
        end
    end
    
end
name = sprintf('Data%s',testParam);
save(name,'TotalRibosomeSS','TotalRibosomeVar','TotalRibosomecv2','GeneArray',...
    'RibosomesInitial','kB','kON','kOFF','alpha');
%%
%plot cv2 vs abundance
hold on
c = colormap(jet(length(ParamArray)));
for n = 1:length(ParamArray)
    for j = 1:length(GeneArray)
        tempSS = TotalRibosomeSS(n,j,:);
        tempcv2 = TotalRibosomecv2(n,j,:);
        plot(tempSS(:),tempcv2(:),'color',c(n,:),...
            'marker','.','markersize',6, 'linestyle','none')
    end
   
end
tempSS = 0;
tempcv2 = 0;
temp2SS = 0;
temp2cv2 = 0;
for n = 1:length(ParamArray)
    name = sprintf('%g %s',ParamArray(n),testParam);
    for j = 1:FullTrials
        tempSS = TotalRibosomeSS(n,j,:);
        tempSS = reshape(tempSS,[1,Runs]);
        tempSS(tempSS == 0) = NaN;
        s = isnan(tempSS);
        temp2SS(j) = geomean(tempSS(~s));
        tempcv2 = TotalRibosomecv2(n,j,:);
        tempcv2 = reshape(tempcv2,[1,Runs]);
        tempcv2(tempcv2 == Inf) = NaN;
        s = isnan(tempcv2);
        temp2cv2(j) = geomean(tempcv2(~s));
    end
    linestore(n) = plot(temp2SS(:),temp2cv2(:),'color',c(n,:),...
        'marker','o','markersize',6, 'markerfacecolor',c(n,:),...
        'markeredgecolor','k','displayname',name);
end
xrange = 100:100000;
yrange = 100./xrange;
plot(xrange,yrange,'color','k')
legend(linestore,'location','northeast')
set(gca,'XScale','log');
set(gca,'YScale','log');
axis([10 100000 10^-4 100]);
xlabel('Ribosome Abundance','FontSize',15)
ylabel('Ribosome cv^2','FontSize',15)
title('Total Ribosome cv^2 vs Abundance')
name = sprintf('cv2AbundanceVary%s.jpg',testParam);
saveas(gcf,name)
%%

figure
hold on
c = colormap(jet(length(ParamArray)));
for n = 1:length(ParamArray)
    for j = 1:length(GeneArray)
        tempBS = TotalRibosomeBS(n,j,:);
        tempBF = TotalRibosomeBF(n,j,:);
        plot(tempBF(:),tempBS(:),'color',c(n,:),...
            'marker','.','markersize',6, 'linestyle','none')
    end
   
end
temp2BS = 0;
temp2BF = 0;
for n = 1:length(ParamArray)
    name = sprintf('%g %s',ParamArray(n),testParam);
      for j = 1:FullTrials
            tempBS = TotalRibosomeBS(n,j,:);
            tempBS = reshape(tempBS,[1,Runs]);
            tempBS(tempBS == 0) = NaN;
            s = isnan(tempBS);
            temp2BS(j) = mean(tempBS(~s));
            tempBF = TotalRibosomeBF(n,j,:);
            tempBF = reshape(tempBF,[1,Runs]);
            tempBF(tempBF == Inf) = NaN;
            s = isnan(tempBF);
            temp2BF(j) = mean(tempBF(~s));
       end
        linestore(n) = plot(temp2BF(:),temp2BS(:),'color',c(n,:),...
        'marker','o','markersize',6, 'markerfacecolor',c(n,:),...
        'markeredgecolor','k','displayname',name);
end
% set(gca,'XScale','log');
% set(gca,'YScale','log');
legend(linestore,'location','northeast')
ylabel('Mean Intensity ON Gene(BS)','FontSize',15)
xlabel('Mean Number ON (BF)','FontSize',15)
title('Burst Frequency Total Ribosome Per Gene ')
name = sprintf('BFvsBSVary%s.jpg',testParam);
saveas(gcf,name)
% %%
% %plot Total Ribosomes over time
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
% %%
% %plot cv2 vs abundance
% figure
% hold on
% c = colormap(jet(FullTrials));
% for k = 1:FullTrials
%     
%     plot(TotalRibosomeSS(k,:),TotalRibosomecv2(k,:),'color',c(k,:),...
%         'linestyle','none','marker','.','markersize',10)
% 
% end
% 
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% axis([.1 100000 10^-6 1]);
% xlabel('Ribosome Abundance','FontSize',15)
% ylabel('Ribosome cv^2','FontSize',15)
% title('Total Ribosome cv^2 vs Abundance')
% saveas(gcf,'Ribsomescv2abundance.jpg')
% %%
% figure
% hold on
% c = colormap(jet(FullTrials));
% for k = 1:FullTrials
%     for j = 1:Runs
%         temp = AllRibosomeRuns(k,j,:,:);
%         temp = reshape(temp,[tspan2,MaxOutput]);
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
% %%
% figure
% hold on
% for k = 1:FullTrials
%     for j = 1:Runs
%         temp = AllRibosomesPerGene(k,j,1:GeneArray(k));
%         temp = temp(:);
%         plot(GeneArray(k),temp,'linestyle','none',...
%             'marker','.','color',c(k,:),'markersize',10)
%     end
% end
% xlabel('Concurrent Genes','FontSize',15)
% ylabel('Ribosomes per Gene','FontSize',15)
% title('Number of Ribsomes per Gene')
% saveas(gcf,'RibsomesperGene.jpg')
% %%
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

% BurstTimesBetweenTrack = diff(BurstTimesTrack,1);
% % GillespieData = GillespieData(200:end,:);
% len = length(GillespieData(:,1));
% GillespieDataAvg = zeros(len,1); %Average curve of all Runs
% mRNADataAvg = zeros(len,1); %Average curve of all Runs
% for h = 1:len
%     GillespieDataAvg(h) = sum(GillespieData(h,:))/Runs;
%     mRNADataAvg(h) = sum(mRNAData(h,:))/Runs;
% end
%     
% GillespieDataA = zeros(len,Runs);
% mRNADataA = zeros(len,Runs);
% 
% for k = 1:Runs
%     disp(k)
%     SSlevelsArray(k) = mean(GillespieData(:,k));
%     SSmRNAArray2(k) = mean(mRNAData(:,k));
%     
%     GillespieDataA(:,k) = GillespieData(:,k) - GillespieDataAvg;
%     mRNADataA(:,k) = mRNAData(:,k) - mRNADataAvg;
%     
%     AutoArrayTemp = xcorr(GillespieDataA(:,k),'unbiased');%,length(tspan2)-1);
%     AutomRNAArrayTemp = xcorr(mRNADataA(:,k),'unbiased');%,length(tspan2)-1);
%     AutoArray(:,k) = AutoArrayTemp;
%     VarienceArray(k) = AutoArrayTemp(len);
%     Avgcv2(k) = AutoArrayTemp(len)/(SSlevelsArray(k))^2;
%     AutomRNAArray(:,k) = AutomRNAArrayTemp;
%     VariencemRNAArray(k) = AutomRNAArrayTemp(len);
%     AvgmRNAcv2(k) = AutomRNAArrayTemp(len)/(SSmRNAArray2(k))^2;
% 
% end
%%

% 
% % 
% % AvgBurstFreq(:) = BurstFreq;
% % AvgBurstSize(:) = BurstSize;
% %%
% % c = colormap(hsv(length(kONArray)));
% % %c = colormap(hsv(length(VarArray)));
% % for i = 1:length(kONArray)
% %     hold on
% %     plot(SSlevelsTot(:,i),Avgcv2Tot(:,i),'linestyle','none','marker','.',...
% %         'markersize',10,'color',c(i,:));
% % end
% % set(gca,'XScale','log');
% % set(gca,'YScale','log');
% save 2StateEncounterData
% %%
% 
% %c = colormap(hsv(length(kONArray)));
% hold on
% for i = 1:7
%     plot(SSmRNAArray(:,i),cv2mRNAArray(:,i),'linestyle','none','marker','.',...
%         'markersize',8,'color','k');
% end
% 
% c = colormap(hsv(length(VarArray)));
% for i = 1:length(VarArray)
%     hold on
%     plot(SSmRNATot(:,i),AvgmRNAcv2Tot(:,i),'linestyle','none','marker','.',...
%         'markersize',10,'color',c(i,:));
% end
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% axis([7 30 .03 1])
% xlabel('mRNA Abundace','FontSize',15)
% ylabel('cv2','FontSize',15)
% title('mRNA cv2 v Abundance From Encounters')
% saveas(gcf,'mRNAcv2vAbundance2StateDirect.jpg')

