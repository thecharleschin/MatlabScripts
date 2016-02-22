clear all
close all
rand('state',sum(100*clock)); %ok<RAND>
 MaxOutput = 1000; %Maximum size expected of the output file


NumGenes = 30;
RibosomesInitial = 300*NumGenes;
Runs = 100;

kB = .001;
kON = .001;
kOFF = 1;
alpha = 5; %20 seconds a transcript, 3 per min
gammam = log(2)/2;

tMax = 600;
startTime = 101;
dt = 1;
tspan = startTime:dt:tMax;
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
BurstSizeRuns = zeros(Runs,MaxOutput,NumGenes);
BurstSizeRRuns = zeros(Runs,MaxOutput,NumGenes);
BurstTimesRuns = zeros(Runs,MaxOutput,NumGenes);
mRNATotalsRuns = zeros(Runs,NumGenes);
RibosomesPerGene = zeros(Runs,NumGenes);

BurstDurationTrack = zeros(MaxOutput,Runs);

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

    disp('run')
    disp(i)
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
    BurstTrack = zeros(1,NumGenes); %track burst index
    BurstSizeTrack = zeros(MaxOutput,NumGenes);%brst size track
    BurstTimes = zeros(MaxOutput,NumGenes);%burst time track
    BurstRTrack = zeros(MaxOutput,NumGenes); %track burst index Ribosomes
    BurstSizeRTrack = zeros(MaxOutput,NumGenes);%brst size track
    BurstRTimes = zeros(MaxOutput,NumGenes);%burst time track
    OnTrack = zeros(1,NumGenes); %Burst state track
    UniqueBurst = 0;

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
        aGamma = xCurrentmRNA(:,1) .* gammam;


        % Compute tau and mu using random variables
        a0 = sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + sum(aGamma);
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
                UniqueBurst = UniqueBurst + 1;
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
            mRNAArray(mRNAnum,3) = BurstTrack(mu);
            mRNAArray(mRNAnum,4) = UniqueBurst;            
            mRNAArray(mRNAnum,5) = T;
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
            mu  = find((sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + cumsum(aGamma) >= r(2)*a0),1,'first');
            %find the next change in state before tau
            T   = T  + tau;
            %xCurrentGenes = xCurrentGenes + RxnMatrixGene(mu,:);
            mRNATotals(mRNAArray(mu,2)) = mRNATotals(mRNAArray(mu,2)) - 1;
            %mRNAnum = mRNAnum - 1;
            xCurrentmRNA(mu) = 0;
            %xCurrentmRNA(end+1) = 0;
            mRNAArray(mu,1) = 0;
            %mRNAArray(end+1,:) = 0;
%                 BurstRTrack(mu) = BurstRTrack(mu)  
%                 if BurstSizeRTrack(BurstRTrack(mu),mu) == 0
%                     BurstRTimes(BurstRTrack(mu),mu) = T;
%                 end
%                 BurstSizeRTrack(BurstRTrack(mu),mu) = BurstSizeRTrack(BurstRTrack(mu),mu) + 1;
%                 
            UnbindingRibosomes = xCurrentRibosomes(mu);
            mRNAArray(mu,6) = UnbindingRibosomes;
            Ribosomes = Ribosomes + UnbindingRibosomes;
            BurstRTrack(mRNAArray(mu,2)) = BurstRTrack(mRNAArray(mu,2))+1;
            BurstSizeRTrack(BurstRTrack(mRNAArray(mu,2)),mRNAArray(mu,2)) = UnbindingRibosomes;
            xCurrentRibosomes(mu) = 0;
            %xCurrentRibosomes(end+1) = 0;

        end



    end


    % Record output
     XGenes(RecordCount,:) = xCurrentGenes;
     XmRNA(RecordCount,:) = xCurrentmRNA(:,1);
     XRibosomes(RecordCount,:) = xCurrentRibosomes;

     mRNATotnum(i) = sum(mRNATotals);
%          for j = 1:mRNATotnum(i)
%             RibosomesPerGene(i,mRNAArray(j,2)) = RibosomesPerGene(i,mRNAArray(j,2)) + xCurrentRibosomes(j);
%          end

    % Record output
    count = count - 1;
    mRNATotalsRuns(i,:) = mRNATotals;
    GeneRuns(i,:,:) = XGenes(startTime:tspan(end),1:NumGenes);
    mRNARuns(i,:,:) = XmRNA(startTime:tspan(end),:);
    RibosomeRuns(i,:,:) = XRibosomes(startTime:tspan(end),:);
    AllmRNAArray(:,:,i) = mRNAArray; 
    BurstTrackRuns(i,:) = BurstTrack;
    BurstSizeRuns(i,:,:) = BurstSizeTrack;
    BurstSizeRRuns(i,:,:) = BurstSizeRTrack;
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
   
save AllData
%%

burstM = zeros(1000,Runs);
burstR = zeros(1000,Runs);

for i = 1:Runs
    temp = AllmRNAArray(:,:,i);
    [Idx, ~] = find(temp(:,5));
    temp = temp(1:length(Idx),:);
	%dont count mRNA that have not decayed
	for j = 1:length(temp(:,1))
		if temp(j,1) == 1
			temp(j,:) = NaN;
		end
	end
	%temp(isnan(temp)) = [];		
	
% 	uniqueBursts = accumarray([temp(:,2),temp(:,3)],1);
% 	pairs = unique([temp(:,2),temp(:,3)],'rows');
	
    count = 1;
    burstIdx1 = temp(1,2);
    burstIdx2 = temp(1,3);
    burstM(count,i) = 1;
    burstR(count,i) = temp(1,6);
    j = 2;
    while j <= length(temp)
        if temp(j,2) == burstIdx1 && temp(j,3) == burstIdx2
            burstM(count,i) = burstM(count,i) + 1;
            burstR(count,i) = burstR(count,i) + temp(j,6);
        else
            count = count + 1;
            burstM(count,i) = burstM(count,i) + 1;
            burstR(count,i) = burstR(count,i) + temp(j,6);
            burstIdx1 = temp(j,2);
            burstIdx2 = temp(j,3);
        end
        
        j = j+ 1;
    end
    burstR(:,i) = burstR(:,i)./burstM(:,i);
end
%%
plot(burstM,burstR,'linestyle','none','marker','.','markersize',8)
set(gca,'YScale','log');

%%
%bin data to take average
binInterval = 0:2:40;
burstMbin = burstM(:);
burstRbin = burstR(:);
burstMbin(burstMbin == 0) = NaN;
%burstR(isnan(burstR)) = [];
binM = zeros(length(binInterval),1);
binR = binM;

for i = 2:length(binInterval)
	for j = 1:length(burstRbin(:))
		if burstMbin(j) >= binInterval(i-1) && burstMbin(j) < binInterval(i)
			binM(i) = burstMbin(j); 
            binR(i) = burstRbin(j);
        end
    end
end

hold on
plot(burstM,burstR,'linestyle','none','marker','.','markersize',10)
set(gca,'YScale','log');
c = colormap(jet(length(binInterval)));
plot(binInterval,binR,'linestyle','none','marker','o','markersize',8,'markerfacecolor','b','markeredgecolor','k')
xlabel('Burst Size (mRNA per Burst)','FontSize',15)
ylabel('Burst Size (Ribosomes per mRNA)','FontSize',15)
title('Ribosomes per mRNA vs mRNA per Burst')
saveas(gcf,'BurstSizemRNARibosome.jpg')