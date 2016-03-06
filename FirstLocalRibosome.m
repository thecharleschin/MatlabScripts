clear all
close all
rand('state',sum(100*clock)); %ok<RAND>
MaxOutput = 1000; %Maximum size expected of the output file


NumGenes = 1;
RibosomesInitial = 300*NumGenes;
Runs = 200;

kON = .02;%gene burst on rate
kOFF = 1; %gene burst off rate
alpha = 2; %mRNA production rate
gammam = log(2)/5; %mRNA decay rate, 5min halflife
kIN = .1; %rate of coming into local ribosome pool
kOUT = 100; %rate of leaving local ribosome pool
kP = 1; %protein Production rate
gammap = log(2)/20;

tMax = 200;
startTime = 1;
dt = 1;
tspan = startTime:dt:tMax;
tspan2 = length(tspan);




Genes = zeros(1,NumGenes*2);
Genes(NumGenes+1:end) = 1;



GeneRuns = zeros(Runs,tspan2,NumGenes);

BurstArray = zeros(MaxOutput,3,NumGenes);
BurstIdx = zeros(NumGenes,1);


%Build Rxn Matrix and Propensity for genes and local pools
%#GenesON #GenesOFF #mRNAPools #LocalRibosomePools #Protein GlobalRibosomePool
%NumGenes NumGenes NumGenes NumGenes NumGenes 1
%#AllBurstON #AllBurstOFF #AllTranscription #AllmRNADecay #AllLocalIN #AllLocalOUT #AllTranslation #AllProteinDecay
%NumGenes NumGenes NumGenes NumGenes NumGenes NumGenes NumGenes NumGenes

RxnMatrix = zeros(NumGenes*8,NumGenes*5+1);

for i = 1:NumGenes
    RxnMatrix(i,i) = 1; %BurstON ONstate
    RxnMatrix(i,i+NumGenes) = -1; %BurstON Offstate
    RxnMatrix(i+NumGenes,i) = -1; %BurstOFF ONstate
    RxnMatrix(i+NumGenes,i+NumGenes) = 1; %BurstOFF OffState
    RxnMatrix(i+NumGenes*2,i+NumGenes*2) = 1; %Transcription
    RxnMatrix(i+NumGenes*3,i+NumGenes*2) = -1; %mRNADecay
    RxnMatrix(i+NumGenes*4,i+NumGenes*3) = 1; %LocalIN
    RxnMatrix(i+NumGenes*5,i+NumGenes*3) = -1; %LocalOUT
    RxnMatrix(i+NumGenes*6,i+NumGenes*4) = 1; %Translation
    RxnMatrix(i+NumGenes*7,i+NumGenes*4) = -1; %ProteinDecay 
   
end

RxnMatrix(NumGenes*4+1:NumGenes*5,NumGenes*5+1) = -1; %LocalIN
RxnMatrix(NumGenes*5+1:NumGenes*6,NumGenes*5+1) = 1; %LocalOUT



x0 = zeros(1,NumGenes*5+1);
x0(1+NumGenes:NumGenes*2) = 1;
x0(NumGenes*5+1) = RibosomesInitial;
aAll = zeros(NumGenes*8,1);
%

parfor i = 1:Runs

    disp('run')
    disp(i)
    
    RxnCount = 1;
    T = 0; %Time
    
    X = zeros(length(tspan), NumGenes*5+1); %Species number tracking
    X(1,:)   = x0;
    xCurrent = x0;

    aAll = zeros(NumGenes*8,1);
    
    RecordTime = dt; %Recording time
    RecordCount = 2;
    mRNAnum = 0;
    StartCount = 0;
    BurstIdx = zeros(NumGenes,1);
    BurstArray = zeros(MaxOutput,3,NumGenes);
    BurstStartCount = BurstIdx;
    count = 1;
    OnDuration = zeros(MaxOutput,1);
    OnTimes = zeros(MaxOutput,1);
    %%%%%%%%%%%%%%%%%%%%%%
    %Gillespie Simulation%
    %%%%%%%%%%%%%%%%%%%%%%

    while T <= tMax

        % Calculate reaction propensities
        x = xCurrent;
        StateON = x(1:NumGenes);
        StateOFF = x(NumGenes+1:NumGenes*2);
        mRNA = x(NumGenes*2+1:NumGenes*3);
        Local = x(NumGenes*3+1:NumGenes*4);
        Protein = x(NumGenes*4+1:NumGenes*5);
        Global = x(NumGenes*5+1);
        
        for j = 1:NumGenes
            aAll(j) = kON*StateOFF(j); %Burst ON
            aAll(j+NumGenes) = kOFF*StateON(j); %Burst OFF
            aAll(j+NumGenes*2) = alpha*StateON(j); %Transcription
            aAll(j+NumGenes*3) = gammam*mRNA(j); %mRNADecay
            aAll(j+NumGenes*4) = Global*kIN;%*mRNA(j);% LocalIN
%             if mRNA(j) == 0
%                 aAll(j+NumGenes*5) = Local(j)*kOUT; 
%             else
            aAll(j+NumGenes*5) = Local(j)*kOUT/(mRNA(j)+1);%LocalOUT
%             end
            aAll(j+NumGenes*6) = kP*mRNA(j)*Local(j); %Translation
            aAll(j+NumGenes*7) = gammap*Protein(j);    %ProteinDecay
        end
        
        

        % Compute tau and mu using random variables
        a0 = sum(aAll);
        r = rand(1,2);
        %tau = -log(r(1))/a0;
        tau = (1/a0)*log(1/r(1));

        %Store information if at time
        tau0 = tau;
        if tau0 > dt
            while tau0 > dt
                X(RecordCount,:) = xCurrent;
                RecordCount = RecordCount + 1;
                RecordTime = RecordTime + dt;
                tau0 = tau0 - dt;
            end
        end

        if T + tau > RecordTime
            X(RecordCount,:) = xCurrent;
            RecordCount = RecordCount + 1;
            RecordTime = RecordTime + dt;
        end
        
        mu  = find((cumsum(aAll) >= r(2)*a0),1,'first');
        T   = T  + tau;
        xCurrent = xCurrent + RxnMatrix(mu,:);
        
        
        if mu <= NumGenes %if burst on
            BurstIdx(mu) = BurstIdx(mu) + 1;
            BurstArray(BurstIdx(mu),3) = T; 
        elseif mu > NumGenes*2 && mu <= NumGenes*3 %if make mRNA
            BurstArray(BurstIdx(mu - NumGenes*2),1,mu - NumGenes*2) = ...
                BurstArray(BurstIdx(mu - NumGenes*2),1,mu - NumGenes*2) +1;
        elseif mu > NumGenes*6 && mu <= NumGenes*7 %if make Protein
            BurstArray(BurstIdx(mu - NumGenes*6),2,mu - NumGenes*6) = ...
                BurstArray(BurstIdx(mu - NumGenes*6),2,mu - NumGenes*6) +1;
        end
        if T >= startTime/dt && StartCount == 0
             BurstStartCount = BurstIdx;
             StartCount = 1;
        end

    end


    % Record output
     X(RecordCount,:) = xCurrent;
   
    
     AllRuns(:,:,i) = X(1:tMax/dt,:);
     AllBurstArray(:,:,:,i) = BurstArray;
     AllBurstIdx(:,i) = BurstIdx;
     AllGene(:,:,i) = X(1:tMax/dt,1:NumGenes*2);
     AllmRNA(:,:,i) = X(1:tMax/dt,NumGenes*2+1:NumGenes*3);
     AllLocal(:,:,i) = X(1:tMax/dt,NumGenes*3+1:NumGenes*4);
     AllProtein(:,:,i) = X(1:tMax/dt,NumGenes*4+1:NumGenes*5);
     AllBurstStartCount(:,i) = BurstStartCount;
    
end

%%
%Pull out burst information
for i = 1:Runs
    Bursttemp = AllBurstArray(:,:,:,i);
    BurstPtemp = zeros(MaxOutput,NumGenes);
    BurstMtemp = zeros(MaxOutput,NumGenes);
    
    for j = 1:NumGenes
        BurstPtemp(:,j) = Bursttemp(:,2,j)./Bursttemp(:,1,j);
        BurstMtemp(:,j) = Bursttemp(:,1,j);
        
        if AllBurstStartCount(j,i) <= 1
        else
            BurstPtemp(1:AllBurstStartCount(j,i)-1,j) = NaN;
            BurstMtemp(1:AllBurstStartCount(j,i)-1,j) = NaN;
        end
        
    end
    BurstMtemp(BurstMtemp == 0) = NaN;
    BurstPtemp(BurstPtemp == Inf) = NaN;
    
    BurstM(i) = nanmean(BurstMtemp(:));
    BurstP(i) = nanmean(BurstPtemp(:));
    BurstMAll(:,i) = BurstMtemp(:);
    BurstPAll(:,i) = BurstPtemp(:);
    
end

save AllData
%%
hold on
xtemp = .1:.1:4;
ytemp = 5.*xtemp.^2;
plot(xtemp,ytemp,'k')

xtemp = .1:.1:4;
ytemp = 5.*xtemp.^3;
plot(xtemp,ytemp,'b')
xtemp = .1:.1:4;
ytemp = 5.*xtemp.^1;
plot(xtemp,ytemp,'r')

plot(BurstM,BurstP,'linestyle','none','marker','.','markersize',10)
set(gca,'YScale','log');
axis([0 Inf .1 Inf])
xlabel('Burst Size (mRNA per Burst)','FontSize',15)
ylabel('Burst Size (Protein per mRNA)','FontSize',15)
name = sprintf('transcriptional vs translational BS');
title(name)
name = sprintf('BSPlot.jpg');
saveas(gcf,name)
%%
figure
hold on
%bin data to take average
binInterval = 0:1:10;
burstMbin = BurstMAll(:);
burstPbin = BurstPAll(:);
burstMbin(burstMbin == 0) = NaN;
%burstR(isnan(burstR)) = [];
binM = zeros(length(binInterval),1);
binP = binM;
binPnum = binP;

for i = 2:length(binInterval)
	for j = 1:length(burstPbin(:))
		if burstMbin(j) >= binInterval(i-1) && burstMbin(j) < binInterval(i)
			binM(i-1) = burstMbin(j); 
            if isnan(burstPbin(j))
            else
                binP(i-1) = binP(i-1) + burstPbin(j);
                binPnum(i-1) = binPnum(i-1) + 1;
            end
        end
    end
end
binP = binP ./ binPnum;
plot(BurstMAll(:),BurstPAll(:),'linestyle','none','marker','.','markersize',10)
plot(binInterval,binP,'linestyle','none','marker','o','markersize',8,...
    'markerfacecolor','b','markeredgecolor','k')
set(gca,'YScale','log');
axis([0 Inf .1 Inf])
xlabel('Burst Size (mRNA per Burst)','FontSize',15)
ylabel('Burst Size (Protein per mRNA)','FontSize',15)
name = sprintf('transcriptional vs translational BS');
title(name)
name = sprintf('BSPlot2.jpg');
saveas(gcf,name)
