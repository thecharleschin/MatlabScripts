%base Gillespie Algorithm for general purpose simulations.

clear all
close all
rand('state',sum(100*clock)); %#ok<RAND>

ParamArray = .01;%logspace(-3,1,10);%[.01,.1,1,10];

name = sprintf('AllDataPower3kP');

for k = 1:length(ParamArray)
%Parameters
RibosomesInitial = 300;
Runs = 200;

kON = .02;%gene burst on rate
kOFF = 1; %gene burst off rate
alpha = 2; %mRNA production rate
gammam = log(2)/5; %mRNA decay rate, 5min halflife
kIN = ParamArray(k); %rate of coming into local ribosome pool
kOUT = 1; %rate of leaving local ribosome pool
kP = .1; %protein Production rate
gammap = log(2)/20;
power = 3;

tMax = 200;
startTime = 1;
dt = 1;
tspan = 0:dt:tMax;
tspan2 = length(tspan);
MaxOutput = 1000;


%GeneON GeneOFF mRNA protein  BoundRibosomes GlobalRibosome
x0 = [0,1,0,0,0,RibosomesInitial];

RxnMatrix = [ 1 -1  0  0  0  0; %Burst ON
             -1  1  0  0  0  0; %Burst OFF
              0  0  1  0  0  0; %transcription
              0  0 -1  0  0  0; %mRNA decay
              0  0  0  1  1 -1; %ribosome in
              0  0  0  0 -1  1; %ribosome out
              0  0  0  1  0  0; %translation
              0  0  0 -1  0  0]; %protein decay


          
parfor i = 1:Runs
    disp(i)
    p = 1;
    Onflag = x0(3);
    mRNAperBurst = zeros(100000,1);

    NumSpecies = size(RxnMatrix, 2);
    RxnCount = 1;
    T = 0; %Time
    Ttrack = zeros(length(tspan), 1); %Time tracking
    X = zeros(length(tspan), NumSpecies); %Species number tracking
    X(1,:)   = x0;
    xCurrent = x0;
    RecordTime = dt; %Recording time
    RecordCount = 2;
    OnTrack = 0;
    count = 1;
    BurstIdx = 0;
    BurstArray = zeros(MaxOutput,3);
    BurstStartCount = BurstIdx;
        StartCount = 0;

    %%%%%%%%%%%%%%%%%%%%%%
    %Gillespie Simulation%
    %%%%%%%%%%%%%%%%%%%%%%

    while T <= tMax


        x = xCurrent;

        StateON = x(1);
        StateOFF = x(2);
        mRNA    = x(3);
        protein = x(4);
        boundR = x(5);
        globalR = x(6);

        % Calculate reaction propensities
        a = [kON*StateOFF;
             kOFF*StateON;
             alpha*StateON; 
             gammam*mRNA;
             (mRNA)*globalR*kIN;
             boundR*kOUT;
             (mRNA^power)*kP*boundR;
             gammap*protein];

        % Compute tau and mu using random variables
        a0 = sum(a);
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
        
        %[~, mu] = histc(r(2)*a0, [0;cumsum(a(:))]);
        mu  = find((cumsum(a) >= r(2)*a0),1,'first');

        %find the next change in state before tau
        T   = T  + tau;
        xCurrent = xCurrent + RxnMatrix(mu,:);
        
        if mu == 1 %if burst on
            BurstIdx = BurstIdx + 1;
            BurstArray(BurstIdx,3) = T; 
        elseif mu == 3 %if make mRNA
            BurstArray(BurstIdx,1) = BurstArray(BurstIdx,1) +1;
        elseif mu == 7 %if make Protein
            BurstArray(BurstIdx,2) = BurstArray(BurstIdx,2) +1;
        end
        if T >= startTime/dt && StartCount == 0
             BurstStartCount = BurstIdx;
             StartCount = 1;
        end


    end


    % Record output
     X(RecordCount,:) = xCurrent;
     AllRuns(:,:,i,k) = X(1:tMax/dt,:);
     AllBurstArray(:,:,i,k) = BurstArray;
     AllBurstStartCount(i,k) = BurstStartCount;
end


%%
%Pull out burst information

for i = 1:Runs
    Bursttemp = AllBurstArray(:,:,i,k);
    BurstPtemp = zeros(MaxOutput,1);
    BurstMtemp = zeros(MaxOutput,1);
    

    BurstPtemp(:) = Bursttemp(:,2)./Bursttemp(:,1);
    BurstMtemp(:) = Bursttemp(:,1);

    if AllBurstStartCount(i,k) <= 1
    else
        BurstPtemp(1:AllBurstStartCount(i,k)-1) = NaN;
        BurstMtemp(1:AllBurstStartCount(i,k)-1) = NaN;
    end

    
    BurstMtemp(BurstMtemp == 0) = NaN;
    BurstPtemp(BurstPtemp == Inf) = NaN;
    
    BurstM(i,k) = nanmean(BurstMtemp(:));
    BurstP(i,k) = nanmean(BurstPtemp(:));
    BurstMAll(:,i,k) = BurstMtemp(:);
    BurstPAll(:,i,k) = BurstPtemp(:);
    

    
end
%%
binInterval = 1:10;
for j = 1:length(binInterval)
   temp = BurstMAll(:,:,k);
   temp2 = find(temp(:) == binInterval(j)); 
   temp = BurstPAll(:,:,k);
   tempP = temp(temp2);
   binP(j,k) = nanmean(tempP);

end
%%
end

save(name)
% %%
% hold on
% xtemp = .1:.1:4;
% ytemp = 5.*xtemp.^2;
% plot(xtemp,ytemp,'k')
% 
% xtemp = .1:.1:4;
% ytemp = 5.*xtemp.^3;
% plot(xtemp,ytemp,'b')
% xtemp = .1:.1:4;
% ytemp = 5.*xtemp.^1;
% plot(xtemp,ytemp,'r')

%%
hold on
c = colormap(jet(length(ParamArray)));
for k = 1:length(ParamArray)
plot(BurstM(:,k),BurstP(:,k),'linestyle','none','marker','.','markersize',10,...
    'color','b')
end
set(gca,'YScale','log');
axis([0 Inf .1 Inf])
xlabel('Burst Size (mRNA per Burst)','FontSize',15)
ylabel('Burst Size (Protein per mRNA)','FontSize',15)
name = sprintf('transcriptional vs translational BS');
title(name)
name = sprintf('BSmeanPlotPower1kP.jpg');
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
name = sprintf('BSPlotPower3kP.jpg');
saveas(gcf,name)

%%
subplot(3,1,1);
hold on
plot(AllRuns(:,3,1,k),'color','b')
plot(AllRuns(:,5,1,k),'color','r')
subplot(3,1,2);
plot(AllRuns(:,6,1,k),'color','b')
subplot(3,1,3);
plot(AllRuns(:,4,1,k),'color','g')
