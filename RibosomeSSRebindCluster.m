function RibosomeSSRebindCluster(jobNum)
%Cluster Variables

%Limit number of cpus this program can utilize
maxNumCompThreads(1);

% seeds = load('seeds.mat');
% seeds = seeds.seeds;
% 
% %Seed random number generator
% rng(seeds(jobNum));

MaxOutput = 1000; %Maximum size expected of the output file

tic

NumGenes = 10;
RibosomesInitial = 200*NumGenes;

kB = .001; %ribosome binding rate from initial pool
kON = .01;%gene burst on rate
kOFF = 1; %gene burst off rate
alpha = 2; %mRNA production rate
gammam = log(2)/10; %mRNA decay rate, 5min halflife
kRb = .1; %rebinding rate from local pool
kRR = .1; %rate to re-randomize, return to large pool
kP = 1;
gammap = log(2)/20;

tMax = 300;
startTime = 1;
dt = .1;
tspan = startTime:dt:tMax;
tspan2 = length(tspan);



a = zeros(NumGenes*2,1);
Genes = zeros(1,NumGenes*2);
Genes(NumGenes+1:end) = 1;
BoundRibosomeArray = zeros(1000,1);
mRNATotals = zeros(1,NumGenes);
aGenes = zeros(NumGenes*2,1);

LocalPools = zeros(1,NumGenes);
aLocalPools = zeros(NumGenes*2,1);


%Build Rxn Matrix and Propensity for genes and local pools
RxnMatrixGene = zeros(NumGenes*2,NumGenes*2);
RxnMatrixLocalPool = zeros(NumGenes*2,NumGenes*2);
for i = 1:NumGenes
    RxnMatrixGene(i,i) = 1;
    RxnMatrixGene(i,NumGenes+i) = -1;
    aGenes(i) = kON;
    RxnMatrixLocalPool(i,i) = -1;
    RxnMatrixLocalPool(i,NumGenes+i) = -1;
end
for i = 1:NumGenes
    RxnMatrixGene(i+NumGenes,i) = -1;
    RxnMatrixGene(i+NumGenes,NumGenes+i) = 1;
    aGenes(i+NumGenes) = kOFF;
end

%

RxnCount = 1;
T = 0; %Time
mRNAArray = zeros(MaxOutput,7);
Ttrack = zeros(length(tspan), 1); %Time tracking
XGenes = zeros(length(tspan), NumGenes*2); %Species number tracking
XmRNA = zeros(length(tspan),length(mRNAArray(:,1)));
XRibosomes = XmRNA;
XProtein = XRibosomes;
Ribosomes = RibosomesInitial;
XGenes(1,:)   = Genes;
xCurrentGenes = Genes;
xCurrentLocalPools = LocalPools;
xCurrentmRNA = zeros(MaxOutput,1);
xCurrentRibosomes = zeros(MaxOutput,1);
xCurrentProtein = zeros(MaxOutput,1);
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
    aGammam = xCurrentmRNA(:,1) .* gammam;
    %protein propensities
    aProtein = xCurrentRibosomes .* kP;
    aProteinDecay = xCurrentProtein .* gammap;

    %local pool propensities
    for j = 1:NumGenes
        aLocalPools(j) = xCurrentLocalPools(j).*kRb.*mRNATotals(j);
        aLocalPools(j+NumGenes) = xCurrentLocalPools(j).*kRR;
    end


    % Compute tau and mu using random variables
    a0 = sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + sum(aGammam) + ...
        sum(aLocalPools) + +sum(aProtein) + sum(aProteinDecay);
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
            XProtein(RecordCount,:) = xCurrentProtein;
            RecordCount = RecordCount + 1;
            RecordTime = RecordTime + dt;
            tau0 = tau0 - dt;
        end
    end

    if T + tau > RecordTime
        XGenes(RecordCount,:) = xCurrentGenes;
        XmRNA(RecordCount,:) = xCurrentmRNA(:,1);
        XRibosomes(RecordCount,:) = xCurrentRibosomes;
        XProtein(RecordCount,:) = xCurrentProtein;
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
    elseif r(2)*a0 < sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + sum(aGammam)
        mu  = find((sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + cumsum(aGammam) >= r(2)*a0),1,'first');
        %find the next change in state before tau
        T   = T  + tau;
        %xCurrentGenes = xCurrentGenes + RxnMatrixGene(mu,:);
        mRNATotals(mRNAArray(mu,2)) = mRNATotals(mRNAArray(mu,2)) - 1;
        %mRNAnum = mRNAnum - 1;
        xCurrentmRNA(mu) = 0;
        %xCurrentmRNA(end+1) = 0;
        mRNAArray(mu,1) = 0;
        %mRNAArray(end+1,:) = 0;
           
        UnbindingRibosomes = xCurrentRibosomes(mu);
        mRNAArray(mu,6) = UnbindingRibosomes;
        mRNAArray(mu,7) = xCurrentProtein(mu);
        xCurrentLocalPools(mRNAArray(mu,2)) = xCurrentLocalPools(mRNAArray(mu,2)) + UnbindingRibosomes;
        BurstRTrack(mRNAArray(mu,2)) = BurstRTrack(mRNAArray(mu,2))+1;
        BurstSizeRTrack(BurstRTrack(mRNAArray(mu,2)),mRNAArray(mu,2)) = UnbindingRibosomes;
        xCurrentRibosomes(mu) = 0;

    elseif r(2)*a0 < sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + sum(aGammam)...
            + sum(aLocalPools)
        mu  = find((sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + sum(aGammam) + cumsum(aLocalPools) >= r(2)*a0),1,'first');
        %find the next change in state before tau
        T   = T  + tau;
        if mu <= NumGenes
            xCurrentLocalPools(mu) = xCurrentLocalPools(mu) - 1;
            randmRNA = randi([1,mRNATotals(mu)],1);
            [idx,~] = find(mRNAArray(:,1) == 1);
            temp = 0;
            count = 0;
            for k = idx
                if mRNAArray(k,2) == mu
                    count  = count + 1;
                end
                if count == randmRNA
                    xCurrentRibosomes(k) = xCurrentRibosomes(k) + 1;
                end
            end  
        else
            xCurrentLocalPools(mu-NumGenes) = xCurrentLocalPools(mu-NumGenes) - 1;
            Ribosomes = Ribosomes + 1;
        end
    elseif r(2)*a0 < sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + sum(aGammam)...
            + sum(aLocalPools) + sum(aProtein)
       mu  = find((sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + sum(aGammam) ...
            + sum(aLocalPools) + cumsum(aProtein) >= r(2)*a0),1,'first');
        %find the next change in state before tau
        T   = T  + tau;
        xCurrentProtein(mu) = xCurrentProtein(mu) + 1;
    else
        mu  = find((sum(aGenes) + sum(aAlpha) + sum(aRibosomes) + cumsum(aGammam)...
            + sum(aLocalPools) + sum(aProtein) + cumsum(aProteinDecay)>= r(2)*a0),1,'first');
        %find the next change in state before tau
        T   = T  + tau;
        xCurrentProtein(mu) = xCurrentProtein(mu) - 1;

    end



end


% Record output
 XGenes(RecordCount,:) = xCurrentGenes;
 XmRNA(RecordCount,:) = xCurrentmRNA(:,1);
 XRibosomes(RecordCount,:) = xCurrentRibosomes;
 XProtein(RecordCount,:) = xCurrentProtein;

 mRNATotnum(i) = sum(mRNATotals);
%          for j = 1:mRNATotnum(i)
%             RibosomesPerGene(i,mRNAArray(j,2)) = RibosomesPerGene(i,mRNAArray(j,2)) + xCurrentRibosomes(j);
%          end

% Record output
count = count - 1;

name = sprintf('DataRibosomeSSRun%g',jobNum);
save(name);
ElapsedTime = toc;
%quit();
end
