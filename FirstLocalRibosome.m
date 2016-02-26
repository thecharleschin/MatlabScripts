clear all
close all
rand('state',sum(100*clock)); %ok<RAND>
MaxOutput = 1000; %Maximum size expected of the output file


NumGenes = 5;
RibosomesInitial = 300*NumGenes;
Runs = 1;

kON = .1;%gene burst on rate
kOFF = 1; %gene burst off rate
alpha = 1; %mRNA production rate
gammam = log(2)/5; %mRNA decay rate, 5min halflife
kIN = 1; %rebinding rate from local pool
kOUT = .1; %rate to re-randomize, return to large pool
kP = .1;
gammap = log(2)/20;

tMax = 100;
startTime = 1;
dt = 1;
tspan = startTime:dt:tMax;
tspan2 = length(tspan);




Genes = zeros(1,NumGenes*2);
Genes(NumGenes+1:end) = 1;



GeneRuns = zeros(Runs,tspan2,NumGenes);


%Build Rxn Matrix and Propensity for genes and local pools
%#GenesON #GenesOFF #mRNAPools #LocalRibosomePools GlobalRibosomePool Protein
%NumGenes NumGenes NumGenes NumGenes 1 1
%#AllBurstON #AllBurstOFF #AllTranscription #AllmRNADecay #AllLocalIN #AllLocalOUT #AllTranslation #AllProteinDecay
%NumGenes NumGenes NumGenes NumGenes NumGenes NumGenes NumGenes NumGenes

RxnMatrix = zeros(NumGenes*8,NumGenes*4+2);

for i = 1:NumGenes
    RxnMatrix(i,i) = 1; %BurstON ONstate
    RxnMatrix(i,i+NumGenes) = -1; %BurstON Offstate
    RxnMatrix(i+NumGenes,i) = -1; %BurstOFF ONstate
    RxnMatrix(i+NumGenes,i+NumGenes) = 1; %BurstOFF OffState
    RxnMatrix(i+NumGenes*2,i+NumGenes*2) = 1; %Transcription
    RxnMatrix(i+NumGenes*3,i+NumGenes*2) = -1; %mRNADecay
    RxnMatrix(i+NumGenes*4,i+NumGenes*3) = 1; %LocalIN
    RxnMatrix(i+NumGenes*5,i+NumGenes*3) = -1; %LocalOUT
   
end

RxnMatrix(NumGenes*4+1:NumGenes*5,NumGenes*4+1) = -1; %LocalIN
RxnMatrix(NumGenes*5+1:NumGenes*6,NumGenes*4+1) = 1; %LocalOUT
RxnMatrix(NumGenes*6+1:NumGenes*7,NumGenes*4+2) = 1; %Translation
RxnMatrix(NumGenes*7+1:NumGenes*8,NumGenes*4+2) = -1; %ProteinDecay 


x0 = zeros(1,NumGenes*4+2);
x0(1+NumGenes:NumGenes*2) = 1;
x0(NumGenes*4+1) = RibosomesInitial;
aAll = zeros(NumGenes*8,1);

%

for i = 1:Runs

    disp('run')
    disp(i)
    
    RxnCount = 1;
    T = 0; %Time
    
    X = zeros(length(tspan), NumGenes*4+2); %Species number tracking
    X(1,:)   = x0;
    xCurrent = x0;

    
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
        x = xCurrent;
        StateON = x(1:NumGenes);
        StateOFF = x(NumGenes+1:NumGenes*2);
        mRNA = x(NumGenes*2+1:NumGenes*3);
        Local = x(NumGenes*3+1:NumGenes*4);
        Global = x(NumGenes*4+1);
        Protein = x(NumGenes*4+2);
        
        for j = 1:NumGenes
            aAll(j) = kON*StateOFF(j); %Burst ON
            aAll(j+NumGenes) = kOFF*StateON(j); %Burst OFF
            aAll(j+NumGenes*2) = alpha*StateON(j); %Transcription
            aAll(j+NumGenes*3) = gammam*mRNA(j); %mRNADecay
            aAll(j+NumGenes*4) = Global*kIN*mRNA(j);% LocalIN
            if mRNA(j) == 0
                aAll(j+NumGenes*5) = Local(j)*kOUT; %LocalOUT
            else
                aAll(j+NumGenes*5) = Local(j)*kOUT/mRNA(j);
            end
            aAll(j+NumGenes*6) = kP*mRNA(j)*Local(j); %Translation
            aAll(j+NumGenes*7) = gammap*Protein;    %ProteinDecay
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

    end


    % Record output
     X(RecordCount,:) = xCurrent;
   

end