function GillespieTransportCluster(jobNum,varNum)
%Cluster Variables

%Limit number of cpus this program can utilize
maxNumCompThreads(1);

seeds = load('seeds.mat');
seeds = seeds.seeds;

%Seed random number generator
rng(seeds(jobNum));

kONArray = [0.001155245,0.001155245,0.001155245,0.001155245,0.001155245,0.005776227,...
    0.005776227,0.005776227,0.005776227,0.005776227,0.011552453,0.011552453,...
    0.011552453,0.011552453,0.011552453,0.057762265,0.057762265,0.057762265,...
    0.057762265,0.057762265,0.11552453,0.11552453,0.11552453,0.11552453,0.11552453];

kOFFArray = [0.000115525,0.000577623,0.001155245,0.005776227,0.011552453,...
    0.000577623,0.002888113,0.005776227,0.028881133,0.057762265,0.001155245,...
    0.005776227,0.011552453,0.057762265,0.11552453,0.005776227,0.028881133,...
    0.057762265,0.288811325,0.57762265,0.011552453,0.057762265,0.11552453,...
    0.57762265,1.155245301];

gammam = log(2)/60; %1 hour halflife
kExport = gammam;
kON = kONArray(varNum);
kOFF = kOFFArray(varNum);
kP = gammam*100;
gammap = log(2)/360;%6 hour halflife
alpha = 5*kOFF;

tMax = 2000;
dt = 1;
tspan = 0:dt:tMax;
tspan2 = length(tspan);

InitialState = [0 1 0 0 0];
%Build Rxn Matrix
%[ON OFF mRNAin mRNAout];
RxnMatrix = [ 1 -1  0  0  0; %burst on
             -1  1  0  0  0; %burst off
              0  0  1  0  0; %transcription
              0  0 -1  1  0; %transport
              0  0  0 -1  0; %mRNA decay
              0  0  0  0  1; %translation
              0  0  0  0 -1];%protein decay



MaxOutput = 100000; %Maximum size expected of the output file
%T = zeros(MaxOutput, 1); %Time tracking
%X = zeros(1, NumSpecies); %Species number tracking
%T(1)     = 0;
%X(1,:)   = x0;
RxnCount = 1;
T = 0; %Time
Ttrack = zeros(length(tspan), 1); %Time tracking
X = zeros(length(tspan),5); %Species number tracking
X(1,:)   = InitialState;
xCurrent = InitialState;
RecordTime = dt; %Recording time
RecordCount = 2;
OnTrack = 0;
count = 1;
OnDuration = zeros(MaxOutput,1);
OnTimes = zeros(MaxOutput,1);
%%%%%%%%%%%%%%%%%%%%%%
%Gillespie Simulation%
%%%%%%%%%%%%%%%%%%%%%%
tic
while T <= tMax

    % Calculate reaction propensities
    x = xCurrent;
    StateON = x(1);
    StateOFF = x(2);
    mRNAin = x(3);
    mRNAout = x(4);
    protein = x(5);

    %generate a
    a = [StateOFF*kON; StateON*kOFF; StateON*alpha; mRNAin*kExport;...
        mRNAout*gammam ;mRNAout*kP; protein*gammap];

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

end


% Record output
 X(RecordCount,:) = xCurrent;

% Record output
count = count - 1;

Traces(:,:) = X(1:tspan2,:);


name = sprintf('DatakONkOFF%gRun%g',varNum,jobNum);
save(name);
ElapsedTime = toc;
quit();
end
