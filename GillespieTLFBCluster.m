function GillespieTLFBCluster(jobNum)

%Cluster Variables

%Limit number of cpus this program can utilize
maxNumCompThreads(1);

seeds = load('seeds.mat');
seeds = seeds.seeds;

%Seed random number generator
rng(seeds(jobNum));

%base Gillespie Algorithm for general purpose simulations.

%Parameters
kONo = .1;
kOFF = 1;
kM = 20;
kP = 100;
gammam = 1;
gammap = .01;
protein0 = 8000;

runs = 1;
tMax = 1000;
dt = 1;
tspan = 0:dt:tMax;
tspan2 = length(tspan);

Name = sprintf('TFFBResults%g.mat',jobNum);

x0 = [0,0,0,1];

RxnMatrix = [ 1  0  0  0; %trascription
              0  1  0  0; %translation
             -1  0  0  0; %mRNA decay
              0 -1  0  0; %protein decay
              0  0 -1  1; %Burst OFF
              0  0  1 -1];%Burst ON

          
          
for i = 1:runs
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

    %%%%%%%%%%%%%%%%%%%%%%
    %Gillespie Simulation%
    %%%%%%%%%%%%%%%%%%%%%%

    while T <= tMax


        x = xCurrent;
        mRNA    = x(1);
        protein = x(2);
        StateON = x(3);
        StateOFF = x(4);
        
        kON = kONo/(1+(protein/protein0));

        % Calculate reaction propensities
        a = [kM*StateON; kP*mRNA; gammam*mRNA; gammap*protein;
            kOFF*StateON; kON*StateOFF];

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
    mRNAData(:,i) = X(1:tMax/dt,1);
    ProteinData(:,i) = X(1:tMax/dt,2);
end
%%
% meanProtein = mean(ProteinData,1);
% varProtein = var(ProteinData,1);
% cv2Protein = varProtein./meanProtein.^2;
% cvProtein = sqrt(varProtein)./meanProtein;
% fanoProtein = varProtein./meanProtein;

save(Name,'ProteinData');
%%
