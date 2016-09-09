function GillespieTLFBTFCluster(jobNum)

%Cluster Variables

%Limit number of cpus this program can utilize
maxNumCompThreads(1);

seeds = load('seeds.mat');
seeds = seeds.seeds;

%Seed random number generator
rng(seeds(jobNum));

%Parameters
kON = .1;
kOFF = 1;
kB = 0.001;
kUb = 8;
kM = 20;
kP = 100;
gammam = 1;
gammap = .01;

runs = 1;
tMax = 1000;
dt = 1;
tspan = 0:dt:tMax;
tspan2 = length(tspan);

%mRNA protein TF geneON geneBound geneFree
x0 = [0,0,1,0,0,1];
Name = sprintf('TLFBTFResults%g.mat',jobNum);

RxnMatrix = [ 1  0  0  0  0  0; %trascription
              0  1  0  0  0  0; %translation
             -1  0  0  0  0  0; %mRNA decay
              0 -1  0  0  0  0; %protein decay
              0  0 -1  1  0 -1; %TF bind to ON
			  0  0  1 -1  0  1; %TF unbind to free
              0 -1  0  0  1 -1; %P bind to OFF
			  0  1  0  0 -1  1];%P Unbind to free

          
          
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

		%mRNA protein TF geneON geneBound geneFree
        x = xCurrent;
        mRNA    = x(1);
        protein = x(2);
		TF = x(3);
        GeneON = x(4);
        ComplexOFF = x(5);
		GeneFree = x(6);
        
        

        % Calculate reaction propensities
        a = [kM*GeneON; kP*mRNA; gammam*mRNA; gammap*protein;
            kON*GeneFree*TF; kOFF*GeneON; kB*protein*GeneFree; kUb*ComplexOFF];

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
    kONTrack(i) = kON;
end
%%
meanProtein = mean(ProteinData,1);
varProtein = var(ProteinData,1);
cv2Protein = varProtein./meanProtein.^2;
cvProtein = sqrt(varProtein)./meanProtein;
fanoProtein = varProtein./meanProtein;


save(Name,'ProteinData');
quit();
end
