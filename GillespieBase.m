%base Gillespie Algorithm for general purpose simulations.

clear all
close all
rand('state',sum(100*clock)); %#ok<RAND>

%Parameters
kON = 1;
kOFF = 1;
kM = 1;
kP = 1;
gammam = .1;
gammap = .1;

tMax = 200;
dt = 1;
tspan = 0:dt:tMax;
tspan2 = length(tspan);

x0 = [0,0,0,1];

RxnMatrix = [ 1  0  0  0; %trascription
              0  1  0  0; %translation
             -1  0  0  0; %mRNA decay
              0 -1  0  0; %protein decay
              0  0 -1  1; %Burst OFF
              0  0  1 -1];%Burst ON

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

    % Calculate reaction propensities
    a = [kM*StateON; kP*mRNA; gammam*mRNA; gammap*protein;
        kOFF*StateON; kON*StateOFF];

    % Compute tau and mu using random variables
    a0 = sum(a);
    r = rand(1,2);
    %tau = -log(r(1))/a0;
    tau = (1/a0)*log(1/r(1));

    %Store information if at time
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
mRNAData = X(1:tMax/dt,1);
ProteinData = X(1:tMax/dt,2);
