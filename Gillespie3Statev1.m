clear all
close all
rand('state',sum(100*clock)); %#ok<RAND>

load SpatialData

% BurstSize = geomean(BurstSizemRNA,1);
% BurstFreq = geomean(BurstFreqmRNA,1);

BurstSize = BurstSizeDirectmean;
BurstFreq = BurstFreqDirectmean;



kM = 40960;
kP = 0;
gammam = 1;
gammap = .1;
lambda = 0;

% kON = 65;
% kOFF = 195130;


temp = BurstSize;%*(10^-6);
kON = 1./(1./BurstFreq-temp);
kOFF = 1./(temp);


%comparable numbers to cube 10% to 50% crowding
varmod = [1.5,1.7,1.7,1.9,2.4,3.2,5];
VarArray = kON.*varmod';%[1,6,25]; %burst durations 
k32o23Array = kON./VarArray;
k23o32Ratio = 1-k32o23Array;
k23factor = kON;
k23Array = k32o23Array.*k23o32Ratio;
k23Array = k23Array.*k23factor;
k32Array = k32o23Array.*(1-k23o32Ratio).*k23factor;


%k23Array(1:3) = 0;
% %timebetween = [0.10594824,0.0866851,0.06852042]; %empty bursts removed
 %timebetween = 1./BurstFreq; %empty bursts included

%comparable numbers for 50% crowding 2High to cube
% VarArray = [8.7,7.7,7.4];
% timebetween = [1/57.9,1/63.3,1/66.2];


%%
Runs = 100;
tMax = 200;
dt = 1/60;
tspan = 0:dt:tMax;
tspan2 = length(tspan);

InitialProtein = 0;

Plasmids = 1;%[1,6,25];
BurstDurationTrack = zeros(100000,Runs);
BurstTimesTrack = zeros(100000,Runs);
SSmRNAArray2 = zeros(Runs,1);

for j = 1:length(VarArray)


    k21 = VarArray(j);
    k12 = kOFF(j);
    k23 = k23Array(j);
    k32 = k32Array(j);


for i = 1:Runs
    disp(j)
    disp(i)
    PercO = kON/(kON + kOFF);
    
    ProteinStart = InitialProtein;%round(InitialProtein + randn*sqrt(InitialProtein));
    
    if Plasmids == 1
        if rand < kON/(kON + kOFF)
            x0 = [0, ProteinStart, 1, 0, 0];
        else
            x0 = [0, ProteinStart, 0, 1, 0];
        end
    else
        Plasmids1 = round(kON/(kON + kOFF) * Plasmids);
        x0 = [0, ProteinStart, Plasmids1, Plasmids - Plasmids1, Plasmids - Plasmids1 - Plasmids2];
    end
    
    
    RxnMatrix = [ 1  0  0  0  0; %trascription
                  0  1  0  0  0; %translation
                 -1  0  0  0  0; %mRNA decay
                  0 -1  0  0  0; %protein decay
                  0  0 -1  1  0; %1 to 2
                  0  0  1 -1  0; %2 to 1
                  0  0  0 -1  1; %2 to 3
                  0  0  0  1 -1];%3 to 2

     p = 1;
        Onflag = x0(3);
        mRNAperBurst = zeros(100000,1);

        MaxOutput = 100000; %Maximum size expected of the output file
        NumSpecies = size(RxnMatrix, 2);
        %T = zeros(MaxOutput, 1); %Time tracking
        %X = zeros(1, NumSpecies); %Species number tracking
        %T(1)     = 0;
        %X(1,:)   = x0;
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
        OnDuration = zeros(MaxOutput,1);
        OnTimes = zeros(MaxOutput,1);
        %%%%%%%%%%%%%%%%%%%%%%
        %Gillespie Simulation%
        %%%%%%%%%%%%%%%%%%%%%%

        while T <= tMax

            % Calculate reaction propensities
            x = xCurrent;
            mRNA    = x(1);
            protein = x(2);
            State1 = x(3);
            State2 = x(4);
            State3 = x(5);
            
            %track burst dynamics
            if State1 == 1
                if OnTrack == 0
                    OnTrack = 1;
                    OnTimes(count) = T;
                    OnDuration(count) = T;
                else
                end
            elseif OnTrack == 1
                OnTrack = 0;
                OnDuration(count) = T - OnDuration(count);
                count = count + 1;
                
            end
                    
            
            a = [kM*State1; kP*mRNA; gammam*mRNA; gammap*protein;
                k12*State1; k21*State2; k23*State2; k32*State3];

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
        mRNAData(:,i) = X(6001:tspan2,1);
        GillespieData(:,i) = X(6001:tspan2,2);
%         BurstDurationTrack(:,i) = OnDuration;
%         BurstTimesTrack(:,i) = OnTimes;
%         BurstNumberTrack(i) = count;
    
end

%BurstTimesBetweenTrack = diff(BurstTimesTrack,1);
% GillespieData = GillespieData(200:end,:);
len = length(GillespieData(:,1));
GillespieDataAvg = zeros(len,1); %Average curve of all Runs
mRNADataAvg = zeros(len,1); %Average curve of all Runs
for h = 1:len
    GillespieDataAvg(h) = sum(GillespieData(h,:))/Runs;
    mRNADataAvg(h) = sum(mRNAData(h,:))/Runs;
end
    
GillespieDataA = zeros(len,Runs);
mRNADataA = zeros(len,Runs);

for k = 1:Runs
    disp(k)
    SSlevelsArray(k) = mean(GillespieData(:,k));
    SSmRNAArray2(k) = mean(mRNAData(:,k));
    
    GillespieDataA(:,k) = GillespieData(:,k) - GillespieDataAvg;
    mRNADataA(:,k) = mRNAData(:,k) - mRNADataAvg;
    
    AutoArrayTemp = xcorr(GillespieDataA(:,k),'unbiased');%,length(tspan2)-1);
    AutomRNAArrayTemp = xcorr(mRNADataA(:,k),'unbiased');%,length(tspan2)-1);
    AutoArray(:,k) = AutoArrayTemp;
    VarienceArray(k) = AutoArrayTemp(len);
    Avgcv2(k) = AutoArrayTemp(len)/(SSlevelsArray(k))^2;
    AutomRNAArray(:,k) = AutomRNAArrayTemp;
    VariencemRNAArray(k) = AutomRNAArrayTemp(len);
    t50mRNAArray(k) = find(AutomRNAArrayTemp(len:end) < .5*AutomRNAArrayTemp(len),1,'first');
    AvgmRNAcv2(k) = AutomRNAArrayTemp(len)/(SSmRNAArray2(k))^2;

end
%%

AutoVar(:,j) = VarienceArray;
AvgSSVar = var(SSlevelsArray);
MeanSS = mean(SSlevelsArray);
Avgcv2Tot(:,j) = Avgcv2;
SSlevelsTot(:,j) = (SSlevelsArray);
t50mRNATot(:,j) = t50mRNAArray;

AutomRNAVar(:,j) = VariencemRNAArray;
AvgmRNAcv2Tot(:,j) = AvgmRNAcv2;
SSmRNATot(:,j) = (SSmRNAArray2);

% BurstDurationTrackTot(:,:,j) = BurstDurationTrack;
% BurstTimesTrackTot(:,:,j) = BurstTimesTrack;
% BurstNumberTrackTot(:,j) = BurstNumberTrack;

end
% 
% AvgBurstFreq(:) = BurstFreq;
% AvgBurstSize(:) = BurstSize;
%%
%c = colormap(hsv(length(kONArray)));
% c = colormap(hsv(length(VarArray)));
% for i = 1:length(VarArray)
%     hold on
%     plot(SSlevelsTot(:,i),Avgcv2Tot(:,i),'linestyle','none','marker','.',...
%         'markersize',10,'color',c(i,:));
% end
% set(gca,'XScale','log');
% set(gca,'YScale','log');
save 3StateData
%%
hold on
% for i = 1:7
%     plot(SSmRNAArray(:,i),cv2mRNAArray(:,i),'linestyle','none','marker','.',...
%         'markersize',8,'color','k');
% end
c = colormap(hsv(length(VarArray)));
%c = colormap(hsv(length(VarArray)));
CrowdArray = [0,10,20,30,40,45,50];
for i = 1:length(VarArray)
    hold on
    name = sprintf('%g%% Crowding',CrowdArray(i));
    linestore(i) = plot(SSmRNATot(:,i),AvgmRNAcv2Tot(:,i),'linestyle','none','marker','o',...
        'markersize',6,'markerfacecolor',c(i,:),'markeredgecolor','k',...
        'displayname',name);
end

legend(linestore,'location','northwest')
set(gca,'XScale','log');
set(gca,'YScale','log');
axis([7 30 .03 1])
xlabel('mRNA Abundance','FontSize',15)
ylabel('cv^2','FontSize',15)
set(gca,'fontsize',15)
title('cv^2 v Abundance 3 State')
saveas(gcf,'mRNAcv2vAbundance3State3.jpg')
saveas(gcf,'mRNAcv2vAbundance3State3.svg')

