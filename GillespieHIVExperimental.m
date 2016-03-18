%base Gillespie Algorithm for general purpose simulations.

clear all
close all
rand('state',sum(100*clock)); %#ok<RAND>

Runs = 50;

paramArray = logspace(0,2,5);%[.1,.5,1,5,10];
paramArray = paramArray./60;

for k = 1:length(paramArray)
    disp('Case')
    disp(k)
    %Parameters
    kON1 = .015; %burst on
    kON2 = .005;
    kOFF = 1/20; %burst off
    kM = 120.98/60;%42.4/60;%.87/60;%120.98/60; %transcription
    kS = paramArray(k); %splicing 
    Sfrac = .046; %fraction of splices that make Rev
    kP = 270/60; %translation
    gammams = .173/60; %mRNAdecay (of spliced)
    gammam = 0;%.173/60; %mRNAdecay (unspliced)
    gammap = .173/60; %protein decay from Wen Thesis

    kR1f = 1.91; %one rev binding to mRNA for 1 bound
    kR1r = log(2)/(4.95); %one rev unbinding to mRNA for 0 bound
    kR2f = 1.01; %one rev binding to mRNA for 2 bound
    kR2r = log(2)/(3.15); %one rev unbinding to mRNA for 1 bound
    kR3f = 1.73; %one rev binding to mRNA for 3 bound
    kR3r = log(2)/(3.64); %one rev unbinding to mRNA for 2 bound
    kR4f = 1.59; %one rev binding to mRNA for 4 bound
    kR4r = log(2)/(3.30); %one rev unbinding to mRNA for 3 bound

    % kR1f = 1.91*60; %one rev binding to mRNA for 1 bound
    % kR1r = log(2)/(4.95/60); %one rev unbinding to mRNA for 0 bound
    % kR2f = 1.01*60; %one rev binding to mRNA for 2 bound
    % kR2r = log(2)/(3.15/60); %one rev unbinding to mRNA for 1 bound
    % kR3f = 1.73*60; %one rev binding to mRNA for 3 bound
    % kR3r = log(2)/(3.64/60); %one rev unbinding to mRNA for 2 bound
    % kR4f = 1.59*60; %one rev binding to mRNA for 4 bound
    % kR4r = log(2)/(3.30/60); %one rev unbinding to mRNA for 3 bound

    kexport = 5.53/60; %export of mRNA and bound rev molecules value from Wen Thesis
    offTime = 5000;

    tMax = 10000;
    dt = 1;
    tspan = 0:dt:tMax;
    tspan2 = length(tspan);

    %GeneON GeneOFF mRNA mRNAs Protein mRNA1R mRNA2R mRNA3R mRNA4R

    x0 = [1,0,0,0,0,0,0,0,0,0];


    RxnMatrix = [ 1 -1  0  0  0  0  0  0  0  0; %Burst ON
                 -1  1  0  0  0  0  0  0  0  0; %Burst OFF
                  0  0  1  0  0  0  0  0  0  0; %transcription
                  0  0 -1  1  0  0  0  0  0  0; %splice
                  0  0  0  0  1  0  0  0  0  0; %translation
                  0  0  0 -1  0  0  0  0  0  0;%mRNAs decay
                  0  0  0  0 -1  0  0  0  0  0;%protein decay
                  0  0 -1  0 -1  1  0  0  0  0;%Rev1 forward
                  0  0  0  0  1 -1  0  0  0  0;%Rev1 back
                  0  0  0  1  1 -1  0  0  0  0;%Rev1 splice
                  0  0  0  0 -1 -1  1  0  0  0;%Rev2 forward
                  0  0  0  0  1  1 -1  0  0  0;%Rev2 back
                  0  0  0  1  2  0 -1  0  0  0;%Rev2 splice
                  0  0  0  0 -1  0 -1  1  0  0;%Rev3 forward
                  0  0  0  0  1  0  1 -1  0  0;%Rev3 back
                  0  0  0  1  3  0  0 -1  0  0;%Rev3 splice
                  0  0  0  0 -1  0  0 -1  1  0;%Rev4 forward
                  0  0  0  0  1  0  0  1 -1  0;%Rev4 back
                  0  0  0  1  4  0  0  0 -1  0;%Rev4 splice
                  0  0  0  0  0  0  0  0 -1  1;%export
                  0  0 -1  0  0  0  0  0  0  0;%mRNAu decay
                  0  0  0  0  0 -1  0  0  0  0;%Rev1 decay
                  0  0  0  0  0  0 -1  0  0  0;%Rev2 decay
                  0  0  0  0  0  0  0 -1  0  0;%Rev3 decay
                  0  0  0  0  0  0  0  0 -1  0];%Rev4 decay

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
        displayTime = 0;
        kON = kON1;

        %%%%%%%%%%%%%%%%%%%%%%
        %Gillespie Simulation%
        %%%%%%%%%%%%%%%%%%%%%%

        while T <= tMax
%             if T > displayTime
% %                 disp(T)
% %                 displayTime = displayTime + 10;
%             end
            if T > offTime 
                kON = kON2;

            end

            x = xCurrent;
            mRNA    = x(3);
            rev = x(5);
            StateON = x(1);
            StateOFF = x(2);
            mRNAs = x(4);
            rev1 = x(6);
            rev2 = x(7);
            rev3 = x(8);
            rev4 = x(9);

            % Calculate reaction propensities
            a = [kON*StateOFF;%Burst ON
                kOFF*StateON;           %Burst OFF
                StateON*kM; %transcription
                      mRNA*kS; %splice
                      mRNAs*kP; %translation
                     mRNAs*gammams;%mRNA decay
                      rev*gammap;%protein decay
                      mRNA*rev*kR1f;%Rev1 forward
                      rev1*kR1r;%Rev1 back
                      rev1*kS;%Rev1 splice
                      rev1*rev*kR2f;%Rev2 forward
                      rev2*kR2r;%Rev2 back
                      rev2*kS;%Rev2 splice
                      rev2*rev*kR3f;%Rev3 forward
                      rev3*kR3r;%Rev3 back
                      rev3*kS;%Rev3 splice
                      rev3*rev*kR4f;%Rev4 forward
                      rev4*kR4r;%Rev4 back
                      rev4*kS;%Rev4 splice
                      rev4*kexport;%export
                      mRNA*gammam;%mRNAu decay
                      rev1*gammam;%rev1 decay
                      rev2*gammam;%rev2 decay
                      rev3*gammam;%rev3 decay
                      rev4*gammam];%rev4 decay

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
            if mu == 4
                if rand() < Sfrac
                     xCurrent = xCurrent + RxnMatrix(mu,:);
                else
                    xCurrent = xCurrent + [0  0 -1  0  0  0  0  0  0  0];
                end
            else
            xCurrent = xCurrent + RxnMatrix(mu,:);
            end

        end


        % Record output
         X(RecordCount,:) = xCurrent;
         AllTraces(:,:,i) = X(1:tMax/dt,:);
    end
    meanRev(:,k) = mean(AllTraces(:,5,:),3);
    temp =  AllTraces(:,5,:);
    temp = reshape(temp,[tMax/dt,Runs]);
    AllRev(:,:,k) = temp;
end
save AllHIVData
%%
hold on
c = colormap(jet(length(paramArray)));
for i = 1:length(paramArray)
    %plot(AllTraces(:,5,i))
    name = sprintf('%.2g kS',paramArray(i));
    linestore(i) = plot(meanRev(:,i),'color',c(i,:),'linewidth',2,...
        'displayname',name);
end
legend(linestore,'location','northeast')
xlabel('Time(min)','FontSize',15)
ylabel('Rev Abundance','FontSize',15)
name = sprintf('Rev Traces Relaxation');
title(name)
name = sprintf('RevTraces.jpg');
saveas(gcf,name)
%%
figure
hold on
for i = 1:length(paramArray)
    normMeanRev(:,i) = meanRev(:,i)./meanRev(offTime,i);
    %plot(AllTraces(:,5,i))
    name = sprintf('%.2g kS',paramArray(i));
    linestore(i)=plot(normMeanRev(offTime:end,i),'color',c(i,:),'linewidth',2,...
        'displayname',name);
end
legend(linestore,'location','northeast')
axis([0 2000 0.3 1.1])
xlabel('Time after off (min)','FontSize',15)
ylabel('Normalized Rev Abundance','FontSize',15)
name = sprintf('Normalized Rev Traces');
title(name)
name = sprintf('RevTracesNorm.jpg');
saveas(gcf,name)
%%
set(gca,'YScale','log');
axis([0 1500 10^-2 10])
name = sprintf('RevTracesNormSemilog.jpg');
saveas(gcf,name)
%plotting
% subplot(3,1,1);
% plot(AllTraces(:,5),'g','linewidth',2)
% subplot(3,1,2);
% hold on
% c = colormap(jet(5));
% plot(AllTraces(:,3),'color',c(1,:))
% plot(AllTraces(:,6),'color',c(2,:))
% plot(AllTraces(:,7),'color',c(3,:))
% plot(AllTraces(:,8),'color',c(4,:))
% plot(AllTraces(:,9),'color',c(5,:))
% subplot(3,1,3);
% plot(AllTraces(:,4),'r','linewidth',2)
% saveas(gcf,'traces.jpg')