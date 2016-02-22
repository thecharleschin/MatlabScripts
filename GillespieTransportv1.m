clear all
close all
rand('state',sum(100*clock)); %#ok<RAND>

gammam = log(2)/60; %1 hour halflife
gammap = log(2)/120;%2 hour halflife
kExport = gammam;
kONArray = [0.1,0.5,1,5,10].*kExport;
kP = gammap*100;



for n = 1:length(kONArray)

    kON = kONArray(n);
    kOFFArray(1) = kON;
    kOFFArray(2) = 10*kON;
    VarArray = [1,2.5,5,7.5,10]*kON;

    for k = 1:length(VarArray)
    disp(k)

    kOFF = VarArray(k);
    alpha = 5*kOFFArray(2);




    Runs = 100;
    tMax = 5000;
    dt = 1;
    tspan = 0:dt:tMax;
    tspan2 = length(tspan);
    TraceStart = 1001;
    
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

    for i = 1:Runs

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

        Traces(i,:,:) = X(TraceStart:tspan2,:);

    end

    AllTraces(k,:,:,:) = Traces;

    end

    %%
    for k = 1:length(VarArray)
        for i = 1:Runs
            InnermRNAVar(k,i) = var(AllTraces(k,i,:,3));
            OutermRNAVar(k,i) = var(AllTraces(k,i,:,4));
            InnermRNAmean(k,i) = mean(AllTraces(k,i,:,3));
            OutermRNAmean(k,i) = mean(AllTraces(k,i,:,4));
            InnermRNAcv2(k,i) = InnermRNAVar(k,i)/ InnermRNAmean(k,i)^2;
            OutermRNAcv2(k,i) = OutermRNAVar(k,i)/OutermRNAmean(k,i)^2;       
        end
    end
    
    name = sprintf('kONData%g',n);
    save(name);

end
%%

hold on
for k = 1:length(VarArray)
    plot(VarArray(k),InnermRNAVar(k,:),'color','b','linestyle','none','marker','.',...
        'markersize',8)
    plot(VarArray(k),OutermRNAVar(k,:),'color','r','linestyle','none','marker','.',...
        'markersize',8)
end
for k = 1:length(VarArray)
    linestore(1) = plot(VarArray(k),mean(InnermRNAVar(k,:)),'markerfacecolor','b','linestyle','none','marker','o',...
        'markersize',8,'markeredgecolor','k','displayname','Inner mRNA');
    linestore(2) = plot(VarArray(k),mean(OutermRNAVar(k,:)),'markerfacecolor','r','linestyle','none','marker','o',...
        'markersize',8,'markeredgecolor','k','displayname','Outer mRNA');
end
legend(linestore,'location','southwest')
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('kOFF Value','FontSize',15)
ylabel('mRNA Variance','FontSize',15)
title('mRNA Variance over kOFF values')
saveas(gcf,'VariancevskOFF.jpg')

figure
hold on
for k = 1:length(VarArray)
    plot(VarArray(k),InnermRNAcv2(k,:),'color','b','linestyle','none','marker','.',...
        'markersize',8)
    plot(VarArray(k),OutermRNAcv2(k,:),'color','r','linestyle','none','marker','.',...
        'markersize',8)
end
for k = 1:length(VarArray)
    linestore(1) = plot(VarArray(k),mean(InnermRNAcv2(k,:)),'markerfacecolor','b','linestyle','none','marker','o',...
        'markersize',8,'markeredgecolor','k','displayname','Inner mRNA');
    linestore(2) = plot(VarArray(k),mean(OutermRNAcv2(k,:)),'markerfacecolor','r','linestyle','none','marker','o',...
        'markersize',8,'markeredgecolor','k','displayname','Outer mRNA');
end
legend(linestore,'location','southwest')
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('kOFF Value','FontSize',15)
ylabel('mRNA cv^2','FontSize',15)
title('mRNA cv^2 over kOFF values')
saveas(gcf,'cv2vskOFF.jpg')

%%
Runnum = 2;
for k = 1:length(VarArray)
    hold off
    temp1 = AllTraces(k,Runnum,:,3);
    linestore(1) = plot(temp1(:),'displayname','Inner mRNA');
    temp2 = AllTraces(k,Runnum,:,4);
    hold on
    linestore(2) = plot(temp2(:),'displayname','Outer mRNA');
    legend(linestore,'location','northeast')
    xlabel('Time (min)','FontSize',15)
    ylabel('Population (mRNA)','FontSize',15)
    name = sprintf('kOFF = %.2g',VarArray(k));
    title(name)
    name = sprintf('kOFF%gTracesRun%g.jpg',k,Runnum);
    saveas(gcf,name)
end

% BurstTimesBetweenTrack = diff(BurstTimesTrack,1);
% % GillespieData = GillespieData(200:end,:);
% len = length(GillespieData(:,1));
% GillespieDataAvg = zeros(len,1); %Average curve of all Runs
% mRNADataAvg = zeros(len,1); %Average curve of all Runs
% for h = 1:len
%     GillespieDataAvg(h) = sum(GillespieData(h,:))/Runs;
%     mRNADataAvg(h) = sum(mRNAData(h,:))/Runs;
% end
%     
% GillespieDataA = zeros(len,Runs);
% mRNADataA = zeros(len,Runs);
% 
% for k = 1:Runs
%     disp(k)
%     SSlevelsArray(k) = mean(GillespieData(:,k));
%     SSmRNAArray2(k) = mean(mRNAData(:,k));
%     
%     GillespieDataA(:,k) = GillespieData(:,k) - GillespieDataAvg;
%     mRNADataA(:,k) = mRNAData(:,k) - mRNADataAvg;
%     
%     AutoArrayTemp = xcorr(GillespieDataA(:,k),'unbiased');%,length(tspan2)-1);
%     AutomRNAArrayTemp = xcorr(mRNADataA(:,k),'unbiased');%,length(tspan2)-1);
%     AutoArray(:,k) = AutoArrayTemp;
%     VarienceArray(k) = AutoArrayTemp(len);
%     Avgcv2(k) = AutoArrayTemp(len)/(SSlevelsArray(k))^2;
%     AutomRNAArray(:,k) = AutomRNAArrayTemp;
%     VariencemRNAArray(k) = AutomRNAArrayTemp(len);
%     AvgmRNAcv2(k) = AutomRNAArrayTemp(len)/(SSmRNAArray2(k))^2;
% 
% end
%%

% 
% % 
% % AvgBurstFreq(:) = BurstFreq;
% % AvgBurstSize(:) = BurstSize;
% %%
% % c = colormap(hsv(length(kONArray)));
% % %c = colormap(hsv(length(VarArray)));
% % for i = 1:length(kONArray)
% %     hold on
% %     plot(SSlevelsTot(:,i),Avgcv2Tot(:,i),'linestyle','none','marker','.',...
% %         'markersize',10,'color',c(i,:));
% % end
% % set(gca,'XScale','log');
% % set(gca,'YScale','log');
% save 2StateEncounterData
% %%
% 
% %c = colormap(hsv(length(kONArray)));
% hold on
% for i = 1:7
%     plot(SSmRNAArray(:,i),cv2mRNAArray(:,i),'linestyle','none','marker','.',...
%         'markersize',8,'color','k');
% end
% 
% c = colormap(hsv(length(VarArray)));
% for i = 1:length(VarArray)
%     hold on
%     plot(SSmRNATot(:,i),AvgmRNAcv2Tot(:,i),'linestyle','none','marker','.',...
%         'markersize',10,'color',c(i,:));
% end
% set(gca,'XScale','log');
% set(gca,'YScale','log');
% axis([7 30 .03 1])
% xlabel('mRNA Abundace','FontSize',15)
% ylabel('cv2','FontSize',15)
% title('mRNA cv2 v Abundance From Encounters')
% saveas(gcf,'mRNAcv2vAbundance2StateDirect.jpg')

