%load data
load LizmRNAData

Trials = length(Data(1,:));
%Start all at 0
for i = 1:Trials
    Data0(:,i) = Data(:,i) - Data(1,i);
end
%General Trend
GenTrend = mean(Data0,2);

len = length(Data(:,1));

%remove general trend by some gain
for i = 1:Trials
    SSlevelsArray(i) = Data0(end,i);
    %find gain that minimizes trace
    %coarse grain
    GainRange = 0:.001:2;
    GainsTemp = zeros(1,length(GainRange));

    for h = 1:length(GainRange)
        DataTemp = Data0(:,i) - GainRange(h).*GenTrend;
        DataTemp2 = xcorr(DataTemp,GenTrend,'biased');

        GainsTemp(h) = abs(DataTemp2(len));
    end

    [~, GainMinIdx] = min(GainsTemp);

    %Finer grain 
    GainRange = GainRange(GainMinIdx) - .002:.00001: GainRange(GainMinIdx) + .002;
    GainsTemp = zeros(1,length(GainRange));

    for h = 1:length(GainRange)
        DataTemp = Data0(:,i) - GainRange(h).*GenTrend;
        DataTemp2 = xcorr(DataTemp,GenTrend,'biased');

        GainsTemp(h) = abs(DataTemp2(len));
    end

    [GainMin, GainMinIdx] = min(GainsTemp);

    %Store the gain values calculated
    GainStore(i) = GainRange(GainMinIdx);
    %Generate the noise traces
    DataA(:,i) = Data0(:,i) - GainRange(GainMinIdx).*GenTrend;
    %generate Autocorrelations
    AutoArrayTemp = xcorr(DataA(:,i),'biased');%,length(tspan2)-1);
    AutoArray(:,i) = AutoArrayTemp;
    VarianceArray(i) = AutoArrayTemp(len);
    cv2Array(i) = AutoArrayTemp(len)/(SSlevelsArray(i))^2;

    %calculate t50
    halfmax = AutoArray(len,i)/2;
    t50Idx = len;
    for j = len:length(AutoArray(:,i))
        if AutoArray(j,i) < halfmax
            t50Idx = j;
            break
        end
    end

    t50Array(i) = t50Idx-len;
end

%%
%plotting
plot(SSlevelsArray,cv2Array,'linestyle','none','marker','.','markersize',10)
set(gca,'XScale','log');
set(gca,'YScale','log');