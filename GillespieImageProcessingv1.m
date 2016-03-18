%Processing from Picture to Abundance Traces
%Outputs - Abundance traces for all days, buckets, and ROIs.

%%%%%%
%NOTE%
%Additional Toolboxes may be required to run all functions. Known required
%toolboxes include:
%Image Processing Toolbox


%Reset workspace state
clear all
close all
rand('state',sum(100*clock)); %#ok<RAND>


%%
%%%%%%%%
%INPUTS%
%%%%%%%%
%Pictures are in same folder as script, named as 'size'um'date'.tif
%example - 10um20140312.tif

%Centers are in same folder as script, named as centers'size'um'date.mat
%example - centers10um20140423.mat

%Dates of files being considered of the form YYMMDD. 20 is appended later
dates = 140416;% [140317,140319,140326,140402,140416,140411,140423];

%size of bucket being analyzed, in um.
bucketsize = 2;

%Range of ROIs to analyze. Radius of ROI, such that diameter = Radius*2 + 1
%example, ROI 10 has radius 10 and diameter 21. start and end range inclusive.
ROIRange = 5;%1:1:70;


%%%%%%%%%%%
%MAIN LOOP%
%%%%%%%%%%%
%For each day being analyzed
for ii = 1:length(dates)
    
    %Grab the date being used
    date = sprintf('20%g',dates(ii)); %04022014, 04162014, 04232014
    disp(date)
    
    %Certain bucket sizes of a certain day only have 22 time points
    if strcmp(date,'20140312') == 1 && bucketsize ~= 5
        tMax = 22;
    else
        tMax = 23;
    end
    

    %[0] grab centers from file, [1] grab centers using algorithm
    GrabCenters = 0;
    
    if GrabCenters == 0
        DynamicFileName = sprintf('centers%gum%s.mat',bucketsize,date);
        load(DynamicFileName)
    end
    count = 1;
    xCenters = 0;
    yCenters = 0;
    BucketIdx = 0;
    BucketGroup = 0; %for use in checking 2um bucket grouping
    %For each bucket, check if center file is 0. if 0, skip, if not, store.
    for i = 1:length(bucketmax(:,1))
        if bucketmax(i,1) == 0
        else
        yCenters(count) = bucketmax(i,1);%Centers are flipped due to indexing
        xCenters(count) = bucketmax(i,2);
        BucketIdx(count) = i;
        BucketGroup(count) = bucketmax(i,3);
        count = count + 1;
        end
    end
    buckets = length(xCenters);
    BucketStart = 0;%last unusable bucket
    BucketEnd = buckets;
    
    %Grab image file
    DynamicFileName = sprintf('%gum%s.tif',bucketsize,date);
    disp(DynamicFileName)
    FileTif= DynamicFileName;
    InfoImage=imfinfo(FileTif);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    FinalImage=zeros(nImage,mImage,NumberImages,'uint16');

    TifLink = Tiff(FileTif, 'r');
    for i=1:NumberImages
       TifLink.setDirectory(i);
       FinalImage(:,:,i)=TifLink.read();
    end
    TifLink.close();
    
    %Convert to double
    Image =double(FinalImage(:,:,1:NumberImages));
    
    %For each ROI value in ROI range
    for n = 1:length(ROIRange)
        ROI = ROIRange(n); %Radius
        disp('ROI')
        disp(ROI)
        ROIArea(n) = ROI^2*pi();
        %%
        %Grab ROI data from image
        pixelTotal = zeros(tMax,buckets);
        pixelsum = zeros(tMax,buckets);
        %for each bucket
        for i = 1:buckets  
            a = xCenters(i);
            b = yCenters(i);
            r = ROI;
            radius = r+.5;
                for slice = 1:tMax
                    B = Image(:,:,slice);
                    for x= a-r:a+r
                        for y= b-r:b+r
                            %if one of these pixel values is within the bounds of this
                             %radius, then it's inside the ROI, and proceed if the value
                               %of the xy box is also within the bounds of the image
                            if sqrt((x-a)^2 + (y-b)^2)<=radius && x>=1 && y>=1 && x<=(nImage-1) && y<=(mImage-1)
                                pixelTotal(slice,i)=pixelTotal(slice,i) +1;
                                pixelsum(slice,i) = pixelsum(slice,i) + B(x,y); % pixel sums of each ROI radius                
                            end

                        end

                    end
                end

        end

       
        GillespieData = pixelsum;
        
        clear FinalImage
        
        DynamicFileName = sprintf('%gum%sROI%gTraces.mat',bucketsize,date,ROI);
        save(DynamicFileName,'GillespieData','BucketIdx','BucketGroup','bucketsize','ROIRange');

    end
end

