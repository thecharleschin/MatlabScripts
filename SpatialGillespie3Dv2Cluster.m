function SpatialGillespie3Dv2Cluster(jobNum)

%Cluster Variables

%Limit number of cpus this program can utilize
maxNumCompThreads(1);

seeds = load('seeds.mat');
seeds = seeds.seeds;

%Seed random number generator
rng(seeds(jobNum));

%%
%%%%%%%%
%INPUTS%
%%%%%%%%

%Turn on commenting display? Mainly for debugging
Comments = 0;
%Turn on gene circuit?
Expression = 1;
%reflective[0] or periodic[1] boundaries
Boundary = 0;
%Wall effects on diffusion?
WallEffects = 0;

%Running variables
tmax = 1000;
dt = 1; %step to record state
tspan = 0:dt:tmax;

%Rate Constants for expression
kM = 10; %k*V
kP = 100;
gammam = 1;
gammap = .1;

%Species
Plasmids = 1; 
Polymerases = 1;
InitialmRNA = 0;
Protein = 0;

%Space
VoxLength = 16; %number of voxels along x direction
VoxWidth = 16; %number of voxles along y direction
VoxHeight = 16;
h = .05; %Side length in um

%Diffusion Constants
Dplas = 60; %Plasmid Bulk Diffusion Constant
Dpoly = 60; %Polymerase Bulk Diffusion Constant
DplasW = 0; %Plasmid Wall Diffusion Constant
DpolyW = 0; %Polymerase Wall Diffusion Constant

%Crowders
Crowding = 1; %0 no crowding, 1 crowding
DynamicCrowding = 1; %0 static crowders, 1 moving crowders
CrowdVol = 50; % percent volume crowded
CrowdD = .006; %Crowder Diffusion Constant
CrowdWallD = .006; %Crowder Diffusion Constant on Wall

Name = sprintf('%gx%gx%gCrowding%gper006DResults%g.mat',VoxLength,VoxWidth,VoxHeight,...
    CrowdVol,jobNum);
%%
tic

%%%%%%%%%%%%%%%
%Batch running%
%%%%%%%%%%%%%%%

%Inital calculations
Volume = (VoxLength*VoxWidth*VoxHeight)*h^3;
kM = kM*(Volume/h^3);
if CrowdVol == 0
    NumCrowders = 0;
else
    NumCrowders = floor(VoxLength*VoxWidth*VoxHeight*(CrowdVol/100));
end
MaxPassage = 500000; %expected maximum number of passages per particle

%kM = kM/h^3; %effective kM second order reaction

% %Initialize seperate species tracking variables
if Plasmids == 1 && Polymerases == 1
    FirstPassageTimes = 0;
    PassageDurations = zeros(MaxPassage,1);
    PassageDurationsLoc = zeros(MaxPassage,1);
    PassageTimes = zeros(MaxPassage,1);
    PassageCounters = 0;
else
%     FirstPassageTimes = zeros(Plasmids,Polymerases,Runs);
%     PassageDurations = zeros(Plasmids,Polymerases,MaxPassage,Runs);
%     PassageTimes = zeros(Plasmids,Polymerases,MaxPassage,Runs);
%     PassageCounters = zeros(Plasmids,Polymerases,Runs);
end

PassageLocations = zeros(VoxLength,VoxWidth,VoxHeight);
PassageLocationsTotal = zeros(VoxLength,VoxWidth,VoxHeight);
ProteinTrackAll = zeros(length(tspan),1);
mRNATrackAll = zeros(length(tspan),1);
mRNATimeTrack = zeros(MaxPassage,1);
mRNAperBurstTrack = zeros(MaxPassage,1);
PlasmidTrackAll = zeros(Plasmids,3,length(tspan));
PolymeraseTrackAll = zeros(Polymerases,3,length(tspan));
CrowdSpaceAll = zeros(VoxLength,VoxWidth,VoxHeight);
CrowdTrackAll = zeros(NumCrowders,3,length(tspan));
TotalEffectiveSpace = 0;
PolymeraseEffectiveSpaceAll = zeros(Polymerases,1);
PlasmidEffectiveSpaceAll = zeros(Plasmids,1);
PlasmidWallHits = zeros(Plasmids,1);
PolymeraseWallHits = zeros(Polymerases,1);
PlasmidCrowderHits = zeros(Plasmids,1);
PolymeraseCrowderHits = zeros(Polymerases,1);
PlasmidVoxDurationTotal = zeros(VoxLength,VoxWidth,VoxHeight);
%disp('initializing')

%%
%%%%%%%%%%%%%%%%%%%%%%
%Initialize Gillespie%
%%%%%%%%%%%%%%%%%%%%%%

%Tracking Arrays for storing data
%Protein is not tracked spatially, can just be placed in single array
ProteinTrack = zeros(1,length(tspan));
ProteinTrack(1) = Protein;%round(Protein + randn*sqrt(Protein)); %First point is initial protein count

mRNATrack = zeros(1,length(tspan));
mRNATrack(1) = InitialmRNA;


%SpatialRxnMatrix
SpatialRxnMatrix = [ 1  0  0; %move up 1 in x
                    -1  0  0; %move down 1 in x
                     0  1  0; %move up 1 in y
                     0 -1  0; %move down 1 in y
                     0  0  1; %move up 1 in z
                     0  0 -1]; %move down 1 in z
%NonSpatial RxnMatrix
%mRNA, Protein
%RxnMatrix = 1;%transcription
 RxnMatrix = [ 1  0; %transcription
               0  1; %translation
              -1  0; %mRNA decay
               0 -1];%Protein decay

%track spatial species on an individual basis
%rows,for each species, xyz location, for all time.
PlasmidTrack = zeros(Plasmids,3,length(tspan));
PolymeraseTrack = zeros(Polymerases,3,length(tspan));
CrowdTrack = zeros(NumCrowders,3,length(tspan));

PassageDurationTrack = zeros(Plasmids,Polymerases,MaxPassage);
PassageDurationLocTrack = zeros(Plasmids,Polymerases,MaxPassage);
PassageTimesTrack = zeros(Plasmids,Polymerases,MaxPassage);
PassageCounterMatrix = ones(Plasmids,Polymerases);
PassageLocationsMatrix = zeros(VoxLength, VoxWidth, VoxHeight);
FirstPassageMatrix = NaN(Plasmids,Polymerases);
MeetMatrix = zeros(Plasmids,Polymerases);
PlasmidWallhit = zeros(Plasmids,1);
PlasmidVoxDuration = zeros(VoxLength,VoxWidth,VoxHeight);
PlasmidVoxCurrentTime = zeros(Plasmids,1);
PolymeraseWallhit = zeros(Polymerases,1);
PlasmidCrowderhit = zeros(Plasmids,1);
PolymeraseCrowderhit = zeros(Polymerases,1);
%Species tracked on a compartment basis. 2D 1to1 map xy, z species, k over
%time
SpeciesSpace = zeros(VoxLength,VoxWidth,VoxHeight,3);
PlasmidSpace = zeros(VoxLength,VoxWidth,VoxHeight);
PolymeraseSpace = zeros(VoxLength,VoxWidth,VoxHeight);

%Track propensities across each compartment space
aSpace = zeros(4,1);

%Distribute crowders in the system
CrowdSpace = zeros(VoxLength,VoxWidth,VoxHeight);
MaxTrials = 500;

if Crowding == 1
    %place crowders randomly in space up to percent volume fill
    failures = 0;
    for i = 1:NumCrowders
        trial = 0;
        while trial < MaxTrials
            comp1 = ceil(rand*VoxLength);
            comp2 = ceil(rand*VoxWidth);
            comp3 = ceil(rand*VoxHeight);
            if CrowdSpace(comp1,comp2,comp3) == 0
                CrowdSpace(comp1,comp2,comp3) = 1;
                CrowdTrack(i,:,1) = [comp1,comp2,comp3];
                break
            end
        end
        if trial == MaxTrials
            failures = failures + 1;
            disp('Crowder placement failure')
            disp(failures)
        end
    end

    CrowdSpaceAll(:,:,:) = CrowdSpace;
    %Distribute spatial species
    failures = 0;
    for i = 1:Plasmids
        trial = 0;
        while trial < MaxTrials
            comp1 = ceil(rand*VoxLength);
            comp2 = ceil(rand*VoxWidth);
            comp3 = ceil(rand*VoxHeight);
            if CrowdSpace(comp1,comp2,comp3) == 0
                PlasmidSpace(comp1,comp2,comp3,1) = PlasmidSpace(comp1,comp2,comp3,1) + 1;
                PlasmidTrack(i,:,1) = [comp1,comp2,comp3];
                break
            end
        end
        if trial == MaxTrials
            failures = failures + 1;
            disp('Plasmid placement failure')
            disp(failures)
        end
    end

    failures = 0;
    for i = 1:Polymerases
        trial = 0;
        while trial < MaxTrials
            comp1 = ceil(rand*VoxLength);
            comp2 = ceil(rand*VoxWidth);
            comp3 = ceil(rand*VoxHeight);
            if CrowdSpace(comp1,comp2,comp3) == 0
                PolymeraseSpace(comp1,comp2,comp3) = PolymeraseSpace(comp1,comp2,comp3) + 1;
                PolymeraseTrack(i,:,1) = [comp1,comp2,comp3];
                break
            end
        end
        if trial == MaxTrials
            failures = failures + 1;
            disp('Polymerase placement failure')
            disp(failures)
        end
    end
end

%Randomly distribute Species into compartments
%Plasmids are on and off according to steady state value
if Crowding == 0
    for i = 1:Plasmids
        comp = ceil([rand*VoxLength, rand*VoxWidth, rand*VoxHeight]);
        PlasmidSpace(comp(1),comp(2),comp(3)) = PlasmidSpace(comp(1),comp(2),comp(3)) + 1;
        PlasmidTrack(i,:,1) = comp;
    end

    for i = 1:Polymerases
        comp = ceil([rand*VoxLength, rand*VoxWidth, rand*VoxHeight]);
        PolymeraseSpace(comp(1),comp(2),comp(3)) = PolymeraseSpace(comp(1),comp(2),comp(3)) + 1;
        PolymeraseTrack(i,:,1) = comp;
    end
end


%%
%Check for explorable space each particle is in.
if Crowding == 1
    TotalEffectiveSpace = VoxLength*VoxWidth*VoxHeight - sum(CrowdSpace(:));
    PlasmidEffectiveSpace = zeros(Plasmids,1);

    PolymeraseEffectiveSpace = zeros(Polymerases,1);

    %randomly explore space from initial condition to see space size
    MaxExplore = 200*TotalEffectiveSpace;
    PlasmidESpace = zeros(VoxLength,VoxWidth,VoxHeight);
    PolymeraseESpace = zeros(VoxLength,VoxWidth,VoxHeight);

    for j = 1:Plasmids
        CurrentSpace = PlasmidTrack(j,:,1);
        PlasmidESpace(CurrentSpace(1),CurrentSpace(2),CurrentSpace(3)) = 1;
        for i = 1:MaxExplore
            RandD = ceil(rand*6);
            NewSpace = CurrentSpace + SpatialRxnMatrix(RandD,:);
            if any(NewSpace == 0)
            elseif NewSpace(1) > VoxLength || NewSpace(2) > VoxWidth || NewSpace(3) > VoxHeight
            elseif CrowdSpace(NewSpace(1),NewSpace(2),NewSpace(3)) == 1 
            else
            CurrentSpace = NewSpace;
            PlasmidESpace(CurrentSpace(1),CurrentSpace(2),CurrentSpace(3)) = 1;
            end

        end
        PlasmidEffectiveSpace(j) = sum(PlasmidESpace(:));
    end
    for j = 1:Polymerases
        CurrentSpace = PolymeraseTrack(j,:,1);
        PolymeraseESpace(CurrentSpace(1),CurrentSpace(2),CurrentSpace(3)) = 1;
        for i = 1:MaxExplore
            RandD = ceil(rand*6);
            NewSpace = CurrentSpace + SpatialRxnMatrix(RandD,:);
            if any(NewSpace == 0)
            elseif NewSpace(1) > VoxLength || NewSpace(2) > VoxWidth || NewSpace(3) > VoxHeight
            elseif CrowdSpace(NewSpace(1),NewSpace(2),NewSpace(3)) == 1 
            else
            CurrentSpace = NewSpace;
            PolymeraseESpace(CurrentSpace(1),CurrentSpace(2),CurrentSpace(3)) = 1;
            end

        end
        PolymeraseEffectiveSpace(j) = sum(PolymeraseESpace(:));
    end

    PolymeraseEffectiveSpaceAll = PolymeraseEffectiveSpace(:);
    PlasmidEffectiveSpaceAll = PlasmidEffectiveSpace(:);

end

SpaceIdx = PlasmidSpace + PolymeraseSpace;
SpaceIdx(SpaceIdx > 0) = 1;

mRNACurrent = mRNATrack(1);
ProteinCurrent = ProteinTrack(1);
PlasmidCurrent = PlasmidTrack(:,:,1);
PolymeraseCurrent = PolymeraseTrack(:,:,1);
CrowdCurrent = CrowdTrack(:,:,1);

%Hop Propensities by individual species
aHopPlasmids = zeros(Plasmids,6);
aHopPolymerase = zeros(Polymerases,6);
aHopPlasmids(:,:) = Dplas/h^2;
aHopPolymerase(:,:) = Dpoly/h^2;

if WallEffects == 0

else
    for i = 1:Plasmids

        if PlasmidCurrent(i,1) == 1 
            if PlasmidCurrent(i,2) == 1 || PlasmidCurrent(i,2) == VoxWidth || ...
                PlasmidCurrent(i,3) == 1 || PlasmidCurrent(i,3) == VoxHeight
            else
            aHopPlasmids(i,1) = DplasW/h^2;
            end
        elseif PlasmidCurrent(i,1) == VoxLength 
            if PlasmidCurrent(i,2) == 1 || PlasmidCurrent(i,2) == VoxWidth || ...
                PlasmidCurrent(i,3) == 1 || PlasmidCurrent(i,3) == VoxHeight
            else
                aHopPlasmids(i,2) = DplasW/h^2;
            end
        elseif PlasmidCurrent(i,2) == 1 
            if PlasmidCurrent(i,1) == 1 || PlasmidCurrent(i,1) == VoxLength || ...
                PlasmidCurrent(i,3) == 1 || PlasmidCurrent(i,3) == VoxHeight
            else
                aHopPlasmids(i,3) = DplasW/h^2;
            end
        elseif PlasmidCurrent(i,2) == VoxWidth 
            if  PlasmidCurrent(i,1) == 1 || PlasmidCurrent(i,1) == VoxLength || ...
                PlasmidCurrent(i,3) == 1 || PlasmidCurrent(i,3) == VoxHeight
            else
                aHopPlasmids(i,4) = DplasW/h^2;
            end
        elseif PlasmidCurrent(i,3) == 1
            if PlasmidCurrent(i,1) == 1 || PlasmidCurrent(i,1) == VoxLength || ...
                PlasmidCurrent(i,2) == 1 || PlasmidCurrent(i,2) == VoxWidth
            else
                aHopPlasmids(i,5) = DplasW/h^2;
            end
        elseif PlasmidCurrent(i,3) == VoxHeight
            if PlasmidCurrent(i,1) == 1 || PlasmidCurrent(i,1) == VoxLength || ...
                PlasmidCurrent(i,3) == 1 || PlasmidCurrent(i,3) == VoxWidth
            else
                aHopPlasmids(i,6) = DplasW/h^2;
            end
        end
    end

    for i = 1:Polymerases
        if PolymeraseCurrent(i,1) == 1 
            if PolymeraseCurrent(i,2) == 1 || PolymeraseCurrent(i,2) == VoxWidth || ...
                PolymeraseCurrent(i,3) == 1 || PolymeraseCurrent(i,3) == VoxHeight
            else
            aHopPolymerase(i,1) = DplasW/h^2;
            end
        elseif PolymeraseCurrent(i,1) == VoxLength 
            if PolymeraseCurrent(i,2) == 1 || PolymeraseCurrent(i,2) == VoxWidth || ...
                PolymeraseCurrent(i,3) == 1 || PolymeraseCurrent(i,3) == VoxHeight
            else
                aHopPolymerase(i,2) = DplasW/h^2;
            end
        elseif PolymeraseCurrent(i,2) == 1 
            if PolymeraseCurrent(i,1) == 1 || PolymeraseCurrent(i,1) == VoxLength || ...
                PolymeraseCurrent(i,3) == 1 || PolymeraseCurrent(i,3) == VoxHeight
            else
                aHopPolymerase(i,3) = DplasW/h^2;
            end
        elseif PolymeraseCurrent(i,2) == VoxWidth 
            if  PolymeraseCurrent(i,1) == 1 || PolymeraseCurrent(i,1) == VoxLength || ...
                PolymeraseCurrent(i,1) == 1 || PolymeraseCurrent(i,3) == VoxHeight
            else
                aHopPolymerase(i,4) = DplasW/h^2;
            end
        elseif PolymeraseCurrent(i,3) == 1
            if PolymeraseCurrent(i,1) == 1 || PolymeraseCurrent(i,1) == VoxLength || ...
                PolymeraseCurrent(i,2) == 1 || PolymeraseCurrent(i,2) == VoxWidth
            else
                aHopPolymerase(i,5) = DplasW/h^2;
            end
        elseif PolymeraseCurrent(i,3) == VoxHeight
            if PolymeraseCurrent(i,1) == 1 || PolymeraseCurrent(i,1) == VoxLength || ...
                PolymeraseCurrent(i,3) == 1 || PolymeraseCurrent(i,3) == VoxWidth
            else
                aHopPolymerase(i,6) = DplasW/h^2;
            end
        end
    end
end


if DynamicCrowding == 0
    if NumCrowders == 0
        aHopCrowders = zeros(1,6);
    else
        aHopCrowders = zeros(NumCrowders,6);
    end
else

    aHopCrowders = zeros(NumCrowders,6);
    aHopCrowders(:,:) = CrowdD/h^2;
    if WallEffects == 1
        for i = 1:NumCrowders
            if CrowdCurrent(i,1) == 1
                aHopCrowders(i,1) = CrowdWallD/h^2;
            elseif CrowdCurrent(i,1) == VoxLength
                aHopCrowders(i,2) = CrowdWallD/h^2;
            elseif CrowdCurrent(i,2) == 1
                aHopCrowders(i,3) = CrowdWallD/h^2;
            elseif CrowdCurrent(i,2) == VoxWidth
                aHopCrowders(i,4) = CrowdWallD/h^2;
            elseif CrowdCurrent(i,3) == 1
                aHopCrowders(i,5) = CrowdWallD/h^2;
            elseif CrowdCurrent(i,3) == VoxHeight
                aHopCrowders(i,6) = CrowdWallD/h^2;
            end
        end
    end
end

aHopPlasmidTot = sum(aHopPlasmids(:));
aHopPolymeraseTot = sum(aHopPolymerase(:));
aHopTot = sum(aHopPlasmids(:)) + sum(aHopPolymerase(:));
aHopCrowdersTot = sum(aHopCrowders(:));

x = [InitialmRNA, Protein];
TTrack = zeros(length(tspan), 1); %Time tracking
T = 0; %Start Time
RxnCount = 1;  %Reaction counter
tempBurstCount = 0;
RecordTime = dt; %Recording time
RecordCount = 2;

if Comments == 1 %Used for debugging
    disp('Plasmids')
    disp(PlasmidSpace)
    disp('Polymerase')
    disp(PolymeraseSpace)
    disp('CrowdSpace')
    disp(CrowdSpace)
    pause
end

%%
%%%%%%%%%%%%%%%%%%%
%Spatial Gillespie%
%%%%%%%%%%%%%%%%%%%
Timer = 0;
%disp('Gillespie Start')
while T < tmax %Run gillespie until time is up

    %Timer for console display
    if floor(T) > Timer
        Timer = Timer + 10;
        disp(Timer)
    end

    %Check for first passage times
    %For each plasmid against each polymerase
    for i = 1:Plasmids
        for j = 1:Polymerases
            %If the two particles share the same voxel
            if PlasmidCurrent(i,:) == PolymeraseCurrent(j,:)
                %If they've never met before, record first passage time
                if isnan(FirstPassageMatrix(i,j))
                    FirstPassageMatrix(i,j) = T;
                end
                %If they were not previously in the same voxel
                %record that they are now and record passage time
                if MeetMatrix(i,j) == 0
                    MeetMatrix(i,j) = 1;
                    PassageDurationTrack(i,j,PassageCounterMatrix(i,j)) = T;
                    PassageTimesTrack(i,j,PassageCounterMatrix(i,j)) = T;
                    hitX = PlasmidCurrent(i,1); hitY = PlasmidCurrent(i,2);
                    hitZ = PlasmidCurrent(i,3);
                    PassageLocationsMatrix(hitX,hitY,hitZ) = PassageLocationsMatrix(hitX,hitY,hitZ) + 1;
                    if hitX == 1 || hitX == VoxLength || hitY == 1 || hitY == VoxWidth ...
                            || hitZ == 1 || hitZ == VoxHeight
                        PassageDurationLocTrack(i,j,PassageCounterMatrix(i,j)) = 1;                           

                    else
                        PassageDurationLocTrack(i,j,PassageCounterMatrix(i,j)) = 0;
                    end
                else

                end
            %If they the two particles are not in the same voxel but 
            %were previously, record the duration of their stay
            elseif MeetMatrix(i,j) == 1
                MeetMatrix(i,j) = 0;
                PassageDurationTrack(i,j,PassageCounterMatrix(i,j)) = T - PassageDurationTrack(i,j,PassageCounterMatrix(i,j));
                PassageCounterMatrix(i,j) = PassageCounterMatrix(i,j) + 1;
                mRNAperBurstTrack(PassageCounterMatrix(i,j)) = tempBurstCount;
                tempBurstCount = 0;

            end
        end
    end

    %Grab species spaces for each species being tracked
%         plasmid = PlasmidSpace(:,:,:);
%         polymerase = PolymeraseSpace(:,:,:);
    mRNA = x(1);%mRNACurrent;
    protein = x(2);%ProteinCurrent;

    %Recalculate all hop propensities for each compartment if needed
    %Hop Propensities by individual species
    if WallEffects == 1
        aHopPlasmids = zeros(Plasmids,6);
        aHopPolymerase = zeros(Polymerases,6);
        aHopPlasmids(:,:) = Dplas/h^2;
        aHopPolymerase(:,:) = Dpoly/h^2;


        for i = 1:Plasmids

            if PlasmidCurrent(i,1) == 1 
                if PlasmidCurrent(i,2) == 1 || PlasmidCurrent(i,2) == VoxWidth || ...
                    PlasmidCurrent(i,3) == 1 || PlasmidCurrent(i,3) == VoxHeight
                else
                aHopPlasmids(i,1) = DplasW/h^2;
                end
            elseif PlasmidCurrent(i,1) == VoxLength 
                if PlasmidCurrent(i,2) == 1 || PlasmidCurrent(i,2) == VoxWidth || ...
                    PlasmidCurrent(i,3) == 1 || PlasmidCurrent(i,3) == VoxHeight
                else
                    aHopPlasmids(i,2) = DplasW/h^2;
                end
            elseif PlasmidCurrent(i,2) == 1 
                if PlasmidCurrent(i,1) == 1 || PlasmidCurrent(i,1) == VoxLength || ...
                    PlasmidCurrent(i,3) == 1 || PlasmidCurrent(i,3) == VoxHeight
                else
                    aHopPlasmids(i,3) = DplasW/h^2;
                end
            elseif PlasmidCurrent(i,2) == VoxWidth 
                if  PlasmidCurrent(i,1) == 1 || PlasmidCurrent(i,1) == VoxLength || ...
                    PlasmidCurrent(i,3) == 1 || PlasmidCurrent(i,3) == VoxHeight
                else
                    aHopPlasmids(i,4) = DplasW/h^2;
                end
            elseif PlasmidCurrent(i,3) == 1
                if PlasmidCurrent(i,1) == 1 || PlasmidCurrent(i,1) == VoxLength || ...
                    PlasmidCurrent(i,2) == 1 || PlasmidCurrent(i,2) == VoxWidth
                else
                    aHopPlasmids(i,5) = DplasW/h^2;
                end
            elseif PlasmidCurrent(i,3) == VoxHeight
                if PlasmidCurrent(i,1) == 1 || PlasmidCurrent(i,1) == VoxLength || ...
                    PlasmidCurrent(i,3) == 1 || PlasmidCurrent(i,3) == VoxWidth
                else
                    aHopPlasmids(i,6) = DplasW/h^2;
                end
            end
        end

        for i = 1:Polymerases
            if PolymeraseCurrent(i,1) == 1 
                if PolymeraseCurrent(i,2) == 1 || PolymeraseCurrent(i,2) == VoxWidth || ...
                    PolymeraseCurrent(i,3) == 1 || PolymeraseCurrent(i,3) == VoxHeight
                else
                aHopPolymerase(i,1) = DplasW/h^2;
                end
            elseif PolymeraseCurrent(i,1) == VoxLength 
                if PolymeraseCurrent(i,2) == 1 || PolymeraseCurrent(i,2) == VoxWidth || ...
                    PolymeraseCurrent(i,3) == 1 || PolymeraseCurrent(i,3) == VoxHeight
                else
                    aHopPolymerase(i,2) = DplasW/h^2;
                end
            elseif PolymeraseCurrent(i,2) == 1 
                if PolymeraseCurrent(i,1) == 1 || PolymeraseCurrent(i,1) == VoxLength || ...
                    PolymeraseCurrent(i,3) == 1 || PolymeraseCurrent(i,3) == VoxHeight
                else
                    aHopPolymerase(i,3) = DplasW/h^2;
                end
            elseif PolymeraseCurrent(i,2) == VoxWidth 
                if  PolymeraseCurrent(i,1) == 1 || PolymeraseCurrent(i,1) == VoxLength || ...
                    PolymeraseCurrent(i,1) == 1 || PolymeraseCurrent(i,3) == VoxHeight
                else
                    aHopPolymerase(i,4) = DplasW/h^2;
                end
            elseif PolymeraseCurrent(i,3) == 1
                if PolymeraseCurrent(i,1) == 1 || PolymeraseCurrent(i,1) == VoxLength || ...
                    PolymeraseCurrent(i,2) == 1 || PolymeraseCurrent(i,2) == VoxWidth
                else
                    aHopPolymerase(i,5) = DplasW/h^2;
                end
            elseif PolymeraseCurrent(i,3) == VoxHeight
                if PolymeraseCurrent(i,1) == 1 || PolymeraseCurrent(i,1) == VoxLength || ...
                    PolymeraseCurrent(i,3) == 1 || PolymeraseCurrent(i,3) == VoxWidth
                else
                    aHopPolymerase(i,6) = DplasW/h^2;
                end
            end
        end

        aHopPlasmidTot = sum(aHopPlasmids(:));
        aHopPolymeraseTot = sum(aHopPolymerase(:));
        aHopTot = sum(aHopPlasmids(:)) + sum(aHopPolymerase(:));
        if DynamicCrowding == 1
            aHopCrowders = zeros(NumCrowders,6);
            aHopCrowders(:,:) = CrowdD/h^2;
            for i = 1:NumCrowders
                if CrowdCurrent(i,1) == 1
                    aHopCrowders(i,1) = CrowdWallD/h^2;
                elseif CrowdCurrent(i,1) == VoxLength
                    aHopCrowders(i,2) = CrowdWallD/h^2;
                elseif CrowdCurrent(i,2) == 1
                    aHopCrowders(i,3) = CrowdWallD/h^2;
                elseif CrowdCurrent(i,2) == VoxWidth
                    aHopCrowders(i,4) = CrowdWallD/h^2;
                elseif CrowdCurrent(i,3) == 1
                    aHopCrowders(i,5) = CrowdWallD/h^2;
                elseif CrowdCurrent(i,3) == VoxHeight
                    aHopCrowders(i,6) = CrowdWallD/h^2;
                end
            end

            aHopCrowdersTot = sum(aHopCrowders(:));
        end
    end

    %gaSpace = gpuArray(aSpace);
    if Expression == 0
        aSpaceTot = 0;
    else
        if Plasmids == 1 && Polymerases == 1
            if PlasmidCurrent == PolymeraseCurrent
                amRNA = kM;
            else
                 amRNA = 0;
            end
            aSpace = [amRNA; kP*mRNA; gammam*mRNA; gammap*protein];
            aSpaceTot = sum(aSpace(:));
        else     
            amRNA = kM.*PlasmidSpace.*PolymeraseSpace;
            aSpace = [sum(amRNA(:)); kP*mRNA; gammam*mRNA; gammap*protein];
            aSpaceTot = sum(aSpace(:));
        end
    end

    % Compute tau and mu using random variables
    a0 = aSpaceTot + aHopTot + aHopCrowdersTot;
    r = rand(1,2);
    %tau = -log(r(1))/a0;
    tau = (1/a0)*log(1/r(1));
    %[~, mu] = histc(r(2)*a0, [0;cumsum(a(:))]);

    %Store information if at time
    if T + tau > RecordTime
        mRNATrack(RecordCount) = x(1);%mRNACurrent;
        ProteinTrack(RecordCount) = x(2);%ProteinCurrent;
        PlasmidTrack(:,:,RecordCount) = PlasmidCurrent;
        PolymeraseTrack(:,:,RecordCount) = PolymeraseCurrent;
        CrowdTrack(:,:,RecordCount) = CrowdCurrent;

        RecordCount = RecordCount + 1;
        RecordTime = RecordTime + dt;

    end

    if Comments == 1
        disp('Random Propensity')
        disp(r(2)*a0)
    end
    if r(2)*a0 <=  aHopCrowdersTot
        if Comments == 1
            disp('CrowderHop')
            disp(aHopCrowdersTot)
            pause
        end
        if WallEffects == 0
        RandCrowder = ceil(rand*NumCrowders);
        RandD = ceil(rand*6);
        else
            [RandCrowder,RandD] = ind2sub([NumCrowders,6],find((cumsum(aHopCrowders(:)) >= r(2)*a0),1,'first'));
        end

        OldPos = CrowdCurrent(RandCrowder,:);
        NewPos = CrowdCurrent(RandCrowder,:) + SpatialRxnMatrix(RandD,:);
        if any(NewPos == 0)
            if Boundary == 0
            else
                if NewPos(1) == 0
                    NewPos(1) = VoxLength;
                elseif NewPos(2) == 0
                    NewPos(2) = VoxWidth;
                else
                    NewPos(3)= VoxHeight;
                end
                CrowdCurrent(RandCrowder,:) = NewPos;
                CrowdSpace(OldPos(1),OldPos(2),OldPos(3)) = 0;
                CrowdSpace(NewPos(1),NewPos(2),NewPos(3)) = 1;
            end
        elseif NewPos(1) > VoxLength || NewPos(2) > VoxWidth || NewPos(3) > VoxHeight
            if Boundary == 0
            else
                if NewPos(1) > VoxLength 
                    NewPos(1) = 1;
                elseif NewPos(2) > VoxWidth
                    NewPos(2) = 1;
                else
                    NewPos(3)= 1;
                end
                CrowdCurrent(RandCrowder,:) = NewPos;
                CrowdSpace(OldPos(1),OldPos(2),OldPos(3)) = 0;
                CrowdSpace(NewPos(1),NewPos(2),NewPos(3)) = 1;
            end
        elseif CrowdSpace(NewPos(1),NewPos(2),NewPos(3)) == 1 || SpaceIdx(NewPos(1),NewPos(2),NewPos(3)) == 1


        else
            CrowdCurrent(RandCrowder,:) = NewPos;
            CrowdSpace(OldPos(1),OldPos(2),OldPos(3)) = 0;
            CrowdSpace(NewPos(1),NewPos(2),NewPos(3)) = 1;

        end

    elseif r(2)*a0 <= aHopPlasmidTot + aHopCrowdersTot %If within plasmid hopping region

        if WallEffects == 0
            HoppingIdx = ceil(rand*(Plasmids));
            RandD = ceil(rand*6);
        else
            [HoppingIdx,RandD] = ind2sub([Plasmids,6],find((cumsum(aHopPlasmids(:)) + aHopCrowdersTot) >= r(2)*a0,1,'first'));
            %[HoppingIdx,RandD] = find((cumsum(aHopPlasmids) + aHopCrowdersTot >= r(2)*a0),1,'first');
        end


        if Comments == 1
        disp('Plasmid Hop')
        disp(RandD)
        pause
        end

        OldPos = PlasmidCurrent(HoppingIdx,:);
        NewPos = PlasmidCurrent(HoppingIdx,:) + SpatialRxnMatrix(RandD,:);
        if any(NewPos == 0)
            if Boundary == 0
                PlasmidWallhit(HoppingIdx,1) = PlasmidWallhit(HoppingIdx,1) + 1;  
            else
                if NewPos(1) == 0
                    NewPos(1) = VoxLength;
                elseif NewPos(2) == 0
                    NewPos(2) = VoxWidth;
                else
                    NewPos(3)= VoxHeight;
                end
                %Index Plasmid to new Position
                PlasmidCurrent(HoppingIdx,:) = NewPos;
                %Move Plasmid in space
                PlasmidSpace(NewPos(1),NewPos(2),NewPos(3)) = PlasmidSpace(NewPos(1),NewPos(2),NewPos(3)) + 1;
                PlasmidSpace(OldPos(1),OldPos(2),OldPos(3)) = PlasmidSpace(OldPos(1),OldPos(2),OldPos(3)) - 1;
                %Track duration Plasmid spent in previous voxel
                PlasmidVoxDuration(OldPos(1),OldPos(2),OldPos(3)) = PlasmidVoxDuration(OldPos(1),OldPos(2),OldPos(3)) + T - PlasmidVoxCurrentTime(HoppingIdx);
                PlasmidVoxCurrentTime(HoppingIdx) = T;
            end

        elseif NewPos(1) > VoxLength || NewPos(2) > VoxWidth || NewPos(3) > VoxHeight
            if Boundary == 0
                PlasmidWallhit(HoppingIdx,1) = PlasmidWallhit(HoppingIdx,1) + 1; 
            else
                if NewPos(1) > VoxLength 
                    NewPos(1) = 1;
                elseif NewPos(2) > VoxWidth
                    NewPos(2) = 1;
                else
                    NewPos(3)= 1;
                end
                PlasmidCurrent(HoppingIdx,:) = NewPos;
                PlasmidSpace(NewPos(1),NewPos(2),NewPos(3)) = PlasmidSpace(NewPos(1),NewPos(2),NewPos(3)) + 1;
                PlasmidSpace(OldPos(1),OldPos(2),OldPos(3)) = PlasmidSpace(OldPos(1),OldPos(2),OldPos(3)) - 1;
                 %Track duration Plasmid spent in previous voxel
                PlasmidVoxDuration(OldPos(1),OldPos(2),OldPos(3)) = PlasmidVoxDuration(OldPos(1),OldPos(2),OldPos(3)) + T - PlasmidVoxCurrentTime(HoppingIdx);
                PlasmidVoxCurrentTime(HoppingIdx) = T;
            end

        elseif CrowdSpace(NewPos(1),NewPos(2),NewPos(3)) == 1
            PlasmidCrowderhit(HoppingIdx,1) = PlasmidCrowderhit(HoppingIdx,1) + 1; 
        else
            PlasmidCurrent(HoppingIdx,:) = NewPos;
            PlasmidSpace(NewPos(1),NewPos(2),NewPos(3)) = PlasmidSpace(NewPos(1),NewPos(2),NewPos(3)) + 1;
            PlasmidSpace(OldPos(1),OldPos(2),OldPos(3)) = PlasmidSpace(OldPos(1),OldPos(2),OldPos(3)) - 1;
            %Track duration Plasmid spent in previous voxel
            PlasmidVoxDuration(OldPos(1),OldPos(2),OldPos(3)) = PlasmidVoxDuration(OldPos(1),OldPos(2),OldPos(3)) + T - PlasmidVoxCurrentTime(HoppingIdx);
            PlasmidVoxCurrentTime(HoppingIdx) = T;
        end

    elseif r(2)*a0 <= aHopPlasmidTot + aHopPolymeraseTot + aHopCrowdersTot %If within plasmid hopping region
        if WallEffects == 0
            HoppingIdx = ceil(rand*(Polymerases));
            RandD = ceil(rand*6);
        else

            [HoppingIdx,RandD] = ind2sub([Plasmids,6],find((cumsum(aHopPolymerase(:)) + aHopPlasmidTot + aHopCrowdersTot >= r(2)*a0),1,'first'));
        end
        if Comments == 1
        disp('Polymerase Hop')
        disp(RandD)
        pause
        end

        OldPos = PolymeraseCurrent(HoppingIdx,:);
        NewPos = PolymeraseCurrent(HoppingIdx,:) + SpatialRxnMatrix(RandD,:);
        if any(NewPos == 0)
            if Boundary == 0
                PolymeraseWallhit(HoppingIdx,1) = PolymeraseWallhit(HoppingIdx,1) + 1;    
            else
                if NewPos(1) == 0
                    NewPos(1) = VoxLength;
                elseif NewPos(2) == 0
                    NewPos(2) = VoxWidth;
                else
                    NewPos(3)= VoxHeight;
                end
                PolymeraseCurrent(HoppingIdx,:) = NewPos;
                PolymeraseSpace(NewPos(1),NewPos(2),NewPos(3)) = PolymeraseSpace(NewPos(1),NewPos(2),NewPos(3)) + 1;
                PolymeraseSpace(OldPos(1),OldPos(2),OldPos(3)) = PolymeraseSpace(OldPos(1),OldPos(2),OldPos(3)) - 1;
            end

        elseif NewPos(1) > VoxLength || NewPos(2) > VoxWidth || NewPos(3) > VoxHeight
            if Boundary == 0
                PolymeraseWallhit(HoppingIdx,1) = PolymeraseWallhit(HoppingIdx,1) + 1;
            else
                if NewPos(1) > VoxLength 
                    NewPos(1) = 1;
                elseif NewPos(2) > VoxWidth
                    NewPos(2) = 1;
                else
                    NewPos(3)= 1;
                end
                PolymeraseCurrent(HoppingIdx,:) = NewPos;
                PolymeraseSpace(NewPos(1),NewPos(2),NewPos(3)) = PolymeraseSpace(NewPos(1),NewPos(2),NewPos(3)) + 1;
                PolymeraseSpace(OldPos(1),OldPos(2),OldPos(3)) = PolymeraseSpace(OldPos(1),OldPos(2),OldPos(3)) - 1;
            end

        elseif CrowdSpace(NewPos(1),NewPos(2),NewPos(3)) == 1
            PolymeraseCrowderhit(HoppingIdx,1) = PolymeraseCrowderhit(HoppingIdx,1) + 1;
        else
            PolymeraseCurrent(HoppingIdx,:) = NewPos;
            PolymeraseSpace(NewPos(1),NewPos(2),NewPos(3)) = PolymeraseSpace(NewPos(1),NewPos(2),NewPos(3)) + 1;
            PolymeraseSpace(OldPos(1),OldPos(2),OldPos(3)) = PolymeraseSpace(OldPos(1),OldPos(2),OldPos(3)) - 1;
        end



    else %if not hop, do reaction

        mu  = find((cumsum(aSpace)+ aHopTot + aHopCrowdersTot >= r(2)*a0),1,'first');

        x = x + RxnMatrix(mu,:);
        
        if mu == 1
            mRNATimeTrack(RxnCount) = T + tau;
            tempBurstCount = tempBurstCount + 1;
            RxnCount = RxnCount + 1;
        end
        if Comments == 1
            disp('Do Rxn')
            disp(mu)
            pause
        end
    end

    T = T + tau;

    if Comments == 1 %Used for debugging
    disp('Plasmids')
    disp(PlasmidSpace)
    disp('Polymerase')
    disp(PolymeraseSpace)
    disp('CrowdSpace')
    disp(CrowdSpace)
    pause
    end
end

ElapsedTime = toc;
save(Name);

quit();
end