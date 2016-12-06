%the interp version of fileConverter saves signals that are interpolated
%with respect to longitudinal position
%run the original file, 'magneto_fileConverter', in order for time-based
%signals

clc
clear
close all

%% Get user's dropbox location
%last update ZS 06.23.16

%get user's dropbox location and prompt to get group file from there
disp('If you have added your Dropbox root already, navigate to the group file and click Open.');
disp('If you have not added your Dropbox root, navigate to that directory and click Open.');
disp('There may be an error displayed when choosing your Dropbox - ignore this.');
DropboxPath = dropbox();

%directory for saving grouped file at end of script
saveDir = [DropboxPath,filesep,'Magneto',filesep,'MagnetoData',filesep];

%directory where individual participant files live
dataDir = uigetdir([filesep,'Users',filesep,'Zach',...
    filesep,'Desktop',filesep,'Work',filesep,'Research',filesep,'DDL',filesep,...
    'Magneto',filesep,'Data',filesep,'_magnetofullset']);


%% Convert .xlsx files to .mat files and save to directory - should need to run one time
%last update ZS 06.14.16
tic
convert = 1; %change to 1 to convert xlxs files for the first time
if convert
    cd(dataDir);
    ExcelFiles = dir('*.xlsx');
    FileBaseNames = {length(ExcelFiles)};
    disp('Consolidating participant drive files & converting xlsx to mat...');
    disp('Sit back, this may take a while...');
    %pNames = [21328 21329 21335 213336 21337 21338 21339 21340 21341 21343];
    for j = 1:(size(ExcelFiles,1)/2)
        [tempName,remain] = strtok(strtok(ExcelFiles(j*2).name,'.'),'Magneto_Sub_');
        FileBaseNames{j,1} = [tempName remain(1:end-1)];
    end
    for i = 1:2:length(ExcelFiles)
        crntFileName = [FileBaseNames{(i+1)/2,1},'.mat'];
        %read in both drive files and concatenate before saving
        [dataMain1,~,~] = xlsread(ExcelFiles(i).name,'','','basic');
        eventCount1{(i+1)/2,1} = dataMain1(:,8);
        [dataMain2,~,~] = xlsread(ExcelFiles(i+1).name,'','','basic');
        eventCount2{(i+1)/2,1} = dataMain2(:,8);
        dataMain = [dataMain1;dataMain2];
        save([FileBaseNames{(i+1)/2,1},'.mat'],'dataMain');
        disp(['Converted files for participant ',crntFileName(1:5)]);
    end
    clear FileBaseNames tempName remain ExcelFiles
end

%% Grab subject names from Magneto .mat files
%last update ZS 05.23.16

cd(dataDir);
MatFileNames = dir('*.mat');
for i = 1:length(MatFileNames)
    fileName = MatFileNames(i).name;
    subjectNames{i} = fileName(1:5);
    fileName_complete{i} = fileName;
end

fileName_complete = fileName_complete';
subjectNames = str2num(cell2mat(subjectNames'));
nFiles = length(subjectNames);
baseNames = unique(subjectNames);
nSubjects = length(baseNames);

%insert relevant data into fields
for i = 1:length(MatFileNames)
    Subjects(i,1).Names = subjectNames(i);
    Subjects(i,1).fileName = fileName_complete{i};
end

%% Allocate memory for a cell array for each signal
%last update ZS 03.15.16

%%%%%%%%%ALLOCATE DATA CELLS & SIGNAL NAMES%%%%%%%%%%%%
signalNames = {'nTrials','nSubjects','nFiles','eventCount','time','DirectionTrigger',...
    'StartTrigger','EndTrigger','ResetTrigger','StartCount','StartEventMat',...
    'DirEventMat','ResetSim','SteeringControl','UserID','SignDirection',...
    'SteeringState','Heading','HeadingError','LaneNum','LaneOffset',...
    'HWPosition','Vx','Vy','PosX','PosY','PosZ','Pitch','Yaw',...
    'subjectNames','cones','Subjects'};

%% Read in .mat files and parse data into vectors
%last update ZS 4.27.16

InterpMethod = 'spline';

%store (time based) measurements from imported spreadsheet in vectors
for i = 1:length(subjectNames)
    disp(['Processing ',fileName_complete{i}]);
    load(fileName_complete{i});
    %eventCount = dataMain(:,8);
    eventCount{i,1} = [eventCount1{i};(eventCount2{i} + max(eventCount1{i}))];
    nTrials_current(i,1) = max(eventCount{i});
    eventCountVec = 1:nTrials_current(i,1);
    
    for j = 1:nTrials_current(i,1)
        TrialIdxRange = find(eventCount{i} == eventCountVec(j));
        %find indices of unique posy values for later signal interpolation
        PosX_whole{i,j} = dataMain(TrialIdxRange,32); %longitudinal position
        StartTrigger_whole{i,j} = dataMain(TrialIdxRange,2);
        
        [posx_unique,PosXUniqueIdx,~] = unique(PosX_whole{i,j},'stable');
        idx_test{i,j} = PosXUniqueIdx;
        
        StartTrigger{i,j} = dataMain(TrialIdxRange,2);
        StartTrigger{i,j} = interp1(posx_unique,StartTrigger{i,j}(PosXUniqueIdx),...
            posx_unique,InterpMethod);
        
        startTrialIdx = find(StartTrigger{i,j} ~= StartTrigger{i,j}(1),1);
        if isempty(startTrialIdx)
            disp(['No start trial idx: Subject ',num2str(i),', Trial ',num2str(j)]);
            if j < 32
                PosX{i,j} = zeros(500,1);
            end
        else
            PosXUniqueIdx = PosXUniqueIdx(startTrialIdx:end);
            
            PosX{i,j} = PosX_whole{i,j}(PosXUniqueIdx);
            StartTrigger{i,j} = StartTrigger_whole{i,j}(PosXUniqueIdx);
            
            PosZ{i,j} = dataMain(TrialIdxRange,33); %vertical
            PosZ{i,j} = ...%PosZ{i,j}(PosYUniqueIdx);
                interp1(PosX{i,j},PosZ{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            
            PosY{i,j} = dataMain(TrialIdxRange,31); %lateral
            PosY{i,j} = ...%= PosX{i,j}(PosYUniqueIdx);
                (interp1(PosX{i,j},PosY{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod));
            PosYOffset(i,j) = PosY{i,j}(1);
            PosY{i,j} = PosY{i,j} - PosYOffset(i,j); %subtract offset (first value of vector)
            time{i,j} = dataMain(TrialIdxRange,1);
            %time{i,j} = time{i,j}(PosYUniqueIdx);
            time{i,j} = (interp1(PosX{i,j},time{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod));
            DirectionTrigger{i,j} = dataMain(TrialIdxRange,3);
            DirectionTrigger{i,j} = interp1(PosX{i,j},DirectionTrigger{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            EndTrigger{i,j} = dataMain(TrialIdxRange,11);
            EndTrigger{i,j} = interp1(PosX{i,j},EndTrigger{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            ResetTrigger{i,j} = dataMain(TrialIdxRange,4);
            ResetTrigger{i,j} = interp1(PosX{i,j},ResetTrigger{i,j}(PosXUniqueIdx),...x
                PosX{i,j},InterpMethod);
            StartCount{i,j} = dataMain(TrialIdxRange,5);
            StartEventMat{i,j} = dataMain(TrialIdxRange,6);
            StartEventMat{i,j} = interp1(PosX{i,j},StartEventMat{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            DirEventMat{i,j} = dataMain(TrialIdxRange,3);
            DirEventMat{i,j} = interp1(PosX{i,j},DirEventMat{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            ResetSim{i,j} = dataMain(TrialIdxRange,7);
            ResetSim{i,j} = interp1(PosX{i,j},ResetSim{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            SteeringControl{i,j} = dataMain(TrialIdxRange,9);
            SteeringControl{i,j} = interp1(PosX{i,j},SteeringControl{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            UserID{i,j} = dataMain(TrialIdxRange,13);
            SignDirection{i,j} = dataMain(TrialIdxRange,15);
            SignDirection{i,j} = interp1(PosX{i,j},SignDirection{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            Heading{i,j} = dataMain(TrialIdxRange,18);
            Heading{i,j} = interp1(PosX{i,j},Heading{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            HeadingError{i,j} = dataMain(TrialIdxRange,19);
            HeadingError{i,j} = interp1(PosX{i,j},HeadingError{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            LaneNum{i,j} = dataMain(TrialIdxRange,22);
            LaneNum{i,j} = interp1(PosX{i,j},LaneNum{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            LaneOffset{i,j} = dataMain(TrialIdxRange,23);
            LaneOffset{i,j} = interp1(PosX{i,j},LaneOffset{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            
            %Convert handwheel position from radians to degrees
            HWPosition{i,j} = rad2deg(dataMain(TrialIdxRange,25));
            HWPosition{i,j} = interp1(PosX{i,j},HWPosition{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            HWPosition{i,j} = HWPosition{i,j} - HWPosition{i,j}(1); %subtract offset (first value of vector)
            
            SteeringState{i,j} = dataMain(TrialIdxRange,16);
            SteeringState{i,j} = interp1(PosX{i,j},SteeringState{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            
            Vx{i,j} = dataMain(TrialIdxRange,28); %longitudinal velocity
            Vx{i,j} = interp1(PosX{i,j},Vx{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            Vy{i,j} = dataMain(TrialIdxRange,29); %lateral velocity
            Vy{i,j} = interp1(PosX{i,j},Vy{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            
            %roll = dataMain(:,33);
            Pitch{i,j} = dataMain(TrialIdxRange,35);
            Pitch{i,j} = interp1(PosX{i,j},Pitch{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            Yaw{i,j} = dataMain(TrialIdxRange,36);
            Yaw{i,j} = interp1(PosX{i,j},Yaw{i,j}(PosXUniqueIdx),...
                PosX{i,j},InterpMethod);
            
            %             SideslipFR{i,j} = dataMain(TrialIdxRange,37);
            %             SideslipFR{i,j} = interp1(PosX{i,j},SideslipFR{i,j}(PosXUniqueIdx),...
            %                 PosX{i,j},InterpMethod);
            %             SideslipFL{i,j} = dataMain(TrialIdxRange,38);
            %             SideslipFL{i,j} = interp1(PosX{i,j},SideslipFL{i,j}(PosXUniqueIdx),...
            %                 PosX{i,j},InterpMethod);
        end
    end
end

%find the largest value of the eventCount vectors that were just processed
nTrials = max(nTrials_current);

%% Input cone locations
%last update ZS 04.01.16

%values taken from ISA markers placed *manually* :(
%IMPORTANT: the cone position values are ordered to take advantage of spline function
%later, do not edit order of cone vectors

cones.X = [
    -175.7,...
    -183.2,...
    -190.7,...
    -220.7,...
    -228.2,...
    -235.7,...
    -265.7,...
    -273.2,...
    -280.7,...
    ] * -1; %needed to invert (check ISA to SAE conversion)

cones.YRightLower = [
    4.910,...
    4.910,...
    4.910,...
    1.550,...
    1.550,...
    1.550,...
    4.910,...
    4.910,...
    4.910,...
    ] + PosYOffset(1,1); %apply offset to match lateral position offset

cones.YRightUpper = [
    8.947,...
    8.947,...
    8.947,...
    5.611,...
    5.611,...
    5.611,...
    8.947,...
    8.947,...
    8.947,...
    ] + PosYOffset(1,1); %apply offset to match lateral position offset

cones.YLeftLower = [
    4.910,...
    4.910,...
    4.910,...
    8.299,...
    8.299,...
    8.299,...
    4.910,...
    4.910,...
    4.910,...
    ] + PosYOffset(1,1); %apply offset to match lateral position offset

cones.YLeftUpper = [
    8.947,...
    8.947,...
    8.947,...
    12.35,...
    12.35,...
    12.35,...
    8.947,...
    8.947,...
    8.947,...
    ] + PosYOffset(1,1); %apply offset to match lateral position offset

cones.LaneWidth = cones.YLeftUpper(1) - cones.YLeftLower(1);


%% Interpolate cone position to create a track
%last update ZS 04.04.16

%Find spline through each set of cones
degree = 7;
x1 = linspace(min(cones.X),max(cones.X),25);

cones.p1 = polyfit(cones.X,cones.YRightLower,degree);
cones.RightLowerPoints = polyval(cones.p1,x1);

cones.p2 = polyfit(cones.X,cones.YRightUpper,degree);
cones.RightUpperPoints = polyval(cones.p2,x1);

cones.p3 = polyfit(cones.X,cones.YLeftLower,degree);
cones.LeftLowerPoints = polyval(cones.p3,x1);

cones.p4 = polyfit(cones.X,cones.YLeftUpper,degree);
cones.LeftUpperPoints = polyval(cones.p4,x1);

%% Save grouped data
%last update ZS 06.23.16
saveName = ['MagnetoGroup_',date];
cd(saveDir);
disp('Saving group data...');
save([saveName,'.mat'],signalNames{:});
disp(['Successfully saved ', saveName,' to ',saveDir]);
toc
