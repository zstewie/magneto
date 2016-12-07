clc
clear
close all

%% Load group file from user's dropbox
%last update ZS 05.17.16

%get user's dropbox location and prompt to get group file from there
disp('If you have added your Dropbox root already, navigate to the group file and click Open.');
disp('If you have not added your Dropbox root, navigate to that directory and click Open.');
disp('There may be an error displayed when choosing your Dropbox - ignore this and run again.');
disp('Locate group .mat file.');
DropboxPath = dropbox();

[GroupFileName,GroupPathName] = uigetfile([DropboxPath,filesep,'Magneto',filesep,...
    'MagnetoData'],'Select group file');
disp('Loading group file...');
load ([GroupPathName GroupFileName]);
disp(['Successfully loaded ',GroupFileName]);

saveDir = [DropboxPath,filesep,'Magneto',filesep,'MagnetoData',filesep];

%% Create settings structure to control processing output of analysis
%last update ZS 06.21.16

settings.includeIncorrectTurns = input('Include incorrect trials [0/1]?: ');
settings.choosePathDev = input('Absolute[0] or relative[1] lateral error?: ');
settings.saveGroupMetrics = input('Save group metrics at end? [0/1]?: ');

%set defaults is no value is entered by user
if isempty(settings.includeIncorrectTurns)
    settings.includeIncorrectTurns = 1; %include incorrect turns by default
end
if isempty(settings.choosePathDev)
    settings.choosePathDev = 1; %relative by default
end
if isempty(settings.saveGroupMetrics)
    settings.saveGroupMetrics = 1; %save group metrics by default
end

%% Preallocate matrices/cells
%last update ZS 2.21.16

%length of snipped signals - all metric signals are first indexed between
%direction and end triggers, then snipped to 400 indices so that all are
%the same length and can take mean, std, etc...
snipLength = 400;

%Preallocate congruence matrix
congruence = NaN(size(PosX));

%Preallocate turn direction matrix
TurnDirection = NaN(size(PosX));

%Trigger event indices
StartEvent.Idx = zeros(size(PosX));
DirEvent.Idx = zeros(size(PosX));
EndEvent.Idx = zeros(size(PosX));
ResetEvent.Idx = zeros(size(PosX));

StartEvent.Pos = zeros(size(PosX));
StartEvent.Time = zeros(size(PosX));

DirEvent.Pos = NaN(size(PosX));
%DirEvent.Idx = NaN(size(PosX));
DirEvent.Time = zeros(size(PosX));

EndEvent.Pos = NaN(size(PosX));
EndEvent.Time = zeros(size(PosX));

ResetEvent.Pos = zeros(size(PosX));
ResetEvent.Time = zeros(size(PosX));

%Handwheel first peak index
HWTurnPkIdx = zeros(size(PosX));
HWTurnPkPos = zeros(size(PosX));
HWTurnPkTime = zeros(size(PosX));
HWTurnPkVal = zeros(size(PosX));
TurnPeakSign = zeros(size(PosX));

%Handwheel metrics
HWPks = cell(size(PosX));
HWLocs = cell(size(PosX));
ReversalsPerSec = NaN(size(PosX));
deltaT = NaN(size(PosX));

%acceleration metrics
accelX = cell(size(PosX));
accelY = cell(size(PosX));

%grouping by trial number (1-16 first block, 17-32 second block)
groupNum = NaN(size(PosX));

%filtering matrices - flags good/bad & correct/incorrect trials
EligibleTrials = ones(size(PosX)); % assume all eligible trials to start (ones) - eligible meaning that all triggers exist in the trial
CorrectTurn = ones(size(PosX)); % assuem all correct to begin - correct meaning participant completed the double lane change maneuver correctly

%matrix for repeat flag
repeatFlag = zeros(size(PosX));

%matrix for output
participant = cell(size(PosX));

nSubjects = size(Subjects,1);

%% Manual eligibility check (based on redcap notes)
%last update ZS 07.22.16

%subject 18554
% EligibleTrials(1,[7,20,23,36,29]) = 0;
% 
% %subject 19119
% EligibleTrials(2,[2,9,13,20,25,32]) = 0;
% 
% %subject 20418 all good
% %subject 20502 all good
% %subject 20588 all good
% %subject 20660 all good
% %subject 20662 all good
% 
% %subject 20679
% EligibleTrials(7,5) = 0;
% 
% %subject 20685
% EligibleTrials(2,13) = 0;
% 
% %subject 20735
% EligibleTrials(10,[10,17,20,28]) = 0;
% 
% %subject 20743
% EligibleTrials(11,[2,27]) = 0;
% 
% %subject 20749 all good
% %subject 20751 all good
% 
% %subject 20752
% EligibleTrials(14,[5,9,14,19,23,28]) = 0;
% 
% %subject 20758
% EligibleTrials(15,[13,16,20]) = 0;
% 
% %subject 20760 - unclear
% %subject 20765 all good
% %subject 20769 all good
% 
% %subject 20770
% EligibleTrials(19,[1,2,5,14,17,20,28]) = 0;
% 
% %subject 20790
% EligibleTrials(20,5) = 0;
% 
% %subject 20791
% EligibleTrials(21,[15,31]) = 0;

%% Identify event triggers, populate congruence/turn direction matrix
%last update ZS 03.18.16

nNominalTrials = 32;
TotalTrialCount = 0;
dt = time{1,1}(2) - time{1,1}(1);

%for each subject, check to see how many trials were completed
for i = 1:nSubjects
    if isempty(find(cellfun(@isempty,PosX(i,:)),1))
        nTrials_current = nTrials - 1;
    else
        nTrials_current = find(cellfun(@isempty,PosX(i,:)),1) - 1;
    end
    Subjects(i).nTrials = nTrials_current;
    EligibleTrials(i,Subjects(i).nTrials+1:end) = 0;
    Subjects(i).nRepeatTrials = nTrials_current - nNominalTrials;
    for j = 1:Subjects(i).nTrials
        if (isempty(StartTrigger{i,j}) || isempty(DirEventMat{i,j}) || isempty(EndTrigger{i,j}) || isempty(ResetTrigger{i,j}))
            %disp(['No trigger position found - subject ',num2str(i),', trial ',num2str(j)]);
            EligibleTrials(i,j) = 0;
            continue
        else
            StartEventIdx_temp = find(StartTrigger{i,j} ~= StartTrigger{i,j}(1),1);
            %DirEventIdx_temp = find(DirEventMat{i,j} ~= DirEventMat{i,j}(1),1);
            %DirEventIdx_temp = find(DirectionTrigger{i,j} ~= DirectionTrigger{i,j}(1),1);
            DirEventIdx_temp = findTrigger(PosX{i,j},'dir');
            %EndEventIdx_temp = find(EndTrigger{i,j} ~= EndTrigger{i,j}(1),1);
            EndEventIdx_temp = findTrigger(PosX{i,j},'end');
            ResetEventIdx_temp = find(ResetTrigger{i,j} ~= ResetTrigger{i,j}(1),1);
        end
        
        if(isempty(DirEventIdx_temp) || isempty(EndEventIdx_temp) || isempty(ResetEventIdx_temp) || isempty(StartEventIdx_temp))
            %disp(['No trigger position found - subject ',num2str(i),', trial ',num2str(j)]);
            EligibleTrials(i,j) = 0;
            StartEventIdx_temp = 0;
            DirEventIdx_temp = 0;
            EndEventIdx_temp = 0;
            ResetEventIdx_temp = 0;
        else
            StartEvent.Idx(i,j) = StartEventIdx_temp;
            DirEvent.Idx(i,j) = DirEventIdx_temp;
            EndEvent.Idx(i,j) = EndEventIdx_temp;
            ResetEvent.Idx(i,j) = ResetEventIdx_temp;
            
            %populate repeatFlag matrix
            if(j > (Subjects(i).nTrials - Subjects(i).nRepeatTrials))
                repeatFlag(i,j) = 1;
            end
            
            %check is this cell is empty, place NaN if so
            if isempty(SteeringState{i,j})
                congruence(i,j) = NaN;
            else
                %populate congruence matrix based on status of
                %SteeringState at the index of the direction trigger
                if SteeringState{i,j}(DirEvent.Idx(i,j)) == 1;
                    congruence(i,j) = 0; %SteeringState = 1 --> incongruent trial
                elseif SteeringState{i,j}(DirEvent.Idx(i,j)) == 0;
                    congruence(i,j) = 1; %SteeringState = 0 --> congruent trial
                end
            end
            
            %use secondary check to confirm steering congruence
            congruence(i,j) = verifySteerState(HWPosition{i,j}(DirEvent.Idx(i,j)),...
                PosY{i,j}(DirEvent.Idx(i,j)));
            
            %populate turn direction matrix based on status of
            %SignDirection at the index of the direction trigger
            if SignDirection{i,j}(DirEvent.Idx(i,j)) == 1
                TurnDirection(i,j) = 1; %TurnDirection = 1 --> right turn signal
            elseif SignDirection{i,j}(DirEvent.Idx(i,j)) == 0
                TurnDirection(i,j) = 0; %TurnDirection = 0 --> left turn signal
            end
            
            StartEvent.Pos(i,j) = PosX{i,j}(StartEvent.Idx(i,j));
            StartEvent.Time(i,j) = time{i,j}(StartEvent.Idx(i,j));
            
            DirEvent.Pos(i,j) = PosX{i,j}(DirEvent.Idx(i,j));
            DirEvent.Time(i,j) = time{i,j}(DirEvent.Idx(i,j));
            
            EndEvent.Pos(i,j) = PosX{i,j}(EndEvent.Idx(i,j));
            EndEvent.Time(i,j) = time{i,j}(EndEvent.Idx(i,j));
            
            ResetEvent.Pos(i,j) = PosX{i,j}(ResetEvent.Idx(i,j));
            ResetEvent.Time(i,j) = time{i,j}(ResetEvent.Idx(i,j));
            
            VyStruct.mean(i,j) = mean(Vy{i,j},'omitnan');
            VyStruct.rms(i,j) = rms(Vy{i,j},'omitnan');
            
            %Calculate acceleration from velocity in x and y
            accelX{i,j} = diff(Vx{i,j}) ./ diff(time{i,j}); %longitudinal
            accelY{i,j} = diff(Vy{i,j}) ./ diff(time{i,j}); %lateral
            
            accelYStruct.mean(i,j) = mean(accelY{i,j},'omitnan');
            accelYStruct.rms(i,j) = rms(accelY{i,j},'omitnan');
            
            %filter hw position, rate, accel, jerk
            wn = 100; % Hz
            [b,a] = butter(3,wn*2/500);
            HWPosition{i,j} = filtfilt(b,a,HWPosition{i,j});
            
            wn = 50; % Hz
            [b,a] = butter(3,wn*2/500);
            
            %calculate 1st, 2nd, 3rd time derivatives of handwheel position
            HWRateStruct.matrix{i,j} = diff(HWPosition{i,j}) ./ diff(time{i,j});
            HWRateStruct.matrix{i,j} = filtfilt(b,a,HWRateStruct.matrix{i,j});
            [HWRateStruct.max(i,j),HWRateStruct.maxIdx(i,j)] = max(HWRateStruct.matrix{i,j});
            HWRateStruct.mean(i,j) = mean(HWRateStruct.matrix{i,j});
            HWRateStruct.rms(i,j) = rms(HWRateStruct.matrix{i,j});
            
            HWAccelStruct.matrix{i,j} = diff(HWRateStruct.matrix{i,j}) ./ diff(time{i,j}(1:end-1));
            HWAccelStruct.matrix{i,j} = filtfilt(b,a,HWAccelStruct.matrix{i,j});
            [HWAccelStruct.max(i,j),HWAccelStruct.maxIdx(i,j)] = max(HWAccelStruct.matrix{i,j});
            HWAccelStruct.mean(i,j) = mean(HWAccelStruct.matrix{i,j});
            HWAccelStruct.rms(i,j) = rms(HWAccelStruct.matrix{i,j});
            
            HWJerkStruct.matrix{i,j} = diff(HWAccelStruct.matrix{i,j}) ./ diff(time{i,j}(1:end-2));
            HWJerkStruct.matrix{i,j} = filtfilt(b,a,HWJerkStruct.matrix{i,j});
            [HWJerkStruct.max(i,j),HWJerkStruct.maxIdx(i,j)] = max(HWJerkStruct.matrix{i,j});
            HWJerkStruct.mean(i,j) = mean(HWJerkStruct.matrix{i,j});
            HWJerkStruct.rms(i,j) = rms(HWJerkStruct.matrix{i,j});
            
            %calculate 1st, 2nd, 3rd time derivatives of yaw
            YawRateStruct.matrix{i,j} = diff(Yaw{i,j}) ./ diff(time{i,j});
            YawRateStruct.matrix{i,j} = filtfilt(b,a,YawRateStruct.matrix{i,j});
            [YawRateStruct.max(i,j),YawRateStruct.maxIdx(i,j)] = max(YawRateStruct.matrix{i,j});
            YawRateStruct.mean(i,j) = mean(YawRateStruct.matrix{i,j});
            YawRateStruct.rms(i,j) = rms(YawRateStruct.matrix{i,j});
            
            YawAccelStruct.matrix{i,j} = diff(YawRateStruct.matrix{i,j}) ./ diff(time{i,j}(1:end-1));
            YawAccelStruct.matrix{i,j} = filtfilt(b,a,YawAccelStruct.matrix{i,j});
            [YawAccelStruct.max(i,j),YawAccelStruct.maxIdx(i,j)] = max(YawAccelStruct.matrix{i,j});
            YawAccelStruct.mean(i,j) = mean(YawAccelStruct.matrix{i,j});
            YawAccelStruct.rms(i,j) = rms(YawAccelStruct.matrix{i,j});
            
            YawJerkStruct.matrix{i,j} = diff(YawAccelStruct.matrix{i,j}) ./ diff(time{i,j}(1:end-2));
            YawJerkStruct.matrix{i,j} = filtfilt(b,a,YawJerkStruct.matrix{i,j});
            [YawJerkStruct.max(i,j),YawJerkStruct.maxIdx(i,j)] = max(YawJerkStruct.matrix{i,j});
            YawJerkStruct.mean(i,j) = mean(YawJerkStruct.matrix{i,j});
            YawJerkStruct.rms(i,j) = rms(YawJerkStruct.matrix{i,j});
            
            %count trials
            TotalTrialCount = TotalTrialCount + 1;
        end
    end
    
    %sort metrics into subject-wise mean
%     crntHWRateRow = HWRateStruct.mean(i,:);
%     crntHWAccelRow = HWAccelStruct.mean(i,:);
%     crntHWJerkRow = HWJerkStruct.mean(i,:);
%     crntVyRow = VyStruct.mean(i,:);
%     crntAccelYRow = accelYStruct.mean(i,:);
%     crntYawRateRow = YawRateStruct.mean(i,:);
%     
%     crntHWRateCongruentMean = crntHWRateRow(congruence(i,:) == 1);
%     crntHWRateIncongruentMean = crntHWRateRow(congruence(i,:) == 0);
%     HWRateStruct.CongruentMean(i,1) = mean(crntHWRateCongruentMean,'omitnan');
%     HWRateStruct.IncongruentMean(i,1) = mean(crntHWRateIncongruentMean,'omitnan');
%     
%     crntHWAccelCongruentMean = crntHWAccelRow(congruence(i,:) == 1);
%     crntHWAccelIncongruentMean = crntHWAccelRow(congruence(i,:) == 0);
%     HWAccelStruct.CongruentMean(i,1) = mean(crntHWAccelCongruentMean,'omitnan');
%     HWAccelStruct.IncongruentMean(i,1) = mean(crntHWAccelIncongruentMean,'omitnan');
%     
%     crntHWJerkCongruentMean = crntHWJerkRow(congruence(i,:) == 1);
%     crntHWJerkIncongruentMean = crntHWJerkRow(congruence(i,:) == 0);
%     HWJerkStruct.CongruentMean(i,1) = mean(crntHWJerkCongruentMean,'omitnan');
%     HWJerkStruct.IncongruentMean(i,1) = mean(crntHWJerkIncongruentMean,'omitnan');
%     
%     crntVyCongruentMean = crntVyRow(congruence(i,:) == 1);
%     crntVyIncongruentMean = crntVyRow(congruence(i,:) == 0);
%     VyStruct.CongruentMean(i,1) = mean(crntVyCongruentMean,'omitnan');
%     VyStruct.IncongruentMean(i,1) = mean(crntVyIncongruentMean,'omitnan');
%     
%     crntAccelYCongruentMean = crntAccelYRow(congruence(i,:) == 1);
%     crntAccelYIncongruentMean = crntAccelYRow(congruence(i,:) == 0);
%     accelYStruct.CongruentMean(i,1) = mean(crntAccelYCongruentMean,'omitnan');
%     accelYStruct.IncongruentMean(i,1) = mean(crntAccelYIncongruentMean,'omitnan');
%     
%     crntYawRateCongruentMean = crntYawRateRow(congruence(i,:) == 1);
%     crntYawRateIncongruentMean = crntYawRateRow(congruence(i,:) == 0);
%     YawRateStruct.CongruentMean(i,1) = mean(crntYawRateCongruentMean,'omitnan');
%     YawRateStruct.IncongruentMean(i,1) = mean(crntYawRateIncongruentMean,'omitnan');
end

[HWRateStruct.IncongruentMean,HWRateStruct.CongruentMean] = ...
    ttestFormat(HWRateStruct.mean,congruence,'congruence');
[HWAccelStruct.IncongruentMean,HWAccelStruct.CongruentMean] = ...
    ttestFormat(HWAccelStruct.mean,congruence,'congruence');
[HWJerkStruct.IncongruentMean,HWJerkStruct.CongruentMean] = ...
    ttestFormat(HWJerkStruct.mean,congruence,'congruence');
[VyStruct.IncongruentMean,VyStruct.CongruentMean] = ...
    ttestFormat(VyStruct.mean,congruence,'congruence');
[accelYStruct.IncongruentMean,accelYStruct.CongruentMean] = ...
    ttestFormat(accelYStruct.mean,congruence,'congruence');
[YawRateStruct.IncongruentMean,YawRateStruct.CongruentMean] = ...
    ttestFormat(YawRateStruct.mean,congruence,'congruence');

[HWRateStruct.IncongruentRMS,HWRateStruct.CongruentRMS] = ...
    ttestFormat(HWRateStruct.rms,congruence,'congruence');
[HWAccelStruct.IncongruentRMS,HWAccelStruct.CongruentRMS] = ...
    ttestFormat(HWAccelStruct.rms,congruence,'congruence');
[HWJerkStruct.IncongruentRMS,HWJerkStruct.CongruentRMS] = ...
    ttestFormat(HWJerkStruct.rms,congruence,'congruence');
[VyStruct.IncongruentRMS,VyStruct.CongruentRMS] = ...
    ttestFormat(VyStruct.rms,congruence,'congruence');
[accelYStruct.IncongruentRMS,accelYStruct.CongruentRMS] = ...
    ttestFormat(accelYStruct.rms,congruence,'congruence');
[YawRateStruct.IncongruentRMS,YawRateStruct.CongruentRMS] = ...
    ttestFormat(YawRateStruct.rms,congruence,'congruence');

%% Find incorrect/correct trials (based on direction trigger)
%last update ZS 06.22.16

IncorrectCount = 0;
for i = 1:size(Subjects,1)
    for j = 1:Subjects(i).nTrials
        if EligibleTrials(i,j)
            %find peaks in lateral position to determine where the first lane
            %change manuever takes place and save this index.
            %access handwheel position at this same index and check sign of
            %that handwheel signal
            
            %figure,plot(PosY{i,j})
            %DirEvent.Idx(i,j)
            %EndEvent.Idx(i,j)
            [AllPks{i,j},AllLocs{i,j}] = findpeaks(abs(PosY{i,j}(DirEvent.Idx(i,j):EndEvent.Idx(i,j))),...
                'MinPeakHeight',1); %find peaks of 1 meter and greater (absolute value used here)
            
            AllLocs{i,j} = AllLocs{i,j} + DirEvent.Idx(i,j); %add direction event idx to account for index offset
            
            %access hwposition vector for each trial using the 1st found location
            %and determine whether negative or positive peak
            
            PosYPkIdx(i,j) = AllLocs{i,j}(1);
            HWTurnPkVal = HWPosition{i,j}(AllLocs{i,j}(1));
            
            if HWTurnPkVal < 0 %if hw position is negative
                TurnPeakSign(i,j)  = 1; %participant turned right
            elseif HWTurnPkVal > 0 %if hw position is positive
                TurnPeakSign(i,j)  = -1; %participant turned left
            end
            
            %Compare hwposition sign, considering congruence and turn
            %direction for all 4 cases
            if congruence(i,j) == 1 && TurnDirection(i,j) == 1 %congruent, right
                %we expect to see TurnPeakSign = 1 and a negative PosYPkVal
                if (TurnPeakSign(i,j) ~= 1)
                    CorrectTurn(i,j) = 0;
                    IncorrectCount = IncorrectCount + 1;
                    %disp(['Incorrect turn - subject ',num2str(i),', trial ',num2str(j)]);
                end
            elseif congruence(i,j) == 1 && TurnDirection(i,j) == 0 %congruent, left
                %we expect to see TurnPeakSign = -1 and positive PosYPkVal
                if (TurnPeakSign(i,j) ~= -1)
                    CorrectTurn(i,j) = 0;
                    IncorrectCount = IncorrectCount + 1;
                    %disp(['Incorrect turn - subject ',num2str(i),', trial ',num2str(j)]);
                end
            elseif congruence(i,j) == 0 && TurnDirection(i,j) == 1 %incongruent, right
                %we expect to see TurnPeakSign = -1 and negative PosYPkVal
                if (TurnPeakSign(i,j) ~= -1)
                    CorrectTurn(i,j) = 0;
                    IncorrectCount = IncorrectCount + 1;
                    %disp(['Incorrect turn - subject ',num2str(i),', trial ',num2str(j)]);
                end
            elseif congruence(i,j) == 0 && TurnDirection(i,j) == 0 %incongruent, left
                %we expect to see TurnPeakSign = 1 and positive PosYPkVal
                if (TurnPeakSign(i,j) ~= 1)
                    CorrectTurn(i,j) = 0;
                    IncorrectCount = IncorrectCount + 1;
                    %disp(['Incorrect turn - subject ',num2str(i),', trial ',num2str(j)]);
                end
            end
            %Make the correct turn check more robust - check for return to
            %single lane
            if(PosY{i,j}(EndEvent.Idx(i,j)) > 2 || PosY{i,j}(EndEvent.Idx(i,j)) < -2)
                CorrectTurn(i,j) = 0;
            end
        end
    end
end

%% Test cell for correct/incorrect trial completion - comment out if not debugging
%last update ZS 03.15.16

% testSubject = 10; %
% Trial = 3;
% figure
% plot(PosX{testSubject,Trial},HWPosition{testSubject,Trial});
% hold on
% scatter(PosX{testSubject,Trial}(AllLocs{testSubject,Trial}),HWPosition{testSubject,Trial}(AllLocs{testSubject,Trial}));
% scatter(PosX{testSubject,Trial}(AllLocs{testSubject,Trial}),PosY{testSubject,Trial}(AllLocs{testSubject,Trial}),'r');
% plot(PosX{testSubject,Trial},PosY{testSubject,Trial});
% plot([PosX{testSubject,Trial}(DirEventIdx(testSubject,Trial)) PosX{testSubject,Trial}(DirEventIdx(testSubject,Trial))],[-5 5],'LineWidth',2);
% plot([PosX{testSubject,Trial}(EndEventIdx(testSubject,Trial)) PosX{testSubject,Trial}(EndEventIdx(testSubject,Trial))],[-5 5],'LineWidth',2);
% plot([PosX{testSubject,Trial}(AllLocs{testSubject,Trial}) PosX{testSubject,Trial}(AllLocs{testSubject,Trial})],[-5 5],'Color','r');
% hold off
% title(['Subject: ',num2str(subjectNames(testSubject)),', Trial: ',num2str(Trial),...
%     ' | Turn Direction: ',num2str(TurnDirection(testSubject,Trial)),...
%     ' | Congruence: ',num2str(congruence(testSubject,Trial)),...
%     ' | FirstPkSign: ',num2str(TurnPeakSign(testSubject,Trial))]);

%% Calculate trial yield after filtering
%last update ZS 06.14.16

BadPathCount = numel(find(CorrectTurn == 0));
IneligibleTrialCount = numel(find(EligibleTrials == 0));
TrialsRemaining = TotalTrialCount - IneligibleTrialCount - BadPathCount;
FractionRemaining = (TrialsRemaining / TotalTrialCount) * 100;

disp('**Correct trial ID Complete**');
disp(['Total Trials = ',num2str(TotalTrialCount)]);
disp(['Ineligible Trials = ',num2str(IneligibleTrialCount)]);
disp(['Trials Omitted from correct path check = ',num2str(BadPathCount)]);
disp(['% Yield = ',num2str(FractionRemaining,3)]);
disp(['Trials Remaining = ',num2str(TrialsRemaining)]);

%find how many bad trials for each subject
for i = 1:nSubjects
    BadTrialsCount4Subject(i,1) = numel(find(EligibleTrials(i,:) == 0));
end

%% Separate trials by congruence, turn direction, and eligibility trial status
%last update ZS 06.06.16

%depending on settings, include incorrect turns in analysis or not
if settings.includeIncorrectTurns
    %for each type of trial, find only those for which GoodTrials(x,y) = 1
    %find congruent, right turn trials
    [CongruentRightTrialsX,CongruentRightTrialsY] = find(congruence == 1 & TurnDirection == 1 & EligibleTrials == 1);
    %find congruent, left turn trials
    [CongruentLeftTrialsX,CongruentLeftTrialsY] = find(congruence == 1 & TurnDirection == 0 & EligibleTrials == 1);
    %find incongruent, right turn trials
    [IncongruentRightTrialsX,IncongruentRightTrialsY] = find(congruence == 0 & TurnDirection == 1 & EligibleTrials == 1);
    %find incongruent, left turn trials
    [IncongruentLeftTrialsX,IncongruentLeftTrialsY] = find(congruence == 0 & TurnDirection == 0 & EligibleTrials == 1);
else
    %for each type of trial, find only those for which GoodTrials(x,y) = 1
    %find congruent, right turn trials
    [CongruentRightTrialsX,CongruentRightTrialsY] = find(congruence == 1 & TurnDirection == 1 & EligibleTrials == 1 & CorrectTurn == 1);
    %find congruent, left turn trials
    [CongruentLeftTrialsX,CongruentLeftTrialsY] = find(congruence == 1 & TurnDirection == 0 & EligibleTrials == 1 & CorrectTurn == 1);
    %find incongruent, right turn trials
    [IncongruentRightTrialsX,IncongruentRightTrialsY] = find(congruence == 0 & TurnDirection == 1 & EligibleTrials == 1 & CorrectTurn == 1);
    %find incongruent, left turn trials
    [IncongruentLeftTrialsX,IncongruentLeftTrialsY] = find(congruence == 0 & TurnDirection == 0 & EligibleTrials == 1 & CorrectTurn == 1);
end

TrialMap = [
    {CongruentRightTrialsX,CongruentRightTrialsY};
    {CongruentLeftTrialsX,CongruentLeftTrialsY};
    {IncongruentRightTrialsX,IncongruentRightTrialsY};
    {IncongruentLeftTrialsX,IncongruentLeftTrialsY};
    ];

%% Find time/dist between trigger and first handwheel peak and lateral vehicle position peak
%last update ZS 06.15.16

%look for peaks that are greater in magnitude than some percentage of
%the max absolute value

%peakscale adjusts how sensitive findpeaks is to the height of the peaks
%it sees. this scaling factor is multiplied by the abs value of the largest
%peak
peakScale = 0.5;

%not sure what to make this parameter at the moment
MinPeakDist = 100;

time2HWPeak_temp = NaN(size(PosX));
dist2HWPeak_temp = NaN(size(PosX));

time2PosYPeak_temp = NaN(size(PosX));
dist2PosYPeak_temp = NaN(size(PosX));

group1Cutoff = 16;

for i = 1:nSubjects
    %counter for assigning group number
    groupCounter = 0;
    for j = 1:Subjects(i).nTrials
        if EligibleTrials(i,j) == 1
            if CorrectTurn(i,j) == 1
                groupCounter = groupCounter + 1;
                if(groupCounter <= group1Cutoff)
                    groupNum(i,j) = 1;
                else
                    groupNum(i,j) = 2;
                end
                
                %search for both positive and negative peaks
                %check if this trial is a left or right turn event
                if (TurnDirection(i,j) == 1) %right turn, check for negative handwheel position peak
                    MinPeakHeight_neg = abs(min(HWPosition{i,j}(DirEvent.Idx(i,j):end))) * peakScale;
                    [HWPks{i,j},HWLocs{i,j}] = findpeaks(abs(HWPosition{i,j}(DirEvent.Idx(i,j):end)),...
                        'MinPeakHeight',MinPeakHeight_neg,'MinPeakDistance',MinPeakDist);
                    HWPks{i,j} = -1 * HWPks{i,j};
                elseif (TurnDirection(i,j) == 0) %left turn, check for positive handwheel position peak
                    MinPeakHeight_pos = max(HWPosition{i,j}(DirEvent.Idx(i,j):end)) * peakScale;
                    [HWPks{i,j},HWLocs{i,j}] = findpeaks(HWPosition{i,j}(DirEvent.Idx(i,j):end),...
                        'MinPeakHeight',MinPeakHeight_pos,'MinPeakDistance',MinPeakDist);
                end
                
                HWLocs{i,j} = HWLocs{i,j} + DirEvent.Idx(i,j); %add findpeaks search index to location vector
                
                %                 if((i == 10 && j == 10) || ((i == 10 && j == 19)))
                %                    figure
                %                    plot(HWPosition{i,j});
                %                    hold on
                %                    plot(abs(HWPosition{i,j}));
                %                    scatter(HWLocs{i,j}(1),HWPks{i,j}(1));
                %                    plot([DirEvent.Idx(i,j) DirEvent.Idx(i,j)],[-5 5]);
                %                    hold off
                %                 end
                
                %mark the index of the first peak in that series
                HWTurnPkIdx(i,j) = HWLocs{i,j}(1);
                HWTurnPkPos(i,j) = PosX{i,j}(HWTurnPkIdx(i,j));
                HWTurnPkTime(i,j) = time{i,j}(HWTurnPkIdx(i,j));
                
                %calculate time between direction trigger and first handwheel peak
                time2HWPeak_temp(i,j) = HWTurnPkTime(i,j) - DirEvent.Time(i,j);
                time2HWPeak.matrix(i,j) = time2HWPeak_temp(i,j);
                
                %calculate distance between direction trigger and first handwheel peak
                dist2HWPeak_temp(i,j) = HWTurnPkPos(i,j) - DirEvent.Pos(i,j);
                dist2HWPeak.matrix(i,j) = dist2HWPeak_temp(i,j);
                
                %calculate time between direction trigger and lateral
                %vehicle position
                time2PosYPeak_temp(i,j) = time{i,j}(PosYPkIdx(i,j)) - DirEvent.Time(i,j);
                time2PosYPeak.matrix(i,j) = time2PosYPeak_temp(i,j);
                
                %calculate (longitudinal) distance between direction trigger and lateral
                %vehicle position
                dist2PosYPeak_temp(i,j) = PosX{i,j}(PosYPkIdx(i,j)) - DirEvent.Pos(i,j);
                dist2PosYPeak.matrix(i,j) = dist2PosYPeak_temp(i,j);
            end
        end
    end
end

time2HWPeak.matrix(time2HWPeak.matrix == 0) = NaN;
dist2HWPeak.matrix(dist2HWPeak.matrix == 0) = NaN;
time2PosYPeak.matrix(time2PosYPeak.matrix == 0) = NaN;
dist2PosYPeak.matrix(dist2PosYPeak.matrix == 0) = NaN;

[time2HWPeak.IncongruentMean,time2HWPeak.CongruentMean] = ...
    ttestFormat(time2HWPeak.matrix,congruence,'congruence');
[dist2HWPeak.IncongruentMean,dist2HWPeak.CongruentMean] = ...
    ttestFormat(dist2HWPeak.matrix,congruence,'congruence');
[time2PosYPeak.IncongruentMean,time2PosYPeak.CongruentMean] = ...
    ttestFormat(time2PosYPeak.matrix,congruence,'congruence');
[dist2PosYPeak.IncongruentMean,dist2PosYPeak.CongruentMean] = ...
    ttestFormat(dist2PosYPeak.matrix,congruence,'congruence');

time2HWPeak.OverallCongruentMean = mean(time2HWPeak.CongruentMean,'omitnan');
dist2HWPeak.OverallCongruentMean = mean(dist2HWPeak.CongruentMean,'omitnan');

time2HWPeak.OverallIncongruentMean = mean(time2HWPeak.IncongruentMean,'omitnan');
dist2HWPeak.OverallIncongruentMean = mean(dist2HWPeak.IncongruentMean,'omitnan');

time2PosYPeak.OverallCongruentMean = mean(time2PosYPeak.CongruentMean,'omitnan');
dist2PosYPeak.OverallCongruentMean = mean(dist2PosYPeak.CongruentMean,'omitnan');

time2PosYPeak.OverallIncongruentMean = mean(time2PosYPeak.IncongruentMean,'omitnan');
dist2PosYPeak.OverallIncongruentMean = mean(dist2PosYPeak.IncongruentMean,'omitnan');

%box plot time between direction trigger and handwheel peak
overallMeans = [time2HWPeak.OverallIncongruentMean time2HWPeak.OverallCongruentMean];
mytitle = 'Time Between Direction Trigger & First Steering Peak';
metricLabel = 'Time (s)';
BoxPlotMetric(time2HWPeak.matrix,overallMeans,mytitle,metricLabel,congruence,[]);

%box plot distance between direction trigger and handwheel peak
overallMeans = [dist2HWPeak.OverallIncongruentMean dist2HWPeak.OverallCongruentMean];
mytitle = 'Distance Between Direction Trigger & First Steering Peak';
metricLabel = 'Distance (m)';
BoxPlotMetric(dist2HWPeak.matrix,overallMeans,mytitle,metricLabel,congruence,[]);

%box plot time between direction trigger and lateral position peak
overallMeans = [time2PosYPeak.OverallIncongruentMean time2PosYPeak.OverallCongruentMean];
mytitle = 'Time Between Direction Trigger & First Vehicle Lateral Position Peak';
metricLabel = 'Time (s)';
BoxPlotMetric(time2PosYPeak.matrix,overallMeans,mytitle,metricLabel,congruence,[]);

%box plot distance between direction trigger and lateral position peak
overallMeans = [dist2PosYPeak.OverallIncongruentMean dist2PosYPeak.OverallCongruentMean];
mytitle = 'Distance Between Direction Trigger & First Vehicle Lateral Position Peak';
metricLabel = 'Distance (m)';
BoxPlotMetric(dist2PosYPeak.matrix,overallMeans,mytitle,metricLabel,congruence,[]);

%split metric based on group number
time2HWPeak.groupSplit = metricSplitter(time2HWPeak.matrix,groupNum);
dist2HWPeak.groupSplit = metricSplitter(time2HWPeak.matrix,groupNum);
time2PosYPeak.groupSplit = metricSplitter(time2PosYPeak.matrix,groupNum);
dist2PosYPeak.groupSplit = metricSplitter(dist2PosYPeak.matrix,groupNum);

%% Metric: Hand Wheel Reversals
% last update ZS 06.15.16
NumHWReversals = NaN(size(PosX));
HWReversal = cell(size(PosX));

%calculate hw reversals for each type of trial
for n = 1:size(TrialMap,1)
    for i = 1:length(TrialMap{n,1})
        x = TrialMap{n,1}(i);
        y = TrialMap{n,2}(i);
        NumHWReversals(x,y) = sum(HWPosition{x,y} .* [0;HWPosition{x,y}(1:end-1)] < 0);
        %calculate number of reversals per second
        deltaT(x,y) = time{x,y}(end) - time{x,y}(1);
        ReversalsPerSec(x,y) = NumHWReversals(x,y) / deltaT(x,y);
    end
end

CongruentReversals_num = NumHWReversals(congruence == 1);
IncongruentReversals_num = NumHWReversals(congruence == 0);
CongruentReversalsPerSec_num = ReversalsPerSec(congruence == 1);
IncongruentReversalsPerSec_num = ReversalsPerSec(congruence == 0);

reversals.num = NumHWReversals;
reversals.numPerSec = ReversalsPerSec;
reversals.OverallCongruentMeanNum = mean(CongruentReversals_num,'omitnan');
reversals.OverallIncongruentMeanNum = mean(IncongruentReversals_num,'omitnan');
reversals.OverallCongruentMedianNum = median(CongruentReversals_num,'omitnan');
reversals.OverallIncongruentMedianNum = median(IncongruentReversals_num,'omitnan');
reversals.OverallCongruentStdNum = std(CongruentReversals_num,1,'omitnan');
reversals.OverallIncongruentStdNum = std(IncongruentReversals_num,1,'omitnan');

reversals.OverallCongruentMeanNumPerSec = mean(CongruentReversalsPerSec_num,'omitnan');
reversals.OverallIncongruentMeanNumPerSec = mean(IncongruentReversalsPerSec_num,'omitnan');
reversals.OverallCongruentMedianNumPerSec = median(CongruentReversalsPerSec_num,'omitnan');
reversals.OverallIncongruentMedianNumPerSec = median(IncongruentReversalsPerSec_num,'omitnan');
reversals.OverallCongruentStdNumPerSec = std(CongruentReversalsPerSec_num,1,'omitnan');
reversals.OverallIncongruentStdNumPerSec = std(IncongruentReversalsPerSec_num,1,'omitnan');

[reversals.IncongruentMeanNum,reversals.CongruentMeanNum] = ...
    ttestFormat(reversals.num,congruence,'congruence');
[reversals.IncongruentMeanNumPerSec,reversals.CongruentMeanNumPerSec] = ...
    ttestFormat(reversals.numPerSec,congruence,'congruence');

[reversals.Group1MeanNum,reversals.Group2MeanNum] = ...
    ttestFormat(reversals.num,groupNum,'group');
[reversals.Group1MeanNumPerSec,reversals.Group2MeanNumPerSec] = ...
    ttestFormat(reversals.numPerSec,groupNum,'group');

%calculate mean, standard deviation, median number of reversals for congruent/incongruent trials
disp('%%%%%%%%%%%%%%%%%');
disp(['Mean # reversals, congruent = ',num2str(reversals.OverallCongruentMeanNum)]);
disp(['Mean # reversals, incongruent = ',num2str(reversals.OverallIncongruentMeanNum)]);
disp(['Median # reversals, congruent = ',num2str(reversals.OverallCongruentMedianNum)]);
disp(['Median # reversals, incongruent = ',num2str(reversals.OverallIncongruentMedianNum)]);
disp(['Reversal standard deviation, congruent = ',num2str(reversals.OverallCongruentStdNum)]);
disp(['Reversal standard deviation, incongruent = ',num2str(reversals.OverallIncongruentStdNum)]);

%split metric based on group number
reversals.groupSplit_num = metricSplitter(reversals.num,groupNum);
reversals.groupSplit_numPerSec = metricSplitter(reversals.numPerSec,groupNum);

%define mean trigger locations for later use
meanDirEventPos = mean(mean(DirEvent.Pos,'omitnan'),'omitnan');
meanDirEventIdx = mean(mean(DirEvent.Idx,'omitnan'),'omitnan');
meanEndEventPos = mean(mean(EndEvent.Pos,'omitnan'),'omitnan');

%% Metric: Vehicle Position
%last update ZS 06.15.16

FontSize = 14;
AxesMag = 7;
TextHeight = 6;

%plot all position vectors
figure
hold all
for i = 1:nSubjects
    for j = 1:Subjects(i).nTrials
        if size(PosX{i,j},1) ~= 500
            if(EligibleTrials(i,j) == 1 && CorrectTurn(i,j) == 1)
                plot(PosX{i,j},PosY{i,j})
            end
        end
    end
end

plot([meanDirEventPos meanDirEventPos],[(-1*AxesMag) AxesMag],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[(-1*AxesMag) AxesMag],'r--','LineWidth',2);

%add labels for direction and end triggers
text(181.5,TextHeight,'Direction Trigger','FontSize',12);
text(264,TextHeight,'End Trigger','FontSize',12);

%plot cones
s = scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
hold off

axis([150 300 (-1*AxesMag) AxesMag]);

title('Vehicle Position, All Participants','FontSize',FontSize)
xlabel('Longitudinal Position (m)','FontSize',FontSize);
ylabel('Lateral Position (m)','FontSize',FontSize);

%%

vehPosLegendLoc = 'southwest';
xSpacing = linspace(cones.X(1),cones.X(end),snipLength);
conePath.rightlower = polyval(cones.p1,xSpacing);
conePath.rightupper = polyval(cones.p2,xSpacing);
conePath.rightMid = (conePath.rightlower + conePath.rightupper) / 2;

conePath.leftlower = polyval(cones.p3,xSpacing);
conePath.leftupper = polyval(cones.p4,xSpacing);
conePath.leftMid = (conePath.leftlower + conePath.leftupper) / 2;

% 'Processed_PosX' and 'Processed_PosY'
% Each row contains three columns of cells for each trial as ordered in 'TrialMap'
% The columns are
% 1 - indexed
% 2 - snipped
% 3 - snipped_concat
Processed_PosX = SignalProcessor(TrialMap, PosX, DirEvent, EndEvent);
Processed_PosY = SignalProcessor(TrialMap, PosY, DirEvent, EndEvent);

% Assigning variables
% CONGRUENT RIGHT
[PosXStruct.CongruentRightSnipped,PosXStruct.CongruentRightConcat] = Processed_PosX{1,2:3};
[PosYStruct.CongruentRightSnipped,PosYStruct.CongruentRightConcat] = Processed_PosY{1,2:3};
PosXStruct.CongruentRightSnippedMean = mean(PosXStruct.CongruentRightConcat,1);
PosYStruct.CongruentRightSnippedMean = mean(PosYStruct.CongruentRightConcat,1);

% CONGRUENT LEFT
[PosXStruct.CongruentLeftSnipped,PosXStruct.CongruentLeftConcat] = Processed_PosX{2,2:3};
[PosYStruct.CongruentLeftSnipped,PosYStruct.CongruentLeftConcat] = Processed_PosY{2,2:3};
PosXStruct.CongruentLeftSnippedMean = mean(PosXStruct.CongruentLeftConcat,1);
PosYStruct.CongruentLeftSnippedMean = mean(PosYStruct.CongruentLeftConcat,1);

%INCONGRUENT RIGHT
[PosXStruct.IncongruentRightSnipped,PosXStruct.IncongruentRightConcat] = Processed_PosX{3,2:3};
[PosYStruct.IncongruentRightSnipped,PosYStruct.IncongruentRightConcat] = Processed_PosY{3,2:3};
PosXStruct.IncongruentRightSnippedMean = mean(PosXStruct.IncongruentRightConcat,1);
PosYStruct.IncongruentRightSnippedMean = mean(PosYStruct.IncongruentRightConcat,1);

%INCONGRUENT LEFT
[PosXStruct.IncongruentLeftSnipped,PosXStruct.IncongruentLeftConcat] = Processed_PosX{4,2:3};
[PosYStruct.IncongruentLeftSnipped,PosYStruct.IncongruentLeftConcat] = Processed_PosY{4,2:3};
PosXStruct.IncongruentLeftSnippedMean = mean(PosXStruct.IncongruentLeftConcat,1);
PosYStruct.IncongruentLeftSnippedMean = mean(PosYStruct.IncongruentLeftConcat,1);

% Plotting the signals
for i = 1 : length( TrialMap(:,1) )
    %plotting the cones
    figure
    hold on
    
    %plot cones
    s = scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
    scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
    scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
    scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
    
    %plot splined path for left turn
    if mod(i,2) == 0
        h_pathEdge = plot(xSpacing,conePath.leftlower,'r','LineWidth',1.5);
        plot(xSpacing,conePath.leftupper,'r','LineWidth',1.5);
    else
        h_pathEdge = plot(xSpacing,conePath.rightlower,'r','LineWidth',1.5);
        plot(xSpacing,conePath.rightupper,'r','LineWidth',1.5);
    end
    
    
    for k = 1:length(TrialMap{i,1})
        % Plotting the 'snipped' signals
        plot(Processed_PosX{i,2}{k,1}, Processed_PosY{i,2}{k,1});
    end
    
    %plot direction event and end event markers
    plot([DirEvent.Pos(TrialMap{i,1}(2),TrialMap{i,2}(2)),...
        DirEvent.Pos(TrialMap{i,1}(2),TrialMap{i,2}(2))],[-10, 10],...
        'LineWidth',2,'Color','k','LineStyle','--');
    plot([EndEvent.Pos(TrialMap{i,1}(2),TrialMap{i,2}(2)), ...
        EndEvent.Pos(TrialMap{i,1}(2),TrialMap{i,2}(2))],[-10, 10],...
        'LineWidth',2,'Color','r','LineStyle','--');
    
    %add labels for direction and end triggers
    text(182,8,'Direction Trigger');
    text(268,8,'End Trigger');
    
    %plot middle of spline (average of upper and lower)
    if mod(i,2) == 0
        h_midLine = plot(xSpacing,conePath.leftMid,'b','LineWidth',1.5);
    else
        h_midLine = plot(xSpacing,conePath.rightMid,'b','LineWidth',1.5);
    end
    
    %base axis frame size on trigger locations
    axis([(DirEvent.Pos(TrialMap{i,1}(1),TrialMap{i,2}(1)) - 20)...
        (EndEvent.Pos(TrialMap{i,1}(1),TrialMap{i,2}(1)) + 20) -10 10]);
    
    %Title the plot
    switch i
        case 1
            title('All Trials Vehicle Position | Congruent, Right Turn','FontSize',FontSize);
        case 2
            title('All Trials Vehicle Position | Congruent, Left Turn','FontSize',FontSize);
        case 3
            title('All Trials Vehicle Position | Incongruent, Right Turn','FontSize',FontSize);
        case 4
            title('All Trials Vehicle Position | Incongruent, Left Turn','FontSize',FontSize);
    end
    
    xlabel('Longitudinal Position (m)','FontSize',FontSize);
    ylabel('Lateral Position (m)','FontSize',FontSize);
    legend([h_pathEdge,h_midLine,s],'Path Edge','Center line','Cone','Location',vehPosLegendLoc);
    hold off;
end

%% Plot mean vehicle position
% last update ZS 06.15.16

%PLOT CONGRUENT VS INCONGRUENT MEAN POSITION FOR LEFT TURN
%find mean path (snipped) and plot using shadederrorbar for congruent left
figure
h_ConLeftPos = shadedErrorBar(PosXStruct.CongruentLeftSnippedMean,PosYStruct.CongruentLeftSnippedMean,...
    std(PosYStruct.CongruentLeftConcat,0,1),'b',1);
hold on

%plot partial cones field
% s = scatter(cones.X(3:9),cones.YRightLower(3:9),[],[1 .35 0],'filled');
% scatter(cones.X(3:9),cones.YRightUpper(3:9),[],[1 .35 0],'filled');
% scatter(cones.X(3:9),cones.YLeftLower(3:9),[],[1 .35 0],'filled');
% scatter(cones.X(3:9),cones.YLeftUpper(3:9),[],[1 .35 0],'filled');
%plot cones
s = scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');

%find mean path (snipped) and plot using shadederrorbar for congruent right
h_InconLeftPos = shadedErrorBar(PosXStruct.CongruentLeftSnippedMean,PosYStruct.IncongruentLeftSnippedMean,...
    std(PosYStruct.IncongruentLeftConcat,0,1),'r',1);

%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-6 6],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-6 6],'b--','LineWidth',1);
text(195,4,'Single Lane');
text(235,4.5,'Double Lane');
text(267,4,'Single Lane');

%add direction trigger
% plot([meanDirEventPos meanDirEventPos],[-6 6],'k--','LineWidth',2);
% plot([meanEndEventPos meanEndEventPos],[-6 6],'r--','LineWidth',2);

%add labels for direction and end triggers
%text(181,5,'Direction Trigger');
%text(272,5,'End Trigger');

% hold off
% title('Mean Vehicle Position | Left Turn','FontSize',FontSize);
% xlabel('Longitudinal Position (m)','FontSize',FontSize);
% ylabel('Lateral Position (m)','FontSize',FontSize);
% legend([h_ConLeftPos.mainLine,h_InconLeftPos.mainLine,s],'Congruent',...
%     'Incongruent','Cone','Location','Southeast');

%PLOT CONGRUENT VS INCONGRUENT MEAN POSITION FOR RIGHT TURN
%find mean path (snipped) and plot using shadederrorbar for congruent left
%figure
h_ConRightPos = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,PosYStruct.CongruentRightSnippedMean,...
    std(PosYStruct.CongruentRightConcat,0,1),'b',1);
%hold on

%plot partial cones
% scatter(cones.X(3:9),cones.YRightLower(3:9),[],[1 .35 0],'filled');
% scatter(cones.X(3:9),cones.YRightUpper(3:9),[],[1 .35 0],'filled');
% scatter(cones.X(3:9),cones.YLeftLower(3:9),[],[1 .35 0],'filled');
% scatter(cones.X(3:9),cones.YLeftUpper(3:9),[],[1 .35 0],'filled');

%plot cones
s = scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');

%find mean path (snipped) and plot using shadederrorbar for congruent right
h_InconRightPos = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,PosYStruct.IncongruentRightSnippedMean,...
    std(PosYStruct.IncongruentRightConcat,0,1),'r',1);

%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-6 6],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-6 6],'b--','LineWidth',1);
text(195,4,'Single Lane');
text(235,4.5,'Double Lane');
text(267,4,'Single Lane');

%add direction trigger
plot([meanDirEventPos meanDirEventPos],[-6 6],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-6 6],'r--','LineWidth',2);

%add labels for direction and end triggers
text(181,5,'Direction Trigger');
text(267,5,'End Trigger');

hold off

axis([162 300 -6 6]);

%title('Mean Vehicle Position | Right Turn','FontSize',FontSize);
title('Mean Vehicle Position','FontSize',FontSize);
xlabel('Longitudinal Position (m)','FontSize',FontSize);
ylabel('Lateral Position (m)','FontSize',FontSize);
% legend([h_ConRightPos.mainLine,h_InconRightPos.mainLine,s],'Congruent',...
%     'Incongruent','Cone','Location','Southeast');
legend([h_ConRightPos.mainLine,h_InconRightPos.mainLine,s],'Congruent',...
    'Incongruent','Cone','Location','Southeast');

%% Find subject-wise congruent paths for relative error metrics (lateral and longitudinal)
%last update ZS 06.13.16

PosYStruct.RelConRightPathMeans = cell(nSubjects,1);
PosXStruct.RelConRightPathMeans = cell(nSubjects,1);
PosYStruct.RelConRightPathUpper = cell(nSubjects,1);
PosYStruct.RelConRightPathLower = cell(nSubjects,1);

PosYStruct.RelConLeftPathMeans = cell(nSubjects,1);
PosXStruct.RelConLeftPathMeans = cell(nSubjects,1);
PosYStruct.RelConLeftPathUpper = cell(nSubjects,1);
PosYStruct.RelConLeftPathLower = cell(nSubjects,1);

for i = 1:nSubjects
    ConRightPosY_temp = [];
    ConLeftPosY_temp = [];
    ConRightPosX_temp = [];
    ConLeftPosX_temp = [];
    
    %we already found all congruent right and left trials
    %find where current subject number is in these trial maps
    trialSearch_ConRight = find(CongruentRightTrialsX == i);
    trialSearch_ConLeft = find(CongruentLeftTrialsX == i);
    
    for j = 1:length(trialSearch_ConRight)
        ConRightTrial_PosY = PosYStruct.CongruentRightSnipped{trialSearch_ConRight(j)};
        ConRightTrial_PosX = PosXStruct.CongruentRightSnipped{trialSearch_ConRight(j)};
        %concatenate these trials
        ConRightPosY_temp = [ConRightPosY_temp;ConRightTrial_PosY'];
        ConRightPosX_temp = [ConRightPosX_temp;ConRightTrial_PosX'];
    end
    
    for j = 1:length(trialSearch_ConLeft)
        ConLeftTrial_PosY = PosYStruct.CongruentLeftSnipped{trialSearch_ConLeft(j)};
        ConLeftTrial_PosX = PosXStruct.CongruentLeftSnipped{trialSearch_ConLeft(j)};
        %concatenate these trials
        ConLeftPosY_temp = [ConLeftPosY_temp;ConLeftTrial_PosY'];
        ConLeftPosX_temp = [ConLeftPosX_temp;ConLeftTrial_PosX'];
    end
    
    %mean these trials and calculate upper and lower path boundaries
    PosYStruct.RelConRightPathMeans{i} = mean(ConRightPosY_temp,1)';
    PosXStruct.RelConRightPathMeans{i} = mean(ConRightPosX_temp,1)';
    PosYStruct.RelConRightPathUpper{i} = PosYStruct.RelConRightPathMeans{i}...
        + cones.LaneWidth / 2;
    PosYStruct.RelConRightPathLower{i} = PosYStruct.RelConRightPathMeans{i}...
        - cones.LaneWidth / 2;
    
    PosYStruct.RelConLeftPathMeans{i} = mean(ConLeftPosY_temp,1)';
    PosXStruct.RelConLeftPathMeans{i} = mean(ConLeftPosX_temp,1)';
    PosYStruct.RelConLeftPathUpper{i} = PosYStruct.RelConLeftPathMeans{i}...
        + cones.LaneWidth / 2;
    PosYStruct.RelConLeftPathLower{i} = PosYStruct.RelConLeftPathMeans{i}...
        - cones.LaneWidth / 2;
end

%package path means and bounds
RelPathPack.means = {
    PosXStruct.RelConRightPathMeans,PosYStruct.RelConRightPathMeans,...
    PosXStruct.RelConLeftPathMeans,PosYStruct.RelConLeftPathMeans
    };
RelPathPack.bounds = {
    PosYStruct.RelConRightPathLower,PosYStruct.RelConRightPathUpper,...
    PosYStruct.RelConLeftPathLower,PosYStruct.RelConLeftPathUpper
    };

%% Metric: Lat/long error between lateral/long position and subject-wise congruent average lat/long position
%last update ZS 06.21.16

disp('Calculating lateral and longitudinal error...');

latError.matrix = cell(size(PosX));
latError.mean = NaN(size(PosX));
longError.matrix = cell(size(PosX));
longError.mean = NaN(size(PosX));
%longError.longLag = cell(size(PosX));

doubleLaneStartIdx = NaN(size(PosX));
singleLaneStartIdx = NaN(size(PosX));

ConLatError_concat = [];
ConLongError_concat = [];
InconLatError_concat = [];
InconLongError_concat = [];

for n = 1:size(TrialMap,1)
    %Congruent, right trials
    for i = 1:length(TrialMap{n,1})
        x = TrialMap{n,1}(i);
        y = TrialMap{n,2}(i);
        
        [~, doubleLaneStartIdx(x,y)] = findNearest(PosX{x,y},cones.X(4));
        [~, singleLaneStartIdx(x,y)] = findNearest(PosX{x,y},cones.X(7));
        
        if settings.choosePathDev %if relative chosen
            if ~mod(n,2) %if this is a left turn
                latError_temp = abs(Processed_PosY{n,2}{i} - PosYStruct.RelConLeftPathMeans{x});
                [longError_temp,halves] = ...
                    findLag(PosXStruct.RelConLeftPathMeans{x},PosYStruct.RelConLeftPathMeans{x},...
                    Processed_PosX{n,2}{i},Processed_PosY{n,2}{i});
            else
                latError_temp = abs(Processed_PosY{n,2}{i} - PosYStruct.RelConRightPathMeans{x});
                [longError_temp,halves] = ...
                    findLag(PosXStruct.RelConRightPathMeans{x},PosYStruct.RelConRightPathMeans{x},...
                    Processed_PosX{n,2}{i},Processed_PosY{n,2}{i});
            end
        else %if absolute chosen (longitudinal computes in same way as relative currently)
            if ~mod(n,2) %if this is a left turn
                latError_temp = abs(Processed_PosY{n,2}{i} - conePath.leftMid');
                [longError_temp,halves] = ...
                    findLag(PosXStruct.RelConLeftPathMeans{x},PosYStruct.RelConLeftPathMeans{x},...
                    Processed_PosX{n,2}{i},Processed_PosY{n,2}{i});
            else
                latError_temp = abs(Processed_PosY{n,2}{i} - conePath.rightMid');
                [longError_temp,halves] = ...
                    findLag(PosXStruct.RelConRightPathMeans{x},PosYStruct.RelConRightPathMeans{x},...
                    Processed_PosX{n,2}{i},Processed_PosY{n,2}{i});
            end
        end
        if ~isempty(longError_temp)
            if n <= 2 %if we are processing the congruent trials
                ConLatError_concat = [ConLatError_concat;latError_temp'];
                ConLongError_concat = [ConLongError_concat;longError_temp'];
            else %if we are processing the incongruent trials
                InconLatError_concat = [InconLatError_concat;latError_temp'];
                InconLongError_concat = [InconLongError_concat;longError_temp'];
            end
            latError.matrix{x,y} = latError_temp;
            latError.mean(x,y) = mean(latError.matrix{x,y},'omitnan');
            longError.matrix{x,y} = longError_temp;
            longError.mean(x,y) = mean(longError.matrix{x,y},'omitnan');
        end
    end
end

%% Split lat/long into subject-wise mean for stats

[latError.IncongruentMean,latError.CongruentMean] = ...
    ttestFormat(latError.mean,congruence,'congruence');
[latError.Group1Mean,latError.Group2Mean] = ...
    ttestFormat(latError.mean,groupNum,'group');
[longError.IncongruentMean,longError.CongruentMean] = ...
    ttestFormat(longError.mean,congruence,'congruence');
[longError.Group1Mean,longError.Group2Mean] = ...
    ttestFormat(longError.mean,groupNum,'group');

%find mean of congruent and incongruent lateral/longitudinal error, mean along 1st
%dimension
latError.OverallCongruentMean = mean(ConLatError_concat,1)';
latError.OverallCongruentStd = std(ConLatError_concat,1)';

longError.OverallCongruentMean = mean(ConLongError_concat,1)';
longError.OverallCongruentStd = std(ConLongError_concat,1)';

latError.OverallIncongruentMean = mean(InconLatError_concat,1)';
latError.OverallIncongruentStd = std(InconLatError_concat,1)';

longError.OverallIncongruentMean = mean(InconLongError_concat,1)';
longError.OverallIncongruentStd = std(InconLongError_concat,1)';

meanDoubleLaneStartIdx = floor(mean(mean(doubleLaneStartIdx,'omitnan')) - meanDirEventIdx);
meanSingleLaneStartIdx = floor(mean(mean(singleLaneStartIdx,'omitnan')) - meanDirEventIdx);

%split metric based on group number
latError.groupSplit = metricSplitter(latError.mean,groupNum);
longError.groupSplit = metricSplitter(longError.mean,groupNum);

%% Split lateral and longitudinal error for each section of the trial
% last update ZS 08.08.16

latError.CongruentMeanSection1 = ...
    mean(latError.OverallCongruentMean(1:meanDoubleLaneStartIdx),'omitnan');
latError.CongruentMeanSection2 = ...
    mean(latError.OverallCongruentMean(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
latError.CongruentMeanSection3 = ...
    mean(latError.OverallCongruentMean(meanSingleLaneStartIdx:end),'omitnan');

latError.CongruentStdSection1 = ...
    mean(latError.OverallCongruentStd(1:meanDoubleLaneStartIdx),'omitnan');
latError.CongruentStdSection2 = ...
    mean(latError.OverallCongruentStd(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
latError.CongruentStdSection3 = ...
    mean(latError.OverallCongruentStd(meanSingleLaneStartIdx:end),'omitnan');

latError.IncongruentMeanSection1 = ...
    mean(latError.OverallIncongruentMean(1:meanDoubleLaneStartIdx),'omitnan');
latError.IncongruentMeanSection2 = ...
    mean(latError.OverallIncongruentMean(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
latError.IncongruentMeanSection3 = ...
    mean(latError.OverallIncongruentMean(meanSingleLaneStartIdx:end),'omitnan');

latError.IncongruentStdSection1 = ...
    mean(latError.OverallIncongruentStd(1:meanDoubleLaneStartIdx),'omitnan');
latError.IncongruentStdSection2 = ...
    mean(latError.OverallIncongruentStd(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
latError.IncongruentStdSection3 = ...
    mean(latError.OverallIncongruentStd(meanSingleLaneStartIdx:end),'omitnan');

%split longitudinal error for each section of the trial
longError.CongruentMeanSection1 = ...
    mean(longError.OverallCongruentMean(1:meanDoubleLaneStartIdx),'omitnan');
longError.CongruentMeanSection2 = ...
    mean(longError.OverallCongruentMean(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
longError.CongruentMeanSection3 = ...
    mean(longError.OverallCongruentMean(meanSingleLaneStartIdx:end),'omitnan');

longError.CongruentStdSection1 = ...
    mean(longError.OverallCongruentStd(1:meanDoubleLaneStartIdx),'omitnan');
longError.CongruentStdSection2 = ...
    mean(longError.OverallCongruentStd(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
longError.CongruentStdSection3 = ...
    mean(longError.OverallCongruentStd(meanSingleLaneStartIdx:end),'omitnan');

longError.IncongruentMeanSection1 = ...
    mean(longError.OverallIncongruentMean(1:meanDoubleLaneStartIdx),'omitnan');
longError.IncongruentMeanSection2 = ...
    mean(longError.OverallIncongruentMean(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
longError.IncongruentMeanSection3 = ...
    mean(longError.OverallIncongruentMean(meanSingleLaneStartIdx:end),'omitnan');

longError.IncongruentStdSection1 = ...
    mean(longError.OverallIncongruentStd(1:meanDoubleLaneStartIdx),'omitnan');
longError.IncongruentStdSection2 = ...
    mean(longError.OverallIncongruentStd(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
longError.IncongruentStdSection3 = ...
    mean(longError.OverallIncongruentStd(meanSingleLaneStartIdx:end),'omitnan');

%plot lateral error (lateral error from congruent mean path)
figure
subplot(2,1,1)
h_ConLatError = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,latError.OverallCongruentMean,...
    std(ConLatError_concat,0,1),'b',1);
hold on
h_InconLatError = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,latError.OverallIncongruentMean,...
    std(InconLatError_concat,0,1),'r',1);

%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-1 3],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-1 3],'b--','LineWidth',1);

%add direction trigger
plot([meanDirEventPos meanDirEventPos],[-1 3],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-1 3],'r--','LineWidth',2);
hold off
title('Mean Lateral Position Error','FontSize',FontSize);
axis([180 300 -1 2]);

ylabel('Lateral Position Error (m)','FontSize',FontSize);
legend([h_ConLatError.mainLine,h_InconLatError.mainLine],'Congruent',...
    'Incongruent','location','southwest');

subplot(2,1,2)
%plot cones
scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
hold on
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
%add direction and end triggers
plot([meanDirEventPos meanDirEventPos],[-15 15],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-15 15],'r--','LineWidth',2);
%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-15 15],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-15 15],'b--','LineWidth',1);
ylabel('Lateral Position (m)','FontSize',FontSize);
xlabel('Longitudinal Position (m)','FontSize',FontSize);
hold off
axis([180 300 -7 7]);
%add labels for direction and end triggers
text(181,6,'Direction Trigger');
text(269,6,'End Trigger');
text(195,4,'Single Lane');
text(235,4,'Double Lane');
text(269,4,'Single Lane');

%plot longitudinal error (error between lateral peak in congruent mean and each trial)
% figure
% h_ConPathLongError = shadedErrorBar([halves{1}.peakPercent_crnt;halves{2}.peakPercent_crnt],longError.CongruentMean,...
%     std(ConLongError_concat,0,1),'b',1);
% hold on
% h_InconPathLongError = shadedErrorBar([halves{1}.peakPercent_crnt;halves{2}.peakPercent_crnt],longError.IncongruentMean,...
%     std(InconLongError_concat,0,1),'r',1);

figure
h_ConPathLongError = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,longError.OverallCongruentMean,...
    std(ConLongError_concat,0,1),'b',1);
hold on
h_InconPathLongError = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,longError.OverallIncongruentMean,...
    std(InconLongError_concat,0,1),'r',1);

%add direction trigger
plot([meanDirEventPos meanDirEventPos],[-20 25],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-20 25],'r--','LineWidth',2);
%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-20 25],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-20 25],'b--','LineWidth',1);
hold off
text(195,18,'Single Lane');
text(235,18,'Double Lane');
text(269,18,'Single Lane');

title('Mean Longitudinal Position Error','FontSize',FontSize);
xlabel('Longitudinal Position (m)','FontSize',FontSize);
ylabel(' Longitudinal Position Error (m)','FontSize',FontSize);
legend([h_ConPathLongError.mainLine,h_InconPathLongError.mainLine],'Congruent',...
    'Incongruent','location','southwest');

%add labels for direction and end triggers
text(181,22,'Direction Trigger');
text(272,22,'End Trigger');

%%
%break down congruent number reversals by trial section
% reversals.CongruentMeanSection1Num = ...
%     mean(reversals.CongruentMeanNum(1:meanDoubleLaneStartIdx),'omitnan');
% reversals.CongruentMeanSection2Num = ...
%     mean(reversals.CongruentMeanNum(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
% reversals.CongruentMeanSection3Num = ...
%     mean(reversals.CongruentMeanNum(meanSingleLaneStartIdx:end),'omitnan');
% 
% %break down congruent reversals per second by trial section
% reversals.CongruentMeanSection1NumPerSec = ...
%     mean(reversals.CongruentMeanNumPerSec(1:meanDoubleLaneStartIdx),'omitnan');
% reversals.CongruentMeanSection2NumPerSec = ...
%     mean(reversals.CongruentMeanNumPerSec(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
% reversals.CongruentMeanSection3NumPerSec = ...
%     mean(reversals.CongruentMeanNumPerSec(meanSingleLaneStartIdx:end),'omitnan');
% 
% %break down incongruent number reversals by trial section
% reversals.IncongruentMeanSection1Num = ...
%     mean(reversals.IncongruentMeanNum(1:meanDoubleLaneStartIdx),'omitnan');
% reversals.IncongruentMeanSection2Num = ...
%     mean(reversals.IncongruentMeanNum(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
% reversals.IncongruentMeanSection3Num = ...
%     mean(reversals.IncongruentMeanNum(meanSingleLaneStartIdx:end),'omitnan');
% 
% %break down incongruent reversals per second by trial section
% reversals.IncongruentMeanSection1NumPerSec = ...
%     mean(reversals.IncongruentMeanNumPerSec(1:meanDoubleLaneStartIdx),'omitnan');
% reversals.IncongruentMeanSection2NumPerSec = ...
%     mean(reversals.IncongruentMeanNumPerSec(meanDoubleLaneStartIdx:meanSingleLaneStartIdx),'omitnan');
% reversals.IncongruentMeanSection3NumPerSec = ...
%     mean(reversals.IncongruentMeanNumPerSec(meanSingleLaneStartIdx:end),'omitnan');

%% Plot cell for all participants' vehicle position and their own congruent mean vehicle position
%last updated ZS 05.31.16

myTitle = cell(21,1);
subjectLatError = cell(21,1);
meanSubjectLatError = cell(21,1);
stdSubjectLatError = cell(21,1);

for n = 1:nSubjects
    %     figure;
    %     hold on
    %     %plot subject's personal congruent baseline path, as well as the
    %     %boundaries of the path by adding/subtracting half of the cones width
    %     plot(PosXStruct.CongruentRightSnippedMean,PosYStruct.RelConRightPathMeans{n},'LineWidth',2);
    %     plot(PosXStruct.CongruentRightSnippedMean,PosYStruct.RelConRightPathMeans{n} + ...
    %         cones.LaneWidth/2,'r','LineWidth',2);
    %     plot(PosXStruct.CongruentRightSnippedMean,PosYStruct.RelConRightPathMeans{n} - ...
    %         cones.LaneWidth/2,'r','LineWidth',2);
    %
    %     for i = 1:length(CongruentRightTrialsX)
    %         x = CongruentRightTrialsX(i);
    %         y = CongruentRightTrialsY(i);
    %         if(CongruentRightTrialsX(i) == n)
    %             plot(PosXStruct.CongruentRightSnipped{i},PosYStruct.CongruentRightSnipped{i});
    %             myTitle{CongruentRightTrialsX(i)} = ...
    %                 ['Congruent Right Vehicle Position, Participant: ',num2str(UserID{x,y}(1))];
    %         end
    %     end
    %
    %     for i = 1:length(IncongruentRightTrialsX)
    %         x = IncongruentRightTrialsX(i);
    %         y = IncongruentRightTrialsY(i);
    %         if(IncongruentRightTrialsX(i) == n)
    %             plot(PosXStruct.IncongruentRightSnipped{i},PosYStruct.IncongruentRightSnipped{i});
    %             myTitle{IncongruentRightTrialsX(i)} = ...
    %                 ['Incongruent Right Vehicle Position, Participant: ',num2str(UserID{x,y}(1))];
    %         end
    %     end
    %
    %     for i = 1:length(CongruentLeftTrialsX)
    %         x = CongruentLeftTrialsX(i);
    %         y = CongruentLeftTrialsY(i);
    %         if(CongruentLeftTrialsX(i) == n)
    %             plot(PosX.CongruentLeftSnipped{i},PosY.CongruentLeftSnipped{i});
    %             myTitle{CongruentLeftTrialsX(i)} = ...
    %                 ['Congruent Left Vehicle Position, Participant: ',num2str(UserID{x,y}(1))];
    %         end
    %     end
    %
    %     for i = 1:length(IncongruentLeftTrialsX)
    %         x = IncongruentLeftTrialsX(i);
    %         y = IncongruentLeftTrialsY(i);
    %         if(IncongruentLeftTrialsX(i) == n)
    %             plot(PosXStruct.IncongruentLeftSnipped{i},PosY.IncongruentLeftSnipped{i});
    %             myTitle{IncongruentLeftTrialsX(i)} = ...
    %                 ['Incongruent Left Vehicle Position, Participant: ',num2str(UserID{x,y}(1))];
    %         end
    %     end
    %
    %     s = scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
    %     scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
    %     scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
    %     scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
    %
    %     hold off
    %     title(myTitle{n});
    %     xlabel('Longitudinal Position (m)');
    %     ylabel('Lateral Position (m)');
    %
    latErrorTemp_con = [];
    latErrorTemp_incon = [];
    for k = 1:Subjects(n).nTrials
        if(~isempty(latError.matrix{n,k}))
            if(congruence(n,k))
                latErrorTemp_con = [latErrorTemp_con latError.matrix{n,k}];
            else
                latErrorTemp_incon = [latErrorTemp_incon latError.matrix{n,k}];
            end
        end
    end
    
    latError.CongruentSubjectPath{n,1} = latErrorTemp_con;
    latError.IncongruentSubjectPath{n,1} = latErrorTemp_incon;
    latError.CongruentSubjectPathMean{n,1} = mean(latError.CongruentSubjectPath{n,1},2);
    latError.IncongruentSubjectPathMean{n,1} = mean(latError.IncongruentSubjectPath{n,1},2);
    latError.CongruentSubjectPathStd{n,1} = std(latError.CongruentSubjectPath{n,1},0,2);
    latError.IncongruentSubjectPathStd{n,1} = std(latError.IncongruentSubjectPath{n,1},0,2);
end

%% shaded error bar plots for lateral error for each subject
% last update 05.23.16

% for i = 1:nSubjects
%     %if(Subjects(i).Names == 20790) %use if statement if focus on one subject
%         figure;
%         h1 = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,latError.CongruentSubjectPathMean{i,1},...
%             latError.CongruentSubjectPathStd{i,1},'b',1);
%         hold on
%         h2 = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,latError.IncongruentSubjectPathMean{i,1},...
%             latError.IncongruentSubjectPathStd{i,1},'r',1);
%         title(['Mean lateral error | Subject ',num2str(Subjects(i).Names)]);
%         xlabel('Longitudinal Vehicle Position (m)');
%         ylabel('Lateral Error (m)');
%
%         %add longitudinal cone section markers
%         plot([cones.X(4) cones.X(4)],[-2 3],'b--','LineWidth',1);
%         plot([cones.X(7) cones.X(7)],[-2 3],'b--','LineWidth',1);
%         text(195,1.5,'Single Lane');
%         text(235,1.5,'Double Lane');
%         text(269,1.5,'Single Lane');
%
%         %add direction trigger
%         plot([meanDirEventPos meanDirEventPos],[-2 3],'k--','LineWidth',2);
%         plot([meanEndEventPos meanEndEventPos],[-2 3],'r--','LineWidth',2);
%         hold off
%
%         %add labels for direction and end triggers
%         text(181,2.25,'Direction Trigger');
%         text(272,2.25,'End Trigger');
%
%         legend([h1.mainLine,h2.mainLine],'Congruent','Incongruent','location',...
%             'southwest');
%     %end
% end

%% Metric: Lateral Velocity
% last update MS 05.16.16

% Process the signals for Lateral Velocity
Processed_Vy = SignalProcessor(TrialMap,Vy,DirEvent,EndEvent);

%Congruent Right
[VyStruct.CongruentRight,VyStruct.CongruentRightSnipped,VyStruct.CongruentRightConcat] = Processed_Vy{1,:};
VyStruct.CongruentRightSnippedMean = mean(VyStruct.CongruentRightConcat,1);

%Congruent Left
[VyStruct.CongruentLeft,VyStruct.CongruentLeftSnipped,VyStruct.CongruentLeftConcat] = Processed_Vy{2,:};
VyStruct.CongruentLeftSnippedMean = mean(VyStruct.CongruentLeftConcat,1);

%Incongruent Right
[VyStruct.IncongruentRight,VyStruct.IncongruentRightSnipped,VyStruct.IncongruentRightConcat] = Processed_Vy{3,:};
VyStruct.IncongruentRightSnippedMean = mean(VyStruct.IncongruentRightConcat,1);

%Incongruent Left
[VyStruct.IncongruentLeft,VyStruct.IncongruentLeftVy_snipped,VyStruct.IncongruentLeftConcat] = Processed_Vy{4,:};
VyStruct.IncongruentLeftSnippedMean = mean(VyStruct.IncongruentLeftConcat,1);

%plot for right turns
figure;
subplot(2,1,1)
h_ConRightVy = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,VyStruct.CongruentRightSnippedMean,...
    std(VyStruct.CongruentRightConcat,0,1),'b',1);
hold on;
h_InconRightVy = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,VyStruct.IncongruentRightSnippedMean,...
    std(VyStruct.IncongruentRightConcat,0,1),'r',1);

%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-0.3 0.3],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-0.3 0.3],'b--','LineWidth',1);
% text(195,0.2,'Single Lane');
% text(235,0.2,'Double Lane');
% text(269,0.2,'Single Lane');

%add direction trigger
plot([meanDirEventPos meanDirEventPos],[-0.3 0.3],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-0.3 0.3],'r--','LineWidth',2);

hold off
title('Mean Lateral Velocity | Right Turn','FontSize',FontSize);
%xlabel('Longitudinal Position (m)');
ylabel('Lateral Velocity (m/s)','FontSize',FontSize);
legend([h_ConRightVy.mainLine,h_InconRightVy.mainLine],'Congruent',...
    'Incongruent','location','southwest');

subplot(2,1,2)
%plot cones
scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
hold on
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
%add direction and end triggers
plot([meanDirEventPos meanDirEventPos],[-15 15],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-15 15],'r--','LineWidth',2);
%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-15 15],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-15 15],'b--','LineWidth',1);
ylabel('Lateral Position (m)','FontSize',FontSize);
xlabel('Longitudinal Position (m)','FontSize',FontSize);
hold off
axis([180 300 -7 7]);

%add labels for direction and end triggers
text(181,6,'Direction Trigger');
text(272,6,'End Trigger');
text(195,4,'Single Lane');
text(235,4,'Double Lane');
text(269,4,'Single Lane');

%plot for left turns
figure;
subplot(2,1,1)
h_ConLeftVy = shadedErrorBar(PosXStruct.CongruentLeftSnippedMean,VyStruct.CongruentLeftSnippedMean,...
    std(VyStruct.CongruentLeftConcat,0,1),'b',1);
hold on
h_InconLeftVy = shadedErrorBar(PosXStruct.CongruentLeftSnippedMean,VyStruct.IncongruentLeftSnippedMean,...
    std(VyStruct.IncongruentLeftConcat,0,1),'r',1);

%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-0.3 0.3],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-0.3 0.3],'b--','LineWidth',1);
% text(195,0.2,'Single Lane');
% text(235,0.2,'Double Lane');
% text(269,0.2,'Single Lane');

%add direction trigger
plot([meanDirEventPos meanDirEventPos],[-0.3 0.3],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-0.3 0.3],'r--','LineWidth',2);

hold off
title('Mean Lateral Velocity | Left Turn','FontSize',FontSize);
%xlabel('Longitudinal Position (m)');
ylabel('Lateral Velocity (m/s)','FontSize',FontSize);
legend([h_ConLeftVy.mainLine,h_InconLeftVy.mainLine],'Congruent',...
    'Incongruent','location','southwest');

subplot(2,1,2)
%plot cones
scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
hold on
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
%add direction and end triggers
plot([meanDirEventPos meanDirEventPos],[-15 15],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-15 15],'r--','LineWidth',2);
%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-15 15],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-15 15],'b--','LineWidth',1);
ylabel('Lateral Position (m)','FontSize',FontSize);
xlabel('Longitudinal Position (m)','FontSize',FontSize);
hold off
axis([180 300 -7 7]);
%add labels for direction and end triggers
text(181,6,'Direction Trigger');
text(272,6,'End Trigger');
text(195,4,'Single Lane');
text(235,4,'Double Lane');
text(269,4,'Single Lane');

%split metric based on group number
%VyStruct.groupSplit = metricSplitter(VyStruct.matrix,groupNum);

%% Metric: Lateral Acceleration
%last update MS 05.16.16

% Process the signals for Lateral Acceleration
Processed_accelY = SignalProcessor(TrialMap,accelY,DirEvent,EndEvent);

%Congruent Right
[accelYStruct.CongruentRight,accelYStruct.CongruentRightSnipped,accelYStruct.CongruentRightConcat] = Processed_accelY{1,:};
accelYStruct.CongruentRightSnippedMean = mean(accelYStruct.CongruentRightConcat,1);

%Congruent Left
[accelYStruct.CongruentLeft,accelYStruct.CongruentLeftSnipped,accelYStruct.CongruentLeftConcat] = Processed_accelY{2,:};
accelYStruct.CongruentLeftSnippedMean = mean(accelYStruct.CongruentLeftConcat,1);

%Incongruent Right
[accelYStruct.IncongruentRight,accelYStruct.IncongruentRightSnipped,accelYStruct.IncongruentRightConcat] = Processed_accelY{3,:};
accelYStruct.IncongruentRightSnippedMean = mean(accelYStruct.IncongruentRightConcat,1);

%Incongruent Left
[accelYStruct.IncongruentLeft,accelYStruct.IncongruentLeftSnipped,accelYStruct.IncongruentLeftConcat] = Processed_accelY{4,:};
accelYStruct.IncongruentLeftSnippedMean = mean(accelYStruct.IncongruentLeftConcat,1);

%UNCOMMENT FOR LATERAL ACCEL PLOTS

%plot for right turns
figure;
subplot(2,1,1)
h_ConRightAccelY = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,accelYStruct.CongruentRightSnippedMean,...
    std(accelYStruct.CongruentRightConcat,0,1),'b',1);
hold on;
h_InconRightAccelY = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,accelYStruct.IncongruentRightSnippedMean,...
    std(accelYStruct.IncongruentRightConcat,0,1),'r',1);
title('Mean Lateral Acceleration | Right Turn','FontSize',FontSize);
%xlabel('Longitudinal Position (m)');
ylabel('Lateral Acceleration (m/s^2)','FontSize',FontSize);
legend([h_ConRightAccelY.mainLine,h_InconRightAccelY.mainLine],'Congruent',...
    'Incongruent','location','southwest');
%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-4 4],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-4 4],'b--','LineWidth',1);

%add direction trigger
plot([meanDirEventPos meanDirEventPos],[-4 4],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-4 4],'r--','LineWidth',2);
hold off

subplot(2,1,2)
%plot cones
scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
hold on
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
%add direction and end triggers
plot([meanDirEventPos meanDirEventPos],[-15 15],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-15 15],'r--','LineWidth',2);
%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-15 15],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-15 15],'b--','LineWidth',1);
ylabel('Lateral Position (m)','FontSize',FontSize);
xlabel('Longitudinal Position (m)','FontSize',FontSize);
hold off
axis([180 300 -7 7]);

%add labels for direction and end triggers
text(181,6,'Direction Trigger');
text(272,6,'End Trigger');
text(195,4,'Single Lane');
text(235,4,'Double Lane');
text(269,4,'Single Lane');

%plot for left turns
figure;
subplot(2,1,1)
h_ConLeftAccelY = shadedErrorBar(PosXStruct.CongruentLeftSnippedMean,accelYStruct.CongruentLeftSnippedMean,...
    std(accelYStruct.CongruentLeftConcat,0,1),'b',1);
hold on;
h_InconLeftAccelY = shadedErrorBar(PosXStruct.CongruentLeftSnippedMean,accelYStruct.IncongruentLeftSnippedMean,...
    std(accelYStruct.IncongruentLeftConcat,0,1),'r',1);

title('Mean Lateral Acceleration | Left Turn','FontSize',FontSize);
%xlabel('Longitudinal Position (m)');
ylabel('Lateral Acceleration (m/s^2)','FontSize',FontSize);
legend([h_ConLeftAccelY.mainLine,h_InconLeftAccelY.mainLine],'Congruent',...
    'Incongruent','location','southwest');

%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-5 5],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-5 5],'b--','LineWidth',1);

%add direction and end triggers
plot([meanDirEventPos meanDirEventPos],[-5 5],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-5 5],'r--','LineWidth',2);
hold off

subplot(2,1,2)
%plot cones
scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
hold on
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
%add direction and end triggers
plot([meanDirEventPos meanDirEventPos],[-15 15],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-15 15],'r--','LineWidth',2);
%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-15 15],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-15 15],'b--','LineWidth',1);
ylabel('Lateral Position (m)','FontSize',FontSize);
xlabel('Longitudinal Position (m)','FontSize',FontSize);
hold off
axis([180 300 -7 7]);

%add labels for direction and end triggers
text(181,6,'Direction Trigger');
text(272,6,'End Trigger');
text(195,4,'Single Lane');
text(235,4,'Double Lane');
text(269,4,'Single Lane');

%% Metric: Yaw & derivatives
%last update MS 05.16.16

% Process the signals for Yaw
Processed_Yaw = SignalProcessor(TrialMap,Yaw,DirEvent,EndEvent);

%Congruent Right
[YawStruct.CongruentRight,YawStruct.CongruentRightSnipped,YawStruct.CongruentRightConcat] = Processed_Yaw{1,:};
YawStruct.CongruentRightSnippedMean = mean(YawStruct.CongruentRightConcat,1);

%Congruent Left
[YawStruct.CongruentLeft,YawStruct.CongruentLeftSnipped,YawStruct.CongruentLeftConcat] = Processed_Yaw{2,:};
YawStruct.CongruentLeftSnippedMean = mean(YawStruct.CongruentLeftConcat,1);

%Incongruent Right
[YawStruct.IncongruentRight,YawStruct.IncongruentRightSnipped,YawStruct.IncongruentRightConcat] = Processed_Yaw{3,:};
YawStruct.IncongruentRightSnippedMean = mean(YawStruct.IncongruentRightConcat,1);

%Incongruent Left
[YawStruct.IncongruentLeft,YawStruct.IncongruentLeftSnipped,YawStruct.IncongruentLeftConcat] = Processed_Yaw{4,:};
YawStruct.IncongruentLeftSnippedMean = mean(YawStruct.IncongruentLeftConcat,1);

%--------------------------%
% 1st derivative - yaw rate

% Process the signals for Yaw Rate
Processed_YawRate = SignalProcessor(TrialMap,YawRateStruct.matrix,DirEvent,EndEvent);

%Congruent Right
[YawRateStruct.CongruentRight,YawRateStruct.CongruentRightSnipped,YawRateStruct.CongruentRightConcat] = Processed_YawRate{1,:};
YawRateStruct.CongruentRightSnippedMean = mean(YawRateStruct.CongruentRightConcat,1);

%Congruent Left
[YawRateStruct.CongruentLeft,YawRateStruct.CongruentLeftSnipped,YawRateStruct.CongruentLeftConcat] = Processed_YawRate{2,:};
YawRateStruct.CongruentLeftSnippedMean = mean(YawRateStruct.CongruentLeftConcat,1);

%Incongruent Right
[YawRateStruct.IncongruentRight,YawRateStruct.IncongruentRightSnipped,YawRateStruct.IncongruentRightConcat] = Processed_YawRate{3,:};
YawRateStruct.IncongruentRightSnippedMean = mean(YawRateStruct.IncongruentRightConcat,1);

%Incongruent Left
[YawRateStruct.IncongruentLeft,YawRateStruct.IncongruentLeftSnipped,YawRateStruct.IncongruentLeftConcat] = Processed_YawRate{4,:};
YawRateStruct.IncongruentLeftSnippedMean = mean(YawRateStruct.IncongruentLeftConcat,1);


%--------------------------%
% 2nd derivative - yaw accel

% Process the signals for Yaw Accel
Processed_YawAccel = SignalProcessor(TrialMap,YawAccelStruct.matrix,DirEvent,EndEvent);

%Congruent Right
[YawAccelStruct.CongruentRight,YawAccelStruct.CongruentRightSnipped,YawAccelStruct.CongruentRightConcat] = Processed_YawAccel{1,:};
YawAccelStruct.CongruentRightSnippedMean = mean(YawAccelStruct.CongruentRightConcat,1);

%Congruent Left
[YawAccelStruct.CongruentLeft,YawAccelStruct.CongruentLeftSnipped,YawAccelStruct.CongruentLeftConcat] = Processed_YawAccel{2,:};
YawAccelStruct.CongruentLeftSnippedMean = mean(YawAccelStruct.CongruentLeftConcat,1);

%Incongruent Right
[YawAccelStruct.IncongruentRight,YawAccelStruct.IncongruentRightSnipped,YawAccelStruct.IncongruentRightConcat] = Processed_YawAccel{3,:};
YawAccelStruct.IncongruentRightSnippedMean = mean(YawAccelStruct.IncongruentRightConcat,1);

%Incongruent Left
[YawAccelStruct.IncongruentLeft,YawAccelStruct.IncongruentLeftSnipped,YawAccelStruct.IncongruentLeftConcat] = Processed_YawAccel{4,:};
YawAccelStruct.IncongruentLeftMeanYawAccel = mean(YawAccelStruct.IncongruentLeftConcat,1);


%% Metric: Handwheel position and derivatives
%last update MS 05.16.16

% Process the signals for HWPosition
Processed_HWPosition = SignalProcessor(TrialMap,HWPosition,DirEvent,EndEvent);

%CONGRUENT, RIGHT HWPOSITION
[HWPositionStruct.CongruentRight,HWPositionStruct.CongruentRightSnipped,HWPositionStruct.CongruentRightConcat] = Processed_HWPosition{1,:};
HWPositionStruct.CongruentRightSnippedMean = mean(HWPositionStruct.CongruentRightConcat,1);

%CONGRUENT, LEFT HWPOSITION
[HWPositionStruct.CongruentLeft,HWPositionStruct.CongruentLeftSnipped,HWPositionStruct.CongruentLeftConcat] = Processed_HWPosition{2,:};
HWPositionStruct.CongruentLeftSnippedMean = mean(HWPositionStruct.CongruentLeftConcat,1);

%INCONGRUENT, RIGHT HWPOSITION
[HWPositionStruct.IncongruentRight,HWPositionStruct.IncongruentRightSnipped,HWPositionStruct.IncongruentRightConcat] = Processed_HWPosition{3,:};
HWPositionStruct.IncongruentRightSnippedMean = mean(HWPositionStruct.IncongruentRightConcat,1);

%INCONGRUENT, LEFT HWPOSITION
[HWPositionStruct.IncongruentLeft,HWPositionStruct.IncongruentLeftSnipped,HWPositionStruct.IncongruentLeftConcat] = Processed_HWPosition{4,:};
HWPositionStruct.IncongruentLeftSnippedMean = mean(HWPositionStruct.IncongruentLeftConcat,1);

%PLOT CONGRUENT VS INCONGRUENT MEAN HW POSITION FOR LEFT TURN
%plot using shadederrorbar for congruent left
figure
subplot(2,1,1)
h_ConLeftHWPos = shadedErrorBar(PosXStruct.CongruentLeftSnippedMean,HWPositionStruct.CongruentLeftSnippedMean,...
    std(HWPositionStruct.CongruentLeftConcat,0,1),'b',1);
hold on
h_InconLeftHWPos = shadedErrorBar(PosXStruct.CongruentLeftSnippedMean,(-1*HWPositionStruct.IncongruentLeftSnippedMean),...
    std(HWPositionStruct.IncongruentLeftConcat,0,1),'r',1);

%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-6 6],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-6 6],'b--','LineWidth',1);

%add direction trigger
plot([meanDirEventPos meanDirEventPos],[-6 6],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-6 6],'r--','LineWidth',2);
hold off

title('Mean Handwheel Position | Left Turn','FontSize',FontSize);
%xlabel('Longitudinal Position (m)');
ylabel('Handwheel Position (deg)','FontSize',FontSize);
legend([h_ConLeftHWPos.mainLine,h_InconLeftHWPos.mainLine],'Congruent','Incongruent','Location','southeast');

subplot(2,1,2)
%plot cones
scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
hold on
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
%add direction and end triggers
plot([meanDirEventPos meanDirEventPos],[-15 15],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-15 15],'r--','LineWidth',2);
%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-15 15],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-15 15],'b--','LineWidth',1);
ylabel('Lateral Position (m)','FontSize',FontSize);
xlabel('Longitudinal Position (m)','FontSize',FontSize);
hold off
axis([180 300 -7 7]);
%add labels for direction and end triggers
text(182,6,'Direction Trigger');
text(271,6,'End Trigger');
text(195,4,'Single Lane');
text(235,4,'Double Lane');
text(269,4,'Single Lane');

%PLOT CONGRUENT VS INCONGRUENT MEAN HW POSITION FOR RIGHT TURN
%find mean path (snipped) and plot using shadederrorbar for congruent right
figure
subplot(2,1,1)
h_ConRightHWPos = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,HWPositionStruct.CongruentRightSnippedMean,...
    std(HWPositionStruct.CongruentRightConcat,0,1),'b',1);
hold on
h_InconRightHWPos = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,(-1*HWPositionStruct.IncongruentRightSnippedMean),...
    std(HWPositionStruct.IncongruentRightConcat,0,1),'r',1);

%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-6 6],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-6 6],'b--','LineWidth',1);

%add direction trigger
plot([meanDirEventPos meanDirEventPos],[-6 6],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-6 6],'r--','LineWidth',2);
hold off

title('Mean Handwheel Position | Right Turn','FontSize',FontSize);
%xlabel('Longitudinal Position (m)');
ylabel('Handwheel Position (deg)','FontSize',FontSize);
legend([h_ConRightHWPos.mainLine,h_InconRightHWPos.mainLine],'Congruent','Incongruent','Location','southeast');

subplot(2,1,2)
%plot cones
scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
hold on
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
%add direction and end triggers
plot([meanDirEventPos meanDirEventPos],[-15 15],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-15 15],'r--','LineWidth',2);
%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-15 15],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-15 15],'b--','LineWidth',1);
ylabel('Lateral Position (m)','FontSize',FontSize);
xlabel('Longitudinal Position (m)','FontSize',FontSize);
hold off
axis([180 300 -7 7]);

%add labels for direction and end triggers
text(182,6,'Direction Trigger');
text(271,6,'End Trigger');
text(195,4,'Single Lane');
text(235,4,'Double Lane');
text(269,4,'Single Lane');

%box plot - reversals & reversals per second
overallMeans = [reversals.IncongruentMeanNum reversals.CongruentMeanNum];
mytitle = 'Number of Handwheel Reversals';
metricLabel = 'Number of Reversals';
BoxPlotMetric(reversals.num,overallMeans,mytitle,metricLabel,congruence,[]);

overallMeans = [reversals.IncongruentMeanNumPerSec reversals.CongruentMeanNumPerSec];
mytitle = 'Number of Handwheel Reversals Per Second';
metricLabel = 'Number of Reversals Per Sec';
BoxPlotMetric(reversals.numPerSec,overallMeans,mytitle,metricLabel,congruence,[]);

%-------------------%
% HW RATE

% Process the signals for HWRate
Processed_HWRate = SignalProcessor(TrialMap,HWRateStruct.matrix,DirEvent,EndEvent);

%CONGRUENT,RIGHT HWRATE (d/dt)
[HWRateStruct.CongruentRight,HWRateStruct.CongruentRightSnipped,HWRateStruct.CongruentRightConcat] = Processed_HWRate{1,:};
HWRateStruct.CongruentRightSnippedMean = mean(HWRateStruct.CongruentRightConcat,1);

%CONGRUENT, LEFT HWRATE (d/dt)
[HWRateStruct.CongruentLeft,HWRateStruct.CongruentLeftSnipped,HWRateStruct.CongruentLeftConcat] = Processed_HWRate{2,:};
HWRateStruct.CongruentLeftSnippedMean = mean(HWRateStruct.CongruentLeftConcat,1);

%INCONGRUENT, RIGHT HWRATE (d/dt)
[HWRateStruct.IncongruentRight,HWRateStruct.IncongruentRightSnipped,HWRateStruct.IncongruentRightConcat] = Processed_HWRate{3,:};
HWRateStruct.IncongruentRightSnippedMean = mean(HWRateStruct.IncongruentRightConcat,1);

%INCONGRUENT, LEFT HWRATE (d/dt)
[HWRateStruct.IncongruentLeft,HWRateStruct.IncongruentLeftSnipped,HWRateStruct.IncongruentLeftConcat] = Processed_HWRate{4,:};
HWRateStruct.IncongruentLeftSnippedMean = mean(HWRateStruct.IncongruentLeftConcat,1);

%PLOT CONGRUENT VS INCONGRUENT MEAN HW RATE FOR LEFT TURN
%plot using shadederrorbar for congruent left
figure
subplot(2,1,1)
h_ConLeftHWRate = shadedErrorBar(PosXStruct.CongruentLeftSnippedMean,HWRateStruct.CongruentLeftSnippedMean,...
    std(HWRateStruct.CongruentLeftConcat,0,1),'b',1);
hold on
h_InconLeftHWRate = shadedErrorBar(PosXStruct.CongruentLeftSnippedMean,(-1*HWRateStruct.IncongruentLeftSnippedMean),...
    std(HWRateStruct.IncongruentLeftConcat,0,1),'r',1);

%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-15 15],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-15 15],'b--','LineWidth',1);

%add direction and end triggers
plot([meanDirEventPos meanDirEventPos],[-15 15],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-15 15],'r--','LineWidth',2);
hold off

title('Mean Handwheel Rate | Left Turn','FontSize',FontSize);
%xlabel('Longitudinal Position (m)');
ylabel('Handwheel Rate (deg/sec)','FontSize',FontSize);
legend([h_ConLeftHWRate.mainLine,h_InconLeftHWRate.mainLine],'Congruent','Incongruent','location','southeast');

subplot(2,1,2)
%plot cones
scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
hold on
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
%add direction and end triggers
plot([meanDirEventPos meanDirEventPos],[-15 15],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-15 15],'r--','LineWidth',2);
%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-15 15],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-15 15],'b--','LineWidth',1);
ylabel('Lateral Position (m)','FontSize',FontSize);
xlabel('Longitudinal Position (m)','FontSize',FontSize);
hold off
axis([180 300 -7 7]);

%add labels for direction and end triggers
text(182,6,'Direction Trigger');
text(271,6,'End Trigger');
text(195,4,'Single Lane');
text(235,4,'Double Lane');
text(269,4,'Single Lane');

%PLOT CONGRUENT VS INCONGRUENT MEAN HW RATE FOR RIGHT TURN
%find mean path (snipped) and plot using shadederrorbar for congruent right
figure
subplot(2,1,1)
h_ConRightHWRate = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,HWRateStruct.CongruentRightSnippedMean,...
    std(HWRateStruct.CongruentRightConcat,0,1),'b',1);
hold on
h_InconRightHWRate = shadedErrorBar(PosXStruct.CongruentRightSnippedMean,(-1*HWRateStruct.IncongruentRightSnippedMean),...
    std(HWRateStruct.IncongruentRightConcat,0,1),'r',1);

%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-15 15],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-15 15],'b--','LineWidth',1);

%add direction trigger
plot([meanDirEventPos meanDirEventPos],[-15 15],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-15 15],'r--','LineWidth',2);
hold off

title('Mean Handwheel Rate | Right Turn','FontSize',FontSize);
%xlabel('Longitudinal Position (m)');
ylabel('Handwheel Rate (deg/sec)','FontSize',FontSize);
legend([h_ConRightHWRate.mainLine,h_InconRightHWRate.mainLine],'Congruent',...
    'Incongruent','location','southeast');

subplot(2,1,2)
%plot cones
s = scatter(cones.X,cones.YRightLower,50,[1 .35 0],'filled');
hold on
scatter(cones.X,cones.YRightUpper,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftLower,50,[1 .35 0],'filled');
scatter(cones.X,cones.YLeftUpper,50,[1 .35 0],'filled');
%add direction and end triggers
plot([meanDirEventPos meanDirEventPos],[-15 15],'k--','LineWidth',2);
plot([meanEndEventPos meanEndEventPos],[-15 15],'r--','LineWidth',2);
%add longitudinal cone section markers
plot([cones.X(4) cones.X(4)],[-15 15],'b--','LineWidth',1);
plot([cones.X(7) cones.X(7)],[-15 15],'b--','LineWidth',1);
ylabel('Lateral Position (m)','FontSize',FontSize);
xlabel('Longitudinal Position (m)','FontSize',FontSize);
hold off
axis([180 300 -7 7]);

%add labels for direction and end triggers
text(182,6,'Direction Trigger');
text(271,6,'End Trigger');
text(195,4,'Single Lane');
text(235,4,'Double Lane');
text(269,4,'Single Lane');

%split metric based on group number
HWRateStruct.groupSplit = metricSplitter(HWRateStruct.mean,groupNum);
HWAccelStruct.groupSplit = metricSplitter(HWAccelStruct.mean,groupNum);
HWJerkStruct.groupSplit = metricSplitter(HWJerkStruct.mean,groupNum);

%find subject-wise section average of metrics
%for all outputs: first column is congruent, second column is incongruent
[HWRateStruct.subjectSection1,HWRateStruct.subjectSection2,HWRateStruct.subjectSection3] = ...
    ttestSectioner(HWRateStruct.matrix,[meanDoubleLaneStartIdx meanSingleLaneStartIdx],congruence);
[HWAccelStruct.subjectSection1,HWAccelStruct.subjectSection2,HWAccelStruct.subjectSection3] = ...
    ttestSectioner(HWAccelStruct.matrix,[meanDoubleLaneStartIdx meanSingleLaneStartIdx],congruence);
[HWJerkStruct.subjectSection1,HWJerkStruct.subjectSection2,HWJerkStruct.subjectSection3] = ...
    ttestSectioner(HWJerkStruct.matrix,[meanDoubleLaneStartIdx meanSingleLaneStartIdx],congruence);
[YawRateStruct.subjectSection1,YawRateStruct.subjectSection2,YawRateStruct.subjectSection3] = ...
    ttestSectioner(YawRateStruct.matrix,[meanDoubleLaneStartIdx meanSingleLaneStartIdx],congruence);

%% Metric: Time spent outside of cones
% last update ZS 06.13.16

disp('Calculating time out of bounds...');

outOfBounds.idx = cell(size(PosX));
outOfBounds.time = zeros(size(PosX));

%Congruent, right trials
%check if any points in the vehicle position are outside of the bounds
%get the indices for which this is true

for n = 1:size(TrialMap,1)
    for i = 1:length(TrialMap{n,1})
        x = TrialMap{n,1}(i);
        y = TrialMap{n,2}(i);
        if mod(n,2) == 0 %if this is a left turn
            absLowerBound = conePath.leftlower';
            absUpperBound = conePath.leftupper';
            relLowerBound = RelPathPack.bounds{3};
            relUpperBound = RelPathPack.bounds{4};
        else
            absLowerBound = conePath.rightlower';
            absUpperBound = conePath.rightupper';
            relLowerBound = RelPathPack.bounds{1};
            relUpperBound = RelPathPack.bounds{2};
        end
        if ~settings.choosePathDev %if the user chose absolute path error
            outOfBounds.idx{x,y} = ((Processed_PosY{n,2}{i} < absLowerBound) | ...
                (CongruentRightPosY_snipped{i} > absUpperBound));
            outOfBounds.time(x,y) = sum((outOfBounds.idx{x,y} == 1) * dt);
        else
            outOfBounds.idx{x,y} = ((Processed_PosY{n,2}{i} < relLowerBound{x}) | ...
                (Processed_PosY{n,2}{i} > relUpperBound{x}));
            outOfBounds.time(x,y) = sum((outOfBounds.idx{x,y} == 1) * dt);
        end
    end
end

[outOfBounds.IncongruentMean,outOfBounds.CongruentMean] = ...
    ttestFormat(outOfBounds.time,congruence,'congruence');

% for i = 1:size(outOfBounds.time,1)
%     crntOutOfBoundsTimeRow = outOfBounds.time(i,:);
%     crntOutOfBoundsTimeCongruent = crntOutOfBoundsTimeRow(congruence(i,:) == 1);
%     crntOutOfBoundsTimeIncongruent = crntOutOfBoundsTimeRow(congruence(i,:) == 0);
%     outOfBounds.CongruentMean(i,1) = mean(crntOutOfBoundsTimeCongruent,'omitnan');
%     outOfBounds.IncongruentMean(i,1) = mean(crntOutOfBoundsTimeIncongruent,'omitnan');
% end

%sum the times from each trial based on congruence
%group 1 includes trial less than trial no. 17
%group 2 includes trials between 17 and 32
%group 3 inlcudes all those greater than 32 (aka the repeats)
outOfBounds.timeCongruent = sum(outOfBounds.time(congruence == 1));
outOfBounds.timeCongruentGroup1 = sum(outOfBounds.time(congruence == 1 & groupNum == 1));
outOfBounds.timeCongruentGroup2 = sum(outOfBounds.time(congruence == 1 & (groupNum == 2 | groupNum == 3)));
outOfBounds.timeIncongruent = sum(outOfBounds.time(congruence == 0));
outOfBounds.timeIncongruentGroup1 = sum(outOfBounds.time(congruence == 0 & groupNum == 1));
outOfBounds.timeIncongruentGroup2 = sum(outOfBounds.time(congruence == 0 & (groupNum == 2 | groupNum == 3)));

disp(['Total time out of bounds, congruent: ',...
    num2str(outOfBounds.timeCongruent),' seconds.']);
disp(['Total time out of bounds, congruent, group 1: ',...
    num2str(outOfBounds.timeCongruentGroup1),' seconds.']);
disp(['Total time out of bounds, congruent, group 2: ',...
    num2str(outOfBounds.timeCongruentGroup2),' seconds.']);
disp(['Total time out of bounds, incongruent: ',...
    num2str(outOfBounds.timeIncongruent),' seconds.']);
disp(['Total time out of bounds, incongruent, group 1: ',...
    num2str(outOfBounds.timeIncongruentGroup1),' seconds.']);
disp(['Total time out of bounds, incongruent, group 2: ',...
    num2str(outOfBounds.timeIncongruentGroup2),' seconds.']);
disp('%%%%%%%%%%%%%%%%%%%%%');

%% Package signals into final matrix of structs
% last update ZS 06.21.16

for k = 1:size(TrialMap,1)
    %Populate congruent,right trials
    for i = 1:length(TrialMap{k,1})
        x = TrialMap{k,1}(i); %get current matrix x index
        y = TrialMap{k,2}(i); %get current matrix y index
        
        %Description
        participant{x,y}.UserID = UserID{x,y}(1); %first entry of vector is sufficient
        participant{x,y}.congruence = congruence(x,y);
        participant{x,y}.TurnDirection = TurnDirection(x,y);
        participant{x,y}.CorrectTurn = CorrectTurn(x,y);
        participant{x,y}.NumBadTrials = BadTrialsCount4Subject(x,1);
        participant{x,y}.NumRepeatTrials = Subjects(x).nRepeatTrials;
        participant{x,y}.isRepeat = repeatFlag(x,y);
        participant{x,y}.groupNum = groupNum(x,y);
        
        %Metrics
        participant{x,y}.HWPosition = Processed_HWPosition{k,1}{i};
        participant{x,y}.HWPosition_snipped = Processed_HWPosition{k,2}{i};
        participant{x,y}.NumHWReversals = reversals.num(x,y);
        participant{x,y}.ReversalsPerSec = reversals.numPerSec(x,y);
        participant{x,y}.HWRate = HWRateStruct.matrix{x,y};
        participant{x,y}.MeanHWRate = HWRateStruct.mean(x,y);
        participant{x,y}.RMSHWRate = HWRateStruct.rms(x,y);
        participant{x,y}.MaxHWRate = HWRateStruct.max(x,y);
        participant{x,y}.HWAccel = HWAccelStruct.matrix{x,y};
        participant{x,y}.MeanHWAccel = HWAccelStruct.mean(x,y);
        participant{x,y}.RMSHWAccel = HWAccelStruct.rms(x,y);
        participant{x,y}.MaxHWAccel = HWAccelStruct.max(x,y);
        participant{x,y}.HWJerk = HWJerkStruct.matrix{x,y};
        participant{x,y}.MeanHWJerk = HWJerkStruct.mean(x,y);
        participant{x,y}.RMSHWJerk = HWJerkStruct.rms(x,y);
        participant{x,y}.MaxHWJerk = HWJerkStruct.max(x,y);
        participant{x,y}.time2HWPeak = time2HWPeak.matrix(x,y);
        participant{x,y}.time2HWPeakCongruentMean = time2HWPeak.CongruentMean;
        participant{x,y}.time2HWPeakIncongruentMean = time2HWPeak.IncongruentMean;
        participant{x,y}.dist2HWPeak = dist2HWPeak.matrix(x,y);
        participant{x,y}.dist2HWPeakCongruentMean = dist2HWPeak.CongruentMean;
        participant{x,y}.dist2HWPeakIncongruentMean = dist2HWPeak.IncongruentMean;
        participant{x,y}.time2PosYPeak = time2PosYPeak.matrix(x,y);
        participant{x,y}.time2PosYPeakCongruentMean = time2PosYPeak.CongruentMean;
        participant{x,y}.time2PosYPeakIncongruentMean = time2PosYPeak.IncongruentMean;
        participant{x,y}.dist2PosYPeak = dist2PosYPeak.matrix(x,y);
        participant{x,y}.dist2PosYCongruentMean = dist2HWPeak.CongruentMean;
        participant{x,y}.dist2PosYIncongruentMean = dist2HWPeak.IncongruentMean;
        participant{x,y}.accelY = Processed_accelY{k,1}{i};
        participant{x,y}.accelY_snipped = Processed_accelY{k,2}{i};
        participant{x,y}.Vy = Processed_Vy{k,1}{i};
        participant{x,y}.Vy_snipped = Processed_Vy{k,2}{i};
        participant{x,y}.Yaw = Processed_YawRate{k,1}{i};
        participant{x,y}.Yaw_snipped = Processed_YawRate{k,2}{i};
        participant{x,y}.latError = latError.matrix{x,y};
        participant{x,y}.latErrorCongruentMean = latError.CongruentMean;
        participant{x,y}.latErrorIncongruentMean = latError.IncongruentMean;
        participant{x,y}.longError = longError.matrix{x,y};
        participant{x,y}.longErrorCongruentMean = longError.CongruentMean;
        participant{x,y}.longErrorIncongruentMean = longError.IncongruentMean;
        participant{x,y}.outOfBounds = outOfBounds.time(x,y);
        participant{x,y}.outOfBoundsTimeSumCongruent = outOfBounds.timeCongruent;
        participant{x,y}.outOfBoundsTimeSumIncongruent = outOfBounds.timeIncongruent;
    end
end

%further packaging
result.metrics = participant;
result.latErrorStruct = latError;
result.longErrorStruct = longError;
result.outOfBoundsStruct = outOfBounds;
result.reversalsStruct = reversals;
result.time2HWPeakStruct = time2HWPeak;
result.dist2HWPeakStruct = dist2HWPeak;
result.time2PosYStruct = time2PosYPeak;
result.dist2PosYStruct = dist2PosYPeak;
result.HWRateStruct = HWRateStruct;
result.HWAccelStruct = HWAccelStruct;
result.YawRateStruct = YawRateStruct;
result.YawAccelStruct = YawAccelStruct;

%concatenate congruent and incongruent means for stats
result.ttestMat_congruent = [reversals.CongruentMeanNum,reversals.CongruentMeanNumPerSec,HWRateStruct.CongruentMean,HWAccelStruct.CongruentMean,HWJerkStruct.CongruentMean,...
    VyStruct.CongruentMean,accelYStruct.CongruentMean,time2HWPeak.CongruentMean,time2PosYPeak.CongruentMean,dist2HWPeak.CongruentMean,dist2PosYPeak.CongruentMean,...
    YawRateStruct.CongruentMean,outOfBounds.CongruentMean,latError.CongruentMean,HWRateStruct.CongruentRMS,HWAccelStruct.CongruentRMS,HWJerkStruct.CongruentRMS,...
    VyStruct.CongruentRMS,accelYStruct.CongruentRMS,YawRateStruct.CongruentRMS,HWRateStruct.subjectSection1.mean(:,1),HWRateStruct.subjectSection2.mean(:,1),...
    HWRateStruct.subjectSection3.mean(:,1),HWAccelStruct.subjectSection1.mean(:,1),HWAccelStruct.subjectSection2.mean(:,1),HWAccelStruct.subjectSection3.mean(:,1),...
    HWJerkStruct.subjectSection1.mean(:,1),HWJerkStruct.subjectSection2.mean(:,1),HWJerkStruct.subjectSection3.mean(:,1),YawRateStruct.subjectSection1.mean(:,1),...
    YawRateStruct.subjectSection2.mean(:,1),YawRateStruct.subjectSection3.mean(:,1)];

result.ttestMat_incongruent = [reversals.IncongruentMeanNum,reversals.IncongruentMeanNumPerSec,HWRateStruct.IncongruentMean,HWAccelStruct.IncongruentMean,HWJerkStruct.IncongruentMean,...
    VyStruct.IncongruentMean,accelYStruct.IncongruentMean,time2HWPeak.IncongruentMean,time2PosYPeak.IncongruentMean,dist2HWPeak.IncongruentMean,dist2PosYPeak.IncongruentMean,...
    YawRateStruct.IncongruentMean,outOfBounds.IncongruentMean,latError.IncongruentMean,HWRateStruct.IncongruentRMS,HWAccelStruct.IncongruentRMS,HWJerkStruct.IncongruentRMS,...
    VyStruct.IncongruentRMS,accelYStruct.IncongruentRMS,YawRateStruct.IncongruentRMS,HWRateStruct.subjectSection1.mean(:,2),HWRateStruct.subjectSection2.mean(:,2),...
    HWRateStruct.subjectSection3.mean(:,2),HWAccelStruct.subjectSection1.mean(:,2),HWAccelStruct.subjectSection2.mean(:,2),HWAccelStruct.subjectSection3.mean(:,2),...
    HWJerkStruct.subjectSection1.mean(:,2),HWJerkStruct.subjectSection2.mean(:,2),HWJerkStruct.subjectSection3.mean(:,2),YawRateStruct.subjectSection1.mean(:,2),...
    YawRateStruct.subjectSection2.mean(:,2),YawRateStruct.subjectSection3.mean(:,2)];

result.ttestMat = [result.ttestMat_congruent,result.ttestMat_incongruent];
%% Call stats function
%last update ZS 05.23.16

%we'll reuse the stats function from profx

% congruenceGroup.name = 'Congruence';
% congruenceGroup.var = congruence;
% turnDirGroup.name = 'TurnDirection';
% turnDirGroup.var = TurnDirection;

%result.stats = profx_stats(result,congruenceGroup);
result = profx_stats(result);

%HW Reversals
% disp(' ');
% disp('Statistical analysis for HW Reversals wrt steering mode');
% stats.reversals.steerMode = profx_stats(reversals.num,congruenceGroup);
% % disp('Statistical analysis for HW Reversals wrt turn direction');
% % stats.reversals.steerMode = profx_stats(reversals.num,turnDirGroup);
%
% %HW Reversals Rate
% disp(' ');
% disp('Statistical analysis for HWRate wrt steering mode');
% stats.hwrate.steerMode = profx_stats(HWRateStruct.mean,congruenceGroup);
% % disp('Statistical analysis for HWRate wrt turn direction');
% % stats.hwrate.steerMode = profx_stats(HWRateStruct.mean,turnDirGroup);
%
% % stats.reversals.turnDir = profx_stats(reversals.num,turnDirGroup);
% % stats.steerMode.reversalsPerSec = profx_stats(reversals.numPerSec,congruenceGroup);
% % stats.turnDir.reversalsPerSec = profx_stats(reversals.numPerSec,turnDirGroup);

%% Save participant.mat
%last update ZS 07.27.16
if settings.saveGroupMetrics
    cd(saveDir);
    saveName = ['MagnetoMetrics_',date];
    save(saveName,'result');
    disp(['Successfully saved ''', saveName,'''.mat to ',saveDir]);
end
