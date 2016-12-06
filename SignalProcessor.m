% Function: SignalProcessor
% -------------------------
% Executes in order to extract different trial types from vehicle signals
function [result] = SignalProcessor(TrialMap,Signal,DirEvent,EndEvent)

snipLength = 400;

mapLength = length(TrialMap(:,1));  % Determines the number of trials

result = cell(mapLength, 3);    % Preallocates result cell

% A loop that runs through all the trials
for i = 1 : mapLength
    
    TMap = TrialMap(i,:);   % Creates a temporary dummy trial map for each trial
    
    % Preallocation of processed signals
    indexed = cell(length(TMap{1}),1); %trials matching the trial type specified by trialmap are pulled out of the matrix
    snipped = cell(length(TMap{1}),1); %the indexed signal is snipped to just the turn event
    snipped_concat = []; %the indexed and snipped signals are concatenated in order to use shaded error bar

    for k = 1:length(TMap{1})
        if (DirEvent.Idx(TMap{1}(k),TMap{2}(k)) ~= 0 && ...
                EndEvent.Idx(TMap{1}(k),TMap{2}(k)) ~= 0)
            indexed{k,1} = Signal{TMap{1}(k),TMap{2}(k)};
            snipped{k,1} = indexed{k,1}(DirEvent.Idx(TMap{1}(k),TMap{2}(k))-75:...
                EndEvent.Idx(TMap{1}(k),TMap{2}(k))+50);
            %snipped{k,1} = snipped{k,1}(1:snipLength)'; %transpose for concatenation
            
%             xSpacing = linspace(DirEvent.Pos(TMap{1}(k),TMap{2}(k)),...
%                 EndEvent.Pos(TMap{1}(k),TMap{2}(k)),snipLength);
            
            %resample data to make all vectors the same length
            %disp(['Sampling ratio: ',num2str(snipLength/length(snipped{k,1}))]);
            snipped{k,1} = resample(snipped{k,1},snipLength,length(snipped{k,1}),0);
            
%             snipped{k,1} = ...
%                 interp1(PosXUnique.pos{TMap{1}(k),TMap{2}(k)}(DirEvent.Idx(TMap{1}(k),...
%                 TMap{2}(k)):EndEvent.Idx(TMap{1}(k),TMap{2}(k))),...
%                 snipped{k,1},xSpacing);


            snipped_concat = [snipped_concat;snipped{k,1}(1:snipLength)'];
            %snipped{k,1} = snipped{k,1}(1:snipLength)'; %transpose back
        end
    end
    
    % Stores the processed signal for each trial in the preallocated results cell
    result{i,1} = indexed;
    result{i,2} = snipped;
    result{i,3} = snipped_concat;
end

end

