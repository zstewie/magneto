function [ crntShift,halves ] = findLag(congruentMeanPosX,congruentMeanPosY,crntTrialPosX,crntTrialPosY)
%FINDLAG Function that computes longitudinal error between current trial
%and congruent mean
%   Inputs: Congruent mean longitudinal position, congruent mean lateral
%   position, current trial longitudinal position, current trial lateral
%   position.
%   Positive lag means that the current trial had a later peak time than the
%   congruent mean. negative means that the current trial had an earlier peak
%   time than the congruent mean.

crntShift = NaN(size(congruentMeanPosX));

%finding lateral peaks and subtracting longitudinal positions
[~,conMeanPeakLoc] = findpeaks(abs(congruentMeanPosY),'MinPeakHeight',...
    0.99 * max(abs(congruentMeanPosY)));
[~,crntPeakLoc] = findpeaks(abs(crntTrialPosY),'MinPeakHeight',...
    0.99 * max(abs(crntTrialPosY)));

%if no peak is found, return empty vars (empty return is handled outside of function)
if isempty(crntPeakLoc)
    crntShift = [];
    halves = {};
    return
else
    conMeanPeakLoc = conMeanPeakLoc(1); 
    crntPeakLoc = crntPeakLoc(1);
end

%grab offset from beginning and end of passed vectors
latPosBaseline_begin = crntTrialPosY(1);
latPosMeanBaseline_begin = congruentMeanPosY(1);
latPosBaseline_end = crntTrialPosY(end);
latPosMeanBaseline_end = congruentMeanPosY(end);

%fill out first half structure - offset to beginning value
firstHalf.meanPosX = congruentMeanPosX(1:conMeanPeakLoc);
firstHalf.crntPosX = crntTrialPosX(1:crntPeakLoc);
firstHalf.meanPosY = congruentMeanPosY(1:conMeanPeakLoc) - latPosMeanBaseline_begin;
firstHalf.crntPosY = crntTrialPosY(1:crntPeakLoc) - latPosBaseline_begin;
firstHalf.peakPercent_crnt = (abs(firstHalf.crntPosY)) / ...
    (max(abs(firstHalf.crntPosY))) * 100;
firstHalf.peakPercent_mean = (abs(firstHalf.meanPosY)) / ...
    (max(abs(firstHalf.meanPosY))) * 100;

%fill out second half structure - offset to end value
secondHalf.meanPosX = congruentMeanPosX(conMeanPeakLoc+1:length(congruentMeanPosX));
secondHalf.crntPosX = crntTrialPosX(crntPeakLoc+1:length(congruentMeanPosX));
secondHalf.meanPosY = congruentMeanPosY(conMeanPeakLoc+1:length(congruentMeanPosX))...
    - latPosMeanBaseline_end;
secondHalf.crntPosY = crntTrialPosY(crntPeakLoc+1:length(congruentMeanPosX))...
    - latPosBaseline_end;
secondHalf.peakPercent_crnt = (abs(secondHalf.crntPosY)) / ...
    (max(abs(secondHalf.crntPosY))) * 100;
secondHalf.peakPercent_mean = (abs(secondHalf.meanPosY)) / ...
    (max(abs(secondHalf.meanPosY))) * 100;

%iterate over first half of vector - interpolate for mean position based on
%peak percent
for i = 1:length(firstHalf.crntPosX)
    crntPercent = firstHalf.peakPercent_crnt(i);
    crntPosition = firstHalf.crntPosX(i);
    crntMeanPosition = interp1(firstHalf.peakPercent_mean,...
        firstHalf.meanPosX,crntPercent);
        if(isnan(crntMeanPosition))
            crntPercent
            disp('Got a nan - 1st half');
        end
    crntShift(i) = crntPosition - crntMeanPosition;
end

%iterate over second half of vector - interpolate for mean position based on
%peak percent
for i = 1:length(secondHalf.crntPosX)
    crntPercent = secondHalf.peakPercent_crnt(i);
    crntPosition = secondHalf.crntPosX(i);
    crntMeanPosition = interp1(secondHalf.peakPercent_mean,...V
        secondHalf.meanPosX,crntPercent);
        if(isnan(crntMeanPosition))
            crntPercent
            disp('Got a nan - 2nd half');
        end
    crntShift(i + crntPeakLoc - 1) = crntPosition - crntMeanPosition;
end

%package structures into halves for use outside of this function
halves = {firstHalf;secondHalf};

end