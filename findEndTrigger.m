function [ idx ] = findEndTrigger( PosX )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%value from fileconverter script (which was taken from simulator)
endTriggerLoc = 280.7;
idx = find(PosX > endTriggerLoc,1);

end

