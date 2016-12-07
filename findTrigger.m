function [ idx ] = findTrigger( PosX,type )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if strcmp(type,'end')
    %value from fileconverter script (which was taken from simulator)
    endTriggerLoc = 280.7;
    idx = find(PosX > endTriggerLoc,1);
elseif strcmp(type,'dir')
    dirTriggerLoc = 180.65;
    idx = find(PosX > dirTriggerLoc,1);
end

end

