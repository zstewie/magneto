function [ steerState ] = verifySteerState(steerInput,PosY,DirIdx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

steerState = 1; %default to congruent
posY_dirTrigger = PosY(DirIdx);

%if signs agree, congruent steering
if(steerInput(DirIdx) < 0 && PosY(DirIdx + 50) > posY_dirTrigger)
    steerState = 1;
elseif (steerInput(DirIdx) > 0 && PosY(DirIdx + 50) < posY_dirTrigger)
    steerState = 1;
else
    steerState = 0;
end

end

