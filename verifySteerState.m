function [ steerState ] = verifySteerState(steerPos,PosY,DirIdx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

steerState = 1; %default to congruent
initPosY = PosY(DirIdx);
delta = 100;

%if signs agree, congruent steering
if(steerPos(DirIdx) > 0)
    if(PosY(DirIdx+delta) > initPosY)
        steerState = 1;
    else
        steerState = 0;
    end
elseif(steerPos(DirIdx) < 0)
    if(PosY(DirIdx+delta) < initPosY)
        steerState = 1;
    else
        steerState = 0;
    end
end

end

