function [ steerState ] = verifySteerState(steerInput,PosY)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

steerState = 1; %default to congruent

[i,j] = find(steerInput > 5,1);

%if signs agree, congruent steering
if((steerInput(i,j) * PosY(i,j)) > 0)
    steerState = 1;
else
    steerState = 0;
end

end

