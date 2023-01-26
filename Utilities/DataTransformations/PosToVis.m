function [outputMazes] = PosToVis(baseMaze,handR,handC,goalR,goalC,thrR,thrC,handVisVal,goalVisVal,thrVisVal,thrFlag)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

outputMazes=repmat(baseMaze,[1 1 length(handR(:))] );

% In case the plotter accidentally gave negative or zero values, set them to 1
goalR(goalR<1)=1;
goalC(goalC<1)=1;
thrR(thrR<1)=1;
thrC(thrC<1)=1;

for kTrial=1:numel(handR)
    outputMazes(handR(kTrial),handC(kTrial),kTrial)=handVisVal;
    outputMazes(goalR(kTrial),goalC(kTrial),kTrial)=goalVisVal;
    if thrFlag==1
        outputMazes(thrR(kTrial),thrC(kTrial),kTrial)=thrVisVal;
    end
end

end

