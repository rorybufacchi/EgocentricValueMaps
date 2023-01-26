function [newR,newC,resetFlag, varargout] = MoveObjDownRand(oldR,oldC,speedR,speedC, varargin)
%MoveObjDownRand Moves any RL pixel downward by one pixel, and possibly
% left or right randomly by one pixel
% oldR is the row where the object starts
% oldC is the columns where the object starts
% speedR is the default speed of the object. downwards is +ve and upwards
% is -ve. Default of all the old versions was 1.
% speedC is the default speed of the object. right is +ve and left
% is -ve. Default of all the old versions was 1.
% ResetFlag indicates whether the object has moved back to the top

global world2D

resetFlag=0;

if numel(varargin)>2
    randomSpread=varargin{3};
else
    randomSpread=2;
end

if oldR > 1 & oldR < size(world2D,1)-1
    %             oldR=oldR+round(rand(1,1))+1;
    newR=oldR+speedR;
    if newR >=size(world2D,1)-1
        newR = size(world2D,1)-1;
    end
    oldCTemp=oldC+round(randomSpread*(rand(1,1)-0.5)) + speedC; % $$$ UNDERSTAND THIS RANDOM MOTION AND PARAMETRISE IT
    while (oldCTemp < 2) | (oldCTemp > size(world2D,2)-1)
        oldCTemp=oldC+round(2*(rand(1,1)-0.5)) + speedC;
        if speedC~=0
            speedC=speedC-sign(speedC);
        end
    end
    newC=oldCTemp;
else
    [newR newC resetFlag]=ResetObjPos();
end

if numel(varargin)>0
    bodySize=varargin{1};
    bodyPos=varargin{2};
    % Store whether and where the body was hit
    bodyTouch=zeros([1 bodySize]);
% %     % This is if bodyPos is the center
% %     if newR==size(world2D,1)-1 & newC >= bodyPos -bodySize/2 & newC <= bodyPos + bodySize/2
% %         bodyTouch(newC-bodyPos + ceil(bodySize/2) )=1;
% %     end  
    if newR==size(world2D,1)-1 & newC >= bodyPos & newC < bodyPos + bodySize
        bodyTouch(newC-bodyPos + 1 )=1;
    end
    varargout{1}=bodyTouch;
    
end