function [newR,newC,resetFlag, varargout] = MoveObjDownRand2(s, oldR,oldC,speedR,speedC, randomSpread)
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

% % %
% % % if strcmp(s.wrld.resetType,'Top')
% % % elseif strcmp(s.wrld.resetType,'Edges')
% % %
% % %     error('Incorrect reset of stimuli. Use one of the proper words')
% % % end

% % if oldR > 1 & oldR < size(world2D,1)-1
% %     %             oldR=oldR+round(rand(1,1))+1;
% %     newR=oldR+speedR;
% %     if newR >=size(world2D,1)-1
% %         newR = size(world2D,1)-1;
% %     end
% %     oldCTemp=oldC+round(randomSpread*(rand(1,1)-0.5)) + speedC; % $$$ UNDERSTAND THIS RANDOM MOTION AND PARAMETRISE IT
% %     while (oldCTemp < 2) | (oldCTemp > size(world2D,2)-1)
% %         oldCTemp=oldC+round(2*(rand(1,1)-0.5)) + speedC;
% %         if speedC~=0
% %             speedC=speedC-sign(speedC);
% %         end
% %     end
% %     newC=oldCTemp;
% % else
% %     [newR newC resetFlag]=ResetObjPos(s);
% % end

if strcmp(s.wrld.resetType,'BottomToTop')
	resetCondition = oldR > 2 & oldR < size(world2D,1)-1;
elseif strcmp(s.wrld.resetType,'Edges')
    resetCondition = (oldR > 2 & oldR < size(world2D,1)-1) & (oldC > 2 & oldC < size(world2D,2)-1);
end

if resetCondition==1
    %             oldR=oldR+round(rand(1,1))+1;
    newR=oldR+speedR + round(randomSpread*(rand(1,1)-0.5));
    if newR >=size(world2D,1)-1
        newR = size(world2D,1)-1;
    end
    if newR <=2
        newR = 2;
    end
    
    oldCTemp=oldC+round(randomSpread*(rand(1,1)-0.5)) + speedC; % $$$ UNDERSTAND THIS RANDOM MOTION AND PARAMETRISE IT
    newC=oldCTemp;
    if strcmp(s.wrld.resetType,'BottomToTop')
        while (oldCTemp < 2) | (oldCTemp > size(world2D,2)-1)
            oldCTemp=oldC+speedC + round(randomSpread*(rand(1,1)-0.5)) ;
            if speedC~=0
                speedC=speedC-sign(speedC);
            end
        end
        newC=oldCTemp;
    elseif strcmp(s.wrld.resetType,'Edges')
        
        if newC >=size(world2D,2)-1
            newC = size(world2D,2)-1;
        end
        if newC <=2
            newC = 2;
        end
        
    else
        error('Incorrect reset of stimuli. Use one of the proper words')
    end

else
    [newR newC resetFlag]=ResetObjPos(s);
end


if s.bdy.size>0
    % Store whether and where the body was hit
    bodyTouch=zeros([1 s.bdy.size]);
    % %     % This is if s.bdy.pos is the center
    % %     if newR==size(world2D,1)-1 & newC >= s.bdy.pos -s.bdy.size/2 & newC <= s.bdy.pos + s.bdy.size/2
    % %         bodyTouch(newC-s.bdy.pos + ceil(s.bdy.size/2) )=1;
    % %     end
    if newR==size(world2D,1)-1 & newC >= s.bdy.pos & newC < s.bdy.pos + s.bdy.size
        bodyTouch(newC-s.bdy.pos + 1 )=1;
    end
    varargout{1}=bodyTouch;
else
    varargout{1}=0;
end