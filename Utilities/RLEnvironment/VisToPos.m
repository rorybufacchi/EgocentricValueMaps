function [handR,handC,goalR,goalC,thrR,thrC] =VisToPos(inpMazes,handVisVal,goalVisVal,thrVisVal)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% $$$ NOT PERFRCT YET, because of occasional overlap!!

for kTrial=1:size(inpMazes,3)
    if rem(kTrial,50000)==0
        disp(kTrial)
    end
    if kTrial==1503
        disp(kTrial)
    end
    try
        try
            [handR(kTrial), handC(kTrial)]=find(inpMazes(:,:,kTrial)==handVisVal);
        catch
            % If can't find the hand prosition, it's probably obscured by the goal
            [handR(kTrial), handC(kTrial)]=find(inpMazes(:,:,kTrial)==goalVisVal);
        end
        try
            [goalR(kTrial), goalC(kTrial)]=find(inpMazes(:,:,kTrial)==goalVisVal);
        catch
            [goalR(kTrial), goalC(kTrial)]=find(inpMazes(:,:,kTrial)==handVisVal);
        end
        if thrVisVal~=50
            [thrR(kTrial), thrC(kTrial)]=find(inpMazes(:,:,kTrial)==thrVisVal);
        end
    catch
        % If through some idiocy on my part, there is no stimulus (i.e.
        % I overwrote all -ve threats with 50 and stuff), then just set
        % the hand in the middle
        if ~sum(sum(inpMazes(:,:,kTrial)~=50 & inpMazes(:,:,kTrial)~=0))
            handR(kTrial)=size(inpMazes,1)-2;
            handC(kTrial)=round(size(inpMazes,2)/2);
            goalR(kTrial)=size(inpMazes,1)-2;
            goalC(kTrial)=round(size(inpMazes,2)/2);
        else
            % $$$ THIS IS TERRIBLE AND I SHOULD DO SOMETHING ABOUT IT
            % (maybe)
% %             error
        end
    end
end

end
