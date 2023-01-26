function [newR,newC,resetFlag, varargout] = ResetObjPos(s,varargin)
%ResetObjPos
% The variable inputs and outputs are for compatibility with older
% versions. They reset the speed of the object

global world2D


% $$ MAYBE I NEED A SETTING FOR THE STANDARD OBJECT RESET POSITION
% % if strcmp(s.wrld.resetType,'BottomToTop')

% % % % newR = 1+ceil(3*rand(1,1)) ;
% % % newR = 2 ;
% % % newC = randsample([2:(size(world2D,2)-1)],1);

newR = randsample([s.wrld.resetPos(1,1):s.wrld.resetPos(1,end)],1);
newC = randsample([s.wrld.resetPos(2,1):s.wrld.resetPos(2,end)],1);

% % elseif strcmp(s.wrld.resetType,'Edges')
% %     
% % end
    
resetFlag=1;


if length(varargin)>0
    allSpR=varargin{1};
    allSpC=varargin{2};
end

if length(varargin)>0
    varargout{1}=randsample(allSpR,1);
    varargout{2}=randsample(allSpC,1);
end

end