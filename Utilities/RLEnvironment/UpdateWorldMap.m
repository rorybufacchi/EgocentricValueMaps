function [] = UpdateWorldMap(w,prvR,prvC,newR,newC,resetFlag,objVal)

global tempworld2D;

% Find the old object position
oldObjPos=find(tempworld2D==objVal);

% reset the old object positions
tempworld2D(oldObjPos)=50;
% Plot all the pixels the threat has passed in this timestep
thrPix=FindPixels(prvR,prvC,newR,newC);
% If the threat is being reset, don't plot a line though
if resetFlag==1 | w.teleportLR ~=0
    tempworld2D(thrPix(1,end),thrPix(2,end))=objVal;
    if w.teleportLR ==-1
        tempworld2D(thrPix(1,1),thrPix(2,1):end-1)=objVal;
        tempworld2D(thrPix(1,end),2:thrPix(2,end))=objVal;
    elseif w.teleportLR ==1
        tempworld2D(thrPix(1,1),2:thrPix(2,1))=objVal;
        tempworld2D(thrPix(1,end),thrPix(2,end):end-1)=objVal;
    end
else
    % If only 1 point, the object is in the same position
    if size(thrPix,2)==1
        tempworld2D(oldObjPos)=objVal;
    else
    % Otherwise, don't use the first point because that is the old point
    for kPix=2:size(thrPix,2)
        tempworld2D(thrPix(1,kPix),thrPix(2,kPix))=objVal;
    end
    end
end

% Draw the borders of the world again
tempworld2D(1,:)=0;tempworld2D(end,:)=0; tempworld2D(:,1)=0; tempworld2D(:,end)=0;