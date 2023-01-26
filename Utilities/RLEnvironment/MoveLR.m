% status = 1 then move ahead successful
% status = 2 then bump into wall or boundary
% status = 3 then goal achieved
% status = 4 then bump into threat
% Move the robot to the next location if no bump
function [row,col,status] = MoveLR(s,w,row,col,action)
global tempworld2D;

% in case there the other body-part moves, use the appropriate move type
% if s.fl.bdyMov==1
if w.mvBdy==1
    % find the underscore Position first to split the name
    usPs = find(s.act.Name{action}=='_');
    actName = s.act.Name{action}(usPs+4:end);
elseif s.fl.bdyMov==1
    usPs = find(s.act.Name{action}=='_');
    actName = s.act.Name{action}(4:usPs-1);
else
    actName = s.act.Name{action};
end

% based on the current direction check whether next location is space or
% bump and get information of use below
[val,valid,w] = LookLR(s,w,row,col,actName);

% check if next location for moving is space
% other wise set the status
% this checks the collision with boundary of world
if valid == 1
    % now check if the next location for space or bump
    % this is for walls inside the world
    if val > 0
        oldRow = row; oldCol = col;
        switch actName
            case {'LEFT','LEFTBUTTON'}
                col = col - 1;
            case {'STAY','BUTTON'}
                col = col;
            case {'RIGHT','RIGHTBUTTON'}
                col = col + 1;
        end
        status = 1;
        
        if w.limbTeleLR==1
            col=size(w.world2D,2)-1;
        elseif w.limbTeleLR==-1
            col=2;
        end
        
        if val == 100
            % goal achieved
            status = 3;
        end
        
        if val ==10
            % Bumped into the threat
            status=4;
        end
        
        % update the current position of the robot in world for display
        tempworld2D(oldRow,oldCol) = 50;
        tempworld2D(row,col) = 60;
    elseif val == 0
        % bump into wall
        status = 2;
    end
else

    % return a bump signal if valid is 0
    status = 2;

end
