function [newPos,s] = TrackDynamics(pos,s)
% Simulates the dynamics of the track spider and butterfly experiment

pos = pos(:);

% Define the points at which the butterflies and spiders can appear,
% relative to the hand
allStimPos = [  0  0 0 ; ... L7, collision course: HAND IS HIT!
               -2  2 0 ; ... L6, collision course: 10 cm from hand in x and y
               -5  2 0 ; ... L5, collision course
               -8  2 0 ; ... L4, collision course
                0 -4 0 ; ... L7, non-coll course: HAND IS MISSED!
               -2 -2 0 ; ... L6, non-coll course: 10 cm from hand in x and y
               -5 -2 0 ; ... L5, non-coll course
               -8 -2 0 ; ... L4, non-collision course
              -11  0 0 ; ... L3
              -14  0 0 ; ... L2
              -17  0 0 ];  % L1

allStimPos = flipud(allStimPos);

allStimPos = allStimPos + s.clc.nearPos';

% If the current position is not on the trajectory, just return the current
% position
if ~any(all( pos == allStimPos' ))
    newPos = pos;

% But if it IS on the trajectory, update the trajectory accordingly
else
    posOnTraj = find(all( pos == allStimPos' ));

    % If on the deterministic parts, just update to the next position
    if posOnTraj ~= 3 
        newPos = allStimPos(posOnTraj + 1 , :);
    else
        % Otherwise update to the collision track, BUT I HAVE TO SET THE
        % RANDOM SPREAD TO ALLOW IT TO GO ONTO THE NON-COLLISION TRACK
        % AS WELL
        newPos = allStimPos(8 , :);
    end
    
end

s.clc.specialTraject = allStimPos;


end