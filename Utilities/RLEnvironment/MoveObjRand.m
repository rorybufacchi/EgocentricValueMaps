function [newX,newY] = MoveObjRand(oldX,oldY)
%MoveObjRand Moves any RL pixel randomly by 1 random x and 1 random y. Does
%not go outside the maze

global maze2D

oldXTemp=oldX+sign((rand(1,1)-0.5));
oldYTemp=oldY+sign((rand(1,1)-0.5));
while (oldXTemp <2) | (oldXTemp > size(maze2D,1)-2)
    oldXTemp=oldX+sign((rand(1,1)-0.5));
end
while (oldYTemp <2) | (oldYTemp > size(maze2D,2)-2)
    oldYTemp=oldY+sign((rand(1,1)-0.5));
end
newX=oldXTemp;
newY=oldYTemp;


end

