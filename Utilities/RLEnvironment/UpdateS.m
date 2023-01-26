function [S] = UpdateS(s, w, tempMaze2D)
%UPDATES Takes current system configuration and updates the state S
%depending on the flags that define the state
%   Detailed explanation placeholder

if s.fl.vis ==1
    if s.fl.vpr==1
        S=[w.lmb.col tempMaze2D(:)'];
    else
        S=tempMaze2D(:)';
    end
elseif s.fl.thr==1
    S=[w.lmb.col w.goal.row w.goal.col w.thr.row w.thr.col];
else
    S=[w.lmb.col w.goal.row w.goal.col];
end

if s.fl.tch==1
    S=[S w.touchV(:)'];
end

if s.fl.bdyMov==1
    S=[S w.bdy.col]; 
end

if s.fl.ToolChange==1
    S=[S w.lmb.ToolPresent];
end

if s.fl.rtTask==1
    S=[S w.rtTask.percieveTouch];
end

if s.fl.extraInput==1
    S = [S w.xtra.currV];
end


end

