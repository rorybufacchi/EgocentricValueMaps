function [rc] = FindPixels(oldR,oldC,newR,newC)
% Outputs the pixels which lie along an interpolated line (in pixel space)

nPoints=200;

rowVec=round(linspace(oldR,newR,nPoints));
colVec=round(linspace(oldC,newC,nPoints));

rc=[rowVec ; colVec];

% all points which are the same as the previous point
rmP=[0 sum(abs(diff(rc,[],2)))==0];

rc(:,find(rmP))=[];