function [] = SquareAxes(h)
%[] = SymAxes(f)
%   Sets the x- and y- axes to be of equal size, but centered around the
%   data
%   

if ~exist('h','var')
    h = gca;
end

xLims = xlim;
yLims = ylim;
xMean = mean(xLims);
yMean = mean(yLims);

maxDiff = max([abs(diff(xLims)) abs(diff(yLims))]);

xlim(xMean + [-maxDiff, maxDiff]./2);
ylim(yMean + [-maxDiff, maxDiff]./2);


end

