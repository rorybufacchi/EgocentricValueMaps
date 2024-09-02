function [] = SquareSymAxes(h)
%[] = SymAxes(f)
%   Sets the x- and y- axes symetrically and equally
%   

if ~exist('h','var')
    h = gca;
end

xLims = xlim;
yLims = ylim;

maxV = max([abs(xLims) abs(yLims)]);

xlim([-max(abs(maxV)) max(abs(maxV))]);
ylim([-max(abs(maxV)) max(abs(maxV))]);
zlim([-max(abs(maxV)) max(abs(maxV))]);


end

