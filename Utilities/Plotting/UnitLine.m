function [] = UnitLine(h)
%[] = UnitLine(f)
%   plots a unit line on the current plot
%   

if ~exist('h','var')
    h = gca;
end

hold on;

xLims = xlim;
yLims = ylim;

maxX = max(abs(xLims));
maxY = max(abs(yLims));

maxVal = max([maxX maxY]);

plot([-maxVal, maxVal], [-maxVal, maxVal],'-.k')

xlim([xLims]);
xlim([yLims]);

end

