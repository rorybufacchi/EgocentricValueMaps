function [] = ShadedPlot(xVals,yValsMean,yValsMin,yValsMax,opts)
%[] = ShadedPlot(xVals,yValsMean,yValsMin,yValsMax,opts)
%   opts.c = colours
%   opts.FaceAlpha is alphaness of the shade
%   opts.PlotMean is a boolean for plotting the mean line

if ~exist('opts','var')
    opts.c = [0.7 0.7 0.7];
end

if ~isfield(opts,'FaceAlpha')
    opts.FaceAlpha = 0.5;
end

if ~isfield(opts,'PlotMean')
    opts.PlotMean = 1;
end


x = [xVals(:); flip(xVals(:))];
y = [yValsMin(:); flip(yValsMax(:))];
patch(x,y,opts.c,'EdgeColor','None','FaceAlpha',opts.FaceAlpha,'HandleVisibility','off'); hold on

switch opts.PlotMean
    case 1
        plot(xVals,yValsMean,'LineWidth',2)
end

end

