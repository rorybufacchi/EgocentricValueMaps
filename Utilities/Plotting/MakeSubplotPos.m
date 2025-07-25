function plotPos    = MakeSubplotPos(fillSize, subSize, subSpace)
plotPos             = (1-fillSize)./2;
for iX              = 2 : length(subSize)
    plotPos             = [plotPos plotPos(end)+subSize(iX-1) + subSpace] ;
end
end
