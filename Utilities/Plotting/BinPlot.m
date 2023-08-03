function [binX, binY, binYsd, fS] = BinPlot(x,y,fS,varargin)
%[] = BinPlot(xDat,yDat,fS)
%
%   Inputs:
%   - fS            figure settings
%   - fS.pltType    can be 'plot','bar','shadeplot'
%   - fS.c          colour of shade; [r g b]
%   - fS.eBars      error bars; 'SD' or 'centile'


if ~isfield(fS,'c')
    fS.c = [0 0 0.7];
end

if ~isfield(fS,'nBins')
    fS.nBins = 5;
end

if ~isfield(fS,'pltType')
%     fS.pltType = 'shadeplot';
    fS.pltType = 'plot';
end

if ~isfield(fS,'eBars')
    fS.eBars = 'SD';
end

% bins and edges
% [bin,edg] = discretize(x,fS.nBins);
[bin,edg] = discretize(x,linspace(min(x(:)),max(x(:)),fS.nBins+1),'IncludedEdge','right');
for iBin=1:fS.nBins
    binYMean(iBin) = nanmean(y(bin==iBin));
    binYMed(iBin) = nanmedian(y(bin==iBin));
    binYsd(iBin) = nanstd(y(bin==iBin));
    binYCentL(iBin) = prctile(y(bin==iBin),25);
    binYCentU(iBin) = prctile(y(bin==iBin),75);
end

binX = (edg(1:end-1)+edg(2:end))./2;

% error bars
switch fS.eBars
    case 'SD'
        binY = binYMean;
        binYlow = binY-binYsd;
        binYup = binY+binYsd;
    case 'centile'
        binY = binYMed;
        binYlow = binYCentL;
        binYup = binYCentU;
end

switch fS.pltType
    case 'plot'
        plot(binX,binY,varargin{:});
    case 'bar'
        bar((edg(1:end-1)+edg(2:end))./2,binY,varargin{:});
        hold on
        for iX = 1:length(binX)
            plot([binX(iX) binX(iX)] , [binYlow(iX) binYup(iX)],'k',varargin{:});
        end
%         hline(
    case 'shadeplot'
%         PlotShade(binX,binY,binYlow,binYup,fS.c,varargin{:})
        ShadedPlot(binX,binY,binYlow,binYup,fS)
end

end

