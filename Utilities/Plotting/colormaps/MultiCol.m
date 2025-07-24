function c = MultiCol(m, colCel)
% c = whitetocol(m, fC, splitLocs)
%  - m is number of divisions
%  - sC is the starting colour
%  - fC is the final colour
%  - splitLocs is a vector defining the locations of color transitions
%
% continuous white to some other color colormap

if ~exist('m','var')
    m=100;
end

nCols = numel(colCel);

c = [];
for iCol = 1:nCols - 1
% % %     c = [c; coltocol(round(m./(nCols - 1)),colCel{iCol},colCel{iCol+1})];

    for iC = 1:3
        tmpMap(:,iC) =  linspace(colCel{iCol}(iC), colCel{iCol+1}(iC), round(m./(nCols - 1)));
    end

    c = [c; tmpMap];
end