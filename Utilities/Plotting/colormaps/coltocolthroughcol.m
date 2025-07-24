function c = coltocolthroughcol(m, sC, mC, fC)
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

c1 = coltocol(round(m./2),sC,mC);
c2 = coltocol(round(m./2),mC,fC);
c = [c1; c2];