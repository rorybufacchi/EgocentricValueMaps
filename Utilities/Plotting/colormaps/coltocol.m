function c = coltocol(m, sC, fC)
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

c1 = coltowhite(m,sC);
c2 = whitetocol(m,fC);
c = [c1; c2];