function c = whitetocol(m,fC,varargin)
% c = whitetocol(m,fC,varargin)
%  - m is number of divisions
%  - fc is colour vector
%
% continuous white to some other colour colourmap
% fC is the final colour


if ~exist('m','var')
    m=100;
end

if ~exist('fC','var')
    fC=[0 1 0];
end

for iC=1:3
    c(:,iC)=linspace(1,fC(iC),m);
end
