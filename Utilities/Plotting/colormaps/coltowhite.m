function c = coltowhite(m,fC,varargin)
% c = coltowhite(m,fC,varargin)
% continuous white to some other colour colourmap
% fC is the final colour


if ~exist('m','var')
    m=100;
end

if ~exist('fC','var')
    fC=[0 1 0];
end

for iC=1:3
c(:,iC)=linspace(fC(iC),1,m);
end
