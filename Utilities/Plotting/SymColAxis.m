function [] = SymColAxis(h)
%[] = SymColAxis(f)
%   Sets the colour axis to be symmetrical
%   
%   - h is an AXIS handle

if ~exist('h','var')
    h = gca;
end

ax = caxis(h);

maxV = max(abs(ax));

newAx = [-maxV maxV];

caxis(h,newAx);


end

