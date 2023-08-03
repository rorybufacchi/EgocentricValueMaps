% -------------------------------------------------------------------------
function [rr cc zz] = FindLinePixels(dR,dC,dZ)
% Finds pixels on a particular line starting at 0

nSamp = 2 .* max ([dR dC dZ]);
rr = round(linspace(0,dR,nSamp));
cc = round(linspace(0,dC,nSamp));
zz = round(linspace(0,dZ,nSamp));

% Find and remove duplicates
duplNs = [0; sum(diff([rr' cc' zz']),2) == 0] == 1;
rr(duplNs) = [];
cc(duplNs) = [];
zz(duplNs) = [];
end
% -------------------------------------------------------------------------