% -------------------------------------------------------------------------
function [rr cc zz] = FindLinePixels(dR,dC,dZ)
% Finds pixels on a particular line starting at 0

nSamp = 2 .* max (abs([dR dC dZ]));
rr = sign(dR) .* round(linspace(0,abs(dR),nSamp));
cc = sign(dC) .* round(linspace(0,abs(dC),nSamp));
zz = sign(dZ) .* round(linspace(0,abs(dZ),nSamp));

% Find and remove duplicates
duplNs = [0; sum(abs(diff([rr' cc' zz'])),2) == 0] == 1;
rr(duplNs) = [];
cc(duplNs) = [];
zz(duplNs) = [];
end
% -------------------------------------------------------------------------