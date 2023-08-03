function [outDat] = NanToZero(inDat)
%[outDat] = NanToZero(inDat)
%   Set all NaNs to zero

outDat               = inDat;
outDat(isnan(inDat)) = 0;

end