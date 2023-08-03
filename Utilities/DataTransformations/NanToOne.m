function [outDat] = NanToOne(inDat,varargin)
%[outDat] = NanToZero(inDat, setVal)
%   Set all NaNs to one, or whatever variable is supplied as setVal

if numel(varargin) == 1
    setVal = varargin{1};
else
    setVal = 1;
end

outDat               = inDat;
outDat(isnan(inDat)) = 1;

end