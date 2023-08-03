function [h] = ImagescInvisNan(alphaVal,varargin)
%[h] = ImagescInvisNan(alphaVal,varargin)
%   Works exaclty like imagesc, but sets NaNs to be invisible by default.
%   alphaval sets the alpha of the non-nan data

if ishandle(varargin{1})
    bI = 1; % base index
else
    bI = 0;
end

if ~isvector(squeeze(varargin{1+bI})) && ndims(squeeze(varargin{1+bI})) == 2
    
    h = imagesc(varargin{:},'AlphaData',~isnan(varargin{1+bI}).*alphaVal);
    
elseif ~isvector(squeeze(varargin{3+bI})) && ndims(squeeze(varargin{3+bI})) == 2
    
    h = imagesc(varargin{:},'AlphaData',~isnan(varargin{3+bI}).*alphaVal);
    
else
    
    error('Can''t find the data to feed to imagesc')
    
end

end

