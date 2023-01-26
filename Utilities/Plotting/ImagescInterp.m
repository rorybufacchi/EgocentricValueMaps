function [New_Image] = ImagescInterp(data,interpFactor,varargin)
%[New_Image] = ImagescInterp(data,interpFactor,varargin)
%   % varargin 1 should be x points. 2 should be y points and 3 should be
%   the plot flag (0 if you don't want to plot)

plotFl = 1;

if nargin==2
    yPlt=1:size(data,2);
    xPlt=1:size(data,1);
elseif numel(varargin) >= 2
    yPlt=varargin{1};
    xPlt=varargin{2};
    if numel(varargin) == 3
        plotFl = varargin{3};
    end
end

class_of_data = class(data);
[x y] = meshgrid([1:size(data,2)],[1:size(data,1)]);
[xi yi] = meshgrid(1:interpFactor:size(data,2), 1:interpFactor:size(data,1));
New_Image = cast(interp2(x,y,double(data),xi,yi,'linear'),class_of_data);
if plotFl ==1
    imagesc(yPlt,xPlt,New_Image);
end

end

