function [splineHandle] = PlotHistSpline(h,varargin)
% [splineHandle] = PlotHistSpline(h,varargin)
%   Plots a spline for the histogram h and returns the handle to the spline
%   varargin functions exaclty like the normal plot options

% put an x-point on either end so that the spline goes to 0
x = [h.BinEdges(1) - h.BinWidth./2 (h.BinEdges + h.BinWidth./2) ];
y = [0 h.Values 0];

% Create the spline using the midpoints and heights
spline_fit = spline(x, y);

% Evaluate the spline on a finer grid for a smooth curve
x_fine = linspace(min(x), max(x), 1000);
y_fine = ppval(spline_fit, x_fine);

% Plot the spline on top of the histogram
splineHandle = plot(x_fine, y_fine, varargin{:});

end