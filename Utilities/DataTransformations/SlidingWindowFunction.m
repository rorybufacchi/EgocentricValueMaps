function [xQuery, yAvg] = SlidingWindowAverage(x, y, varargin)
%
% [xQuery, yAvg] = SlidingWindowAverage(x, y, function, varargin)
%
% Description:
% Computes a function over a y-sliding window for specified x-query points.
%
% Inputs:
% x          - The original x-coordinate data points.
% y          - The original y-coordinate data points.
% slFun      - The function that will be applied
% 
% xQuery     - The x-coordinates at which to calculate the averaged y-values.
% window_size- The size of the sliding window on the x-axis.
%
% Outputs:
% xQuery     - The x-coordinates at which the averaged y-values are calculated.
% yAvg       - The averaged y-values corresponding to xQuery.
%
% Example usage:
% [xQuery, yAvg] = sliding_window_average(x, y, xQuery, window_size);

minX = nanmin(x(~isinf(x)));
maxX = nanmax(x(~isinf(x)));

% Check for optional xQuery and window_size
if nargin < 3
    slFun = @(x) mean(x);
else
    slFun = varargin{1};
end

% Check for optional xQuery and window_size
if nargin < 4
    xQuery = linspace(minX, maxX, 100);
else
    xQuery = varargin{1};
end

if nargin < 5
    window_size = (maxX - minX) ./4;  % Default window size
else
    window_size = varargin{2};
end

% Sort the x and y arrays
[x_sorted, sort_idx] = sort(x);
y_sorted = y(sort_idx);

% Initialize output
yAvg = zeros(size(xQuery));

% Slide the window across the specified xQuery points and compute the average
for i = 1:length(xQuery)
    x_center = xQuery(i);
    x_start = x_center - window_size / 2;
    x_end = x_center + window_size / 2;

    % Find points that fall within the window
    idx = find(x_sorted >= x_start & x_sorted <= x_end);

    % If there are no points within this window, set the average to NaN
    if isempty(idx)
        yAvg(i) = NaN;
        continue;
    end

    % Compute the average of y-values within the window
    yAvg(i) = slFun(y_sorted(idx));
end

end