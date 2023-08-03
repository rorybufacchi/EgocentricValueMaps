function h = hline(y, varargin)

ax = axis;
h = plot(ax(1:2), [y y], varargin{:});