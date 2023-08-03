function h = vline(x, varargin)

ax = axis;
h = plot([x x], ax(3:4), varargin{:});