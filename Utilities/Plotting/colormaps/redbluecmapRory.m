function c = redbluecmap(m,n,varargin)
%REDBLUECMAPRORY creates a red and blue colormap.
%
%   REDBLUECMAPRORY(M,N) returns an N-by-3 matrix containing a red and blue
%   diverging color palette. Uses the values from the matlab function but
%   interpolates N times between each
%
%   Example:
% 
%       % Reset the colormap of the current figure, type
%             colormap(redbluecmap)
%
%   See also CLUSTERGRAM, COLORMAP, COLORMAPEDITOR, REDGREENCMAP.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2009/01/30 14:41:04 $

% Reference: 
% http://colorbrewer.org.

%== Setting default
if nargin < 1 || isempty(m) || ~isnumeric(m)
    m = 11; n = 10;
end

if nargin < 2 || isempty(n) || ~isnumeric(n)
    n = 10;
end

if ~isscalar(m)
    m = m(:);
end

m = max(abs(fix(m)), 3);
m = min(m, 11);

switch (m)
    case 3
        c = [239	138     98;
             247	247     247;
             103	169     207];
    case 4
        c = [202	0       32;
             244	165     130;
             146	197     222;
             5      113     176];
    case 5
        c = [202	0       32;
             244	165     130;
             247	247     247;
             146	197     222;
             5      113     176];
    case 6
        c = [178	24      43;
             239	138     98;
             253	219     199;
             209	229     240;
             103	169     207;
             33     102     172];
    case 7
        c = [178	24      43;
             239	138     98;
             253	219     199;
            247     247     247;
            209     229     240;
            103     169     207;
            33      102     172];
    case 8
        c = [178	24      43;
             214	96      77;
             244	165     130;
             253	219     199;
             209	229     240;
             146	197     222;
             67     147     195;
             33     102     172];
    case 9
        c = [178	24      43;
             214	96      77;
             244	165     130;
             253	219     199;
             247	247     247;
             209	229     240;
             146	197     222;
             67     147     195;
             33     102     172];
    case 10
        c = [103	0       31;
            178     24      43;
            214     96      77;
            244     165     130;
            253     219     199;
            209     229     240;
            146     197     222;
            67	    147     195;
            33      102     172;
            5       48      97];
    case 11
        c = [103    0       31;
            178     24      43;
            214     96      77;
            244     165     130;
            253     219     199;
            247     247     247;
            209     229     240;
            146     197     222;
            67      147     195;
            33      102     172;
            5       48      97];
end

for iM = 1:m-1
    for iC = 1:3
        c2((iM-1).*n +[1:n],iC) = linspace(c(iM,iC),c(iM+1,iC),n);
    end
end



c = flipud(c2/255);