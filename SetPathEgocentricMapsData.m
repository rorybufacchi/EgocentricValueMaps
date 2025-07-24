function [basePath] = SetPathEgocentricMapsData()
% [foldName] = SetPathEgocentricMaps()
%   Adds the data path chosen by the user.

addedPaths  = {'\Utilities','\GeneratedData','\Stats','\Figures','\Main'};

msgbox(     'Please select folder for the (pre-)generated data');
pause(      0.5);
basePath    = uigetdir([],'Please select folder to set path');

addpath(    genpath(basePath));

end