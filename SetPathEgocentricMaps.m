function [basePath addedPaths] = SetPathEgocentricMaps()
% [foldName addedPaths] = SetPathEgocentricMaps()
%   Adds the necessary paths. Returns the name of the base folder, and the
%   list of path names

addedPaths = {'\Utilities','\GeneratedData','\Stats','\Figures','\Main'};

msgbox('Please select folder that will serve as the base file path');
pause(0.5)
basePath = uigetdir([],'Please select BASE FOLDER to set path');


for iPath = 1:length(addedPaths)
    addedPaths{iPath} = [basePath addedPaths{iPath}];
    addpath(genpath( addedPaths{iPath} ));
end

end