function [map_interp] = BlueWhiteRedDavide3(nPoints)
% continuous blue to red colourmap, copied from Davide

if ~exist('nPoints','var')
    nPoints = 400;
end

currCols = caxis;

% intrPs = 1:0.01:4.99;

map=[0 0 0.5; 0.1 0.1 0.8; 1 1 1; 1 1 1; 1 1 1; 0.8 0.1 0.1; 0.5 0 0]; % colorbar blu--->bianco--->rosso

mapS = size(map,1);


intrPs = linspace(1,mapS - 0.001,nPoints);

map_interp=[interp1(1:mapS,map(:,1),intrPs)',interp1(1:1:mapS,map(:,2),intrPs)',...
    interp1(1:1:mapS,map(:,3),intrPs)'];


