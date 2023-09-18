function c = coltocol(m, sC, fC, splitLocs)
% c = whitetocol(m, fC, splitLocs)
%  - m is number of divisions
%  - sC is the starting colour
%  - fC is the final colour
%  - splitLocs is a vector defining the locations of color transitions
%
% continuous white to some other color colormap

if ~exist('m', 'var')
    m = 100;
end

if ~exist('splitLocs', 'var')
    splitLocs = 1;
end

if ~exist('sC', 'var')
    fC = [1 1 1];
end

if ~exist('fC', 'var')
    fC = [0 1 0];
end


% sort splitLocs in ascending order and ensure they are in the range [0, 1]
splitLocs = unique(max(min(splitLocs, 1), 0));
% append 1 at the end, if is is not there
splitLocs = unique([splitLocs(:); 1]);


% Update the colours linearly to match the number of splitting locations
if size(fC,1) ~= numel(splitLocs)
    fCold = fC;
    for iC = 1:3

        % Create new vector of intermediate colours
        tmpFC = linspace(sC(:,iC),fCold(end,iC),numel(splitLocs) + 1);

        % Put the correct colour in the correct place
        for iSL = 2:numel(tmpFC)
            fC(iSL - 1, iC) = tmpFC(iSL);
        end

    end
end


c = zeros(m,3);  % initialize color map

% Set number of points inbetween each colour
splitM = round(m .* splitLocs);
% ensure that the total number of colour points is correct (robust to
% rounding errors)
splitM(end) = m;


% create colour map
for iC=1:3
    currCol = sC(:,iC);
    prvM    = 0;
    for iSL = 1:length(splitLocs)


        cPoints = (1:splitM(iSL)) + prvM;
        c(cPoints,iC) = linspace(currCol,fC(iSL,iC),splitM(iSL));

        currCol = fC(iSL,iC);
        prvM    = splitM(iSL);
    end
end

c


% % % function c = whitetocol(m,fC,splitLocs)
% % % % c = whitetocol(m,fC,varargin)
% % % %  - m is number of divisions
% % % %  - fc is colour vector
% % % %  - 
% % % %
% % % % continuous white to some other colour colourmap
% % % % fC is the final colour
% % % 
% % % 
% % % if ~exist('m','var')
% % %     m=100;
% % % end
% % % 
% % % if ~exist('fC','var')
% % %     fC=[0 1 0];
% % % end
% % % 
% % % for iC=1:3
% % %     c(:,iC)=linspace(1,fC(iC),m);
% % % end



% % % % Put the correct colour in the correct place
% % % currCol = 1;
% % % for iSL = 1:numel(splitLocs)
% % % 
% % %     % Create new vector of intermediate colours
% % %     endCol = 1 - ((1 - fCold(end,iC)) ./ numel(splitLocs) ) .* iSL;
% % %     tmpFC = linspace(1,endCol,numel(splitLocs));
% % % 
% % % 
% % %     fC(iSL,iC) = tmpFC(iSL);
% % % 
% % %     currCol = endCol;
% % % end
