function [newQ] = CalcQDirect(s, varargin)
% -------------------------------------------------------------------------
% CalcQDirect(s)
% Calculate Q values, or Hit probabilities

s = DefaultSettings(s);

nA = size(s.clc.actConsequence,1); % Number of actions

if isempty(varargin)
    varargin = {[],[]};
end

% Allow loading of previously calculated Q values
if isempty(varargin{1})
    newQ = zeros([nA s.wrld.size]);
else
    newQ = varargin{1};
    
    % Take into account the possibility of adding actions
    if size(newQ,1) < nA
        actDiff = nA - size(newQ,1);
        newQ(end+ [1 : actDiff],:,:,:) = zeros([actDiff size(newQ,[2 3 4]) ] );
    end
end


% Provide alternative Q values for following a different policy
if isempty(varargin{2})
    altQ = zeros([nA s.wrld.size]);
else
    altQ = varargin{2};
    % Take into account the possibility of adding actions
    if size(altQ,1) < nA
        actDiff = nA - size(altQ,1);
        altQ(end+ [1 : actDiff],:,:,:) = zeros([actDiff size(altQ,[2 3 4]) ] );
    end
end

% first initialise with the contact reward. Also behind the 'skin surface'
% [hence the startSR:end in the 2nd dimension]
for iVol = 1:length(s.clc.startSZ)
    if numel(s.clc.startRew) == 1
        if s.clc.thinLimbsFl == 1
            newQ(:,s.clc.startSR(iVol),s.clc.startSC(iVol),s.clc.startSZ(iVol)) = s.clc.startRew;
        else
            newQ(:,s.clc.startSR(iVol):end,s.clc.startSC(iVol),s.clc.startSZ(iVol)) = s.clc.startRew;
        end
    else
        % If there are multiple rewards, it is the hand sliding over the
        % body situation, so I shouls split the rewards accordingly
        for iSplitInd = 1:length(s.clc.rewSplitInd)
            if iVol >= s.clc.rewSplitInd(iSplitInd) & iVol < s.clc.rewSplitInd(iSplitInd+1)
                if s.clc.thinLimbsFl == 1
                    newQ(:,s.clc.startSR(iVol),s.clc.startSC(iVol),s.clc.startSZ(iVol)) = s.clc.startRew(iSplitInd);
                else
                    newQ(:,s.clc.startSR(iVol):end,s.clc.startSC(iVol),s.clc.startSZ(iVol)) = s.clc.startRew(iSplitInd);
                end
            end
        end
    end
end

% Also initialise the tool rewards, IF there is a tool
if isfield(s.clc,'toolSR')
for iVol = 1:length(s.clc.toolSR)
    newQ(:,s.clc.toolSR(iVol),s.clc.toolSC(iVol),s.clc.toolSZ(iVol)) = s.clc.startRew;
end
end

% Create a variable that doesn't get updated, to be used for checking
% collision only
checkCollision = newQ ~= 0;

% Figure out hwo far from the edges the calculation should stop, based on
% the available actions
maxActDists = max(abs(s.clc.actConsequence),[],1) + 1;
% Add in the effects of world dynamics: Don't want to go so far 'forward in
% the world', that we look for q-values in the columns.
% NOTE: to take into account dynamics that vary by stimulus position, check
% for the displacements in a bunch of different positions
maxActDists = maxActDists + ...
    max(cell2mat(arrayfun(...
    @(x) abs(s.clc.stimDynams(x + [1 2 3]) - (x + [1 2 3]) ) , ...
    [1:15]  , 'UniformOutput' , false)'));

% Take into account random spread of the stimulus
maxActDists = maxActDists + cellfun(@(dim) max(abs(dim)), s.clc.randSpread);

% Take into account perceptual errors
maxActDists = maxActDists + cellfun(@(dim) max(abs(dim)), s.clc.sensSpread);

% next, loop through rows backwards and add sums
for iIt    = 1:s.clc.nReps
for iSweeps = 1:s.clc.nSteps % this can be used to sweep through all values first, before using the updates values to continue Q-value calculation
oldQ = newQ;
for iSR     = s.wrld.size(1) - maxActDists(1) : -1 : maxActDists(1) +1

    % reverse direction to get symmetric result if doing more than 1
    % repetition
    if mod(iIt,2) == 0
        colOrder = [(1 + maxActDists(2)) : (s.wrld.size(2) - maxActDists(2))];
    else
        colOrder = [(s.wrld.size(2) - maxActDists(2)) : -1 : (1 + maxActDists(2))] ;
    end
    
    for iSC = colOrder
% % %         for iSZ = 1 + maxActDists(3) : s.wrld.size(3) - maxActDists(3)
            for iSZ = 1 + (maxActDists(3)-1) : s.wrld.size(3) - (maxActDists(3)-1)

            % Only update Q if it is not an 'originally rewarded' location
            if s.clc.thinLimbsFl == 0
                checkRows = all( [iSC iSZ]' == [s.clc.startSC ; s.clc.startSZ]);
            else
                checkRows = any(all( [iSR iSC iSZ]' == [s.clc.startSR ; s.clc.startSC ; s.clc.startSZ]));
            end
            if iSR >= s.clc.startSR(checkRows)
                % Do nothing
                calcFl = 0; % flag for calculating Q
            elseif isfield(s.clc,'toolSR')
               % Also check whether the current position is a tool position
               if any(all( [iSR iSC iSZ]' == [s.clc.toolSR ; s.clc.toolSC ; s.clc.toolSZ]))
                   calcFl = 0;
               else
                   calcFl = 1;
               end
            else
                calcFl = 1;
            end

            if calcFl == 1;
                nextPos = s.clc.stimDynams([iSR iSC iSZ]);

                % Define random spread of stimulus
                rndSprR = s.clc.randSpread{1};
                rndSprC = s.clc.randSpread{2};
                rndSprZ = s.clc.randSpread{3};

                sprPrR  = s.clc.spreadProb{1};
                sprPrC  = s.clc.spreadProb{2};
                sprPrZ  = s.clc.spreadProb{3};


                % Define sensory spread
                snsSprR = s.clc.sensSpread{1};
                snsSprC = s.clc.sensSpread{2};
                snsSprZ = s.clc.sensSpread{3};

                snsPrR  = s.clc.sensProb{1};
                snsPrC  = s.clc.sensProb{2};
                snsPrZ  = s.clc.sensProb{3};

                % Make exception for the track situation: needs to not have
                % random variance EXCEPT when stimulus is at splitting
                % location
                if isfield(s.clc,'specialTraject')
                    if all([iSR iSC iSZ] == s.clc.specialTraject(3,:))
                        rndSprR = 0;
                        rndSprC = [0 -4];
                        rndSprZ = 0;

                        sprPrR  = 1;
                        sprPrC  = [.5 .5];
                        sprPrZ  = 1;
                    end
                end


                for iAct = 1:nA

                    % DEBUGGING
                    if iAct == 2 ...
                        % & ...
                        %    iSR == 12 & ... 11 & ...
                        %    iSC == 7 & ... 11
                        %    iSZ == 1 %13
                        % any(s.clc.startSR - 1 == iSR) & ...
                        %    any(s.clc.startSC)     == iSC & ...
                        %    any(s.clc.startSZ)     == iSZ,

                        % disp('test')
                    end

                    % Find expected value across future possible states: loop
                    % through all possible next positions and calculate value
                    % times probability of being in that state
                    nextQ = [];
                    for iSprR = 1:length(rndSprR) % row spread
                        for iSnsR = 1:length(snsSprR)
                            nextR = rndSprR(iSprR) + nextPos(1);
                            nextR = nextR + snsSprR(iSnsR);
                        for iSprC = 1:length(rndSprC) % column spread
                            for iSnsC = 1:length(snsSprC)
                                nextC = rndSprC(iSprC) + nextPos(2);
                                nextC = nextC + snsSprC(iSnsC);
                            for iSprZ = 1:length(rndSprZ) % height spread
                                for iSnsZ = 1:length(snsSprZ)
                                    nextZ = rndSprZ(iSprZ) + nextPos(3);
                                    nextZ = nextZ + snsSprZ(iSnsZ);
             
                                % add to possible future qs
                                actNextR = nextR + s.clc.actConsequence(iAct,1);
                                actNextC = nextC + s.clc.actConsequence(iAct,2);
                                actNextZ = nextZ + s.clc.actConsequence(iAct,3);

                                % DO it with collision:
                                % First calculate which blocks are in the line  
                                dR = actNextR - iSR;
                                dC = actNextC - iSC;
                                dZ = actNextZ - iSZ;

                                % Then find all blocks on that line
                                [rR cC zZ] = FindLinePixels(dR,dC,dZ);
                                rR = rR + iSR;
                                cC = cC + iSC;
                                zZ = zZ + iSZ;
                                % And check the Q value at ANY of those
                                % points is equal the starting reward: that
                                % indicates that a collision has happened,
                                % if 'interception' is rewarded, rather
                                % than 'landing on', or 'teleporting
                                % through'
                                tmpQ = NaN;
                                if s.clc.rewardInterceptFl == 1
                                    for iRR = 1:numel(rR)
                                    % Allow for using alternative policy
                                    if s.clc.useAltPolicyFl == 1
                                        [dmyQ, nextAct] = max( altQ(:,...
                                        actNextR , ...
                                        actNextC , ...
                                        actNextZ  ) );
                                        if squeeze( checkCollision(nextAct,...
                                                rR(iRR) , ...
                                                cC(iRR) , ...
                                                zZ(iRR)  ) ) == 1
                                            tmpQ = s.clc.startRew;
                                        end
                                    else
                                        if squeeze( max( checkCollision(:,...
                                                rR(iRR) , ...
                                                cC(iRR) , ...
                                                zZ(iRR)  ) )) == 1
                                            tmpQ = s.clc.startRew;
                                        end
                                    end 
                                    end
                                end

                                % SO if there is any possible collision,
                                % then use that Q value. Otherwise use the
                                % q-value of the place the stimulus ends up
                                % in
                                if isnan(tmpQ)
                                    % Allow for using alternative policy
                                    if s.clc.useAltPolicyFl == 1
                                        [~, nextAct] = max( altQ(:,...
                                        actNextR , ...
                                        actNextC , ...
                                        actNextZ  ) );
                                    else
                                        nextAct = 1:size(oldQ,1);
                                    end
                                    switch s.lp.alg
                                        case 'Q'
                                            tmpQ = squeeze( max( oldQ(nextAct,...
                                                actNextR , ...
                                                actNextC , ...
                                                actNextZ  ) ));
                                        case 'SARSA'
                                            tmpQ = squeeze( mean( oldQ(nextAct,...
                                                actNextR , ...
                                                actNextC , ...
                                                actNextZ  ) ));
                                    end
                                end

                                % Then update the next possible Q value
                                nextQ = [nextQ ; ...
                                    sprPrR(iSprR) .* sprPrC(iSprC) .* sprPrZ(iSprZ) .* ...
                                    snsPrR(iSnsR) .* snsPrC(iSnsC) .* snsPrZ(iSnsZ) .* ...
                                    s.clc.gammaVal .* ...
                                    tmpQ];
 
                            end
                            end
                        end
                        end
                    end
                    end
                    % Sum all Qs into new Q
                    newQ(iAct,iSR,iSC,iSZ) = sum(nextQ);

                    % For faster updates,calculate the Q-values from latest
                    % approximation.For unbiased updates, instead set the
                    % flag below to 1
                    if s.clc.stepUpdateFl == 0
                        oldQ = newQ;
                    end
                end

            end

        end
    end
end
end
iIt
end


% Then set the value of the touch condition back to 0
for iVol = 1:length(s.clc.startSZ)
    if s.clc.thinLimbsFl
        newQ(:,s.clc.startSR(iVol),s.clc.startSC(iVol),s.clc.startSZ(iVol)) = 0;
    else
        newQ(:,s.clc.startSR(iVol):end,s.clc.startSC(iVol),s.clc.startSZ(iVol)) = 0;
    end
end

% Or set tool space to the slightly nearer position
if isfield(s.clc,'toolSR')
    for iVol = 1:length(s.clc.toolSR)
        newQ(:,s.clc.toolSR(iVol),s.clc.toolSC(iVol),s.clc.toolSZ(iVol)) = ...
            newQ(:,s.clc.toolSR(iVol)+1,s.clc.toolSC(iVol),s.clc.toolSZ(iVol));
    end
end


end
% -------------------------------------------------------------------------