tic


% Random stimulus dynamics
% % % rSpr = [-1 0 1];
% % % rSprPr = gaussmf(rSpr,[1 0]) ./ sum(gaussmf(rSpr,[1 0]));
% % % cSpr = [-1 0 1]; %%% cSpr = [ -1 0 1 ];
% % % cSprPr = gaussmf(cSpr,[1 0]) ./ sum(gaussmf(cSpr,[1 0]));
% % % zSpr = [-1 0 1]; %%% zSpr = [ -1 0 1 ];
% % % zSprPr = gaussmf(zSpr,[1 0]) ./ sum(gaussmf(zSpr,[1 0]));
rSpr = [0];
rSprPr = gaussmf(rSpr,[1 0]) ./ sum(gaussmf(rSpr,[1 0]));
cSpr = [0]; 
cSprPr = gaussmf(cSpr,[1 0]) ./ sum(gaussmf(cSpr,[1 0]));
zSpr = [0];
zSprPr = gaussmf(zSpr,[1 0]) ./ sum(gaussmf(zSpr,[1 0]));



% Sensory stimulus spreading: from Straka et al ( * 1/2 because their units
% were seconds, and ours are half-seconds
rSnsSpr = [0]; 
cSnsSpr = [0]; 
zSnsSpr = [0]; 
% Set uncertainties to Straka Noel Hoffmann values
rSnsSprPr   = 1; 
cSnsSprPr   = 1; %gaussmf(cSnsSpr,[czSigma 0]) ./ sum(gaussmf(cSnsSpr,[czSigma 0]));
zSnsSprPr   = 1; %gaussmf(zSnsSpr,[czSigma 0]) ./ sum(gaussmf(zSnsSpr,[czSigma 0]));

doubleUnc = 0;


% -------------------------------------------------------------------------
sHnd2                 = sHnd;
sHnd2.clc.baseVel     = [0 0 0];
sHnd2.clc.thinLimbsFl = 1;
% #################################################################
% =========================================================
% HAND by chest, TOWARDS
sHnd2.lp.alg = 'Q';
sHnd2.clc.stimDynams = @(pos) pos + sHnd2.clc.baseVel;

sHnd2.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sHnd2.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

sHnd2.clc.randSpread = {rSpr   cSpr   zSpr}; 
sHnd2.clc.spreadProb = {rSprPr cSprPr zSprPr};


% % % sHnd2.wrld.size     = [61 27 1];
% % % sHnd2.clc.nearPos   = [49 11 1];
% % % % % % sHnd2.clc.startSR   = [49 49 49 49 49 49 49 49 49];
% % % % % % sHnd2.clc.startSC   = [10 11 12 10 11 12 10 11 12];
% % % % % % sHnd2.clc.startSZ   = [ 1  1  1  2  2  2  3  3  3];
% % % sHnd2.clc.startSR   = [49 49 49 ];
% % % sHnd2.clc.startSC   = [10 11 12 ];
% % % sHnd2.clc.startSZ   = [ 1  1  1 ];


% % % sHnd2.wrld.size     = [61 27 3];
% % % sHnd2.clc.nearPos   = [30 15 2];
% % % sHnd2.clc.startSR   = [30 30 30];
% % % sHnd2.clc.startSC   = [14 15 16];
% % % sHnd2.clc.startSZ   = [ 2  2  2];

% % % sHnd2.wrld.size     = [61 27 3];
% % % sHnd2.clc.nearPos   = [30 15 2];
% % % sHnd2.clc.startSR   = [29 29 29 30 30 30 31 31 31];
% % % sHnd2.clc.startSC   = [14 15 16 14 15 16 14 15 16];
% % % sHnd2.clc.startSZ   = [ 2  2  2  2  2  2  2  2  2];


sHnd2.wrld.size     = [61 61 3];
sHnd2.clc.nearPos   = [30 30 2];
sHnd2.clc.startSR   = [29 29 29 30 30 30 31 31 31];
sHnd2.clc.startSC   = [29 30 31 29 30 31 29 30 31];
sHnd2.clc.startSZ   = [ 2  2  2  2  2  2  2  2  2];



% % % sHnd2.clc.actConsequence = [0  0  0 ; ... % action 1 stay
% % %                             0  1  0 ; ... % action 2 left SLOW
% % %                             0 -1  0 ; ... % action 3 right SLOW
% % %                             0  0  1 ; ... % action 4 up SLOW
% % %                             0  0 -1 ; ... % action 5 down SLOW
% % %                            -1  0  0 ; ... % action 6 forward SLOW
% % %                             1  0  0 ; ... % action 7 backward SLOW
% % %                             0  2  0 ; ... % action 2 left MEDIUM
% % %                             0 -2  0 ; ... % action 3 right MEDIUM
% % %                             0  0  2 ; ... % action 4 up MEDIUM
% % %                             0  0 -2 ; ... % action 5 down MEDIUM
% % %                            -2  0  0 ; ... % action 12 forward MEDIUM
% % %                             2  0  0 ; ... % action 13 backward MEDIUM
% % %                             0  3  0 ; ... % action 2 left MEDIUM
% % %                             0 -3  0 ; ... % action 3 right MEDIUM
% % %                             0  0  3 ; ... % action 4 up MEDIUM
% % %                             0  0 -3 ; ... % action 5 down MEDIUM
% % %                            -3  0  0 ; ... % action 18 forward FAST
% % %                             3  0  0 ; ... % action 19 backward FAST
% % %                             0  4  0 ; ... % action 2 left MEDIUM
% % %                             0 -4  0 ; ... % action 3 right MEDIUM
% % %                             0  0  4 ; ... % action 4 up MEDIUM
% % %                             0  0 -4 ; ... % action 5 down MEDIUM
% % %                            -4  0  0 ; ... % action 18 forward FAST
% % %                             4  0  0 ; ... % action 19 backward FAST
% % %                             0  5  0 ; ... % action 2 left MEDIUM
% % %                             0 -5  0 ; ... % action 3 right MEDIUM
% % %                             0  0  5 ; ... % action 4 up MEDIUM
% % %                             0  0 -5 ; ...
% % %                            -5  0  0 ; ... % action 18 forward FAST
% % %                             5  0  0 ;     % action 19 backward FAST
                            

% % % sHnd2.clc.actConsequence = [0  0  0 ; ... % action 1 stay
% % %                             0  1  0 ; ... % action 2 left SLOW
% % %                             0 -1  0 ; ... % action 3 right SLOW
% % %                            -1  0  0 ; ... % action 6 forward SLOW
% % %                             1  0  0 ; ... % action 7 backward SLOW
% % %                             0  2  0 ; ... % action 2 left MEDIUM
% % %                             0 -2  0 ; ... % action 3 right MEDIUM
% % %                            -2  0  0 ; ... % action 12 forward MEDIUM
% % %                             2  0  0 ; ... % action 13 backward MEDIUM
% % %                             0  3  0 ; ... % action 2 left MEDIUM
% % %                             0 -3  0 ; ... % action 3 right MEDIUM
% % %                            -3  0  0 ; ... % action 18 forward FAST
% % %                             3  0  0 ; ... % action 19 backward FAST
% % %                             0  4  0 ; ... % action 2 left MEDIUM
% % %                             0 -4  0 ; ... % action 3 right MEDIUM
% % %                            -4  0  0 ; ... % action 18 forward FAST
% % %                             4  0  0 ; ... % action 19 backward FAST
% % %                             0  5  0 ; ... % action 2 left MEDIUM
% % %                             0 -5  0 ; ... % action 3 right MEDIUM
% % %                            -5  0  0 ; ... % action 18 forward FAST
% % %                             5  0  0 ];    % action 19 backward FAST


% With Diagonal actions
sHnd2.clc.actConsequence = [0  0  0 ; ... % action 1 stay
                            0  1  0 ; ... % action 2 left SLOW
                            0 -1  0 ; ... % action 3 right SLOW
                           -1  0  0 ; ... % action 6 forward SLOW
                            1  0  0 ; ... % action 7 backward SLOW
                            1  1  0 ; ... % action 2 left SLOW ------- DIAGONAL ------
                           -1 -1  0 ; ... % action 3 right SLOW ------- DIAGONAL ------
                           -1  1  0 ; ... % action 6 forward SLOW ------- DIAGONAL ------
                            1 -1  0 ; ... % action 7 backward SLOW ------- DIAGONAL ------
                            0  2  0 ; ... % action 2 left MEDIUM
                            0 -2  0 ; ... % action 3 right MEDIUM
                           -2  0  0 ; ... % action 12 forward MEDIUM
                            2  0  0 ; ... % action 13 backward MEDIUM
                            2  2  0 ; ... % action 2 left SLOW ------- DIAGONAL ------
                           -2 -2  0 ; ... % action 3 right SLOW ------- DIAGONAL ------
                           -2  2  0 ; ... % action 6 forward SLOW ------- DIAGONAL ------
                            2 -2  0 ; ... % action 7 backward SLOW ------- DIAGONAL ------
                            0  3  0 ; ... % action 2 left MEDIUM
                            0 -3  0 ; ... % action 3 right MEDIUM
                           -3  0  0 ; ... % action 18 forward FAST
                            3  0  0 ; ... % action 19 backward FAST
                            3  3  0 ; ... % action 2 left SLOW ------- DIAGONAL ------
                           -3 -3  0 ; ... % action 3 right SLOW ------- DIAGONAL ------
                           -3  3  0 ; ... % action 6 forward SLOW ------- DIAGONAL ------
                            3 -3  0 ; ... % action 7 backward SLOW ------- DIAGONAL ------
                            0  4  0 ; ... % action 2 left MEDIUM
                            0 -4  0 ; ... % action 3 right MEDIUM
                           -4  0  0 ; ... % action 18 forward FAST
                            4  0  0 ; ... % action 19 backward FAST
                            4  4  0 ; ... % action 2 left SLOW ------- DIAGONAL ------
                           -4 -4  0 ; ... % action 3 right SLOW ------- DIAGONAL ------
                           -4  4  0 ; ... % action 6 forward SLOW ------- DIAGONAL ------
                            4 -4  0 ; ... % action 7 backward SLOW ------- DIAGONAL ------
                            0  5  0 ; ... % action 2 left MEDIUM
                            0 -5  0 ; ... % action 3 right MEDIUM
                           -5  0  0 ; ... % action 18 forward FAST
                            5  0  0 ];    % action 19 backward FAST

                            
                            

%-----------------------------
% POS REW
% REWARD SIZE
sHnd2.clc.startRew          = 1;
sHnd2.clc.nSteps            = 2;
sHnd2.clc.nReps             = 2;
sHnd2.clc.stepUpdateFl      = 0;
sHnd2.clc.rewardInterceptFl = 0;


[newQ ] = CalcQDirect(sHnd2);

% store q values and attributes
cQ = 1;
allQ2.qVals{cQ}       = newQ;
allQ2.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ2.rew(cQ)         = sHnd2.clc.startRew; % towards
allQ2.bodyPart{cQ}    = 'HandByChestSARSA'; 
allQ2.centerpos{cQ}   = sHnd2.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% Set the stay q-values to the reward under the hand
for iVol = 1:length(sHnd2.clc.startSZ)
    allQ2.qVals{cQ}(1,sHnd2.clc.startSR(iVol),sHnd2.clc.startSC(iVol),sHnd2.clc.startSZ(iVol)) = sHnd2.clc.startRew;
end

timSpent = toc

% % % % Hopefully final save
% % % save('Data\qValsForTableExpt.mat','allQ2')

figure; imagesc(squeeze(max(allQ2(1,:).qVals{1}(:,:,:,2))))

figure; imagesc(squeeze(mean(allQ2(1,:).qVals{1}(:,:,:,2))))
