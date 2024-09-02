%% -------------------------------------------------------------------------
cQ = 142;


% Sensory stimulus spreading: from Straka et al ( * 1/2 because their units
% were seconds, and ours are half-seconds
rSnsSpr = 0; 
cSnsSpr = 0; 
zSnsSpr = 0; 

% Set uncertainties to Straka Noel Hoffmann values
rSnsSprPr   = 1;
cSnsSprPr   = 1;
zSnsSprPr   = 1;

doubleUnc = 0;

sBack = s;
% Because the body is so small compared to the whole modelled volume, don't
% set all the points behind the body to be rewarding: that would simulate a
% human train, not just a human
sBack.clc.thinLimbsFl = 1;

sBack.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sBack.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

sBack.wrld.size = [91 90 31];
sBack.clc.nearPos = [sBack.wrld.size(1)-18 17 13]';

sBack.clc.nReps = 1;


sBack.clc.gammaVal   = 0.8;
% Speed is 1.667 m/s, but direction is towards the body, at 45 degrees
% That's almost 12 blocks / second in each direction: ((5/3) .* cos(pi/4)) ./ 0.1
sBack.clc.baseVel    = [12 -12 0]; 


% Random stimulus dynamics
rSpr = [-1 0 1];
rSprPr = gaussmf(rSpr,[0.5 0]) ./ sum(gaussmf(rSpr,[0.5 0]));
cSpr = [-1 0 1]; 
cSprPr = gaussmf(cSpr,[0.5 0]) ./ sum(gaussmf(cSpr,[0.5 0]));
zSpr = [-1 0 1]; 
zSprPr = gaussmf(zSpr,[0.5 0]) ./ sum(gaussmf(zSpr,[0.5 0]));

% x y z, Deterministic stimulus dynamics
sBack.clc.stimDynams =     @(pos) pos + sBack.clc.baseVel; % For approaching, set speed positive
sBack.clc.randSpread =     {rSpr cSpr zSpr}; % Put a little biit of x and z variability in? Kind of arbitrary
sBack.clc.spreadProb =     {rSprPr cSprPr  zSprPr}; % x y z, probabilities of spread


% NOTE: If the limb can move at different speeds, BOTH speeds have to be
% put in as potential actions, because the model doesn't have any collision
% calculation that takes into account overshoot

% FIRST CALCULATE FOR BODY - moves more slowly than limb, let's say. Also
% ONLY has negative potential rewards
sBack.clc.actConsequence = [ 0  0  0 ; ... % action 1 stay
    0  1  0 ; ... % action 2 left
    0 -1  0 ; ... % action 3 right
    0  0  1 ; ... % action 4 up
    0  0 -1 ; ... % action 5 down
    -1  0  0 ; ... % action 6 forward
    1  0  0 ];    % action 7 backward-

% Location and size of limb, and where the Q-values are calculated FROM
sBack.clc.startSR = []; sBack.clc.startSC = []; sBack.clc.startSZ = [];


bdyRs = sBack.clc.nearPos(1) + [1 0 0 0 1];
bdyCs = sBack.clc.nearPos(2) + [-2:2];
bdyZs = sBack.clc.nearPos(3) + [-2:2];
iRew = 0; % Initialise rewarded block counter
for  iZ = 1:length(bdyZs)
    for  iC = 1:length(bdyCs)
        iRew = iRew + 1;
        sBack.clc.startSR(iRew) =  bdyRs(iC);
        sBack.clc.startSC(iRew) =  bdyCs(iC);
        sBack.clc.startSZ(iRew) =  bdyZs(iZ);
    end
end

%% -----------------------------
% POS REW
% REWARD SIZE
sBack.clc.startRew   =  1;
[newQ ] = CalcQDirect(sBack);
% store q values and attributes
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %sBack.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sBack.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'Back'; 
allQ.centerpos{cQ}   = sBack.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central



% -----------------------------
% NEG REW
% REWARD SIZE
cQ = cQ + 1;
sBack.clc.startRew   =  -1;
[newQ ] = CalcQDirect(sBack);
% store q values and attributes
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %sBack.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sBack.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'Back'; 
allQ.centerpos{cQ}   = sBack.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%% Also calculate SARSA HP, and MI for the back ================================


% ================================= SARSA =================================
cQ = 145;

sBack.lp.alg = 'SARSA';
sBack.clc.nReps = 2; 

% POS REW 
% REWARD SIZE
cQ = cQ + 1;
sBack.clc.startRew   =  1;
[newQ ] = CalcQDirect(sBack);
% store q values and attributes
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %sBack.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sBack.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'BackSARSA'; 
allQ.centerpos{cQ}   = sBack.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW
% REWARD SIZE
cQ = cQ + 1;
sBack.clc.startRew   =  -1;
[newQ ] = CalcQDirect(sBack);
% store q values and attributes
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %sBack.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sBack.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'BackSARSA'; 
allQ.centerpos{cQ}   = sBack.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central



%% ================================ HITPROB ================================
cQ = 147;

sBackHP                     = sBack;
sBackHP.clc.actConsequence  = [0 0 0]; 
sBackHP.clc.gammaVal        = 1;
sBackHP.clc.startRew        = 1;

sBackHP.clc.nReps = 1; 

% Hit probability - body
% -----------------------------
% NOTE: deltaT is 0.5, so I have to halve all the velocities and spreads?


% Random stimulus dynamics
rSpr = [-1 0 1];
rSprPr = gaussmf(rSpr,[1 0]) ./ sum(gaussmf(rSpr,[0.5 0]));
cSpr = [-1 0 1]; %%% cSpr = [ -1 0 1 ];
cSprPr = gaussmf(cSpr,[1 0]) ./ sum(gaussmf(cSpr,[0.5 0]));
zSpr = [-1 0 1]; %%% zSpr = [ -1 0 1 ];
zSprPr = gaussmf(zSpr,[1 0]) ./ sum(gaussmf(zSpr,[0.5 0]));

% divide probabilities over 2 because deltaT = 1/2
rSprPr = rSprPr./([2 1 2]) + rSprPr(1).*([0 1 0]); 
cSprPx = cSprPr./([2 1 2]) + cSprPr(1).*([0 1 0]); 
zSprPr = zSprPr./([2 1 2]) + zSprPr(1).*([0 1 0]); 

% x y z, Deterministic stimulus dynamics
sBackHP.clc.stimDynams =     @(pos) pos + s.clc.baseVel./2; % For approaching, set speed positive
sBackHP.clc.randSpread =     {rSpr cSpr zSpr}; % Put a little biit of x and z variability in? Kind of arbitrary
sBackHP.clc.spreadProb =     {rSprPr cSprPr  zSprPr}; %


% Sensory stimulus spreading: from Straka et al ( * 1/2 because their units
% were seconds, and ours are half-seconds, and * 1/2 again because the back is in
% voxels of 10cm, whereas all other bodyparts have voxels of 5cm
rSnsSpr = [-7:7]; 
cSnsSpr = [-2:2]; 
zSnsSpr = [-2:2]; 

% Set uncertainties to Straka Noel Hoffmann values
rSigma      = sqrt(2.5.^2 + 30.^2) ./ sBackHP.clc.binW ./4;
rSnsSprPr   = gaussmf(rSnsSpr,[rSigma 0]) ./ sum(gaussmf(rSnsSpr,[rSigma 0])); %sqrt(2.5.^2 + 30.^2)./s.clc.binW
czSigma     = sqrt(5.^2 + 5.^2) ./ sBackHP.clc.binW ./4;
cSnsSprPr   = gaussmf(cSnsSpr,[czSigma 0]) ./ sum(gaussmf(cSnsSpr,[czSigma 0]));
zSnsSprPr   = gaussmf(zSnsSpr,[czSigma 0]) ./ sum(gaussmf(zSnsSpr,[czSigma 0]));

sBackHP.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sBackHP.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};


sBackHP.clc.baseVel    = [6 -6 0]; % Same velocity as Q values, but now in units of half-seconds
sBackHP.clc.stimDynams = @(pos) pos + sBackHP.clc.baseVel;


sBackHP.clc.nReps = 1; 

sBackHP.clc.stepUpdateFl = 1; % whether to run a full sweep update or not
sBackHP.clc.nSteps = 1;

% Set the importance of false positives and false negatives
sForUtility.clc.FN = 5;
sForUtility.clc.FP = 1;

backHP = CalcQDirect(sBackHP);
cQ = cQ + 1;
allQ.qVals{cQ}       = HitProbToUtility(backHP,sForUtility);
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sBackHP.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'BackHP'; 
allQ.centerpos{cQ}   = sBackHP.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central



%% ======================== Multisensory Integration =======================
% Parameters are based of Bertoni et al 2021
cQ = 148;

% -----------------------------
sBackMI          = sBackHP;
sBackMI.clc.binW = 10;

% Stimulus should not move at all
rSpr    = 0;
rSprPr  = 1;
cSpr    = 0;
cSprPr  = 1;
zSpr    = 0;
zSprPr  = 1;

% x y z, Deterministic stimulus dynamics
sBackMI.clc.stimDynams =     @(pos) pos; % No speed
sBackMI.clc.randSpread =     {rSpr cSpr zSpr}; % Put a little biit of x and z variability in? Kind of arbitrary
sBackMI.clc.spreadProb =     {rSprPr cSprPr  zSprPr}; %

% Sensory uncertainty
rSnsSprMI = [-7:7]; 
cSnsSprMI = [-2:2]; 
zSnsSprMI = [-2:2]; 

% Set uncertainties to Bertoni 2021 or Noel 2018 values
% Bertoni 2021: 13cm for proprioception and 11cm for visual
% Noel 2018: 36 mm --> that's too small, it will never fit any of the real data
% (But divide by 2, because the 
rSigmaMI      = sqrt(13.^2 + 11.^2) ./ sBackMI.clc.binW ;
rSnsSprMIPr   = gaussmf(rSnsSprMI,[rSigmaMI 0]) ./ sum(gaussmf(rSnsSprMI,[rSigmaMI 0])); %sqrt(2.5.^2 + 30.^2)./s.clc.binW
czSigmaMI     = sqrt(13.^2 + 11.^2) ./ sBackMI.clc.binW ;
cSnsSprMIPr   = gaussmf(cSnsSprMI,[czSigmaMI 0]) ./ sum(gaussmf(cSnsSprMI,[czSigmaMI 0]));
zSnsSprMIPr   = gaussmf(zSnsSprMI,[czSigmaMI 0]) ./ sum(gaussmf(zSnsSprMI,[czSigmaMI 0]));


sBackMI.clc.sensSpread = {rSnsSprMI cSnsSprMI zSnsSprMI};
sBackMI.clc.sensProb   = {rSnsSprMIPr cSnsSprMIPr zSnsSprMIPr};


sBackMI.clc.baseVel    = [0 0 0]; % Same velocity as Q values, but now in units of half-seconds
sBackMI.clc.stimDynams = @(pos) pos + sBackMI.clc.baseVel;


sBackMI.clc.nReps = 1; 

sBackMI.clc.stepUpdateFl = 1; % whether to run a full sweep update or not
sBackMI.clc.nSteps = 1;

cQ = cQ + 1;
allQ.qVals{cQ}       = CalcQDirect(sBackMI);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sBackMI.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'BackMI'; 
allQ.centerpos{cQ}   = sBackMI.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%% =============================== Distance ===============================
cQ = 150;
allQ.qVals{cQ}       = permute(CalcDistToBody(sBackMI),[4 1 2 3]);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sBackMI.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'BackDist'; 
allQ.centerpos{cQ}   = sBackMI.clc.nearPos;


%%

% Hopefully final save
% save('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Data\NewWithHitProbs_MultSensInt_AndDist_AndUncertainQ_andSARSAQ_V3_WITHBACKDAT_ThinBack.mat','allQ','-v7.3')
save('Data\NewWithHitProbs_MultSensInt_AndDist_AndUncertainQ_AndSARSA_FINAL_I_hope_v4.mat','-v7.3')

%% Functions 

function observerOutput = ObserverPerspective(s,qVals)
% observerOutput = ObserverPerspective(s,qVals)
% Inverts the observerdynamics, so that it calculates Q or hitprob from the
% perspective of an observer delivering the stimuli
% inverse dynamics for P(o|s)

sobs = s;
sobsBack.clc.sensSpread = cellfun(@(x) -x, sBack.clc.sensSpread, 'UniformOutput', false); 
% no movement, only observation uncertainty allowed
sobsBack.clc.baseVel           = [0 0 0]; 
sobsBack.clc.actConsequence    = [0 0 0]; 
sobsBack.clc.randSpread        = {[0] [0] [0]};
sobsBack.clc.spreadProb        = {[1] [1] [1]};
sobsBack.clc.stimDynams        = @(pos) pos + sobsBack.clc.baseVel;
sobsBack.clc.nSteps            = 1;
sobsBack.clc.nReps             = 1; 
sobsBack.clc.stepUpdateFl      = 1; 

% hit encoding from observer's perspective
observerOutput = CalcQDirect(sobs,qVals);
end