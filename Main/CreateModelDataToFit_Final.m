% -----
% USE:
% 
% Run this script to generate the Q-values and hit probabilities used to
% fit empirical data from scrath
% -----
%
% =========================================================================
%                              !!! WARNING !!!
% =========================================================================
% Running this whole script uo to a day to run
% It is !!! HIGHLY RECOMMENDED !!! to either: 
%       --------------------------   
% A) Skip using this script, and instead just load the pre-computed 
%    outputs , which are available $$$ in the 'GeneratedData' sub-directory
% OR      -------------
% B) Run each section separately:
%    some sections can take a long time to run
% =========================================================================
%

% First make sure the path is set right though
cFileName       = matlab.desktop.editor.getActiveFilename;
codePath        = cFileName(1:end-33);
addpath(        genpath(codePath));


%% -------------------------------------------------------------------------
% Set base settings

% First run a quick dummy agent to get the world-structure 'w' 
% (Needed for plotting later)
s.rp.maxActions                 = 10; 
s.wrld.size                     = [14 15];
s                               = DefaultSettings(s);
[~, w]                          = RunRLfun(s);
s                               = DefaultSettings(s);

% Set bin width
s.clc.binW                      = 5;



%% -------------------------------------------------------------------------
%  Create theoretical Q-values

% Because (A) the furthest distance for most expts ~= 100cm, (B) the distance 
% between stimulus locations tends to be between 15 and 22 cm, and (C) to 
% allow for velocity effects, we opted to used a 5cm^3 bin spatial size.
% We also made the world 305cm by 135cm by 155cm --> allows hand to be placed 
% to the side, and head to be placed above body
% NOTE: this involves leaving 20cm behind the participants' back


% #################################################################
% Baseline Q-values
% =========================================================
% STIM MOVING TOWARDS TORSO

s.wrld.size     = [61 27 31];
s.clc.nearPos   = [s.wrld.size(1)-12 11 13]';
s.clc.nReps     = 1;
allQ            = table();
s.clc.gammaVal  = 0.8;
s.clc.baseVel   = [4 0 0]; 

% Random stimulus dynamics
rSpr            = [-1 0 1];
rSprPr          = gaussmf(rSpr,[1 0]) ./ sum(gaussmf(rSpr,[1 0]));
cSpr            = [-1 0 1];
cSprPr          = gaussmf(cSpr,[1 0]) ./ sum(gaussmf(cSpr,[1 0]));
zSpr            = [-1 0 1];
zSprPr          = gaussmf(zSpr,[1 0]) ./ sum(gaussmf(zSpr,[1 0]));

% x y z, Deterministic stimulus dynamics
s.clc.stimDynams =     @(pos) pos + s.clc.baseVel; % For approaching, set speed positive
s.clc.randSpread =     {rSpr cSpr zSpr}; % Put a little biit of x and z variability in? Kind of arbitrary
s.clc.spreadProb =     {rSprPr cSprPr  zSprPr}; % x y z, probabilities of spread

% CREATE ACTIONS AVAILABLE TO BODY. It moves more slowly than limb, max 10cm/s.
% NOTE: If the limb can move at different speeds, BOTH speeds have to be
% put in as potential actions, because the model doesn't have any collision
% calculation that takes into account overshoot
sBdy = s;
sBdy.clc.actConsequence = [ 
     0  0  0 ; ... % action 1 stay
     0  1  0 ; ... % action 2 left
     0 -1  0 ; ... % action 3 right
     0  0  1 ; ... % action 4 up
     0  0 -1 ; ... % action 5 down
    -1  0  0 ; ... % action 6 forward
     1  0  0 ];    % action 7 backward-

% Location and size of body part, and where the Q-values are calculated FROM
sBdy.clc.startSR = []; sBdy.clc.startSC = []; sBdy.clc.startSZ = [];
bdyRs               = s.clc.nearPos(1) + [2 1 1 0 0 0 1 1 2];
bdyCs               = s.clc.nearPos(2) + [-4:4];
bdyZs               = s.clc.nearPos(3) + [-5:5];
iRew                = 0; % Initialise rewarded block counter
for  iZ = 1:length(bdyZs)
    for  iC = 1:length(bdyCs)
        iRew = iRew + 1;
        sBdy.clc.startSR(iRew) =  bdyRs(iC);
        sBdy.clc.startSC(iRew) =  bdyCs(iC);
        sBdy.clc.startSZ(iRew) =  bdyZs(iZ);
    end
end

% -----------------------------
% CALC Q FOR POS REW
% REWARD SIZE
sBdy.clc.startRew   =  1;
[newQ ]             = CalcQDirect(sBdy);
% store q values and attributes
allQ.qVals{1}       = newQ;
allQ.dir(1)         = 1; % towards
allQ.rew(1)         = sBdy.clc.startRew;
allQ.bodyPart{1}    = 'Trunk'; 
allQ.centerpos{1}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% CALC Q FOR NEG REW
% REWARD SIZE
sBdy.clc.startRew   =  -1;
[newQ ]             = CalcQDirect(sBdy);
% store q values and attributes
allQ.qVals{2}       = newQ;
allQ.dir(2)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(2)         = sBdy.clc.startRew; % towards
allQ.bodyPart{2}    = 'Trunk'; 
allQ.centerpos{2}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% =========================================================
% STIM MOVING AWAY FROM TORSO
sBdy.clc.stimDynams =     @(pos) pos - s.clc.baseVel; % For receding, set speed negative

% -----------------------------
% POS REW
% REWARD SIZE
sBdy.clc.startRew   = 1;
[newQ ]             = CalcQDirect(sBdy);
% store q values and attributes
allQ.qVals{3}       = newQ;
allQ.dir(3)         = -1; 
allQ.rew(3)         = sBdy.clc.startRew;
allQ.bodyPart{3}    = 'Trunk'; 
allQ.centerpos{3}   = sBdy.clc.nearPos; 

% -----------------------------
% NEG REW
% REWARD SIZE
sBdy.clc.startRew   = -1;
[newQ ] = CalcQDirect(sBdy);
% store q values and attributes
allQ.qVals{4}       = newQ;
allQ.dir(4)         = -1; 
allQ.rew(4)         = sBdy.clc.startRew;
allQ.bodyPart{4}    = 'Trunk'; 
allQ.centerpos{4}   = sBdy.clc.nearPos; 


% #################################################################
% =========================================================
% STIM MOVING TOWARDS HEAD

sHed                    = s;
sHed.clc.nearPos        = s.clc.nearPos + [0 0 11]';
sHed.clc.actConsequence = sBdy.clc.actConsequence;
% Location and size of limb, and where the Q-values are calculated FROM
sHed.clc.startSR = []; sHed.clc.startSC = []; sHed.clc.startSZ = [];
hedRs                   = sHed.clc.nearPos(1) + [1 0 1];
hedCs                   = sHed.clc.nearPos(2) + [-1:1];
hedZs                   = sHed.clc.nearPos(3) + [-2:2];
iRew = 0; % Initialise rewarded block counter
for  iZ = 1:length(hedZs)
    for  iC = 1:length(hedCs)
        iRew = iRew + 1;
        sHed.clc.startSR(iRew) =  hedRs(iC);
        sHed.clc.startSC(iRew) =  hedCs(iC);
        sHed.clc.startSZ(iRew) =  hedZs(iZ);
    end
end

% -----------------------------
% POS REW
% REWARD SIZE
sHed.clc.startRew   = 1;
[newQ ]             = CalcQDirect(sHed);
% store q values and attributes
allQ.qVals{5}       = newQ;
allQ.dir(5)         = 1;
allQ.rew(5)         = sHed.clc.startRew;
allQ.bodyPart{5}    = 'Head'; 
allQ.centerpos{5}   = sBdy.clc.nearPos;


% -----------------------------
% NEG REW
% REWARD SIZE
sHed.clc.startRew   = -1;
[newQ ]             = CalcQDirect(sHed);
% store q values and attributes
allQ.qVals{6}       = newQ;
allQ.dir(6)         = 1; 
allQ.rew(6)         = sHed.clc.startRew;
allQ.bodyPart{6}    = 'Head'; 
allQ.centerpos{6}   = sBdy.clc.nearPos;


% =========================================================
% STIM MOVING AWAY FROM HEAD
sHed.clc.stimDynams =     @(pos) pos - s.clc.baseVel; % For receding, set speed negative

% -----------------------------
% POS REW
% REWARD SIZE
sHed.clc.startRew   = 1;
[newQ ]             = CalcQDirect(sHed);
% store q values and attributes
allQ.qVals{7}       = newQ;
allQ.dir(7)         = -1; 
allQ.rew(7)         = sHed.clc.startRew;
allQ.bodyPart{7}    = 'Head'; 
allQ.centerpos{7}   = sBdy.clc.nearPos;

% -----------------------------
% NEG REW
% REWARD SIZE
sHed.clc.startRew   =  -1;
[newQ ]             = CalcQDirect(sHed);
% store q values and attributes
allQ.qVals{8}       = newQ;
allQ.dir(8)         = -1; 
allQ.rew(8)         = sHed.clc.startRew;
allQ.bodyPart{8}    = 'Head'; 
allQ.centerpos{8}   = sBdy.clc.nearPos; 


%% #################################################################
% ==========================================================
% HAND by chest, STIMULUS MOVING TOWARDS HAND

sHnd = s;
sHnd.clc.nearPos = sHnd.clc.nearPos;
% Hand movementrs need to be allowed to be faster
sHnd.clc.actConsequence = [ ...
     0  0  0 ; ... % action 1  stay
     0  1  0 ; ... % action 2  left
     0 -1  0 ; ... % action 3  right
     0  0  1 ; ... % action 4  up
     0  0 -1 ; ... % action 5  down
    -1  0  0 ; ... % action 6  forward
     1  0  0 ; ... % action 7  backward
     0  2  0 ; ... % action 8  left
     0 -2  0 ; ... % action 9  right
     0  0  2 ; ... % action 10 up
     0  0 -2 ; ... % action 11 down
    -2  0  0 ; ... % action 12 forward
     2  0  0 ; ... % action 13 backward
    -3  0  0 ; ... % action 14 forward
     3  0  0 ; ... % action 15 backward
    -4  0  0 ; ... % action 16 forward 
     4  0  0 ; ... % action 17 backward
    -5  0  0 ; ... % action 18 forward
     5  0  0 ; ... % action 19 backward
     0  3  0 ; ... % action 20 left
     0 -3  0 ; ... % action 21 right
     0  0  3 ; ... % action 22 up
     0  0 -3 ; ... % action 23 down
     0  4  0 ; ... % action 24 left
     0 -4  0 ; ... % action 24 right
     0  0  4 ; ... % action 26 up
     0  0 -4 ; ... % action 27 down
     0  5  0 ; ... % action 28 left
     0 -5  0 ; ... % action 29 right
     0  0  5 ; ... % action 30 up
     0  0 -5 ; ... % action 31 down
     1  3  0 ; ... % action 32 left
     1 -3  0 ; ... % action 33 right
     1  0  3 ; ... % action 34 up
     1  0 -3 ; ... % action 35 down
     1  4  0 ; ... % action 36 left
     1 -4  0 ; ... % action 37 right
     1  0  4 ; ... % action 38 up
     1  0 -4 ; ... % action 39 down
     1  5  0 ; ... % action 40 left
     1 -5  0 ; ... % action 41 right
     1  0  5 ; ... % action 42 up
     1  0 -5 ];    % action 43 down


% Location and size of limb, and where the Q-values are calculated FROM
sHnd.clc.startSR = []; sHnd.clc.startSC = []; sHnd.clc.startSZ = [];
hndRs   = sHnd.clc.nearPos(1) + [0 0 0];
hndCs   = sHnd.clc.nearPos(2) + [-1:1];
hndZs   = sHnd.clc.nearPos(3) + [-1:1];
iRew    = 0; % Initialise rewarded block counter
for  iZ = 1:length(hndZs)
    for  iC = 1:length(hndCs)
        iRew = iRew + 1;
        sHnd.clc.startSR(iRew) =  hndRs(iC);
        sHnd.clc.startSC(iRew) =  hndCs(iC);
        sHnd.clc.startSZ(iRew) =  hndZs(iZ);
    end
end

%-----------------------------
% POS REW
% REWARD SIZE
sHnd.clc.startRew   = 1;
[newQ ]             = CalcQDirect(sHnd);
% store q values and attributes
allQ.qVals{9}       = newQ;
allQ.dir(9)         = 1;
allQ.rew(9)         = sHnd.clc.startRew;
allQ.bodyPart{9}    = 'HandByChest'; 
allQ.centerpos{9}   = sHnd.clc.nearPos; 

% -----------------------------
% NEG REW
% REWARD SIZE
sHnd.clc.startRew   = -1;
[newQ ]             = CalcQDirect(sHnd);
% store q values and attributes
allQ.qVals{10}      = newQ;
allQ.dir(10)        = 1;
allQ.rew(10)        = sHnd.clc.startRew;
allQ.bodyPart{10}   = 'HandByChest'; 
allQ.centerpos{10}  = sHnd.clc.nearPos;


% =========================================================
% STIMULUS MOVING AWAY FROM HAND
sHnd.clc.stimDynams =     @(pos) pos - s.clc.baseVel; % For receding, set speed negative

% -----------------------------
% POS REW
% REWARD SIZE
sHnd.clc.startRew   = 1;
[newQ ]             = CalcQDirect(sHnd);
% store q values and attributes
allQ.qVals{11}      = newQ;
allQ.dir(11)        = -1;
allQ.rew(11)        = sHnd.clc.startRew;
allQ.bodyPart{11}   = 'HandByChest'; 
allQ.centerpos{11}  = sHnd.clc.nearPos; 

% -----------------------------
% NEG REW
% REWARD SIZE
sHnd.clc.startRew   =  -1;
[newQ ]             = CalcQDirect(sHnd);
% store q values and attributes
allQ.qVals{12}      = newQ;
allQ.dir(12)        = -1; 
allQ.rew(12)        = sHnd.clc.startRew; 
allQ.bodyPart{12}   = 'HandByChest'; 
allQ.centerpos{12}  = sHnd.clc.nearPos;


%% #################################################################
% ==========================================================
% HAND by side, stimulus all directions 

sHndSide = sHnd;
sHndSide.clc.nearPos = sHnd.clc.nearPos + [0 10 0]';
sHndSide.clc.startSC = sHnd.clc.startSC + 10;

% Move the data, then shift the Q values next - don't need to recompute
allQ(13:16,:)           = allQ(9:12,:);
allQ.bodyPart(13:16)    = repmat({'HandBySide'},[1 4])';

handMidPos      = sHnd.clc.nearPos(2);
handNewPos      = sHndSide.clc.nearPos(2);
sideHandShift   = handNewPos - handMidPos;

% Find transposable width on the narrower side
moveWidth       = s.wrld.size(2) - handNewPos;

% Shift the data
for iPsi = 13:16
    allQ.qVals{iPsi} = zeros(size(allQ.qVals{iPsi}));

    % Flip the Q-vals because it's symetric anyway and the extent on the
    % right is bigger than on the left
    leftCopyDist    = min([handMidPos - 1 , size(    allQ.qVals{iPsi},3) - handNewPos]);
    rightCopyDist   = min([handNewPos - 1 , size(    allQ.qVals{iPsi},3) - handMidPos]);

    % Move the data from the right of the central hand to the left of the
    % lateral hand, and flip it
    allQ.qVals{iPsi}(:, :, handNewPos + [-rightCopyDist:0]   ,:) = ...
        flip(allQ.qVals{iPsi-4}(:, :, handMidPos + [0:rightCopyDist] ,:),3);

    % Move the data from the left of the central hand to the right of the
    % lateral hand, and flip it
    allQ.qVals{iPsi}(:, :, handNewPos + [1:leftCopyDist],:) = ...
        flip(allQ.qVals{iPsi-4}(:, :, handMidPos + [-leftCopyDist:-1] ,:),3);
    allQ.centerpos{iPsi} =  sHndSide.clc.nearPos;
        
end


%% #################################################################
% TOOL type 1: stick
% ==========================================================
% TOOL by HAND, STIMULUS MOVING TOWARDS

sHndTool = sHnd;
sHndTool.clc.stimDynams =  @(pos) pos + s.clc.baseVel;

% Location and size of Tool tip, and where the Q-values are calculated FROM
sHndToolTool.clc.toolSR = []; sHndTool.clc.toolSC = []; sHnd.clc.toolSZ = [];
tlRs = sHndTool.clc.nearPos(1) - 13; % 60 cm away
tlCs = sHndTool.clc.nearPos(2) + 0;
tlZs = sHndTool.clc.nearPos(3) + 0;
iRew = 0; % Initialise rewarded block counter
for  iC = 1:length(tlCs)
        iRew = iRew + 1;
        sHndTool.clc.toolSR(iRew) =  tlRs(iC);
        sHndTool.clc.toolSC(iRew) =  tlCs(iC);
        sHndTool.clc.toolSZ(iRew) =  tlZs(iC);
end

%-----------------------------
% POS REW
% REWARD SIZE
sHndTool.clc.stimDynams =     @(pos) pos + s.clc.baseVel;

sHndTool.clc.startRew   = 1;
[newQ ]                 = CalcQDirect(sHndTool);
% store q values and attributes
allQ.qVals{17}          = newQ;
allQ.dir(17)            = 1;
allQ.rew(17)            = sHndTool.clc.startRew;
allQ.bodyPart{17}       = 'HandByChestPlusTool'; 
allQ.centerpos{17}      = sHndTool.clc.nearPos;


% -----------------------------
% NEG REW
% REWARD SIZE
sHndTool.clc.stimDynams =     @(pos) pos + s.clc.baseVel;

sHndTool.clc.startRew   = -1;
[newQ ]                 = CalcQDirect(sHndTool);
% store q values and attributes
allQ.qVals{18}          = newQ;
allQ.dir(18)            = 1; 
allQ.rew(18)            = sHndTool.clc.startRew;
allQ.bodyPart{18}       = 'HandByChestPlusTool'; 
allQ.centerpos{18}      = sHndTool.clc.nearPos; 

% -----------------------------
% POS REW, AWAY
% REWARD SIZE
sHndTool.clc.stimDynams =     @(pos) pos - s.clc.baseVel;

sHndTool.clc.startRew   = 1;
[newQ ]                 = CalcQDirect(sHndTool);
% store q values and attributes
allQ.qVals{19}          = newQ;
allQ.dir(19)            = -1; 
allQ.rew(19)            = sHndTool.clc.startRew;
allQ.bodyPart{19}       = 'HandByChestPlusTool'; 
allQ.centerpos{19}      = sHndTool.clc.nearPos;


% -----------------------------
% NEG REW, AWAY
% REWARD SIZE
sHndTool.clc.stimDynams =     @(pos) pos - s.clc.baseVel;

sHndTool.clc.startRew   = -1;
[newQ ]                 = CalcQDirect(sHndTool);
% store q values and attributes
allQ.qVals{20}          = newQ;
allQ.dir(20)            = -1;
allQ.rew(20)            = sHndTool.clc.startRew;
allQ.bodyPart{20}       = 'HandByChestPlusTool'; 
allQ.centerpos{20}      = sHndTool.clc.nearPos;


%% #################################################################
% Spider Butterfly Track
% Legacy: before revisions, here we created one positive reward, and TWO 
% negative rewards, I think. One where the negative is equal reward to the 
% positive, and the other where negative is 1.25x the absolute reward.
% HOWEVER, after during the revision we found a better way: set the
% negative reward as a fitted variable later down the line, and ignore the
% 1.25x
% ==========================================================


% Create separate Q-values that follow the experiment's 'stimulus on a
% fixed track' dynamics
sTrack                      = sHnd;
sTrack.clc.nReps            = 1;
% Hand can't move for the track experiment
sTrack.clc.actConsequence   = [ 0  0  0 ];    % action 1 stay 

% Random stimulus dynamics
rSprT   = [0];
rSprPrT = [1];
cSprT   = [0];
cSprPrT = [1];
zSprT   = [0];
zSprPrT = [1];

% Change settings to inlcude this experiment's dynamics
[newPos,sTrack]         = TrackDynamics([0 0 0],sTrack);
sTrack.clc.stimDynams   = @(pos) TrackDynamics(pos,sTrack); % Use the track dynamics
sTrack.clc.randSpread   = {rSprT cSprT zSprT}; % Put a little biit of x and z variability in? Kind of arbitrary
sTrack.clc.spreadProb   = {rSprPrT cSprPrT  zSprPrT}; % x y z, probabilities of spread

%-----------------------------
% POS REW
sTrack.clc.startRew = 1;
[newQ ]             = CalcQDirect(sTrack);
% store q values and attributes
allQ.qVals{21}      = newQ;
allQ.dir(21)        = 1; 
allQ.rew(21)        = sTrack.clc.startRew;
allQ.bodyPart{21}   = 'HandByChestTrack'; 
allQ.centerpos{21}  = sTrack.clc.nearPos;


%-----------------------------
% NEG REW
sTrack.clc.startRew = -1;
[newQ ]             = CalcQDirect(sTrack);
% store q values and attributes
allQ.qVals{22}      = newQ;
allQ.dir(22)        = 1;
allQ.rew(22)        = sTrack.clc.startRew;
allQ.bodyPart{22}   = 'HandByChestTrack'; 
allQ.centerpos{22}  = sTrack.clc.nearPos; 

%-----------------------------
% NEG REW, bigger
sTrack.clc.startRew = -1.25;
[newQ ]             = CalcQDirect(sTrack);
% store q values and attributes
allQ.qVals{23}      = newQ;
allQ.dir(23)        = 1;
allQ.rew(23)        = sTrack.clc.startRew; 
allQ.bodyPart{23}   = 'HandByChestTrack'; 
allQ.centerpos{23}  = sTrack.clc.nearPos; 


% #################################################################
% Make extra q with reward -1.25 for the more threatening situation for body as well
allQ(24,:)          = allQ(2,:);
allQ.qVals{24}      = 1.25 .* allQ.qVals{2};
allQ.rew(24)        = 1.25 .* allQ.rew(24);
allQ(25,:)          = allQ(6,:);
allQ.qVals{25}      = 1.25 .* allQ.qVals{6};
allQ.rew(25)        = 1.25 .* allQ.rew(6);

% #################################################################
% Make extra q with NO VALUES because nothing is active
allQ(26,:)          = allQ(1,:);
allQ.qVals{26}      = zeros(size(allQ.qVals{1}));
allQ.rew(26)        = 0;
allQ.bodyPart{26}   = 'None'; 


%% #################################################################
% TOOL type 2: rake
% ==========================================================
% TOOL by HAND, TOWARDS

sHndToolRake = sHnd;
sHndToolRake.clc.stimDynams =  @(pos) pos + s.clc.baseVel;

% Location and size of Tool tip, and where the Q-values are calculated FROM
sHndToolRake.clc.toolSR = []; sHndToolRake.clc.toolSC = []; sHnd.clc.toolSZ = [];
tlRs = sHndToolRake.clc.nearPos(1) - [29 29 29 29 29]; % 145 cm away
tlCs = sHndToolRake.clc.nearPos(2) + [2  -1  0  1  2];
tlZs = sHndToolRake.clc.nearPos(3) + [0   0  0  0  0];
iRew = 0; % Initialise rewarded block counter
for  iC = 1:length(tlCs)
        iRew = iRew + 1;
        sHndToolRake.clc.toolSR(iRew) =  tlRs(iC);
        sHndToolRake.clc.toolSC(iRew) =  tlCs(iC);
        sHndToolRake.clc.toolSZ(iRew) =  tlZs(iC);
end

%-----------------------------
% POS REW
sHndToolRake.clc.stimDynams =     @(pos) pos + s.clc.baseVel;

sHndToolRake.clc.startRew   =  1;
[newQ ]                     = CalcQDirect(sHndToolRake);
% store q values and attributes
allQ.qVals{27}              = newQ;
allQ.dir(27)                = 1; 
allQ.rew(27)                = sHndToolRake.clc.startRew;
allQ.bodyPart{27}           = 'HandByChestPlusToolRake'; 
allQ.centerpos{27}          = sHndToolRake.clc.nearPos; 


% -----------------------------
% NEG REW
sHndToolRake.clc.stimDynams =     @(pos) pos + s.clc.baseVel;

sHndToolRake.clc.startRew   = -1;
[newQ ]                     = CalcQDirect(sHndToolRake);
% store q values and attributes
allQ.qVals{28}              = newQ;
allQ.dir(28)                = 1; 
allQ.rew(28)                = sHndToolRake.clc.startRew;
allQ.bodyPart{28}           = 'HandByChestPlusToolRake'; 
allQ.centerpos{28}          = sHndToolRake.clc.nearPos; 

%  -----------------------------
% POS REW, AWAY
sHndToolRake.clc.stimDynams =     @(pos) pos - s.clc.baseVel;

sHndToolRake.clc.startRew   = 1;
[newQ ]                     = CalcQDirect(sHndToolRake);
% store q values and attributes
allQ.qVals{29}              = newQ;
allQ.dir(29)                = -1; 
allQ.rew(29)                = sHndToolRake.clc.startRew; 
allQ.bodyPart{29}           = 'HandByChestPlusToolRake'; 
allQ.centerpos{29}          = sHndToolRake.clc.nearPos;


% -----------------------------
% NEG REW, AWAY
% REWARD SIZE
sHndToolRake.clc.stimDynams =     @(pos) pos - s.clc.baseVel;

sHndToolRake.clc.startRew   = -1;
[newQ ]                     = CalcQDirect(sHndToolRake);
% store q values and attributes
allQ.qVals{30}              = newQ;
allQ.dir(30)                = -1;
allQ.rew(30)                = sHndToolRake.clc.startRew; 
allQ.bodyPart{30}           = 'HandByChestPlusToolRake'; 
allQ.centerpos{30}          = sHndToolRake.clc.nearPos;



%% #################################################################
% Monkey arm

% ==========================================================
% ARM FORWARD ['right']
sBdy.clc.nReps = 2;
sArm = sBdy;
sArm.clc.stimDynams =  @(pos) pos + s.clc.baseVel;

% Location and size of limb, and where the Q-values are calculated FROM
% [treat arm like tool so that space behind it is not set to 1]
sArm.clc.toolSR = []; sArm.clc.toolSC = []; sArm.clc.toolSZ = [];
armRs = sArm.clc.nearPos(1) - [14 14 ...
                               13 13 ...
                               12 12 ...
                               11 11 ...
                               10 10 ...
                                9  9 ...
                                8  8 ...
                                7  7 ...
                                6  6 ...
                                5  5 ...
                                4  4 ...
                                3  3 ...
                                2  2 2 ...
                                1  1 ...
                                0  0 ];
armCs = sArm.clc.nearPos(2) + [ 5  6 ...
                                5  6 ...
                                5  6 ...
                                5  6 ...
                                5  6 ...
                                5  6 ...
                                5  6 ...
                                5  6 ...
                                5  6 ...
                                5  6 ...
                                5  6 ...
                                5  6 ...
                                4  5  6 ...
                                4  5 ...
                                4  5];
armZs = sArm.clc.nearPos(3) + [ 0  0 ...
                                0  0 ...
                                0  0 ...
                                0  0 ...
                                0  0 ...
                                0  0 ...
                                0  0 ...
                                0  0 ...
                                0  0 ...
                                0  0 ...
                                0  0 ...
                                1  1 ...
                                1  2  2 ...
                                2  2 ...
                                3  3];
iRew = 0; % Initialise rewarded block counter
for  iZ = 1:length(armZs)
        sArm.clc.toolSR(iZ) =  armRs(iZ);
        sArm.clc.toolSC(iZ) =  armCs(iZ);
        sArm.clc.toolSZ(iZ) =  armZs(iZ);
end

sArm.clc.startSR =  60;
sArm.clc.startSC =  sArm.clc.nearPos(2);
sArm.clc.startSZ =  sArm.clc.nearPos(3);

% -----------------------------
% POS REW
sArm.clc.stimDynams =     @(pos) pos + s.clc.baseVel;

sArm.clc.startRew   = 1;
[newQ ]             = CalcQDirect(sArm);
% store q values and attributes
allQ.qVals{31}      = newQ;
allQ.dir(31)        = 1; 
allQ.rew(31)        = sArm.clc.startRew;
allQ.bodyPart{31}   = 'ArmForward'; 
allQ.centerpos{31}  = sArm.clc.nearPos;


% -----------------------------
% NEG REW
% REWARD SIZE
sArm.clc.startRew   = -1;
[newQ ]             = CalcQDirect(sArm);
% store q values and attributes
allQ.qVals{32}      = newQ;
allQ.dir(32)        = 1;
allQ.rew(32)        = sArm.clc.startRew;
allQ.bodyPart{32}   = 'ArmForward'; 
allQ.centerpos{32}  = sArm.clc.nearPos;


% ==========================================================
% ARM LEFT
sArmL = sArm;

% Location and size of limb, and where the Q-values are calculated FROM
sArmL.clc.toolSR = []; sArmL.clc.toolSC = []; sArmL.clc.toolSZ = [];
armRs = sArmL.clc.nearPos(1) - [12 12 ...
                               11 11 11 ...
                               10 10 10 ...
                                9  9  9 ...
                                8  8  8 ...
                                7  7  7 ...
                                6  6  6 ...
                                5  5 ...
                                4  4  4 ...
                                3  3 ...
                                2  2 2 ...
                                1  1 ...
                                0  0 ];
armCs = sArmL.clc.nearPos(2) +  [-6,-5, ...
                                -6,-5,-4, ...
                                -5,-4,-3, ...
                                -4,-3,-2, ...
                                -3,-2,-1, ...
                                -2,-1,-0, ...
                                -1, 0, 1, ...
                                 0, 1, ...
                                 0, 1, 2, ...
                                 1, 2, ...
                                 1, 2, 3, ...
                                 2, 3, ...
                                 2, 3];
armZs = sArmL.clc.nearPos(3) + [ 0  0 ...
                                0  0  0 ...
                                0  0  0 ...
                                0  0  0 ...
                                0  0  0 ...
                                0  0  0 ...
                                0  0  0 ...
                                0  0 ...
                                0  0  0 ...
                                1  1 ...
                                2  2  2 ...
                                2  2 ...
                                3  3];
iRew = 0; % Initialise rewarded block counter
for  iZ = 1:length(armZs)
        sArmL.clc.toolSR(iZ) =  armRs(iZ);
        sArmL.clc.toolSC(iZ) =  armCs(iZ);
        sArmL.clc.toolSZ(iZ) =  armZs(iZ);
end


%-----------------------------
% POS REW
sArmL.clc.stimDynams =     @(pos) pos + s.clc.baseVel;

sArmL.clc.startRew  = 1;
[newQ ]             = CalcQDirect(sArmL);
% store q values and attributes
allQ.qVals{33}      = newQ;
allQ.dir(33)        = 1; 
allQ.rew(33)        = sArmL.clc.startRew;
allQ.bodyPart{33}   = 'ArmLeft'; 
allQ.centerpos{33}  = sArmL.clc.nearPos;


% -----------------------------
% NEG REW
sArmL.clc.startRew  = -1;
[newQ ]             = CalcQDirect(sArmL);
% store q values and attributes
allQ.qVals{34}      = newQ;
allQ.dir(34)        = 1; 
allQ.rew(34)        = sArmL.clc.startRew;
allQ.bodyPart{34}   = 'ArmLeft'; 
allQ.centerpos{34}  = sArmL.clc.nearPos;



%% #################################################################
% =========================================================
% MONKEY HEAD - CONSTRAINED FORWARD

sHedConstr = sHed;
sHedConstr.clc.actConsequence = [ 0  0  0 ] ; % action 1 stay ONLY

sHedConstr.clc.stimDynams =  @(pos) pos + sHedConstr.clc.baseVel;

% -----------------------------
% POS REW
sHedConstr.clc.startRew =  1;
[newQ ]                 = CalcQDirect(sHedConstr);
% store q values and attributes
allQ.qVals{35}          = newQ;
allQ.dir(35)            = 1; 
allQ.rew(35)            = sHedConstr.clc.startRew;
allQ.bodyPart{35}       = 'HeadConstr'; 
allQ.centerpos{35}      = sHedConstr.clc.nearPos; 

% -----------------------------
% NEG REW
sHedConstr.clc.startRew =  -1;
[newQ ]                 = CalcQDirect(sHedConstr);
% store q values and attributes
allQ.qVals{36}          = newQ;
allQ.dir(36)            = 1;
allQ.rew(36)            = sHedConstr.clc.startRew; 
allQ.bodyPart{36}       = 'HeadConstr'; 
allQ.centerpos{36}      = sHedConstr.clc.nearPos; 


%% #################################################################
% =========================================================
% MONKEY HEAD 15 degree angle
sHed15 = sHedConstr;

% Move the data, then shift the Q values next
allQ(37:38,:)           = allQ(35:36,:);
allQ.bodyPart(37:38)    = repmat({'RotatedHead'},[1 2])';

% Shift the data
for iPsi = 37:38
    % shift the q values to the side by 15 degree, i.e. once for every 4
    % rows, starting at the face
    shiftAmount = 1;
    for iR = (sHed15.clc.nearPos(1) + 2) : -1 : 1
        allQ.qVals{iPsi}(:,iR, 1:end-shiftAmount ,:) = allQ.qVals{iPsi-2}(:,iR, (shiftAmount+1):end ,:);
        if mod(iR - (sHed15.clc.nearPos(1) + 3),4) == 0
            shiftAmount = shiftAmount + 1;
        end
    end
end


%% #################################################################
% =========================================================
% Hit probability - body
% -----------------------------
% NOTE: deltaT is 0.5, so I have to halve all the velocities and spreads?

sBdyHP = sBdy;
sBdyHP.clc.actConsequence = [0 0 0]; 
sBdyHP.clc.gammaVal = 1;
sBdyHP.clc.startRew = 1;

% Random stimulus dynamics
rSpr    = [-1 0 1];
rSprPr  = gaussmf(rSpr,[1 0]) ./ sum(gaussmf(rSpr,[1 0]));
cSpr    = [-1 0 1];
cSprPr  = gaussmf(cSpr,[1 0]) ./ sum(gaussmf(cSpr,[1 0]));
zSpr    = [-1 0 1];
zSprPr  = gaussmf(zSpr,[1 0]) ./ sum(gaussmf(zSpr,[1 0]));

% divide probabilities over 2 because deltaT = 1/2
rSprPr = rSprPr./([2 1 2]) + rSprPr(1).*([0 1 0]); 
cSprPx = cSprPr./([2 1 2]) + cSprPr(1).*([0 1 0]); 
zSprPr = zSprPr./([2 1 2]) + zSprPr(1).*([0 1 0]); 

% x y z, Deterministic stimulus dynamics
sBdyHP.clc.stimDynams =     @(pos) pos + s.clc.baseVel./2; % For approaching, set speed positive
sBdyHP.clc.randSpread =     {rSpr cSpr zSpr}; % Put a little biit of x and z variability in (kind of arbitrary)
sBdyHP.clc.spreadProb =     {rSprPr cSprPr  zSprPr};

% Sensory stimulus spreading: from Straka et al (note: 1/2 because their 
% units were seconds, and ours are half-seconds
rSnsSpr = [-7:7]; 
cSnsSpr = [-2:2]; 
zSnsSpr = [-2:2]; 

% Set uncertainties to Straka Noel Hoffmann values
rSigma      = sqrt(2.5.^2 + 30.^2) ./ sBdyHP.clc.binW ./2;
rSnsSprPr   = gaussmf(rSnsSpr,[rSigma 0]) ./ sum(gaussmf(rSnsSpr,[rSigma 0]));
czSigma     = sqrt(5.^2 + 5.^2) ./ sBdyHP.clc.binW ./2;
cSnsSprPr   = gaussmf(cSnsSpr,[czSigma 0]) ./ sum(gaussmf(cSnsSpr,[czSigma 0]));
zSnsSprPr   = gaussmf(zSnsSpr,[czSigma 0]) ./ sum(gaussmf(zSnsSpr,[czSigma 0]));

sBdyHP.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sBdyHP.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

sBdyHP.clc.baseVel    = [2 0 0]; % Same velocity as Q values, but now in units of half-seconds
sBdyHP.clc.stimDynams = @(pos) pos + sBdyHP.clc.baseVel;

sBdyHP.clc.nReps        = 1; 
sBdyHP.clc.stepUpdateFl = 1; % whether to run a full sweep update or not
sBdyHP.clc.nSteps       = 1;

% Set the importance of false positives and false negatives
sForUtility.clc.FN = 5;
sForUtility.clc.FP = 1;

bodyHP                  = calcQDirect(sBdyHP);
cQ                      = 39;
allQ.qVals{cQ}          = HitProbToUtility(bodyHP,sForUtility);
allQ.dir(cQ)            = 1; 
allQ.rew(cQ)            = sBdyHP.clc.startRew; 
allQ.bodyPart{cQ}       = 'TrunkHP'; 
allQ.centerpos{cQ}      = sBdyHP.clc.nearPos; 

sBdyHP.clc.baseVel      = [-2 0 0]; 
sBdyHP.clc.stimDynams   = @(pos) pos + sBdyHP.clc.baseVel;
bodyHPaway              = calcQDirect(sBdyHP);
cQ                      = 40;
allQ.qVals{cQ}          = HitProbToUtility(bodyHPaway,sForUtility);
allQ.dir(cQ)            = -1;
allQ.rew(cQ)            = sBdyHP.clc.startRew;
allQ.bodyPart{cQ}       = 'TrunkHP'; 
allQ.centerpos{cQ}      = sBdyHP.clc.nearPos; 


% =========================================================
% Make Head/face hitprob
% -----------------------------
sHedHP = sHed;
% Variables to copy over
varsForHP = {'actConsequence','gammaVal','startRew','stimDynams',...
             'randSpread','spreadProb','sensSpread','sensProb',...
             'baseVel','stimDynams','nReps','stepUpdateFl','nSteps'};
for iV = 1:numel(varsForHP)
    sHedHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
headHP                  = calcQDirect(sHedHP);
cQ                      = 41;
allQ.qVals{cQ}          = HitProbToUtility(headHP,sForUtility);
allQ.dir(cQ)            = 1;
allQ.rew(cQ)            = sHedHP.clc.startRew;
allQ.bodyPart{cQ}       = 'HeadHP'; 
allQ.centerpos{cQ}      = sHedHP.clc.nearPos; 


sHedHP.clc.baseVel      = [-2 0 0]; 
sHedHP.clc.stimDynams   = @(pos) pos + sHedHP.clc.baseVel;
hedHPaway               = calcQDirect(sHedHP);
cQ                      = 42;
allQ.qVals{cQ}          = HitProbToUtility(headHPaway,sForUtility);;
allQ.dir(cQ)            = -1;
allQ.rew(cQ)            = sHedHP.clc.startRew;
allQ.bodyPart{cQ}       = 'HeadHP'; 
allQ.centerpos{cQ}      = sHedHP.clc.nearPos; 


% =========================================================
% Make hand hitprob
% -----------------------------
sHndHP = sHnd;
% Variables to copy over
for iV = 1:numel(varsForHP)
    sHndHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
handHP                  = calcQDirect(sHndHP);
cQ                      = 43;
allQ.qVals{cQ}          = HitProbToUtility(handHP,sForUtility);
allQ.dir(cQ)            = 1;
allQ.rew(cQ)            = sHndHP.clc.startRew;
allQ.bodyPart{cQ}       = 'HandByChestHP'; 
allQ.centerpos{cQ}      = sHndHP.clc.nearPos; 


sHndHP.clc.baseVel      = [-2 0 0];
sHndHP.clc.stimDynams   = @(pos) pos + sHndHP.clc.baseVel;
handHPaway              = calcQDirect(sHndHP);
cQ                      = 44;
allQ.qVals{cQ}          = HitProbToUtility(handHPaway,sForUtility);
allQ.dir(cQ)            = -1;
allQ.rew(cQ)            = sHndHP.clc.startRew;
allQ.bodyPart{cQ}       = 'HandByChestHP'; 
allQ.centerpos{cQ}      = sHndHP.clc.nearPos; 


% =========================================================
% ------------------------------------------------
% HAND BY SIDE HITPROB
sHndSideHP              = sHndHP;
sHndSideHP.clc.nearPos  = sHndSide.clc.nearPos;
sHndSideHP.clc.startSC  = sHndSide.clc.startSC;

% Move the data, then shift the Q values next
allQ(45:46,:)           = allQ(43:44,:);
allQ.bodyPart(45:46)    = repmat({'HandBySideHP'},[1 2])';

% Shift the data
for iPsi = 45:46
    allQ.qVals{iPsi}    = zeros(size(allQ.qVals{iPsi}));

    % Flip the data because it's symetric anyway and the extent on the
    % right is bigger than on the left
    leftCopyDist  = min([handMidPos - 1 , size(    allQ.qVals{iPsi},3) - handNewPos]);
    rightCopyDist = min([handNewPos - 1 , size(    allQ.qVals{iPsi},3) - handMidPos]);
    

    % Move the data from the right of the central hand to the left of the
    % lateral hand, and flip it
    allQ.qVals{iPsi}(:, :, handNewPos + [-rightCopyDist:0]   ,:) = ...
        flip(allQ.qVals{iPsi-2}(:, :, handMidPos + [0:rightCopyDist] ,:),3);

    % Move the data from the left of the central hand to the right of the
    % lateral hand, and flip it
    allQ.qVals{iPsi}(:, :, handNewPos + [1:leftCopyDist],:) = ...
        flip(allQ.qVals{iPsi-2}(:, :, handMidPos + [-leftCopyDist:-1] ,:),3);
    allQ.centerpos{iPsi}   = sHnd.clc.nearPos;

    allQ.centerpos{iPsi} =  sHndSide.clc.nearPos;   
end
% ------------------------------------------------



%% #################################################################
% Make Hand tool hitprobs 
% -----------------------------
sHndToolHP = sHndTool;
for iV = 1:numel(varsForHP)
    sHndToolHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
handToolHP                  = calcQDirect(sHndToolHP);
cQ                          = 47;
allQ.qVals{cQ}              = HitProbToUtility(handToolHP,sForUtility);
allQ.dir(cQ)                = 1;
allQ.rew(cQ)                = sHndToolHP.clc.startRew;
allQ.bodyPart{cQ}           = 'HandByChestPlusToolHP'; 
allQ.centerpos{cQ}          = sHndToolHP.clc.nearPos; 


sHndToolHP.clc.baseVel      = [-2 0 0];
sHndToolHP.clc.stimDynams   = @(pos) pos + sHndToolHP.clc.baseVel;
handToolHPaway              = calcQDirect(sHndToolHP);
cQ                          = 48;
allQ.qVals{cQ}              = HitProbToUtility(handToolHPaway,sForUtility);
allQ.dir(cQ)                = -1;
allQ.rew(cQ)                = sHndToolHP.clc.startRew;
allQ.bodyPart{cQ}           = 'HandByChestPlusToolHP'; 
allQ.centerpos{cQ}          = sHndToolHP.clc.nearPos; 


sHndToolRakeHP              = sHndToolRake;
for iV = 1:numel(varsForHP)
    sHndToolRakeHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
handToolRakeHP              = calcQDirect(sHndToolHP);
cQ                          = 49;
allQ.qVals{cQ}              = HitProbToUtility(handToolRakeHP,sForUtility);
allQ.dir(cQ)                = 1;
allQ.rew(cQ)                = sHndToolRakeHP.clc.startRew;
allQ.bodyPart{cQ}           = 'HandByChestPlusToolRakeHP'; 
allQ.centerpos{cQ}          = sHndToolRakeHP.clc.nearPos; 


sHndToolRakeHP.clc.baseVel      = [-2 0 0];
sHndToolRakeHP.clc.stimDynams   = @(pos) pos + sHndToolRakeHP.clc.baseVel;
handToolRakeHPaway              = calcQDirect(sHndToolRakeHP);
cQ                      = 50;
allQ.qVals{cQ}          = HitProbToUtility(handToolRakeHPaway,sForUtility);
allQ.dir(cQ)            = -1;
allQ.rew(cQ)            = sHndToolRakeHP.clc.startRew;
allQ.bodyPart{cQ}       = 'HandByChestPlusToolRakeHP'; 
allQ.centerpos{cQ}      = sHndToolRakeHP.clc.nearPos; 


%% #################################################################
% Make Hand track hitprob
% -----------------------------
sTrackHP            = sTrack;
% Variables to copy over
varsForTrackHP      = {'actConsequence','gammaVal','startRew',...
                       'sensSpread','sensProb',...
                       'baseVel','nReps','stepUpdateFl','nSteps'};
for iV = 1:numel(varsForTrackHP)
    sTrackHP.clc.(varsForTrackHP{iV}) = sBdyHP.clc.(varsForTrackHP{iV});
end
trackHP             = calcQDirect(sTrackHP);
cQ                  = 51;
allQ.qVals{cQ}      = trackHP;
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sTrackHP.clc.startRew; 
allQ.bodyPart{cQ}   = 'HandByChestTrackHP'; 
allQ.centerpos{cQ}  = sTrackHP.clc.nearPos; 


%% #################################################################
%  Make Monkey hitprobs 
%  -----------------------------
sArmHP              = sArm;
% For some reason, I made the monkey's arm a tool, but I can't remember why. 
% Oh well, I should take that into account
sArmHP.clc.startSC  = sArmHP.clc.toolSC;
sArmHP.clc.startSR  = sArmHP.clc.toolSR;
sArmHP.clc.startSZ  = sArmHP.clc.toolSZ;
for iV = 1:numel(varsForHP)
    sArmHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
cQ                  = 52;
allQ.qVals{cQ}      = calcQDirect(sArmHP);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sArmHP.clc.startRew; 
allQ.bodyPart{cQ}   = 'ArmForwardHP'; 
allQ.centerpos{cQ}  = sArmHP.clc.nearPos; 


sArmLHP             = sArmL;
% For some reason, I made the monkey's arm a tool, but I can't remember why. 
% Oh well, I should take that into account
sArmLHP.clc.startSC = sArmLHP.clc.toolSC;
sArmLHP.clc.startSR = sArmLHP.clc.toolSR;
sArmLHP.clc.startSZ = sArmLHP.clc.toolSZ;
for iV = 1:numel(varsForHP)
    sArmLHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
cQ                  = 53;
allQ.qVals{cQ}      = calcQDirect(sArmLHP);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sArmLHP.clc.startRew; 
allQ.bodyPart{cQ}   = 'ArmLeftHP'; 
allQ.centerpos{cQ}  = sArmLHP.clc.nearPos; 


sHedConstrHP        = sHedConstr;
for iV = 1:numel(varsForHP)
    sHedConstrHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
cQ                  = 54;
allQ.qVals{cQ}      = calcQDirect(sHedConstrHP);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sHedConstrHP.clc.startRew; 
allQ.bodyPart{cQ}   = 'HeadConstrHP'; 
allQ.centerpos{cQ}  = sHedConstrHP.clc.nearPos; 


sHed15HP            = sHed15;
for iV = 1:numel(varsForHP)
    sHed15HP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
cQ                  = 55;
allQ.qVals{cQ}      = calcQDirect(sHed15HP);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sHed15HP.clc.startRew; 
allQ.bodyPart{cQ}   = 'RotatedHeadHP'; 
allQ.centerpos{cQ}  = sHed15HP.clc.nearPos; 


%% #################################################################
% Make multisensory integration model data
% Parameters are based of Bertoni et al 2021

% -----------------------------
sBdyMI = sBdyHP;

% Stimulus should not move at all
rSpr    = 0;
rSprPr  = 1;
cSpr    = 0;
cSprPr  = 1;
zSpr    = 0;
zSprPr  = 1;

% x y z, Deterministic stimulus dynamics
sBdyMI.clc.stimDynams =     @(pos) pos; % No speed
sBdyMI.clc.randSpread =     {rSpr cSpr zSpr}; % Put a little biit of x and z variability in? Kind of arbitrary
sBdyMI.clc.spreadProb =     {rSprPr cSprPr  zSprPr}; %

% Sensory uncertainty
rSnsSprMI = [-7:7]; 
cSnsSprMI = [-2:2]; 
zSnsSprMI = [-2:2]; 

% Set uncertainties to Bertoni 2021 or Noel 2018 values
% Bertoni 2021: 13cm for proprioception and 11cm for visual
% Noel 2018: 36 mm --> that's too small, it will never fit any of the real data
rSigmaMI      = sqrt(13.^2 + 11.^2) ./ sBdyMI.clc.binW ;
rSnsSprMIPr   = gaussmf(rSnsSprMI,[rSigmaMI 0]) ./ sum(gaussmf(rSnsSprMI,[rSigmaMI 0]));
czSigmaMI     = sqrt(13.^2 + 11.^2) ./ sBdyMI.clc.binW ./2;
cSnsSprMIPr   = gaussmf(cSnsSprMI,[czSigmaMI 0]) ./ sum(gaussmf(cSnsSprMI,[czSigmaMI 0]));
zSnsSprMIPr   = gaussmf(zSnsSprMI,[czSigmaMI 0]) ./ sum(gaussmf(zSnsSprMI,[czSigmaMI 0]));

sBdyMI.clc.sensSpread   = {rSnsSprMI cSnsSprMI zSnsSprMI};
sBdyMI.clc.sensProb     = {rSnsSprMIPr cSnsSprMIPr zSnsSprMIPr};

sBdyMI.clc.baseVel      = [0 0 0]; % Same velocity as Q values, but now in units of half-seconds
sBdyMI.clc.stimDynams   = @(pos) pos + sBdyMI.clc.baseVel;

sBdyMI.clc.nReps        = 1;
sBdyMI.clc.stepUpdateFl = 1; % whether to run a full sweep update or not
sBdyMI.clc.nSteps       = 1;

cQ = 56;
allQ.qVals{cQ}       = calcQDirect(sBdyMI);
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sBdyMI.clc.startRew;
allQ.bodyPart{cQ}    = 'TrunkMI'; 
allQ.centerpos{cQ}   = sBdyMI.clc.nearPos;

% =========================================================
% Make Head/face multsensint
% -----------------------------
sHedMI              = sHedHP;
varsForMI = {'stimDynams',...
             'randSpread','spreadProb','sensSpread','sensProb',...
             'baseVel','stimDynams','nReps','stepUpdateFl','nSteps'}
for iV = 1:numel(varsForMI)
    sHedMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ                  = 57;
allQ.qVals{cQ}      = calcQDirect(sHedMI);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sHedMI.clc.startRew;
allQ.bodyPart{cQ}   = 'HeadMI'; 
allQ.centerpos{cQ}  = sHedMI.clc.nearPos; 

% =========================================================
% Make hand multsensint
% -----------------------------
sHndMI              = sHndHP;
for iV = 1:numel(varsForMI)
    sHndMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ                  = 58;
allQ.qVals{cQ}      = calcQDirect(sHndMI);
allQ.dir(cQ)        = 1;
allQ.rew(cQ)        = sHndMI.clc.startRew;
allQ.bodyPart{cQ}   = 'HandByChestMI';
allQ.centerpos{cQ}  = sHndMI.clc.nearPos; 

% =========================================================
% ------------------------------------------------
% Hand by side multisens integration 
sHndSideMI              = sHndSideHP;
sHndSideMI.clc.nearPos  = sHndSideHP.clc.nearPos;
sHndSideMI.clc.startSC  = sHndSideHP.clc.startSC;


% Move the data, then shift the Q values next
allQ(59,:)          = allQ(58,:);
allQ.bodyPart(59)   = repmat({'HandBySideMI'},[1 1])';

% Shift the data
for iPsi = 59
    allQ.qVals{iPsi} = zeros(size(allQ.qVals{iPsi}));

    % FLip the data because it's symetric anyway and the extent on the
    % right is bigger than on the left
    leftCopyDist  = min([handMidPos - 1 , size(    allQ.qVals{iPsi},3) - handNewPos]);
    rightCopyDist = min([handNewPos - 1 , size(    allQ.qVals{iPsi},3) - handMidPos]);
    

    % Move the data from the right of the central hand to the left of the
    % lateral hand, and flip it
    allQ.qVals{iPsi}(:, :, handNewPos + [-rightCopyDist:0]   ,:) = ...
        flip(allQ.qVals{iPsi-1}(:, :, handMidPos + [0:rightCopyDist] ,:),3);

    % Move the data from the left of the central hand to the right of the
    % lateral hand, and flip it
    allQ.qVals{iPsi}(:, :, handNewPos + [1:leftCopyDist],:) = ...
        flip(allQ.qVals{iPsi-1}(:, :, handMidPos + [-leftCopyDist:-1] ,:),3);
    allQ.centerpos{iPsi}   = sHnd.clc.nearPos;

    allQ.centerpos{iPsi} =  sHndSide.clc.nearPos;   
end
% ------------------------------------------------


% =========================================================
% Make Hand tool multsensints 
% -----------------------------
sHndToolMI          = sHndToolHP;
for iV = 1:numel(varsForMI)
    sHndToolMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ                  = 60;
allQ.qVals{cQ}      = calcQDirect(sHndToolMI);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sHndToolMI.clc.startRew;
allQ.bodyPart{cQ}   = 'HandByChestPlusToolMI';
allQ.centerpos{cQ}  = sHndToolMI.clc.nearPos; 

% -----------------------------
sHndToolRakeMI      = sHndToolRakeHP;
for iV = 1:numel(varsForMI)
    sHndToolRakeMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ                  = 61;
allQ.qVals{cQ}      = calcQDirect(sHndToolRakeMI);
allQ.dir(cQ)        = 1;
allQ.rew(cQ)        = sHndToolRakeMI.clc.startRew;
allQ.bodyPart{cQ}   = 'HandByChestPlusToolRakeMI';
allQ.centerpos{cQ}  = sHndToolRakeMI.clc.nearPos;


% =========================================================
% Make Hand track multsensint
% -----------------------------
sTrackMI            = sTrackHP;
for iV = 1:numel(varsForMI)
    sTrackMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ                  = 62;
allQ.qVals{cQ}      = calcQDirect(sTrackMI);
allQ.dir(cQ)        = 1;
allQ.rew(cQ)        = sTrackMI.clc.startRew;
allQ.bodyPart{cQ}   = 'HandByChestTrackMI';
allQ.centerpos{cQ}  = sTrackMI.clc.nearPos;

% =========================================================
% Make Monkey multsensint 
% -----------------------------
sArmMI              = sArmHP;
for iV = 1:numel(varsForMI)
    sArmMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ                  = 63;
allQ.qVals{cQ}      = calcQDirect(sArmMI);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sArmMI.clc.startRew; 
allQ.bodyPart{cQ}   = 'ArmForwardMI'; 
allQ.centerpos{cQ}  = sArmMI.clc.nearPos; 


sArmLMI             = sArmLHP;
for iV = 1:numel(varsForMI)
    sArmLMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ                  = 64; 
allQ.qVals{cQ}      = calcQDirect(sArmLMI);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sArmLMI.clc.startRew; 
allQ.bodyPart{cQ}   = 'ArmLeftMI'; 
allQ.centerpos{cQ}  = sArmLMI.clc.nearPos; 


sHedConstrMI        = sHedConstrHP;
for iV = 1:numel(varsForMI)
    sHedConstrMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ                  = 65; 
allQ.qVals{cQ}      = calcQDirect(sHedConstrMI);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sHedConstrMI.clc.startRew; 
allQ.bodyPart{cQ}   = 'HeadConstrMI'; 
allQ.centerpos{cQ}  = sHedConstrMI.clc.nearPos; 


sHed15MI            = sHed15HP;
for iV = 1:numel(varsForMI)
    sHed15MI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ                  = 66;
allQ.qVals{cQ}      = calcQDirect(sHed15MI);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sHed15MI.clc.startRew; 
allQ.bodyPart{cQ}   = 'RotatedHeadMI'; 
allQ.centerpos{cQ}  = sHed15MI.clc.nearPos; 


%% #################################################################
%  Make a distance metric, so that I can use it to create exponential and sigmoid functions


% =========================================================
% Make Body distance
% -----------------------------
cQ                  = 67;
allQ.qVals{cQ}      = permute(CalcDistToBody(sBdyMI),[4 1 2 3]);
allQ.dir(cQ)        = 1;
allQ.rew(cQ)        = sBdyMI.clc.startRew;
allQ.bodyPart{cQ}   = 'TrunkDist'; 
allQ.centerpos{cQ}  = sBdyMI.clc.nearPos;

% =========================================================
% Make Head distance
% -----------------------------
cQ                  = 68;
allQ.qVals{cQ}      = permute(CalcDistToBody(sHedMI),[4 1 2 3]);
allQ.dir(cQ)        = 1;
allQ.rew(cQ)        = sHedMI.clc.startRew;
allQ.bodyPart{cQ}   = 'HeadDist'; 
allQ.centerpos{cQ}  = sHedMI.clc.nearPos;

% =========================================================
% Make Hand Chest distance
% -----------------------------
cQ                  = 69;
allQ.qVals{cQ}      = permute(CalcDistToBody(sHndMI),[4 1 2 3]);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sHndMI.clc.startRew;
allQ.bodyPart{cQ}   = 'HandByChestDist';
allQ.centerpos{cQ}  = sHndMI.clc.nearPos; 

% =========================================================
% Make Hand Side distance
% -----------------------------
cQ                      = 70;
sHndSideMI.clc.startSC  = sHndMI.clc.startSC + (sHndSideMI.clc.nearPos(2) - sHndMI.clc.nearPos(2));
allQ.qVals{cQ}          = permute(CalcDistToBody(sHndSideMI),[4 1 2 3]);
allQ.dir(cQ)            = 1;
allQ.rew(cQ)            = sHndSideMI.clc.startRew;
allQ.bodyPart{cQ}       = 'HandBySideDist';
allQ.centerpos{cQ}      = sHndSideMI.clc.nearPos; 

% =========================================================
% Make Hand tool distance
% -----------------------------
cQ                          = 71;
% First include the tool locations in the locations from which to calculate distance
sHndToolDist                = sHndToolMI; 
sHndToolDist.clc.startSC    = [sHndToolMI.clc.startSC , sHndToolMI.clc.toolSC];
sHndToolDist.clc.startSR    = [sHndToolMI.clc.startSR , sHndToolMI.clc.toolSR];
sHndToolDist.clc.startSZ    = [sHndToolMI.clc.startSZ , sHndToolMI.clc.toolSZ];
% Then continue as usual
allQ.qVals{cQ}       = permute(CalcDistToBody(sHndToolDist),[4 1 2 3]);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sHndToolDist.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolDist';
allQ.centerpos{cQ}   = sHndToolDist.clc.nearPos; 


% First include the tool locations in the locations from which to calculate distance
sHndToolRakeDist             = sHndToolRakeMI; 
sHndToolRakeDist.clc.startSC = [sHndToolRakeMI.clc.startSC , sHndToolRakeMI.clc.toolSC];
sHndToolRakeDist.clc.startSR = [sHndToolRakeMI.clc.startSR , sHndToolRakeMI.clc.toolSR];
sHndToolRakeDist.clc.startSZ = [sHndToolRakeMI.clc.startSZ , sHndToolRakeMI.clc.toolSZ];
% Then continue as usual
cQ                  = 72;
allQ.qVals{cQ}      = permute(CalcDistToBody(sHndToolRakeDist),[4 1 2 3]);
allQ.dir(cQ)        = 1;
allQ.rew(cQ)        = sHndToolRakeDist.clc.startRew;
allQ.bodyPart{cQ}   = 'HandByChestPlusToolRakeDist';
allQ.centerpos{cQ}  = sHndToolRakeDist.clc.nearPos; 

% =========================================================
% Make Hand track distance
% -----------------------------
cQ                  = 73;
allQ.qVals{cQ}      = permute(CalcDistToBody(sTrackMI),[4 1 2 3]);
allQ.dir(cQ)        = 1;
allQ.rew(cQ)        = sTrackMI.clc.startRew;
allQ.bodyPart{cQ}   = 'HandByChestTrackDist';
allQ.centerpos{cQ}  = sTrackMI.clc.nearPos;

% =========================================================
% Make Monkey distance
% -----------------------------
cQ                  = 74;
allQ.qVals{cQ}      = permute(CalcDistToBody(sArmMI),[4 1 2 3]);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sArmMI.clc.startRew; 
allQ.bodyPart{cQ}   = 'ArmForwardDist'; 
allQ.centerpos{cQ}  = sArmMI.clc.nearPos; 

cQ                  = 75; 
allQ.qVals{cQ}      = permute(CalcDistToBody(sArmLMI),[4 1 2 3]);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sArmLMI.clc.startRew; 
allQ.bodyPart{cQ}   = 'ArmLeftDist'; 
allQ.centerpos{cQ}  = sArmLMI.clc.nearPos; 

cQ                  = 76; 
allQ.qVals{cQ}      = permute(CalcDistToBody(sHedConstrMI),[4 1 2 3]);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sHedConstrMI.clc.startRew; 
allQ.bodyPart{cQ}   = 'HeadConstrDist'; 
allQ.centerpos{cQ}  = sHedConstrMI.clc.nearPos; 

cQ                  = 77;
allQ.qVals{cQ}      = permute(CalcDistToBody(sHed15MI),[4 1 2 3]);
allQ.dir(cQ)        = 1; 
allQ.rew(cQ)        = sHed15MI.clc.startRew; 
allQ.bodyPart{cQ}   = 'RotatedHeadDist'; 
allQ.centerpos{cQ}  = sHed15MI.clc.nearPos;




%% ------------------------------------------------------------------------
% #################################################################
%  Create Q-values using SARSA, instead of Q-learning

% Sensory stimulus spreading
rSnsSpr = [0]; 
cSnsSpr = [0]; 
zSnsSpr = [0]; 
% Set uncertainties 
rSigma      = sqrt(5.^2 + 5.^2) ./ sBdyHP.clc.binW ./2;
rSnsSprPr   = 1;
czSigma     = sqrt(5.^2 + 5.^2) ./ sBdyHP.clc.binW ./2;
cSnsSprPr   = 1;
zSnsSprPr   = 1;

% =========================================================
% BODY TOWARDS
sBdy.lp.alg = 'SARSA';
sBdy.clc.stimDynams     = @(pos) pos + sBdy.clc.baseVel;
sBdy.clc.sensSpread     = {rSnsSpr cSnsSpr zSnsSpr};
sBdy.clc.sensProb       = {rSnsSprPr cSnsSprPr zSnsSprPr};
sBdy.clc.nReps          = 2; 
sBdy.clc.stepUpdateFl   = 0; % whether to run a full sweep update or not
sBdy.clc.nSteps         = 1;

% -----------------------------
% POS REW
sBdy.clc.startRew   =  1;
[newQ ] = CalcQDirect(sBdy);
% store q values and attributes
cQ = 110;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sBdy.clc.startRew;
allQ.bodyPart{cQ}    = 'TrunkSARSA'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos; 

% -----------------------------
% NEG REW
sBdy.clc.startRew =  -1;
[newQ ] = CalcQDirect(sBdy);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sBdy.clc.startRew; 
allQ.bodyPart{cQ}    = 'TrunkSARSA'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos; 

% =========================================================
% BODY AWAY
sBdy.clc.stimDynams =     @(pos) pos - sBdy.clc.baseVel; % For receding, set speed negative

% -----------------------------
% POS REW
sBdy.clc.startRew =  1;
[newQ ] = CalcQDirect(sBdy);
% store q values and attributes
cQ = cQ +1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1;
allQ.rew(cQ)         = sBdy.clc.startRew;
allQ.bodyPart{cQ}    = 'TrunkSARSA'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos;

% -----------------------------
% NEG REW
sBdy.clc.startRew =  -1;
[newQ ] = CalcQDirect(sBdy);
% store q values and attributes
cQ = cQ + 1
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1;
allQ.rew(cQ)         = sBdy.clc.startRew;
allQ.bodyPart{cQ}    = 'TrunkSARSA'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos;

% #################################################################
% =========================================================
% HEAD TOWARDS
sHed.lp.alg = 'SARSA';
sHed.clc.stimDynams =     @(pos) pos + sHed.clc.baseVel;
sHed.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sHed.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

% -----------------------------
% POS REW
sHed.clc.startRew =  1;
[newQ ] = CalcQDirect(sHed);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}        = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sHed.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HeadSARSA'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos;

% -----------------------------
% NEG REW
sHed.clc.startRew =  -1;
[newQ ] = CalcQDirect(sHed);
cQ = cQ + 1;
% store q values and attributes
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sHed.clc.startRew;
allQ.bodyPart{cQ}    = 'HeadSARSA'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos;


% =========================================================
% HEAD AWAY
sHed.clc.stimDynams  =     @(pos) pos - sHed.clc.baseVel; % For receding, set speed negative

% -----------------------------
% POS REW
sHed.clc.startRew =  1;
[newQ ] = CalcQDirect(sHed); 
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1;
allQ.rew(cQ)         = sHed.clc.startRew;
allQ.bodyPart{cQ}    = 'HeadSARSA';
allQ.centerpos{cQ}   = sBdy.clc.nearPos;

% -----------------------------
% NEG REW
% REWARD SIZE
sHed.clc.startRew   =  -1;
[newQ ] = CalcQDirect(sHed);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}      = newQ;
allQ.dir(cQ)        = -1;
allQ.rew(cQ)        = sHed.clc.startRew;
allQ.bodyPart{cQ}   = 'HeadSARSA'; 
allQ.centerpos{cQ}  = sBdy.clc.nearPos; 


% #################################################################
% =========================================================
% HAND by chest, TOWARDS
sHnd.lp.alg = 'SARSA';
sHnd.clc.stimDynams = @(pos) pos + sHnd.clc.baseVel;
sHnd.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sHnd.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

%-----------------------------
% POS REW
sHnd.clc.startRew =  1;
[newQ ] = CalcQDirect(sHnd);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sHnd.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestSARSA'; 
allQ.centerpos{cQ}   = sHnd.clc.nearPos; 

% -----------------------------
% NEG REW
sHnd.clc.startRew    =  -1;
[newQ ] = CalcQDirect(sHnd);
% store q values and attributes
cQ = cQ +1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sHnd.clc.startRew; 
allQ.bodyPart{cQ}    = 'HandByChestSARSA'; 
allQ.centerpos{cQ}   = sHnd.clc.nearPos; 


% =========================================================
% HAND AWAY
sHnd.clc.stimDynams =     @(pos) pos - sHnd.clc.baseVel; 


% -----------------------------
% POS REW
sHnd.clc.startRew =  1;
[newQ ] = CalcQDirect(sHnd);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1;
allQ.rew(cQ)         = sHnd.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestSARSA'; 
allQ.centerpos{cQ}   = sHnd.clc.nearPos; 

% -----------------------------
% NEG REW
sHnd.clc.startRew =  -1;
[newQ ] = CalcQDirect(sHnd);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1;
allQ.rew(cQ)         = sHnd.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestSARSA'; 
allQ.centerpos{cQ}   = sHnd.clc.nearPos;


% #################################################################
% ==========================================================
% HAND by side, all directions 

% Move the data, then shift the Q values next
cQ = cQ + [1:4];
allQ(cQ,:)       = allQ(cQ - 4,:);
allQ.bodyPart(cQ)    = repmat({'HandBySideSARSA'},[1 4])';

% Shift the data
for iPsi = cQ
    allQ.qVals{iPsi} = zeros(size(allQ.qVals{iPsi}));

    % Flip the data because it's symetric anyway and the extent on the
    % right is bigger than on the left
    leftCopyDist  = min([handMidPos - 1 , size(    allQ.qVals{iPsi},3) - handNewPos]);
    rightCopyDist = min([handNewPos - 1 , size(    allQ.qVals{iPsi},3) - handMidPos]);
    
    % Move the data from the right of the central hand to the left of the
    % lateral hand, and flip it
    allQ.qVals{iPsi}(:, :, handNewPos + [-rightCopyDist:0]   ,:) = ...
        flip(allQ.qVals{iPsi-4}(:, :, handMidPos + [0:rightCopyDist] ,:),3);

    % Move the data from the left of the central hand to the right of the
    % lateral hand, and flip it
    allQ.qVals{iPsi}(:, :, handNewPos + [1:leftCopyDist],:) = ...
        flip(allQ.qVals{iPsi-4}(:, :, handMidPos + [-leftCopyDist:-1] ,:),3);
    allQ.centerpos{iPsi}   = sHnd.clc.nearPos;
    allQ.centerpos{iPsi} =  sHndSide.clc.nearPos;
        
end

% #################################################################
% TOOL USE
% =========================================================
% TOOL by HAND, TOWARDS

sHndTool.clc.stimDynams = @(pos) pos + sHndTool.clc.baseVel;

sHndTool.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sHndTool.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

%-----------------------------
% POS REW
% REWARD SIZE
sHndTool.lp.alg = 'SARSA';
sHndTool.clc.startRew =  1;
[newQ] = CalcQDirect(sHndTool);
% store q values and attributes
cQ = cQ(end) + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sHndTool.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolSARSA'; 
allQ.centerpos{cQ}   = sHndTool.clc.nearPos;

% -----------------------------
% NEG REW
sHndTool.clc.startRew =  -1;
[newQ] = CalcQDirect(sHndTool);
% store q values and attributes
cQ = cQ + 1
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sHndTool.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolSARSA'; 
allQ.centerpos{cQ}   = sHndTool.clc.nearPos; 

% -----------------------------
% POS REW, AWAY
sHndTool.clc.stimDynams =     @(pos) pos - sHndTool.clc.baseVel;
sHndTool.clc.startRew =  1;
[newQ] = CalcQDirect(sHndTool);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1;
allQ.rew(cQ)         = sHndTool.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolSARSA'; 
allQ.centerpos{cQ}   = sHndTool.clc.nearPos;

% -----------------------------
% NEG REW, AWAY
sHndTool.clc.startRew =  -1;
[newQ] = CalcQDirect(sHndTool);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1;
allQ.rew(cQ)         = sHndTool.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolSARSA'; 
allQ.centerpos{cQ}   = sHndTool.clc.nearPos;


% #################################################################
% Spider Butterfly Track; assume the uncertainties are too small, 
% therefore no new Q-value calculation is needed
% ==========================================================


% #################################################################
% TOOL version 2: rake
% ==========================================================
% TOOL by HAND, TOWARDS
sHndToolRake.lp.alg = 'SARSA';
sHndToolRake.clc.stimDynams = @(pos) pos + sHndToolRake.clc.baseVel;
sHndToolRake.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sHndToolRake.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

%-----------------------------
% POS REW
sHndToolRake.clc.startRew =  1;
[newQ ] = CalcQDirect(sHndToolRake);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sHndToolRake.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeSARSA'; 
allQ.centerpos{cQ}   = sHndToolRake.clc.nearPos; 

% -----------------------------
% NEG REW
sHndToolRake.clc.startRew =  -1;
[newQ ] = CalcQDirect(sHndToolRake);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sHndToolRake.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeSARSA'; 
allQ.centerpos{cQ}   = sHndToolRake.clc.nearPos;

%  -----------------------------
% POS REW, AWAY
sHndToolRake.clc.stimDynams = @(pos) pos - sHndToolRake.clc.baseVel;
sHndToolRake.clc.startRew =  1;
[newQ ] = CalcQDirect(sHndToolRake);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1;
allQ.rew(cQ)         = sHndToolRake.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeSARSA'; 
allQ.centerpos{cQ}   = sHndToolRake.clc.nearPos;

% -----------------------------
% NEG REW, AWAY
sHndToolRake.clc.startRew =  -1;
[newQ ] = CalcQDirect(sHndToolRake);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1;
allQ.rew(cQ)         = sHndToolRake.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeSARSA'; 
allQ.centerpos{cQ}   = sHndToolRake.clc.nearPos;


% #################################################################
% Monkey arm
% ==========================================================
% ARM FORWARD ['right']
sArm.lp.alg = 'SARSA';
sArm.clc.stimDynams = @(pos) pos + sArm.clc.baseVel;
sArm.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sArm.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

% -----------------------------
% POS REW
sArm.clc.startRew =  1;
[newQ ] = CalcQDirect(sArm);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sArm.clc.startRew;
allQ.bodyPart{cQ}    = 'ArmForwardSARSA'; 
allQ.centerpos{cQ}   = sArm.clc.nearPos;


% -----------------------------
% NEG REW
sArm.clc.startRew =  -1;
[newQ ] = CalcQDirect(sArm);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sArm.clc.startRew;
allQ.bodyPart{cQ}    = 'ArmForwardSARSA'; 
allQ.centerpos{cQ}   = sArm.clc.nearPos;


% ==========================================================
% ARM LEFT
sArm.lp.alg = 'SARSA';
sArmL.clc.stimDynams = @(pos) pos + sArmL.clc.baseVel;
sArmL.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sArmL.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

%-----------------------------
% POS REW
sArmL.lp.alg = 'SARSA';
sArmL.clc.startRew =  1;
[newQ ] = CalcQDirect(sArmL);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sArmL.clc.startRew;
allQ.bodyPart{cQ}    = 'ArmLeftSARSA'; 
allQ.centerpos{cQ}   = sArmL.clc.nearPos;


% -----------------------------
% NEG REW
% REWARD SIZE
sArmL.clc.startRew =  -1;
[newQ ] = CalcQDirect(sArmL);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sArmL.clc.startRew;
allQ.bodyPart{cQ}    = 'ArmLeftSARSA'; 
allQ.centerpos{cQ}   = sArmL.clc.nearPos;



% #################################################################
% =========================================================
% MONKEY HEAD TOWARDS - CONSTRAINED

sHedConstr.lp.alg = 'SARSA';
sHedConstr.clc.stimDynams = @(pos) pos + sHedConstr.clc.baseVel;

sHedConstr.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sHedConstr.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

% -----------------------------
% POS REW
% REWARD SIZE
sHedConstr.clc.startRew =  1;
[newQ ] = CalcQDirect(sHedConstr);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sHedConstr.clc.startRew;
allQ.bodyPart{cQ}    = 'HeadConstrSARSA'; 
allQ.centerpos{cQ}   = sHedConstr.clc.nearPos;

% -----------------------------
% NEG REW
% REWARD SIZE
sHedConstr.clc.startRew =  -1;
[newQ ] = CalcQDirect(sHedConstr);
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sHedConstr.clc.startRew;
allQ.bodyPart{cQ}    = 'HeadConstrSARSA'; 
allQ.centerpos{cQ}   = sHedConstr.clc.nearPos; 


% #################################################################
% =========================================================
% 'MONKEY' HEAD 15 degree angle
%  Move the straight-on monkey head data, then shift the Q values

cQ= cQ + [1:2];
allQ(cQ,:)       = allQ(cQ - 2,:);
allQ.bodyPart(cQ)    = repmat({'RotatedHeadSARSA'},[1 2])';
% Shift the data
for iPsi = cQ
    % shift the q values to the side by 15 degree, i.e. once for every 4
    % rows, starting at the face
    shiftAmount = 1;
    for iR = (sHed15.clc.nearPos(1) + 2) : -1 : 1
        allQ.qVals{iPsi}(:,iR, 1:end-shiftAmount ,:) = allQ.qVals{iPsi-2}(:,iR, (shiftAmount+1):end ,:);

        if mod(iR - (sHed15.clc.nearPos(1) + 3),4) == 0
            shiftAmount = shiftAmount + 1;
        end
    end
end


%% ------------------------------------------------------------------------
% #################################################################
% Create data simulating stimulation in the rear space -- necessary for
% fitting the Taffou 2014 data added during revision

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

sBack          = s;
sBack.clc.binW = 10; % Bins are bigger cause space for this experiment is bigger

% Because the body is so small compared to the whole modelled volume, don't
% set all the points behind the body to be rewarding
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
sBack.clc.actConsequence = [ ...
     0  0  0 ; ... % action 1 stay
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

% -----------------------------
% POS REW
sBack.clc.startRew   =  1;
[newQ ] = CalcQDirect(sBack);
% store q values and attributes
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sBack.clc.startRew;
allQ.bodyPart{cQ}    = 'Back'; 
allQ.centerpos{cQ}   = sBack.clc.nearPos; 


% -----------------------------
% NEG REW
cQ = cQ + 1;
sBack.clc.startRew   =  -1;
[newQ ] = CalcQDirect(sBack);
% store q values and attributes
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sBack.clc.startRew; 
allQ.bodyPart{cQ}    = 'Back'; 
allQ.centerpos{cQ}   = sBack.clc.nearPos; 


% Also calculate SARSA HP, and MI for the back ================================
% ================================= SARSA =================================
cQ = 145;

sBack.lp.alg = 'SARSA';
sBack.clc.nReps = 2; 

% -----------------------------
% POS REW 
cQ = cQ + 1;
sBack.clc.startRew   =  1;
[newQ ] = CalcQDirect(sBack);
% store q values and attributes
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sBack.clc.startRew;
allQ.bodyPart{cQ}    = 'BackSARSA'; 
allQ.centerpos{cQ}   = sBack.clc.nearPos; 


% -----------------------------
% NEG REW
cQ = cQ + 1;
sBack.clc.startRew   =  -1;
[newQ ] = CalcQDirect(sBack);
% store q values and attributes
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sBack.clc.startRew;
allQ.bodyPart{cQ}    = 'BackSARSA'; 
allQ.centerpos{cQ}   = sBack.clc.nearPos; 


% ================================ HITPROB ================================
cQ = 147;

sBackHP                     = sBack;
sBackHP.clc.actConsequence  = [0 0 0]; 
sBackHP.clc.gammaVal        = 1;
sBackHP.clc.startRew        = 1;
sBackHP.clc.nReps           = 1; 

% Hit probability - body
% -----------------------------
% NOTE: deltaT is 0.5, so I have to halve all the velocities and spreads

% Random stimulus dynamics
rSpr = [-1 0 1];
rSprPr = gaussmf(rSpr,[1 0]) ./ sum(gaussmf(rSpr,[0.5 0]));
cSpr = [-1 0 1];
cSprPr = gaussmf(cSpr,[1 0]) ./ sum(gaussmf(cSpr,[0.5 0]));
zSpr = [-1 0 1];
zSprPr = gaussmf(zSpr,[1 0]) ./ sum(gaussmf(zSpr,[0.5 0]));

% divide probabilities over 2 because deltaT = 1/2
rSprPr = rSprPr./([2 1 2]) + rSprPr(1).*([0 1 0]); 
cSprPx = cSprPr./([2 1 2]) + cSprPr(1).*([0 1 0]); 
zSprPr = zSprPr./([2 1 2]) + zSprPr(1).*([0 1 0]); 

% x y z, Deterministic stimulus dynamics
sBackHP.clc.stimDynams =     @(pos) pos + s.clc.baseVel./2; % For approaching, set speed positive
sBackHP.clc.randSpread =     {rSpr cSpr zSpr}; % Put a little biit of x and z variability in? Kind of arbitrary
sBackHP.clc.spreadProb =     {rSprPr cSprPr  zSprPr};

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


% ======================== Multisensory Integration =======================
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
rSigmaMI      = sqrt(13.^2 + 11.^2) ./ sBackMI.clc.binW ;
rSnsSprMIPr   = gaussmf(rSnsSprMI,[rSigmaMI 0]) ./ sum(gaussmf(rSnsSprMI,[rSigmaMI 0])); 
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
allQ.dir(cQ)         = 1;
allQ.rew(cQ)         = sBackMI.clc.startRew;
allQ.bodyPart{cQ}    = 'BackMI'; 
allQ.centerpos{cQ}   = sBackMI.clc.nearPos;

% =============================== Distance ===============================
cQ = 150;
allQ.qVals{cQ}       = permute(CalcDistToBody(sBackMI),[4 1 2 3]);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sBackMI.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'BackDist'; 
allQ.centerpos{cQ}   = sBackMI.clc.nearPos;


%% #################################################################
%  Create Q-values with added uncertainty, by blurring existing values
% (i.e. assuming that the agent has learned Q under less uncertainty, and
% is now simply trying to apply that learned Q-value)

% Sensory stimulus spreading
rSnsSpr = [-3:3]; 
cSnsSpr = [-3:3]; 
zSnsSpr = [-3:3]; 

% Set uncertainties to near-optimal values
rSigma      = sqrt(2.5.^2 + 2.5.^2) ./ sBdyHP.clc.binW ./2;
rSnsSprPr   = gaussmf(rSnsSpr,[rSigma 0]) ./ sum(gaussmf(rSnsSpr,[rSigma 0]));
czSigma     = sqrt(1.^2 + 1.^2) ./ sBdyHP.clc.binW ./2;
cSnsSprPr   = gaussmf(cSnsSpr,[czSigma 0]) ./ sum(gaussmf(cSnsSpr,[czSigma 0]));
zSigma      = sqrt(2.5.^2 + 2.5.^2) ./ sBdyHP.clc.binW ./2;
zSnsSprPr   = gaussmf(zSnsSpr,[zSigma 0]) ./ sum(gaussmf(zSnsSpr,[zSigma 0]));

s                   = sBdy;
s.clc.sensSpread    = {rSnsSpr cSnsSpr zSnsSpr};
s.clc.sensProb      = {rSnsSprPr cSnsSprPr zSnsSprPr};

% Table indices of the q-values which are to be blurred
baseQs = [1:20 27:38 142 143]; % 

% Indices of the table which will store the blurred Qs
uncQs = [78:109 144 145];

% First just copy the original data
allQ(uncQs,:) = allQ(baseQs,:);

for iQ = 1:length(baseQs)
    disp(['Current iteration ' num2str(iQ) '/' num2str(length(baseQs)) ]);
    currQ = allQ.qVals{uncQs(iQ)};

    newQ = currQ;
    % Loop through y, then x, then z, and calculate the new blurred q
    % values
    maxActDists = cellfun(@(dim) max(abs(dim)), s.clc.sensSpread);
    for iSR     = s.wrld.size(1) - maxActDists(1) : -1 : maxActDists(1) +1
        for iSC = (1 + maxActDists(2)) : (s.wrld.size(2) - maxActDists(2))
            for iSZ = 1 + maxActDists(3) : s.wrld.size(3) - maxActDists(3)

                snsSprR = iSR + s.clc.sensSpread{1};
                snsSprC = iSC + s.clc.sensSpread{2};
                snsSprZ = iSZ + s.clc.sensSpread{3};

                % Blur the Q-values as dictated by the uncertainty, and put
                % them into the new Q
                tmpQ = zeros([size(newQ,1) 1]);
                for iSnsR = 1:length(snsSprR)
                for iSnsC = 1:length(snsSprC)
                for iSnsZ = 1:length(snsSprZ)
                    tmpQ = tmpQ + ...
                        s.clc.sensProb{1}(iSnsR) .* s.clc.sensProb{2}(iSnsR) .* s.clc.sensProb{3}(iSnsR) .* ...
                        currQ(:,snsSprR(iSnsR),snsSprC(iSnsC),snsSprZ(iSnsZ));
                end
                end
                end
            end
        end
    end
    allQ.qVals{uncQs(iQ)}    =  newQ;
    allQ.bodyPart{uncQs(iQ)} = [allQ.bodyPart{baseQs(iQ)} 'Unc'];  
end

