addpath('F:\Programs\Matlab\Utilities\')
addpath('F:\Programs\Matlab\Utilities\plotting\colormaps\')
addpath('F:\Programs\Matlab\Utilities\plotting\Plot_VoxelSurf_V1_11')

addpath('F:\Programs\Matlab\Utilities\')
addpath('F:\Programs\Matlab\Utilities\plotting\');
addpath('F:\Programs\Matlab\Toolboxes\CircStat2012a');
addpath('F:\Programs\Matlab\Utilities\plotting\colormaps');


load('Results\ForFigures\Fig1_Results_v3')
s=rS(end).s;
w=rS(end).w;
addpath(genpath('F:\Projects\DPPS\DefenseAgent\Scripts\rlsimplepps'));
addpath('F:\Programs\Matlab\Utilities\plotting\colormaps');

% Set gammas and stuff
allGammas = 0.7;

% settings for plot
sFP=s;
sFP.plt.lmbRow = s.wrld.size(1)-2;
sFP.plt.rowLims=[6.5 s.wrld.size(1)-1.5];
sFP.plt.colLims=[3.5 s.wrld.size(2)-3.5];
sFP.plt.cBarFl=0;
sFP.plt.meanLimbCols=1;
sFP=DefaultSettings(sFP);
sFP.plt.axesVis=0;

% figure settings
fS.gridXstart = -4.5;
fS.gridXstep = 1;
fS.gridYstart = 3.5;
fS.gridYstep = 1;

% Load existing calculations
% load('F:\Projects\DPPS\DefenseAgent\Data\FittingWorkspace_10iterations_WIthTool_WIthExpandedHand_2Iterations_PLUSDIAGONALS_2MoreIterations_GOODFIT_2_PLusRake_PlusArm_asifTool_plusRotHead_CONSTR.mat')
% load('F:\Projects\DPPS\DefenseAgent\Data\NewWithHitProbs_MultSensInt_AndDist_AndUncertainQ.mat')
% load('F:\Projects\DPPS\DefenseAgent\Data\NewWithHitProbs_MultSensInt_AndDist_AndUncertainQ_FromObserverPerspective.mat')
load('F:\Projects\DPPS\DefenseAgent\Data\NewWithHitProbs_MultSensInt_AndDist_AndUncertainQ_FromObserverPerspective_v2.mat')

%% -------------------------------------------------------------------------
% NOEL 2015 PPS SIZE

% #################################################################
% =========================================================
% BODY TOWARDS

% Furthest distance = 100cm. Inter-stim dist between 15 and 22 cm
% --> to allow for speed effects, maybe use 5cm bin size?
% Then make world 305cm by 135cm by 155cm? --> allow hand to be placed to the
% side, and head to be placed above body
% Then I guess we can put the hand volume to 2 blocks? let's try
s.wrld.size = [61 27 31];
s.clc.nearPos = [s.wrld.size(1)-12 11 13]';

s.clc.nReps = 1; % 

% NOTE this involves leaving 20cm behind the participants' back


allQ = table();

s.clc.gammaVal   = 0.8;
s.clc.baseVel    = [4 0 0]; 

% % % % x y z, Deterministic stimulus dynamics
% % % s.clc.stimDynams =     @(pos) pos + [0 0 0]; % For approaching, set speed positive
% % % s.clc.randSpread =     {[0 1] [-1   0  1] [-1  0  1]}; % Put a little biit of x and z variability in? Kind of arbitrary
% % % s.clc.spreadProb =     {[.5 .5] [0.2 .6 0.2]  [ .2 .6 .2]}; % x y z, probabilities of spread

% Random stimulus dynamics
rSpr = [-1 0 1];
rSprPr = gaussmf(rSpr,[1 0]) ./ sum(gaussmf(rSpr,[1 0]));
cSpr = [-1 0 1]; %%% cSpr = [ -1 0 1 ];
cSprPr = gaussmf(cSpr,[1 0]) ./ sum(gaussmf(cSpr,[1 0]));
zSpr = [-1 0 1]; %%% zSpr = [ -1 0 1 ];
zSprPr = gaussmf(zSpr,[1 0]) ./ sum(gaussmf(zSpr,[1 0]));

% x y z, Deterministic stimulus dynamics
s.clc.stimDynams =     @(pos) pos + s.clc.baseVel; % For approaching, set speed positive
s.clc.randSpread =     {rSpr cSpr zSpr}; % Put a little biit of x and z variability in? Kind of arbitrary
s.clc.spreadProb =     {rSprPr cSprPr  zSprPr}; % x y z, probabilities of spread


% NOTE: If the limb can move at different speeds, BOTH speeds have to be
% put in as potential actions, because the model doesn't have any collision
% calculation that takes into account overshoot

% FIRST CALCULATE FOR BODY - moves more slowly than limb, let's say. Also
% ONLY has negative potential rewards
sBdy = s;
sBdy.clc.actConsequence = [ 0  0  0 ; ... % action 1 stay
    0  1  0 ; ... % action 2 left
    0 -1  0 ; ... % action 3 right
    0  0  1 ; ... % action 4 up
    0  0 -1 ; ... % action 5 down
    -1  0  0 ; ... % action 6 forward
    1  0  0 ;...    % action 7 backward-
    0  2  0 ; ... % action 2 left
    0 -2  0 ; ... % action 3 right
    0  0  2 ; ... % action 4 up
    0  0 -2 ; ... % action 5 down
    -2  0  0 ; ... % action 6 forward
    2  0  0 ];    % action 7 backward-

% Location and size of limb, and where the Q-values are calculated FROM
sBdy.clc.startSR = []; sBdy.clc.startSC = []; sBdy.clc.startSZ = [];


bdyRs = s.clc.nearPos(1) + [2 1 1 0 0 0 1 1 2];
bdyCs = s.clc.nearPos(2) + [-4:4];
bdyZs = s.clc.nearPos(3) + [-5:5];
iRew = 0; % Initialise rewarded block counter
for  iZ = 1:length(bdyZs)
    for  iC = 1:length(bdyCs)
        iRew = iRew + 1;
        sBdy.clc.startSR(iRew) =  bdyRs(iC);
        sBdy.clc.startSC(iRew) =  bdyCs(iC);
        sBdy.clc.startSZ(iRew) =  bdyZs(iZ);
    end
end

% % % % -----------------------------
% % % % POS REW
% % % % REWARD SIZE
% % % sBdy.clc.startRew =  1;
% % % % % % [newQ ] = CalcQDirect(sBdy);
% % % [newQ ] = CalcQDirect(sBdy);
% % % % store q values and attributes
% % % allQ.qVals{1}        = newQ;
% % % allQ.dir(1)         = 1; %s.clc.stimDynams([0 0 0]); % towards
% % % allQ.rew(1)         = sBdy.clc.startRew; % towards
% % % allQ.bodyPart{1}        = 'Trunk'; 
% % % allQ.centerpos{1}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

%%


% -----------------------------
% NEG REW
% REWARD SIZE
sBdy.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sBdy);
[newQ ] = CalcQDirect(sBdy,allQ.qVals{2});
% store q values and attributes
allQ.qVals{2}        = newQ;
allQ.dir(2)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(2)         = sBdy.clc.startRew; % towards
allQ.bodyPart{2}        = 'Trunk'; 
allQ.centerpos{2}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% =========================================================
% BODY AWAY
sBdy.clc.stimDynams =     @(pos) pos - s.clc.baseVel; % For receding, set speed negative

% -----------------------------
% POS REW
% REWARD SIZE
sBdy.clc.startRew =  1;
% % % [newQ ] = CalcQDirect(sBdy);
[newQ ] = CalcQDirect(sBdy,allQ.qVals{3});
% store q values and attributes
allQ.qVals{3}        = newQ;
allQ.dir(3)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(3)         = sBdy.clc.startRew; % towards
allQ.bodyPart{3}        = 'Trunk'; 
allQ.centerpos{3}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

% -----------------------------
% NEG REW
% REWARD SIZE
sBdy.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sBdy);
[newQ ] = CalcQDirect(sBdy,allQ.qVals{4});
% store q values and attributes
allQ.qVals{4}        = newQ;
allQ.dir(4)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(4)         = sBdy.clc.startRew; % towards
allQ.bodyPart{4}        = 'Trunk'; 
allQ.centerpos{4}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% #################################################################
% =========================================================
% HEAD TOWARDS

sHed = s;
sHed.clc.nearPos = s.clc.nearPos + [0 0 11]';
sHed.clc.actConsequence = sBdy.clc.actConsequence;
% % % [ 0  0  0 ; ... % action 1 stay
% % %     0  1  0 ; ... % action 2 left
% % %     0 -1  0 ; ... % action 3 right
% % %     0  0  1 ; ... % action 4 up
% % %     0  0 -1 ; ... % action 5 down
% % %     -1  0  0 ; ... % action 6 forward
% % %     1  0  0 ];    % action 7 backward-
% Location and size of limb, and where the Q-values are calculated FROM
sHed.clc.startSR = []; sHed.clc.startSC = []; sHed.clc.startSZ = [];
hedRs = sHed.clc.nearPos(1) + [1 0 1];
hedCs = sHed.clc.nearPos(2) + [-1:1];
hedZs = sHed.clc.nearPos(3) + [-2:2];
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
sHed.clc.startRew =  1;
% % % [newQ ] = CalcQDirect(sHed);
[newQ ] = CalcQDirect(sHed,allQ.qVals{5});
% store q values and attributes
allQ.qVals{5}        = newQ;
allQ.dir(5)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(5)         = sHed.clc.startRew; % towards
allQ.bodyPart{5}        = 'Head'; 
allQ.centerpos{5}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW
% REWARD SIZE
sHed.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sHed);
[newQ ] = CalcQDirect(sHed,allQ.qVals{6});
% store q values and attributes
allQ.qVals{6}        = newQ;
allQ.dir(6)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(6)         = sHed.clc.startRew; % towards
allQ.bodyPart{6}        = 'Head'; 
allQ.centerpos{6}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% =========================================================
% HEAD AWAY
sHed.clc.stimDynams =     @(pos) pos - s.clc.baseVel; % For receding, set speed negative

% -----------------------------
% POS REW
% REWARD SIZE
sHed.clc.startRew =  1;
% % % [newQ ] = CalcQDirect(sHed);
[newQ ] = CalcQDirect(sHed,allQ.qVals{7});
% store q values and attributes
allQ.qVals{7}        = newQ;
allQ.dir(7)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(7)         = sHed.clc.startRew; % towards
allQ.bodyPart{7}        = 'Head'; 
allQ.centerpos{7}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

% -----------------------------
% NEG REW
% REWARD SIZE
sHed.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sHed);
[newQ ] = CalcQDirect(sHed,allQ.qVals{8});
% store q values and attributes
allQ.qVals{8}        = newQ;
allQ.dir(8)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(8)         = sHed.clc.startRew; % towards
allQ.bodyPart{8}        = 'Head'; 
allQ.centerpos{8}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%% #################################################################
% ===============================================s==========
% HAND by chest, TOWARDS

sHnd = s;
sHnd.clc.nearPos = sHnd.clc.nearPos;
% Hand movementrs need to be allowed to be faster
sHnd.clc.actConsequence = [ 0  0  0 ; ... % action 1 stay
    0  1  0 ; ... % action 2 left SLOW
    0 -1  0 ; ... % action 3 right SLOW
    0  0  1 ; ... % action 4 up SLOW
    0  0 -1 ; ... % action 5 down SLOW
    -1  0  0 ; ... % action 6 forward SLOW
    1  0  0 ; ... % action 7 backward SLOW
    0  2  0 ; ... % action 2 left MEDIUM
    0 -2  0 ; ... % action 3 right MEDIUM
    0  0  2 ; ... % action 4 up MEDIUM
    0  0 -2 ; ... % action 5 down MEDIUM
    -2  0  0 ; ... % action 12 forward MEDIUM
    2  0  0 ; ... % action 13 backward MEDIUM
    -3  0  0 ; ... % action 18 forward FAST
    3  0  0 ; ... % action 19 backward FAST
    -4  0  0 ; ... % action 18 forward FAST
    4  0  0 ; ... % action 19 backward FAST
    -5  0  0 ; ... % action 18 forward FAST
    5  0  0 ; ... % action 19 backward FAST
    0  3  0 ; ... % action 2 left MEDIUM
    0 -3  0 ; ... % action 3 right MEDIUM
    0  0  3 ; ... % action 4 up MEDIUM
    0  0 -3 ; ... % action 5 down MEDIUM
    0  4  0 ; ... % action 2 left MEDIUM
    0 -4  0 ; ... % action 3 right MEDIUM
    0  0  4 ; ... % action 4 up MEDIUM
    0  0 -4 ; ... % action 5 down MEDIUM
    0  5  0 ; ... % action 2 left MEDIUM
    0 -5  0 ; ... % action 3 right MEDIUM
    0  0  5 ; ... % action 4 up MEDIUM
    0  0 -5 ; ...
    1  3  0 ; ... % action 2 left MEDIUM
    1 -3  0 ; ... % action 3 right MEDIUM
    1  0  3 ; ... % action 4 up MEDIUM
    1  0 -3 ; ... % action 5 down MEDIUM
    1  4  0 ; ... % action 2 left MEDIUM
    1 -4  0 ; ... % action 3 right MEDIUM
    1  0  4 ; ... % action 4 up MEDIUM
    1  0 -4 ; ... % action 5 down MEDIUM
    1  5  0 ; ... % action 2 left MEDIUM
    1 -5  0 ; ... % action 3 right MEDIUM
    1  0  5 ; ... % action 4 up MEDIUM
    1  0 -5 ; ...
    6  0  0 ; ... Forward utramegafast
    7  0  0 ] ; %Forward utramegafast



% Location and size of limb, and where the Q-values are calculated FROM
sHnd.clc.startSR = []; sHnd.clc.startSC = []; sHnd.clc.startSZ = [];
hndRs = sHnd.clc.nearPos(1) + [0 0 0];
hndCs = sHnd.clc.nearPos(2) + [-1:1];
hndZs = sHnd.clc.nearPos(3) + [-1:1];
iRew = 0; % Initialise rewarded block counter
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
sHnd.clc.startRew =  1;
% % % [newQ ] = CalcQDirect(sHnd,);
[newQ ] = CalcQDirect(sHnd,allQ.qVals{9});
% store q values and attributes
allQ.qVals{9}        = newQ;
allQ.dir(9)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(9)         = sHnd.clc.startRew; % towards
allQ.bodyPart{9}        = 'HandByChest'; 
allQ.centerpos{9}   = sHnd.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW
% REWARD SIZE
sHnd.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sHnd);
[newQ ] = CalcQDirect(sHnd,allQ.qVals{10});
% store q values and attributes
allQ.qVals{10}        = newQ;
allQ.dir(10)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(10)         = sHnd.clc.startRew; % towards
allQ.bodyPart{10}        = 'HandByChest'; 
allQ.centerpos{10}   = sHnd.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% =========================================================
% HAND AWAY
sHnd.clc.stimDynams =     @(pos) pos - s.clc.baseVel; % For receding, set speed negative


% -----------------------------
% POS REW
% REWARD SIZE
sHnd.clc.startRew =  1;
% % % [newQ ] = CalcQDirect(sHnd);
[newQ ] = CalcQDirect(sHnd,allQ.qVals{11});
% store q values and attributes
allQ.qVals{11}        = newQ;
allQ.dir(11)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(11)         = sHnd.clc.startRew; % towards
allQ.bodyPart{11}        = 'HandByChest'; 
allQ.centerpos{11}   = sHnd.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

% -----------------------------
% NEG REW
% REWARD SIZE
sHnd.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sHnd);
[newQ ] = CalcQDirect(sHnd,allQ.qVals{12});
% store q values and attributes
allQ.qVals{12}        = newQ;
allQ.dir(12)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(12)         = sHnd.clc.startRew; % towards
allQ.bodyPart{12}        = 'HandByChest'; 
allQ.centerpos{12}   = sHnd.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%% #################################################################
% ===============================================s==========
% HAND by side, all directions 

% $$$ --> NEXT FIGURE OUT HOW TO FLIP THIS BADBOY so that the useful data
% ends up on the LEFT of the hand --> bigger space to use :D

sHndSide = sHnd;
sHndSide.clc.nearPos = sHnd.clc.nearPos + [0 10 0]';
sHndSide.clc.startSC = sHnd.clc.startSC + 10;


% Move the data, then shift the Q values next
allQ(13:16,:)       = allQ(9:12,:);
allQ.bodyPart(13:16)    = repmat({'HandBySide'},[1 4])';

% $$$ CLEAN THIS UP w.r.t. new hndside.clc.nearpos defined above
handMidPos      = sHnd.clc.nearPos(2);
handNewPos      = sHndSide.clc.nearPos(2);
sideHandShift   = handNewPos - handMidPos;

% Find transposable width on the narrower side
moveWidth = s.wrld.size(2) - handNewPos;

% Shift the data
for iPsi = 13:16
    allQ.qVals{iPsi} = zeros(size(allQ.qVals{iPsi}));
% % %     allQ.qVals{iPsi}(:, :, sideHandShift : end,:) = ...
% % %         allQ.qVals{iPsi-4}(:, :, 1 : moveWidth + handMidPos + 1,:);
% % % 
% % %     allQ.centerpos{iPsi}   = sHnd.clc.nearPos;

    % FLip the data because it's symetric anyway and the extent on the
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


%% #################################################################
% TOOL TIME!
% ===============================================s==========
% TOOL by HAND, TOWARDS

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

sHndTool.clc.startRew =  1;
% % % [newQ ] = CalcQDirect(sHndTool);
[newQ ] = CalcQDirect(sHndTool,allQ.qVals{17});
% store q values and attributes
allQ.qVals{17}       = newQ;
allQ.dir(17)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(17)         = sHndTool.clc.startRew; % towards
allQ.bodyPart{17}    = 'HandByChestPlusTool'; 
allQ.centerpos{17}   = sHndTool.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW
% REWARD SIZE
sHndTool.clc.stimDynams =     @(pos) pos + s.clc.baseVel;

sHndTool.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sHndTool);
[newQ ] = CalcQDirect(sHndTool,allQ.qVals{18});
% store q values and attributes
allQ.qVals{18}       = newQ;
allQ.dir(18)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(18)         = sHndTool.clc.startRew; % towards
allQ.bodyPart{18}    = 'HandByChestPlusTool'; 
allQ.centerpos{18}   = sHndTool.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

% -----------------------------
% POS REW, AWAY
% REWARD SIZE
sHndTool.clc.stimDynams =     @(pos) pos - s.clc.baseVel;

sHndTool.clc.startRew =  1;
% % % [newQ ] = CalcQDirect(sHndTool);
[newQ ] = CalcQDirect(sHndTool,allQ.qVals{19});
% store q values and attributes
allQ.qVals{19}       = newQ;
allQ.dir(19)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(19)         = sHndTool.clc.startRew; % towards
allQ.bodyPart{19}    = 'HandByChestPlusTool'; 
allQ.centerpos{19}   = sHndTool.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW, AWAY
% REWARD SIZE
sHndTool.clc.stimDynams =     @(pos) pos - s.clc.baseVel;

sHndTool.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sHndTool);
[newQ ] = CalcQDirect(sHndTool,allQ.qVals{20});
% store q values and attributes
allQ.qVals{20}       = newQ;
allQ.dir(20)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(20)         = sHndTool.clc.startRew; % towards
allQ.bodyPart{20}    = 'HandByChestPlusTool'; 
allQ.centerpos{20}   = sHndTool.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%% #################################################################
% Spider Butterfly Track
% Here I'll crete one positive reward, and TWO negative rewards, I think.
% One where the negative is equal reward to the positive, and the other
% where negative is double the absolute reward.
% Then I combine those with the STANDARD hand-field, and maybe also the
% body field
% ===============================================s==========



sTrack = sHnd;
sTrack.clc.nReps = 1;
% Hand can't move for the track experiment
sTrack.clc.actConsequence = [ 0  0  0 ];    % action 1 stay 


% Random stimulus dynamics
rSprT = [0];
rSprPrT = [1];
cSprT = [0]; %%% cSpr = [ -1 0 1 ];
cSprPrT = [1];
zSprT = [0]; %%% zSpr = [ -1 0 1 ];
zSprPrT = [1];

% x y z, Deterministic stimulus dynamics
[newPos,sTrack] = TrackDynamics(pos,sTrack)

sTrack.clc.stimDynams =     @(pos) TrackDynamics(pos,sTrack); % Use the track dynamics
sTrack.clc.randSpread =     {rSprT cSprT zSprT}; % Put a little biit of x and z variability in? Kind of arbitrary
sTrack.clc.spreadProb =     {rSprPrT cSprPrT  zSprPrT}; % x y z, probabilities of spread

%-----------------------------
% POS REW
% REWARD SIZE 
sTrack.clc.startRew =  1;
[newQ ] = CalcQDirect(sTrack);
% store q values and attributes
allQ.qVals{21}       = newQ;
allQ.dir(21)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(21)         = sTrack.clc.startRew; % towards
allQ.bodyPart{21}    = 'HandByChestTrack'; 
allQ.centerpos{21}   = sTrack.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%-----------------------------
% NEG REW
% REWARD SIZE 
sTrack.clc.startRew =  -1;
[newQ ] = CalcQDirect(sTrack);
% store q values and attributes
allQ.qVals{22}       = newQ;
allQ.dir(22)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(22)         = sTrack.clc.startRew; % towards
allQ.bodyPart{22}    = 'HandByChestTrack'; 
allQ.centerpos{22}   = sTrack.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

%-----------------------------
% NEG REW, ~![BIG]!~
% REWARD SIZE 
sTrack.clc.startRew =  -1.25;
[newQ ] = CalcQDirect(sTrack);
% store q values and attributes
allQ.qVals{23}       = newQ;
allQ.dir(23)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(23)         = sTrack.clc.startRew; % towards
allQ.bodyPart{23}    = 'HandByChestTrack'; 
allQ.centerpos{23}   = sTrack.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%% #################################################################
% Make extra q with reward -2 for the more threatening situation for body as well
allQ(24,:) = allQ(2,:);
allQ.qVals{24} = 1.25 .* allQ.qVals{2};
allQ.rew(24)   = 1.25 .* allQ.rew(24);
allQ(25,:)     = allQ(6,:);
allQ.qVals{25} = 1.25 .* allQ.qVals{6};
allQ.rew(25)   = 1.25 .* allQ.rew(6);

% #################################################################
% Make extra q with NO VALUES because nothing is active
allQ(26,:)          = allQ(1,:);
allQ.qVals{26}      = zeros(size(allQ.qVals{1}));
allQ.rew(26)        = 0;
allQ.bodyPart{26}   = 'None'; 

%% #################################################################
% TOOL TIME version 2: rake
% ===============================================s==========
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
% REWARD SIZE
sHndToolRake.clc.stimDynams =     @(pos) pos + s.clc.baseVel;

sHndToolRake.clc.startRew =  1;
% % % [newQ ] = CalcQDirect(sHndToolRake);
[newQ ] = CalcQDirect(sHndToolRake,allQ.qVals{27});
% store q values and attributes
allQ.qVals{27}       = newQ;
allQ.dir(27)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(27)         = sHndToolRake.clc.startRew; % towards
allQ.bodyPart{27}    = 'HandByChestPlusToolRake'; 
allQ.centerpos{27}   = sHndToolRake.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW
% REWARD SIZE
sHndToolRake.clc.stimDynams =     @(pos) pos + s.clc.baseVel;

sHndToolRake.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sHndToolRake);
[newQ ] = CalcQDirect(sHndToolRake,allQ.qVals{28});
% store q values and attributes
allQ.qVals{28}       = newQ;
allQ.dir(28)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(28)         = sHndToolRake.clc.startRew; % towards
allQ.bodyPart{28}    = 'HandByChestPlusToolRake'; 
allQ.centerpos{28}   = sHndToolRake.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

%  -----------------------------
% POS REW, AWAY
% REWARD SIZE
sHndToolRake.clc.stimDynams =     @(pos) pos - s.clc.baseVel;

sHndToolRake.clc.startRew =  1;
% % % [newQ ] = CalcQDirect(sHndToolRake);
[newQ ] = CalcQDirect(sHndToolRake,allQ.qVals{29});
% store q values and attributes
allQ.qVals{29}       = newQ;
allQ.dir(29)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(29)         = sHndToolRake.clc.startRew; % towards
allQ.bodyPart{29}    = 'HandByChestPlusToolRake'; 
allQ.centerpos{29}   = sHndToolRake.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW, AWAY
% REWARD SIZE
sHndToolRake.clc.stimDynams =     @(pos) pos - s.clc.baseVel;

sHndToolRake.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sHndToolRake);
[newQ ] = CalcQDirect(sHndToolRake,allQ.qVals{30});
% store q values and attributes
allQ.qVals{30}       = newQ;
allQ.dir(30)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(30)         = sHndToolRake.clc.startRew; % towards
allQ.bodyPart{30}    = 'HandByChestPlusToolRake'; 
allQ.centerpos{30}   = sHndToolRake.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%% #################################################################
% Monkey arm
% ===============================================s==========
% ARM FORWARD ['right']

sBdy.clc.nReps = 2;

sHnd.clc.nReps = 1;
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
% REWARD SIZE
sArm.clc.stimDynams =     @(pos) pos + s.clc.baseVel;

sArm.clc.startRew =  1;
[newQ ] = CalcQDirect(sArm);
% % % [newQ ] = CalcQDirect(sArm,allQ.qVals{31});
% store q values and attributes
allQ.qVals{31}       = newQ;
allQ.dir(31)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(31)         = sArm.clc.startRew; % towards
allQ.bodyPart{31}    = 'ArmForward'; 
allQ.centerpos{31}   = sArm.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central



% -----------------------------
% NEG REW
% REWARD SIZE
sArm.clc.startRew =  -1;
[newQ ] = CalcQDirect(sArm);
% % % [newQ ] = CalcQDirect(sArm,allQ.qVals{32});
% store q values and attributes
allQ.qVals{32}       = newQ;
allQ.dir(32)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(32)         = sArm.clc.startRew; % towards
allQ.bodyPart{32}    = 'ArmForward'; 
allQ.centerpos{32}   = sArm.clc.nearPos;



%% ===============================================s==========
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
% REWARD SIZE
sArmL.clc.stimDynams =     @(pos) pos + s.clc.baseVel;

sArmL.clc.startRew =  1;
[newQ ] = CalcQDirect(sArmL);
% % % [newQ ] = CalcQDirect(sArmL,allQ.qVals{31});
% store q values and attributes
allQ.qVals{33}       = newQ;
allQ.dir(33)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(33)         = sArmL.clc.startRew; % towards
allQ.bodyPart{33}    = 'ArmLeft'; 
allQ.centerpos{33}   = sArmL.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central



% -----------------------------
% NEG REW
% REWARD SIZE
sArmL.clc.startRew =  -1;
[newQ ] = CalcQDirect(sArmL);
% % % [newQ ] = CalcQDirect(sArmL,allQ.qVals{34});
% store q values and attributes
allQ.qVals{34}       = newQ;
allQ.dir(34)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(34)         = sArmL.clc.startRew; % towards
allQ.bodyPart{34}    = 'ArmLeft'; 
allQ.centerpos{34}   = sArmL.clc.nearPos;



%% #################################################################
% =========================================================
% MONKEY HEAD TOWARDS - CONSTRAINED

sHedConstr = sHed;
sHedConstr.clc.actConsequence = [ 0  0  0 ] ; % action 1 stay ONLY

sHedConstr.clc.stimDynams =  @(pos) pos + sHedConstr.clc.baseVel;

% -----------------------------
% POS REW
% REWARD SIZE
sHedConstr.clc.startRew =  1;
[newQ ] = CalcQDirect(sHedConstr);
% % % [newQ ] = CalcQDirect(sHedConstr,allQ.qVals{35});
% store q values and attributes
allQ.qVals{35}        = newQ;
allQ.dir(35)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(35)         = sHedConstr.clc.startRew; % towards
allQ.bodyPart{35}    = 'HeadConstr'; 
allQ.centerpos{35}   = sHedConstr.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW
% REWARD SIZE
sHedConstr.clc.startRew =  -1;
[newQ ] = CalcQDirect(sHedConstr);
% % % [newQ ] = CalcQDirect(sHedConstr,allQ.qVals{36});
% store q values and attributes
allQ.qVals{36}        = newQ;
allQ.dir(36)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(36)         = sHedConstr.clc.startRew; % towards
allQ.bodyPart{36}    = 'HeadConstr'; 
allQ.centerpos{36}   = sHedConstr.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central



%% #################################################################
% =========================================================
% 'MONKEY' HEAD 15 degree angle

sHed15 = sHedConstr;

% Move the data, then shift the Q values next
allQ(37:38,:)       = allQ(35:36,:);


allQ.bodyPart(37:38)    = repmat({'RotatedHead'},[1 2])';

% Shift the data
for iPsi = 37:38

    % shift the q values to the side by 15 degree, i.e. once for every 4
    % rows, starting at the face
    shiftAmount = 1;
% % %     for iR = (sHed15.clc.nearPos(1) + 2) : -3 : 4
% % %         allQ.qVals{iPsi}(:,iR:-1:iR-2, 1:end-shiftAmount ,:) = allQ.qVals{iPsi-30}(:,iR:-1:iR-2, (shiftAmount+1):end ,:);
% % %         shiftAmount = shiftAmount + 1;
% % %     end


for iR = (sHed15.clc.nearPos(1) + 2) : -1 : 1
    allQ.qVals{iPsi}(:,iR, 1:end-shiftAmount ,:) = allQ.qVals{iPsi-2}(:,iR, (shiftAmount+1):end ,:);
    
    if mod(iR - (sHed15.clc.nearPos(1) + 3),4) == 0
        shiftAmount = shiftAmount + 1;
    end
end

            
end



%% Hit probability - body
% -----------------------------
% NOTE: deltaT is 0.5, so I have to halve all the velocities and spreads?
sBdyHP = sBdy;
sBdyHP.clc.actConsequence = [0 0 0]; 
sBdyHP.clc.gammaVal = 1;
sBdyHP.clc.startRew = 1;

% Random stimulus dynamics
rSpr = [-1 0 1];
rSprPr = gaussmf(rSpr,[1 0]) ./ sum(gaussmf(rSpr,[1 0]));
cSpr = [-1 0 1]; %%% cSpr = [ -1 0 1 ];
cSprPr = gaussmf(cSpr,[1 0]) ./ sum(gaussmf(cSpr,[1 0]));
zSpr = [-1 0 1]; %%% zSpr = [ -1 0 1 ];
zSprPr = gaussmf(zSpr,[1 0]) ./ sum(gaussmf(zSpr,[1 0]));

% divide probabilities over 2 because deltaT = 1/2
rSprPr = rSprPr./([2 1 2]) + rSprPr(1).*([0 1 0]); 
cSprPx = cSprPr./([2 1 2]) + cSprPr(1).*([0 1 0]); 
zSprPr = zSprPr./([2 1 2]) + zSprPr(1).*([0 1 0]); 

% x y z, Deterministic stimulus dynamics
sBdyHP.clc.stimDynams =     @(pos) pos + s.clc.baseVel./2; % For approaching, set speed positive
sBdyHP.clc.randSpread =     {rSpr cSpr zSpr}; % Put a little biit of x and z variability in? Kind of arbitrary
sBdyHP.clc.spreadProb =     {rSprPr cSprPr  zSprPr}; %


% Sensory stimulus spreading: from Straka et al ( * 1/2 because their units
% were seconds, and ours are half-seconds
rSnsSpr = [-7:7]; 
cSnsSpr = [-2:2]; 
zSnsSpr = [-2:2]; 

% Set uncertainties to Straka Noel Hoffmann values
rSigma      = sqrt(2.5.^2 + 30.^2) ./ sBdyHP.clc.binW ./2;
rSnsSprPr   = gaussmf(rSnsSpr,[rSigma 0]) ./ sum(gaussmf(rSnsSpr,[rSigma 0])); %sqrt(2.5.^2 + 30.^2)./s.clc.binW
czSigma     = sqrt(5.^2 + 5.^2) ./ sBdyHP.clc.binW ./2;
cSnsSprPr   = gaussmf(cSnsSpr,[czSigma 0]) ./ sum(gaussmf(cSnsSpr,[czSigma 0]));
zSnsSprPr   = gaussmf(zSnsSpr,[czSigma 0]) ./ sum(gaussmf(zSnsSpr,[czSigma 0]));

sBdyHP.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sBdyHP.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};


sBdyHP.clc.baseVel    = [2 0 0]; % Same velocity as Q values, but now in units of half-seconds
sBdyHP.clc.stimDynams = @(pos) pos + sBdyHP.clc.baseVel;


sBdyHP.clc.nReps = 1; 

sBdyHP.clc.stepUpdateFl = 1; % whether to run a full sweep update or not
sBdyHP.clc.nSteps = 1;

% Set the importance of false positives and false negatives
sForUtility.clc.FN = 5;
sForUtility.clc.FP = 1;

bodyHP = calcQDirect(sBdyHP);
cQ = 39;
allQ.qVals{cQ}       = HitProbToUtility(bodyHP,sForUtility);
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sBdyHP.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'TrunkHP'; 
allQ.centerpos{cQ}   = sBdyHP.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


sBdyHP.clc.baseVel    = [-2 0 0]; % Same velocity as Q values, but now in units of half-seconds
sBdyHP.clc.stimDynams = @(pos) pos + sBdyHP.clc.baseVel;
bodyHPaway = calcQDirect(sBdyHP);
cQ = 40;
allQ.qVals{cQ}       = HitProbToUtility(bodyHPaway,sForUtility);
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sBdyHP.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'TrunkHP'; 
allQ.centerpos{cQ}   = sBdyHP.clc.nearPos; 





%% Calculate the field fromthe observer's perspective
% % % % inverse dynamics for P(o|s)
% % % sBdyHPobs = sBdyHP;
% % % sBdyHPobs.clc.sensSpread = cellfun(@(x) -x, sBdyHP.clc.sensSpread, 'UniformOutput', false); 
% % % % no movement, only observation uncertainty allowed
% % % sBdyHPobs.clc.baseVel           = [0 0 0]; 
% % % sBdyHPobs.clc.actConsequence    = [0 0 0]; 
% % % sBdyHPobs.clc.randSpread        = {[0] [0] [0]};
% % % sBdyHPobs.clc.spreadProb        = {[1] [1] [1]};
% % % sBdyHPobs.clc.stimDynams        = @(pos) pos + sBdyHPobs.clc.baseVel;
% % % sBdyHPobs.clc.nSteps            = 1;
% % % 
% % % % hit encoding from observer's perspective
% % % bodyHPobs = calcQDirect(sBdyHPobs,bodyHP);


% % % %% Make Body hitprob, receding stim
% % % % -----------------------------
% % % sBdyHPaway = sBdy;
% % % 
% % % varsForHP = {'actConsequence','gammaVal','startRew','stimDynams',...
% % %              'randSpread','spreadProb','sensSpread','sensProb',...
% % %              'baseVel','stimDynams','nReps','stepUpdateFl','nSteps'}
% % % for iV = 1:numel(varsForHP)
% % %     sBdyHPaway.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
% % % end
% % % 
% % % sBdyHPaway.clc.baseVel    = [-2 0 0]; % Same velocity as Q values, but now in units of half-seconds
% % % sBdyHPaway.clc.stimDynams = @(pos) pos + sBdyHP.clc.baseVel;
% % % 
% % % bodyAway = calcQDirect(sBdyHPaway);



%% Make Head/face hitprob
% -----------------------------
sHedHP = sHed;
varsForHP = {'actConsequence','gammaVal','startRew','stimDynams',...
             'randSpread','spreadProb','sensSpread','sensProb',...
             'baseVel','stimDynams','nReps','stepUpdateFl','nSteps'}
for iV = 1:numel(varsForHP)
    sHedHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
headHP = calcQDirect(sHedHP);
cQ = 41;
allQ.qVals{cQ}       = HitProbToUtility(headHP,sForUtility);
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHedHP.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HeadHP'; 
allQ.centerpos{cQ}   = sHedHP.clc.nearPos; 


sHedHP.clc.baseVel    = [-2 0 0]; % Same velocity as Q values, but now in units of half-seconds
sHedHP.clc.stimDynams = @(pos) pos + sHedHP.clc.baseVel;
hedHPaway = calcQDirect(sHedHP);
cQ = 42;
allQ.qVals{cQ}       = HitProbToUtility(headHPaway,sForUtility);;
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHedHP.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HeadHP'; 
allQ.centerpos{cQ}   = sHedHP.clc.nearPos; 


%% Make hand hitprob
% -----------------------------
sHndHP = sHnd;
for iV = 1:numel(varsForHP)
    sHndHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
handHP = calcQDirect(sHndHP);
cQ = 43;
allQ.qVals{cQ}       = HitProbToUtility(handHP,sForUtility);
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndHP.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestHP'; 
allQ.centerpos{cQ}   = sHndHP.clc.nearPos; 


sHndHP.clc.baseVel    = [-2 0 0]; % Same velocity as Q values, but now in units of half-seconds
sHndHP.clc.stimDynams = @(pos) pos + sHndHP.clc.baseVel;
handHPaway = calcQDirect(sHndHP);
cQ = 44;
allQ.qVals{cQ}       = HitProbToUtility(handHPaway,sForUtility);
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndHP.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestHP'; 
allQ.centerpos{cQ}   = sHndHP.clc.nearPos; 


% ------------------------------------------------
% HAND BY sidE HitPROB
sHndSideHP = sHndHP;
sHndSideHP.clc.nearPos = sHndSide.clc.nearPos;
sHndSideHP.clc.startSC = sHndSide.clc.startSC;


% Move the data, then shift the Q values next
allQ(45:46,:)       = allQ(43:44,:);
allQ.bodyPart(45:46)    = repmat({'HandBySideHP'},[1 2])';

% Shift the data
for iPsi = 45:46
    allQ.qVals{iPsi} = zeros(size(allQ.qVals{iPsi}));

    % FLip the data because it's symetric anyway and the extent on the
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



%% Make Hand tool hitprobs 
% -----------------------------
sHndToolHP = sHndTool;
for iV = 1:numel(varsForHP)
    sHndToolHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
handToolHP = calcQDirect(sHndToolHP);
cQ = 47;
allQ.qVals{cQ}       = HitProbToUtility(handToolHP,sForUtility);
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndToolHP.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolHP'; 
allQ.centerpos{cQ}   = sHndToolHP.clc.nearPos; 


sHndToolHP.clc.baseVel    = [-2 0 0]; % Same velocity as Q values, but now in units of half-seconds
sHndToolHP.clc.stimDynams = @(pos) pos + sHndToolHP.clc.baseVel;
handToolHPaway = calcQDirect(sHndToolHP);
cQ = 48;
allQ.qVals{cQ}       = HitProbToUtility(handToolHPaway,sForUtility);
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndToolHP.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolHP'; 
allQ.centerpos{cQ}   = sHndToolHP.clc.nearPos; 


sHndToolRakeHP = sHndToolRake;
for iV = 1:numel(varsForHP)
    sHndToolRakeHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
handToolRakeHP = calcQDirect(sHndToolHP);
cQ = 49;
allQ.qVals{cQ}       = HitProbToUtility(handToolRakeHP,sForUtility);
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndToolRakeHP.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeHP'; 
allQ.centerpos{cQ}   = sHndToolRakeHP.clc.nearPos; 


sHndToolRakeHP.clc.baseVel    = [-2 0 0]; % Same velocity as Q values, but now in units of half-seconds
sHndToolRakeHP.clc.stimDynams = @(pos) pos + sHndToolRakeHP.clc.baseVel;
handToolRakeHPaway = calcQDirect(sHndToolRakeHP);
cQ = 50;
allQ.qVals{cQ}       = HitProbToUtility(handToolRakeHPaway,sForUtility);
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndToolRakeHP.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeHP'; 
allQ.centerpos{cQ}   = sHndToolRakeHP.clc.nearPos; 



%% Make Hand track hitprob
%-----------------------------
sTrackHP = sTrack;
varsForTrackHP = {'actConsequence','gammaVal','startRew',...
             'sensSpread','sensProb',...
             'baseVel','nReps','stepUpdateFl','nSteps'}
for iV = 1:numel(varsForTrackHP)
    sTrackHP.clc.(varsForTrackHP{iV}) = sBdyHP.clc.(varsForTrackHP{iV});
end
trackHP = calcQDirect(sTrackHP);
cQ = 51;
allQ.qVals{cQ}       = trackHP;
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sTrackHP.clc.startRew; 
allQ.bodyPart{cQ}    = 'HandByChestTrackHP'; 
allQ.centerpos{cQ}   = sTrackHP.clc.nearPos; 


%% Make Monkey hitprobs 
% -----------------------------
sArmHP = sArm;
% For some reason, I made the monkey's arm a tool, but I can't remember why. 
% Oh well, I should take that into account
sArmHP.clc.startSC = sArmHP.clc.toolSC;
sArmHP.clc.startSR = sArmHP.clc.toolSR;
sArmHP.clc.startSZ = sArmHP.clc.toolSZ;
for iV = 1:numel(varsForHP)
    sArmHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
cQ = 52;
allQ.qVals{cQ}       = calcQDirect(sArmHP);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sArmHP.clc.startRew; 
allQ.bodyPart{cQ}    = 'ArmForwardHP'; 
allQ.centerpos{cQ}   = sArmHP.clc.nearPos; 


sArmLHP = sArmL;
% For some reason, I made the monkey's arm a tool, but I can't remember why. 
% Oh well, I should take that into account
sArmLHP.clc.startSC = sArmLHP.clc.toolSC;
sArmLHP.clc.startSR = sArmLHP.clc.toolSR;
sArmLHP.clc.startSZ = sArmLHP.clc.toolSZ;
for iV = 1:numel(varsForHP)
    sArmLHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
cQ = 53;
allQ.qVals{cQ}       = calcQDirect(sArmLHP);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sArmLHP.clc.startRew; 
allQ.bodyPart{cQ}    = 'ArmLeftHP'; 
allQ.centerpos{cQ}   = sArmLHP.clc.nearPos; 


sHedConstrHP = sHedConstr;
for iV = 1:numel(varsForHP)
    sHedConstrHP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
cQ = 54;
allQ.qVals{cQ}       = calcQDirect(sHedConstrHP);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sHedConstrHP.clc.startRew; 
allQ.bodyPart{cQ}    = 'HeadConstrHP'; 
allQ.centerpos{cQ}   = sHedConstrHP.clc.nearPos; 


sHed15HP = sHed15;
for iV = 1:numel(varsForHP)
    sHed15HP.clc.(varsForHP{iV}) = sBdyHP.clc.(varsForHP{iV});
end
cQ = 55;
allQ.qVals{cQ}       = calcQDirect(sHed15HP);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sHed15HP.clc.startRew; 
allQ.bodyPart{cQ}    = 'RotatedHeadHP'; 
allQ.centerpos{cQ}   = sHed15HP.clc.nearPos; 


%% Make multisensory integration model
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
rSnsSprMIPr   = gaussmf(rSnsSprMI,[rSigmaMI 0]) ./ sum(gaussmf(rSnsSprMI,[rSigmaMI 0])); %sqrt(2.5.^2 + 30.^2)./s.clc.binW
czSigmaMI     = sqrt(13.^2 + 11.^2) ./ sBdyMI.clc.binW ./2;
cSnsSprMIPr   = gaussmf(cSnsSprMI,[czSigmaMI 0]) ./ sum(gaussmf(cSnsSprMI,[czSigmaMI 0]));
zSnsSprMIPr   = gaussmf(zSnsSprMI,[czSigmaMI 0]) ./ sum(gaussmf(zSnsSprMI,[czSigmaMI 0]));


sBdyMI.clc.sensSpread = {rSnsSprMI cSnsSprMI zSnsSprMI};
sBdyMI.clc.sensProb   = {rSnsSprMIPr cSnsSprMIPr zSnsSprMIPr};


sBdyMI.clc.baseVel    = [0 0 0]; % Same velocity as Q values, but now in units of half-seconds
sBdyMI.clc.stimDynams = @(pos) pos + sBdyMI.clc.baseVel;


sBdyMI.clc.nReps = 1; 

sBdyMI.clc.stepUpdateFl = 1; % whether to run a full sweep update or not
sBdyMI.clc.nSteps = 1;

cQ = 56;
allQ.qVals{cQ}       = calcQDirect(sBdyMI);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sBdyMI.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'TrunkMI'; 
allQ.centerpos{cQ}   = sBdyMI.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% Make Head/face multsensint
% -----------------------------
sHedMI = sHedHP;
varsForMI = {'stimDynams',...
             'randSpread','spreadProb','sensSpread','sensProb',...
             'baseVel','stimDynams','nReps','stepUpdateFl','nSteps'}
for iV = 1:numel(varsForMI)
    sHedMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ = 57;
allQ.qVals{cQ}       = calcQDirect(sHedMI);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sHedMI.clc.startRew;
allQ.bodyPart{cQ}    = 'HeadMI'; 
allQ.centerpos{cQ}   = sHedMI.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% Make hand multsensint
% -----------------------------
sHndMI = sHndHP;
for iV = 1:numel(varsForMI)
    sHndMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ = 58;
allQ.qVals{cQ}       = calcQDirect(sHndMI);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sHndMI.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestMI';
allQ.centerpos{cQ}   = sHndMI.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% ------------------------------------------------
% Hand by side multisens integration 
sHndSideMI = sHndSideHP;
sHndSideMI.clc.nearPos = sHndSideHP.clc.nearPos;
sHndSideMI.clc.startSC = sHndSideHP.clc.startSC;


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


% Make Hand tool multsensints 
% -----------------------------
sHndToolMI = sHndToolHP;
for iV = 1:numel(varsForMI)
    sHndToolMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ = 60;
allQ.qVals{cQ}       = calcQDirect(sHndToolMI);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sHndToolMI.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolMI';
allQ.centerpos{cQ}   = sHndToolMI.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
sHndToolRakeMI = sHndToolRakeHP;
for iV = 1:numel(varsForMI)
    sHndToolRakeMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ = 61;
allQ.qVals{cQ}       = calcQDirect(sHndToolRakeMI);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sHndToolRakeMI.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeMI';
allQ.centerpos{cQ}   = sHndToolRakeMI.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% Make Hand track multsensint
% -----------------------------
sTrackMI = sTrackHP;
for iV = 1:numel(varsForMI)
    sTrackMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ = 62;
allQ.qVals{cQ}       = calcQDirect(sTrackMI);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sTrackMI.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestTrackMI';
allQ.centerpos{cQ}   = sTrackMI.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% Make Monkey multsensint 
% -----------------------------
sArmMI = sArmHP;
for iV = 1:numel(varsForMI)
    sArmMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ = 63;
allQ.qVals{cQ}       = calcQDirect(sArmMI);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sArmMI.clc.startRew; 
allQ.bodyPart{cQ}    = 'ArmForwardMI'; 
allQ.centerpos{cQ}   = sArmMI.clc.nearPos; 


sArmLMI = sArmLHP;
for iV = 1:numel(varsForMI)
    sArmLMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ = 64; 
allQ.qVals{cQ}       = calcQDirect(sArmLMI);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sArmLMI.clc.startRew; 
allQ.bodyPart{cQ}    = 'ArmLeftMI'; 
allQ.centerpos{cQ}   = sArmLMI.clc.nearPos; 


sHedConstrMI = sHedConstrHP;
for iV = 1:numel(varsForMI)
    sHedConstrMI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ = 65; 
allQ.qVals{cQ}       = calcQDirect(sHedConstrMI);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sHedConstrMI.clc.startRew; 
allQ.bodyPart{cQ}    = 'HeadConstrMI'; 
allQ.centerpos{cQ}   = sHedConstrMI.clc.nearPos; 


sHed15MI = sHed15HP;
for iV = 1:numel(varsForMI)
    sHed15MI.clc.(varsForMI{iV}) = sBdyMI.clc.(varsForMI{iV});
end
cQ = 66;
allQ.qVals{cQ}       = calcQDirect(sHed15MI);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sHed15MI.clc.startRew; 
allQ.bodyPart{cQ}    = 'RotatedHeadMI'; 
allQ.centerpos{cQ}   = sHed15MI.clc.nearPos; 




save('F:\Projects\DPPS\DefenseAgent\Data\NewWithHitProbs_MultSensInt_V3.mat')



%% Make a distance metric, so that I can use it to create exponential and sigmoid functions


% Make Body distance
% -----------------------------
cQ = 67;
allQ.qVals{cQ}       = permute(CalcDistToBody(sBdyMI),[4 1 2 3]);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sBdyMI.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'TrunkDist'; 
allQ.centerpos{cQ}   = sBdyMI.clc.nearPos;

% Make Head distance
% -----------------------------
cQ = 68;
allQ.qVals{cQ}       = permute(CalcDistToBody(sHedMI),[4 1 2 3]);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sHedMI.clc.startRew;
allQ.bodyPart{cQ}    = 'HeadDist'; 
allQ.centerpos{cQ}   = sHedMI.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% Make Hand Chest distance
% -----------------------------
cQ = 69;
allQ.qVals{cQ}       = permute(CalcDistToBody(sHndMI),[4 1 2 3]);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sHndMI.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestDist';
allQ.centerpos{cQ}   = sHndMI.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% Make Hand Side distance
% -----------------------------
cQ = 70;
% I didn't update the starting stim cols, because I calculated Qvalue and
% hitprob differently. But distance is quick, so I can just do it directly
% here
sHndSideMI.clc.startSC = sHndMI.clc.startSC + (sHndSideMI.clc.nearPos(2) - sHndMI.clc.nearPos(2));
allQ.qVals{cQ}       = permute(CalcDistToBody(sHndSideMI),[4 1 2 3]);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sHndSideMI.clc.startRew;
allQ.bodyPart{cQ}    = 'HandBySideDist';
allQ.centerpos{cQ}   = sHndSideMI.clc.nearPos; 


% Make Hand tool distance
% -----------------------------
cQ = 71;
% First include the tool locations in the locations from which to calculate distance
sHndToolDist             = sHndToolMI; 
sHndToolDist.clc.startSC = [sHndToolMI.clc.startSC , sHndToolMI.clc.toolSC];
sHndToolDist.clc.startSR = [sHndToolMI.clc.startSR , sHndToolMI.clc.toolSR];
sHndToolDist.clc.startSZ = [sHndToolMI.clc.startSZ , sHndToolMI.clc.toolSZ];
% Then continue as usual
allQ.qVals{cQ}       = permute(CalcDistToBody(sHndToolDist),[4 1 2 3]);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sHndToolDist.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolDist';
allQ.centerpos{cQ}   = sHndToolDist.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% First include the tool locations in the locations from which to calculate distance
sHndToolRakeDist             = sHndToolRakeMI; 
sHndToolRakeDist.clc.startSC = [sHndToolRakeMI.clc.startSC , sHndToolRakeMI.clc.toolSC];
sHndToolRakeDist.clc.startSR = [sHndToolRakeMI.clc.startSR , sHndToolRakeMI.clc.toolSR];
sHndToolRakeDist.clc.startSZ = [sHndToolRakeMI.clc.startSZ , sHndToolRakeMI.clc.toolSZ];
% Then continue as usual
cQ = 72;
allQ.qVals{cQ}       = permute(CalcDistToBody(sHndToolRakeDist),[4 1 2 3]);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sHndToolRakeDist.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeDist';
allQ.centerpos{cQ}   = sHndToolRakeDist.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% Make Hand track distance
% -----------------------------
cQ = 73;
allQ.qVals{cQ}       = permute(CalcDistToBody(sTrackMI),[4 1 2 3]);
allQ.dir(cQ)         = 1; % towards
allQ.rew(cQ)         = sTrackMI.clc.startRew;
allQ.bodyPart{cQ}    = 'HandByChestTrackDist';
allQ.centerpos{cQ}   = sTrackMI.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% Make Monkey distance
% -----------------------------
cQ = 74;
allQ.qVals{cQ}       = permute(CalcDistToBody(sArmMI),[4 1 2 3]);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sArmMI.clc.startRew; 
allQ.bodyPart{cQ}    = 'ArmForwardDist'; 
allQ.centerpos{cQ}   = sArmMI.clc.nearPos; 

cQ = 75; 
allQ.qVals{cQ}       = permute(CalcDistToBody(sArmLMI),[4 1 2 3]);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sArmLMI.clc.startRew; 
allQ.bodyPart{cQ}    = 'ArmLeftDist'; 
allQ.centerpos{cQ}   = sArmLMI.clc.nearPos; 

cQ = 76; 
allQ.qVals{cQ}       = permute(CalcDistToBody(sHedConstrMI),[4 1 2 3]);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sHedConstrMI.clc.startRew; 
allQ.bodyPart{cQ}    = 'HeadConstrDist'; 
allQ.centerpos{cQ}   = sHedConstrMI.clc.nearPos; 

cQ = 77;
allQ.qVals{cQ}       = permute(CalcDistToBody(sHed15MI),[4 1 2 3]);
allQ.dir(cQ)         = 1; 
allQ.rew(cQ)         = sHed15MI.clc.startRew; 
allQ.bodyPart{cQ}    = 'RotatedHeadDist'; 
allQ.centerpos{cQ}   = sHed15MI.clc.nearPos; 


save('F:\Projects\DPPS\DefenseAgent\Data\NewWithHitProbs_MultSensInt_AndDist.mat')


%%
% % % %% #################################################################
% % % % ===============================================s==========
% % % % HAND by NEAR side [resting position], all directions 
% % % 
% % % % $$$ --> NEXT FIGURE OUT HOW TO FLIP THIS BADBOY so that the useful data
% % % % ends up on the LEFT of the hand --> bigger space to use :D
% % % 
% % % sHndRest = sHnd;
% % % sHndRest.clc.nearPos = sHnd.clc.nearPos + [0, 4, -5]';
% % % 
% % % 
% % % % Move the data, then shift the Q values next
% % % allQ(25:28,:)       = allQ(9:12,:);
% % % allQ.bodyPart(25:28)    = repmat({'HandRest'},[1 4])';
% % % 
% % % % $$$ CLEAN THIS UP w.r.t. new hndside.clc.nearpos defined above
% % % handMidPos      = sHnd.clc.nearPos(2);
% % % handNewPos      = sHndRest.clc.nearPos(2);
% % % 
% % % handMidPosUD    = sHnd.clc.nearPos(3);
% % % handNewPosUD    = sHndRest.clc.nearPos(3);
% % % 
% % % sideHandShift   = handNewPos - handMidPos;
% % % udHandShift     = handNewPosUD - handMidPosUD;
% % % 
% % % % Find transposable width on the narrower side
% % % moveWidth   = s.wrld.size(2) - handNewPos;
% % % moveHeight  = s.wrld.size(3) - handNewPosUD;
% % % 
% % % % Shift the data
% % % for iPsi = 25:28
% % %     allQ.qVals{iPsi} = zeros(size(allQ.qVals{iPsi}));
% % % % % %     allQ.qVals{iPsi}(:, :, sideHandShift : end,:) = ...
% % % % % %         allQ.qVals{iPsi-4}(:, :, 1 : moveWidth + handMidPos + 1,:);
% % % % % % 
% % % % % %     allQ.centerpos{iPsi}   = sHnd.clc.nearPos;
% % % 
% % %     % FLip the data left rightbecause it's symetric anyway and the extent on the
% % %     % right is bigger than on the left
% % %     leftCopyDist  = min([handMidPos - 1 , size(    allQ.qVals{iPsi},3) - handNewPos]);
% % %     rightCopyDist = min([handNewPos - 1 , size(    allQ.qVals{iPsi},3) - handMidPos]);
% % % 
% % %     upCopyDist    = min([size(allQ.qVals{iPsi},4) - handMidPosUD , size(allQ.qVals{iPsi},4) - handNewPosUD]);
% % %     downCopyDist  = min([handMidPosUD - 1, handNewPosUD - 1]);
% % % 
% % %     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % %     % FIRST DATA BELOW HAND
% % %     % Move the data from the right of the central hand to the left of the
% % %     % lateral hand, and flip it
% % %     allQ.qVals{iPsi}(:, :, handNewPos + [-rightCopyDist:0]   , [-downCopyDist:0] + handNewPosUD ) = ...
% % %         flip(allQ.qVals{iPsi-16}(:, :, handMidPos + [0:rightCopyDist] ,[-downCopyDist:0] + handMidPosUD ),3);
% % % 
% % %     % Move the data from the left of the central hand to the right of the
% % %     % lateral hand, and flip it
% % %     allQ.qVals{iPsi}(:, :, handNewPos + [1:leftCopyDist],[-downCopyDist:0] + handNewPosUD) = ...
% % %         flip(allQ.qVals{iPsi-16}(:, :, handMidPos + [-leftCopyDist:-1] ,[-downCopyDist:0] + handMidPosUD ),3);
% % %     
% % %     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % %     % THEN DATA ABOVE HAND
% % %     % Move the data from the right of the central hand to the left of the
% % %     % lateral hand, and flip it
% % %     allQ.qVals{iPsi}(:, :, handNewPos + [-rightCopyDist:0]   , [1:upCopyDist] + handNewPosUD ) = ...
% % %         flip(allQ.qVals{iPsi-16}(:, :, handMidPos + [0:rightCopyDist] ,[1:upCopyDist]+ handMidPosUD ),3);
% % % 
% % %     % Move the data from the left of the central hand to the right of the
% % %     % lateral hand, and flip it
% % %     allQ.qVals{iPsi}(:, :, handNewPos + [1:leftCopyDist],[1:upCopyDist] + handNewPosUD) = ...
% % %         flip(allQ.qVals{iPsi-16}(:, :, handMidPos + [-leftCopyDist:-1] ,[1:upCopyDist] + handMidPosUD ),3);
% % %     
% % %     
% % %     
% % %     allQ.centerpos{iPsi}   = sHnd.clc.nearPos;
% % % 
% % %     allQ.centerpos{iPsi} =  sHndRest.clc.nearPos;
% % %         
% % % end


% % % %% -------------------------------------------------------------------------
% % % % FOR LUCA effects of movement freedom
% % % 
% % % % =========================================================
% % % % BODY TOWARDS
% % % 
% % % s.wrld.size = [27 27 31];
% % % 
% % % 
% % % s.clc.gammaVal   =     0.8;
% % % % x y z, Deterministic stimulus dynamics
% % % s.clc.stimDynams =     @(pos) pos + [0 0 0]; % For approaching, set speed positive
% % % s.clc.randSpread =     {[0 1] [-1   0  1] [-1  0  1]}; % Put a little biit of x and z variability in? Kind of arbitrary
% % % s.clc.spreadProb =     {[.5 .5] [0.3 .4 0.3]  [ .3 .4 .3]}; % x y z, probabilities of spread
% % % % NOTE: If the limb can move at different speeds, BOTH speeds have to be
% % % % put in as potential actions, because the model doesn't have any collision
% % % % calculation that takes into account overshoot
% % % 
% % % % FIRST CALCULATE FOR BODY - moves more slowly than limb, let's say. Also
% % % % ONLY has negative potential rewards
% % % sBdy = s;
% % % sBdy.clc.actConsequence = [ 0  0  0 ; ... % action 1 stay
% % %     0  1  0 ; ... % action 2 left
% % %     0 -1  0 ; ... % action 3 right
% % %     0  0  1 ; ... % action 4 up
% % %     0  0 -1 ; ... % action 5 down
% % %     -1  0  0 ; ... % action 6 forward
% % %     1  0  0 ];    % action 7 backward-
% % % 
% % % % Location and size of limb, and where the Q-values are calculated FROM
% % % sBdy.clc.startSR = []; sBdy.clc.startSC = []; sBdy.clc.startSZ = [];
% % % bdyZs = [8:18];
% % % bdyCs = [7:15];
% % % bdyRs = [21 20 20 19 19 19 20 20 21];
% % % iRew = 0; % Initialise rewarded block counter
% % % for  iZ = 1:length(bdyZs)
% % %     for  iC = 1:length(bdyCs)
% % %         iRew = iRew + 1;
% % %         sBdy.clc.startSR(iRew) =  bdyRs(iC);
% % %         sBdy.clc.startSC(iRew) =  bdyCs(iC);
% % %         sBdy.clc.startSZ(iRew) =  bdyZs(iZ);
% % %     end
% % % end
% % % 
% % % % -----------------------------
% % % % MANY movements, NEG reward
% % % sBdy.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sBdy);
% % % % store q values and attributes
% % % allQ.qVals{20}        = newQ;
% % % allQ.dir(20)         = 1; %s.clc.stimDynams([0 0 0]); % towards
% % % allQ.rew(20)         = sBdy.clc.startRew; % towards
% % % allQ.bodyPart{20}        = 'TrunkLUCA'; 
% % % allQ.centerpos{20}   = [0 0 0]; % x y z position of this plot in the OVERALL space --> trunk is central
% % % allQ.moveOpts{20}    = sBdy.clc.actConsequence;
% % % 
% % % 
% % % % -----------------------------
% % % % ONE movement, NEG reward
% % % sBdy.clc.actConsequence = [ 0  0  0  ];    % STAY
% % % [newQ ] = CalcQDirect(sBdy);
% % % % store q values and attributes
% % % allQ(21,:)        = allQ(20,:);
% % % allQ.qVals{21}    = newQ;
% % % allQ.moveOpts{21}    = sBdy.clc.actConsequence;
% % % 
% % % 
% % % 
% % % % =========================================================
% % % % BODY TOWARDS
% % % 
% % % sHed = sBdy;
% % % sHed.clc.actConsequence = [ 0  0  0 ; ... % action 1 stay
% % %     0  1  0 ; ... % action 2 left
% % %     0 -1  0 ; ... % action 3 right
% % %     0  0  1 ; ... % action 4 up
% % %     0  0 -1 ; ... % action 5 down
% % %     -1  0  0 ; ... % action 6 forward
% % %     1  0  0 ];    % action 7 backward-
% % % % Location and size of limb, and where the Q-values are calculated FROM
% % % sHed.clc.startSR = []; sHed.clc.startSC = []; sHed.clc.startSZ = [];
% % % hedRs = [20 19 20];
% % % hedCs = [10:12];
% % % hedZs = [22:26];
% % % iRew = 0; % Initialise rewarded block counter
% % % for  iZ = 1:length(hedZs)
% % %     for  iC = 1:length(hedCs)
% % %         iRew = iRew + 1;
% % %         sHed.clc.startSR(iRew) =  hedRs(iC);
% % %         sHed.clc.startSC(iRew) =  hedCs(iC);
% % %         sHed.clc.startSZ(iRew) =  hedZs(iZ);
% % %     end
% % % end
% % % 
% % % % -----------------------------
% % % % MANY movements, NEG reward
% % % sHed.clc.startRew =  -1;
% % % [newQ ] = CalcQDirect(sHed);
% % % % store q values and attributes
% % % allQ.qVals{22}        = newQ;
% % % allQ.dir(22)         = 1; %s.clc.stimDynams([0 0 0]); % towards
% % % allQ.rew(22)         = sHed.clc.startRew; % towards
% % % allQ.bodyPart{22}        = 'HeadLUCA'; 
% % % allQ.centerpos{22}   = [0 0 0]; % x y z position of this plot in the OVERALL space --> trunk is central
% % % allQ.moveOpts{22}    = sHed.clc.actConsequence;
% % % 
% % % % -----------------------------
% % % % ONE movement, NEG reward
% % % sHe d.clc.actConsequence = [ 0  0  0  ];    % STAY
% % % [newQ ] = CalcQDirect(sHed);
% % % % store q values and attributes
% % % allQ(23,:)        = allQ(22,:);
% % % allQ.qVals{23}    = newQ;
% % % allQ.moveOpts{23}    = sHed.clc.actConsequence;

%% Plot 2D 

% iAct = 1:13;


% iAct = 1;
% newQ = allQ.qVals{1} ;%+ allQ.qVals{9}; % Just body

% newQ = allQ.qVals{1}  + allQ.qVals{5}; % Head plus body, positive reward, towards
% newQ = allQ.qVals{2}  + allQ.qVals{6}; % Head plus body, negative reward, towards
% newQ = newQ(iAct,:,:,:);

% newQ = allQ.qVals{1} + allQ.qVals{5} + allQ.qVals{9}(1:7,:,:,:); % Body plus Head plus hand, positive reward, towards

% newQ = allQ.qVals{2} + allQ.qVals{6} + allQ.qVals{15}(1:7,:,:,:); % Body plus Head NEGATIVE plus hand onside POSITIVE

% newQ = allQ.qVals{3} + allQ.qVals{7} + allQ.qVals{9}(13:end,:,:,:); % AWAY, hand by side


% newQ = allQ.qVals{13} ;%+ allQ.qVals{9}; % The two hand positions
% sTmp = sHndSide;

% newQ = allQ.qVals{11} ;%+ allQ.qVals{9}; % The two hand positions
% sTmp = sHnd;

% newQ = allQ.qVals{25} ;%+ allQ.qVals{9}; % The two hand positions
% sTmp = sHndRest;
% iAct = 1:43;

% newQ = allQ.qVals{6} ;%+ allQ.qVals{9}; % The two hand positions
% sTmp = sHed;
% iAct = 1:7;


% newQ = allQ.qVals{35}   ;
% sTmp = sHedConstr;
% iAct = 1;

% % % % % newQ = allQ.qVals{37}  - allQ.qVals{35} ;
% % % newQ = allQ.qVals{38}   ;
% % % sTmp = sHed15;
% % % iAct = 1;


% % newQ = allQ.qVals{37}  - allQ.qVals{35} ;
% % newQ = allQ.qVals{38}   ;


% -----
% newQ = bodyHP;
% sTmp = sBdyHP;
% iAct = 1;

% newQ = bodyUtilityHP;
% sTmp = sBdyHP;
% iAct = 1;


% newQ = allQ.qVals{82}   ;
% sTmp = sHed;
% iAct = 1;


% newQ = allQ.qVals{90}   ;
% sTmp = sHndSide;
% iAct = 1;

newQ = allQ.qVals{86}   ;
sTmp = sHnd;
iAct = 1;

% newQ = allQ.qVals{43}   ;
% sTmp = sHndHP;
% iAct = 1;

% newQ = allQ.qVals{54}   ;
% sTmp = sHedConstrHP;
% iAct = 1;


% newQ = allQ.qVals{65}   ;
% sTmp = sHedConstrMI;
% iAct = 1;


% % newQ = allQ.qVals{37};
% newQ = allQ.qVals{101};
% sTmp = sHedConstr;
% iAct = 1;


% newQ = permute(handSideDist,[4 1 2 3]) ;
% sTmp = sHndSide;


% newQ = bodyHPaway;
% sTmp = sBdyHP;
% iAct = 1;

% newQ = handHP;
% sTmp = sHndHP;
% iAct = 1;

% 
% newQ = headHP;
% sTmp = sHedHP;
% iAct = 1;

% newQ = bodyHPobs;
% sTmp = sBdyHPobs;
% iAct = 1;
% -----



% newQ = allQ.qVals{2} ;%+ allQ.qVals{9}; % The two hand positions
% sTmp = sBdy;
% iAct = 1:7;

% newQ = allQ.qVals{27} ;%+ allQ.qVals{9}; % The two hand positions
% sTmp = sHndToolRake;
% iAct = 1:43;

% newQ = allQ.qVals{31} ;
% sTmp = sArm;
% iAct = 1:7;

% newQ = allQ.qVals{33} ;
% sTmp = sArmL;
% iAct = 1:7;

newQ2 = max(newQ(iAct,:,:,:),[],1);
% newQ2 = mean(newQ(iAct,:,:,:),1);

% % % % FOR LUCA: body and head, all movements available
% % % newQ = allQ.qVals{20}(1:7,:,:,:) + allQ.qVals{22}(1:7,:,:,:); 
% % % % % % % FOR LUCA: body and head, NO movements available
% % % % % % newQ = allQ.qVals{21} + allQ.qVals{23}; 


f.TwoDimPlot.f        = figure('Position',[20 20 1200 600]);
f.TwoDimPlot.ax{1}    = axes('Position',[.05 .1 .4 .8]);

imagesc(squeeze(-newQ2(1,:,:,sTmp.clc.nearPos(3) + 1 )) );
% colormap(whitetocol(100,[0 0 0.7]))
colormap(redbluecmapRory)
GridOverImage(fS,f.TwoDimPlot.ax{1});
caxis([-1 1])
colorbar


f.TwoDimPlot.ax{2}    = axes('Position',[.55 .1 .4 .8]);
imagesc(squeeze(-newQ2(1,:,sTmp.clc.nearPos(2),:))' );
hold on; axis xy
% colormap(whitetocol(100,[0 0 0.7]))
colormap(redbluecmapRory)
GridOverImage(fS,f.TwoDimPlot.ax{2});
caxis([-1 1])


%  3D plot
figure('Position',[500 20 600 600])
tmpQ = -newQ2;
% % % tmpQ(abs(tmpQ) < prctile(abs(newQ(:)),90) ) = 0;
tmpQ(abs(tmpQ) < 0.2 ) = 0;
voxelSurf(squeeze(tmpQ(1,:,:,:)),true,[0.1 0.8 0.1 0.8 0.1 0.8],.3);
% voxelSurf(squeeze(tmpQ(1,:,:,:)),true,[0.8 0.1 0.8 0.1 0.8 0.1 ],.4);
% colormap(whitetocol(100,[0 0 0.7]))
colormap(redbluecmapRory)
view([-2 2 2])
caxis([-1 1])
colorbar


figure,plot(squeeze(-newQ2(:,:,sTmp.clc.nearPos(2),sTmp.clc.nearPos(3)))' )
% figure,plot(squeeze(-newQ2(:,:,sHnd.clc.nearPos(2),sHnd.clc.nearPos(3)))' )
% figure,plot(squeeze(-newQ2(:,:,sHndSide.clc.nearPos(2),sHndSide.clc.nearPos(3)))' )

%% #################################################################
% Fit all data: [NOTE: make sure that data from the same experiment is
% always sequential; i.e. d, = exp1, d , = exp 2 etc]
% CREATE data table

sBdy.clc.binW    = 5;
sHnd.clc.binW    = 5;
sHndSide.clc.binW= 5;
sHed.clc.binW    = 5;


% % % sBdy.clc.binW    = 4.5;
% % % sHnd.clc.binW    = 4.5;
% % % sHndSide.clc.binW= 4.5;
% % % sHed.clc.binW    = 4.5;

d = table;

% Put row counter first so it's easier to eyeball things
d.iD(1) = 1;


% % % d.psiSettings{1,1} = 'max_all'; %'average'; 'average_all';
% % % d.psiSettings{1,2} = 'no_norm'; %'diff_from_mean_divmean_norm'; %'diff_from_mean_norm'; % 'diff_from_mean_norm' ; %'sum_norm'; %'no_norm';

d.psiSettings{1,1} = 'max_then_avg'; %'average'; 'average_all';
d.psiSettings{1,2} = 'no_norm'; %'diff_from_mean_divmean_norm'; %'diff_from_mean_norm'; % 'diff_from_mean_norm' ; %'sum_norm'; %'no_norm';

% % % d.psiSettings{1,1} = 'average_all'; %'average'; 'average_all';
% % % d.psiSettings{1,2} = 'no_norm'; %'diff_from_mean_divmean_norm'; %'diff_from_mean_norm'; % 'diff_from_mean_norm' ; %'sum_norm'; %'no_norm';




% d.psiSettings{1,1} = 'average_all'; %'average'; 'average_all';
% d.psiSettings{1,2} = 'diff_norm'; %'diff_from_mean_divmean_norm'; %'diff_from_mean_norm'; % 'diff_from_mean_norm' ; %'sum_norm'; %'no_norm';

% d.psiSettings{1,1} = 'max_then_avg'; %'average'; 'average_all';
% d.psiSettings{1,2} = 'rel_still_norm'; %'diff_from_mean_divmean_norm'; %'diff_from_mean_norm'; % 'diff_from_mean_norm' ; %'sum_norm'; %'no_norm';

% d.psiSettings{1,1} = 'average'; %'average'; 'average_all';
% d.psiSettings{1,2} = 'no_norm'; %'diff_from_mean_divmean_norm'; %'diff_from_mean_norm'; % 'diff_from_mean_norm' ; %'sum_norm'; %'no_norm';


d.psiSettings{1,3} = '';



% ==============================
% APPROACHING TRUNK 1
cD              = 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';
d.cmPos{cD}      = [ 5, 24, 43, 62, 81, 100 ; ...
                    0,  0,  0,  0,  0,   0 ; ...
                    0,  0,  0,  0,  0,   0  ];
d = AddbinPos(d,sBdy,cD); % approaching body

% Create real values
d.realDat{cD} = [0.646 0.755 0.568 0.313 0.028 -0.213] .* 73.39449541;
d.realSTE{cD} = abs([0.474 0.530 0.396 0.038 -0.227 -0.484] .* 73.39449541 - d.realDat{cD})  ;
d.exp(cD)     = 1;

% Create 'successor representation'
% For modelling trunk towards, we take the trunk, head, and hand-by-side
% Average over body-part related q values
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{cD}}}; 

% ==============================
% RECEDING TRUNK 2
cD               = cD + 1;
d(cD,:)          = d(cD-1,:);
d.dir(cD)        = -1;
d.realDat{cD}    = [0.183 0.211 0.364 0.479 0.394 0.359] .* 73.39449541;
d.realSTE{cD}    = abs([0.069 -0.039 0.131 0.252 0.207 0.174] .* 73.39449541 - d.realDat{cD})  ;
d.exp(cD)        = 1;

d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{-1}}};

% ==============================
% APPROACHING HAND ONLY 3
cD               = cD + 1;
d(cD,:)          = d(cD-1,:);
d.dir(cD)        = 1;
d.tactloc{cD}    = 'HandSide';
d.tact{cD}       = 'Hand';
d.cmPos{cD}      = [ 5  27, 50, 71, 93 ; ...
                    0,  0,  0,  0,  0, ; ...
                    0,  0,  0,  0,  0, ];
d = AddbinPos(d,sHndSide,cD); % approaching hand on side

% Create real values
d.realDat{cD} = [1.044 0.969 0.547 0.285 -0.234] .* 45.62737643;
d.realSTE{cD} = abs([0.831 0.711 0.231 -0.021 -0.516] .* 45.62737643 - d.realDat{cD})  ;
d.exp(cD)     = 2;


d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{1}}}; 

% ==============================
% RECEDING HAND ONLY 4
cD            = cD + 1;
d(cD,:)       = d(cD - 1,:);
d.dir(cD)     = -1;
d.realDat{cD} = [0.662 0.624 0.360 0.175 0.208] .* 45.62737643;
d.realSTE{cD} = abs([0.409 0.372 0.078 0.014 0.012] .* 45.62737643 - d.realDat{cD})  ;
d.exp(cD)     = 2;

% Settings for 'successor representation'
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{-1}}};





% #################################################################
% ==============================
% APPROACHING TRUNK WITH HAND ON TRUNK, TRUNK STIMULATED 5
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';
d.cmPos{cD}      = [ 5, 24, 43, 62, 81, 100 ; ...
                    0,  0,  0,  0,  0,   0 ; ...
                    0,  0,  0,  0,  0,   0  ];
d = AddbinPos(d,sBdy,cD); % approaching body

% Create real values
d.realDat{cD} = [0.983 0.865 0.751 0.490 0.440 0.046] .* 56.02240896;
d.realSTE{cD} = abs([0.777 0.508 0.404 0.095 0.182 0.334] .* 56.02240896 - d.realDat{cD})  ;
d.exp(cD)     = 3;

d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChest'}}, ...
    {'dir' ,{1}}}; 

% ==============================
% APPROACHING HAND WITH HAND ON TRUNK, HAND STIMULATED 6
cD               = cD + 1; handOnTrunkApproachHandStimD = cD;
d(cD,:)          = d(cD - 1,:);
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'Hand';

% ==============================
% RECEDING TRUNK WITH HAND ON TRUNK, TRUNK STIMULATED 7
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.dir(cD)        = -1;
d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';

% Create real values
d.realDat{cD} = [0.425 0.532 0.579 0.581 0.702 0.362] .* 56.02240896;
d.realSTE{cD} = abs([0.148 0.200 0.219 0.184  0.344 0.012] .* 56.02240896 - d.realDat{cD})  ;
d.exp(cD)     = 3;

d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChest'}}, ...
    {'dir' ,{-1}}}; 

% ==============================
% RECEDING HAND WITH HAND ON TRUNK, HAND STIMULATED 8
cD               = cD + 1;
d(cD,:)          = d(cD-1,:);
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'Hand';

% #################################################################
% ==============================
% APPROACHING BETWEEN HAND AND TRUNK WITH HAND STIMULATED AND HAND AWAY
% FROM CHEST 9

cD = 8

cD               = cD + 1;
d(cD,:)          = d(handOnTrunkApproachHandStimD,:);

% auditory stimulus position offset to the participants' right

% 15 Probably works, BUT because the hand and body OVERLAP when the hand is
% on the chest, while they don't overlap when it isn't --> that's kinda the
% opposite to what NOEL et al's explanation is.
d.cmPos{cD}      = [ 5, 24, 43, 62, 81, 100 ; ...
                     0,  0,  0,  0,  0,   0 ; ...
                     0,  0,  0,  0,  0,   0  ] + [0 20 0]';

% This one could work, BUT it requires either more expansive lateral hand
% movement, OR a change in tactile contribution to the psi multipliers
% between hand at chest and hand by side
% d.cmPos{cD}      = [ 5, 24, 43, 62, 81, 100 ; ...
%                      0,  0,  0,  0,  0,   0 ; ...
%                      0,  0,  0,  0,  0,   0  ] + [0 25 0]';

d = AddbinPos(d,sBdy,cD); % approaching body [ but to the side, as defined above]

d.tactloc{cD}    = 'HandSide';
d.tact{cD}       = 'Hand';

d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{1}}}; 

% Create real values <-- DO THIS DO THIS
d.realDat{cD} = [0.566 0.665 0.046 0.089 -0.303 -0.647] .* 56.11672278;
d.realSTE{cD} = abs([0.415 0.536 -0.143 -0.035 -0.447 -0.839] .* 56.11672278 - d.realDat{cD})  ;
d.exp(cD)     = 4;

% ==============================
% APPROACHING BETWEEN HAND AND TRUNK WITH HAND STIMULATED AND HAND NEAR/ON
% CHEST 10
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);

d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'Hand';

d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChest'}}, ...
    {'dir' ,{1}}}; 

% Create real values <-- DO THIS DO THIS
d.realDat{cD} = [0.730 0.694 0.412 -0.006 -0.227 -0.343] .* 56.11672278;
d.realSTE{cD} = abs([0.583 0.538 0.240 -0.182 -0.472 -0.536] .* 56.11672278 - d.realDat{cD})  ;
d.exp(cD)     = 4;

% #################################################################
% ==============================
% APPROACHING FACE WITH FACE STIMULATED 11
cD               = cD + 1;  faceApproachFaceStimD = cD;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'Head';
d.cmPos{cD}      = [ 5, 37, 69, 101, 133, 165, 197 ; ...
                     0,  0,  0,   0,   0,   0,   0 ; ...
                     0,  0,  0,   0,   0,   0,   0];

d = AddbinPos(d,sHed,cD); % approaching head

% Create real values 
d.realDat{cD} = [1.043 0.785 0.566 0.599 0.493 0.543 0.424] .* 63.84676776;
d.realSTE{cD} = abs([0.855 0.629 0.409 0.450 0.331 0.372 0.242] .* 63.84676776 - d.realDat{cD})  ;

d.exp(cD)     = 5;

d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{1}}}; 

% ============================== 12
% APPROACHING BODY WITH FACE STIMULATED
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'Head';
d = AddbinPos(d,sBdy,cD); % approaching body

% Create real values
d.realDat{cD} = [0.38 0.366 0.463 0.314 0.446 0.314 0.394] .* 63.84676776;
d.realSTE{cD} = abs([0.131 0.129 0.252 0.106 0.214 0.084 0.190] .* 63.84676776 - d.realDat{cD}) ;

% ==============================
% APPROACHING BODY WITH BODY STIMULATED 13
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';
d = AddbinPos(d,sBdy,cD); % approaching body

% Create real values
d.realDat{cD} = [0.959 0.801 0.628 0.433 0.440 0.438 0.358] .* 51.02040816;
d.realSTE{cD} = abs([0.818 0.649 0.467 0.283 0.243 0.259 0.215] .* 51.02040816 - d.realDat{cD})  ;

d.exp(cD)     = 6;

% ==============================
% APPROACHING FACE WITH BODY STIMULATED 14
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';
d = AddbinPos(d,sHed,cD); % approaching Head

% Create real values
d.realDat{cD} = [0.609 0.539 0.354 0.376 0.189 0.251 0.194] .* 51.02040816;
d.realSTE{cD} = abs([0.351 0.386 0.117 0.120 -0.033 0.006 -0.020] .* 51.02040816 - d.realDat{cD})  ;


% #################################################################
% ==============================
% VISUAL APPROACHING FACE WITH TRUNK STIMULATED 15
cD               = cD + 1;
d(cD,:)          = d(faceApproachFaceStimD ,:); 

% Put position a bit more inbetween head and body for the visual one
d.cmPos{cD}      = [ 5, 37, 69, 101, 133, 165, 197 ; ...
                     0,  0,  0,   0,   0,   0,   0 ; ...
                   -10,-10,-10, -10, -10, -10, -10];

d.tactloc{cD}    = 'Trunk';
d.tact{cD}       = 'Trunk';
d = AddbinPos(d,sHed,cD); % approaching Head [ but slightly lower because of cmpos above]

% Create real values
d.realDat{cD} = [1.278 1.000 0.639 0.531 0.380 0.295 -0.073] .* 48.95104895;
d.realSTE{cD} = abs([1.118 0.838 0.363 0.335 0.222 0.085 -0.271] .* 48.95104895 - d.realDat{cD})  ;

d.exp(cD)     = 7;

d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandBySide'}}, ...
    {'dir' ,{1}}}; 

% ==============================
% VISUAL APPROACHING FACE WITH FACE STIMULATED 16
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'Head';

% Create real values
d.realDat{cD} = [1.096 0.948 0.092 0.356 0.023 -0.083 -0.114] .* 48.95104895;
d.realSTE{cD} = abs([0.944 0.794 -0.184 0.169 -0.135 -0.300 -0.305] .* 48.95104895 - d.realDat{cD})  ;


% #################################################################
% Fit holmes data: [NOTE: make sure that data from the same experiment is
% always sequential; i.e. d, = exp1, d , = exp 2 etc]
% CREATE data table

% ==============================
% NO TOOL for TOOL EXPERIMENT 17
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'HandCCE';
d.cmPos{cD}      = [5, 45, 70 ; ...
                    0,  0,  0 ; ...
                    0,  0,  0 ];
d = AddbinPos(d,sBdy,cD); % approaching body

% Create real values
d.realDat{cD} = [1.714 0.641 0.342] .* 87.95074758;
% Transform from 95% confidence interval to SE
d.realSTE{cD} = ((abs([2.104 0.972 0.677] .* 87.95074758 - d.realDat{cD})) ./ 47.5) .* 34.1 ;
d.exp(cD)     = 8;

d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChest'}}}; 



% ==============================
% YES TOOL for TOOL EXPERIMENT 18
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);

% Create real values
d.realDat{cD} = [1.587 0.410 0.799] .* 87.95074758;
% Transform from 95% confidence interval to SE
d.realSTE{cD} = ((abs([1.911 0.675 1.071] .* 87.95074758 - d.realDat{cD})) ./ 47.5) .* 34.1  ;

d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChestPlusTool'}}}; 

% #################################################################
% Fit Huijsmans data: spider on a track
% CREATE data table

% ==============================
% BUTTERFLY on track towards hand 19
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'HandSpider';
d.cmPos{cD}      = [ 5, 20, 35, 50, 65, 80 ; ...
                    10, 10, 10,  0,  0,  0 ; ...
                     0,  0,  0, 0,  0,  0 ] + [0 0 0]';
d = AddbinPos(d,sTrack,cD); % approaching hand on chest

% Create real values
d.realDat{cD} = [1.487 1.830 2.014 1.924 2.754 2.935].* 19.17545542;
d.realSTE{cD} = abs([0.735 1.047 1.256 1.148 1.961 2.172] .* 19.17545542 - d.realDat{cD}) ;
d.realDat{cD} = fliplr(d.realDat{cD}); d.realSTE{cD} = fliplr(d.realSTE{cD}); % accidentally pasted them the wrong way round...
d.exp(cD)     = 9;


d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChestTrack'}} , ...
    {'rew',{1,-1}} }; 

% ==============================
% BUTTERFLY on track away from hand 20
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.cmPos{cD}      = [ 5, 20, 35, 50, 65, 80 ; ...
                    -10, -10, -10,  0,  0,  0 ; ...
                     0,  0,  0, 0,  0,  0 ] + [0 0 0]';
d = AddbinPos(d,sTrack,cD); % approaching hand on chest

% Create real values
d.realDat{cD} = [1.487 1.830 2.014 1.896 2.292 2.427].* 19.17545542;
d.realSTE{cD} = abs([0.735 1.047 1.256 1.092 1.470 1.637] .* 19.17545542 - d.realDat{cD}) ;
d.realDat{cD} = fliplr(d.realDat{cD}); d.realSTE{cD} = fliplr(d.realSTE{cD}); % accidentally pasted them the wrong way round...
d.exp(cD)     = 9;


% ==============================
% SPIDER on track towards hand 21
cD               = cD + 1;
d(cD,:)          = d(cD - 2,:);
d.psiSplit{cD} = { {'bodyPart',{'Trunk','Head','HandByChestTrack'}} , ...
    {'rew',{1,-1.25}} }; 

% Create real values
d.realDat{cD} = [1.105 1.730 2.129 2.219 3.282 3.164] .* 19.17545542;
d.realSTE{cD} = abs([0.322 0.926 1.391 1.478 2.533 2.27] .* 19.17545542 - d.realDat{cD}) ;
d.realDat{cD} = fliplr(d.realDat{cD}); d.realSTE{cD} = fliplr(d.realSTE{cD}); % accidentally pasted them the wrong way round...
d.exp(cD)     = 9;


% ==============================
% SPIDER on track away from hand 22
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
d.cmPos{cD}      = [ 5, 20, 35, 50, 65, 80 ; ...
                    -10, -10, -10,  0,  0,  0 ; ...
                     0,  0,  0, 0,  0,  0 ] + [0 0 0]';
d = AddbinPos(d,sTrack,cD); % approaching hand on chest

% Create real values
d.realDat{cD} = [1.105 1.730 2.129 1.986 2.608 3.045] .* 19.17545542;
d.realSTE{cD} = abs([0.322 0.926 1.391 1.224 1.797 2.161] .* 19.17545542 - d.realDat{cD}) ;
d.realDat{cD} = fliplr(d.realDat{cD}); d.realSTE{cD} = fliplr(d.realSTE{cD}); % accidentally pasted them the wrong way round...
d.exp(cD)     = 9;





% #################################################################
% Fit Wamain data: mu power during VR object judgements
% CREATE data table

cD = 22
% % % % ==============================
% % % % PROTOTYPICAL OBJECT judgements 23
% % % cD               = cD + 1;
% % % d.dir(cD)        = 1;
% % % d.tactloc{cD}    = 'HandSide';
% % % d.tact{cD}       = 'HandReachEEG';
% % % d.cmPos{cD}      = [20 87.5 160 ; ...
% % %                      0,  0,   0 ; ...
% % %                      0,  0,   0 ];
% % % 
% % % d = AddbinPos(d,sHndSide,cD); % approaching hand on chest
% % % 
% % % % Create real values
% % % tmpD          = [0.926 0.617 0.112; 0.284 0.000 0.201];
% % % d.realDat{cD} = (tmpD(1,:) - tmpD(2,:)) .* 0.138657793;
% % % tmpErr        = abs([0.494 0.153 -0.260 ; -0.068 -0.363 -0.169] - tmpD);
% % % d.realSTE{cD} = (sqrt(sum(tmpErr.^2)) .* 0.138657793) ./sqrt(17);
% % % d.exp(cD)     = 10;
% % % 
% % % 
% % % d.psiSplit{cD} = { {'bodyPart',{'HandBySide'}} , ...
% % %     {'rew',{1}} }; 
% % % 
% % % % ==============================
% % % % SCRAMBLED OBJECT on track away from hand 24
% % % cD               = cD + 1;
% % % d(cD,:)          = d(cD - 1,:);
% % % 
% % % % Create real values
% % % tmpD          = [0.520 0.481 0.700; 0.666 0.444 0.553];
% % % d.realDat{cD} = (tmpD(1,:) - tmpD(2,:)) .* 0.138657793;
% % % tmpErr        = abs([0.119 0.094 0.292 ; 0.357 0.109 0.251] - tmpD);
% % % d.realSTE{cD} = (sqrt(sum(tmpErr.^2)) .* 0.138657793) ./sqrt(17);
% % % d.exp(cD)     = 10;
% % % 
% % % 
% % % d.psiSplit{cD} = { {'bodyPart',{'None'}} , ...
% % %     {'rew',{1}} }; 


% ==============================
% PROTOTYPICAL OBJECT judgements - SCRAMBLED OBJECT judgements 23
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'HandSide';
d.tact{cD}       = 'HandReachEEG';
d.cmPos{cD}      = [20 87.5 160 ; ...
                     0,  0,   0 ; ...
                     0,  0,   0 ];

d = AddbinPos(d,sHndSide,cD); % approaching hand on chest

% Create real values
tmpD          = [0.926 0.617 0.112; 0.284 0.000 0.201 ; ...
                 0.520 0.481 0.700; 0.666 0.444 0.553];
d.realDat{cD} = (tmpD(1,:) - tmpD(2,:) - ( tmpD(3,:) - tmpD(4,:) )) .* 0.138657793;
tmpErr        = abs([0.494 0.153 -0.260 ; -0.068 -0.363 -0.169; ...
                     0.119 0.094 0.292 ; 0.357 0.109 0.251] - tmpD);
d.realSTE{cD} = (sqrt(sum(tmpErr.^2)) .* 0.138657793) ./sqrt(17);
d.exp(cD)     = 10;


d.psiSplit{cD} = { {'bodyPart',{'HandBySide'}} , ...
    {'rew',{1}} }; 


% #################################################################
% Fit Ronga data: alpha power during tool use
% CREATE data table

cD = 23;

% ==============================
% AFETR COGNITIVE TRAINING: alpha div beta, to baseline because beta is unchanged
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'HandChest';
d.tact{cD}       = 'HandToolEEG';
d.cmPos{cD}      = [5, 145 ; ...
                    0,  0 ; ...
                    0,  0 ];

d = AddbinPos(d,sHndToolRake,cD); % approaching hand on chest

% % % % Create real values: subtractive normalisation
% % % tmpD          = [1.748 1.491; 0.557 0.555];
% % % d.realDat{cD} = (tmpD(1,:) - tmpD(2,:) ) .* 0.343524562;
% % % tmpErr        = abs([1.946 1.689 ; 0.654 0.640 ] - tmpD);
% % % d.realSTE{cD} = (sqrt(sum(tmpErr.^2)) .* 0.343524562) ;

% Create real values: divisive normalisation: Account for an OVERALL effect
% of training type on ALL power
tmpD          = [1.748 1.491] ;
d.realDat{cD} = tmpD(1,:) .* 0.343524562  ./ 0.556;
tmpErr        = abs([1.946 1.689] - tmpD)  ;
d.realSTE{cD} = tmpErr .* 0.343524562 ./ 0.556 ;%./sqrt(18);

d.exp(cD)     = 11;

d.psiSplit{cD} = { {'bodyPart',{'HandByChest'}} , ...
    {'rew',{1,-1}} }; 


% ==============================
% AFETR TOOL TRAINING: alpha div beta, to baseline because beta is unchanged
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);

% Create real values: divisive normalisation: Account for an OVERALL effect
% of training type on ALL power
tmpD          = [2.407 2.403] ;
d.realDat{cD} = tmpD(1,:) .* 0.343524562 ./ 0.661;
tmpErr        = abs([2.933 2.981] - tmpD) ;
d.realSTE{cD} = tmpErr .* 0.343524562 ./ 0.661 ;%./sqrt(18);

d.exp(cD)     = 11;


d.psiSplit{cD} = { {'bodyPart',{'HandByChestPlusTool'}} , ...
    {'rew',{1,-1}} }; 


% #################################################################
% Fit Quinlan dPOS data: fMRI moving vs stationary
% CREATE data table

cD = 25;

% ==============================
% AFETR COGNITIVE TRAINING: alpha div beta, to baseline because beta is unchanged
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'HeadMRIdPOS';
d.cmPos{cD}      = [15, 38, 84 ; ...
                     0,  0,  0 ; ...
                     0,  0,  0 ];

d = AddbinPos(d,sHed,cD); % approaching hand on chest


% Create real values: divisive normalisation: Account for an OVERALL effect
% of training type on ALL power
tmpD          = [2.237 1.661 1.036] ;
d.realDat{cD} = tmpD(1,:) .* 0.375335121 ;
tmpErr        = abs([3.557 2.675 1.696] - tmpD)  ;


% Transform from 95% confidence interval to SE
d.realSTE{cD} = ((abs([2.104 0.972 0.677] .* 87.95074758 - d.realDat{cD})) ./ 47.5) .* 34.1 ;
d.realSTE{cD} = tmpErr .* 0.375335121 .* (68.2/99.17) ./sqrt(18);

d.exp(cD)     = 12;

d.psiSplit{cD} = { {'bodyPart',{'Head'}} , ...
    {'rew',{1,-1}} }; 


% #################################################################
% Fit Holt DIPS PMV data: fMRI approach - withdrawal difference, faces cars spheres
% CREATE data table

cD = 26;

% ==============================
% LOOMING - RECEDING FACE, CAR, SPHERE DIPS
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'HeadMRIdDIPS';
d.cmPos{cD}      = [20, 150, 150; ...
                     0,   0,   0; ...
                     0,   0,   0];

d = AddbinPos(d,sHed,cD); % approaching hand on chest

% Create real values
tmpD          = [1.866 1.616 1.262; 0.207 1.367 1.315];
d.realDat{cD} = (tmpD(1,:) - tmpD(2,:)) .* 0.049979175;
tmpErr        = abs([2.231 1.911  1.543 ; 0.421  1.456 1.421] - tmpD);
d.realSTE{cD} = (sqrt(sum(tmpErr.^2)) .* 0.049979175) ;
d.exp(cD)     = 10;

d.exp(cD)     = 13;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'Head'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 



% ==============================
% LOOMING - RECEDING FACE, CAR, SPHERE PMv
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'HeadMRIdDIPS';
d.cmPos{cD}      = [20, 150, 150; ...
                     0,   0,   0; ...
                     0,   0,   0];

d = AddbinPos(d,sHed,cD); % approaching hand on chest

% Create real values
tmpD          = [1.174 0.558 0.661; -0.056 0.548 0.652];
d.realDat{cD} = (tmpD(1,:) - tmpD(2,:)) .* 0.049979175;
tmpErr        = abs([1.546 0.914 1.016; 0.126 0.796 0.935] - tmpD);
d.realSTE{cD} = (sqrt(sum(tmpErr.^2)) .* 0.049979175) ;
d.exp(cD)     = 10;

d.exp(cD)     = 13;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'Head'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 


% #################################################################
% Fit Graziano single neuron responses: ARM
% CREATE data table

cD = 28;

% ==============================
% ARM STRAIGHT, ARM RESPONSES
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Arm';
d.tact{cD}       = 'Macaque';
d.cmPos{cD}      = [ 20,  20,  20,  20; ...
                    -30,   0,  30,  60; ...
                     10,  10,  10,  10];

d = AddbinPos(d,sArm,cD); % approaching arm

% Create real values
tmpD          = [1.077 1.525 2.593 1.683];
d.realDat{cD} = (tmpD(1,:)) .* 36.954915;
tmpErr        = abs([1.244 1.709 2.652 1.900] - tmpD);
d.realSTE{cD} = tmpErr .* 36.954915 ;

d.exp(cD)     = 14;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'ArmForward'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 


% ==============================
% ARM LEFT, ARM RESPONSES
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
% Create real values
tmpD          = [1.575 2.069 1.682 1.219];
d.realDat{cD} = (tmpD(1,:)) .* 36.954915;
tmpErr        = abs([1.773 2.215 1.825 1.458] - tmpD);
d.realSTE{cD} = tmpErr .* 36.954915 ;

d.exp(cD)     = 14;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'ArmLeft'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 


% #################################################################
% Fit Graziano single neuron responses: HEAD
% CREATE data table

cD = 30;

% ==============================
% HEAD STRAIGHT, HEAD RESPONSES
cD               = cD + 1;
d.dir(cD)        = 1;
d.tactloc{cD}    = 'Head';
d.tact{cD}       = 'MacaqueHead';
% % % d.cmPos{cD}      = { 5.*[1:17], 5.*[1:19], 5.*[1:20], 5.*[1:19], 5.*[1:17]; ...
% % %                      5.*-[2 2 3 4 4 5 6 6 7 8 8 9 10 10 11 12 12], 5.*-[1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 6] , 5.*zeros([20 1]) , 5.*[1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 6] , 5.*[2 2 3 4 4 5 6 6 7 8 8 9 10 10 11 12 12]; ... 
% % %                      zeros([17 1]) ,  zeros([19 1]) ,  zeros([20 1]) ,  zeros([19 1]) , zeros([17 1]) };
% % % % Shorten slightly so it all fits in calculated area
% % % d.cmPos{cD}      = { 5.*[3:13], 5.*[3:15], 5.*[3:16], 5.*[3:15], 5.*[3:13]; ...
% % %                      5.*[3 4 5 5 6 7 7 8 9 9 10 ], 5.*[ 2 2 3 3 3 4 4 4 5 5 5 6 6] , 5.*zeros([1 14]) , 5.*-[ 2 2 3 3 3 4 4 4 5 5 5 6 6 ] , 5.*-[3 4 5 5 6 7 7 8 9 9 10 ]; ... 
% % %                      zeros([1 11]) ,  zeros([1 13]) ,  zeros([1 14]) ,  zeros([1 13]) , zeros([1 11]) };

% % 50cm 
% d.cmPos{cD}      = { 5.*[3:12], 5.*[3:12], 5.*[3:12], 5.*[3:12], 5.*[3:12]; ...
%                      5.*[3 4 5 5 6 7 7 8 9 9], 5.*[ 2 2 3 3 3 4 4 4 5 5 ] , 5.*zeros([1 10]) , 5.*-[ 2 2 3 3 3 4 4 4 5 5 ] , 5.*-[3 4 5 5 6 7 7 8 9 9 ]; ... 
%                      zeros([1 10]) ,  zeros([1 10]) ,  zeros([1 10]) ,  zeros([1 10]) , zeros([1 10]) };

% Start even further from face
d.cmPos{cD}      = { 5.*[4:13], 5.*[4:15], 5.*[4:16], 5.*[4:15], 5.*[4:13]; ...
                     5.*[4 5 5 6 7 7 8 9 9 10 ], 5.*[ 2 3 3 3 4 4 4 5 5 5 6 6] , 5.*zeros([1 13]) , 5.*-[ 2 3 3 3 4 4 4 5 5 5 6 6 ] , 5.*-[ 4 5 5 6 7 7 8 9 9 10 ]; ... 
                     zeros([1 10]) ,  zeros([1 12]) ,  zeros([1 13]) ,  zeros([1 12]) , zeros([1 10]) };



d = AddbinPos(d,sHed,cD); % approaching arm

% Create real values
tmpD          = [2.187 2.315 4.551 2.742 1.919];
d.realDat{cD} = (tmpD(1,:)) .* 20.55921053
tmpErr        = abs([2.763 2.635 4.756 2.971 2.184] - tmpD);
d.realSTE{cD} = tmpErr .* 20.55921053 ;

d.exp(cD)     = 15;
% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'HeadConstr'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 



% ==============================
% HEAD ROTATED, HEAD RESPONSES
cD               = cD + 1;
d(cD,:)          = d(cD - 1,:);
% Create real values
tmpD          = [1.528 1.910 2.816 3.868 2.358];
d.realDat{cD} = (tmpD(1,:)) * 20.55921053;
tmpErr        = abs([1.757 2.146 3.070 4.153 2.580] - tmpD);
d.realSTE{cD} = tmpErr * 20.55921053;

% Create psisplit
d.psiSplit{cD} = { {'bodyPart',{'RotatedHead'}} , ...
    {'rew',{1,-1}}, {'dir',{1}} }; 


% Use the same settings for calculating all psi [just with different
% features]
for iD = 2:size(d,1)
    d.psiSettings(iD,:) = d.psiSettings(1,:);
end

% Set the psis which need to average over rows
% % % d.psiSettings{cD-5,3} = 'avOverRows'; % fMRI approach - withdrawal difference, faces cars spheres
% % % d.psiSettings{cD-4,3} = 'avOverRows';
d.psiSettings{cD-3,3} = 'avOverRows'; % Macaque single neurons
d.psiSettings{cD-2,3} = 'avOverRows';
d.psiSettings{cD-1,3} = 'avOverRows'; 
d.psiSettings{cD,3} = 'avOverRows';


% add d-counter to d
for iD = 1:size(d,1)
    d.iD(iD) = iD;
end


dBackup = d;

%% Data normalisation options

% % Shift to zero mean for each default defined experiment
% allExps = unique(d.exp);
% for iExp = 1:numel(allExps)
%     inclD   = find(d.exp == iExp);
%     tmpD    = [d.realDat{inclD}];
%     tmpMean = mean(tmpD);
%     
%     % loop through all d that belong to experiment and remove mean
%     for cD = inclD'
%         d.realDat{cD} = d.realDat{cD} - tmpMean
%     end
% end


% % z score for each default defined experiment
% allExps = unique(d.exp);
% for iExp = 1:numel(allExps)
%     inclD   = find(d.exp == iExp);
%     tmpD    = [d.realDat{inclD}];
%     tmpMean = mean(tmpD);
%     tmpSD = std(tmpD);
%     
%     % loop through all d that belong to experiment and remove mean
%     for cD = inclD'
%         d.realDat{cD} = (d.realDat{cD} - tmpMean) ./ tmpSD;
%         d.realSTE{cD} = (d.realSTE{cD} ) ./ tmpSD;
%     end
% end




% scale between 0 and 1 for each default defined experiment
allExps = unique(d.exp);
for iExp = 1:numel(allExps)
    inclD   = find(d.exp == iExp);
    tmpD    = [d.realDat{inclD}];
    tmpMin  = min(tmpD);
    tmpMax  = max(tmpD);
    
    % loop through all d that belong to experiment and remove mean
    for cD = inclD'
        d.realDat{cD} = (d.realDat{cD} - tmpMin) ./ (tmpMax - tmpMin);
        d.realSTE{cD} = (d.realSTE{cD} ) ./ (tmpMax - tmpMin);
    end
end



% % Shift to zero minimum for each default defined experiment
% allExps = unique(d.exp);
% for iExp = 1:numel(allExps)
%     inclD   = find(d.exp == iExp);
%     tmpD    = [d.realDat{inclD}];
%     tmpMin  = min(tmpD);
%     
%     % loop through all d that belong to experiment and remove mean
%     for cD = inclD'
%         d.realDat{cD} = d.realDat{cD} - tmpMin;
%     end
% end




%% Other offset options

% Set single offset per paper
d.exp( 1:end) = 1;

% % Set single offset per paper
% d.exp( 1:16) = 1;
% d.exp(17:18) = 2;
% d.exp(19:22) = 3;
% d.exp(23) = 4;
% d.exp(24:25) = 5;
% d.exp(24:25) = 6;
% d.exp(26) = 7;
% d.exp(27:28) = 8;
% d.exp(29:32) = 8;


% % Set single offset per modality and paper
% d.exp( 1:14) = 1;
% d.exp(15:16) = 2;
% d.exp(17:18) = 3;
% d.exp(19:22) = 4;

% Set offset for each response modality
for iD = [1:16 19:22]
    d.exp(iD) = 1;
end
for iD = 17:18
    d.exp(iD) = 2;
end
for iD = 23 % EEG percent task a vs task b
    d.exp(iD) = 3;
end
for iD = 24:25 % EEG overall task c vs task d
    d.exp(iD) = 4;
end
for iD = 26 % fMRI PERCENT move vs stat
    d.exp(iD) = 5;
end
for iD = 27:28 % fMRI PERCENT Approach vs withdraw
    d.exp(iD) = 6;
end
for iD = 29:32  % Macaques single neur
    d.exp(iD) = 7;
end


% % Set single offset per modality and tactile location
% d.exp([1:2 5 7 13 14])  = 1; % Trunk audiotact NOEL
% d.exp([3 4 6 8:10])     = 2; % Hand audiotact NOEL
% d.exp([11:12])          = 3; % Head audiotact NOEL
% d.exp([15])             = 4; % Trunk visuotact NOEL
% d.exp([16])             = 5; % Head visuotact NOEL
% d.exp([17:18])          = 6; % Hand CCE HOLMES
% d.exp([19:22])          = 7; % Hand visuotact HUIJSMANS


% % THIS WORKS WHEN SETTING EXPERIMENT MEAN TO ZERO
% % Set single offset per tactile location and response type
% d.exp([1:2 5 7 13 14 15])  = 1; % Trunk NOEL
% d.exp([3 4 6 8:10])        = 2; % Hand NOEL
% d.exp([11:12 16])          = 3; % Head  NOEL
% d.exp([17:18])             = 4; % Hand CCE HOLMES
% d.exp([19:22])             = 5; % Hand visuotact HUIJSMANS

% % Maximal intercepts
% for iD = 1:18
%     d.exp(iD) = ceil(iD./2);
% end
% d.exp([19:22])  = ceil((iD+1)./2); 


% % % % Set offset for every condition
% % % for iD = 1:size(d,1)
% % %     d.exp(iD) = iD;
% % % end

%% Other multiplier options

% Set mutliplier to be the same everywhere 
for iD = 1:size(d,1)
    d.tact{iD} = 'Default';
end

% % Set mutliplier to be the same for each paper
% for iD = 1:16
%     d.tact{iD} = ['NOEL'];
% end
% for iD = 17:18
%     d.tact{iD} = ['HOLMES'];
% end
% for iD = 19:22
%     d.tact{iD} = ['HUIJSMANS'];
% end
% for iD = 23
%     d.tact{iD} = ['WAIMAN'];
% end
% for iD = 24:25
%     d.tact{iD} = ['RONGA'];
% end
% for iD = 26
%     d.tact{iD} = ['QUINLAN'];
% end
% for iD = 27:28
%     d.tact{iD} = ['HOLT'];
% end
% for iD = 29:32
%     d.tact{iD} = ['GRAZIANO'];
% end



% Set mutliplier to be the same for each response modality
for iD = [1:16 19:22]
    d.tact{iD} = ['VisTat'];
end
for iD = 17:18
    d.tact{iD} = ['CCE'];
end
for iD = 23
    d.tact{iD} = ['EEG perc task a vs task b'];
end
for iD = 24:25
    d.tact{iD} = ['EEG overall task c vs task d'];
end
for iD = 26
    d.tact{iD} = ['fMRI percent move vs stat'];
end
for iD = 27:28
    d.tact{iD} = ['fMRI perent appraoch vs withdraw'];
end
for iD = 29:32
    d.tact{iD} = ['NEUR'];
end



% % Set mutliplier to be the same for each tactile location by response
% % modality
% for iD = 19:22
%     d.tact{iD} = ['Hand'];
% end

%% Other stimulus posision options
% 
% % Move everything one bin away so that near pos is 5-10 cm instead of 0-5
% % Set offset for every condition
% for iD = 1:size(d,1)
%     if ~strcmp(d.tact{iD},'HandSpider')
%         d.binPos{iD}(1,:) = d.binPos{iD}(1,:) - 1;
%     end
% end



% Move everything one bin away so that near pos is 5-10 cm instead of 0-5
% Set offset for every condition
for iD = 1:size(d,1)
    if ~iscell(d.cmPos{iD})
        if iD <= 28 % don't shift monkey data
            
            d.binPos{iD}(1,:) = d.binPos{iD}(1,:) - 1;
            
            % don't shift position 9 or earlier, because things will go out
            % of bounds, and responses are basically the that far away
            % anyway
            d.binPos{iD}(1,d.binPos{iD}(1,:) < 9) = 9;

        end
    end
end


%% Other options for psi: 



% % % % ======================= THIS ONE IS GOOD ================================
% % % % % Using rule of Trunk-primacy 
% % % d.psiSplit{1}  = { {'bodyPart',{'Trunk','Head','HandByChest'}},  {'dir',{1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{2}  = { {'bodyPart',{'Trunk','Head','HandByChest'}},  {'dir',{-1}} , {'rew',{1,-1}}};
% % % d.psiSplit{3}  = { {'bodyPart',{'HandBySide'}},         {'dir',{1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{4}  = { {'bodyPart',{'HandBySide'}},         {'dir',{-1}} , {'rew',{1,-1}}};
% % % d.psiSplit{5}  = { {'bodyPart',{'Trunk','Head','HandByChest'}}, {'dir',{1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{6}  = { {'bodyPart',{'Trunk','HandByChest'}},        {'dir',{1}}  , {'rew',{1,-1}}};
% % % d.psiSplit{7}  = { {'bodyPart',{'Trunk','Head','HandByChest'}}, {'dir',{-1}} , {'rew',{1,-1}}}; 
% % % d.psiSplit{8}  = { {'bodyPart',{'Trunk','HandByChest'}},        {'dir',{-1}} , {'rew',{1,-1}}};
% % % d.psiSplit{9}  = { {'bodyPart',{'HandBySide'}},         {'dir',{1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{10} = { {'bodyPart',{'Trunk','HandByChest'}},        {'dir',{1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{11}  = { {'bodyPart',{'Head'}},                      {'dir',{1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{12} = { {'bodyPart',{'Head'}},                       {'dir',{1}}  , {'rew',{1,-1}}};
% % % d.psiSplit{13} = { {'bodyPart',{'Head','Trunk'}},               {'dir',{1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{14} = { {'bodyPart',{'Head','Trunk'}},               {'dir',{1}}  , {'rew',{1,-1}}};
% % % d.psiSplit{15} = { {'bodyPart',{'Head','Trunk'}},               {'dir',{1}}  , {'rew',{1,-1}}};
% % % d.psiSplit{16} = { {'bodyPart',{'Head'}},                       {'dir',{1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{17} = { {'bodyPart',{'HandByChest'}} ,       {'dir',{1}}  , {'rew',{1,-1}}};
% % % d.psiSplit{18} = { {'bodyPart',{'Trunk','HandByChestPlusTool'}},{'dir',{1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{19} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,  {'dir',{1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{20} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,  {'dir',{1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{21} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,  {'dir',{1}}  , {'rew',{1,-1.25}}}; 
% % % d.psiSplit{22} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,  {'dir',{1}}  , {'rew',{1,-1.25}}}; 
% % % d.psiSplit{23} = { {'bodyPart',{'HandBySide'}}               ,  {'dir',{1}}  , {'rew',{1}} }; % This only positive rew is justified by the 'reachability' task 
% % % % % % d.psiSplit{23} = { {'bodyPart',{'HandBySide','Head','Trunk'}}               ,  {'dir',{1}}  , {'rew',{1}} }; % This only positive rew is justified by the 'reachability' task 
% % % % % % d.psiSplit{24} = { {'bodyPart',{'None'}}                     ,  {'dir',{1}}  , {'rew',{0}} }; % This only positive rew is justified by the 'reachability' task 
% % % d.psiSplit{24} = { {'bodyPart',{'Trunk','HandByChest'}} ,       {'dir',{1,-1}}  , {'rew',{1,-1}}};
% % % d.psiSplit{25} = { {'bodyPart',{'Trunk','HandByChestPlusToolRake'}},{'dir',{1,-1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{26} = { {'bodyPart',{'Head'}},{'dir',{1,-1}}  , {'rew',{1,-1}}}; 
% % % d.psiSplit{27} = { {'bodyPart',{'Head'}},{'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
% % % d.psiSplit{28} = { {'bodyPart',{'Head'}},{'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
% % % d.psiSplit{29} = { {'bodyPart',{'ArmForward'}},{'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
% % % d.psiSplit{30} = { {'bodyPart',{'ArmLeft'}},{'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
% % % d.psiSplit{31} = { {'bodyPart',{'HeadConstr'}},{'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head forward
% % % d.psiSplit{32} = { {'bodyPart',{'RotatedHead'}},{'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head to side



dQ = d;
% ======================= THIS ONE IS GOOD ================================
% % Using rule of Trunk-primacy 
dQ.psiSplit{1}  = { {'bodyPart',{'Trunk','Head','HandByChest'}}, {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{2}  = { {'bodyPart',{'Trunk','Head','HandByChest'}}, {'dir',{-1}} , {'rew',{1,-1}}};
dQ.psiSplit{3}  = { {'bodyPart',{'HandBySide'}},                 {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{4}  = { {'bodyPart',{'HandBySide'}},                 {'dir',{-1}} , {'rew',{1,-1}}};
dQ.psiSplit{5}  = { {'bodyPart',{'Trunk','Head','HandByChest'}}, {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{6}  = { {'bodyPart',{'Trunk','HandByChest'}},        {'dir',{1}}  , {'rew',{1,-1}}};
dQ.psiSplit{7}  = { {'bodyPart',{'Trunk','Head','HandByChest'}}, {'dir',{-1}} , {'rew',{1,-1}}}; 
dQ.psiSplit{8}  = { {'bodyPart',{'Trunk','HandByChest'}},        {'dir',{-1}} , {'rew',{1,-1}}};
dQ.psiSplit{9}  = { {'bodyPart',{'HandBySide'}},                 {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{10} = { {'bodyPart',{'Trunk','HandByChest'}},        {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{11} = { {'bodyPart',{'Head'}},                       {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{12} = { {'bodyPart',{'Head'}},                       {'dir',{1}}  , {'rew',{1,-1}}};
dQ.psiSplit{13} = { {'bodyPart',{'Head','Trunk'}},               {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{14} = { {'bodyPart',{'Head','Trunk'}},               {'dir',{1}}  , {'rew',{1,-1}}};
dQ.psiSplit{15} = { {'bodyPart',{'Head','Trunk'}},               {'dir',{1}}  , {'rew',{1,-1}}};
dQ.psiSplit{16} = { {'bodyPart',{'Head'}},                       {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{17} = { {'bodyPart',{'HandByChest'}} ,               {'dir',{1}}  , {'rew',{1,-1}}};
dQ.psiSplit{18} = { {'bodyPart',{'Trunk','HandByChestPlusTool'}},{'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{19} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,  {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{20} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,  {'dir',{1}}  , {'rew',{1,-1}}}; 
dQ.psiSplit{21} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,  {'dir',{1}}  , {'rew',{1,-1.25}}}; 
dQ.psiSplit{22} = { {'bodyPart',{'Trunk','HandByChestTrack'}} ,  {'dir',{1}}  , {'rew',{1,-1.25}}}; 
dQ.psiSplit{23} = { {'bodyPart',{'HandBySide','Head','Trunk'}},  {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dQ.psiSplit{24} = { {'bodyPart',{'Trunk','Head','HandByChest'}}, {'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL
dQ.psiSplit{25} = { {'bodyPart',{'Trunk','Head','HandByChestPlusToolRake'}},{'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL 
dQ.psiSplit{26} = { {'bodyPart',{'Trunk','Head','HandBySide'}},  {'dir',{1,-1}}  , {'rew',{1,-1}}};  % fMRI moving vs stationary dPOS
dQ.psiSplit{27} = { {'bodyPart',{'Trunk','Head','HandBySide'}},  {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dQ.psiSplit{28} = { {'bodyPart',{'Trunk','Head','HandBySide'}},  {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dQ.psiSplit{29} = { {'bodyPart',{'ArmForward'}},                 {'dir',{ 1}}  ,  {'rew',{1,-1}}}; % Single neuron - arm forward
dQ.psiSplit{30} = { {'bodyPart',{'ArmLeft'}},                    {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
dQ.psiSplit{31} = { {'bodyPart',{'HeadConstr'}},                 {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head forward
dQ.psiSplit{32} = { {'bodyPart',{'RotatedHead'}},                {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head to side




% % NO negative valence
% for iD = 1:size(d,2)
%     d.psiSplit{iD}{3} = {'rew',{1}};
% end


dUn = d;
% ======================= Uncertain Q-values ================================
% % Using rule of Trunk-primacy 
dUn.psiSplit{1}  = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestUnc'}},  {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{2}  = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestUnc'}},  {'dir',{-1}} , {'rew',{1,-1}}};
dUn.psiSplit{3}  = { {'bodyPart',{'HandBySideUnc'}},                        {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{4}  = { {'bodyPart',{'HandBySideUnc'}},                        {'dir',{-1}} , {'rew',{1,-1}}};
dUn.psiSplit{5}  = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestUnc'}},  {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{6}  = { {'bodyPart',{'TrunkUnc','HandByChestUnc'}},            {'dir',{1}}  , {'rew',{1,-1}}};
dUn.psiSplit{7}  = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestUnc'}},  {'dir',{-1}} , {'rew',{1,-1}}}; 
dUn.psiSplit{8}  = { {'bodyPart',{'TrunkUnc','HandByChestUnc'}},            {'dir',{-1}} , {'rew',{1,-1}}};
dUn.psiSplit{9}  = { {'bodyPart',{'HandBySideUnc'}},                        {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{10} = { {'bodyPart',{'TrunkUnc','HandByChestUnc'}},            {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{11} = { {'bodyPart',{'HeadUnc'}},                              {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{12} = { {'bodyPart',{'HeadUnc'}},                              {'dir',{1}}  , {'rew',{1,-1}}};
dUn.psiSplit{13} = { {'bodyPart',{'HeadUnc','TrunkUnc'}},                   {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{14} = { {'bodyPart',{'HeadUnc','TrunkUnc'}},                   {'dir',{1}}  , {'rew',{1,-1}}};
dUn.psiSplit{15} = { {'bodyPart',{'HeadUnc','TrunkUnc'}},                   {'dir',{1}}  , {'rew',{1,-1}}};
dUn.psiSplit{16} = { {'bodyPart',{'HeadUnc'}},                              {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{17} = { {'bodyPart',{'HandByChestUnc'}} ,                      {'dir',{1}}  , {'rew',{1,-1}}};
dUn.psiSplit{18} = { {'bodyPart',{'TrunkUnc','HandByChestPlusToolUnc'}},    {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{19} = { {'bodyPart',{'TrunkUnc','HandByChestTrack'}} ,      {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{20} = { {'bodyPart',{'TrunkUnc','HandByChestTrack'}} ,      {'dir',{1}}  , {'rew',{1,-1}}}; 
dUn.psiSplit{21} = { {'bodyPart',{'TrunkUnc','HandByChestTrack'}} ,      {'dir',{1}}  , {'rew',{1,-1.25}}}; 
dUn.psiSplit{22} = { {'bodyPart',{'TrunkUnc','HandByChestTrack'}} ,      {'dir',{1}}  , {'rew',{1,-1.25}}}; 
dUn.psiSplit{23} = { {'bodyPart',{'HandBySideUnc','HeadUnc','TrunkUnc'}},   {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dUn.psiSplit{24} = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestUnc'}},  {'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL
dUn.psiSplit{25} = { {'bodyPart',{'TrunkUnc','HeadUnc','HandByChestPlusToolRakeUnc'}},{'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL 
dUn.psiSplit{26} = { {'bodyPart',{'TrunkUnc','HeadUnc','HandBySideUnc'}},   {'dir',{1,-1}}  , {'rew',{1,-1}}};  % fMRI moving vs stationary dPOS
dUn.psiSplit{27} = { {'bodyPart',{'TrunkUnc','HeadUnc','HandBySideUnc'}},   {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dUn.psiSplit{28} = { {'bodyPart',{'TrunkUnc','HeadUnc','HandBySideUnc'}},   {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dUn.psiSplit{29} = { {'bodyPart',{'ArmForwardUnc'}},                        {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
dUn.psiSplit{30} = { {'bodyPart',{'ArmLeftUnc'}},                           {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
dUn.psiSplit{31} = { {'bodyPart',{'HeadConstrUnc'}},                        {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head forward
dUn.psiSplit{32} = { {'bodyPart',{'RotatedHeadUnc'}},                       {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head to side



dSr = d;
% ======================= Uncertain Q-values ================================
% % Using rule of Trunk-primacy 
dSr.psiSplit{1}  = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestSARSA'}},  {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{2}  = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestSARSA'}},  {'dir',{-1}} , {'rew',{1,-1}}};
dSr.psiSplit{3}  = { {'bodyPart',{'HandBySideSARSA'}},                        {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{4}  = { {'bodyPart',{'HandBySideSARSA'}},                        {'dir',{-1}} , {'rew',{1,-1}}};
dSr.psiSplit{5}  = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestSARSA'}},  {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{6}  = { {'bodyPart',{'TrunkSARSA','HandByChestSARSA'}},            {'dir',{1}}  , {'rew',{1,-1}}};
dSr.psiSplit{7}  = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestSARSA'}},  {'dir',{-1}} , {'rew',{1,-1}}}; 
dSr.psiSplit{8}  = { {'bodyPart',{'TrunkSARSA','HandByChestSARSA'}},            {'dir',{-1}} , {'rew',{1,-1}}};
dSr.psiSplit{9}  = { {'bodyPart',{'HandBySideSARSA'}},                        {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{10} = { {'bodyPart',{'TrunkSARSA','HandByChestSARSA'}},            {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{11} = { {'bodyPart',{'HeadSARSA'}},                              {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{12} = { {'bodyPart',{'HeadSARSA'}},                              {'dir',{1}}  , {'rew',{1,-1}}};
dSr.psiSplit{13} = { {'bodyPart',{'HeadSARSA','TrunkSARSA'}},                   {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{14} = { {'bodyPart',{'HeadSARSA','TrunkSARSA'}},                   {'dir',{1}}  , {'rew',{1,-1}}};
dSr.psiSplit{15} = { {'bodyPart',{'HeadSARSA','TrunkSARSA'}},                   {'dir',{1}}  , {'rew',{1,-1}}};
dSr.psiSplit{16} = { {'bodyPart',{'HeadSARSA'}},                              {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{17} = { {'bodyPart',{'HandByChestSARSA'}} ,                      {'dir',{1}}  , {'rew',{1,-1}}};
dSr.psiSplit{18} = { {'bodyPart',{'TrunkSARSA','HandByChestPlusToolSARSA'}},    {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{19} = { {'bodyPart',{'TrunkSARSA','HandByChestTrack'}} ,      {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{20} = { {'bodyPart',{'TrunkSARSA','HandByChestTrack'}} ,      {'dir',{1}}  , {'rew',{1,-1}}}; 
dSr.psiSplit{21} = { {'bodyPart',{'TrunkSARSA','HandByChestTrack'}} ,      {'dir',{1}}  , {'rew',{1,-1.25}}}; 
dSr.psiSplit{22} = { {'bodyPart',{'TrunkSARSA','HandByChestTrack'}} ,      {'dir',{1}}  , {'rew',{1,-1.25}}}; 
dSr.psiSplit{23} = { {'bodyPart',{'HandBySideSARSA','HeadSARSA','TrunkSARSA'}},   {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dSr.psiSplit{24} = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestSARSA'}},  {'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL
dSr.psiSplit{25} = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandByChestPlusToolRakeSARSA'}},{'dir',{1,-1}}  , {'rew',{1,-1}}}; % EEG RONGA TOOL 
dSr.psiSplit{26} = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandBySideSARSA'}},   {'dir',{1,-1}}  , {'rew',{1,-1}}};  % fMRI moving vs stationary dPOS
dSr.psiSplit{27} = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandBySideSARSA'}},   {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dSr.psiSplit{28} = { {'bodyPart',{'TrunkSARSA','HeadSARSA','HandBySideSARSA'}},   {'dir',{ 1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dSr.psiSplit{29} = { {'bodyPart',{'ArmForwardSARSA'}},                        {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
dSr.psiSplit{30} = { {'bodyPart',{'ArmLeftSARSA'}},                           {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - arm forward
dSr.psiSplit{31} = { {'bodyPart',{'HeadConstrSARSA'}},                        {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head forward
dSr.psiSplit{32} = { {'bodyPart',{'RotatedHeadSARSA'}},                       {'dir',{ 1}}  , {'rew',{1,-1}}}; % Single neuron - head to side





dHP = d;
% ======================= HITPROB VERSION: ================================
% % Using rule of Trunk-primacy 
dHP.psiSplit{1}  = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestHP'}},   {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{2}  = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestHP'}},   {'dir',{-1}} , {'rew',{1}}};
dHP.psiSplit{3}  = { {'bodyPart',{'HandBySideHP'}},                       {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{4}  = { {'bodyPart',{'HandBySideHP'}},                       {'dir',{-1}} , {'rew',{1}}};
dHP.psiSplit{5}  = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestHP'}},   {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{6}  = { {'bodyPart',{'TrunkHP','HandByChestHP'}},            {'dir',{1}}  , {'rew',{1}}};
dHP.psiSplit{7}  = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestHP'}},   {'dir',{-1}} , {'rew',{1}}}; 
dHP.psiSplit{8}  = { {'bodyPart',{'TrunkHP','HandByChestHP'}},            {'dir',{-1}} , {'rew',{1}}};
dHP.psiSplit{9}  = { {'bodyPart',{'HandBySideHP'}},                       {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{10} = { {'bodyPart',{'TrunkHP','HandByChestHP'}},            {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{11} = { {'bodyPart',{'HeadHP'}},                             {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{12} = { {'bodyPart',{'HeadHP'}},                             {'dir',{1}}  , {'rew',{1}}};
dHP.psiSplit{13} = { {'bodyPart',{'HeadHP','TrunkHP'}},                   {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{14} = { {'bodyPart',{'HeadHP','TrunkHP'}},                   {'dir',{1}}  , {'rew',{1}}};
dHP.psiSplit{15} = { {'bodyPart',{'HeadHP','TrunkHP'}},                   {'dir',{1}}  , {'rew',{1}}};
dHP.psiSplit{16} = { {'bodyPart',{'HeadHP'}},                             {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{17} = { {'bodyPart',{'HandByChestHP'}} ,                     {'dir',{1}}  , {'rew',{1}}};
dHP.psiSplit{18} = { {'bodyPart',{'TrunkHP','HandByChestPlusToolHP'}},    {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{19} = { {'bodyPart',{'TrunkHP','HandByChestTrackHP'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{20} = { {'bodyPart',{'TrunkHP','HandByChestTrackHP'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{21} = { {'bodyPart',{'TrunkHP','HandByChestTrackHP'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{22} = { {'bodyPart',{'TrunkHP','HandByChestTrackHP'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dHP.psiSplit{23} = { {'bodyPart',{'HandBySideHP','HeadHP','TrunkHP'}} ,   {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dHP.psiSplit{24} = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestHP'}} ,  {'dir',{1,-1}}  , {'rew',{1}}}; % EEG RONGA TOOL
dHP.psiSplit{25} = { {'bodyPart',{'TrunkHP','HeadHP','HandByChestPlusToolRakeHP'}},{'dir',{1,-1}}  , {'rew',{1}}}; % EEG RONGA TOOL 
dHP.psiSplit{26} = { {'bodyPart',{'TrunkHP','HeadHP','HandBySideHP'}},    {'dir',{1,-1}}  , {'rew',{1}}};  % fMRI moving vs stationary dPOS
dHP.psiSplit{27} = { {'bodyPart',{'TrunkHP','HeadHP','HandBySideHP'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dHP.psiSplit{28} = { {'bodyPart',{'TrunkHP','HeadHP','HandBySideHP'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dHP.psiSplit{29} = { {'bodyPart',{'ArmForwardHP'}},                       {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dHP.psiSplit{30} = { {'bodyPart',{'ArmLeftHP'}},                          {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dHP.psiSplit{31} = { {'bodyPart',{'HeadConstrHP'}},                       {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head forward
dHP.psiSplit{32} = { {'bodyPart',{'RotatedHeadHP'}},                      {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head to side


dMI = d;
% ======================= MultiSENSORY INTEGRATION VERSION: ================================
% % Using rule of Trunk-primacy 
dMI.psiSplit{1}  = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestMI'}},   {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{2}  = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestMI'}},   {'dir',{1}} , {'rew',{1}}};
dMI.psiSplit{3}  = { {'bodyPart',{'HandBySideMI'}},                       {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{4}  = { {'bodyPart',{'HandBySideMI'}},                       {'dir',{1}} , {'rew',{1}}};
dMI.psiSplit{5}  = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestMI'}},   {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{6}  = { {'bodyPart',{'TrunkMI','HandByChestMI'}},            {'dir',{1}}  , {'rew',{1}}};
dMI.psiSplit{7}  = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestMI'}},   {'dir',{1}} , {'rew',{1}}}; 
dMI.psiSplit{8}  = { {'bodyPart',{'TrunkMI','HandByChestMI'}},            {'dir',{1}} , {'rew',{1}}};
dMI.psiSplit{9}  = { {'bodyPart',{'HandBySideMI'}},                       {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{10} = { {'bodyPart',{'TrunkMI','HandByChestMI'}},            {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{11} = { {'bodyPart',{'HeadMI'}},                             {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{12} = { {'bodyPart',{'HeadMI'}},                             {'dir',{1}}  , {'rew',{1}}};
dMI.psiSplit{13} = { {'bodyPart',{'HeadMI','TrunkMI'}},                   {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{14} = { {'bodyPart',{'HeadMI','TrunkMI'}},                   {'dir',{1}}  , {'rew',{1}}};
dMI.psiSplit{15} = { {'bodyPart',{'HeadMI','TrunkMI'}},                   {'dir',{1}}  , {'rew',{1}}};
dMI.psiSplit{16} = { {'bodyPart',{'HeadMI'}},                             {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{17} = { {'bodyPart',{'HandByChestMI'}} ,                     {'dir',{1}}  , {'rew',{1}}};
dMI.psiSplit{18} = { {'bodyPart',{'TrunkMI','HandByChestPlusToolMI'}},    {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{19} = { {'bodyPart',{'TrunkMI','HandByChestTrackMI'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{20} = { {'bodyPart',{'TrunkMI','HandByChestTrackMI'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{21} = { {'bodyPart',{'TrunkMI','HandByChestTrackMI'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{22} = { {'bodyPart',{'TrunkMI','HandByChestTrackMI'}} ,      {'dir',{1}}  , {'rew',{1}}}; 
dMI.psiSplit{23} = { {'bodyPart',{'HandBySideMI','HeadMI','TrunkMI'}} ,   {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dMI.psiSplit{24} = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestMI'}} ,  {'dir',{1}}  , {'rew',{1}}}; % EEG RONGA TOOL
dMI.psiSplit{25} = { {'bodyPart',{'TrunkMI','HeadMI','HandByChestPlusToolRakeMI'}},{'dir',{1}}  , {'rew',{1}}}; % EEG RONGA TOOL 
dMI.psiSplit{26} = { {'bodyPart',{'TrunkMI','HeadMI','HandBySideMI'}},    {'dir',{1}}  , {'rew',{1}}};  % fMRI moving vs stationary dPOS
dMI.psiSplit{27} = { {'bodyPart',{'TrunkMI','HeadMI','HandBySideMI'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dMI.psiSplit{28} = { {'bodyPart',{'TrunkMI','HeadMI','HandBySideMI'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dMI.psiSplit{29} = { {'bodyPart',{'ArmForwardMI'}},                       {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dMI.psiSplit{30} = { {'bodyPart',{'ArmLeftMI'}},                          {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dMI.psiSplit{31} = { {'bodyPart',{'HeadConstrMI'}},                       {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head forward
dMI.psiSplit{32} = { {'bodyPart',{'RotatedHeadMI'}},                      {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head to side


dDi = d;
% ======================= DISTANCE VERSION: ================================
% % Using rule of Trunk-primacy 
dDi.psiSplit{1}  = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestDist'}},   {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{2}  = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestDist'}},   {'dir',{1}} , {'rew',{1}}};
dDi.psiSplit{3}  = { {'bodyPart',{'HandBySideDist'}},                           {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{4}  = { {'bodyPart',{'HandBySideDist'}},                           {'dir',{1}} , {'rew',{1}}};
dDi.psiSplit{5}  = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestDist'}},   {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{6}  = { {'bodyPart',{'TrunkDist','HandByChestDist'}},              {'dir',{1}}  , {'rew',{1}}};
dDi.psiSplit{7}  = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestDist'}},   {'dir',{1}} , {'rew',{1}}}; 
dDi.psiSplit{8}  = { {'bodyPart',{'TrunkDist','HandByChestDist'}},              {'dir',{1}} , {'rew',{1}}};
dDi.psiSplit{9}  = { {'bodyPart',{'HandBySideDist'}},                           {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{10} = { {'bodyPart',{'TrunkDist','HandByChestDist'}},              {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{11} = { {'bodyPart',{'HeadDist'}},                                 {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{12} = { {'bodyPart',{'HeadDist'}},                                 {'dir',{1}}  , {'rew',{1}}};
dDi.psiSplit{13} = { {'bodyPart',{'HeadDist','TrunkDist'}},                     {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{14} = { {'bodyPart',{'HeadDist','TrunkDist'}},                     {'dir',{1}}  , {'rew',{1}}};
dDi.psiSplit{15} = { {'bodyPart',{'HeadDist','TrunkDist'}},                     {'dir',{1}}  , {'rew',{1}}};
dDi.psiSplit{16} = { {'bodyPart',{'HeadDist'}},                                 {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{17} = { {'bodyPart',{'HandByChestDist'}} ,                         {'dir',{1}}  , {'rew',{1}}};
dDi.psiSplit{18} = { {'bodyPart',{'TrunkDist','HandByChestPlusToolDist'}},      {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{19} = { {'bodyPart',{'TrunkDist','HandByChestTrackDist'}} ,        {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{20} = { {'bodyPart',{'TrunkDist','HandByChestTrackDist'}} ,        {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{21} = { {'bodyPart',{'TrunkDist','HandByChestTrackDist'}} ,        {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{22} = { {'bodyPart',{'TrunkDist','HandByChestTrackDist'}} ,        {'dir',{1}}  , {'rew',{1}}}; 
dDi.psiSplit{23} = { {'bodyPart',{'HandBySideDist','HeadDist','TrunkDist'}} ,   {'dir',{1}}  , {'rew',{1}} }; % EEG WAMAIN This only positive rew is justified by the 'reachability' task 
dDi.psiSplit{24} = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestDist'}} ,  {'dir',{1}}  , {'rew',{1}}}; % EEG RONGA TOOL
dDi.psiSplit{25} = { {'bodyPart',{'TrunkDist','HeadDist','HandByChestPlusToolRakeDist'}},{'dir',{1}}  , {'rew',{1}}}; % EEG RONGA TOOL 
dDi.psiSplit{26} = { {'bodyPart',{'TrunkDist','HeadDist','HandBySideDist'}},    {'dir',{1}}  , {'rew',{1}}};  % fMRI moving vs stationary dPOS
dDi.psiSplit{27} = { {'bodyPart',{'TrunkDist','HeadDist','HandBySideDist'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding DIPS
dDi.psiSplit{28} = { {'bodyPart',{'TrunkDist','HeadDist','HandBySideDist'}},    {'dir',{1}}  , {'rew',{1}}}; % fMRI looming - receding PMV
dDi.psiSplit{29} = { {'bodyPart',{'ArmForwardDist'}},                           {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dDi.psiSplit{30} = { {'bodyPart',{'ArmLeftDist'}},                              {'dir',{1}}  , {'rew',{1}}}; % Single neuron - arm forward
dDi.psiSplit{31} = { {'bodyPart',{'HeadConstrDist'}},                           {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head forward
dDi.psiSplit{32} = { {'bodyPart',{'RotatedHeadDist'}},                          {'dir',{1}}  , {'rew',{1}}}; % Single neuron - head to side


%% #################################################################
% Create 'successor representations'
for iD = 1:size(d,1) 
    [dQ.psi{iD},  dQ.psiLines{iD},  dQ.psiLineDescr{iD}]  = MakePsi(allQ, dQ(iD,:));
    [dHP.psi{iD}, dHP.psiLines{iD}, dHP.psiLineDescr{iD}] = MakePsi(allQ, dHP(iD,:));
    [dMI.psi{iD}, dMI.psiLines{iD}, dMI.psiLineDescr{iD}] = MakePsi(allQ, dMI(iD,:));
    [dDi.psi{iD}, dDi.psiLines{iD}, dDi.psiLineDescr{iD}] = MakePsi(allQ, dDi(iD,:));
    [dUn.psi{iD}, dUn.psiLines{iD}, dUn.psiLineDescr{iD}] = MakePsi(allQ, dUn(iD,:));
    [dSr.psi{iD}, dSr.psiLines{iD}, dSr.psiLineDescr{iD}] = MakePsi(allQ, dSr(iD,:));

    % Add Experiment information to psilinedescr
    for iDescr = 1:numel(dQ.psiLineDescr{iD})
        dQ.psiLineDescr{iD}{iDescr}   = [dQ.psiLineDescr{iD}{iDescr}   'Exp == ' num2str(dQ.exp(iD))];
    end
    for iDescr = 1:numel(dHP.psiLineDescr{iD})
        dHP.psiLineDescr{iD}{iDescr} = [dHP.psiLineDescr{iD}{iDescr} 'Exp == ' num2str(dHP.exp(iD))];
    end
    for iDescr = 1:numel(dMI.psiLineDescr{iD})
        dMI.psiLineDescr{iD}{iDescr} = [dMI.psiLineDescr{iD}{iDescr} 'Exp == ' num2str(dMI.exp(iD))];
    end
    for iDescr = 1:numel(dDi.psiLineDescr{iD})
        dDi.psiLineDescr{iD}{iDescr} = [dDi.psiLineDescr{iD}{iDescr} 'Exp == ' num2str(dDi.exp(iD))];
    end
    for iDescr = 1:numel(dUn.psiLineDescr{iD})
        dUn.psiLineDescr{iD}{iDescr} = [dUn.psiLineDescr{iD}{iDescr} 'Exp == ' num2str(dUn.exp(iD))];
    end
    for iDescr = 1:numel(dSr.psiLineDescr{iD})
        dSr.psiLineDescr{iD}{iDescr} = [dSr.psiLineDescr{iD}{iDescr} 'Exp == ' num2str(dSr.exp(iD))];
    end
end

% Make the distances negative, so that values are maximal near the limbs,
% and positive slopes can fit the data properly
dNegDi = dDi;
dNegDi.psi = cellfun(@(psi) -psi, dDi.psi, 'UniformOutput', false);

%% ========================================================================
% Fit all data: 

% Fit the main model:
[fitRes.Q]      = Fit3Ddat(dQ,'Linear');
% Fit the uncertain model:
[fitRes.QUN]    = Fit3Ddat(dUn,'Linear');
% Fit the SARSA model:
[fitRes.SRS]    = Fit3Ddat(dSr,'Linear');
% Fit the hit-probability model:
[fitRes.HP]     = Fit3Ddat(dHP,'Linear');
% Fit the mutlisensory integration model:
[fitRes.MI]     = Fit3Ddat(dMI,'Linear');
% Fit the exponential model:
[fitRes.EXP]    = Fit3Ddat(dDi,'Exponential');
% Fit the sigmoid model:
[fitRes.SIG]    = Fit3Ddat(dDi,'Sigmoid');
% Fit the linear model:
[fitRes.LIN]    = Fit3Ddat(dNegDi,'Linear');


% Save fitted data
tic
save('Data\LatestFittedData.mat','-v7.3')
toc

% $$$ HERE HERE
%% $$$ CHECK what's going wrong:


p0 = ([150 500 .5  120 0.6 120 0.5 6 200 800 150 9 0.5  13  12.5  22 12  9.5 -37 17.9 0.1 1 -0.2 -.1 12 -70]') ; %ones([ 1, numel(unique(d.tact)).*size(d.psi{1},1) + numel(unique(d.exp))])';
lb = zeros(size(p0)); lb(1+end-numel(unique(d.exp)):end) = -100; % allow the offsets to be -ve
ub = ones(size(p0)) .* 1000;


% % % Define Function to be optimised
% % fitType = 'Linear';

% Extra multiplier for exponential effect --> make it inverse exponential
fitType = 'Exponential';
p0 = [p0;   -.16 ;    -10];
lb = [lb;    -20 ; -1000];
ub = [ub;      -.0001 ;  -0.00001];


% % % % Extra multiplier for exponential effect --> make it inverse exponential
% % % fitType = 'Sigmoid';
% % % p0 = [p0;     -8 ;     -8];
% % % lb = [lb; -1000 ; -1000];
% % % ub = [ub;     0 ;  1000];


% OK AWAY OFFSET IS JUSTIFIED BY AN EXPECTATION/ATTENTION EFFECT: all aways seem to
% have a higher baseline
% Extra multiplier for direction
p0 = [p0; 20];
lb = [lb; -100];
ub = [ub; 1000];

% Create successor matrix with separate features that are  depending on 
% touch location, for faster optimisation
[psiMat allDat psiLabels] = MakePsiMat(d);


FunToOpt = @(p) ErrorFun(psiMat,allDat,d,p,fitType);

% Set optimisation options - keep it simple
A = []; b = []; Aeq = []; beq = []; a = tic;

% Run optimisation
OPTIONS = optimset('TolCon',1e-10);
[p,sSqErr,exitflag,output,lambda,grad,hessian] = fmincon(FunToOpt,p0,A,b,Aeq,beq,lb,ub,[],OPTIONS);
optTime = toc(a)

% Extract optimised data
[sSqErrFinal, dFitFinal] = ErrorFun(psiMat,allDat,d,p,fitType);
d = InsertModelledDat(d,dFitFinal,allDat);


%% Plot the qualities of the fit against eachother somehow?

f.STATS.f = figure('Position',[20 20 1200 400]);


allF = fields(fitRes);
fieldNames = {'Q Value' , 'HitProb', 'MultInt', 'Uncert Q', 'SARSA', 'Exponential', 'Sigmoid', 'Linear'};

% % % for iF = 1:length(allF)
% % % 
% % %     sprintf(' %s model. GoF: %.2f. AIC: %.2f. BIC: %.2f', allF{iF}, fitRes.(allF{iF}).gofScore, fitRes.(allF{iF}).AIC, fitRes.(allF{iF}).BIC)
% % % %     fitRes.(allF{iF}).gofScore
% % % 
% % % end


for iF = 1:length(allF)

%     cAx = subplot(4, length(allF) , [iF,   iF + length(allF)] );
    cAx = subplot(3, length(allF) , iF );


    plot(chi2pdf(0:500 ,k),'LineWidth',2); hold on
    yLims = ylim;
    plot([fitRes.(allF{iF}).chiSq fitRes.(allF{iF}).chiSq], yLims,'r','LineWidth',2);

    xlim([0 500])

    xlabel('total error (ChiSquare)')
    if iF == 1
        ylabel('probability');
    end
%     title(['Rejection of model (null) hypothesis: p = ' num2str(fitRes.(allF{iF}).pVal1) ')'])
    title(['p = ' num2str(fitRes.(allF{iF}).pVal1) ])


    %     sprintf(' %s model. GoF: %.2f. AIC: %.2f. BIC: %.2f', allF{iF}, fitRes.(allF{iF}).gofScore, fitRes.(allF{iF}).AIC, fitRes.(allF{iF}).BIC)
    %     fitRes.(allF{iF}).gofScore

end

legend('Expected distribution of total error if data is generated by a process like the model', ...
    'Observed total error', 'Location','South')


tmpGoF = arrayfun(@(iF) fitRes.(allF{iF}).gofScore ,1:length(allF));
tmpAIC = arrayfun(@(iF) fitRes.(allF{iF}).AIC ,1:length(allF));
tmpBIC = arrayfun(@(iF) fitRes.(allF{iF}).BIC ,1:length(allF));

% % % % % cAx = subplot(4,7, 15:21);
% % % % cAx = subplot(3,7, (15:21) - 7);
% % % % bar(tmpGoF)
% % % % xlim([0.5 7.5])
% % % % ylabel('Goodness of Fit')
% % % % set(cAx,'xTicklabels',[]);
% % % % cAx.XAxis.Visible = 'off';
% % % % box off
% % % % 
% % % % % cAx = subplot(4,7, 22:28);
% % % % cAx = subplot(3,7, (22:28) - 7);
% % % % bar(( [tmpAIC; tmpBIC] - [tmpAIC(1) ; tmpBIC(1)]  )' )
% % % % xlabel('Model type')
% % % % ylabel('Delta Information Criterion')
% % % % legend('AIC','BIC')
% % % % set(cAx,'xTicklabels',fieldNames );
% % % % box off


% cAx = subplot(4,7, 15:21);
cAx = subplot(3,7, (8:21) );
yyaxis left
bar((1:length(tmpGoF)) - 0.25 , tmpGoF , .2, 'FaceColor', [0.1 0.1 0.8])
xlim([0.6 8.4])
ylabel('Goodness of Fit')
box off

% cAx = subplot(4,7, 22:28);
% cAx = subplot(3,7, (22:28) - 7);
yyaxis right
% $$$ HERE HERE
bar( (1:length(tmpGoF))       ,  ( tmpAIC - tmpAIC(1)  )' , .2 , 'FaceColor', [0.7 0 0.2] ); hold on
bar( (1:length(tmpGoF)) + 0.25,  ( tmpBIC - tmpBIC(1)  )' , .2 , 'FaceColor', [0.85 0.4 0])
% bar( (1:length(tmpGoF)) ,  ( [tmpAIC; tmpBIC] - [tmpAIC(1) ; tmpBIC(1)]  )' , 1/3 )
xlabel('Model type')
ylabel('Delta Information Criterion')
legend('GoF','AIC','BIC')
ylim([-28 140]);
set(cAx,'xTicklabels',fieldNames );
box off




% figure
% plotyy([1:7], ([tmpAIC; tmpBIC] - [tmpAIC(1) ; tmpBIC(1)]  )' , ...
%        [1:7],  tmpGoF , 'bar', 'bar')

% bar(tmpAIC)
% bar(tmpBIC)


%% Plot the qualities of the fit against eachother VERSION 2

f.STATS.f = figure('Position',[20 -20 400 1200]);


allF = fields(fitRes);
fieldNames = {'Q Value' , 'HitProb', 'MultInt', 'Uncert Q', 'SARSA', 'Exponential', 'Sigmoid', 'Linear'};

% % % for iF = 1:length(allF)
% % % 
% % %     sprintf(' %s model. GoF: %.2f. AIC: %.2f. BIC: %.2f', allF{iF}, fitRes.(allF{iF}).gofScore, fitRes.(allF{iF}).AIC, fitRes.(allF{iF}).BIC)
% % % %     fitRes.(allF{iF}).gofScore
% % % 
% % % end

tmpGoFAll = arrayfun(@(iF) fitRes.(allF{iF}).gofScore ,1:length(allF));
tmpAICAll = arrayfun(@(iF) fitRes.(allF{iF}).AIC ,1:length(allF));
tmpBICAll = arrayfun(@(iF) fitRes.(allF{iF}).BIC ,1:length(allF));

for iF = 1:length(allF)

%     cAx = subplot(4, length(allF) , [iF,   iF + length(allF)] );
    cAx = subplot(length(allF) , 2 , 1 + (iF-1) .* 2 );


    plot(chi2pdf(0:500 ,k),'LineWidth',2); hold on
    yLims = ylim;
    plot([fitRes.(allF{iF}).chiSq fitRes.(allF{iF}).chiSq], yLims,'r','LineWidth',2);

    xlim([0 500])

    xlabel('total error (ChiSquare)')
    if iF == length(allF)
        ylabel('probability');
    end
%     title(['Rejection of model (null) hypothesis: p = ' num2str(fitRes.(allF{iF}).pVal1) ')'])
    title(['p = ' num2str(fitRes.(allF{iF}).pVal1) ])


    %     sprintf(' %s model. GoF: %.2f. AIC: %.2f. BIC: %.2f', allF{iF}, fitRes.(allF{iF}).gofScore, fitRes.(allF{iF}).AIC, fitRes.(allF{iF}).BIC)
    %     fitRes.(allF{iF}).gofScore

    if iF == length(allF)
        legend('Expected distribution of total error if data is generated by a process like the model', ...
               'Observed total error', 'Location','South');
    end

    cAx = subplot(length(allF) , 2 , iF .* 2 );


    tmpGoF = tmpGoFAll(iF);
    tmpAIC = tmpAICAll(iF);
    tmpBIC = tmpBICAll(iF);

    yyaxis left
    bar((1:length(tmpGoF)) - 0.25 , tmpGoF , .2, 'FaceColor', [0.1 0.1 0.8])
%     xlim([0.6 8.4])
    ylim([0 24])
    ylabel('Normalized Error')
    box off

    % cAx = subplot(4,7, 22:28);
    % cAx = subplot(3,7, (22:28) - 7);
    yyaxis right
    % $$$ HERE HERE
    bar( (1:length(tmpGoF))       ,  ( tmpAIC - tmpAICAll(1)  )' , .2 , 'FaceColor', [0.7 0 0.2] ); hold on
    bar( (1:length(tmpGoF)) + 0.25,  ( tmpBIC - tmpBICAll(1)  )' , .2 , 'FaceColor', [0.85 0.4 0])
    % bar( (1:length(tmpGoF)) ,  ( [tmpAIC; tmpBIC] - [tmpAIC(1) ; tmpBIC(1)]  )' , 1/3 )
    xlabel('')
    ylabel('Delta IC')
    if iF == length(allF)
        legend('GoF','AIC','BIC')
    end
%     ylim([-28 140]);
    ylim([0 140]);
    set(cAx,'xTicklabels',fieldNames );
    box off

end





% figure
% plotyy([1:7], ([tmpAIC; tmpBIC] - [tmpAIC(1) ; tmpBIC(1)]  )' , ...
%        [1:7],  tmpGoF , 'bar', 'bar')

% bar(tmpAIC)
% bar(tmpBIC)





%% Plot fitted data


f.PREDVSREAL.f = figure('Position',[20 20 900 900]);
% % % errorbar(allDat, dFitFinal, [d.realSTE{:}],'.k','MarkerSize',10,'LineWidth',1,'CapSize',3);

% fType = 'EXP';
% fType = 'LIN';
fType = 'QUN';
% fType = 'HP';
% fType = 'Q';

dFitFinal   = fitRes.(fType).dFitFinal;
dToPlt      = fitRes.(fType).d;


plot(allDat, dFitFinal, '.b','MarkerSize',25);
hold on;
xlabel('Original Data (a.u.)');
ylabel('Fitted Data (a.u.)');
xlim([-20 85]);
ylim([-20 85]);
axis square
UnitLine;
title('Model fits the data well overall');

psiLabels = [psiLabels {'Exp1','Exp2','Exp3','Exp4','Exp5','Exp5'}];
f.PARAMSTRENGTH.f = figure('Position',[20 20 600 600]);
% plot(p, 'x')
bar(abs(p));
set(gca,'XTick',[1:numel(p)]);
set(gca,'XTickLabels',psiLabels);



% Figure with data from each condition separately
f.AllPLOTS.f = figure('Position',[20 20 1200 1200])
nPls = size(d,1); % number of plots
spW = ceil(sqrt(nPls)); % subplot width
for iSpl = 1:nPls
    subplot(spW,spW,iSpl)


    if iSpl == 31 | iSpl == 32
        plD = [0 0 0 0 0 ;  -30 -15 0 15 30 ; 0 0 0 0 0];
    else
        plD = dToPlt.cmPos{iSpl};
    end
    
    if iSpl >= 29
        plDim = 2;
    else
        plDim = 1;
    end

    errorbar(plD(plDim,:),dToPlt.realDat{iSpl},dToPlt.realSTE{iSpl},'.','MarkerSize',15,'LineWidth',1.5);
    hold on;
    plot(plD(plDim,:),dToPlt.fittedDat{iSpl},'b');
    plot(plD(plDim,:),zeros(size(dToPlt.cmPos{iSpl}(1,:))),'-.k');

% % %     if iSpl <= 16
% % %         ylim([-30 90])
% % %     end
% % % 
% % %     if iSpl >= 19 & iSpl <= 22
% % %         ylim([-30 80])
% % %     end 

    title(num2str(iSpl))

end





%% Figures showing data in more interpretable manners

f.NEATDATA.f        = figure('Position',[20 20 1400 1200])
f.ALLDATABYEXPTS.f  = figure('Position',[20 20 900 900])

plOpts.shadeFl  = 1;
plOpts.errBFl   = 1;
% plOpts.errBopts = {'LineWidth',2,'CapSize',3};
plOpts.errBopts = {'LineWidth',.5,'CapSize',1};

plotTitles = {'Art 1, Exp 1: AudTac Trunk Stim','Art 1, Exp 2: AudTac Hand Stim',... 1
              'Art 1, Exp 3: AudTac Trunk/Hand Stim, Hand on Chest', ... 2
              'Art 1, Exp 4: AudTac Between Hand and Chest, Hand Stim',... 3
              'Art 1, Exp 5: AudTact conr/incongr Head Stim', ... 4
              'Art 1, Exp 6: AudTact conr/incongr Chest Stim', ... 5
              'Art 1, Exp 7: VisTact Head/Chest Stim', ... 6
              'Art 2, Exps 1-7 (N=180): Tool Use', ... 7
              'Art 3: Spider vs Butterfly approaching',... 8
              'Art 4: EEG Mu during reachability judgements', ... 9
              'Art 5: EEG Alpha during rake tool use',... 10
              'Art 6: fMRI dPOS', ... 11
              'Art 7: fMRI dIPS PMv', ... 12
              'Art 8: Macaque single arm-senstive neurons', ... 13
              'Art 8: Macaque single face-senstive neurons'}; % 14

allConds = { {'Towards','Away'}, {'Towards','Away'}, ... 1
               {'Towards','Away'}, ... 2
               {'Hand by Side', 'Hand by Chest'}, ... 3 
               {'Towards Head', 'Towards Body'}, ... 4 
               {'Towards Body', 'Towards Head'}, ... 5 
               {'Stim Chest', 'Stim Head'}, ... 6 
               {'No Tool', 'Tool'}, ... 7 
               {'Butterfly Coll', 'Butterfly Non-Coll', 'Spider Coll', 'Spider Non-Coll'}, ... 8
               {''}, ... 9  'Reachability mu decrease'
               {'Post cogn training','Post tool training'},... 10
               {'Moving vs stationary difference'}, ... 11
               {'dIPS','PMv'}, ... 12
               {'Arm Right','arm Left'}, ... 13
               {'Look Ahead','Look Right'}}; % 14

plType  =    {'Line', 'Line',... 1sHndToolRake
              'Line', ... 2
              'Line', ... 3
              'Line', ... 4
              'Line', ... 5
              'Line', ... 6
              'Line', ... 7
              'Line', ... 8
              'Bar', ... 9
              'Bar', ... 10
              'Line', ... 11
              'Bar', ... 12
              'Line', ... 13
              'Line' }; % 14


yLabs  =      {'Reaction time speeding (ms)', 'Reaction time speeding (ms)',... 1
              'Reaction time speeding (ms)', ... 2
              'Reaction time speeding (ms)', ... 3
              'Reaction time speeding (ms)', ... 4
              'Reaction time speeding (ms)', ... 5
              'Reaction time speeding (ms)', ... 6
              'Reaction time speeding (ms)', ... 7
              'Reaction time speeding (ms)', ... 8
              'Norm mu power (a.u.)', ... 9
              'Norm alpha power (a.u.)', ... 10
              'Bold difference (%)', ... 11
              'Bold difference (%)', ... 12
              'Neuronal response (% of max)', ... 13
              'Neuronal response (% of max)'}; % 14


tickDescr =  { {'',''}, {'',''}, ... 1
               {'',''}, ... 2
               {'', ''}, ... 3 
               {'', ''}, ... 4 
               {'', ''}, ... 5 
               {'', ''}, ... 6 
               {'', ''}, ... 7 
               {'', '', '', ''}, ... 8
               {'Near','max reach','far'}, ... 9  'Reachability mu decrease'
               {'',''},... 10
               {''}, ... 11
               {'Faces','Cars','Spheres'}, ... 12
               {'Arm Right','Arm Left'}, ... 13
               {'Look Ahead','Look Right'}}; % 14

xLabs =  { {'Dist from body part (cm)'}, {'Dist from body part (cm)'}, ... 1
               {'Dist from body part (cm)'}, ... 2
               {'Dist from body part (cm)'}, ... 3 
               {'Dist from body part (cm)'}, ... 4 
               {'Dist from body part (cm)'}, ... 5 
               {'Dist from body part (cm)'}, ... 6 
               {'Dist from body part (cm)'}, ... 7 
               {'Dist from body part (cm)'}, ... 8
               {'Stimulus position'}, ... 9  'Reachability mu decrease'
               {'Dist from body part (cm)'},... 10
               {'Dist from body part (cm)'}, ... 11
               {'Stimulus Type'}, ... 12
               {'Dist from body part (cm)'}, ... 13
               {'Angle of approach (deg)'}}; % 14



% The lines from the data table to include
inclDs = {[1 2], [3 4], ... 1
          [5 7], ... 2
          [9 10],... 3
          [11 12], ... 4
          [13 14], ... 5
          [15 16], ... 6
          [17 18], ... 7
          [19:22], ... 8
          [23], ... 9
          [24 25], ... 10
          [26], ... 11
          [27 28] ... 12
          [29 30] ... 13
          [31 32]}; % 14

nPls = numel(inclDs); % number of plots
spW  = ceil(sqrt(nPls)); % subplot width

% colOrd = {0,0,1; 1,0,0 ; 0,1,0 ; 0,0,1};
colOrd = {[0,0,1]; [1,0,0]; [0,0,.5]; [.5,0,0]};


for iSpl = 1:nPls

    figure(f.NEATDATA.f)

    subplot(spW,spW,iSpl)

    inclRows = inclDs{iSpl};

    cPlType = plType{iSpl};

    plDat = [];
    ftDat = [];
    plSTE = [];
    plX   = [];

    for iDat = 1:numel(inclDs{iSpl})
        plDat(iDat,:) = d.realDat{inclRows(iDat)};
        plSTE(iDat,:) = d.realSTE{inclRows(iDat)};
        if iSpl == 14
            plX(iDat,:)   = d.cmPos{inclRows(iDat)}(2,:,:);
        elseif iSpl == 15
            plX(iDat,:) = [-30 -15 0 15 30];
        else
            plX(iDat,:)   = d.cmPos{inclRows(iDat)}(1,:,:);
        end
        ftDat(iDat,:) = d.fittedDat{inclRows(iDat)};
    end

    switch cPlType

        case 'Line'
            
            h = plot(plX',ftDat','LineWidth',1.5);
            set(h,{'color'},colOrd(1:numel(h),:));
            hold on;
            
            if plOpts.errBFl == 1
                h = errorbar(plX',plDat',plSTE','.','MarkerSize',20,plOpts.errBopts{:});
                set(h,{'color'},colOrd(1:numel(h),:));
            else
                h = plot(plX',plDat','.','MarkerSize',20,'LineWidth',.5);
                set(h,{'color'},colOrd(1:numel(h),:));
            end

            plot(plX',plDat','ok','MarkerSize',7,'LineWidth',1);
            
            if plOpts.shadeFl == 1
                for iD = 1:size(plX,1)
                    plOpts.c = colOrd{iD};
                    plOpts.PlotMean = 0;
                    plOpts.FaceAlpha = 0.2;
                    ShadedPlot(plX(iD,:),plDat(iD,:),plDat(iD,:) - plSTE(iD,:),plDat(iD,:) + plSTE(iD,:),plOpts);
                end
            end


        case 'Bar'
% % %             % For the fMRI data, split it in particular way
% % %             if strcmp(plotTitles{iSpl},'Art 7: fMRI dIPS')
% % %                 plDat(1,:) = [d.realDat{inclRows(1:)}];
% % %                 plDat(2,:) = [d.realDat{inclRows(4:6)}];
% % %                 plSTE(1,:) = [d.realSTE{inclRows(1:3)}];
% % %                 plSTE(2,:) = [d.realSTE{inclRows(4:6)}];
% % %                 ftDat(1,:) = [d.fittedDat{inclRows(1:3)}];
% % %                 ftDat(2,:) = [d.fittedDat{inclRows(4:6)}];
% % %             end
            
% % %             b = bar(plDat'); hold on
            b = bar(ftDat'); hold on
            for iDat = 1:size(plDat,1)
% % % 
% % %                 errorbar(b(iDat).XEndPoints, plDat(iDat,:), plSTE(iDat,:),'.k','MarkerSize',0.001,plOpts.errBopts{:});
% % % 
% % %                 h = plot(b(iDat).XEndPoints, ftDat(iDat,:),'.k','MarkerSize',20);
% % %                 plot(b(iDat).XEndPoints, ftDat(iDat,:),'ok','MarkerSize',7);
% % %                 set(h,{'color'},colOrd(iDat,:));


                h = errorbar(b(iDat).XEndPoints, plDat(iDat,:), plSTE(iDat,:),'.k','MarkerSize',0.001,plOpts.errBopts{:});
                set(h,{'color'},colOrd(iDat,:));

                h = plot(b(iDat).XEndPoints, plDat(iDat,:),'.k','MarkerSize',20);
                plot(b(iDat).XEndPoints, plDat(iDat,:),'ok','MarkerSize',7);
                set(h,{'color'},colOrd(iDat,:));


            end
           
            set(gca,'XTickLabels',tickDescr{iSpl});
            
    end

    title(plotTitles{iSpl})

    legend(allConds{iSpl})

    xlabel(xLabs{iSpl})

    ylabel(yLabs{iSpl})

    box off




% Figure with data from each condition separately

figure(f.ALLDATABYEXPTS.f );

h = plot(plDat(:), ftDat(:),'.','MarkerSize',25);
hold on;
xlabel('Original Data (a.u.)')
ylabel('Fitted Data (a.u.)')
xlim([-20 85])
ylim([-20 85])
axis square
UnitLine
title('Model fits the data well overall')


end

extraCOlOrd = [linspace(0,.8,nSpl); ...
               [linspace(.2,.8,floor(nSpl/2)), linspace(.8,.2,ceil(nSpl/2)) ] ; ...
               linspace(.2,.8,nSpl)]';
set(gca, 'ColorOrder','default')



%% Save figures


allFields = fields(f);
for iF = 1:length(allFields)
    cF = allFields{iF};

% % %     set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\' cF 'V4.tif'] , 'tif')
% % %     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\' cF 'V4.pdf'] , 'pdf')
    
%     print(f.(cF).f,'-vector','-dsvg',['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\' cF 'V4.epsc']) % svg
%     print(f.(cF).f,'-vector','-dsvg',['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\' cF 'V4.svg']) % svg
%     print(f.(cF).f,'-vector','-dpdf',['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\' cF 'V4.pdf']) % pdf


% % %     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\' cF 'V2.eps'] , 'epsc')
% % %     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\' cF 'V2.pdf'] , 'pdf')

%     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\' cF '.tif'] , 'tif')
end


%% Make figures showing body part surfaces

f.Model3DBase.f = figure('Position',[20 20 1400 1200]);

lims3D   = [105 365 ; -20 120 ; 0 150];
sizePlot = [0 300 0 130 0 150]; % [0.1 0.8 0.1 0.8 0.1 0.8]

s.plt.vidFl = 0;

% ------------------------------------------------------
% First skin for all body parts

% HAND
sFPl = sHndSide ;
plQ = allQ(13,:);
sFPl.clc.startRew       = 1;
sFPl.clc.plS.iAct       = 1:43;
sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
sFPl.clc.plS.plDim      = 3;
sFPl.clc.plS.offset     = 2; % semi-arbitrary offset to get around how voxelsurf.m deals with negative numbers
sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < 0.3;
sFPl.clc.plS.plSkin     = 'Y';
sFPl.clc.plS.plField    = 'N'; %

% [1 sY 1 sX 1 sZ]
sFPl.clc.plS.volSettings= {true,sizePlot,1};

[newQ, f] = PlotQMaps(sFPl,plQ,f); hold on

ylim(lims3D(1,:));
xlim(lims3D(2,:));
zlim(lims3D(3,:));

% set(gca,'XTick',[0.2 0.4 0.6 0.8])
% set(gca,'YTick',[0.2 0.4 0.6 0.8])
% set(gca,'ZTick',[0.2 0.4 0.6 0.8])
% set(gca,'XTickLabels',{''})
% set(gca,'YTickLabels',{''})
% set(gca,'ZTickLabels',{''})
caxis([1 4])


% BODY
sFPl = sBdy;
plQ = allQ(2,:);
sFPl.clc.startRew       = -1;
sFPl.clc.plS.iAct       = 1:7;
sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
sFPl.clc.plS.plDim      = 3;
sFPl.clc.plS.offset     = 3; 
sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < 0.3;
sFPl.clc.plS.plSkin     = 'Y';
sFPl.clc.plS.plField    = 'N';

sFPl.clc.plS.volSettings= {true,sizePlot,1};
[newQ, f] = PlotQMaps(sFPl,plQ,f);

caxis([1 4])

% HEAD
sFPl = sHed;
plQ = allQ(6,:);
sFPl.clc.startRew       = -1;
sFPl.clc.plS.iAct       = 1:7;
sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
sFPl.clc.plS.plDim      = 3;
sFPl.clc.plS.offset     = 3; 
sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < 0.3;
sFPl.clc.plS.plSkin     = 'Y';
sFPl.clc.plS.plField    = 'N';

sFPl.clc.plS.volSettings= {true,sizePlot,1};
[newQ, f] = PlotQMaps(sFPl,plQ,f);

caxis([1 4])

% MapsColourBar();

view([30 25])
axis off


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if s.plt.vidFl == 1
    % Video without fields

    s.plt.vidFileName = 'F:\Projects\DPPS\DefenseAgent\Results\Videos\RotateWithWIthoutFields.avi';

    s.plt.vidFR = 36;
    s.plt.addVidAngs = linspace(0,360 .* 1.5,s.plt.vidFR * 6)

    v = RotateAndFilm(gcf,s,[]);
end

title('Body Surface')



% ------------------------------------------------------
% Next FIELDS for all body parts

% Make copy of figure
cAx = gca;
f.Model3DFields.f = figure('Position',[20 20 1400 1200]);
cAx = copyobj(cAx,f.Model3DFields.f);

title('')


% BODY
sFPl = sBdy;
plQ = allQ(2,:);
sFPl.clc.startRew       = -1;
sFPl.clc.plS.iAct       = 1:7;
sFPl.clc.plS.ActFun     = @(x) x(1,:,:,:,:);
% sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
% sFPl.clc.plS.ActFun     = @(x) median(x,1);
sFPl.clc.plS.plDim      = 3;
sFPl.clc.plS.offset     = 3; 
sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < 0.3;
sFPl.clc.plS.plSkin     = 'N';
sFPl.clc.plS.plField    = 'Y';

sFPl.clc.plS.volSettings= {true,sizePlot,.25};
[newQ, f] = PlotQMaps(sFPl,plQ,f);

% MapsColourBar();
% title('Trunk PPS')

ylim(lims3D(1,:));
xlim(lims3D(2,:));
zlim(lims3D(3,:));


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if s.plt.vidFl == 1
    % Video with body field
    s.plt.addVidAngs = linspace(0,360 .* 1,s.plt.vidFR * 4)
    v = RotateAndFilm(gcf,s,v);
end



% HEAD
sFPl = sHed;
plQ = allQ(6,:);
sFPl.clc.startRew       = -1;
sFPl.clc.plS.iAct       = 1:7;
sFPl.clc.plS.ActFun     = @(x) x(1,:,:,:,:);
% sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
% sFPl.clc.plS.ActFun     = @(x) median(x,1);
sFPl.clc.plS.plDim      = 3;
sFPl.clc.plS.offset     = 3; 
sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < 0.3;
sFPl.clc.plS.plSkin     = 'N';
sFPl.clc.plS.plField    = 'Y';

sFPl.clc.plS.volSettings= {true,sizePlot,.25};
[newQ, f] = PlotQMaps(sFPl,plQ,f);

% MapsColourBar();

% title('Trunk + Head PPS')

ylim(lims3D(1,:));
xlim(lims3D(2,:));
zlim(lims3D(3,:));


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if s.plt.vidFl == 1
    s.plt.addVidAngs = linspace(0,360 .* 1,s.plt.vidFR * 4)
    v = RotateAndFilm(gcf,s,v);
end


% HAND
sFPl = sHndSide ;
plQ = allQ(13,:);
sFPl.clc.startRew       = 1;
sFPl.clc.plS.iAct       = 1:43;
sFPl.clc.plS.ActFun     = @(x) x(1,:,:,:,:);
% sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
% sFPl.clc.plS.ActFun     = @(x) median(x,1);
sFPl.clc.plS.plDim      = 3;
sFPl.clc.plS.offset     = 2; 
sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < 0.45
sFPl.clc.plS.plSkin     = 'N';
sFPl.clc.plS.plField    = 'Y'; %

sFPl.clc.plS.volSettings= {true,sizePlot,.25};
[newQ, f] = PlotQMaps(sFPl,plQ,f); hold on

caxis([1 4])
colormap(BlueWhiteRedDavide3)

% MapsColourBar();

ylim(lims3D(1,:));
xlim(lims3D(2,:));
zlim(lims3D(3,:));


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if s.plt.vidFl == 1
    s.plt.addVidAngs = linspace(0,360 .* 2.5 ,s.plt.vidFR * 10)
    v = RotateAndFilm(gcf,s,v);
    close(v)
end

title('Trunk + Head + Hand PPS')


%% Make figures showing body part movements for theoretical fig


lims3D   = [105 365 ; -20 120 ; 0 150];
sizePlot = [0 300 0 130 0 150]; % [0.1 0.8 0.1 0.8 0.1 0.8]

s.plt.vidFl = 0;

actNames    = {'Stay','Up','Down','Left','Right','Forward','Back'};
actCriteria = {@(x) sum(x==0,2)==3, ... STAY
               @(x) x(:,3) > 0, ...     UP
               @(x) x(:,3) < 0, ...     DOWN
               @(x) x(:,2) > 0, ...     LEFT
               @(x) x(:,2) < 0, ...     RIGHT
               @(x) x(:,1) < 0, ...     FORWARD
               @(x) x(:,1) > 0}; %      BACK

taskNames = {'Goal','Threat'};
taskDefs  = [1, -1];

bdyPartNames = {'Trunk','Head','Hand'};
bdyPartTags  = {'sBdy', 'sHed','sHnd'}'

qsToPlot = [ 1  2; ... Body goal threat (towards)
             5  6; ... Head goal threat (towards)
            13 14];  % Hand goal threat (towards)


% ------------------------------------------------------
% First skin for all body parts

f.Model3DBase.f = figure('Position',[20 20 1400 1200]);

% HAND
% Check whether it's got the wrong starting column
if all(sHndSide.clc.startSC < 15)
    sHndSide.clc.startSC = sHndSide.clc.startSC + 10
end
sFPl = sHndSide ;
plQ = allQ(14,:);
sFPl.clc.startRew       = 0.1;
sFPl.clc.plS.iAct       = 1:43;
sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
sFPl.clc.plS.plDim      = 3;
sFPl.clc.plS.offset     = 0; % semi-arbitrary offset to get around how voxelsurf.m deals with negative numbers
sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < -1;
sFPl.clc.plS.plSkin     = 'Y';
sFPl.clc.plS.plField    = 'N'; %

% [1 sY 1 sX 1 sZ]
sFPl.clc.plS.volSettings= {true,sizePlot,1};

[newQ, f] = PlotQMaps(sFPl,plQ,f); hold on

ylim(lims3D(1,:));
xlim(lims3D(2,:));
zlim(lims3D(3,:));

% set(gca,'XTick',[0.2 0.4 0.6 0.8])
% set(gca,'YTick',[0.2 0.4 0.6 0.8])
% set(gca,'ZTick',[0.2 0.4 0.6 0.8])
% set(gca,'XTickLabels',{''})
% set(gca,'YTickLabels',{''})
% set(gca,'ZTickLabels',{''})
caxis([1 4]);
hold on


% BODY
sFPl = sBdy;
plQ = allQ(2,:);
sFPl.clc.startRew       = 0.1;
sFPl.clc.plS.iAct       = 1:7;
sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
sFPl.clc.plS.plDim      = 3;
sFPl.clc.plS.offset     = 0; 
sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < -1;
sFPl.clc.plS.plSkin     = 'Y';
sFPl.clc.plS.plField    = 'N';

sFPl.clc.plS.volSettings= {true,sizePlot,1};
[newQ, f] = PlotQMaps(sFPl,plQ,f);

caxis([1 4])

% HEAD
sFPl = sHed;
plQ = allQ(6,:);
sFPl.clc.startRew       = 0.1;
sFPl.clc.plS.iAct       = 1:7;
sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
sFPl.clc.plS.plDim      = 3;
sFPl.clc.plS.offset     = 0; 
sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < -1;
sFPl.clc.plS.plSkin     = 'Y';
sFPl.clc.plS.plField    = 'N';

sFPl.clc.plS.volSettings= {true,sizePlot,1};
[newQ, f] = PlotQMaps(sFPl,plQ,f);

caxis([1 4])

% MapsColourBar();

view([-28 22.5])
axis off


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if s.plt.vidFl == 1
    % Video without fields

    s.plt.vidFileName = 'F:\Projects\DPPS\DefenseAgent\Results\Videos\RotateWithWIthoutFields.avi';

    s.plt.vidFR = 36;
    s.plt.addVidAngs = linspace(0,360 .* 1.5,s.plt.vidFR * 6)

    v = RotateAndFilm(gcf,s,[]);
end

title('Body Surface')


for iBP = 1:length(bdyPartNames)

% Select body part
eval(['sFPl = ' bdyPartTags{iBP} ';']);

for iTask = 1:length(taskDefs)

% Select data for body part and task
plQ = allQ(qsToPlot(iBP,iTask),:);

for iAct = 1:length(actNames)

% ------------------------------------------------------
% Next FIELDS for all body parts


% % % % Create new figure or switch to existing figure
% % % if iTask == 1 & iAct == 1
% % % 
% % %     figure(f.Model3DBase.f)
% % % 
% % %     currFigName = ['Model3DForTheory_' actNames{iAct} 'BodyPart_' bdyPartNames{iBP}];
% % %     currFigName = ['Model3DForTheory_' actNames{iAct} '_BodyPart_' bdyPartNames{iBP} '_Act_' num2str(iAct)];
% % % 
% % %     % Make copy of figure
% % %     cAx = gca;
% % % % % %     f.(currFigName).f = figure('Position',[20 20 2400 1200]);
% % %     f.(currFigName).f = figure('Position',[20 20 1400 1200]);
% % %     cAx = copyobj(cAx,f.(currFigName).f);
% % % else
% % %     figure(f.(currFigName).f)
% % % end

% Create new figure or switch to existing figure
figure(f.Model3DBase.f)

currFigName = ['Model3DForTheory_' actNames{iAct} 'BodyPart_' bdyPartNames{iBP}];
currFigName = ['Model3DForTheory_' actNames{iAct} '_BodyPart_' bdyPartNames{iBP} '_Act_' num2str(iAct)];

% Make copy of figure
cAx = gca;
% % %     f.(currFigName).f = figure('Position',[20 20 2400 1200]);
f.(currFigName).f = figure('Position',[20 20 1400 1200]);
cAx = copyobj(cAx,f.(currFigName).f);





sFPl.clc.plS.iAct       = actCriteria{iAct}(sFPl.clc.actConsequence);
% sFPl.clc.plS.ActFun     = @(x) x(1,:,:,:,:);
% sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
sFPl.clc.plS.ActFun     = @(x) max(abs(x),[],1);
% sFPl.clc.plS.ActFun     = @(x) median(x,1);
sFPl.clc.plS.plDim      = 3;
sFPl.clc.plS.offset     = 3; 
sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < 0.3;
sFPl.clc.plS.plSkin     = 'N';
sFPl.clc.plS.plField    = 'Y';

sFPl.clc.plS.volSettings= {true,sizePlot,.25};


% subplot(length(actNames),2,iTask + (iAct-1) .* length(taskDefs))
% subplot(2,length(actNames),iAct + (iTask-1) .* length(actNames))
% % % axS.xOffset = 0.05;
% % % axS.yOffset = 0.05
% % % axes('Position',[axS.xOffset + (iAct-1) .* ((1 - 2.*axS.xOffset) ./ length(actNames)) , ...
% % %                  axS.yOffset + (iTask-1) .* ((1 - 2.*axS.yOffset) ./ length(taskDefs)), ... 
% % %                  ((1 - 2.*axS.xOffset) ./ length(actNames)), ...
% % %                  ((1 - 2.*axS.yOffset) ./ length(taskDefs))]);

% % % axes('Position',[axS.xOffset + (iTask-1) .* ((1 - 2.*axS.xOffset) ./ length(taskDefs)) , ...
% % %                  axS.yOffset + (iAct-1) .* ((1 - 2.*axS.yOffset) ./ length(actNames)), ... 
% % %                  ((1 - 2.*axS.xOffset) ./ length(taskDefs)), ...
% % %                  ((1 - 2.*axS.yOffset) ./ length(actNames))]);
sFPl.clc.startRew   = taskDefs(iTask);
[newQ, f]           = PlotQMaps(sFPl,plQ,f);
    caxis([0 1])
% % % cLims = caxis;
% % % caxis(max(abs(cLims)) .* [-1 1] );


ylim(lims3D(1,:));
xlim(lims3D(2,:));
zlim(lims3D(3,:));

hold on

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if s.plt.vidFl == 1
    % Video with body field
    s.plt.addVidAngs = linspace(0,360 .* 1,s.plt.vidFR * 4)
    v = RotateAndFilm(gcf,s,v);
end

set(gca,'Visible','off')

title([ actNames{iAct} ' ' taskNames(iTask)])


% Adapt colourmap to task
if iTask == 1 
%     colormap(gca, whitetocol(100, [0 0 0.7]) );
    colormap(coltocol(100,[0.3 0.3 0.3],[0 0 0.7]))
else
%     colormap(gca, whitetocol(100, [0.7 0 0]) );
    colormap(coltocol(100,[0.3 0.3 0.3],[0.7 0 0]))
end



% % % % HEAD
% % % sFPl = sHed;
% % % plQ = allQ(6,:);
% % % sFPl.clc.plS.iAct       = actCriteria{iAct}(sHed.clc.actConsequence);
% % % % sFPl.clc.plS.ActFun     = @(x) x(1,:,:,:,:);
% % % % sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
% % % sFPl.clc.plS.ActFun     = @(x) max(abs(x),[],1);
% % % % sFPl.clc.plS.ActFun     = @(x) median(x,1);
% % % sFPl.clc.plS.plDim      = 3;
% % % sFPl.clc.plS.offset     = 3; 
% % % sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < 0.3;
% % % sFPl.clc.plS.plSkin     = 'N';
% % % sFPl.clc.plS.plField    = 'Y';
% % % 
% % % sFPl.clc.plS.volSettings= {true,sizePlot,.25};
% % % 
% % % for iTask = 1:length(taskNames)
% % %     subplot(length(actNames),2,iTask + (iAct-1) .* length(taskDefs))
% % %     sFPl.clc.startRew   = taskDefs(iTask);
% % %     [newQ, f] = PlotQMaps(sFPl,plQ,f);
% % % 
% % %     % MapsColourBar();
% % % 
% % %     % title('Trunk + Head PPS')
% % % 
% % %     ylim(lims3D(1,:));
% % %     xlim(lims3D(2,:));
% % %     zlim(lims3D(3,:));
% % % 
% % % 
% % %     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % %     if s.plt.vidFl == 1
% % %         s.plt.addVidAngs = linspace(0,360 .* 1,s.plt.vidFR * 4)
% % %         v = RotateAndFilm(gcf,s,v);
% % %     end
% % % end
% % % 
% % % % HAND
% % % 
% % % sFPl = sHndSide ;
% % % plQ = allQ(13,:);
% % % sFPl.clc.startRew       = 1;
% % % sFPl.clc.plS.iAct       = actCriteria{iAct}(sHnd.clc.actConsequence);
% % % % sFPl.clc.plS.ActFun     = @(x) x(1,:,:,:,:);
% % % % sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
% % % sFPl.clc.plS.ActFun     = @(x) max(abs(x),[],1);
% % % % sFPl.clc.plS.ActFun     = @(x) median(x,1);
% % % sFPl.clc.plS.plDim      = 3;
% % % sFPl.clc.plS.offset     = 2; 
% % % sFPl.clc.plS.excl3D     = @(x) abs(x-sFPl.clc.plS.offset) < 0.45
% % % sFPl.clc.plS.plSkin     = 'N';
% % % sFPl.clc.plS.plField    = 'Y'; %
% % % 
% % % sFPl.clc.plS.volSettings= {true,sizePlot,.25};
% % % 
% % % for iTask = 1:length(taskNames)
% % %     subplot(length(actNames),2,iTask + (iAct-1) .* length(taskDefs))
% % %     sFPl.clc.startRew   = taskDefs(iTask);
% % %     [newQ, f] = PlotQMaps(sFPl,plQ,f); hold on
% % % 
% % %     caxis([1 4])
% % %     colormap(BlueWhiteRedDavide3)
% % % 
% % %     % MapsColourBar();
% % % 
% % %     ylim(lims3D(1,:));
% % %     xlim(lims3D(2,:));
% % %     zlim(lims3D(3,:));
% % % 
% % % 
% % %     % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % %     if s.plt.vidFl == 1
% % %         s.plt.addVidAngs = linspace(0,360 .* 2.5 ,s.plt.vidFR * 10)
% % %         v = RotateAndFilm(gcf,s,v);
% % %         close(v)
% % %     end
% % % 
% % % end

end
end
sgtitle(bdyPartNames{iBP});
end


%%


sFPl = sHnd;
plQ = allQ(9,:);
sFPl.clc.plS.iAct       = 1:43;

% sFPl = sHndSide ;
% plQ = allQ(13,:);
% sFPl.clc.plS.iAct       = 1:43;


% sFPl = sBdy;
% plQ = allQ(2,:);
% sFPl.clc.plS.iAct       = 1:7;


sFPl = sHed;
plQ = allQ(6,:);
sFPl.clc.plS.iAct       = 1:7;

% plot settings
sFPl.clc.plS.ActFun     = @(x) max(x,[],1);
% sFPl.clc.plS.ActFun = @(x) mean(x,1);
sFPl.clc.plS.plDim      = 3;
% sFPl.clc.plS.excl3D     = @(x) abs(x) < 00;
sFPl.clc.plS.excl3D     = @(x) abs(x) < 0.1;
sFPl.clc.plS.plSkin     = 'Y';
sFPl.clc.plS.plField    = 'Y'; % If this is no, do not plot ANY of the field

% sFPl.clc.plS.volSettings= {true,[0.1 0.8 0.1 0.8 0.1 0.8],.3};
% sFPl.clc.plS.volSettings= {true,[0.7 0.8 0.7 0.8 0.7 0.8],.3};
% sFPl.clc.plS.volSettings= {true};

sFPl.clc.plS.volSettings= {true,[0.1 0.8 0.1 0.8 0.1 0.8],.5};

[newQ, f] = PlotQMaps(sFPl,plQ,f);



function [cb] = MapsColourBar()
    cLims   = caxis;
    cb      = colorbar
    cb.Ticks= cLims;
    cb.TickLabels = {'Pos','Neg'};
end

% -------------------------------------------------------------------------
function [v] = RotateAndFilm(h,s,v)
% Rotates and films a 3d plot, given by handle
    
    if isempty(v)
        v = VideoWriter(s.plt.vidFileName);
        v.FrameRate=s.plt.vidFR;
        v.Quality = 100;
        open(v)
    end
    
    [baseAz baseEl] = view;

    axis vis3d

    for iAng = 1:numel(s.plt.addVidAngs)

        view([baseAz + s.plt.addVidAngs(iAng), baseEl])

        frame=getframe(h);

        writeVideo(v,frame);
    end

end

% -------------------------------------------------------------------------
function [allQ, f] = PlotQMaps(s,allQ,f)
% Plots q values. Best not to feed the entirety of allQ at the same time,
% but rather specific rows

iAct = s.clc.plS.iAct;

for iQ = 1:size(allQ,1)
    cQ = allQ(iQ,:);

    plQ = s.clc.plS.ActFun(cQ.qVals{1}(iAct,:,:,:));

    % Switch between 1D, 2D and 3D plots
    switch s.clc.plS.plDim
        case 1
            plot(squeeze(-plQ(:,:,s.clc.nearPos(2),s.clc.nearPos(3)))' )

        case 2
            % 2D figure settings
            fS.gridXstart = -4.5;
            fS.gridXstep = 1;
            fS.gridYstart = 3.5;
            fS.gridYstep = 1;


            ax{1}    = axes('Position',[.05 .1 .4 .8]);
            imagesc(squeeze(-plQ(1,:,:,s.clc.nearPos(3))) );
            % colormap(whitetocol(100,[0 0 0.7]))
            colormap(redbluecmapRory)
            GridOverImage(fS,ax{1});
            caxis([-1 1])
            colorbar


            ax{2}    = axes('Position',[.55 .1 .4 .8]);
            imagesc(squeeze(-plQ(1,:,s.clc.nearPos(2),:))' );
            hold on; axis xy
            % colormap(whitetocol(100,[0 0 0.7]))
            colormap(redbluecmapRory)
            GridOverImage(fS,ax{2});
            caxis([-1 1])


        case 3
            %  3D plot
            switch s.clc.plS.plField
                case 'Y'
                    tmpQ = -plQ;
                case 'N'
                    tmpQ = zeros(size(plQ));
            end

            switch s.clc.plS.plSkin
                case 'Y'
                    % Set the the 'skin surface' to a particular value
                    % [hence the startSR:end in the 2nd dimension]
                    for iVol = 1:length(s.clc.startSZ)
                        tmpQ(:,s.clc.startSR(iVol),s.clc.startSC(iVol),s.clc.startSZ(iVol)) = -s.clc.startRew;
                    end

                    if isfield(s.clc,'toolSR')
                        for iVol = 1:length(s.clc.toolSR)
                            newQ(:,s.clc.toolSR(iVol),s.clc.toolSC(iVol),s.clc.toolSZ(iVol)) = -s.clc.startRew;
                        end
                    end
            end



            % Have to set tmpQ to all positive, becasue voxelsurf is weird
            % with negative numbers
            tmpQ = tmpQ + s.clc.plS.offset;
            
            % Make the excluded points invisible
            tmpQ( s.clc.plS.excl3D(tmpQ) ) = 0;

            voxelSurf(squeeze(tmpQ(1,:,:,:)),s.clc.plS.volSettings{:});

            % colormap(whitetocol(100,[0 0 0.7]))
            % colormap(redbluecmapRory)
            % colormap(BlueWhiteRedDavide3)

            

    end

    % --> AND add option for video
end

end
% -------------------------------------------------------------------------




%% FUNCTIONS


% -------------------------------------------------------------------------
function [psi psiLines lineDescr] = MakePsi(allQ, d)
% Return a successor representation specific to the dataset

[psiLines lineDescr]   = FindPsiLines(allQ, d.psiSplit{1});

% Add to psi in manner dependent on the settings
psi =[];
for iPsi = 1:numel(psiLines)


    % q values for all actions of a particular [iPsi] bodypart
    if ~iscell(d.cmPos{1})
        tmpPsi         = abs(ExtractQ(allQ.qVals(psiLines{iPsi} ),d.binPos{1}));
    end

    switch d.psiSettings{1,3}
        case 'avOverRows'
    

            cBinPos = d.binPos{1};
            

            % If stimulus poisiont is defined as a cell, average over each
            % entry within a give cell
            if iscell(d.cmPos{1})

                % Initialise to NaN so that different lengths can be used
                tmpTmpPsi = nan([max(cellfun(@(x) numel(x), d.binPos{1}),[],'all')  ...
                    size(allQ.qVals{psiLines{iPsi} },1) size(d.cmPos{1},2)]);

                for iTraj = 1:size(d.cmPos{1},2) % loop through trajectories


                    avRows = d.binPos{1}{1,iTraj};

                    for iRow = 1:numel(avRows)
                        cRow = avRows(iRow);

                        cBinPos = [d.binPos{1}{1,iTraj}(iRow) d.binPos{1}{2,iTraj}(iRow) d.binPos{1}{3,iTraj}(iRow)]';


                        % rows within trajectory, then actions, then positions,
                        % i.e. trajectories
                        tmpTmpPsi(iRow,:,iTraj)         = abs(ExtractQ(allQ.qVals(psiLines{iPsi} ),cBinPos));


                    end
                end

            % If no further indication is give, loop over 1m [100cm], stricly along
            % the row dimension
            else
                avRows = min(d.binPos{1}(1,:)) : -1 : (min(d.binPos{1}(1,:)) - 20);
                for iRow = 1:numel(avRows)
                    cRow = avRows(iRow);

                    cBinPos(1,:) = cRow;

                    tmpTmpPsi(iRow,:,:,:,:,:,:)         = abs(ExtractQ(allQ.qVals(psiLines{iPsi} ),cBinPos));


                end

                
            end

            % Take into account possibility of only 1 possible action
            if size(tmpTmpPsi,2) == 1
                tmpPsi(1,:) = squeeze(nanmean(tmpTmpPsi,1));
            else
                tmpPsi = squeeze(nanmean(tmpTmpPsi,1));
            end

            
    end
    

    
    % Placeholder for relative action probability to 'stay still' action
    switch d.psiSettings{1,2}
        case 'rel_still_norm'
            stayPsi = tmpPsi(1,:);
            tmpPsi = tmpPsi - stayPsi;

    end


    switch d.psiSettings{1,1}
        case {'average', 'average_all'}
            tmpPsi = mean(tmpPsi,1);
        case {'max','max_all','max_then_avg'}
            tmpPsi = max(tmpPsi,[],1);
        case 'raw'
    end
    psi = [psi ; tmpPsi];
end



switch d.psiSettings{1,1}
    case {'average_all','max_then_avg'}
        psi = mean(psi,1);
    case 'max_all'
        psi = mean(psi,1);
end

% Normalise across successor features --> simulates mutual inhibition, for
% example
switch d.psiSettings{1,2}
    case 'sum_norm'
        psi = psi ./ (sum(psi,1));
    case 'diff_from_mean_norm'
        psi = psi - mean(psi,1);
    case 'diff_from_mean_divmean_norm'
        psi = (psi - mean(psi,1))./ mean(psi,1);
    case 'no_norm'
end
% Adjust for case where everything dividing by zero
psi(isnan(psi)) = 0;

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function [allLines lineDescr] = FindPsiLines(allQ, psiSplit)
% Finds the lines in allQ of the Psi features that will be included
% Also returns a description of what those lines entail
elements = {};
% Create list of elements: cell array with N vectors to combine
for iSC = 1:numel(psiSplit) % loop through split conditions
    elements = [elements, {psiSplit{iSC}{2}}];
end

combinations        = cell(1, numel(elements)); %set up the varargout result
[combinations{:}]   = ndgrid(elements{:}); % feed each element as separate input
combinations        = cellfun(@(x) x(:), combinations,'uniformoutput',false); %Transform into vectors

allSplits = arrayfun(@(ind) psiSplit{ind}{1},1:numel(psiSplit),'UniformOutput',false);

% Loop through combinations to create logical indexes
allLines = [];
for iSS = 1:numel(combinations{1}) % sub split index

    % Make a Boolean specifying which lines of allQ correspond to this
    % combination of conditions
    inclLines = ones([size(allQ,1) 1]);
    lineDescr{iSS} = '';
    for iSC = 1:numel(allSplits)
        currSpl = allSplits{iSC}; % current split

        % Find the lines in allQ that should be included
        inclLines = inclLines & IsEqual(allQ.(currSpl),{combinations{iSC}{iSS}});
        % Only add to the description if the current split also has
        % alternative options. Otherwise just consider it an effect that
        % Psi naturally takes into account - i.e.e there ren't multiple
        % psis for that effect
        if numel(psiSplit{iSC}{2}) == 1
        else
            lineDescr{iSS} = [lineDescr{iSS} currSpl ' == ' num2str(combinations{iSC}{iSS}) '. '];
        end
    end

    allLines = [allLines; {find( inclLines )}];
%     lineDescr{iSS} = 
end
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function [eqBool] = IsEqual(d1,d2)
% Compares any input, whether a string or a number. Takes CELL input. Note
% at least one of the entries must only contain one value

% Convert EVERYTHING to strings inside cells
if ~ iscell(d1)
    d1 = num2cell(d1);
end
if ~ iscell(d2)
    d1 = num2cell(d2);
end
d1 = cellfun(@(subD) num2str(subD), d1, 'UniformOutput', false);
d2 = cellfun(@(subD) num2str(subD), d2, 'UniformOutput', false);

% Make sure they are the same shape
d1 = d1(:);
d2 = d2(:);

% Compare the strings
eqBool = strcmp(d1,d2);
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function [d] = InsertModelledDat(d,dFitFinal,allDat)
% Takes modelled data and puts it back into the data structure for easier
% plotting
    cP = 0; % current point
    for iD = 1:size(d,1)
        nP = numel(d.realDat{iD}); % number of points
        d.fittedDat{iD} = dFitFinal(cP + [1:nP]);
        cP = cP + nP;
    end
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function [psiMat allDat psiLabels] = MakePsiMat(d)
% Makes a successor state matrix that allows faster fitting of the data,
% assuming different multipliers for each tactile location

allTacts = unique(d.tact);
numTacts = numel(allTacts);
numFeats = size(d.psi{1},1);

allDat   = [d.realDat{:}];
numDats  = numel(allDat);

% % % psiMat = zeros([numFeats .* numTacts, numDats]); % $$$ HERE!!!
psiMat = nan([numFeats .* numTacts, numDats]);

% Current Row
cRow  = 0;
cRow2 = 0;
doneTacts = [];
for iD = 1:size(d,1)

    cTact = d.tact{iD};

    % Stick the psi in the appropriate tactile zone [the rest is 0]
    tactRow = find(strcmp(allTacts,cTact));
    tactRow = ((tactRow - 1) .* numFeats) + 1;

    % Current data
    cDat    = d.realDat{iD};
    cDLen   = numel(cDat);

    if iD == 3
        disp('Test')
    end
    psiMat(tactRow:tactRow + numFeats - 1, [1:cDLen] + cRow) = d.psi{iD};
    
    % Store the description of the features that the weight refers to
    if ~ismember(cTact, doneTacts)
        for iPsi = 1:numFeats
            doneTacts = [doneTacts {cTact}];
            % --> DO I NEED TO PUT SOMETHING IN HERE THAT TAKES INTO
            % ACCOUNT HOW MANY SUBFEATURES THERE ARE?? OR do I do it above,
            % in FindPsiLines?
            psiLabels{cRow2 + iPsi} = ['Tact Loc: ' cTact '. ' d.psiLineDescr{iD}{iPsi}];
        end
        cRow2 = cRow2 + iPsi;
    end

    cRow = cRow + cDLen;
end
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Calculate fit quality and return best modelled data
function [sSqErr, dFit, residuals] = ErrorFun(psiMat,allDat,d,p,fitType)
% % %     % 'fitted' data: estimated values
% % %     dFit = sum( psiMat .* p(1:size(psiMat,1)) ) + p(size(psiMat,1) + 1);

    % With an offset for each experiment
    currInd = 0;
% % %     allExps = unique(d.exp);
% % %     expCol = [d.exp{:}];
% % %     for iExp = 1:numel(allExps)
% % %         expDat = [d.realDat{d.exp == iExp}];
% % %         expOffset(currInd + [1:numel(expDat)]) = p(size(psiMat,1) + iExp);
% % %         currInd = numel(expOffset);
% % %     end

    allExps = unique(d.exp);
    % Make a list of experiment numbers for each entry in the d table
    expCol = [];
    for iD = 1:size(d,1)
        expCol = [expCol, d.exp(iD) .* ones(size(d.realDat{iD})) ];
    end
    for iExp = 1:numel(allExps)
        expPos = expCol == iExp;
        expOffset(expPos) = p(size(psiMat,1) + iExp);
    end

    % $$$ Do non-linear fitting here? --> give it a setting for fittype
    switch fitType
        case 'Exponential'
            % Parameter controlling decay rate of the exponential
% % %             psiMat  = exp( psiMat .* p(size(psiMat,1) + iExp + 1) ); 
            psiMat  = exp( (psiMat - p(size(psiMat,1) + iExp + 2) ) .* (p(size(psiMat,1) + iExp + 1))  ); 
            iExp    = iExp + 2;
        case 'Sigmoid'
            psiMat  = 1 ./ (  1 + exp(-(psiMat -  p(size(psiMat,1) + iExp + 2)     )   ./   (p(size(psiMat,1) + iExp + 1))   )  ) ;
    end

    % Add a separate offset for stimuli moving away direction
    if numel(p) > size(psiMat,1) + iExp
        currInd = 0;
        for iD = 1:size(d,1)
            dDat = [d.realDat{iD}];
            if d.dir(iD) == -1
                dirOffset(currInd + [1:numel(dDat)]) = p(size(psiMat,1) + iExp + 1);
            else
                dirOffset(currInd + [1:numel(dDat)]) = 0;
            end
            currInd = numel(dirOffset);
        end
    else
        dirOffset = 0;
    end

    % 'fitted' data: estimated values
    dFit = nansum( psiMat .* p(1:size(psiMat,1)) ) + expOffset + dirOffset;

    residuals = (dFit - allDat).^2;

    sSqErr = nansum(residuals);
 

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function d = AddbinPos(d,s,dNums)
for iFitD = dNums

    % If there are multiple entries per 'location' adjust accordingly
    if iscell(d.cmPos{iFitD})

        d.binPos{iFitD} = cell(size(d.cmPos{iFitD}));
        tmpCmPos = d.cmPos{iFitD};


        for iTraj = 1:size(d.cmPos{iFitD},2) % loop through trajectories

            % The rows have to be flipped, because small == far
            tmpCmPos{1,iTraj} = -d.cmPos{iFitD}{1,iTraj};
            for iDim = 1:size(d.cmPos{iFitD},1)
                % Convert real positions to voxels
                d.binPos{iFitD}{iDim,iTraj} = round(s.clc.nearPos(iDim) + tmpCmPos{iDim,iTraj} ./ s.clc.binW);
            end

        end
    else

        tmpCmPos = d.cmPos{iFitD};

        % The rows have to be flipped, because small == far
        tmpCmPos(1,:) = -d.cmPos{iFitD}(1,:);
        % Convert real positions to voxels
        d.binPos{iFitD} = round(s.clc.nearPos + tmpCmPos ./ s.clc.binW);
        % But for rows, use FLOOR so that the closest position isn't IN the bodypart
        d.binPos{iFitD}(1,:) = floor(s.clc.nearPos(1,:) + tmpCmPos(1,:) ./ s.clc.binW);
    end
end
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Extract appropriate q values using the d.bin positions
function outQ = ExtractQ(allqVals, binPos)
outQ = [];
for iQ = 1:numel(allqVals)
% % %     if iQ == 9
% % %         disp('test')
% % %     end
    qVals = allqVals{iQ};
    
    outQ = [outQ ; cell2mat(arrayfun(@(pos)  qVals(:,binPos(1,pos),binPos(2,pos),binPos(3,pos)),...
    [1:size(binPos,2)], 'UniformOutput', false))];
end
end
% -------------------------------------------------------------------------


function [ gofScore chiSq pVal1 pVal2 ] = GoFFun( errSq, stdDevSq, k )
%GoFFun Calculates the Goodness of Fit scores

chiSq = nansum(errSq(:)./stdDevSq(:));

gofScore = (chiSq - k) ./sqrt(2.*k);

pVal1=1-chi2cdf(chiSq,k);
pVal2=1-normcdf(gofScore);
end


% -------------------------------------------------------------------------
function [fitRes] = Fit3Ddat(d, fitType)
% Calculates the optimal fit for data andmodel defined in d. 
% Outputs the fitting results


% base parameters
% one per psi part, per tactile loc,and then an offset per expt
p0 = ([150 500 .5  120 0.6 120 0.5 6 200 800 150 9 0.5  13  12.5  22 12  9.5 -37 17.9 0.1 1 -0.2 -.1 12 -70]') ; %ones([ 1, numel(unique(d.tact)).*size(d.psi{1},1) + numel(unique(d.exp))])';
lb = zeros(size(p0)); lb(1+end-numel(unique(d.exp)):end) = -100; % allow the offsets to be -ve
ub = ones(size(p0)) .* 1000;


% Set optimisation options - keep it simple
A = []; b = []; Aeq = []; beq = []; a = tic;


switch fitType
    case 'Linear'
    case 'Exponential'
        p0 = [p0;   -.06 ;    -5];
        lb = [lb;    -20 ; -1000];
        ub = [ub;      0 ;  1000];
    case 'Sigmoid'
        p0 = [p0;     -8 ;     -8];
        lb = [lb; -1000 ; -1000];
        ub = [ub;     0 ;  1000];
end

% Extra multiplier for direction [justified by expectation effec: all
% away-moving seems to have a higher baseline]
p0 = [p0; 20];
lb = [lb; -100];
ub = [ub; 1000];


% OK AWAY OFFSET IS JUSTIFIED BY AN EXPECTATION/ATTENTION EFFECT: all aways seem to
% have a higher baseline
% Extra multiplier for direction
p0 = [p0; 20];
lb = [lb; -100];
ub = [ub; 1000];

% Create successor matrix with separate features that are  depending on 
% touch location, for faster optimisation
[psiMat allDat psiLabels] = MakePsiMat(d);

% Specify optimised function
FunToOpt = @(p) ErrorFun(psiMat,allDat,d,p,fitType);

% Run optimisation
OPTIONS = optimset('TolCon',1e-10);
[p,sSqErr,exitflag,output,lambda,grad,hessian] = fmincon(FunToOpt,p0,A,b,Aeq,beq,lb,ub,[],OPTIONS);
optTime = toc(a)

% Extract optimised data
[sSqErrFinal, dFitFinal] = ErrorFun(psiMat,allDat,d,p,fitType);

% Store data in d
d = InsertModelledDat(d,dFitFinal,allDat);


% -------------------
% Make some statistical measure of goodness of fit

% First adjust all the input variance to standard error
dataVar = [d.realSTE{:}].^2;
errSq   = (allDat - dFitFinal).^2;
chiSq   = nansum(errSq(:)./dataVar(:))


% Calculate the chi square confidence interval using bootstrapping
chiSqDistr  = arrayfun(@(x) nansum(randsample(errSq(:)./dataVar(:),numel(errSq(:)),1)) ,[1:100000] );
chiSqCI     = prctile(chiSqDistr,[2.5 97.5])

N   = numel(errSq(:)); 
dof = numel(p);
k   = N - dof;

% Calculate goodness of fit
[ gofScore chiSqCheck pVal1 pVal2 ] = GoFFun( errSq, dataVar, k );


% calculate variance of residuals
resVar  = nansum(errSq(:)./ k )
% calculate log likelihood of data
logLike = -(  N                    .* log(2 .*pi .* resVar) ./ 2 )  - ...
             (1 ./ (2.* resVar))   .* nansum(errSq(:)) ;
% calculate AIC and BIC
AIC     = -2 .* logLike + dof .* 2;
BIC     = -2 .* logLike + dof .* log(N);



f.STATS.f = figure('Position',[20 20 600 300]);

plot(chi2pdf(0:250 ,k),'LineWidth',2); hold on
yLims = ylim;
plot([chiSq chiSq], yLims,'r','LineWidth',2);

xlim([0 250])

xlabel('total error (ChiSquare)')
ylabel('probability')

legend('Expected distribution of total error if data is generated by a process like the model', ...
    'Observed total error', 'Location','South')

title(['Action value is a satisfactory explanation of the data ' ...
       ' (Hypothesis cant be rejected: p = ' num2str(pVal1) ')'])

fitRes.d            = d;
fitRes.sSqErrFinal  = sSqErrFinal;
fitRes.dFitFinal    = dFitFinal;
fitRes.chiSq        = chiSq;
fitRes.chiSqDistr   = chiSqDistr;
fitRes.gofScore     = gofScore;
fitRes.chiSqCheck   = chiSqCheck;
fitRes.pVal1        = pVal1;
fitRes.pVal2        = pVal2;
fitRes.resVar       = resVar; 
fitRes.logLike      = logLike;
fitRes.AIC          = AIC;
fitRes.BIC          = BIC;

end
% -------------------------------------------------------------------------



