function [sNew] = DefaultSettings(varargin)
%DefaultSettings Sests all settings that haven't been specified to default
%values
% First input should be the new settings

% -------------------------------------------------------------------------
% Flags
s.fl.thr = 0; % threat flag
s.fl.onlyThr = 0; % only threat flag: changes value of goal to be negative
s.fl.vis = 0; % vision flag
s.fl.vpr = 0; % vision + proprioception flag
s.fl.hist = 0; % history flag
s.fl.tch = 0; % touch flag
% Flag for multiple outputs. If 1, it gives 1 output for each action. If 0,
% actions are inputs
s.fl.mo = 1;
s.fl.vid = 0; % $$$ NEeed to neaten this up!!! - now it might go crazy - don't know
s.fl.dspm = 0; % 1 means display world and 0 means no display
% This is if the relearning should only use more recent data
s.fl.rdat = 1;
% Oversample high and low reward cases
s.fl.os = 1;
s.fl.newNet = 1; % start a new network when running a simulation or relearning
s.fl.newTable = 1; % start a new table when running a simulation or relearning
s.fl.trainNet = 1;
s.fl.trainTable = 1;
s.fl.bdyMov = 0; % If this is set to 1, the body can move as well
s.fl.ToolChange = 0; % If this is 1, the tool will change size/existence from (goal) reset to reset
s.fl.rtTask = 0; % If this is 1, the agent will perform a "Reaction Time" task instead
s.fl.extraInput = 0; % If == 1, the extra input is used as input to the network
s.fl.perfWhileLearn = 0; % Set to calculate network performance after every s.prf.skipBatches batches
s.fl.deepNet = 0; % set to 1 if training to be done using deep learning toolbox
s.fl.grabAndEat = 0; % set to 1 if the goal needs to be moved to a specific location to receive reward
s.fl.defendZone = 0; % set to 1 if the agent needs to prevent the stimulus from touching a particular block

% -------------------------------------------------------------------------
% 'Stimulus' Properties
% Set possible and current threat speeds in Rows (R) and Columns (C)
s.thr.alSpR = 1;%[1 2];
s.thr.alSpC = [-1 0 1];
s.thr.randSpr = [2 2]; % random spread - if 2, the threat has a small chance to move by 1 pixel. Dimension 1 is uo/down (rows), dimension 2 is left/right (columns)
% Set possible and current goal speeds in Rows (R) and Columns (C)
s.gol.alSpR = 1; %[1 2];
s.gol.alSpC  = [-1 0 1];
s.gol.randSpr = [2 2]; % random spread - if 2, the goal has a small chance to move by 1 pixel left or right
% This is for the Reaction Time Task
s.rtt.threshold = 2;
s.rtt.mu = 0;
s.rtt.std = 1;
% Alternative option for NoiseFun: @(mu,std) (2.*std.*rand([1,1])-std+mu)
s.rtt.NoiseFun = @(mu,std) normrnd(mu,std,1,1); % This is the noise in the 'tactile' input. 
% This is for the extra info fed into the neural network
s.xtra.useThr = 0;


% -------------------------------------------------------------------------
% World Properties
s.wrld.size = [8 15];
s.wrld.resetType = 'BottomTop_InfLR'; % manner in which to reset goal/threats after they exit the world. OPTIONS: 'BottomTop','BottomTop_InfLR','Edges'
s.wrld.resetPos(1,:) =  [2 3];
try
    s.wrld.resetPos(2,:) = [2 (varargin{1}.wrld.size(2)-1)];
catch
    s.wrld.resetPos(2,:) = [2 (s.wrld.size(2)-1)];
end

% -------------------------------------------------------------------------
% Body propeties
s.bdy.size = 0;
s.bdy.pos = 0;
s.bdy.startCol = 2;
s.bdy.startRow = [];
% s.bdy.pos = floor(s.wrld.size(2)/2)+1 - floor(s.bdy.size/2);

% -------------------------------------------------------------------------
% Limb properties
s.lmb.startRow = [];
s.lmb.startCol = [];
s.lmb.size = 1;
s.lmb.ToolRows = [0];
s.lmb.ToolProb = 0.5; % Probability of tool on each trial, if toolsswitching present


% -------------------------------------------------------------------------
% Learning parameters
s.lp.alg = 'Q'; %Options: 'Q','SARSA'
s.lp.epsilon = 0.1;
s.lp.alpha = 0.6;%; 0.3;
s.lp.gamma = 0.5;
% Network learning parameters................................
s.lp.netS = [6 5 4 3];
s.lp.TrnFn = 'trainlm'; % Training function
s.lp.TrnEp = 100; % Number of training epochs
s.lp.ShwW = 0; % Show training window
s.lp.mxF = 5; % Maximum number of fails (per epoch?)
s.lp.n1stTr = 1000; % Number of gradient descents to do 1st time round
s.lp.neurTypes = 'tansig'; % types of neurons for use in default Matlab network
                 % Options: purelin, poslin [i.e. relu], tansig, logsig, softmax, hardlim, tribas, radbas
% Batch learning parameters..................................
s.lp.bSiz = 7e4; % Normal batch size
s.lp.b1Siz =  7e4; % Size of first batch
s.lp.rcSiz = s.lp.bSiz; % If using only recent data to learn, set that batch size
s.lp.retr = 2e4;% Number of actions to take before network is retrained (1 for every step)
s.lp.nOs = 30; % Number of times to oversample, if flag is set
s.lp.nOsBatch = 1; % Number of batches to oversample
s.lp.hiRewDef = 1.1; % Definition of what a 'high reward to oversample' is
% Running learning parameters..................................
s.lp.minExp = 2e6; % Minimum number of epochs before first training starts $$$ 2e6 is def better
s.lp.sWs = s.lp.minExp*2;%s.lp.minExp+1; % number of steps to take before saving workspace
s.lp.dispA = 2e4; % number of steps to take before displaying how many steps have been taken
s.lp.dispW = 1;% number of steps to take before displaying the world
s.lp.plQ = s.lp.minExp+2; % number of steps to take before plotting Q values
s.lp.maxRetrainSteps = 1000; % number of times to retrani the network while learing


% -------------------------------------------------------------------------
% Action parameters
% Policy: epsilon-greedy
s.act.pi = @(optimal_action, n_actions, randNum) (randNum > s.lp.epsilon) * optimal_action + (randNum <= s.lp.epsilon) * randi(n_actions);
% Action types
s.act.Name = {'LEFT','STAY','RIGHT'}; %  names of actions
% Set status message for goal and bump
s.act.BUMP = 2;
s.act.GOAL = 3;
s.act.THREAT = 4;
% Rewards for the various outcomes
s.act.BumpRew = -1;
s.act.GoalRew = 2;
s.act.ThreatRew = -2;
s.act.DefaultMoveRew =  -0.001;
s.act.DefaultStayRew =  0;
s.act.bdyBumpRew = -1;
s.act.bdyGoalRew = 4;
s.act.bdyThreatRew = -4;
s.act.bdyMoveRew =  -0.001;
s.act.bdyStayRew =  0;
s.act.rtRew = 0;
s.act.eatRew = 0;
s.act.shieldRew = 0;
% Value approximation function. 
% If this is table, it will use the table to find the optimum action. If it
% is 'Net', it will use a network
s.act.ValAprFun = 'Net' ; 
s.act.bdyMovesLimb = 0; % Defines whether the limb is constrained by the body's location
s.act.bdyLimbProx  = 2; % max distance between body and limb

% -------------------------------------------------------------------------
% Runtime parameters
s.rp.numIter = 1; % change this value to set max iterations
s.rp.maxActions = 2e6 ; % Maximum number of action taken before instance ends. 2e6 2e5 forr quick

% -------------------------------------------------------------------------
% Re-learning parameters
s.rl.maxRetr = 5;

% -------------------------------------------------------------------------
% Plot settings
s.plt.ON = 1; % set to zero if the polot functions shouldn't plot figures $$ incomplete in most functions so far
s.plt.meanOSVfl = 0; % set to 1 if averaging all possible other variables
s.plt.plotThrFl = 0; % 1 if plotting the effect of the threat, 0 if plotting the effect of the stimulus
s.plt.intrpFact = 1.0; % Interpolates the images
s.plt.plAct = 2; % Default action number to plot
s.plt.rowLag = 1; % lag of the visible object
s.plt.colLag = 0;
s.plt.nvRl = 1; % non-visible object lag (i.e. thr if goal is vis and vice-versa)
s.plt.nvCl = 0;
s.plt.lmbVisVal = 60; % Visual value of limb (and below, goal and threat)
s.plt.golVisVal = 100;
s.plt.thrVisVal = 10;
s.plt.lmbRow = [];
s.plt.lmbCol = 3;
s.plt.bdyCol = 4;
s.plt.OneActFl = 1; % If only the value of 1 action should be plotted by the DisplActValsFun
s.plt.vidFl = 0;
s.plt.vidFR = 1.5; % video framerate 
s.plt.startLayer  = 1 ; % First layer to plot multiple neurons from
try
    s.plt.stopLayer = length(varargin{1}.lp.netS) ; % Last layer to plot multiple neurons from
catch
    s.plt.stopLayer = length(s.lp.netS) ;
end
s.plt.vidFileName = 'D:\DPPS\DefenseAgent\Results\Videos\Test.avi';
s.plt.cBarFl = 1; % set to 0 if you don't want to plot colorbars
s.plt.pltType = 'Imagesc'; % Changes which type of plot happens. OPTIONS: 'Imagesc' 'Binned'
s.plt.distanceType = 'Absolute'; %  OPTIONS: 'Absolute' 'Row' 'Column' 'AbsRow' 'AbsColumn'
s.plt.sequentialLimbCols = 1;
s.plt.nBins = 10;
s.plt.separateLimbCols = 0; % plots output for each possible limb column if set to 1. 
s.plt.meanLimbCols = 0; % If set to 0, plots summary stats across hand columns
s.plt.rowLims = [1.5 s.wrld.size(1)-0.5]; % limits of rows to make Imagesc of
s.plt.colLims = [1.5 s.wrld.size(2)-1.5]; % limits of rows to make Imagesc of
s.plt.ToolPresent = 0;
s.plt.rtTouchVal = 0; % touch value for the reaction time task
s.plt.holdingRew = 0; % holding object value for rewarded eating
s.plt.axesVis = 1; % If zero, don't show axes
s.plt.fitSigmoid = 0; % Fit a sigmoid to Q- or Neural activity values
s.plt.DistFromTool = 0; % if 1, all correlations and plots will be done in a tool-tip centric manner
s.plt.fitSigmLow = [-Inf,0, 0, 0]; % Lower values of parameters for the sigmoid fitting function
s.plt.fitSigmUp = [Inf, Inf, Inf, Inf]; % Upper values of parameters for the sigmoid fitting function
s.plt.fitSigmStart = [-1, .5, 0, 1]; % Starting values of parameters for the sigmoid fitting function
s.plt.fitSigmDistTypes = {'aD','rD','cD','rDabs','cDabs'};
s.plt.nPm = 100; % number of perumtations for fpermutation tests

% -------------------------------------------------------------------------
% For newer figure settings
s.plt.varargs           = {};
s.plt.fS.c              = [0 0 0.7];
s.plt.fS.nBins          = 5;
s.plt.fS.pltType        = 'shadeplot';
s.plt.grid              = 0; % if 1, plots a grid over any imagesc image
s.plt.showLimb          = 0; % if is 1, plots a circle at the position of the limb
s.plt.showBodyCol       = []; % the number of the column where to plot the body
s.plt.fancyWorld        = 0; % For plotting of the actual world state
s.plt.otherStateVars    = [s.plt.lmbRow-6 5]; % For calculating neural activation with the threats/goals in additional positions
s.plt.contours          = 0;
s.plt.contourVals       = [-0.3, 0.3];

% -------------------------------------------------------------------------
% Figure settins I'm not sure about $$$
s.plt.nSubCols = 2;

% -------------------------------------------------------------------------
% Settings for graph/network analysis
s.nta.useNetWeight  = 'Norm';  % options: 'Norm','Raw','Ranked'
s.nta.comparison    = 'Valence'; % options: 'Valence','Bodypart'
s.nta.compScale     = 'Graded'; % options: 'Graded','Absolute'

% -------------------------------------------------------------------------
% Performance testing parameters
s.prf.maxActions    = 1e3;
s.prf.newNet        = 0;
s.prf.newTable      = 1;
s.prf.epsilon       = 0;
s.prf.nRep          = 100;
s.prf.skipBatches   = 1;

% -------------------------------------------------------------------------
% Video plot variables
% % % s.vid.fancyPlot = 0;


% -------------------------------------------------------------------------
% Data fitting variables
s.clc.startRew              =  1;
s.clc.nearPos               = [s.wrld.size(1)-12 11 13]';
s.clc.startSR               = 12; % Limb row
s.clc.startSC               =  8; % Limb column
s.clc.startSZ               =  1; % Limb hieght (not important in 2D)
s.clc.nReps                 = 1; 
s.clc.stepUpdateFl          = 0; % Whether to update in timesteps - especially important for hitprob and multisens integration
s.clc.RewardBehindSurfaceFl = 1; % Whether or not to have more depth to the body
s.clc.checkCollisionFl      = 1; % Whether to check if a collision has happened anywhere along the line
s.clc.useAltPolicyFl        = 0; % Whether to use the alternate Q-values supplied as the policy
s.clc.nSteps                = 1; 
s.clc.gammaVal              = 0.8;
s.clc.baseVel               = [4  0  0]; 
s.clc.actConsequence        = [0  0  0];
s.clc.maximiseSimilarityType = 'OverallQ'; % Alternatives: 'ChosenAction', For Policy improvement (e.g. egocentric map use).  

% Random stimulus dynamics
rSpr    = [-1 0 1];
rSprPr  = gaussmf(rSpr,[1 0]) ./ sum(gaussmf(rSpr,[1 0]));
cSpr    = [-1 0 1]; %%% cSpr = [ -1 0 1 ];
cSprPr  = gaussmf(cSpr,[1 0]) ./ sum(gaussmf(cSpr,[1 0]));
zSpr    = [-1 0 1]; %%% zSpr = [ -1 0 1 ];
zSprPr  = gaussmf(zSpr,[1 0]) ./ sum(gaussmf(zSpr,[1 0]));
% x y z, Deterministic stimulus dynamics
s.clc.stimDynams =     @(pos) pos + s.clc.baseVel; % For approaching, set speed positive
s.clc.randSpread =     {rSpr cSpr zSpr}; % Put a little biit of x and z variability in? Kind of arbitrary
s.clc.spreadProb =     {rSprPr cSprPr  zSprPr}; % x y z, probabilities of spread
% Random sensory uncertainties
s.clc.sensSpread = {[0] , [0] , [0]}; 
s.clc.sensProb   = {[1] , [1] , [1]}; 






% $$ STILL EXPERIMENTAL
s.plt.stimRow = [];
s.plt.stimCol = [];

% =========================================================================
% Substitute the newly specified settings for the default ones
if numel(varargin)>0
    sNew = varargin{1};
    
    sFields = fieldnames(s);
    
    for iFN = 1:length(sFields)  %field name counter
        
        % If thesubstrucure doesn't exist, put it there in its entirety
        if ~isfield(sNew,sFields{iFN})
            sNew.(sFields{iFN}) = s.(sFields{iFN});
                      
        else
            % If the substructure does exist, try looping through the
            % sub-structure, to check whether there are any missing bits
            sSubFields = fieldnames(s.(sFields{iFN}));
            
            for iSFN = 1:length(sSubFields)

                if ~isfield(sNew.(sFields{iFN}),sSubFields{iSFN})
                    
                    sNew.(sFields{iFN}).(sSubFields{iSFN}) = ...
                        s.(sFields{iFN}).(sSubFields{iSFN});
                end
            end
        end
        
    end
    
else
    sNew = s;
end



end

