% Sensory stimulus spreading: from Straka et al ( * 1/2 because their units
% were seconds, and ours are half-seconds
rSnsSpr = [-2:2]; 
cSnsSpr = [-2:2]; 
zSnsSpr = [-2:2]; 

% Set uncertainties to Straka Noel Hoffmann values
rSigma      = sqrt(5.^2 + 5.^2) ./ sBdyHP.clc.binW ./2;
rSnsSprPr   = gaussmf(rSnsSpr,[rSigma 0]) ./ sum(gaussmf(rSnsSpr,[rSigma 0])); %sqrt(2.5.^2 + 30.^2)./s.clc.binW
czSigma     = sqrt(5.^2 + 5.^2) ./ sBdyHP.clc.binW ./2;
cSnsSprPr   = gaussmf(cSnsSpr,[czSigma 0]) ./ sum(gaussmf(cSnsSpr,[czSigma 0]));
zSnsSprPr   = gaussmf(zSnsSpr,[czSigma 0]) ./ sum(gaussmf(zSnsSpr,[czSigma 0]));

doubleUnc = 0;


% -------------------------------------------------------------------------

% #################################################################
% =========================================================
% BODY TOWARDS
sBdy.clc.stimDynams =     @(pos) pos + sBdy.clc.baseVel;

sBdy.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sBdy.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

sBdy.clc.nReps = 5; 
sBdy.clc.stepUpdateFl = 0; % whether to run a full sweep update or not
sBdy.clc.nSteps = 1;



% -----------------------------
% POS REW
% REWARD SIZE
sBdy.clc.startRew   =  1;
% % % [newQ ] = CalcQDirect(sBdy);
[newQ ] = CalcQDirect(sBdy);

if doubleUnc == 1
[newQ ] = ObserverPerspective(sBdy,newQ);
end
% store q values and attributes
cQ = 78;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sBdy.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'TrunkUnc'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central



% -----------------------------
% NEG REW
% REWARD SIZE
sBdy.clc.startRew =  -1;
[newQ ] = CalcQDirect(sBdy);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sBdy,newQ);
end
% % % [newQ ] = CalcQDirect(sBdy,allQ.qVals{2});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sBdy.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'TrunkUnc'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% =========================================================
% BODY AWAY
sBdy.clc.stimDynams =     @(pos) pos - sBdy.clc.baseVel; % For receding, set speed negative

% -----------------------------
% POS REW
% REWARD SIZE
sBdy.clc.startRew =  1;
[newQ ] = CalcQDirect(sBdy);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sBdy,newQ);
end
% % % [newQ ] = CalcQDirect(sBdy,allQ.qVals{3});
% store q values and attributes
cQ = cQ +1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sBdy.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'TrunkUnc'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

% -----------------------------
% NEG REW
% REWARD SIZE
sBdy.clc.startRew =  -1;
[newQ ] = CalcQDirect(sBdy);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sBdy,newQ);
end
% % % [newQ ] = CalcQDirect(sBdy,allQ.qVals{4});
% store q values and attributes
cQ = cQ + 1
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sBdy.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'TrunkUnc'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

%Intermediate save
% save('F:\Projects\DPPS\DefenseAgent\Data\NewWithHitProbs_MultSensInt_AndDist_AndUncertainQ.mat')


% #################################################################
% =========================================================
% HEAD TOWARDS
sHed.clc.stimDynams =     @(pos) pos + sHed.clc.baseVel;

sHed.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sHed.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

%% -----------------------------
% POS REW
% REWARD SIZE
sHed.clc.startRew =  1;
[newQ ] = CalcQDirect(sHed);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHed,newQ);
end
% % % [newQ ] = CalcQDirect(sHed,allQ.qVals{5});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}        = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHed.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HeadUnc'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW
% REWARD SIZE
sHed.clc.startRew =  -1;
[newQ ] = CalcQDirect(sHed);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHed,newQ);
end
% % % [newQ ] = CalcQDirect(sHed,allQ.qVals{6});
cQ = cQ + 1;
% store q values and attributes
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHed.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HeadUnc'; 
allQ.centerpos{cQ}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% =========================================================
% HEAD AWAY
sHed.clc.stimDynams  =     @(pos) pos - sHed.clc.baseVel; % For receding, set speed negative

% -----------------------------
% POS REW $$$
% REWARD SIZE
sHed.clc.startRew =  1;
[newQ ] = CalcQDirect(sHed); 
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHed,newQ); % $$$ UNCOMMENT THIS
end
% % % [newQ ] = CalcQDirect(sHed,allQ.qVals{7}); % $$$ COMMENT THIS
% $$$ ADD CQ IS CQ PLUS ONE
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ; % $$$ CQ HERE AND BELOW
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHed.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HeadUnc';  % $$$ HEAD UNCERTAIN
allQ.centerpos{cQ}   = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

% -----------------------------
% NEG REW
% REWARD SIZE
sHed.clc.startRew   =  -1;
[newQ ] = CalcQDirect(sHed);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHed,newQ);
end
% % % [newQ ] = CalcQDirect(sHed,allQ.qVals{8});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}      = newQ;
allQ.dir(cQ)        = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)        = sHed.clc.startRew; % towards
allQ.bodyPart{cQ}   = 'HeadUnc'; 
allQ.centerpos{cQ}  = sBdy.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%Intermediate save
% save('F:\Projects\DPPS\DefenseAgent\Data\NewWithHitProbs_MultSensInt_AndDist_AndUncertainQ.mat')


% #################################################################
% =========================================================
% HAND by chest, TOWARDS

sHnd.clc.stimDynams = @(pos) pos + sHnd.clc.baseVel;

sHnd.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sHnd.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

%-----------------------------
% POS REW
% REWARD SIZE
sHnd.clc.startRew =  1;
[newQ ] = CalcQDirect(sHnd);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHnd,newQ);
end
% % % [newQ ] = CalcQDirect(sHnd,allQ.qVals{9});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHnd.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestUnc'; 
allQ.centerpos{cQ}   = sHnd.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW
% REWARD SIZE
sHnd.clc.startRew    =  -1;
[newQ ] = CalcQDirect(sHnd);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHnd,newQ);
end
% % % [newQ ] = CalcQDirect(sHnd,allQ.qVals{10});
% store q values and attributes
cQ = cQ +1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHnd.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestUnc'; 
allQ.centerpos{cQ}   = sHnd.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% =========================================================
% HAND AWAY
sHnd.clc.stimDynams =     @(pos) pos - sHnd.clc.baseVel; % For receding, set speed negative


% -----------------------------
% POS REW
% REWARD SIZE
sHnd.clc.startRew =  1;
[newQ ] = CalcQDirect(sHnd);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHnd,newQ);
end
% % % [newQ ] = CalcQDirect(sHnd,allQ.qVals{11});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHnd.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestUnc'; 
allQ.centerpos{cQ}   = sHnd.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

% -----------------------------
% NEG REW
% REWARD SIZE
sHnd.clc.startRew =  -1;
[newQ ] = CalcQDirect(sHnd);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHnd,newQ);
end
% % % [newQ ] = CalcQDirect(sHnd,allQ.qVals{12});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHnd.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestUnc'; 
allQ.centerpos{cQ}   = sHnd.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% #################################################################
% ===============================================s==========
% HAND by side, all directions 

% Move the data, then shift the Q values next
cQ = cQ + [1:4];
allQ(cQ,:)       = allQ(cQ - 4,:);
allQ.bodyPart(cQ)    = repmat({'HandBySideUnc'},[1 4])';



% Shift the data
for iPsi = cQ
    allQ.qVals{iPsi} = zeros(size(allQ.qVals{iPsi}));

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



%Intermediate save
% save('F:\Projects\DPPS\DefenseAgent\Data\NewWithHitProbs_MultSensInt_AndDist_AndUncertainQ.mat')



%% #################################################################
% TOOL TIME!
% =========================================================
% TOOL by HAND, TOWARDS

sHndTool.clc.stimDynams = @(pos) pos + sHndTool.clc.baseVel;

sHndTool.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sHndTool.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

%-----------------------------
% POS REW
% REWARD SIZE

sHndTool.clc.startRew =  1;
[newQ] = CalcQDirect(sHndTool);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHndTool,newQ);
end
% % % [newQ ] = CalcQDirect(sHndTool,allQ.qVals{17});
% store q values and attributes
% cQ = cQ(end) + 1; % $$$ HERE it goes wrong for some reason @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
cQ = 94;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndTool.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolUnc'; 
allQ.centerpos{cQ}   = sHndTool.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW
% REWARD SIZE

sHndTool.clc.startRew =  -1;
[newQ] = CalcQDirect(sHndTool);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHndTool,newQ);
end
% % % [newQ ] = CalcQDirect(sHndTool,allQ.qVals{18});
% store q values and attributes
cQ = cQ + 1
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndTool.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolUnc'; 
allQ.centerpos{cQ}   = sHndTool.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

% -----------------------------
% POS REW, AWAY
% REWARD SIZE
sHndTool.clc.stimDynams =     @(pos) pos - sHndTool.clc.baseVel;

sHndTool.clc.startRew =  1;
[newQ] = CalcQDirect(sHndTool);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHndTool,newQ);
end
% % % [newQ ] = CalcQDirect(sHndTool,allQ.qVals{19});
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndTool.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolUnc'; 
allQ.centerpos{cQ}   = sHndTool.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW, AWAY
% REWARD SIZE

sHndTool.clc.startRew =  -1;
[newQ] = CalcQDirect(sHndTool);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHndTool,newQ);
end
% % % [newQ ] = CalcQDirect(sHndTool,allQ.qVals{20});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndTool.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolUnc'; 
allQ.centerpos{cQ}   = sHndTool.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%Intermediate save
% save('F:\Projects\DPPS\DefenseAgent\Data\NewWithHitProbs_MultSensInt_AndDist_AndUncertainQ.mat')



% #################################################################
% Spider Butterfly Track; let's just assume the uncertainties are too small
% to affect this one, so I don't need a new Q-value calculation for it
% ===============================================s==========



% #################################################################
% TOOL TIME version 2: rake
% ===============================================s==========
% TOOL by HAND, TOWARDS

sHndToolRake.clc.stimDynams = @(pos) pos + sHndToolRake.clc.baseVel;

sHndToolRake.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sHndToolRake.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

%-----------------------------
% POS REW
% REWARD SIZE

sHndToolRake.clc.startRew =  1;
[newQ ] = CalcQDirect(sHndToolRake);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHndToolRake,newQ);
end
% % % [newQ ] = CalcQDirect(sHndToolRake,allQ.qVals{27});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndToolRake.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeUnc'; 
allQ.centerpos{cQ}   = sHndToolRake.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW
% REWARD SIZE

sHndToolRake.clc.startRew =  -1;
[newQ ] = CalcQDirect(sHndToolRake);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHndToolRake,newQ);
end
% % % [newQ ] = CalcQDirect(sHndToolRake,allQ.qVals{28});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndToolRake.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeUnc'; 
allQ.centerpos{cQ}   = sHndToolRake.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central

%  -----------------------------
% POS REW, AWAY
% REWARD SIZE

sHndToolRake.clc.stimDynams = @(pos) pos - sHndToolRake.clc.baseVel;

sHndToolRake.clc.startRew =  1;
[newQ ] = CalcQDirect(sHndToolRake);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHndToolRake,newQ);
end
% % % [newQ ] = CalcQDirect(sHndToolRake,allQ.qVals{29});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndToolRake.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeUnc'; 
allQ.centerpos{cQ}   = sHndToolRake.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%% -----------------------------
% NEG REW, AWAY
% REWARD SIZE

sHndToolRake.clc.startRew =  -1;
[newQ ] = CalcQDirect(sHndToolRake);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHndToolRake,newQ);
end
% % % [newQ ] = CalcQDirect(sHndToolRake,allQ.qVals{30});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = -1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHndToolRake.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HandByChestPlusToolRakeUnc'; 
allQ.centerpos{cQ}   = sHndToolRake.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


%Intermediate save
% save('F:\Projects\DPPS\DefenseAgent\Data\NewWithHitProbs_MultSensInt_AndDist_AndUncertainQ.mat')


% #################################################################
% Monkey arm
% ===============================================s==========
% ARM FORWARD ['right']

sArm.clc.stimDynams = @(pos) pos + sArm.clc.baseVel;

sArm.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sArm.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};


% -----------------------------
% POS REW
% REWARD SIZE

sArm.clc.startRew =  1;
[newQ ] = CalcQDirect(sArm);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sArm,newQ);
end
% % % [newQ ] = CalcQDirect(sArm,allQ.qVals{31});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sArm.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'ArmForwardUnc'; 
allQ.centerpos{cQ}   = sArm.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central



% -----------------------------
% NEG REW
% REWARD SIZE
sArm.clc.startRew =  -1;
[newQ ] = CalcQDirect(sArm);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sArm,newQ);
end
% % % [newQ ] = CalcQDirect(sArm,allQ.qVals{32});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sArm.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'ArmForwardUnc'; 
allQ.centerpos{cQ}   = sArm.clc.nearPos;



% ===============================================s==========
% ARM LEFT

sArmL.clc.stimDynams = @(pos) pos + sArmL.clc.baseVel;

sArmL.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sArmL.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};


%-----------------------------
% POS REW
% REWARD SIZE

sArmL.clc.startRew =  1;
[newQ ] = CalcQDirect(sArmL);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sArmL,newQ);
end
% % % [newQ ] = CalcQDirect(sArmL,allQ.qVals{31});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sArmL.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'ArmLeftUnc'; 
allQ.centerpos{cQ}   = sArmL.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central



% -----------------------------
% NEG REW
% REWARD SIZE
sArmL.clc.startRew =  -1;
[newQ ] = CalcQDirect(sArmL);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sArmL,newQ);
end
% % % [newQ ] = CalcQDirect(sArmL,allQ.qVals{34});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sArmL.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'ArmLeftUnc'; 
allQ.centerpos{cQ}   = sArmL.clc.nearPos;



% #################################################################
% =========================================================
% MONKEY HEAD TOWARDS - CONSTRAINED

sHedConstr.clc.stimDynams = @(pos) pos + sHedConstr.clc.baseVel;

sHedConstr.clc.sensSpread = {rSnsSpr cSnsSpr zSnsSpr};
sHedConstr.clc.sensProb   = {rSnsSprPr cSnsSprPr zSnsSprPr};

% -----------------------------
% POS REW
% REWARD SIZE
sHedConstr.clc.startRew =  1;
[newQ ] = CalcQDirect(sHedConstr);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHedConstr,newQ);
end
% % % [newQ ] = CalcQDirect(sHedConstr,allQ.qVals{35});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHedConstr.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HeadConstrUnc'; 
allQ.centerpos{cQ}   = sHedConstr.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central


% -----------------------------
% NEG REW
% REWARD SIZE
sHedConstr.clc.startRew =  -1;
[newQ ] = CalcQDirect(sHedConstr);
if doubleUnc == 1
[newQ ] = ObserverPerspective(sHedConstr,newQ);
end
% % % [newQ ] = CalcQDirect(sHedConstr,allQ.qVals{36});
% store q values and attributes
cQ = cQ + 1;
allQ.qVals{cQ}       = newQ;
allQ.dir(cQ)         = 1; %s.clc.stimDynams([0 0 0]); % towards
allQ.rew(cQ)         = sHedConstr.clc.startRew; % towards
allQ.bodyPart{cQ}    = 'HeadConstrUnc'; 
allQ.centerpos{cQ}   = sHedConstr.clc.nearPos; % x y z position of this plot in the OVERALL space --> trunk is central



% #################################################################
% =========================================================
% 'MONKEY' HEAD 15 degree angle

cQ= cQ + [1:2];

% Move the data, then shift the Q values next
allQ(cQ,:)       = allQ(cQ - 2,:);


allQ.bodyPart(cQ)    = repmat({'RotatedHeadUnc'},[1 2])';

% $$$ HERE 

% Shift the data
for iPsi = cQ

    % shift the q values to the side by 15 degree, i.e. once for every 4
    % rows, starting at the face
    shiftAmount = 1;
% %     for iR = (sHed15.clc.nearPos(1) + 2) : -3 : 4
% %         allQ.qVals{iPsi}(:,iR:-1:iR-2, 1:end-shiftAmount ,:) = allQ.qVals{iPsi-30}(:,iR:-1:iR-2, (shiftAmount+1):end ,:);
% %         shiftAmount = shiftAmount + 1;
% %     end

for iR = (sHed15.clc.nearPos(1) + 2) : -1 : 1
    allQ.qVals{iPsi}(:,iR, 1:end-shiftAmount ,:) = allQ.qVals{iPsi-2}(:,iR, (shiftAmount+1):end ,:);
    
    if mod(iR - (sHed15.clc.nearPos(1) + 3),4) == 0
        shiftAmount = shiftAmount + 1;
    end
end

            
end


% Hopefully final save
save('F:\Projects\DPPS\DefenseAgent\Data\NewWithHitProbs_MultSensInt_AndDist_AndUncertainQ_andSARSAQ_V3.mat')


%% Functions 

function observerOutput = ObserverPerspective(s,qVals)
% observerOutput = ObserverPerspective(s,qVals)
% Inverts the observerdynamics, so that it calculates Q or hitprob from the
% perspective of an observer delivering the stimuli
% inverse dynamics for P(o|s)

sobs = s;
sobs.clc.sensSpread = cellfun(@(x) -x, s.clc.sensSpread, 'UniformOutput', false); 
% no movement, only observation uncertainty allowed
sobs.clc.baseVel           = [0 0 0]; 
sobs.clc.actConsequence    = [0 0 0]; 
sobs.clc.randSpread        = {[0] [0] [0]};
sobs.clc.spreadProb        = {[1] [1] [1]};
sobs.clc.stimDynams        = @(pos) pos + sobs.clc.baseVel;
sobs.clc.nSteps            = 1;
sobs.clc.nReps             = 1; 
sobs.clc.stepUpdateFl      = 1; 

% hit encoding from observer's perspective
observerOutput = CalcQDirect(sobs,qVals);
end