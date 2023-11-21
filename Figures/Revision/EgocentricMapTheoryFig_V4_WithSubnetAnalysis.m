
load('Results\ForFigures\Fig1_Results_v3')
s = DefaultSettings(rS(end).s);
% w = rS(end).w;
%
%% Load subnets


% % % load('Results\ForFigures\SuccessorState\NetworkAnalysisWorkSpace2.mat')
% % % clear allBinNPFpm allBinNPFpm

% all Threat Presences to test the extents with. 0 is just goal, 1 is just
% 'threat' (which could also have positive valence)
aTP = [0 1];


squishAct = nan([14 15 7 18 2 (length(rS)-1)]);

for iM = 1:1 %length(rS)
    for iV = 1:2 %length(aTP)
    
    
        sFP = DefaultSettings(rS(iM).s); 
        sFP.plt.otherStateVars = 3;
        w = rS(iM).w;
        net = rS(iM).net;
        Qtable = rS(iM).Qtable;
        
        sFP.plt.plotThrFl = aTP(iV);
           
        sFP.plt.rowLims = [1.5 sFP.wrld.size(1)-0.5];
        sFP.plt.stimRow=[3:size(w.world2D,1)-3];
        
        % -------------------------------------------------------------------------
        % Calculate the correlati-ness of each velocity, for Q-values and neural
        % activity
        sFP.plt.meanLimbCols = 1;
        sFP.plt.fitSigmoid = 0;
        sFP.plt.lmbCol = [2:(sFP.wrld.size(2)-1)];
        sFP.plt.stimCol= [2:(sFP.wrld.size(2)-1)] ;
        [Q,allNeurAct] = CalcNetOutput(sFP,w,net);

        sFP.plt.plAct = 1;
        sFP.plt.ON = 0;
        for iL = 1:size(allNeurAct,5)
            for iN = 1:size(allNeurAct,6)
                [dmy, squishAct(:,:,iL,iN,iV,iM)] = DisplActValsFun(sFP,w,permute(allNeurAct(:,:,:,:,iL,iN), [3 4 1 2 5 6]));
            end
        end

        % Divide the neurons up in different ways
        % (first make sure I'm only taking the 2nd half of the network
        netLastHalf = ceil(size(allNeurAct,5)./2) : size(allNeurAct,5);
        tmpABrat    = allNetAR(iM).A_B_rat(netLastHalf,:);
        
        % General selectivity rankings
        tmpNeurs            = allNetAR(iM).A_B_rat(netLastHalf,:);
        toSort = abs(tmpNeurs - 0.5);
        [dmy sortInds]       = sort(toSort(:));
        [dmy selectRanks]       = sort(sortInds);
        selectRanks = reshape(selectRanks,size(tmpNeurs));

        % Most unselective
        % $$$ HERE DO THIS NEXT



    end
end

allNetAR(iM).A_B_rat
squishAct(:,:,iL,iN,iV,iM)


%%


allGammas = [0.7];

iGamma = 1;


% New settings - let's see if this make things look neater
s.clc.RewardBehindSurfaceFl = 0;
s.clc.checkCollisionFl      = 0;

s.clc.maximiseSimilarityType = 'OverallQ'; %'WinningQ' ; %'ChosenAction'; % 'OverallQ'




% settings for plot
sFP=s;
sFP.plt.lmbRow = s.wrld.size(1)-2;
sFP.plt.rowLims=[6.5 s.wrld.size(1)-1.5];
sFP.plt.colLims=[3.5 s.wrld.size(2)-3.5];
sFP.plt.cBarFl=0;
sFP.plt.meanLimbCols=1;
sFP=DefaultSettings(sFP);
sFP.plt.axesVis=1;

fS.gridXstart = -4.5;
fS.gridXstep = 1;
fS.gridYstart = 3.5;
fS.gridYstep = 1;



sFP.plt.lmbCol=8;


% Data fitting variables
s.clc.startRew =  1;
s.clc.startSR = 12;
s.clc.startSC =  8;
s.clc.startSZ =  1;
s.clc.nearPos = [s.wrld.size(1)-0 8 1]';
s.clc.nReps = 1;

s.clc.stepUpdateFl = 0; % Whether to update in timesteps - especially important for hitprob and multisens integration
s.clc.nSteps = 1;

s.clc.gammaVal   = allGammas(iGamma);
s.clc.baseVel    = [1 0 0];

% Random stimulus dynamics
rSpr = 0;
rSprPr = 1;
% % % rSpr = [-1 0 1];
% % % rSprPr = 1./3 .* [1 1 1];
% % % % rSprPr = [.1 .8 .1];
% % % rSpr = [0 1];
% % % rSprPr = 1./2 .* [1 1];


cSpr = 0;
cSprPr = 1;
% % % cSpr = [-1 0 1];
% % % % cSprPr = 1./3 .* [1 1 1];
% % % cSprPr = [.2 .6 .2];

zSpr = 0; %%% zSpr = [ -1 0 1 ];
zSprPr = 1;
% x y z, Deterministic stimulus dynamics
s.clc.stimDynams =     @(pos) pos + s.clc.baseVel; % For approaching, set speed positive
s.clc.randSpread =     {rSpr cSpr zSpr}; % Put a little biit of x and z variability in? Kind of arbitrary
s.clc.spreadProb =     {rSprPr cSprPr  zSprPr}; % x y z, probabilities of spread
% Random sensory uncertainties
s.clc.sensSpread = {[0] , [0] , [0]};
s.clc.sensProb   = {[1] , [1] , [1]};

% % % % s.clc.actConsequence = ...
% % % %    [0  0  0 ; ... % action 1 stay
% % % %     0  1  0 ; ... % action 2 left
% % % %     0 -1  0 ; ... % action 3 right
% % % %     1  0  0 ; ... % action 4 up
% % % %    -1  0  0 ; ... % action 5 down
% % % %     1  1  0 ; ... % action 6 left up
% % % %     1 -1  0 ; ... % action 7 right up
% % % %    -1  1  0 ; ... % action 8 left down
% % % %    -1 -1  0 ]; % action 4 right down


% % % s.clc.actConsequence = ...
% % %    [0  0  0 ; ... % action 1 stay
% % %     0  1  0 ; ... % action 2 left
% % %     0 -1  0 ; ... % action 3 right
% % %     1  0  0 ; ... % action 4 up
% % %    -1  0  0 ];    % action 5 down

s.clc.actConsequence = ...
    [0  0  0 ; ... % action 1 stay
    0  1  0 ; ... % action 2 left
    0 -1  0];     % action 3 right


fS.cAxis = allGammas(iGamma) .* [-1 1] .* max(cSprPr);
% fS.cAxis = allGammas(iGamma) .* [-1 1] ./3;
% fS.cAxis = allGammas(iGamma) .* [-.7 .7];


s.wrld.size = [14 15 1];

clear baseTasks

% POLICY: grab goal
% --------------------------------
% Task: grab goal
s.clc.useAltPolicyFl = 0;
s.clc.startRew = 1;
polGtaskGQtmp = CalcHPDirect(s);
polGtaskGQ = repmat(permute(polGtaskGQtmp,[4 5 2 3 1]),[14 15 1 1]);
baseTasks{1,1}.allQ = polGtaskGQ;
baseTasks{1,1}.name = 'Pol: Goal, Task: Goal';

% Task: avoid threat
polGtaskTQ = -polGtaskGQ;
baseTasks{1,2}.allQ = polGtaskTQ;
baseTasks{1,2}.name = 'Pol: Goal, Task: Threat';

% Task: wide goal
sTmp = s;
sTmp.clc.useAltPolicyFl = 1;
sTmp.clc.startSC =  [ 7  8  9];
sTmp.clc.startSR =  [12 12 12];
sTmp.clc.startSZ =  [ 1  1  1];
sTmp.clc.startRew = 1;
polGtaskWideGQtmp = CalcHPDirect(sTmp, [], polGtaskGQtmp);
polGtaskWideGQ = repmat(permute(polGtaskWideGQtmp,[4 5 2 3 1]),[14 15 1 1]);
baseTasks{1,3}.allQ = polGtaskWideGQ;
baseTasks{1,3}.name = 'Pol: Goal,Task: Wide Goal';


% POLICY: avoid threat
% --------------------------------
% Task: avoid threat
s.clc.startRew = -1;
polTtaskTQtmp = CalcHPDirect(s);
polTtaskTQ = repmat(permute(polTtaskTQtmp,[4 5 2 3 1]),[14 15 1 1]);
baseTasks{2,2}.allQ = polTtaskTQ;
baseTasks{2,2}.name = 'Pol: Threat, Task: Goal';

% Task: grab goal
polTtaskGQ = -polTtaskTQ;
baseTasks{2,1}.allQ = polTtaskGQ;
baseTasks{2,1}.name = 'Pol: Threat, Task: Threat';

% Task: wide goal
sTmp = s;
sTmp.clc.useAltPolicyFl = 1;
sTmp.clc.startSC =  [ 7  8  9];
sTmp.clc.startSR =  [12 12 12];
sTmp.clc.startSZ =  [ 1  1  1];
sTmp.clc.startRew = 1;
polTtaskWideGQtmp = CalcHPDirect(sTmp, [], polTtaskTQtmp);
polTtaskWideGQ = repmat(permute(polTtaskWideGQtmp,[4 5 2 3 1]),[14 15 1 1]);
baseTasks{2,3}.allQ = polTtaskWideGQ ; %widepolGtaskTQ(:,:,:,:,1);
baseTasks{2,3}.name = 'Pol: Threat,Task: Wide Goal';



% POLICY: grab WIDE goal
% --------------------------------
% Task: wide goal
sTmp = s;
sTmp.clc.useAltPolicyFl = 0;
sTmp.clc.startSC =  [ 7  8  9];
sTmp.clc.startSR =  [12 12 12];
sTmp.clc.startSZ =  [ 1  1  1];
sTmp.clc.startRew = 1;
polWideGtaskWideGQtmp = CalcHPDirect(sTmp);
polWideGtaskWideGQ = repmat(permute(polWideGtaskWideGQtmp,[4 5 2 3 1]),[14 15 1 1]);
baseTasks{3,3}.allQ = polWideGtaskWideGQ;
baseTasks{3,3}.name = 'Pol: Wide Goal,Task: WideGoal';



% Task: grab goal
s.clc.useAltPolicyFl = 1;
s.clc.startRew = 1;
polWideGtaskGQtmp = CalcHPDirect(s,[], polWideGtaskWideGQtmp);
polWideGtaskGQ = repmat(permute(polWideGtaskGQtmp,[4 5 2 3 1]),[14 15 1 1]);
baseTasks{3,1}.allQ = polWideGtaskGQ;
baseTasks{3,1}.name = 'Pol: Wide Goal, Task: Goal';

% Task: avoid threat
s.clc.useAltPolicyFl = 1;
s.clc.startRew = -1;
polWideGtaskTQtmp = CalcHPDirect(s,[], polWideGtaskWideGQtmp);
polWideGtaskTQ = repmat(permute(polWideGtaskTQtmp,[4 5 2 3 1]),[14 15 1 1]);
baseTasks{3,2}.allQ = polWideGtaskTQ;
baseTasks{3,2}.name = 'Pol: Wide Goal, Task: Threat';




% Reset alternative policy use
s.clc.useAltPolicyFl = 0;


% Create baseQ (matrix of basis tasks)
clear baseQ
for iPol = 1:size(baseTasks,1)
    qTmpForPol = arrayfun(@(iTask) baseTasks{iPol,iTask}.allQ, 1:size(baseTasks,2), 'UniformOutput', false);
    % row col row col act pol task
    baseQ(:,:,:,:,:,iPol,:) = cat(6,qTmpForPol{:});
end

% % % % Create matrix of basis tasks
% % % baseQ = arrayfun(@(iTask) baseTasks{iTask}.allQ, 1:numel(baseTasks), 'UniformOutput', false);
% % % baseQ = cat(6,baseQ{:});


% Define new Q-values to be created
% -------------------------------------------------------

% % % % Task: 'stay in place for 2 timesteps, then grab'
% % % stationaryQ              = zeros(size(polGtaskGQ));
% % % stationaryQ(:,:,1:12,:)  = allGammas(iGamma).^2 .* polGtaskGQ(:,:,3:14,:);
% % % stationaryQ(:,:,7:11,8) = polGtaskGQ(:,:,7:11,8);
% % % allQtoRecreate = stationaryQ(:,:,:,:,1);


% Task: 'stay in place for 1 timestep, then grab'
stationaryQ              = zeros(size(polGtaskGQ(:,:,:,:,:)));
stationaryQ(:,:,1:13,:,:)  = allGammas(iGamma).^1 .* polGtaskGQ(:,:,2:14,:,:);
stationaryQ(:,:,6:11,8,:)  = polGtaskGQ(:,:,6:11,8,:);
tasksToRecreate{1}.allQtoRecreate = stationaryQ; %stationaryQ(:,:,:,:,1);
tasksToRecreate{1}.name  = 'Stay1ThenGrab';

% % % % Task: 'stay in place for 2 timesteps, then grab'
% % % stationaryQ              = zeros(size(polGtaskGQ(:,:,:,:,:)));
% % % stationaryQ(:,:,1:12,:,:)  = allGammas(iGamma).^1 .* polGtaskGQ(:,:,3:14,:,:);
% % % stationaryQ(:,:,6:11,8,:)  = polGtaskGQ(:,:,6:11,8,:);
% % % tasksToRecreate{1}.allQtoRecreate = stationaryQ; %stationaryQ(:,:,:,:,1);
% % % tasksToRecreate{1}.name  = 'Stay2ThenGrab';


% Task: Avoid WIDE stimulus
s.clc.startSC =  [ 7  8  9];
s.clc.startSR =  [12 12 12];
s.clc.startSZ =  [ 1  1  1];

s.clc.startRew = -1;
widepolGtaskTQ = CalcHPDirect(s);
widepolGtaskTQ = repmat(permute(widepolGtaskTQ,[4 5 2 3 1]),[14 15 1 1]);
tasksToRecreate{2}.allQtoRecreate = widepolGtaskTQ; %widepolGtaskTQ(:,:,:,:,1);
tasksToRecreate{2}.name = 'WideThreat';
s.clc.startSC =  8;

% Task: Pass 2 blocks by the right of the stimulus
stationaryQ              = zeros(size(polGtaskGQ(:,:,:,:,:)));
stationaryQ(:,:,:,3:15,:)  = polGtaskGQ(:,:,:,1:13,:);
stationaryQ(:,:,:,1,:)     = polGtaskGQ(:,:,:,15,:);
tasksToRecreate{3}.allQtoRecreate = stationaryQ;  %stationaryQ(:,:,:,:,1);
tasksToRecreate{3}.name  = 'PassByRight2blocks';

% Task: flickering stimulus
flickerQ                 = polGtaskGQ(:,:,:,:,1);
% flickerQ(:,:,1:2:end,:)  = widepolGtaskTQ(:,:,1:2:end,:,1); % polGtaskTQ(:,:,1:2:end,:,1);
flickerQ(:,:,1:2:end,:)  = polTtaskTQ(:,:,1:2:end,:,1); % polGtaskTQ(:,:,1:2:end,:,1);
tasksToRecreate{4}.allQtoRecreate = flickerQ; %flickerQ(:,:,:,:,1);
tasksToRecreate{4}.name     = 'FlickerStim';



% Plot the original Q-values
f.EgocentricTheory.f = figure('Position',[20 -20 1800 1800]),
sFP.act.Name = {'Stay','Left','Right','Up','Down','LeftUp','RightUp','LeftDown','RightDown'};
sFP.nCols = 3 + size(baseTasks,1) .* size(s.clc.actConsequence,1) ;
sFP.nRows = 5;
iCol = 1;
for iPol = 1:size(baseTasks,1)
    for iTask = 1:size(baseTasks,2)
        currQ = baseTasks{iPol,iTask}.allQ;
        for iAct = 1:size(currQ,5)
            subplot(sFP.nRows, sFP.nCols, sFP.nCols.*(iAct - 1) + iCol);
            sFP.plt.plAct = iAct;
            DisplActValsFun(sFP,w,currQ); hold on
            GridOverImage(fS,gca);
            colormap(flip(redbluecmapRory));
            caxis(fS.cAxis);
            title(['Q: ' sFP.act.Name{iAct}]);
            set(gca,'box','off');
            axis square
            title([ baseTasks{iPol,iTask}.name ]);
            if iCol == 1
                ylabel([ sFP.act.Name{iAct} ]);
            end
        end
        iCol = iCol + 1;
    end
end


sFP.plt.plAct = 1;
cAct = sFP.plt.plAct;

% Plot the Q values to be fitted
for iFig = 1:numel(tasksToRecreate)

    allQtoRecreate  = tasksToRecreate{iFig}.allQtoRecreate; %polGtaskTQ(:,:,:,:,1);
    allQtoPlot      = allQtoRecreate(:,:,:,:,cAct);

    subplot(sFP.nRows, sFP.nCols, sFP.nCols.*(iFig - 1) + size(baseTasks,1) .* size(s.clc.actConsequence,1) + 3 );
    DisplActValsFun(sFP,w,allQtoPlot); hold on
    GridOverImage(fS,gca);
    colormap(flip(redbluecmapRory));
    caxis(fS.cAxis);
    title(['Optimal Q for task: '  tasksToRecreate{iFig}.name])
    set(gca,'box','off')
    axis square
    
% % %     % ==================================================================
% % %     % FIT THE DATA BY JUST LINEARLY WEIGHTING ACTIONS (not identical to
% % %     % Barreto, but simple)
% % %     % Reshape polGtaskGQ for fitting
% % %     baseQtmp  = baseQ(:,:,:,:,:); % This version uses both tasks and individual actions as elements
% % %     baseQFlat = permute(baseQtmp, [5 1 2 3 4]);
% % %     baseQFlat = reshape(baseQFlat, size(baseQtmp,5), []);
% % % 
% % %     % Reshape allQtoRecreate
% % %     allQtoRecreateFlat = allQtoRecreate(:);
% % % 
% % %     corrFacts = (baseQFlat'\allQtoRecreateFlat);
% % %     fittedDataFlat = baseQFlat' * corrFacts;
% % % 
% % %     % Reshape the fitted data back to its original shape
% % %     fittedData = reshape(fittedDataFlat, 14, 15, 14, 15);
% % %     % ==================================================================

    % ==================================================================
    % FIT THE DATA MORE CLOSE TO THE BARRETO METHOD (i.e. don't use actions
    % as separate features)
    % Define function to be optimised
    FunToOpt = @(p) ErrFun(baseQ, allQtoRecreate, p, cAct, sFP);

    % Set optimisation options - keep it simple
    A = []; b = []; Aeq = []; beq = []; a = tic;
    w0 = ones([size(baseQ,7), 1]);
    lb = -[Inf Inf]';
    ub = [Inf Inf]';
    % Run optimisation
    OPTIONS = optimset('TolCon',1e-10);
    [p,sSqErr,exitflag,output,lambda,grad,hessian] = fmincon(FunToOpt,w0,A,b,Aeq,beq,lb,ub,[],OPTIONS);
    optTime = toc(a)

    % Extract optimised data
    [sSqErrFinal, bestPsi] = ErrFun(baseQ, allQtoRecreate, p, cAct, sFP);
    % fittedData =

    fittedData   = zeros(size(baseQ(:,:,:,:,1)));
    weightedData = sum(baseQ(:,:,:,:,cAct,:,:,:) .* permute(p,[2 3 4 5 6 7 1]),7);
    for r = 1:size(bestPsi,1)
        for c = 1:size(bestPsi,2)
            for rr = 1:size(bestPsi,1)
                for cc = 1:size(bestPsi,2)

                    % Select the correct policy for each condition
                    fittedData(r,c,rr,cc) = weightedData(r,c,rr,cc,bestPsi(r,c,rr,cc));

                end
            end
        end
    end
    % ==================================================================


    subplot(sFP.nRows, sFP.nCols, sFP.nCols.*(iFig - 1) + size(baseTasks,1) .* size(s.clc.actConsequence,1) + 2  );
    DisplActValsFun(sFP,w,fittedData); hold on
    GridOverImage(fS,gca);
    colormap(flip(redbluecmapRory));
    caxis(fS.cAxis)
    title(['Approximation :' sprintf('%.2f  ' , p)] );
    set(gca,'box','off')
    axis square

    % Plot the chosen policy
    subplot(sFP.nRows, sFP.nCols, sFP.nCols.*(iFig - 1) + size(baseTasks,1) .* size(s.clc.actConsequence,1) + 1 );
    DisplActValsFun(sFP,w,bestPsi); hold on
    GridOverImage(fS,gca);
%     caxis([1 2])
    title(['Optimal Psi'] );
    set(gca,'box','off')
    axis square

end

% % % % Show the weighted sum of the successor representation
% % % subplot(sFP.nRows, sFP.nCols, sFP.nCols.*(numel(tasksToRecreate)) + 6  );
% % % psiVals = max(abs(baseQtmp),[],5);
% % % sFP.plt.plAct = 1;
% % % DisplActValsFun(sFP,w,psiVals); hold on
% % % GridOverImage(fS,gca);
% % % colormap(flip(redbluecmapRory));
% % % caxis(fS.cAxis);
% % % title(['PSI']);
% % % set(gca,'box','off')
% % % axis square


% % % % Show 'data'
% % % subplot(sFP.nRows, sFP.nCols, sFP.nCols.*(numel(tasksToRecreate)) + 7  );
% % % rng('default');
% % % psiVals = max(abs(baseQtmp),[],5) + (rand(size(psiVals))-0.5)./4;
% % % sFP.plt.plAct = 1;
% % % DisplActValsFun(sFP,w,psiVals); hold on
% % % GridOverImage(fS,gca);
% % % colormap(flip(redbluecmapRory));
% % % caxis(fS.cAxis);
% % % title(['empirical data']);
% % % set(gca,'box','off')
% % % axis square

colormap(flip(redbluecmapRory));

% Go through the optimal policies and change colour maps locally
for iFig = 1:numel(tasksToRecreate)
            subplot(sFP.nRows, sFP.nCols, sFP.nCols.*(iFig - 1) + size(baseTasks,1) .* size(s.clc.actConsequence,1) + 1  );
%             colormap(gca,coltocol(100,[0.8 0.8 0.8],[0.2 0.2 0.2]));
%             caxis([.5 3])
            colormap(gca,[.3 .3 .8 ; .8 .3 .3; 0 0 .6])
end

%% Save figures
allFields = fields(f);
for iF = 1:length(allFields)
    cF = allFields{iF};
    set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['Results\ForFigures\TheoreticalFigMap\' cF 'BIGGER.eps'] , 'epsc')
    saveas(f.(cF).f,['Results\ForFigures\TheoreticalFigMap\' cF 'BIGGER.pdf'] , 'pdf')
end


%% FUNCTIONS

function plotPos = MakeSubplotPos(fillSize, subSize, subSpace)
plotPos=(1-fillSize)./2;
for iX=2:length(subSize)
    plotPos = [plotPos plotPos(end)+subSize(iX-1) + subSpace] ;
end
end


function [newQ optQ] = CalcQDirect(s)
newQ = zeros([s.wrld.size s.wrld.size]);
optQ = zeros([s.wrld.size s.wrld.size]);

% first initialise with the contact reward
for iLC = 1:length(newQ)
    newQ(:,iLC,s.clc.startSC,iLC) = s.clc.startRew;
    optQ(:,iLC,s.clc.startSC,iLC) = s.clc.startRew;
end

% next, loop through rows backwards and add sums
for iLC = 3:length(newQ)-2
    for iSR = s.clc.startSC - 1 : -1 : 1

        for iSC = 3:size(newQ,4) - 2
            colsBelow = iSC + s.clc.randSpread;

            % For each action, take the maximum expected value that it could lead
            % to (i.e. accounting for random spread after action is taken)
            newQ(:,iLC,iSR,iSC) = s.clc.gammaVal .* ...
                sum(s.clc.spreadProb .* ...
                squeeze(optQ(:,iLC,iSR + 1,colsBelow )),2);

            for iAct = 1:length(s.clc.actConsequence)
                optQ(:,iLC,iSR,iSC) = max(optQ(:,iLC,iSR,iSC), s.clc.gammaVal .* ...
                    sum(s.clc.spreadProb .* ...
                    squeeze(optQ(:,iLC,iSR + 1,colsBelow + s.clc.actConsequence(iAct) )),2));
            end

        end

    end
end

% Then set the value of the touch condition back to 0
for iLC = 1:length(newQ)
    newQ(:,iLC,s.clc.startSC,iLC) = 0;
end

end


% Calculate error w.r.t. optimal data
function [sSqErr bestPsi] = ErrFun(baseQ, allQtoRecreate, w, cAct, s)

% if sFP.clc.useActsAsBases
% % % allSSqErr = (sum(baseQ(:,:,:,:,cAct,:,:,:) .* permute(w,[2 3 4 5 6 7 1]),7) - allQtoRecreate) .^ 2;

newQ = sum(baseQ .* permute(w,[2 3 4 5 6 7 1]),7);

switch s.clc.maximiseSimilarityType
    case 'OverallQ'
        % Optimise weights over all actions, so sum over the actions and the tasks
        allSSqErr = sum( (newQ - allQtoRecreate) .^ 2 , 5);
    case 'ChosenAction'
        [~, chosenAct]           = max(newQ, [] ,5);
        [~, chosenActToRecreate] = max(allQtoRecreate , [], 5);

        allSSqErr = sum( ((chosenAct - chosenActToRecreate) ~= 0) , 5);

    case 'WinningQ' 
        [maxQs, chosenAct]   = max( newQ, [] ,5);
        [maxQToRecreate, chosenActToRecreate] = max(allQtoRecreate  , [], 5);
        
        % Compare to most valuable alternative
        winningQ            = zeros(size(maxQs));
        winningQToRecreate  = zeros(size(maxQToRecreate));
        for iR = 1:size(newQ,1)
            for iC = 1:size(newQ,2)
                for iRR = 1:size(newQ,3)
                    for iCC  = 1:size(newQ,4)
                        for iPol = 1:size(newQ,6)
                            otherActs = (1:size(newQ,5) ~= chosenAct(iR,iC,iRR,iCC,:,iPol));
                            winningQ(iR, iC, iRR, iCC,:,iPol)  = ...
                                maxQs(iR, iC, iRR, iCC,:,iPol) - ...
                                max(newQ(iR, iC, iRR, iCC,otherActs,iPol), [], 5);


                            otherActsToRecreate = (1:size(newQ,5) ~= chosenActToRecreate(iR,iC,iRR,iCC));
                            winningQToRecreate(iR, iC, iRR, iCC,:) = ...
                                maxQToRecreate(iR, iC, iRR, iCC,:) - ...
                                max(allQtoRecreate(iR, iC, iRR, iCC,otherActsToRecreate), [], 5);
                        end
                    end
                end
            end
        end
        
     allSSqErr = sum( (winningQ - winningQToRecreate).^2 , 5);
end

[minErrs, bestPsi] =  min( allSSqErr , [] ,6);

sSqErr = sum( minErrs , 'all');

% % % [bestPsi{1}, bestPsi{2}, bestPsi{3}, bestPsi{4}, bestPsi{5}] = ind2sub([size(minErrs) , size(allSSqErr , 6) ], minLinInd);

end

