
% Initial Novelty value. $$$ ASSUME THIS RESETS EVERY TIME THE HAND IS MOVED?
vNovInit = 2 .*[-1 -1]';
% Initial Position value
vPosInit = [-.8 -0]';


% Rewards for 'novelty', i.e. having experienced a stimulus before recently
rNov           = [ 0  0]';
% rewards for each position
rPos           = [-2 -0]';
rPosDuringTest = [ 0  0]';

% Learning rate
alphaPos = 0.2;
alphaNov = 0.5;

% Discount factor (make zero for 'end of simulation is now' situation
% (Note I'm assuming a lack of movement here, I guess)
gamma = 0;


% Learning timepoints
learningSteps = 1:14; % 12 + 2
startLearnStep = 3;
testingSteps  = 1:6;

[vPosDuringLearn1, vNovDuringLearn1, ...
          vPosDuringTest1,  vNovDuringTest1,  ...
          vPosDuringLearn2, vNovDuringLearn2, ...
          vPosDuringTest2,  vNovDuringTest2 ] = ...
              CalculateZanini2021SCR( vPosInit, alphaPos, gamma, rPos, ...
                                      vNovInit, alphaNov, rNov);


figure('Position', [ 50, 50, 400, 900]),
subplot(4,1,1)
plot(1:11 , abs(vPosDuringLearn1 + vNovDuringLearn1)','-o');
title('Learning phase');
legend('near','far')

subplot(4,1,2)
plot(testingSteps , abs(vPosDuringTest1 + vNovDuringTest1)','-o');
title('Testing phase1');
legend('near','far')

subplot(4,1,3)
plot(1:11 , abs(vPosDuringLearn2 + vNovDuringLearn2)','-o');
title('Learning phase2');
legend('near','far')

subplot(4,1,4)
plot(testingSteps , abs(vPosDuringTest2 + vNovDuringTest2)','-o');
title('Testing phase2');
legend('near','far')

%% $$$ Try to find the optimal parameters


avOverPhaseFl = 1;


realDat = {[1.029 0.814 0.732 0.671 0.719 0.741 0.662 0.645 0.577 0.633 0.665 ; ... Learning 1, Near
            0.745 0.612 0.478 0.498 0.382 0.342 0.458 0.448 0.472 0.404 0.394], ... Learning 1, Far
           [1.214 1.055 0.809 0.632 0.689 0.795 ; ... Testing 1, Near
            0.518 0.415 0.287 0.371 0.318 0.386], ... Testing 1, Far
           [0.779 0.560 0.769 0.699 0.679 0.642 0.633 0.751 0.672 0.799 0.552 ; ... Learning 2, Near
            0.602 0.485 0.450 0.483 0.452 0.378 0.382 0.459 0.393 0.366 0.582], ... Learning 2, Far
           [1.239 0.899 0.671 0.619 0.492 0.400 ; ... Testing 2, Near
            0.649 0.345 0.410 0.300 0.382 0.354]};... Testing 2, Far

% Highest estimate
realStd = {[1.090 0.868 0.790 0.723 0.780 0.828 0.713 0.694 0.643 0.680 0.755 ; ... Learning 1, Near
            0.809 0.667 0.530 0.556 0.450 0.394 0.492 0.485 0.523 0.457 0.450], ... Learning 1, Far
           [1.325 1.132 0.924 0.727 0.795 0.926 ; ... Testing 1, Near
            0.594 0.470 0.329 0.448 0.380 0.443], ... Testing 1, Far
           [0.860 0.616 0.836 0.785 0.752 0.716 0.711 0.834 0.739 0.921 0.588 ; ... Learning 2, Near
            0.675 0.551 0.504 0.529 0.498 0.434 0.441 0.512 0.456 0.421 0.652], ... Learning 2, Far
           [1.367 1.001 0.760 0.731 0.573 0.487 ; ... Testing 2, Near
            0.729 0.415 0.500 0.368 0.434 0.411]};... Testing 2, Far
% Actual standard deviation
for iD = 1 : numel(realDat)
    realDat{iD} = realDat{iD} .* 0.3660 ;
    realStd{iD} = realStd{iD} .* 0.3660 ; % .* sqrt(41)
end
% Correct the values to the right units
for iD = 1 : numel(realDat)
    realStd{iD} = realStd{iD} - realDat{iD};
end

% Take the average from the two phases
if avOverPhaseFl == 1
    for iD = 1 : numel(realDat) / 2
        realDat{iD} =     (realDat{iD}    + realDat{iD+2})./2;
%         realStd{iD} = sqrt(realStd{iD}.^2 + realStd{iD+2}.^2 + ((realDat{iD} + realDat{iD+2}).^2)./4 ); % .* sqrt(41)
        % More conservative std estimatesion
        realStd{iD} = sqrt(realStd{iD}.^2 + realStd{iD+2}.^2 ); % .* sqrt(41)
    end
    realDat = realDat(1:2);
    realStd = realStd(1:2);
end

p0 = [vPosInit(:); rPos(:); alphaPos; gamma; ... [1 2] [3 4] [5] [6]
      vNovInit(:); rNov(:); alphaNov(:); ... [7 8] [9 10] [11]
      0; 1]; % offset, mult [12 13]
lb = [[-10 -10]' ; [-10 -10]' ; 0 ; 0; ...
      [-10 -10]' ; [-10 -10]' ; 0 ; ...
      -10; -10];
ub = [[   0    0]' ; [   0    0]' ; 1 ; 1; ...
      [   0    0]' ; [   0    0]' ; 1 ; ...
      100; 100];



% % % FunToOpt = @(params) CalcErrSCR(realDat,[ [-1 -.5]';  p0]);

% % % % This is just optimising everything
% % % FunToOpt = @(params) CalcErrSCR(realDat,params);

% Select jparameters for fitting
fitParamIDs = [3 4 5 7 11 12 13];
% Repeating params 7 means that the novelty score is the same for both position
FunToOpt = @(params) CalcErrSCR(realDat,[p0(1:2); params(1:3); p0(6) ; params([4 4]); p0(9:10); params(5:7)]);
% % % fitParamIDs = [1 2 3 4 5 7 11 12 13];
% % % % Repeating params 7 means that the novelty score is the same for both position
% % % FunToOpt = @(params) CalcErrSCR(realDat,[params(1:5); p0(6) ; params([6 6]); p0(9:10); params(7:9)]);

% Set optimisation options - keep it simple
A = []; b = []; Aeq = []; beq = []; a = tic;
% Run optimisation
tic
OPTIONS = optimset('TolCon',1e-10);
[params,sSqErr,exitflag,output,lambda,grad,hessian] = fmincon(FunToOpt,p0(fitParamIDs),A,b,Aeq,beq,lb(fitParamIDs),ub(fitParamIDs),[],OPTIONS)
toc

% % % [sSqErr, SqErr, modelDat] = CalcErrSCR(realDat,p);
[sSqErr, SqErr, modelDat] = CalcErrSCR(realDat,[p0(1:2); params(1:3); p0(6) ; params([4 4]); p0(9:10); params(5:7)]);
% % % [sSqErr, SqErr, modelDat] = CalcErrSCR(realDat,[params(1:5); p0(6) ; params([6 6]); p0(9:10); params(7:9)]);

fitRes = CalcFitMetrics(realDat,realStd,modelDat,params)

params

% plot the data
allTitles = {'Learning phase 1','Testing phase1','Learning phase2','Testing phase2'};
colours = [0.8 0.3 0.3 ; 0.3 0.3 0.8];
figure('Position', [ 50, 50, 400, 900]),
for iD = 1:numel(realDat)

    subplot(numel(realDat),1,iD)
%     plot(realDat{iD}','-o');
    errorbar(1:size(realDat{iD},2) , realDat{iD}', realStd{iD}', 'o');
    hold on;
    plot(modelDat{iD}','-','LineWidth',2);
    title(allTitles{iD});
    legend('near','far')

    colororder(colours)
end


%% $$$ NOW it's time to do the reachability paper
load('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Data\allQ.mat')
load('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Data\sHnd.mat')
%%

% Find the values of hand- actions at the 'height of the table', using the 
% Qs averaged from approaching and retreating
% % % qOnTable = squeeze((allQ(13,:).qVals{1} + allQ(15,:).qVals{1}) ./ 2);
% % % qOnTable = squeeze(qOnTable(:,:,allQ(13,:).centerpos{1}(2),:));

qOnTable = allQ2(1,:).qVals{1}(:,:,:,2);

% Define the positions of the real stimuli
% Assume hand is just at the base of the central block
blockXsize = 88.56/7;
blockYsize = 49.81/6;
stimPosX   = blockXsize .* -4   + blockXsize .* (1:7);
stimPosY   = blockYsize .* -0.5 + blockYsize .* (1:6); 
[sX sY]    = meshgrid(stimPosX,stimPosY);


onTableQPosY = 62 - (1:size(qOnTable,2));
onTableQPosY = 5.* (  onTableQPosY - ((size(qOnTable,2)+1) - allQ2(1,:).centerpos{1}(1))  );
onTableQPosX = 5.* (  (1:size(qOnTable,3)) - allQ2(1,:).centerpos{1}(2)  );
[qX qY]      = meshgrid(onTableQPosX,onTableQPosY);

% % % figure,imagesc(qX(:),qY(:),sqrt(qX.^2 + qY.^2));
figure,imagesc(qX(:),qY(:),squeeze(max(qOnTable)) );
% figure,imagesc(qX(:),qY(:),squeeze(qOnTable(19,:,:)) ); % Hand forward
hold on
scatter(sX,sY,20,[0 0 0],'filled');
axis xy

% Find action values at positions in question
clear stimPosQ
for iSX = 1:numel(stimPosX)
    xBin = findnearest(stimPosX(iSX), onTableQPosX );
    for iSY = 1:numel(stimPosY)
        yBin = findnearest(stimPosY(iSY), onTableQPosY );
        stimPosQ(:,iSY,iSX) = qOnTable(:,yBin,xBin);
    end    
end

figure,imagesc(squeeze(max(stimPosQ)))
axis xy



%% Calculate probability of performing reaching actions

sHnd2.clc.baseVel = [0 0 0];

sHnd2.rch.temp = 0.07;
sHnd2.rch.nTP  =  20;
% =========================================================================
% CASE: 50-50 chance of reward everywhere (just a flat line basically)
sHnd2.rch.farVal = 1.0;
sHnd2.rch.nearVal = 1.0;


sHnd2.rch.tableXMin = -88.56/2 ;
sHnd2.rch.tableXMax =  88.56/2 ;
sHnd2.rch.tableYMin =   8      ; % Width of half of hand
sHnd2.rch.tableYMax =  49.81   ;

[handLocProbs, transProbs, onTableQPosX, onTableQPosY] = CalcReachProbabilities(sHnd2,allQ2);

figure,imagesc(onTableQPosX, onTableQPosY, handLocProbs(:,:,sHnd2.rch.nTP)); axis xy
xlim([-50 50]);
ylim([0 100]);
SquareAxes;

% =========================================================================
% CASE: 75% chance of reward far, 25% of reward near
sHnd2.rch.farVal  = 1.50;
sHnd2.rch.nearVal = 0.75;

[handLocProbs_FarRew, transProbs_FarRew] = CalcReachProbabilities(sHnd2,allQ2);

figure,imagesc(onTableQPosX, onTableQPosY, handLocProbs_FarRew(:,:,sHnd2.rch.nTP)); axis xy
xlim([-50 50]);
ylim([0 100]);
SquareAxes;


% =========================================================================
% CASE: 25% chance of reward far, 75% of reward near
sHnd2.rch.farVal  = 0.75;
sHnd2.rch.nearVal = 1.50;

[handLocProbs_NearRew, transProbs_NearRew] = CalcReachProbabilities(sHnd2,allQ2);

figure,imagesc(onTableQPosX, onTableQPosY, handLocProbs_NearRew(:,:,sHnd2.rch.nTP)); axis xy
xlim([-50 50]);
ylim([0 100]);
SquareAxes;



%% Calculate reaching action proababilities for a variety of different learning rates, temperatures, and blocks

s2.nSteps           = 39;
s2.initV            = [1 1]';
s2.alpha            = 0.4;
s2.rew              = [1.5 0.75]';
s2.gamma            = 0;
valsForReach        = applyTDLearn(s2); % Familiarisation

% Make a plot to test what kind of alpha we need
figure,plot([1, 1 ; valsForReach'])

allAlphas = 0.1:0.1:0.9;

allTemps = [ 0.05 0.07 0.1 0.14 0.2 0.3];

% $$$ GOT HERE: iTemp 1, iAlpha 7, iLEarnTP 7

% $$$ In the meantime, for model fitting and stuff, just set
% $$$ iTemp 2, iAlpha 4, iLearnTP 1:40
for iTemp = 3:numel(allTemps)

    sHnd2.rch.temp = allTemps(iTemp);

    for iAlpha = 4 %1:numel(allAlphas)

        % Set learning rate Alpha
        s2.alpha = allAlphas(iAlpha);

        % Update values of states as a function of time spent leanring
        valsForReach       = applyTDLearn(s2);

        % Insert the 'no learning yet' stage
        valsForReach       = [s2.initV, valsForReach];

        for iLearnTP = 1:40

            sHnd2.rch.nTP  = 20;
            % =========================================================================
            % CASE: 50-50 chance of reward everywhere (just a flat line basically)
            sHnd2.rch.farVal = 1.0;
            sHnd2.rch.nearVal = 1.0;

            sHnd2.rch.tableXMin = -88.56/2 ;
            sHnd2.rch.tableXMax =  88.56/2 ;
            sHnd2.rch.tableYMin =   8      ; % Width of half of hand
            sHnd2.rch.tableYMax =  49.81   ;

            if iAlpha == 4 & iLearnTP == 1
                % This only needs to happen once because learning changes ABOLUTELY NOTHING
                [handLocProbs_Learn(:,:,:,iLearnTP,iAlpha,iTemp), ~, onTableQPosX, onTableQPosY] = CalcReachProbabilities(sHnd2,allQ2);
            else
                handLocProbs_Learn(:,:,:,iLearnTP,iAlpha,iTemp) = handLocProbs_Learn(:,:,:,1,1,iTemp);
            end

            % =========================================================================
            % CASE: 75% chance of reward far, 25% of reward near
            sHnd2.rch.farVal  = valsForReach(1,iLearnTP);
            sHnd2.rch.nearVal = valsForReach(2,iLearnTP);

            [handLocProbs_Learn_FarRew(:,:,:,iLearnTP,iAlpha,iTemp)] = CalcReachProbabilities(sHnd2,allQ2);


            % =========================================================================
            % CASE: 25% chance of reward far, 75% of reward near
            sHnd2.rch.farVal  = valsForReach(2,iLearnTP);
            sHnd2.rch.nearVal = valsForReach(1,iLearnTP);

            [handLocProbs_Learn_NearRew(:,:,:,iLearnTP,iAlpha,iTemp)] = CalcReachProbabilities(sHnd2,allQ2);

            disp(['alpha: '   num2str(iAlpha)]);
            disp(['learnTP: ' num2str(iLearnTP)]);
        end
    end
    save('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Data\Fitting_Learning_Workspace_v2.mat','-v7.3')
end


%% Reproduce the findings from the paper
%


% $$$ NEXT I need to calculate the sumsquare error and then optimise for
% sum square error using different values of cTemp. 

% $$$ NOTE ALSO: maybe I can run some more values of cTemp while I'm having
% dinner

sHnd2.rch.plotFl     =    0;

sHnd2.rch.cDiffuseTP =   20;
sHnd2.rch.cAlpha     =    4;
sHnd2.rch.cTemp      =    2;
sHnd2.rch.cLearnTP   = 1:40; 
sHnd2.rch.gaussKern  =    1; 
sHnd2.rch.gaussSigm  =    1;


sHnd2.rch.nFreeParams = 0;

[fitRes, normSSqErr, realDat, realDatStd, modelDat] = ReproduceReachingExperiment(sHnd2,...
    handLocProbs_Learn,handLocProbs_Learn_FarRew,handLocProbs_Learn_NearRew,...
    qX,qY,onTableQPosX,onTableQPosY);

fitRes



%% $$$ OK, I think next is just to add a little bit of RANDOM NOISE IN, then
% everything should fit? OR do I have to go through the whole fitting
% shebang with all the different temperature values, while NOT adding in
% additional noise? Hmm Hmm indeed



%% Functions


function [fitRes, errSq, realDat, realDatStd, modelDat] = ReproduceReachingExperiment(s,handLocProbs_Learn,handLocProbs_Learn_FarRew,handLocProbs_Learn_NearRew,...
    qX,qY,onTableQPosX,onTableQPosY)

% =========================================================================
% Make gaussian blurred version

handLocProbs_Learn_G         = GaussBlur(handLocProbs_Learn,        s.rch.gaussKern,s.rch.gaussSigm );
handLocProbs_Learn_FarRew_G  = GaussBlur(handLocProbs_Learn_FarRew, s.rch.gaussKern,s.rch.gaussSigm );
handLocProbs_Learn_NearRew_G = GaussBlur(handLocProbs_Learn_NearRew,s.rch.gaussKern,s.rch.gaussSigm );


% =========================================================================
% Plot and reproduce 2D reaching results

% Control
pC = squeeze(mean(handLocProbs_Learn_G(:,:,s.rch.cDiffuseTP,s.rch.cLearnTP,s.rch.cAlpha,s.rch.cTemp), 4));
% Far
pF = squeeze(mean(handLocProbs_Learn_FarRew_G(:,:,s.rch.cDiffuseTP,s.rch.cLearnTP,s.rch.cAlpha,s.rch.cTemp), 4));
% Near
pN = squeeze(mean(handLocProbs_Learn_NearRew_G(:,:,s.rch.cDiffuseTP,s.rch.cLearnTP,s.rch.cAlpha,s.rch.cTemp), 4));

blockXsize = 88.56/7;
blockYsize = 49.81/6;
stimPosX   = blockXsize .* -4   + blockXsize .* (1:7);
stimPosY   = 2.5 + blockYsize .* -0.5 + blockYsize .* (1:6); 

% Find action values at positions in question
clear endProbs
for iSX = 1:numel(stimPosX)
    xBinMin = stimPosX(iSX) - blockXsize./2 ;
    xBinMax = stimPosX(iSX) + blockXsize./2 ;
    for iSY = 1:numel(stimPosY)
        if iSY == 4 & iSX == 6
            disp('test')
        end
        yBinMin = stimPosY(iSY) - blockYsize./2 ;
        yBinMax = stimPosY(iSY) + blockYsize./2 ;
        inThisBin = qY >= yBinMin & qY <= yBinMax & ...
                    qX >= xBinMin & qX <= xBinMax;
        [rr cc] = find(inThisBin);

        sumProbsC = 0;
        sumProbsF = 0;
        sumProbsN = 0;
        for iSmallerBin = 1:numel(rr)
             sumProbsC = sumProbsC + sum(pC(rr(iSmallerBin), cc(iSmallerBin)));
             sumProbsF = sumProbsF + sum(pF(rr(iSmallerBin), cc(iSmallerBin)));
             sumProbsN = sumProbsN + sum(pN(rr(iSmallerBin), cc(iSmallerBin)));
        end
        % calculate the probabilities of ending up in a given block:
        % Multiply by the surface ratio between the experiment and modelled
        % blocks, and divide by the number of samples that fell within the
        % experimental block
        pResC(iSY,iSX) = ((blockXsize .* blockYsize) ./ 25) .* (sumProbsC ./ numel(rr));
        pResF(iSY,iSX) = ((blockXsize .* blockYsize) ./ 25) .* (sumProbsF ./ numel(rr));
        pResN(iSY,iSX) = ((blockXsize .* blockYsize) ./ 25) .* (sumProbsN ./ numel(rr));
    end
end

if s.rch.plotFl == 1
    figure('Position' , [100, 100, 900, 400]);
    subplot(1,3,1)
    imagesc(stimPosX,stimPosY,pResF)
    title('Far group')
    axis xy
    colorbar
    clim([0 0.07])

    subplot(1,3,2)
    imagesc(stimPosX,stimPosY,pResC)
    title('Control group')
    axis xy
    colorbar
    clim([0 0.07])

    subplot(1,3,3)
    imagesc(stimPosX,stimPosY,pResN)
    title('Near group')
    axis xy
    colorbar
    clim([0 0.07])

    colormap jet
end

% =========================================================================
% Plot and reproduce 1D reaching distances results

earlyBlocks =  1:3;
lateBlocks  = 31:33;
% Control
pEarlyC = squeeze(handLocProbs_Learn_G(:,:,s.rch.cDiffuseTP,earlyBlocks,s.rch.cAlpha,s.rch.cTemp));
pLateC  = squeeze(handLocProbs_Learn_G(:,:,s.rch.cDiffuseTP,lateBlocks, s.rch.cAlpha,s.rch.cTemp));
% Far
pEarlyF = squeeze(handLocProbs_Learn_FarRew_G(:,:,s.rch.cDiffuseTP,earlyBlocks,s.rch.cAlpha,s.rch.cTemp));
pLateF  = squeeze(handLocProbs_Learn_FarRew_G(:,:,s.rch.cDiffuseTP,lateBlocks, s.rch.cAlpha,s.rch.cTemp));
% Near
pEarlyN = squeeze(handLocProbs_Learn_NearRew_G(:,:,s.rch.cDiffuseTP,earlyBlocks,s.rch.cAlpha,s.rch.cTemp));
pLateN  = squeeze(handLocProbs_Learn_NearRew_G(:,:,s.rch.cDiffuseTP,lateBlocks, s.rch.cAlpha,s.rch.cTemp));


realDat = {[0.619 0.535 0.481 ; ... Reaching distance during 1st blocks
           0.698 0.528 0.420]};... Reaching distance during last blocks

% Highest estimate
realDatStd = {[0.740 0.632 0.586 ; ... Reaching distance during 1st blocks
               0.783 0.620 0.520]};... Reaching distance during last blocks
% Actual proper units
for iD = 1 : size(realDat{1},1)
    realDat{1}(iD,:)    = realDat{1}(iD,:)    .* 45.50 ;
    realDatStd{1}(iD,:) = realDatStd{1}(iD,:) .* 45.50 ; 
end
% Correct the values to the right units
for iD = 1 : size(realDat{1},1)
    realDatStd{1}(iD,:) = realDatStd{1}(iD,:) - realDat{1}(iD,:);
end

% Calculate the mean reach distances and standard deviations
[meanEarlyReachC, stdEarlyReachC] = CalcMeanReachDist(qY, pEarlyC);
[meanLateReachC,  stdLateReachC ] = CalcMeanReachDist(qY, pLateC);
[meanEarlyReachF, stdEarlyReachF] = CalcMeanReachDist(qY, pEarlyF);
[meanLateReachF,  stdLateReachF ] = CalcMeanReachDist(qY, pLateF);
[meanEarlyReachN, stdEarlyReachN] = CalcMeanReachDist(qY, pEarlyN);
[meanLateReachN,  stdLateReachN ] = CalcMeanReachDist(qY, pLateN);

modelDat{1}    = [meanEarlyReachF meanEarlyReachC meanEarlyReachN ; ...
                  meanLateReachF  meanLateReachC  meanLateReachN ];
modelDatStd{1} = [stdEarlyReachF  stdEarlyReachC  stdEarlyReachN ; ...
                  stdLateReachF   stdLateReachC   stdLateReachN ];

if s.rch.plotFl == 1
    figure,
    subplot(1,2,1);
    errorbar([0.85 1.15 ; 1.85 2.15 ; 2.85 3.15] , realDat{1}', realDatStd{1}','.k');
    hold on;
    bar(realDat{1}');
    ylim([0 40])

    subplot(1,2,2);
    errorbar([0.85 1.15 ; 1.85 2.15 ; 2.85 3.15], modelDat{1}', modelDatStd{1}','.k');
    hold on;
    bar(modelDat{1}');
    ylim([0 40])
end


% =========================================================================
% Plot and reproduce the 'proximal selection' as a function of block

distPositions  = onTableQPosY >= ( blockYsize .* 4 + 2.5);
proxPositions  = onTableQPosY <  ( blockYsize .* 2 + 2.5);

% Control
pProxC = squeeze(sum(handLocProbs_Learn_G(proxPositions,:,s.rch.cDiffuseTP,:,s.rch.cAlpha,s.rch.cTemp),[1 2]));
pDistC = squeeze(sum(handLocProbs_Learn_G(distPositions,:,s.rch.cDiffuseTP,:,s.rch.cAlpha,s.rch.cTemp),[1 2]));
% Far Rewards
pProxF = squeeze(sum(handLocProbs_Learn_FarRew_G(proxPositions,:,s.rch.cDiffuseTP,:,s.rch.cAlpha,s.rch.cTemp),[1 2]));
pDistF = squeeze(sum(handLocProbs_Learn_FarRew_G(distPositions,:,s.rch.cDiffuseTP,:,s.rch.cAlpha,s.rch.cTemp),[1 2]));
% Near Rewards
pProxN = squeeze(sum(handLocProbs_Learn_NearRew_G(proxPositions,:,s.rch.cDiffuseTP,:,s.rch.cAlpha,s.rch.cTemp),[1 2]));
pDistN = squeeze(sum(handLocProbs_Learn_NearRew_G(distPositions,:,s.rch.cDiffuseTP,:,s.rch.cAlpha,s.rch.cTemp),[1 2]));

% Store real data
pProxRealN = [0.555 0.520 0.537 0.521 0.520 0.550 0.546 0.631 0.645 0.600 0.566 0.575 0.606 0.604 0.573 0.556 0.586 0.559 0.615 0.593 0.601 0.642 0.615 0.674 0.546 0.584 0.549 0.574 0.616 0.636 0.678 0.587 0.600 0.603 0.532 0.575 0.588 0.612 0.665 0.575] .* 124.5330012 ./ 100;
pProxRealC = [0.492 0.459 0.461 0.418 0.441 0.403 0.468 0.385 0.378 0.488 0.384 0.435 0.379 0.449 0.389 0.445 0.532 0.448 0.466 0.475 0.489 0.509 0.467 0.468 0.379 0.425 0.508 0.443 0.475 0.396 0.467 0.457 0.459 0.389 0.444 0.449 0.462 0.500 0.530 0.413] .* 124.5330012 ./ 100;
pProxRealF = [0.445 0.363 0.287 0.257 0.302 0.282 0.309 0.280 0.246 0.279 0.202 0.248 0.309 0.221 0.208 0.207 0.218 0.234 0.257 0.309 0.239 0.246 0.308 0.231 0.209 0.229 0.218 0.253 0.195 0.215 0.270 0.237 0.239 0.220 0.224 0.246 0.220 0.209 0.274 0.178] .* 124.5330012 ./ 100;

% Note: shift x- axis by half a block, so that we're simulating the average
% learnign that happened during the block itself, rather than the result at
% the end of the block
if s.rch.plotFl == 1
    blockXAx = s.rch.cLearnTP;
    figure,
    hold off
    plot(blockXAx - 0.5, pProxN,'b');
    hold on
    plot(blockXAx, pProxRealN ,'.b','MarkerSize',20);
    plot(blockXAx - 0.5, pProxF, 'r');
    plot(blockXAx, pProxRealF ,'.r','MarkerSize',20);
    plot(blockXAx - 0.5, pProxC, 'k');
    plot(blockXAx, pProxRealC ,'.k','MarkerSize',20);
    ylim([0 1]);
end


% =========================================================================
% Plot and reproduce the reachability thresholds

% First, use relations from the paper to estimate modelled reachability
% ------------------------------------------------------------------
% Relationship with distal targets selected

% Coordinate data: figure size
sl.xMin = 208.485;
sl.xMax = 209.496;
sl.yMin =   0.634;
sl.yMax =   1.571;
% Real data: figure size
sl.xRMin = -10;
sl.xRMax =  15;
sl.yRMin =   0;
sl.yRMax = 350;
% Coordinate data: slope points
sl.x0    = 208.888;
sl.y0    =   0.634;
sl.xSMin = 208.485;
sl.xSMax = 209.360;
sl.ySMin =   0.691;
sl.ySMax =   1.222;
% Calculate the relationship
[rDistSlope rDistInt] = FindLinEqFromFig(sl);

% ------------------------------------------------------------------
% Relationship with proximal targets selected

% Coordinate data: figure size
sl.xMin = 208.491;
sl.xMax = 209.511;
sl.yMin =  -0.439;
sl.yMax =   0.514;
% Real data: figure size
sl.xRMin = -10;
sl.xRMax =  15;
sl.yRMin =   0;
sl.yRMax = 350;
% Coordinate data: slope points
sl.x0    = 208.899;
sl.y0    =  -0.439;
sl.xSMin = 208.485;
sl.xSMax = 209.378;
sl.ySMin =   0.127;
sl.ySMax =  -0.439;
% Calculate the relationship
[rProxSlope rProxInt] = FindLinEqFromFig(sl);

nDistF = mean(pDistF .* 400);
nProxF = mean(pProxF .* 400);
nDistC = mean(pDistC .* 400);
nProxC = mean(pProxC .* 400);
nDistN = mean(pDistN .* 400);
nProxN = mean(pProxN .* 400);

% Now find the predicted 'boundary' changes, and average them
bDistF = (nDistF - rDistInt) ./ rDistSlope;
bProxF = (nProxF - rProxInt) ./ rProxSlope;
bF     = (bProxF + bDistF)   ./ 2

bDistC = (nDistC - rDistInt) ./ rDistSlope;
bProxC = (nProxC - rProxInt) ./ rProxSlope;
bC     = (bProxC + bDistF) ./ 2

bDistN = (nDistN - rDistInt) ./ rDistSlope;
bProxN = (nProxN - rProxInt) ./ rProxSlope;
bN     = (bProxN + bDistF)   ./ 2

realDat{2}    = [0.305 -0.051 -0.302] .* 7.886;
realDatStd{2} = [0.435 (- 0.051 + (0.139 - 0.051))  (- 0.302 + (0.421 - 0.302)) ] .* 7.886;
realDatStd{2} = realDatStd{2} - realDat{2};

modelDat{2}    = [bF bC bN];

% =======================================================================
% Normalise and calculate error squared 
% (use the same scale to normalise the modelled data as the real data, ofc)
realM  = cellfun(@(cDat) mean(cDat(:)) , realDat , 'UniformOutput' , false);
realSD = cellfun(@(cDat) std(cDat(:))  , realDat , 'UniformOutput' , false);

realDatForErrorCalc    = arrayfun(@(iD) (realDat{iD} - realM{iD}) ./ realSD{iD} , 1:numel(realDat) , 'UniformOutput' , false);
realDatStdForErrorCalc = arrayfun(@(iD) (realDatStd{iD}         ) ./ realSD{iD} , 1:numel(realDat) , 'UniformOutput' , false);

modelDatForErrorCalc   = arrayfun(@(iD) (modelDat{iD} - realM{iD}) ./ realSD{iD} , 1:numel(modelDat) , 'UniformOutput' , false);

% Squared error
errSq = cell2mat( ...
    arrayfun(@(iCond) (modelDatForErrorCalc{iCond}(:) - realDatForErrorCalc{iCond}(:))' .^ 2 , ...
    (1 : numel(realDat)) , 'UniformOutput' , false)  );

% Sum Square error
errSq = sqrt(sum(errSq(:)));

dataVar = cell2mat(...
    arrayfun(@(iCond) realDatStdForErrorCalc{iCond}(:)' .^2 , ...
    (1 : numel(realDat))  , 'UniformOutput' , false));
chiSq   = nansum(errSq(:)./dataVar(:));


% Calculate the chi square confidence interval using bootstrapping
chiSqDistr     = arrayfun(@(x) nansum(randsample(errSq(:)./dataVar(:),numel(errSq(:)),1)) ,[1:100000] );
fitRes.chiSqCI = prctile(chiSqDistr,[2.5 97.5])

N   = numel(errSq(:));
dof = numel(s.rch.nFreeParams);
k   = N - dof;

% Calculate goodness of fit
[ fitRes.gofScore chiSqCheck fitRes.pVal1 fitRes.pVal2 ] = GoFFun( errSq, dataVar, k );


end



%%
function [normalized_blurred_matrix] = GaussBlur(full_prob_matrix, kernel_size, sigma)

% Make it work for >2-dim arrays
normalized_blurred_matrix = full_prob_matrix;

full_prob_matrix = full_prob_matrix(:,:,:);

for iFurtherDim = 1:size(full_prob_matrix,3)
    prob_matrix = full_prob_matrix(:,:,iFurtherDim);
    % Create a Gaussian kernel
    [x, y] = meshgrid(-floor(kernel_size/2):floor(kernel_size/2), -floor(kernel_size/2):floor(kernel_size/2));
    gaussian_kernel = exp(-(x.^2 + y.^2) / (2 * sigma^2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:)); % Normalize the kernel

    % Store the total probability of the original matrix
    total_prob = sum(prob_matrix(:));

    % Perform the convolution (2D Gaussian blur)
    blurred_matrix = conv2(prob_matrix, gaussian_kernel, 'same');

    % Normalize the blurred matrix to maintain the total probability
    normalized_blurred_matrix(:,:,iFurtherDim) = (blurred_matrix / sum(blurred_matrix(:))) * total_prob;
end

end



function [rSlope rInt] = FindLinEqFromFig(sl)

% Scaling to real units from figure size units
multFactY = (sl.yRMax - sl.yRMin) ./ (sl.yMax - sl.yMin);
multFactX = (sl.xRMax - sl.xRMin) ./ (sl.xMax - sl.xMin);

% Find slope
rSlope = ((sl.ySMax - sl.ySMin) .* multFactY )  ./ ((sl.xSMax - sl.xSMin) .* multFactX);
% Find intercept
rInt = ((sl.ySMax   - sl.y0)   .* multFactY )  - ((sl.xSMax - sl.x0)     .* multFactX) .* rSlope;
end


function [meanDists, stdDists] = CalcMeanReachDist(allDistances, distanceProbs)

weightDists = allDistances .* distanceProbs;
weightDists = weightDists(:); % flatten
meanDists   = sum(weightDists./3);
stdDists    = sqrt(sum(((allDistances - meanDists).^2) .* distanceProbs./3, 'all') ) ./ sqrt(20);

end


function [handLocProbs, transProbs, onTableQPosX, onTableQPosY] = CalcReachProbabilities(s,allQ)


qOnTable = allQ(1,:).qVals{1}(:,:,:,2);


onTableQPosY = 62 - (1:size(qOnTable,2));
onTableQPosY = 5.* (  onTableQPosY - ((size(qOnTable,2)+1) - allQ(1,:).centerpos{1}(1))  );
onTableQPosX = 5.* (  (1:size(qOnTable,3)) - allQ(1,:).centerpos{1}(2)  );
[qX, qY]      = meshgrid(onTableQPosX,onTableQPosY);

handLocProbs = zeros(size(qX));
% Starting hand position is always at 0
handLocProbs(qX == 0 & qY == 0) = 1;
repmat(handLocProbs,[1 1 s.rch.nTP]);

% Also store transition probabilities
transProbs = zeros([size(qX), size(qX), s.rch.nTP]);


for iTP = 1:s.rch.nTP

    % loop through possible hand locations
    for iHLX = 1 : size(handLocProbs,2)
        for iHLY = 1 : size(handLocProbs,1)

            % Only calculate next possible hand locations if there is some
            % probability of being at the current hand location
            if handLocProbs(iHLY, iHLX, iTP) > 0

                % Reset the Q-values around the hand
                qAroundHand = zeros(size(qOnTable));

                % Find which positions fall 'on the table'
                nonZeroXVals = onTableQPosX >= s.rch.tableXMin - onTableQPosX(iHLX) & onTableQPosX <= s.rch.tableXMax - onTableQPosX(iHLX);
                nonZeroYVals = onTableQPosY >= s.rch.tableYMin - onTableQPosY(iHLY) & onTableQPosY <= s.rch.tableYMax - onTableQPosY(iHLY);
                qAroundHand(:,nonZeroYVals,nonZeroXVals) = qOnTable(:,nonZeroYVals,nonZeroXVals);

                % Change the values of the positions on the table to
                % reflect the experimental condition $$$
                farPositions  = onTableQPosY  > ((s.rch.tableYMax ./ 2) + s.rch.tableYMin - onTableQPosY(iHLY));
                nearPositions = onTableQPosY <= ((s.rch.tableYMax ./ 2) + s.rch.tableYMin - onTableQPosY(iHLY));
                qAroundHand(:, farPositions  ,:,:) =  qAroundHand(:, farPositions  ,:,:) .* s.rch.farVal ;
                qAroundHand(:, nearPositions ,:,:) =  qAroundHand(:, nearPositions ,:,:) .* s.rch.nearVal;

                % Remove the options to move in the z-plane
                qAroundHand(s.clc.actConsequence(:,3) ~= 0,:,:) = 0;

                % Pick the action stochastically
                % Mx across target positions first, then stochastic action selection?
                maxValAllLocs = max(qAroundHand,[],[2 3]);
                actProbs      = exp(maxValAllLocs ./ s.rch.temp) ./  sum( exp(maxValAllLocs ./ s.rch.temp) ,'all' );

                % Update future positions based on current actions
                for iAct = 1:size(actProbs,1)
                    % Find the effect of a given action
                    nextHandLocInd = [iHLY iHLX] - s.clc.actConsequence(iAct,1:2) - s.clc.baseVel(1:2);
                    % Deal with going out of bounds
                    nextHandLocInd(nextHandLocInd < 1) = 1;
                    if nextHandLocInd(1) > size(qAroundHand,2)
                        nextHandLocInd(1) = size(qAroundHand,2);
                    end
                    if nextHandLocInd(2) > size(qAroundHand,3)
                        nextHandLocInd(2) = size(qAroundHand,3);
                    end

                    % Multiply the action probability by the effect, and update the probability of ending
                    % up at a given location.
                    % Basically, this is future location probabilities, conditional on current location
                    transProbs(nextHandLocInd(1), nextHandLocInd(2), iHLY, iHLX, iTP)   =  actProbs(iAct)  ; % additional probabilities to end up in this point dure to the currently calculated action
                end

                %                 figure, imagesc(transProbs(:, :, iHLY, iHLX, iTP))
            end


        end
    end

    % Sum over all previous states, to find all future states

    handLocProbs(:,:,iTP + 1) = sum(transProbs(:,:,:,:, iTP) .* permute(handLocProbs(:, :, iTP),[4 5 1 2 3]) ,[3 4]);

end

end

function [fitRes] = CalcFitMetrics(realDat,realStd,modelDat,params)
    % ---------------------------------------------------------------------
    % Make some additional statistical measure of goodness of fit
    % First adjust all the input variance to standard error
    dataVar = [realStd{:}].^2;
    errSq   = ([realDat{:}] - [modelDat{:}]).^2;
    chiSq   = nansum(errSq(:)./dataVar(:));
    
    
    % Calculate the chi square confidence interval using bootstrapping
    chiSqDistr     = arrayfun(@(x) nansum(randsample(errSq(:)./dataVar(:),numel(errSq(:)),1)) ,[1:100000] );
    fitRes.chiSqCI = prctile(chiSqDistr,[2.5 97.5])
    
    N   = numel(errSq(:)); 
    dof = numel(params);
    k   = N - dof;
    
    % Calculate goodness of fit
    [ fitRes.gofScore chiSqCheck fitRes.pVal1 fitRes.pVal2 ] = GoFFun( errSq, dataVar, k );
    
    
    % calculate variance of residuals
    resVar     = nansum(errSq(:)./ k );
    % calculate log likelihood of data
    logLike    = -(  N                 .* log(2 .*pi .* resVar) ./ 2 )  - ...
                  (1 ./ (2.* resVar))  .* nansum(errSq(:)) ;
    % calculate AIC and BIC
    fitRes.AIC = -2 .* logLike + dof .* 2;
    fitRes.BIC = -2 .* logLike + dof .* log(N);

end


function [ gofScore chiSq pVal1 pVal2 ] = GoFFun( errSq, stdDevSq, k )
    %GoFFun Calculates the Goodness of Fit scores
    
    chiSq = nansum(errSq(:)./stdDevSq(:));
    
    gofScore = (chiSq - k) ./sqrt(2.*k);
    
    pVal1=1-chi2cdf(chiSq,k);
    pVal2=1-normcdf(gofScore);
end


function [sSqErr, sqErr, modelDat] = CalcErrSCR(realDat,params)

    % Ensure parameters are in row format first
    params  = params(:);

    % Assign the parameters to variables for easier understanding
    vPosInit = params([1,2]);
    rPos     = params([3,4]);
    alphaPos = params(5);
    gamma    = params(6);
    vNovInit = params([7,8]);
    rNov     = params([9,10]);
    alphaNov = params(11);
    offset   = params(12);
    mult     = params(13);

    % Calculate predicted values
    [vPosDuringLearn1, vNovDuringLearn1, ...
     vPosDuringTest1,  vNovDuringTest1,  ...
     vPosDuringLearn2, vNovDuringLearn2, ...
     vPosDuringTest2,  vNovDuringTest2 ] = ...
              CalculateZanini2021SCR( vPosInit, alphaPos, gamma, rPos, ...
                                      vNovInit, alphaNov, rNov);

    % Put modelled data into a cell so that it can be directly compared to
    % the real data
    modelDat = {offset + abs(vPosDuringLearn1 + vNovDuringLearn1) .* mult , ...
                offset + abs(vPosDuringTest1  +  vNovDuringTest1) .* mult ,  ...
                offset + abs(vPosDuringLearn2 + vNovDuringLearn2) .* mult , ...
                offset + abs(vPosDuringTest2  +  vNovDuringTest2) .* mult };

    if numel(realDat) == 2
        modelDat = modelDat(1:2);
    end

    % Squared error
    sqErr = cell2mat( ...
            arrayfun(@(iCond) (modelDat{iCond} - realDat{iCond}) .^ 2 , ...
                              (1 : numel(realDat)) , 'UniformOutput' , false)  );

    % Sum Square error
    sSqErr = sqrt(sum(sqErr(:)));

end


function [vPosDuringLearn1, vNovDuringLearn1, ...
          vPosDuringTest1,  vNovDuringTest1,  ...
          vPosDuringLearn2, vNovDuringLearn2, ...
          vPosDuringTest2,  vNovDuringTest2 ] = ...
              CalculateZanini2021SCR( vPosInit, alphaPos, gamma, rPos, ...
                                      vNovInit, alphaNov, rNov)
    
    % POSITION value - familiarisation
    s.nSteps            = 2;
    s.initV             = vPosInit;
    s.alpha             = alphaPos;
    s.rew               = [0 0]';
    s.gamma             = gamma;
    vPosDuringFamiliar  = applyTDLearn(s); % Familiarisation
    % POSITION value - learning 1
    s.nSteps            = 11;
    s.initV             = vPosDuringFamiliar(:,end);
    s.rew               = rPos;
    vPosDuringLearn1    = applyTDLearn(s);
    % POSITION value - testing 1
    s.nSteps            = 6;
    s.initV             = vPosDuringLearn1(:,end);
    s.alpha             = alphaPos;
    s.rew               = [0 0]';
    s.gamma             = gamma;
    vPosDuringTest1     = applyTDLearn(s); % Testing 1
    % POSITION value - learning 2
    s.nSteps            = 11;
    s.initV             = vPosDuringTest1(:,end);
    s.rew               = rPos;
    vPosDuringLearn2    = applyTDLearn(s);
    % POSITION value - testing 2
    s.nSteps            = 6;
    s.initV             = vPosDuringLearn2(:,end);
    s.alpha             = alphaPos;
    s.rew               = [0 0]';
    s.gamma             = gamma;
    vPosDuringTest2     = applyTDLearn(s); % Testing 1
    
    
    % NOVELTY value - familiarisation
    s.nSteps            = 2;
    s.initV             = vNovInit;
    s.alpha             = alphaNov;
    s.rew               = [0 0]';
    s.gamma             = gamma;
    vNovDuringFamiliar  = applyTDLearn(s); % Familiarisation
    % NOVELTY value - learning
    s.nSteps            = 11;
    s.initV             = vNovInit; vNovDuringFamiliar(:,end); % $$$ OR do I reset this to the starting value??
%     s.initV             = vNovDuringFamiliar(:,end); % $$$ OR do I reset this to the starting value??
    s.rew               = rNov;
    vNovDuringLearn1    = applyTDLearn(s);
    % NOVELTY value - testing 1
    s.nSteps            = 6;
    s.initV             = vNovInit; vNovDuringLearn1(:,end); % $$$ OR do I reset this to the starting value??
%     s.initV             = vNovDuringLearn1(:,end); % $$$ OR do I reset this to the starting value??
    s.alpha             = alphaNov;
    s.rew               = [0 0]';
    s.gamma             = gamma;
    vNovDuringTest1     = applyTDLearn(s); % Testing 1
    % NOVELTY value - learning
    s.nSteps            = 11;
    s.initV             = vNovInit; vNovDuringTest1(:,end); % $$$ OR do I reset this to the starting value??
%     s.initV             = vNovDuringTest1(:,end); % $$$ OR do I reset this to the starting value??
    s.rew               = rNov;
    vNovDuringLearn2    = applyTDLearn(s);
    % NOVELTY value - testing 2
    s.nSteps            = 6;
    s.initV             = vNovInit; vNovDuringLearn2(:,end); % $$$ OR do I reset this to the starting value??
%     s.initV             = vNovDuringLearn2(:,end); % $$$ OR do I reset this to the starting value??
    s.alpha             = alphaNov;
    s.rew               = [0 0]';
    s.gamma             = gamma;
    vNovDuringTest2     = applyTDLearn(s); % Testing 1

end


function [outVals] = applyTDLearn(s)

    for iStep = 1 : s.nSteps + 1
        if iStep == 1
            outVals(:,iStep) = s.initV;
        else
            outVals(:,iStep) = (1 - s.alpha) .*                    outVals(:,iStep-1) + ...
                                    s.alpha .* (s.rew + s.gamma .* outVals(:,1));
        end
    end
    % Remove the initial value
    outVals(:,1) = [];

end