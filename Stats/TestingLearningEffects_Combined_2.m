
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

%% Try to find the optimal parameters


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

% Select parameters for fitting
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

%% plot the data
allTitles = {'Learning phase 1','Testing phase1','Learning phase2','Testing phase2'};
colours = [0.8 0.3 0.3 ; 0.3 0.3 0.8];
f.LearnSCR.f = figure('Position', [ 50, 50, 400, 900]),

colorOrd = [ 1 0.2 0.2; 0.2 0.2 1 ];
for iD = 1:numel(realDat)

%     subplot(numel(realDat),1,iD)
    subplot(numel(realDat),1,iD)
    errorbar(1:size(realDat{iD},2) , realDat{iD}', realStd{iD}', '-o');

% % %     for iCond = 1:2
% % %         opts.c = colorOrd(iCond,:);
% % %         opts.plotSpec = {'--','Color',opts.c};
% % %         ShadedPlot(1:size(realDat{iD},2), realDat{iD}(iCond,:) , realDat{iD}(iCond,:) - realStd{iD}(iCond,:), realDat{iD}(iCond,:) + realStd{iD}(iCond,:), opts);
% % %         hold on;
% % %     end


    hold on;
    plot(modelDat{iD}','-o','LineWidth',2);
    title(allTitles{iD});
    legend('near','far')

    colororder(colours)

    legend('Real Near','Modelled Near','Real Far','Modelled Far')
    ylabel('SCR response')
    xlabel('trial number')
end


%% Save figures if called for
allFields = fields(f);
for iF = 1:length(allFields)
    cF = allFields{iF};
    set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\LearningEffects\' cF '.eps'] , 'epsc')
    saveas(f.(cF).f,['C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\LearningEffects\' cF '.pdf'] , 'pdf')
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

allAlphas = [0.1:0.1:0.9 , 0.07 0.05 0.03];

allTemps = [ 0.05 0.07 0.1 0.14 0.2 0.3];

% $$$ GOT HERE: iTemp 1, iAlpha 7, iLEarnTP 7

% $$$ In the meantime, for model fitting and stuff, just set
% $$$ iTemp 2, iAlpha 4, iLearnTP 1:40
for iTemp = 1:numel(allTemps)

    sHnd2.rch.temp = allTemps(iTemp);

    for iAlpha = [10 11 12] %1:numel(allAlphas)

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
    save('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Data\Fitting_Learning_Workspace_v4.mat','-v7.3')
end


%% $$$ Just make sure the Control values are properly filled out
for iTemp = 1:6

    sHnd2.rch.temp = allTemps(iTemp);
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

    handLocProbs_Learn(:,:,:,:,:,iTemp) =  repmat(handLocProbs,[1 1 1 size(handLocProbs_Learn_FarRew,4) size(handLocProbs_Learn_FarRew,5)]);

end
%% Reproduce the findings from the paper:
%  Find the optimal parameters

% sHnd2.rch.fitReachLinearly          = 1;
sHnd2.rch.fitReachLinearly          = 0;
sHnd2.rch.renormaliseOnTableSurface = 1;

sHnd2.rch.plotFl     =    0;

sHnd2.rch.cDiffuseTP =   20;
sHnd2.rch.cAlpha     =    4;
sHnd2.rch.cTemp      =    2;
sHnd2.rch.cLearnTP   = 1:40; 
sHnd2.rch.gaussKern  =    6; 
sHnd2.rch.gaussSigm  =    1;


sHnd2.rch.nFreeParams = 0;

% % % [fitRes, normSSqErr, realDat, realDatStd, modelDat] = ReproduceReachingExperiment(sHnd2,...
% % %     handLocProbs_Learn,handLocProbs_Learn_FarRew,handLocProbs_Learn_NearRew,...
% % %     qX,qY,onTableQPosX,onTableQPosY);


% Loop through the fitting parameters

% % % clear normSSqErrAll gofScoreAll pValAll



for iGaussSigm = 1:5

    sHnd2.rch.gaussSigm  = iGaussSigm - 1;
    % ------------------------------------------------
    % Make gaussian blurred version
    handLocProbs_Learn_G         = GaussBlur(handLocProbs_Learn,        sHnd2.rch.gaussKern,sHnd2.rch.gaussSigm );
    handLocProbs_Learn_FarRew_G  = GaussBlur(handLocProbs_Learn_FarRew, sHnd2.rch.gaussKern,sHnd2.rch.gaussSigm );
    handLocProbs_Learn_NearRew_G = GaussBlur(handLocProbs_Learn_NearRew,sHnd2.rch.gaussKern,sHnd2.rch.gaussSigm );

    for iAlpha = 1:size(handLocProbs_Learn_FarRew,5)
        for iDiffuseTP = 1:21
            for iTemp  = 1:size(handLocProbs_Learn_FarRew,6)

                sHnd2.rch.cAlpha     = iAlpha ;
                sHnd2.rch.cDiffuseTP = iDiffuseTP;
                sHnd2.rch.cTemp      = iTemp;


                [fitRes, normSSqErrAll(iAlpha,iDiffuseTP, iTemp, iGaussSigm), realDat, realDatStd, modelDat] = ...
                    ReproduceReachingExperiment(sHnd2,...
                    handLocProbs_Learn_G,handLocProbs_Learn_FarRew_G,handLocProbs_Learn_NearRew_G,...
                    qX,qY,onTableQPosX,onTableQPosY);


                gofScoreAll(iAlpha,iDiffuseTP, iTemp, iGaussSigm) = fitRes.gofScore;
                pValAll(iAlpha,iDiffuseTP, iTemp, iGaussSigm)     = fitRes.pVal1;

            end
        end
        
    end

    % save('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Data\Fitting_Learning_Workspace_CheckedFitting.mat_V3','-v7.3')
    save('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Data\Fitting_Learning_Workspace_CheckedFitting_NOTfittingReachabilityDirectly.mat','-v7.3')

    iAlpha
    iDiffuseTP
    iTemp
    iGaussSigm
end


%% $$$ Try to look at the best results

% normSSqErrAll(:,1,:) = NaN;
normSSqErrAll(normSSqErrAll==0) = NaN;


sHnd2.rch.svPlFl                    = 1;

sHnd2.rch.fitReachLinearly          = 1;
sHnd2.rch.renormaliseOnTableSurface = 1;

[dmyMin, minInd] = min(normSSqErrAll,[],'all')
[dmyMin, minInd] = max(pValAll,[],'all')

[rr, cc, hh1, hh2] = ind2sub(size(pValAll), minInd)

sHnd2.rch.plotFl     =    1;

sHnd2.rch.cAlpha     = rr;
sHnd2.rch.cDiffuseTP = cc;
sHnd2.rch.cTemp      = hh1;
sHnd2.rch.gaussSigm  = hh2 - 1; % Because there's also zero gaussian

% pValTmp = pValAll;
% % % % pValTmp(:,:,:,2:end) = 0;
% pValTmp(:,:,:,[1 3:end]) = 0;
% [dmyMin, minInd] = max(normSSqErrAll,[],'all')
% % [dmyMin, minInd] = max(pValTmp,[],'all')
% [rr, cc, hh1, hh2] = ind2sub(size(pValTmp), minInd)
% 
% sHnd2.rch.plotFl     =    1;
% 
% sHnd2.rch.cAlpha     = rr;
% sHnd2.rch.cDiffuseTP = cc;
% sHnd2.rch.cTemp      = hh1;
% sHnd2.rch.gaussSigm  = hh2 - 1; % Because there's also zero gaussian


% % % rr = 1;
% % % cc = 15;
% % % hh1 = 1;
% % % hh2 = 5;


% ------------------------------------------------
% Make gaussian blurred version
handLocProbs_Learn_G         = GaussBlur(handLocProbs_Learn,        sHnd2.rch.gaussKern,sHnd2.rch.gaussSigm );
handLocProbs_Learn_FarRew_G  = GaussBlur(handLocProbs_Learn_FarRew, sHnd2.rch.gaussKern,sHnd2.rch.gaussSigm );
handLocProbs_Learn_NearRew_G = GaussBlur(handLocProbs_Learn_NearRew,sHnd2.rch.gaussKern,sHnd2.rch.gaussSigm );

%% ------------------------------------------------
% Check the model fits
[fitRes, sSqNormErr, realDat, realDatStd, modelDat] = ...
    ReproduceReachingExperiment(sHnd2,...
    handLocProbs_Learn_G,handLocProbs_Learn_FarRew_G,handLocProbs_Learn_NearRew_G,...
    qX,qY,onTableQPosX,onTableQPosY);


sSqNormErr

fitRes.pVal1



%% Functions


function [fitRes, sumErrSq, realDat, realDatStd, modelDat] = ReproduceReachingExperiment(s,handLocProbs_Learn,handLocProbs_Learn_FarRew,handLocProbs_Learn_NearRew,...
    qX,qY,onTableQPosX,onTableQPosY)

% =========================================================================
% Some default settings
blockXsize = 88.56/7;
blockYsize = 49.81/6;
stimPosX   = blockXsize .* -4   + blockXsize .* (1:7);
stimPosY   = 2.5 + blockYsize .* -0.5 + blockYsize .* (1:6); 


handLocProbs_Learn_G = handLocProbs_Learn;
handLocProbs_Learn_FarRew_G = handLocProbs_Learn_FarRew;
handLocProbs_Learn_NearRew_G = handLocProbs_Learn_NearRew;

% =========================================================================
% Renormalise probability to just the table's surface area
if s.rch.renormaliseOnTableSurface == 1

    onTable = qY >= stimPosY(1) - blockYsize./2 & qY <= stimPosY(end) + blockYsize./2 & ...
              qX >= stimPosX(1) - blockXsize./2 & qX <= stimPosX(end) + blockXsize./2;
    
    % When normalising, add a tiny tiny offset so that we don't get NaNs
    % --> just allow probability to be 0 instead
    handLocProbs_Learn_G         = handLocProbs_Learn_G         .* onTable;
    handLocProbs_Learn_G         = handLocProbs_Learn_G         ./ (sum(handLocProbs_Learn_G,[1 2]) + 10e-100) ;
    handLocProbs_Learn_FarRew_G  = handLocProbs_Learn_FarRew_G  .* onTable;
    handLocProbs_Learn_FarRew_G  = handLocProbs_Learn_FarRew_G  ./ (sum(handLocProbs_Learn_FarRew_G,[1 2]) + 10e-100) ;
    handLocProbs_Learn_NearRew_G = handLocProbs_Learn_NearRew_G .* onTable;
    handLocProbs_Learn_NearRew_G = handLocProbs_Learn_NearRew_G ./ (sum(handLocProbs_Learn_NearRew_G,[1 2]) + 10e-100) ;
end

% =========================================================================
% Plot and reproduce 2D reaching results

% Control
pC = squeeze(mean(handLocProbs_Learn_G(:,:,s.rch.cDiffuseTP,s.rch.cLearnTP,s.rch.cAlpha,s.rch.cTemp), 4));
% Far
pF = squeeze(mean(handLocProbs_Learn_FarRew_G(:,:,s.rch.cDiffuseTP,s.rch.cLearnTP,s.rch.cAlpha,s.rch.cTemp), 4));
% Near
pN = squeeze(mean(handLocProbs_Learn_NearRew_G(:,:,s.rch.cDiffuseTP,s.rch.cLearnTP,s.rch.cAlpha,s.rch.cTemp), 4));


% Find probabilities of ending up at certain positions on the table
clear endProbs
for iSX = 1:numel(stimPosX)
    xBinMin = stimPosX(iSX) - blockXsize./2 ;
    xBinMax = stimPosX(iSX) + blockXsize./2 ;
    for iSY = 1:numel(stimPosY)
%         if iSY == 4 & iSX == 6
%             disp('test')
%         end
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

% Extra backup for small rounding errors, btu it shouldn't do much
if s.rch.renormaliseOnTableSurface == 1
    pResC = pResC ./ sum(pResC(:));
    pResF = pResF ./ sum(pResF(:));
    pResN = pResN ./ sum(pResN(:));
end

if s.rch.plotFl == 1
    f.reachHeatMaps.f =figure('Position' , [100, 100, 1200, 250]);
    subplot(1,3,1)
    imagesc(stimPosX,stimPosY,pResF)
    title('Far group')
    axis xy
    colorbar
    clim([0 0.075])

    subplot(1,3,2)
    imagesc(stimPosX,stimPosY,pResC)
    title('Control group')
    axis xy
    colorbar
    clim([0 0.075])

    subplot(1,3,3)
    imagesc(stimPosX,stimPosY,pResN)
    title('Near group')
    axis xy
    colorbar
    clim([0 0.075])

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

colours = [0.8 0.3 0.3 ; 0.3 0.3 0.8;  0.8 0.3 0.3 ; 0.3 0.3 0.8];
if s.rch.plotFl == 1
    f.reachDistances.f = figure('Position',[100, 300, 1000 400]);
    subplot(1,2,1);
    errorbar([0.85 1.15 ; 1.85 2.15 ; 2.85 3.15] , realDat{1}', realDatStd{1}','.k');
    hold on;
    bar(realDat{1}');
    ylim([0 40])
    title('Real Data')
    set(gca,'XTick',[1 2 3]);
    set(gca,'XTickLabel',{'Far Rew','Control','Near Rew'});

    subplot(1,2,2);
    errorbar([0.85 1.15 ; 1.85 2.15 ; 2.85 3.15], modelDat{1}', modelDatStd{1}','.k');
    hold on;
    set(gca,'XTick',[1 2 3]);
    set(gca,'XTickLabel',{'Far Rew','Control','Near Rew'});
    h = bar(modelDat{1}');
    ylim([0 40])
    title('Modelled Data');

    colororder(colours)
    legend(h,'Early Blocks','Late Blocks');
    
end

% =========================================================================
% Plot and reproduce the 'proximal selection' as a function of block

% % % distPositions  = onTableQPosY >= ( blockYsize .* 4 + 2.5);
% % % proxPositions  = onTableQPosY <= ( blockYsize .* 2 + 2.5);
% $$$ MAYBE THIS WORKS?.. let's see...
distPositions  = onTableQPosY > ( blockYsize .* 4 + 2.5);
proxPositions  = onTableQPosY <= ( blockYsize .* 3 + 2.5);

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


realDat{2}    = [pProxRealN ; pProxRealC; pProxRealF];
modelDat{2}   = [pProxN'    ; pProxC'   ; pProxF'];
% Estimate the variability of the data 
realDatStd{2} = [repmat(std(pProxRealN(1:end)) , size(pProxRealN))  ; ...
                 repmat(std(pProxRealC(1:end)) , size(pProxRealC))  ; ...
                 repmat(std(pProxRealF(1:end)) , size(pProxRealF)) ];



% Note: shift x- axis by half a block, so that we're simulating the average
% learnign that happened during the block itself, rather than the result at
% the end of the block
% % % if s.rch.plotFl == 1
% % %     blockXAx = s.rch.cLearnTP;
% % %     figure,
% % %     hold off
% % %     plot(blockXAx - 0.5, pProxN,'b','LineWidth',2);
% % %     hold on
% % % %     plot(blockXAx, pProxRealN ,'.b','MarkerSize',20);
% % %     errorbar(blockXAx, pProxRealN ,realDatStd{2}(1,:), '.b','MarkerSize',15, 'CapSize',0);
% % %     plot(blockXAx - 0.5, pProxF, 'r','LineWidth',2);
% % % %     plot(blockXAx, pProxRealF ,'.r','MarkerSize',20);
% % %     errorbar(blockXAx, pProxRealF ,realDatStd{2}(3,:), '.r','MarkerSize',15, 'CapSize',0);
% % %     plot(blockXAx - 0.5, pProxC, 'k','LineWidth',2);
% % % %     plot(blockXAx, pProxRealC ,'.k','MarkerSize',20);
% % %     errorbar(blockXAx, pProxRealC ,realDatStd{2}(2,:), '.k','MarkerSize',15, 'CapSize',0);
% % %     ylim([0 1]);
% % % end

if s.rch.plotFl == 1
    f.PropNearResponses.f = figure('Position',[500 300 800 500]);
    colorOrd = [0.2 0.2 1 ; 0 0 0 ; 1 0.2 0.2 ];
    for iCond = 1:size(realDat{2},1)
        opts.c = colorOrd(iCond,:);
        opts.plotSpec = {'--','Color',opts.c};
        blockXAx = s.rch.cLearnTP;
%         plot(blockXAx - 0.5, modelDat{2}(iCond,:),'Color',opts.c,'LineWidth',2);
        plot(blockXAx, modelDat{2}(iCond,:),'Color',opts.c,'LineWidth',2);
        hold on;
        ShadedPlot(blockXAx, realDat{2}(iCond,:) , realDat{2}(iCond,:) - realDatStd{2}(iCond,:), realDat{2}(iCond,:) + realDatStd{2}(iCond,:), opts);
    end
    ylim([0 1]);
end
legend('Modelled Near','Real Near','Modelled Control','Real Control','Modelled Far','Real Far');


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
sl.xSMin = 208.49;
sl.xSMax = 209.378;
sl.ySMin =   0.127;
sl.ySMax =  -0.439;
% Calculate the relationship
[rProxSlope rProxInt] = FindLinEqFromFig(sl);


% ------------------------------------------------------------------
% Relationship with change of movement amplitude

% Coordinate data: figure size
sl.xMin = 209.841;
sl.xMax = 210.853;
sl.yMin =  -0.434;
sl.yMax =   0.519;
% Real data: figure size
sl.xRMin = -10;
sl.xRMax =  15;
sl.yRMin =-200;
sl.yRMax = 300;
% Coordinate data: slope points
sl.x0    = 210.244;
sl.y0    =  -0.054;
sl.xSMin = 209.841;
sl.xSMax = 210.720;
sl.ySMin =  -0.104;
sl.ySMax =   0.021;
% Calculate the relationship
[rMoveSlope rMoveInt] = FindLinEqFromFig(sl);


nDistF = mean(pDistF .* 400);
nProxF = mean(pProxF .* 400);
nDistC = mean(pDistC .* 400);
nProxC = mean(pProxC .* 400);
nDistN = mean(pDistN .* 400);
nProxN = mean(pProxN .* 400);

nMoveF = modelDat{1}(2,1) - modelDat{1}(1,1);
nMoveC = modelDat{1}(2,2) - modelDat{1}(1,2);
nMoveN = modelDat{1}(2,3) - modelDat{1}(1,3);

% Now find the predicted 'boundary' changes, and average them
bDistF = (nDistF - rDistInt) ./ rDistSlope;
bProxF = (nProxF - rProxInt) ./ rProxSlope;
bMoveF = (nMoveF - rMoveInt) ./ rMoveSlope;
bF     = (bProxF + bDistF + bMoveF)   ./ 3;
% bF     = median(bProxF + bDistF + bMoveF);

bDistC = (nDistC - rDistInt) ./ rDistSlope;
bProxC = (nProxC - rProxInt) ./ rProxSlope;
bMoveC = (nMoveC - rMoveInt) ./ rMoveSlope;
bC     = (bProxC + bDistC + bMoveC) ./ 3;
% bC     = median(bMoveC + bDistF + bMoveC) ;

bDistN = (nDistN - rDistInt) ./ rDistSlope;
bProxN = (nProxN - rProxInt) ./ rProxSlope;
bMoveN = (nMoveN - rMoveInt) ./ rMoveSlope;
bN     = (bProxN + bDistN + bMoveN)   ./ 3;
% bN     = median(bProxN + bDistF + bMoveN)  ;

realDat{3}    = [0.305 -0.051 -0.302] .* 7.886;
realDatStd{3} = [0.435 (- 0.051 + (0.139 - 0.051))  (- 0.302 + (0.421 - 0.302)) ] .* 7.886;
realDatStd{3} = realDatStd{3} - realDat{3};

modelDat{3}    = [bF bC bN];
% Either take the values directly, or fit them
if s.rch.fitReachLinearly == 1
    % Formulate the design matrix
    X            = [ones([numel(modelDat{3}) , 1]), modelDat{3}'];
    % Solve for coefficients using the normal equations
    coefficients = (X' * X) \ (X' * realDat{3}');
    % Extract the intercept and slope
    intrcpt      = coefficients(1);
    slp          = coefficients(2);

    modelDat{3} = modelDat{3} .* slp + intrcpt;
end

if s.rch.plotFl == 1
    f.ReachEstimates.f = figure('Position',[1000 300 800 500]);
    subplot(1,2,1)
    bar(realDat{3});
    hold on;
    errorbar(realDat{3}, realDatStd{3},'.k');
    title('real reaching estimate change')
    set(gca,'XTickLabels',{'Far','Control','Near'})
    ylim([-4 4])

    subplot(1,2,2)
    bar(modelDat{3});
    set(gca,'XTickLabels',{'Far','Control','Near'})
    title('modelled reaching estimate change')
    ylim([-4 4])
end

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
sumErrSq = sqrt(sum(errSq(:)));

dataVar = cell2mat(...
    arrayfun(@(iCond) realDatStdForErrorCalc{iCond}(:)' .^2 , ...
    (1 : numel(realDat))  , 'UniformOutput' , false));
chiSq   = nansum(errSq(:)./dataVar(:));

% Calculate the chi square confidence interval using bootstrapping
chiSqDistr     = arrayfun(@(x) nansum(randsample(chiSq,numel(errSq(:)),1)) ,[1:100000] );
fitRes.chiSqCI = prctile(chiSqDistr,[2.5 97.5]);

N   = numel(errSq(:));
dof = numel(s.rch.nFreeParams) + s.rch.fitReachLinearly.*2 ;
k   = N - dof;

% Calculate goodness of fit
[ fitRes.gofScore chiSqCheck fitRes.pVal1 fitRes.pVal2 ] = GoFFun( errSq, dataVar, k );


% Save figures if called for
if s.rch.svPlFl == 1
    allFields = fields(f);
    for iF = 1:length(allFields)
        cF = allFields{iF};
        set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
        saveas(f.(cF).f,['C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\LearningEffects\' cF 'V6.eps'] , 'epsc')
        saveas(f.(cF).f,['C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\ModelEmpirical\BitsAndPieces\FromMatlab\LearningEffects\' cF 'V6.pdf'] , 'pdf')
    end
end

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
    if sigma == 0
        % Allow for there to be no smoothing, i.e. 0 sigma
        % (If we don't do this, we get NaN in the center of the kernel
        gaussian_kernel = double(x == 0 & y == 0);
    else
        gaussian_kernel = exp(-(x.^2 + y.^2) / (2 * sigma^2));
    end
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:)); % Normalize the kernel

    % Store the total probability of the original matrix
    total_prob = sum(prob_matrix(:));

    % Perform the convolution (2D Gaussian blur)
    blurred_matrix = conv2(prob_matrix, gaussian_kernel, 'same');

    % Normalize the blurred matrix to maintain the total probability
    if sum(blurred_matrix(:)) == 0 | total_prob == 0
        normalized_blurred_matrix(:,:,iFurtherDim) = 0;
    else
        normalized_blurred_matrix(:,:,iFurtherDim) = (blurred_matrix / sum(blurred_matrix(:))) * total_prob;
    end
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