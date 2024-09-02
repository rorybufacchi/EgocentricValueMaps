
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
    for iD = 1 : numel(realDat)./2
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



% % % FunToOpt = @(params) CalcErr(realDat,[ [-1 -.5]';  p0]);

% This is just optimising everything
FunToOpt = @(params) CalcErr(realDat,params);

% Select jparameters for fitting
fitParamIDs = [3 4 5 7 11 12 13];
% Repeating params 7 means that the novelty score is the same for both position
FunToOpt = @(params) CalcErr(realDat,[p0(1:2); params(1:3); p0(6) ; params([4 4]); p0(9:10); params(5:7)]);
% % % fitParamIDs = [1 2 3 4 5 7 11 12 13];
% % % % Repeating params 7 means that the novelty score is the same for both position
% % % FunToOpt = @(params) CalcErr(realDat,[params(1:5); p0(6) ; params([6 6]); p0(9:10); params(7:9)]);

% Set optimisation options - keep it simple
A = []; b = []; Aeq = []; beq = []; a = tic;
% Run optimisation
tic
OPTIONS = optimset('TolCon',1e-10);
[params,sSqErr,exitflag,output,lambda,grad,hessian] = fmincon(FunToOpt,p0(fitParamIDs),A,b,Aeq,beq,lb(fitParamIDs),ub(fitParamIDs),[],OPTIONS)
toc

% % % [sSqErr, SqErr, modelDat] = CalcErr(realDat,p);
[sSqErr, SqErr, modelDat] = CalcErr(realDat,[p0(1:2); params(1:3); p0(6) ; params([4 4]); p0(9:10); params(5:7)]);
% % % [sSqErr, SqErr, modelDat] = CalcErr(realDat,[params(1:5); p0(6) ; params([6 6]); p0(9:10); params(7:9)]);

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
onTableQPosY = 5.* (  onTableQPosY - (62 - allQ(9,:).centerpos{1}(1))  );
onTableQPosX = 5.* (  (1:size(qOnTable,3)) - allQ(9,:).centerpos{1}(3)  );
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


%% $$$ Let's try making it a step-by-step thing?
% % % 
% % % % $$$ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % % % $$$ MAYBE LATER STORE ALL STIMPOSQ FOR EACH TIMEPOINT AND POSITION?
% % % % $$$ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % % 
% % % temp = 0.01;
% % % nTP  = 2;
% % % 
% % % handLocProbs = zeros(size(qX));
% % % % Starting hand position is always at 0
% % % handLocProbs(qX == 0 & qY == 0) = 1;
% % % repmat(handLocProbs,[1 1 nTP]);
% % % 
% % % % Also store transition probabilities
% % % transProbs = zeros([size(qX), size(qX), nTP]);
% % % 
% % % % Possibly make q-values with reward directly under the hand
% % % qAroundHand = qOnTable;
% % % 
% % % 
% % % for iTP = 1:nTP
% % % 
% % %     % loop through possible hand locations
% % %     for iHLX = 1 : size(handLocProbs,2)
% % %         for iHLY = 1 : size(handLocProbs,1)
% % % 
% % %             % Only calculate next possible hand locations if there is some
% % %             % probability of being at the current hand location
% % %             if handLocProbs(iHLY, iHLX, iTP) > 0
% % % 
% % %                 curHandLoc = [onTableQPosX(iHLX) onTableQPosY(iHLY)];
% % % 
% % %                 % Find action values at positions in question, relative to current hand position
% % %                 stimPosQ = zeros([43,6,7]);
% % %                 for iSX = 1:numel(stimPosX)
% % %                     xBin = findnearest(stimPosX(iSX) - onTableQPosX(iHLX), onTableQPosX ); % $$$ HERE I have to subtract or add the current position, but not sure which
% % %                     for iSY = 1:numel(stimPosY)
% % %                         yBin = findnearest(stimPosY(iSY) - onTableQPosY(iHLY) , onTableQPosY );
% % %                         stimPosQ(:,iSY,iSX) = qAroundHand(:,yBin,xBin);
% % %                     end
% % %                 end
% % % 
% % %                 % Remove the options to move in the z-plane
% % %                 stimPosQ(sHnd.clc.actConsequence(:,3) ~= 0,:,:) = 0;
% % % 
% % %                 % Pick the action stochastically
% % %                 actProbs                     = exp(stimPosQ./temp) ./ sum( exp(stimPosQ./temp) );
% % %                 % Update future positions based on current actions
% % %                 for iSX = 1:numel(stimPosX)
% % %                     for iSY = 1:numel(stimPosY)
% % %                         for iAct = 1:size(actProbs,1)
% % %                             % Find the effect of a given action
% % %                             nextHandLocInd = [iHLY iHLX] + sHnd.clc.actConsequence(iAct,1:2) - sHnd.clc.baseVel(1:2);
% % %                             % Multiply the action probability by the effect, and update the probability of ending
% % %                             % up at a given location.
% % %                             % Basically, this is future location probabilities, conditional on current location
% % %                             transProbs(nextHandLocInd(1), nextHandLocInd(2), iHLY, iHLX, iTP)   =   ...
% % %                                 transProbs(nextHandLocInd(1), nextHandLocInd(2), iHLY, iHLX, iTP) + ...
% % %                                 actProbs(iAct,iSY,iSX)  ;
% % %                         end
% % %                     end
% % %                 end
% % %                 % renormalize the probabilities (because we had numel(stimPosX) * numel(stimPosY) separate scenarios, with uniform probability)
% % %                 transProbs(:, :, iHLY, iHLX, iTP) = transProbs(:, :, iHLY, iHLX, iTP) ./ (numel(stimPosX) * numel(stimPosY));
% % % 
% % %                 figure,imagesc(squeeze(max(stimPosQ)))
% % % 
% % %                 figure, imagesc(transProbs(:, :, iHLY, iHLX, iTP))
% % %             end
% % % 
% % % 
% % %         end
% % %     end
% % % 
% % %     % Sum over all previous states, to find all future states
% % % 
% % %     handLocProbs(:,:,iTP + 1) = sum(transProbs(:,:,:,:, iTP) .* permute(handLocProbs(:, :, iTP),[4 5 1 2 3]) ,[3 4])
% % % 
% % % end


%% $$$ NEW step-by-step version, SIMPLIFIED

% $$$ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% $$$ MAYBE LATER STORE ALL STIMPOSQ FOR EACH TIMEPOINT AND POSITION?
% $$$ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

temp = 0.7;
nTP  = 5;

tableXMin = -88.56/2 ;
tableXMax =  88.56/2 ;
tableYMin =   0      ;
tableYMax =  49.81   ;


handBasePos = allQ2(1,:).centerpos{1}(1:2)';

handLocProbs = zeros(size(qX));
% Starting hand position is always at 0
handLocProbs(qX == 0 & qY == 0) = 1;
repmat(handLocProbs,[1 1 nTP]);

% Also store transition probabilities
transProbs = zeros([size(qX), size(qX), nTP]);


for iTP = 1:nTP

    % loop through possible hand locations
    for iHLX = 1 : size(handLocProbs,2)
        for iHLY = 1 : size(handLocProbs,1)

            % Only calculate next possible hand locations if there is some
            % probability of being at the current hand location
            if handLocProbs(iHLY, iHLX, iTP) > 0

                % Reset the Q-values around the hand
                qAroundHand = zeros(size(qOnTable));

                curHandLoc = [onTableQPosX(iHLX) onTableQPosY(iHLY)];

                % Find which positions fall 'on the table'
% % %                 yOffset = handBasePos(1) - iHLY;
% % %                 xOffset = handBasePos(1) - iHLX;

                % $$$ THIS NEEDS TO BE UPDATED TO TAKE THE CORRECT ZONE
                % when the current hand location changes
                nonZeroXVals = onTableQPosX >= tableXMin - onTableQPosX(iHLX) & onTableQPosX <= tableXMax - onTableQPosX(iHLX);
                nonZeroYVals = onTableQPosY >= tableYMin - onTableQPosY(iHLY) & onTableQPosY <= tableYMax - onTableQPosY(iHLY);
                qAroundHand(:,nonZeroYVals,nonZeroXVals) = qOnTable(:,nonZeroYVals,nonZeroXVals);

                % Remove the options to move in the z-plane
                qAroundHand(sHnd.clc.actConsequence(:,3) ~= 0,:,:) = 0;

                % Pick the action stochastically
% % %                 % Either stochastically only within target position
% % %                 actProbs = exp(qAroundHand./temp) ./ sum( exp(qAroundHand./temp) );
% % %                 % OR stochastically across ALL target positions too!
% % %                 actProbs = exp(qAroundHand./temp) ./  sum( exp(qAroundHand./temp) ,'all' );
% % %                 actProbs = sum(actProbs,[2 3]);
                % $$$ OR do I do some kind of max across target positions
                % first, and then perform the stochastic action selection?
                % GAH SO MANY OPTIONS!!
                % --> I guess I just try --> let's stick with this for now
                maxValAllLocs = max(qAroundHand,[],[2 3]);
                actProbs      = exp(maxValAllLocs./temp) ./  sum( exp(maxValAllLocs./temp) ,'all' );

                % Update future positions based on current actions
                for iAct = 1:size(actProbs,1)
                    % Find the effect of a given action
                    nextHandLocInd = [iHLY iHLX] - sHnd.clc.actConsequence(iAct,1:2) - sHnd.clc.baseVel(1:2);
                    % Deal with going out of bounds
                    nextHandLocInd(nextHandLocInd < 1) = 1;
                    if nextHandLocInd(1) > 61
                        nextHandLocInd(1) = 61;
                    end
                    if nextHandLocInd(2) > 31
                        nextHandLocInd(2) = 31;
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

    handLocProbs(:,:,iTP + 1) = sum(transProbs(:,:,:,:, iTP) .* permute(handLocProbs(:, :, iTP),[4 5 1 2 3]) ,[3 4])

end


%%

% $$$ Consider adding an option to put a function handle in to "CalcErr", so that I
% can re-use the CalcErrAndGoF for the other dataset
fitRes = CalcFitMetrics(realDat,realStd,modelDat,p0)




%% Functions


% function Bol

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


function [sSqErr, sqErr, modelDat] = CalcErr(realDat,params)

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