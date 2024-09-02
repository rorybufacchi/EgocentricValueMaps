
% Initial Novelty value. $$$ ASSUME THIS RESETS EVERY TIME THE HAND IS MOVED?
vNovInit = 5 .*[-1 -1]';
% Initial Position value
vPosInit = [-1 -0.5]';


% Rewards for 'novelty', i.e. having experienced a stimulus before recently
rNov           = [ 0  0]';
% rewards for each position
rPos           = [-2 -0]';
rPosDuringTest = [ 0  0]';

% Learning rate
alphaPos = 0.2;
alphaNov = 0.6;

% Discount factor (make zero for 'end of simulation is now' situation
% (Note I'm assuming a lack of movement here, I guess)
gamma = 0;

% Habituation size. $$$ THIS NEEDS TO BE SOMETHING THAT EXPONENTIALLY
% APPROACHES A CERTAIN VALUE: like starts off at 2 and then slowly
% approaches 1 or something
% habMag = 
% $$$ OR I COULD DO IT BY LEARNING A SEPARATE EFFECT OR SOMETHING?...



% Learning timepoints
learningSteps = 1:14; % 12 + 2
startLearnStep = 3;

testingSteps  = 1:6;



% Calculate value for the familiarization and learning steps
vPosDuringLearn = zeros([numel(vPosInit) , numel(learningSteps)]);
vNovDuringLearn = zeros([numel(vNovInit) , numel(learningSteps)]);
for iStep = learningSteps
    if iStep == 1
        vPosDuringLearn(:,iStep) = vPosInit;
        vNovDuringLearn(:,iStep) = vNovInit;
    else
        if iStep < startLearnStep
            rPosTmp = [0 0]'; % Don't offer a reward for the habituation session
        else
            rPosTmp = rPos; 
        end
        vPosDuringLearn(:,iStep) = (1-alpha) .* vPosDuringLearn(:,iStep-1) + alphaPos .* (rPosTmp + gamma .* vPosDuringLearn(:,1));
        vNovDuringLearn(:,iStep) = (1-alpha) .* vNovDuringLearn(:,iStep-1) + alphaNov .* (rNov    + gamma .* vNovDuringLearn(:,1));
    end
end

% POSITION value - familiarisation
s.nSteps            = 2;
s.initV             = vPosInit;
s.alpha             = alphaPos;
s.rew               = [0 0]';
s.gamma             = gamma;
vPosDuringFamiliar  = applyTDLearn(s); % Familiarisation

% POSITION value - learning
s.nSteps            = 12;
s.initV             = vPosDuringFamiliar(:,end);
s.rew               = rPos;
vPosDuringLearn     = applyTDLearn(s);

%%



% Calculate value for the first testing phase
vPosDuringTest1 = zeros([numel(vPosInit) , numel(testingSteps)]);
vNovDuringTest1 = zeros([numel(vNovInit) , numel(testingSteps)]);
for iStep = testingSteps
    if iStep == 1
        vPosDuringTest1(:,iStep) = vPosDuringLearn(:,end); % Position effect carries over
        vNovDuringTest1(:,iStep) = vNovInit; % Novelty effect restarts?
    else
        vPosDuringTest1(:,iStep) = (1-alpha) .* vPosDuringTest1(:,iStep-1) + alphaPos .* (rPosDuringTest + gamma .* vPosDuringTest1(:,1));
        vNovDuringTest1(:,iStep) = (1-alpha) .* vNovDuringTest1(:,iStep-1) + alphaNov .* (rNov           + gamma .* vNovDuringTest1(:,1));
    end
end

figure,
subplot(2,1,1)
plot(1:12 , abs(vPosDuringLearn + vNovDuringLearn)','-o');
title('Learning phase');
legend('near','far')

subplot(2,1,2)
plot(testingSteps , abs(vPosDuringTest1 + vNovDuringTest1)','-o');
title('Testing phase');
legend('near','far')


% $$$ NEXT put the learning thing in a function, and redo the learning
% thing but for zero reward. PLUS add in a 'habituation effect

% valDuringLearn = vPosInit

%% Functions

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