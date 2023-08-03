function [utilHP] = HitProbToUtility(objHP,s)
% [utilHP] = HitProbToUtility(objHP)
% Calculate 'utility' hit-probs, following Straka's method
% objHP:    Objective hit probability
% utilHP:   Hitprobability for maximum utility
% s.clc.FP: Weighting of false positives 
% s.clc.FN: Weighting of false negatives

objNoHP = 1 - objHP;

lossFun = @(yPred,y) s.clc.FP .* (  max([zeros(size(yPred)) , yPred-y],[],2).^2  ) + ...
                     s.clc.FN .* (  max([zeros(size(yPred)) , y-yPred],[],2).^2  );

% Find predicted hit probabilities (y) which maximise 'utility'(based on
% the loss function that arises from the the wrighted importances of false
% positives and false negatives
yPreds      = [0:0.01:1]';
yPredMat    = repmat(yPreds, size(objHP));
yPredMatHit = repmat(lossFun(yPreds,1), size(objHP));
yPredMatMiss= repmat(lossFun(yPreds,0), size(objHP));
utilityLoss = objHP    .* yPredMatHit + ...
              objNoHP  .* yPredMatMiss ; %lossFun(ones(size(body)))


[minHP minInd] = min(utilityLoss);

utilHP = yPredMat(minInd);



end