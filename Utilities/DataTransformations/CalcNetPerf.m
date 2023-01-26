function [rewPerAct, rewPerActSD, rewSums, sE] = CalcNetPerf(s,net,extraInput)
%CalcNetPerf Calculates network performance
%
%   outputs:
%
%   rewPerAct(:,1): total positive rewards per actions (timestep)
%   rewPerAct(:,2): total negative rewards per actions (timestep)
%   rewSums(:,1):   total positive rewards
%   rewSums(:,2):   total negative rewards
%   sE:             storedEpochs of the last-performed net performance calculation

% For older version, load possibly absent settings
s = DefaultSettings(s);

% new settings for running the performance checks
s2=s;

s2.rp.maxActions=s.prf.maxActions;
s2.fl.newNet=s.prf.newNet;
s2.fl.newTable=s.prf.newTable;
s2.lp.epsilon=s.prf.epsilon;


for iRep=1:s.prf.nRep
% %     if s2.fl.extraInput == 1
% %         [~, ~, sE, ~, ~] = RunRLfun(s2,net,[],[],extraInput);
% %     else
% %         [~, ~, sE, ~, ~] = RunRLfun(s2,net,[],[],[]);
% %     end

    [~, ~, sE, ~, ~] = RunRLfun(s2,net,[],[],extraInput);
    
    rewSums(iRep,1)=sum(sE.R(sE.R>1)); % sum of all positive rewards when changing goal neurons
    rewSums(iRep,2)=sum(sE.R(sE.R<-1));% sum of all negative rewards when changing goal neurons
end


rewPerAct = nanmean(rewSums)./s2.rp.maxActions;
rewPerActSD = nanstd(rewSums)./s2.rp.maxActions;



end

