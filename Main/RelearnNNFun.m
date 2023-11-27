function [net, Qtable, perf] = RelearnNNFun(s,w,storedEps,net,Qtable,extraInput)
% [net, Qtable, perf] = RelearnNNFun(s,w,storedEps,net,Qtable,extraInput)
% 
% Learns or relearns value Q from a set of stored [actionsa & state transitions]
% 
% -------------------------------------------------------------------------
% INPUTS:
% ----------
% s:          Settings
% w:          Description of agent's environment (world)
% storedEps:  History of actions and results --> used to train 'net' and/or
%             update Qtable
% net:        Neural network to be trained (can be pre-trained)
% Qtable:     [OPTIONAL] Provides q values in tabular format. Only used if
%             appropriate setting is chosen in 's' input
% extraInput: [OPTIONAL] Stores the value of additional environmental objects
%
% -------------------------------------------------------------------------
% OUTPUTS:
% ----------
% net:        Neural network after training
% Qtable:     Table of learned q values, relating states and actions
% perf:       [OPTIONAL] Performance of network on the trained task
%
% ----------
% NOTE:
% If s.fl.trainNet   == 1, it trains  'net'
% if s.fl.trainTable == 1, it updates 'Qtable'
%
%


tic

if ~exist('extraInput','var')
    extraInput = [];
end

if s.fl.newNet==1
    net=feedforwardnet(s.lp.netS);
    net=ChangeNeuronType(net,s); % Change neuron types
    net.performFcn= 'mse';
    net=init(net);
    net.trainFcn = s.lp.TrnFn;
    net.trainParam.epochs=s.lp.TrnEp;
    net.trainParam.showWindow = s.lp.ShwW;
    net.trainParam.max_fail=s.lp.mxF;
    firstPass=1;
    firstTrainNum=s.lp.n1stTr;
else
    firstPass=0;
end
retrEps=1000;

% Define learning rate in detail
if s.fl.mo==1
    alphaStart=0.5;
else
    alphaStart=0.3;
end
alpha=alphaStart;
% adaptive alpha function
alphaAdptPar=100;
% alphaAdptFn=@(x,y) exp(-(x)./y);
alphaAdptFn=@(x,y) 1;

parFlag=0;

eligibleEps=find(~isnan(storedEps.A(:,1)));

perf.rewPerAct = []; % This is to initialise the performance structure

clear y yy

%%
tic
for kBatch=1:s.rl.maxRetr 
    
    if rem(kBatch,1)==0 & s.fl.trainNet==1
        kBatch
    end
    
    if kBatch==1 & firstPass==1
        batchSize=s.lp.b1Siz;
    else
        batchSize=s.lp.bSiz;
    end
    
    % batch Trials
    bTrls=datasample(eligibleEps,batchSize);
    % If tool presence is being investigated, ensure that the first batch
    % includes at least 1 epoch with a tool present
    if kBatch == 1 && s.fl.ToolChange == 1 && s.lmb.ToolProb > 0 && s.fl.hist == 0
        tmpS=cell2mat(storedEps.S);
        toolInd = find(tmpS(:,4) == 1,1);
        bTrls(1) = toolInd;
    end
    % Likewise, ensure that the first batch includes a held object,if
    % appropriate
    if kBatch == 1 && s.act.eatRew == 1 && s.fl.hist == 0
        tmpS=cell2mat(storedEps.S);
        holdInd = find(tmpS(:,4 + s.fl.ToolChange) == 1,1);
        bTrls(2) = holdInd;
    end
    batchEps.prvS = storedEps.prvS(bTrls,: );
    batchEps.S = storedEps.S(bTrls,: );
    batchEps.A = storedEps.A(bTrls,: );
    batchEps.R = storedEps.R(bTrls,: );

    % Store the actions in the next state, in case of on-policy (e.g. SARSA) learning
    if strcmp(s.lp.alg,'SARSA')
        nextTrl         = min(bTrls + 1,size(storedEps.prvS,1));
        batchEps.nextA  = storedEps.A(nextTrl); 
    end
    
    if s.fl.os==1
        if firstPass==1
            osTimes=s.lp.nOs;
        else
            osTimes=s.lp.nOsBatch;
        end
        for kOS=1:s.lp.nOs
            highRewSamps=abs(batchEps.R)>s.lp.hiRewDef;
            % Oversample rewarding cases for the first pass
            bsTmp=batchSize-sum(highRewSamps);
            
            bTrls=datasample(eligibleEps,bsTmp);
            batchEps.prvS(~highRewSamps) = storedEps.prvS(bTrls,: );
            batchEps.S(~highRewSamps) = storedEps.S(bTrls,: );
            batchEps.A(~highRewSamps) = storedEps.A(bTrls,: );
            batchEps.R(~highRewSamps) = storedEps.R(bTrls,: );

            if strcmp(s.lp.alg,'SARSA')
                nextTrl                        = min(bTrls + 1,size(storedEps.prvS,1));
                batchEps.nextA(~highRewSamps)  = storedEps.A(nextTrl);
            end

        end
    end
    
    
    % flatten states for input to neural network
    x = cell2mat(batchEps.prvS); % pre states
    xx = cell2mat(batchEps.S); % post states
    A=batchEps.A; % Action
    R=batchEps.R; % Reward
    
    if s.fl.trainNet==1
        
        if s.fl.mo==1
            tmpR=nan(size(R));
            tmpR=repmat(tmpR,[1 s.act.numA]);
            
            % -------------------------------------------------------------
            % update conflicting responses:
            [inpTypes inpPos]=unique([x ],'rows');
            for kAct=1:s.act.numA
                tmpR((A==kAct),kAct)=R(A==kAct);
            end
            
            R=tmpR;
        end

        if firstPass==1
            
            net.trainParam.epochs=firstTrainNum;
            
            if s.fl.mo==1
                if parFlag==1
                    net=train(net,x', R','useParallel','yes');
                else
                    net=train(net,x', R' );
                end
            else
                if parFlag==1
                    net=train(net,[x  A]', R','useParallel','yes');
                else
                    net=train(net,[x  A]', R' );
                end
            end
            firstPass=0;
            %     net.trainFcn = 'traingd';
            net.trainFcn = 'trainlm';
            net.trainParam.showWindow = 0;
            net.trainParam.epochs=retrEps; %net.trainParam.epochs=10;
            
            % Make it able to learn a bit better for vision, which is much
            % tougher (I think)
            if s.fl.vis == 1
                net.trainParam.epochs=50;
                net.trainParam.mu=10;
                %             net.trainParam.max_fail=20;
            end
        end
        
        % predict Q-values of prestates
        if s.fl.mo==1
            y = net(x')';
        else
            y = net([x A]')';
        end

        clear yy
        
        % Update learning rate
        alpha=alphaStart*alphaAdptFn(kBatch,alphaAdptPar);
        
        if alpha>0
            % predict Q-values of poststates
            if s.fl.mo==1
                yy = net(xx')';
            else
                for kAction=1:s.act.numA
                    yy(:,kAction) = net([xx kAction*ones(size(xx,1),1)]')';
                end
            end
            
            % Training algorithm. Options: 'Q','SARSA'
            switch s.lp.alg
                case 'Q'
                    % maximum Q-value for each poststate
                    yynext = max(yy, [], 2);
                case 'SARSA'
                    % Q-value of next action that was taken
                    yynext = yy(sub2ind(size(yy), (1:size(yy, 1))', batchEps.nextA));
            end
        end
        
        % -----------------------------------------------------------------
        % calculate discounted future reward for each state transition in batch
        if w.runEnd==1 || s.lp.gamma == 0
            y = R;
        else
            if s.fl.mo==1
                for kAct=1:s.act.numA
                    y(~isnan(R(:,kAct)),kAct) = y(~isnan(R(:,kAct)),kAct) + ...
                        alpha.* ...
                        (R(~isnan(R(:,kAct)),kAct) + s.lp.gamma.* yynext(~isnan(R(:,kAct))) - y(~isnan(R(:,kAct)),kAct) );
                end
            else
                y = y + alpha.*(R + s.lp.gamma.* yynext -y );
            end
        end
        

        % -----------------------------------------------------------------
        % train the network (copied from nntrain())
        if s.fl.mo==1
            net=train(net,x',y');
            y2 = net(x')'; % $$ For debugging
        else
            net=train(net,[x A]',y');
            y2 = net([x A]')';
        end
       
    % Leave output network empty if not training
    else
        net=[];
    end
    
    if s.fl.trainTable==1
        
        
        if isempty(s.lmb.startRow)
            prvRow=size(w.world2D,1)-2; 
            row=size(w.world2D,1)-2; 
        else
            row = s.lmb.startRow;
            prvRow=s.lmb.startRow;
        end
        
        try
            for iAct=1:length(A)

                % [tbld: double checl the indexing here]
                col=batchEps.S{iAct}(1);
                goalR=batchEps.S{iAct}(2);
                goalC=batchEps.S{iAct}(3);
                prvCol=batchEps.prvS{iAct}(1);
                prvgoalR=batchEps.prvS{iAct}(2);
                prvgoalC=batchEps.prvS{iAct}(3);
                action=batchEps.A(iAct);
                rewardVal=batchEps.R(iAct);
                
                
                % update Q table to compare to NN Q for later use
                Qtable(prvRow,prvCol,prvgoalR,prvgoalC,action) = ...
                    Qtable(prvRow,prvCol,prvgoalR,prvgoalC,action) + s.lp.alpha* ...
                    (rewardVal + s.lp.gamma*max(Qtable(row,col,goalR,goalC,:)) - ...
                    Qtable(prvRow,prvCol,prvgoalR,prvgoalC,action));
            end
        catch
            warning('For now, to plot the learning table, disable all special options')
            Qtable=[];
        end
    else
        Qtable=[];
    end
    
    
    if s.fl.perfWhileLearn == 1 && mod(kBatch-1,s.prf.skipBatches) == 0
        [perf.rewPerAct(kBatch,:),perf.rewPerActSD(kBatch,:),perf.rewSums(:,:,kBatch),~] = CalcNetPerf(s,net,extraInput);
    end
    
end
toc

end


