function [net, Qtable, perf] = RelearnNNFunSC(s,w,storedEps,net,Qtable,extraInput)
% If s.fl.trainNet==1, use a network
% if s.fl.trainTable =1, use a table
% perf returns measures of network performance


% $$$$$$$$$$$$ I THINK TRY: really big body, small hand - should mainly
% avoid
% $$$$$$$$ NEXT THING TO DO:::::!!! check whether there are any
% simultaneous hits of hand and body - if so, figure out how to change
% that! $$$$$


% $$ THINGS TO TRY:
% NOW: 1) Add in touch. 2) Add in body. 3) Do reaction time 'study' (visual towards and away from the body - training with positive and negative valence)
% Can also add tool use. Even washout effects of tool use
% Try 'stimulating' neurons and seeing what effect is on a) behaviour and b) other neurons
% vision with different pixels (i.e. 'eye')
% Fix the display vision thing for single neurons
% Add extra limb/body with more/less reward - would have to destroy the
%   'threats/rewards' when they are captured
% Think about how to add in touch
% Using LSTM instead of history - i.e. 'memory' . Need to check out deep
%   learning toolbox for this
% Successor state representation??
% Make it a full limb that can move forward/backward as well
% With vision try adding in multiple objects

% $$ NEEED to double check the stored epochs when using the body - do it
% with only goal for simplicity first

% $$$
% One thing that makes a huge difference is the amount of stored epochs -
% the more the better definitely! With 2 speeds, and threat and goal, 4 mil
% is pretty good. For vision that's too much for the RAM though.

% Some things I can try:
% THen I will try runnnin it as an actual reinforcement learning algirthm

% $$ NEED TO find out why the plotting thing doesn't care about rowlag and
% collag


tic

if ~exist('extraInput','var')
    extraInput = [];
end

if s.fl.newNet==1
    % netSize=[40 12 12 8 8 6 6 5]; % This one for threats with different speeds (and 50 batches)
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

% $$ NOTE I SET: net.trainParam.goal=1e-6

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
    if kBatch == 1 && s.fl.ToolChange == 1 && s.lmb.ToolProb > 0 && s.fl.hist == 0 % $$$ THE HIST IS EXPERIMENTAL!!
        tmpS=cell2mat(storedEps.S);
        toolInd = find(tmpS(:,4) == 1,1);
        bTrls(1) = toolInd;
    end
    batchEps.prvS = storedEps.prvS(bTrls,: );
    batchEps.S = storedEps.S(bTrls,: );
    batchEps.A = storedEps.A(bTrls,: );
    batchEps.R = storedEps.R(bTrls,: );
    
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
        end
    end
    
    
    % flatten states for input to neural network
    x = cell2mat(batchEps.prvS); % pre states $$$$$$$$$$$$
    xx = cell2mat(batchEps.S); % post states
    A=batchEps.A; % Action
    R=batchEps.R; % Reward
    
    if s.fl.trainNet==1
        
        if s.fl.mo==1
            tmpR=nan(size(R));
            tmpR=repmat(tmpR,[1 s.act.numA]);
            
            [inpTypes inpPos]=unique([x ],'rows');
            for kAct=1:s.act.numA
                tmpR((A==kAct),kAct)=R(A==kAct);
            end
            
            % % %         % --------------------------------------------------------------
            % % %         % update conflicting responses:
            % % %         % Find events where the state is the same. Note which actions
            % % %         % were ones that were actually taken. Replace the R values of
            % % %         % the actions that weren't taken with randomly sampled ones
            % % %         % that were taken. If none of the memory has a a particular
            % % %         % state-action pair, leave the resulting R as nan, and
            % % %         % ultimately keep the updated y as itself
            % % %         for kInp=1:numel(inpPos)
            % % %             sameInp=find(min(x==inpTypes(kInp,:),[],2));
            % % %             for kAct=1:s.act.numA
            % % %                 sameStateR=tmpR(sameInp,kAct);
            % % %                 % Replace nans with actual observed values.
            % % %                 if sum(~isnan(sameStateR))==0
            % % %                     % If no values observed, do what - LEAVE NAN?? I think so
            % % %                 elseif sum(~isnan(sameStateR))==1
            % % %                     sameStateR(isnan(sameStateR))=sameStateR(~isnan(sameStateR));
            % % %                 else
            % % %                     sameStateR(isnan(sameStateR))=randsample(sameStateR(~isnan(sameStateR)),sum(isnan(sameStateR)),1);
            % % %                 end
            % % %                 tmpR(sameInp,kAct)=sameStateR;
            % % %             end
            % % %         end
            % % %         % --------------------------------------------------------------
            
            R=tmpR;
            % % %         R(isnan(R))=-0.001;
        end
        % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
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
        y1=y; % $$ For debugging
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
            % maximum Q-value for each poststate
            yymax = max(yy, [], 2);
        end
        
        % calculate discounted future reward for each state transition in batch
        if w.runEnd==1 || s.lp.gamma == 0
            y = R;
        else
            if s.fl.mo==1
                for kAct=1:s.act.numA
                    % $$$ HERE NEED TO CHANGE TO NOT ISNAN I think
                    % % %                 y(A==kAct,kAct) = y(A==kAct,kAct) + ...
                    % % %                     alpha.* ...
                    % % %                     (R(A==kAct,kAct) + s.lp.gamma.* yymax(A==kAct) - y(A==kAct,kAct) );
                    y(~isnan(R(:,kAct)),kAct) = y(~isnan(R(:,kAct)),kAct) + ...
                        alpha.* ...
                        (R(~isnan(R(:,kAct)),kAct) + s.lp.gamma.* yymax(~isnan(R(:,kAct))) - y(~isnan(R(:,kAct)),kAct) );
                end
            else
                y = y + alpha.*(R + s.lp.gamma.* yymax -y );
            end
        end
        
        % % %     for kMem = 1:size(x,1)
        % % %         % only change one action, other Q-values stay the same
        % % %         if w.runEnd==1 || s.lp.gamma == 0
        % % %         else
        % % %             % $$$$ MAYBE I need to change the reward function
        % % %             %                     y(kMem) = reward + s.lp.gamma * yymax(kMem);
        % % %             y(kMem) =y(kMem) + alpha*(R(kMem) + s.lp.gamma * yymax(kMem) -y(kMem) );
        % % %             % %                     Q(prvRow,prvCol,prvgoalC,prvgoalC,action) = Q(prvRow) + alpha*(rewardVal+s.lp.gamma*max(Q(row,:)) - Q(prvRow));
        % % %         end
        % % %     end
        
        % train the network (copied from nntrain())
        if s.fl.mo==1
            net=train(net,x',y');
            y2 = net(x')'; % $$ For debugging
        else
            net=train(net,[x A]',y');
            y2 = net([x A]')';
        end
        
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
        
%         warning('For now, to plot the learning table, disable all special options')
        try
            for iAct=1:length(A)

                % $$ NEED TO GET THE INDEXING RIGHT HERE
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
    
    %     if rem(kBatch,1)==0
    %         save('C:\Users\Rory\Documents\DPPS\DefenseAgent\Scripts\rlsimplepps\testWorkSpace_History_MultSpeedsDirs_VisPropr_V3')
    %     end
    
    if s.fl.perfWhileLearn == 1 && mod(kBatch-1,s.prf.skipBatches) == 0
        [perf.rewPerAct(kBatch,:),perf.rewPerActSD(kBatch,:),perf.rewSums(:,:,kBatch),~] = CalcNetPerf(s,net,extraInput);
    end
    
end
toc

end


