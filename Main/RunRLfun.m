function [s, w, storedEps, net, Qtable] = RunRLfun(s,net,Qtable,storedEps,extraInput)
% [s, w, storedEps, net, Qtable] = RunRLfun(s,net,Qtable,storedEps,extraInput)
% 
% Runs an environment in which an agent is rewarded (or punished) for
% coming into contact with environmental objects
% 
% -------------------------------------------------------------------------
% INPUTS:
% ----------
% s:          settings
% net:        neural network to be used and trained (can be pre-trained)
% Qtable:     [OPTIONAL] Provides q values in tabular format. Only used if
%             appropriate setting is chosen in 's' input
% storedEps:  [OPTIONAL] History of actions and results from a previous run
%             Will be    used to train the network
% extraInput: [OPTIONAL] Stores the value of additional environmental objects
%
% -------------------------------------------------------------------------
% OUTPUTS:
% ----------
% s:          settings
% w:          world configuration at end of run
% storedEps:  History of actions taken and results encountered during run
%             Will be used to train the network. 
% net:        Neural network that was used and trained
% Qtable:     Table of learned q values, relating states and actions
%
%

global world2D;
global tempworld2D;

% =========================================================================
% Complete the settings using defaults
s=DefaultSettings(s);
% If the body can move, add body movement actions: add  1 'action' for 
% each possible combination of action (not the most realisitic but is simple)
if s.fl.bdyMov==1
    if length(s.act.Name)==3
        tempActs2=[];
        for iAct=1:length(s.act.Name)
            tempActs{1}=['lmb' s.act.Name{iAct} '_bdyLEFT'];
            tempActs{2}=['lmb' s.act.Name{iAct} '_bdySTAY'];
            tempActs{3}=['lmb' s.act.Name{iAct} '_bdyRIGHT'];
            tempActs2=[tempActs2 tempActs];
        end
        s.act.Name=tempActs2;
        clear tempActs tempActs2
    elseif length(s.act.Name)==9
    else
        error('Wrong action names supplied!')
    end
end


% =========================================================================
% Initialise the algorithm

% -------------------------------------------------------------------------
% Set the number of actions
s.act.numA=length(s.act.Name); % number of possible actions

% -------------------------------------------------------------------------
% Initialise world
[allCol,allgoalR,allgoalC]=meshgrid(1:s.wrld.size(2),1:s.wrld.size(1),1:s.wrld.size(2));
world2D=50*ones(s.wrld.size);
world2D(1,:)=0;world2D(end,:)=0; world2D(:,1)=0; world2D(:,end)=0;

% -------------------------------------------------------------------------
% Set threat and goal position and create threat
w.thr.row   =round(s.wrld.size(1)/2); w.thr.col     =round(s.wrld.size(2)/2);
w.goal.row  =round(s.wrld.size(1)/2); w.goal.col    =round(s.wrld.size(2)/2);

% -------------------------------------------------------------------------
% Initialise the row and column to plot/put the limb and body in
if isempty(s.lmb.startRow)
    w.lmb.row=size(world2D,1)-2;
else
    w.lmb.row = s.lmb.startRow;
end
if isempty(s.lmb.startCol)
    w.lmb.col=2;
else
    w.lmb.col = s.lmb.startCol;
end

if isempty(s.bdy.startRow)
    w.bdy.row=size(world2D,1)-1;
else
    w.bdy.row = s.bdy.startRow;
end
if isempty(s.bdy.startCol)
    w.bdy.col=2;
else
    w.bdy.col = s.bdy.startCol;
end
startR = w.lmb.row;
startC = w.lmb.col;
startBdyRow=w.bdy.row;
startBdyCol=w.bdy.col;
w.world2D=world2D;
rFl=0;

% -------------------------------------------------------------------------
% Initialise tool
w.lmb.ToolRows=s.lmb.ToolRows;
if max(w.lmb.ToolRows~=0)
    w.lmb.ToolPresent=1;
else
    w.lmb.ToolPresent=0;
end

% -------------------------------------------------------------------------
% Initialise the holding status
% This indicates whether the limb is holding a reward, specifically for the
% environment in which the agent needs to 'eat' to receive reward.
w.lmb.holdingReward = 0; 

% -------------------------------------------------------------------------
% Set speeds of 'stimuli' - sample randomly,
thrSpeedR=randsample(s.thr.alSpR,1);
thrSpeedC=randsample(s.thr.alSpR,1);
goalSpeedR=randsample(s.gol.alSpR,1);
goalSpeedC=randsample(s.gol.alSpC,1);

% -------------------------------------------------------------------------
% Start video
if s.fl.vid==1
    v = VideoWriter(s.plt.vidFileName);
    v.FrameRate=s.plt.vidFR;
    v.Quality = 100;
    open(v)
end

% -------------------------------------------------------------------------
% Check for problems in settings
if s.lp.b1Siz>s.lp.minExp
    error('The First batch should be equal to or larger than the minimum history size!')
end

% -------------------------------------------------------------------------
%Initialise stored states, actions and rewards
if ~exist('storedEps','var') | isempty(storedEps)
    % The number of actions already taken on previous runs
    startAction=0;
    
    storedEps.prvS = cell([s.rp.numIter*s.rp.maxActions,1]);
    storedEps.S = cell([s.rp.numIter*s.rp.maxActions,1]);
    storedEps.A = nan([s.rp.numIter*s.rp.maxActions,1]);
    storedEps.R = nan([s.rp.numIter*s.rp.maxActions,1]);
else
    % The number of actions already taken on previous runs
    startAction=length(storedEps.A);
    
    storedEps.prvS = [ storedEps.prvS ; cell([s.rp.numIter*s.rp.maxActions,1])  ];
    storedEps.S = [ storedEps.S ; cell([s.rp.numIter*s.rp.maxActions,1])  ];
    storedEps.A = [ storedEps.A ; nan([s.rp.numIter*s.rp.maxActions,1])  ];
    storedEps.R = [ storedEps.R ; nan([s.rp.numIter*s.rp.maxActions,1])  ];
end


% -------------------------------------------------------------------------
% make some copies of world to use later for display
orgworld2D = world2D;
orgworld2D(w.lmb.row,w.lmb.col) = 50;
tempworld2D = orgworld2D;

% -------------------------------------------------------------------------
% Initialise extra input if necessary
if s.xtra.useThr==1 && s.fl.extraInput==1
    w.xtra.currV = squeeze(extraInput(w.lmb.row,w.lmb.col,w.thr.row,w.thr.col,:))'; %Current value
elseif s.fl.extraInput==1
    w.xtra.currV = squeeze(extraInput(w.lmb.row,w.lmb.col,w.goal.row,w.goal.col,:))'; %Current value
end


% -------------------------------------------------------------------------
% initialise action count properly
countActions = 0;

% -------------------------------------------------------------------------
% Initialise touch value (start not touched)
s.lmb.size=1;
w.touchV=zeros([1 s.bdy.size+s.lmb.size]); % This is actually the body touch
w.rtTask.percieveTouch=0; % The touch value for the RT task

% -------------------------------------------------------------------------
% build a state action matrix by finding all valid states from world
% we have s.act.numA actions for each state (move left, move right). So the limb
% can only move on 1 dimension, meaning that the 2nd dimension is the
% position of the stimulus. For now I will only give x position
if s.fl.newTable==1
    Q = zeros(size(world2D,1),size(world2D,2),size(world2D,1),size(world2D,2),s.act.numA);
    Qtable=Q;
%     firstPass=1;
else
%     firstPass=0;
end
if ~exist('Qtable','var') && s.fl.newNet==0
    error('No Q table supplied, and no new Qtable requested')
end

% -------------------------------------------------------------------------
% Set status message for goal and bump
BUMP = s.act.BUMP ;
GOAL = s.act.GOAL;
THREAT = s.act.THREAT;

% -------------------------------------------------------------------------
% % % Initialise videos
if  s.fl.dspm == 1
    screen = get(0,'ScreenSize');
    hworld=figure('Units', 'pixels','Position',screen*2);

    set(hworld, 'Renderer', 'painters');
    opengl hardware
end

% -------------------------------------------------------------------------
% Run the first State update
S=UpdateS(s, w, tempworld2D); 

% -------------------------------------------------------------------------
% Expand the state S if 1 frame of history is being stored
if s.fl.hist==1
    S=repmat(S,[1 2]);
end

% -------------------------------------------------------------------------
% Initialise net
if s.fl.newNet==1
    if s.fl.deepNet==1
        net.Layers = [
            imageInputLayer(size(S'),'Normalization','none')
            fullyConnectedLayer(200,"Name","fc_1")
            fullyConnectedLayer(100,"Name","fc_2")
            fullyConnectedLayer(50,"Name","fc_3")
            fullyConnectedLayer(s.act.numA,"Name","fc_4")
            regressionLayer("Name","regressionoutput")];
        
        firstTrainNum=s.lp.n1stTr;
        
        % ---------------------------------------------------------------------
    else
        net=feedforwardnet(s.lp.netS);
        net=ChangeNeuronType(net,s);
        net=init(net);
        net.trainFcn = s.lp.TrnFn;
        net.trainParam.epochs=s.lp.TrnEp;
        net.trainParam.showWindow = s.lp.ShwW;
        net.trainParam.max_fail=s.lp.mxF;
        firstTrainNum=s.lp.n1stTr; % Number of gradient descents to do 1st time round
    end
    
    if s.fl.newTable==1
        firstPass=1; % this is to initialise the network properly before the first batch
    elseif s.fl.newTable==0
        firstPass=0;
    end
else
    firstPass=0;
end
if ~exist('net','var') && s.fl.newNet==0
    error('No network supplied, and no new network requested')
end
w.runEnd=0; % This is to be used for a game with separate runs

%%
for i=1:s.rp.numIter
% % %     tempworld2D(w.thr.row,w.thr.col) = 10;
% % %     tempworld2D(w.goal.row,w.goal.col) = 100;
    w.lmb.row = startR; w.lmb.col = startC;
    w.bdy.col=startBdyCol; w.bdy.row=startBdyRow;
    status = -1;
    actionsWithinRun =0;
    
    while actionsWithinRun < s.rp.maxActions %status ~= GOAL
        
        % -----------------------------------------------------------------
        % record the current state variables for use later
        w.lmb.prvRow  = w.lmb.row;  w.lmb.prvCol  = w.lmb.col;
        w.goal.prvRow = w.goal.row; w.goal.prvCol = w.goal.col;
        w.thr.prvRow  = w.thr.row;  w.thr.prvCol  = w.thr.col;
        w.lmb.prvToolRows=w.lmb.ToolRows;
        w.lmb.PrvToolPresent=w.lmb.ToolPresent;
        w.lmb.prvHoldingReward = w.lmb.holdingReward;

        w.bdy.prvCol=w.bdy.col;
        prvtouchV = w.touchV;
        w.rtTask.prvTouch= w.rtTask.percieveTouch;
        if s.fl.extraInput==1, w.xtra.prvV = w.xtra.currV; end

        
        % -----------------------------------------------------------------
        % If the stimulus was reset, change the tool length/existence
        if s.fl.ToolChange==1 & rFl==1
            if rand(1)<s.lmb.ToolProb
                w.lmb.ToolRows=s.lmb.ToolRows;
                w.lmb.ToolPresent=1;
            else
                % First clear the tool from the world map, then clear the tool
                tempworld2D(w.lmb.row-w.lmb.ToolRows,w.lmb.col) = 50;
                w.lmb.ToolRows=[0];
                w.lmb.ToolPresent=0;
            end
        end
        
        % -----------------------------------------------------------------
        % Store the previous state
        prvS=S;
        
        % -----------------------------------------------------------------
        % select an action value i.e. Direction
        % which has the maximum value of Q in it
        % if more than one actions has same value then select randomly from them
        if firstPass==1
            action=randperm(s.act.numA);
            action=action(1);
            bestAction=action;
        else
            clear Qest;
            % If using network, calculate action values using  net
            if strcmp(s.act.ValAprFun,'Net')
                if s.fl.mo==1
                    if s.fl.deepNet==1
                        Qest = predict(net,reshape(prvS',[size(prvS,2) 1 1 size(prvS,1)]));
                    else
                        Qest=net(prvS');
                    end
                else
                    for kAct=1:s.act.numA
                        Qest(kAct)=net([prvS  kAct]');
                    end
                end
                % else caulculate action values using the table
            elseif strcmp(s.act.ValAprFun,'Table')
                for kAct=1:s.act.numA
                    Qest(kAct)=Qtable(w.lmb.prvRow,w.lmb.prvCol,w.goal.prvRow,w.goal.prvCol,kAct);
                end
            else
                warning('The function to calculate value should be "Net" or "Table" ')
            end
            
            [val,index] = max(Qest);
            indexAll = find(Qest == val);
            action = indexAll(1);
            % If more than one action has optimal value, pick randomly
            if size(indexAll,1) > 1
                index = 1+round(rand*(length(indexAll)-1));
                action = indexAll(index,1);
            end
            try
                bestAction=action;
            catch
                index = 1+round(rand*(length(indexAll)-1));
                action = indexAll(1,index);
                bestAction=action;
            end
            % Select action dpending on policy
            action = s.act.pi(action,s.act.numA,rand());

        end
        
        % -----------------------------------------------------------------
        % display the world after some steps
        if rem(countActions,s.lp.dispW) == 0 & s.fl.dspm == 1 & exist('action')
            DisplayWorld(s, w, action, hworld,storedEps,countActions)
            % Make it possible to record video
            if s.fl.vid==1
                frame=getframe(gcf);
                writeVideo(v,frame);
            end
        end
        
        % -----------------------------------------------------------------
        % Move the Threat [TBD: make a function to move both threats and goals]
        if s.fl.thr==1
            [w.thr.row,w.thr.col,rFl,bodyTouch,w.teleportLR]=MoveObj(s,w.bdy.col,w.thr.row,w.thr.col,thrSpeedR,thrSpeedC, s.thr.randSpr);
            w.touchV(s.lmb.size+1 : end)=bodyTouch;
            bodyStatus=-bodyTouch;
            % If the threat is being reset, give it a different speed
            if rFl==1
                thrSpeedR=randsample(s.thr.alSpR,1);
                thrSpeedC=randsample(s.thr.alSpC,1);
            end
            UpdateWorldMap(w,w.thr.prvRow,w.thr.prvCol,w.thr.row,w.thr.col,rFl,10)
        end
        
        % -----------------------------------------------------------------
        % Move the Goal
        [w.goal.row,w.goal.col,rFl,bodyTouch,w.teleportLR]=MoveObj(s,w.bdy.col,w.goal.row,w.goal.col,goalSpeedR,goalSpeedC, s.gol.randSpr);
        w.touchV(s.lmb.size+1 : end)=bodyTouch;
        if s.fl.onlyThr==1
            bodyStatus=-bodyTouch;
        else
            bodyStatus=bodyTouch;
        end
        % If the goal is being reset, give it a different speed
        if rFl==1
            goalSpeedR=randsample(s.gol.alSpR,1); % [TBD: put this and other instances into the reset obj position function]
            goalSpeedC=randsample(s.gol.alSpC,1);
        end
        UpdateWorldMap(w,w.goal.prvRow,w.goal.prvCol,w.goal.row,w.goal.col,rFl,100)
        
        % -------------------------------------------------------------------------
        % Update extra input if necessary
        if s.xtra.useThr==1 && s.fl.extraInput==1
            w.xtra.currV = squeeze(extraInput(w.lmb.row,w.lmb.col,w.thr.row,w.thr.col,:))'; %Current value
        elseif s.fl.extraInput==1
            w.xtra.currV = squeeze(extraInput(w.lmb.row,w.lmb.col,w.goal.row,w.goal.col,:))'; %Current value
        end
        
        
        % -----------------------------------------------------------------
        % Perform the selected action i.e. MoveLeft or MoveRight (the 'mod'
        % is in case the body is also used for actions
        if s.fl.bdyMov==1 
            %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
            w.mvBdy=1;
            lmbAction=ceil(action/sqrt(length(s.act.Name)));
            [w.bdy.row,w.bdy.col,statusBdy] = MoveLR(s,w,w.bdy.row,w.bdy.col,action);
            %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
        else
            lmbAction=action;
        end
        % If there is a tool, find out what its status is. Do it before
        % updating the limbCol of course
        clear toolStatus % <-- This is because the length of the tool could change
        if max(w.lmb.ToolRows~=0) 
            for iToolRow=1:length(w.lmb.ToolRows)
                %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
                w.mvBdy = 0;
                [dmyR,dmC,toolStatus(iToolRow)] = ...
                    MoveLR(s,w,w.lmb.row-w.lmb.ToolRows(iToolRow),w.lmb.col,action);
                %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
            end
        end
        %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬
        w.mvBdy=0;
        [w.lmb.row,w.lmb.col,status] = MoveLR(s,w,w.lmb.row,w.lmb.col,action);
        %¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬

        % -----------------------------------------------------------------
        % IF the reaction time task is active, update the touch values
        if mod(countActions,30)>0 & mod(countActions,30)<5
            w.rtTask.realTouch=s.rtt.threshold;
            w.rtTask.percieveTouch=s.rtt.threshold+s.rtt.NoiseFun(s.rtt.mu,s.rtt.std);
        else
            w.rtTask.realTouch=0;
            w.rtTask.percieveTouch=w.rtTask.realTouch+s.rtt.NoiseFun(s.rtt.mu,s.rtt.std);
        end
        
        
        
        % -----------------------------------------------------------------
        % If the limb/tool/body intercepts an object, reset the position
        [w, rFl] = TouchObjReset(status,GOAL, THREAT, rFl, s,w);
        if s.fl.bdyMov==1, [w, rFl] = TouchObjReset(statusBdy,GOAL, THREAT, rFl, s,w); end
        if max(w.lmb.ToolRows~=0), for iToolRow=1:length(w.lmb.ToolRows)
                [w, rFl] = TouchObjReset(toolStatus(iToolRow),GOAL, THREAT, rFl, s,w);
        end;end
        
        % -----------------------------------------------------------------
        % count the actions required to reach the goal
        countActions = countActions + 1;
        
        
        % -----------------------------------------------------------------
        % Get the reward values for the base task
        if s.fl.bdyMov==1
            allRewards=s.act; 
            rewardVal = StatusConsequence(s,status,ismember(action,[4 5 6]),GOAL,BUMP,THREAT,allRewards);
            % Add the reward values for events regarding the body
            clear allRewards
            allRewards.BumpRew=s.act.bdyBumpRew;
            allRewards.GoalRew=s.act.bdyGoalRew;
            allRewards.ThreatRew=s.act.bdyThreatRew;
            allRewards.DefaultMoveRew=s.act.bdyMoveRew;
            allRewards.DefaultStayRew=s.act.bdyStayRew;
            rewardVal = rewardVal + StatusConsequence(s,statusBdy,ismember(action,[2 5 8]),GOAL,BUMP,THREAT,allRewards);
        else
            allRewards=s.act;
            rewardVal = StatusConsequence(s,status,ismember(action,2),GOAL,BUMP,THREAT,allRewards);
        end
        
        % -----------------------------------------------------------------
        % If there is a tool, add the rewards for tool movement/touch
        if max(w.lmb.ToolRows~=0), for iToolRow=1:length(w.lmb.ToolRows)
            % Assume identical rewards for tool and for limb contact
            allRewards=s.act;
            % Assume that moving the stick doesn't cost any extra
            % effort over moving the limb (i.e. set stationaryFlag = 1)
            rewardVal = rewardVal + StatusConsequence(s,toolStatus(iToolRow),1,GOAL,BUMP,THREAT,allRewards);
            end;end
        
        % -----------------------------------------------------------------
        % If the task is a reaction-time task to tactile stimuli, add the relevant rewards
        if s.fl.rtTask==1
            if contains(s.act.Name{action},'BUTTON')
                if w.rtTask.realTouch>=s.rtt.threshold
                    rewardVal=rewardVal + s.act.rtRew;
                else
                    rewardVal=rewardVal - s.act.rtRew;
                end
            end
        end

        % -----------------------------------------------------------------
        % If the task involves moving 'caught food' to the 'mouth', add the 
        % relevant rewards
        if s.fl.grabAndEat == 1
            % Only if the agent is actually holding something, and if the
            % limb is in the correct position
            if w.lmb.holdingReward == 1 && w.lmb.col == ceil(s.wrld.size(2)./2)
                w.lmb.holdingReward = 0; % Reward is eaten, so hand is empty
                rewardVal = rewardVal + s.act.eatRew;
            end
        end

        % -----------------------------------------------------------------
        % If the task involves defending a particular area, add the
        % relevant rewards
        if s.fl.defendZone == 1
            if w.goal.row == 12 & w.goal.col == ceil(s.wrld.size(2)./2)
                rewardVal = rewardVal + s.act.shieldRew;
            end
        end
        
        % -----------------------------------------------------------------
        % Other consequences of touch $$ might have to think more here $$
        % For now, this is to allow for touch of the body IF the body isn't
        % moving
        if s.fl.bdyMov==0
            rewardVal=rewardVal+4.*sum(bodyStatus);
        end
        bodyStatus=0;
        
        
       
        % -----------------------------------------------------------------
        % Update S --------------------------------------------------------
        if s.fl.hist==1
            presStart=numel(S)/2 +1; % start of present values
            S(1:end/2)=S(presStart:end);
        else
            presStart=1;
        end
        S(presStart:end)=UpdateS(s, w, tempworld2D);
        
        % -----------------------------------------------------------------
        % Store state transitions, actions and rewards
        % prev state, next state, action, reward
        storedEps.prvS{startAction+countActions,1}=prvS;
        storedEps.S{startAction+countActions,1}=S;
        storedEps.A(startAction+countActions,1)=action;
        storedEps.R(startAction+countActions,1)=rewardVal;
        
        % -----------------------------------------------------------------
        % Show the number of actions taken and save the workspace
        if rem(countActions,s.lp.dispA)==0
            countActions
        end
        
        
        % -----------------------------------------------------------------
        % Perform Q-learning
        if countActions>=s.lp.minExp & rem(countActions,s.lp.retr)==0
            
            % Possibly only use the events from the s.lp.rcSiz epochs
            if s.fl.rdat==1
                eligibleEps=find(~isnan(storedEps.A(:,1))' &  ...
                    [1:numel(storedEps.A(:,1))] > (startAction+countActions-s.lp.rcSiz) );
                if firstPass==1
                    batchSize=s.lp.b1Siz;
                else
                    batchSize=s.lp.bSiz;
                end
                bTrls=datasample(eligibleEps,batchSize);
                batchEps.prvS = storedEps.prvS(bTrls,: );
                batchEps.S = storedEps.S(bTrls,: );
                batchEps.A = storedEps.A(bTrls,: );
                batchEps.R = storedEps.R(bTrls,: );
                [net Qtable] = RelearnNNFunSC(s,w,batchEps,net,Qtable);
                
            % OR just learn using the full amount of stored epochs
            else
                [net Qtable] = RelearnNNFunSC(s,w,storedEps,net,Qtable);
            end
            s.fl.newNet=0;
            s.fl.newTable=0;
        end
        
        actionsWithinRun=actionsWithinRun+1;
    end
    
    iterationCount(i,1) = countActions;
    
    
end

w.world2D=world2D;
w.tempworld2D=tempworld2D;

% Ensure that the important parts of the stored states are easy to access
[w] = StoreStateAttributes(storedEps,w);

end

function [w, rFl] = TouchObjReset(status,GOAL, THREAT, rFl, s,w)
if status == GOAL | status == THREAT
    % Update the touch value
    w.touchV(s.lmb.size)=1; % [TBD: This needs to be changed if the limb ends up being bigger]
    % Put the threat/goal back up to a starting position
    if status == GOAL
        [w.goal.row w.goal.col, rFl, goalSpeedR, goalSpeedC]=ResetObjPos(s,s.gol.alSpR, s.gol.alSpC);
        UpdateWorldMap(w,w.goal.prvRow,w.goal.prvCol,w.goal.row,w.goal.col,rFl,100);

        % Set 'grabbed status', in case the food-to-mouth scenario applies
        w.lmb.holdingReward = 1;
    elseif status == THREAT
        [w.thr.row w.thr.col, rFl, thrSpeedR, thrSpeedC]=ResetObjPos(s,s.thr.alSpR, s.thr.alSpC);
        UpdateWorldMap(w,w.thr.prvRow,w.thr.prvCol,w.thr.row,w.thr.col,rFl,10);
    end
else
    % otherwise there was no touch
    w.touchV(s.lmb.size)=0;
end
end

%%
function DisplayWorld(s, w, action, hworld, storedEps, countActions)

global tempworld2D

if s.fl.rtTask==1
    subplot(2,1,1)
end

set(0,'CurrentFigure',hworld)

if s.fl.rtTask==1
    subplot(2,1,1);
end


if s.plt.fancyWorld == 1
    plotWrld                                = zeros(size(w.world2D) - 1);
    plotWrld(w.goal.row - 1,w.goal.col - 1) = - 1;
    if s.fl.thr == 1
        plotWrld(w.thr.row - 1,w.thr.col - 1)   =   1;
    end
    imagesc(plotWrld);
    colormap(redbluecmapRory(100,20));
    hold on
    
    if s.plt.grid == 1    
        fS.gridXstart   = -4.5;
        fS.gridXstep    = 1;
        fS.gridYstart   = 0.5;
        fS.gridYstep    = 1;
        fS.color        = [.8 .8 .8];
        GridOverImage(fS,gca);
    end

    if s.fl.bdyMov == 1
        plot(w.bdy.col - 1, w.bdy.row - 1, '^' ,'MarkerFaceColor','w', ...
            'MarkerEdgeColor','k', 'MarkerSize',20,'LineWidth',5);
    end

    % Plot hand
    tmpFig = imread('F:\Projects\DPPS\DefenseAgent\Documentation\Figures\Finalising_All\OnlyHand.png');
    alphachannel = 1-double(all(tmpFig == 255, 3));
    hold on
    h = imshow(tmpFig,'XData',w.lmb.col + [-1.5 -.5],'YData',w.lmb.row+ [-1.5 -.5]);
    set(h, 'AlphaData', alphachannel);

    xlim([.5 s.wrld.size(2)-.5])
    caxis([-1.5 1.5])
    axis off
    hold off
else
    imagesc(tempworld2D);%,colorbar;
end

bA=s.act.Name{action};
title(['Best Action: ' bA]);
drawnow

if s.fl.rtTask==1
    subplot(2,1,2)
    if mod(countActions,10)==0 & sum(~isnan(storedEps.R))>100
        tmpS=cell2mat(storedEps.S);
        plotEnd=sum(~isnan(storedEps.R));
        if s.fl.extraInput==1
            plot(tmpS(plotEnd-100:plotEnd,end-length(w.xtra.currV)));
        else
            plot(tmpS(plotEnd-100:plotEnd,end));
        end
        hold on
        plot(storedEps.A(plotEnd-100:plotEnd)==4 & storedEps.R(plotEnd-100:plotEnd)>0);
        title('touch perception');
        hold off
        drawnow
    end
end

end

%%
function rewardVal = StatusConsequence(s,status,stationaryFl,GOAL,BUMP,THREAT,allRewards)


    if status == BUMP
        rewardVal = allRewards.BumpRew;
    elseif status == GOAL
        if s.fl.onlyThr==1
            rewardVal = -allRewards.GoalRew;
        else
            rewardVal = allRewards.GoalRew;
        end
    elseif status == THREAT
        rewardVal = allRewards.ThreatRew;
    else
        if stationaryFl == 1
            rewardVal = allRewards.DefaultStayRew;
        else
            rewardVal = allRewards.DefaultMoveRew;
        end
    end

end
