%
% -----
% USE:
% 
% Run this script to generate the data from scratch. 
% -----
%
% =========================================================================
%                              !!! WARNING !!!
% =========================================================================
% It is !!! HIGHLY RECOMMENDED !!! to either: 
%       --------------------------   
% A) Just load the existing outputs from this script, which are available
%    in the 'GeneratedData' sub-directory
% OR      -------------
% B) Run each section separately:
%    some sections can take a VERY VERY long time to run
%    (Yes,this script is highly un-optimised in terms of runtime)
% =========================================================================
%
%

[basePath addedPaths] = SetPathEgocentricMaps();


%% $$$ HERE HERE
foldName(1:end-10)



s.rp.maxActions = 2e6 ;
s.lp.minExp = 2e7;


%% Data for Model 1: Artificial agents produce response fields anchored to the limb
clear s rS ntRS rtRS ytRS olRS

disp('Here 1')

% THIS CREATES THE FILE BELOW:
% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig1_Results_v3.mat')

% -------------------------------------------------------------------------
% Base settings
s.gol.alSpR = 1;
s.gol.alSpC = [0 0];
s.gol.randSpr = [0 0];
s.rp.maxActions = 2e6 ;
s.lp.minExp = 2e7;
% s.rp.maxActions = 2e3 ;
% s.lp.minExp = 2e4;
s.lp.b1Siz = 1e4;
s.wrld.size = [14 15];
s.wrld.resetType = 'BottomTop_InfLR';
s.lp.epsilon = 0.3;
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=0;
% s.rl.maxRetr=100;
s.rl.maxRetr=200;

% -------------------------------------------------------------------------
% Settings for different plots
gammas=[0.7 0.3 0.7 0.3 ...
    0.7 0.3 0.7 0.3 ];
actionNames={{'STAY'},{'STAY'},{'LEFT','STAY','RIGHT'},{'LEFT','STAY','RIGHT'} ...
    {'STAY'},{'STAY'},{'LEFT','STAY','RIGHT'},{'LEFT','STAY','RIGHT'}};
goalDynamics=[[0 0] ; [0 0] ; [0 0] ; [0 0] ; ...
    [0 3] ; [0 3] ; [0 3] ; [0 3]] ;

% -------------------------------------------------------------------------
% Run each type of model
for iM=1:length(gammas)
    tic
    % special case, to improve accuracy of estimates when necessary
    if iM == 6 | iM==5
        tmps=s;
        s.lp.netS=[12 12 12 12];
        s.fl.trainNet=1;
        s.fl.newNet=1;
        s.fl.newTable=1;
        s.rl.maxRetr=200;
    elseif iM == 7
        s=tmps;
    end
    
    % Change specific settings of the model
    s.lp.gamma = gammas(iM);
    s.act.Name = actionNames{iM};
    s.gol.randSpr = goalDynamics(iM,:);
    % Run the model and learn Q (table)
    [s w storedEps dmyNet Qtable] = RunRLfun(s);
    [net Qtable1] = RelearnNNFun(s,w,storedEps,dmyNet,Qtable);
    % Store results Structure
    rS(iM).s=s; rS(iM).Qtable=Qtable1; rS(iM).w=w; rS(iM).net=net;
    toc
    iM
end

svNm=[foldName 'Z\ForFigures\Fig1_Results.mat'];
try
    save(svNm,'rS','-v7.3')
catch
    save(svNm,'rS')
end

%% Find performance of different network sizes for proximity and position dependence
clear s rS ntRS rtRS ytRS olRS

disp('Here 2')

% -------------------------------------------------------------------------
% Base settings
s.gol.alSpR = 1;
s.gol.alSpC = [0 0];
s.gol.randSpr = [2 2];
s.lp.b1Siz = 1e4;
s.wrld.size = [14 15];
s.wrld.resetType = 'BottomTop_InfLR';
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr=50;
% % s.rl.maxRetr=2;

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=2;

netSizes={[9],[9 9],[9 9 9],[9 9 9 9],[9 9 9 9 9]};

% Run the model and learn Q (table)
[s w storedEps net Qtable] = RunRLfun(s);

% Run each type of model
for iM = 1:length(netSizes)
    tic
    
    s.fl.newNet = 1;
    s.fl.newTable = 1;
    
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};
    
    % Relearn Q with different netsizes
    [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
    % Store network results Structure
    ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
    ntRS(iM).perf = perf;
    
    bFld = 'F:\Projects\DPPS\DefenseAgent\Results\Performance\ProximityPosition\NetSizes\';
    svNm = [bFld 'Mults_of_9_v3'];
    save(svNm,'ntRS','-v7.3');
    
    iM
    toc
end

%% Find performance of different network sizes for MULTIPLE LIMBS and proximity position dependence
clear s rS ntRS rtRS ytRS olRS

disp('Here 3')

% -------------------------------------------------------------------------
% Base settings
s.gol.alSpR = 1;
s.gol.alSpC = [0 0];
s.gol.randSpr = [2 2];
s.lp.b1Siz = 1e4;
s.lp.minExp = 4e7; % minumum number of actions before training sarts (big here so no training)
s.lp.dispA = 1e6; % only show action count very infreequently
% % % s.wrld.size = [14 15];
s.wrld.size = [10 15]; % $$$ DOING THIS ONE for simplicity
s.wrld.resetType = 'BottomTop_InfLR';
s.rp.maxActions = 4e6 ;
s.act.bdyGoalRew = 2;
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr=100;
s.fl.bdyMov = 1;

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=5;

netSizes={[9],[9 9],[9 9 9],[9 9 9 9],[9 9 9 9 9],[9 9 9 9 9], ...
    9*ones(1,6),9*ones(1,7),9*ones(1,8),9*ones(1,9),9*ones(1,10),9*ones(1,11)};

% -------------------------------------------------------------------------
% Run the model
[s w storedEps net Qtable] = RunRLfun(s);
% -------------------------------------------------------------------------
% Run each type of model
for iM = 1:length(netSizes)
    
    s.fl.newNet = 1;
    s.fl.newTable = 1;
    
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};

    % Relearn Q with different netsizes
    [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
    % Store network results Structure
    ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
    ntRS(iM).perf = perf;
    
    bFld = 'F:\Projects\DPPS\DefenseAgent\Results\Performance\ProximityPosition\NetSizes\';
    svNm = [bFld '2LIMBS_Mults_of_9_V4'];
    save(svNm,'ntRS','-v7.3');
    
end

%% Data for plot $2$ about distance and position dependence, and corresponding stats
clear s rS ntRS rtRS ytRS olRS

disp('Here 4')

% -------------------------------------------------------------------------
% Base settings
s.gol.alSpR = 1;
s.gol.alSpC = [0 0];
s.gol.randSpr = [2 2];
s.lp.b1Siz = 1e4;
s.wrld.size = [14 15];
s.wrld.resetType = 'BottomTop_InfLR';
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr=100;

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=5;

netSizes={[12 10 8 6],[9 9 9 9],[6 8 10 12]};

% -------------------------------------------------------------------------
% Run the model
 [s w storedEps net Qtable] = RunRLfun(s);
% -------------------------------------------------------------------------
% Run each type of model
for iM=1:length(netSizes)
    
    s.fl.newNet = 1;
    s.fl.newTable = 1;
    
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};
    
    
    [net Qtable] = RelearnNNFun(s,w,storedEps,net,Qtable);
    % Store results Structure
    rS(iM).s=s; rS(iM).Qtable=Qtable; rS(iM).w=w; rS(iM).net=net;
    
    svNm='F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_v5_randRC.mat';
    try
        save(svNm,'rS','-v7.3')
    catch
        save(svNm,'rS')
    end
    iM
end

%% Data for plot $2$, but with MULTIPLE LIMBS, about distance and position dependence, and corresponding stats
clear s rS ntRS rtRS ytRS olRS

disp('Here 5')

% -------------------------------------------------------------------------
% Base settings
s.gol.alSpR = 1;
s.gol.alSpC = [0 0];
s.gol.randSpr = [2 2];
s.lp.b1Siz = 1e4;
s.lp.minExp = 4e7; % minumum number of actions before training sarts (big here so no training)
s.lp.dispA = 1e6; % only show action count very infreequently
% % % s.wrld.size = [14 15];
s.wrld.size = [10 15]; % $$$ DOING THIS ONE for simplicity
s.wrld.resetType = 'BottomTop_InfLR';
s.rp.maxActions = 4e6 ;
s.act.bdyGoalRew = 2;
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr=100;
s.fl.bdyMov = 1;

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=5;

netSizes={[13 11 9 7 5],[9 9 9 9 9],[5 7 9 11 13]};

% -------------------------------------------------------------------------
% Run the model
[s w storedEps net Qtable] = RunRLfun(s);
% -------------------------------------------------------------------------
% Run each type of model
for iM = 1:length(netSizes)
    
    s.fl.newNet = 1;
    s.fl.newTable = 1;
    
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};

    % Relearn Q with different netsizes
    [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
    % Store network results Structure
    ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
    ntRS(iM).perf = perf;
    
    bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\';
    svNm = [bFld 'Fig_Dist_Pos_Dependence_2LIMBS_v2.mat'];
    save(svNm,'ntRS','-v7.3');
    iM
end






%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% % % % % % % % % % % % %% Run a RT scenario with and without loaded QPPS.
% % % % % % % % % % % % clear s rS ntRS rtRS ytRS olRS
% % % % % % % % % % % % 
% % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % % First load normal Q-values
% % % % % % % % % % % % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_V4.mat');
% % % % % % % % % % % % 
% % % % % % % % % % % % iM=2;
% % % % % % % % % % % % 
% % % % % % % % % % % % s = rS(iM).s;
% % % % % % % % % % % % w = rS(iM).w;
% % % % % % % % % % % % net = rS(iM).net;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.plt.rowLims=[1.5 s.wrld.size(1)-1.5];
% % % % % % % % % % % % s.plt.colLims=[1.5 s.wrld.size(2)-2.5];
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % % s.plt.stimRow = [3:size(w.world2D,1)-1];
% % % % % % % % % % % % % s.plt.stimCol = [2:size(w.world2D,2)-1];
% % % % % % % % % % % % s.plt.lmbCol = 8 %[2:size(w.world2D,2)-1];
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % [QPPS,allNeurAct] = CalcNetOutput(s,w,net);
% % % % % % % % % % % % 
% % % % % % % % % % % % DisplActValsFun(s,w,QPPS);
% % % % % % % % % % % % 
% % % % % % % % % % % % caxis([0 2]);
% % % % % % % % % % % % colormap(whitetocol(100,[0 0 0.7]));
% % % % % % % % % % % % 
% % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % % Then run the reaction time scenario
% % % % % % % % % % % % % $$$ I'm redoing it witht the proper world size
% % % % % % % % % % % % 
% % % % % % % % % % % % % clear s rS ntRS rtRS ytRS olRS net storedEps w
% % % % % % % % % % % % 
% % % % % % % % % % % % clear rS storedEps;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.lp.netS = [9 9];
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.newNet = 1;
% % % % % % % % % % % % s.fl.newTable = 1;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.rp.maxActions = 2e5 ;
% % % % % % % % % % % % s.lp.minExp = 2e6 ;
% % % % % % % % % % % % 
% % % % % % % % % % % % % s.rp.maxActions = 2e6 ;
% % % % % % % % % % % % % s.lp.minExp = 2e7 ;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.rtTask=1;
% % % % % % % % % % % % s.fl.hist=0;
% % % % % % % % % % % % s.act.Name={'LEFT','STAY','RIGHT','BUTTON'};
% % % % % % % % % % % % s.fl.dspm=0;
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % s.act.rtRew=1;
% % % % % % % % % % % % s.rtt.threshold=3;
% % % % % % % % % % % % 
% % % % % % % % % % % % for defaultRew = [0, -0.01]
% % % % % % % % % % % %     
% % % % % % % % % % % %     clear rS
% % % % % % % % % % % %     
% % % % % % % % % % % %     s.act.BumpRew = defaultRew;
% % % % % % % % % % % %     s.act.GoalRew = defaultRew;
% % % % % % % % % % % %     s.act.ThreatRew = defaultRew;
% % % % % % % % % % % %     s.act.DefaultMoveRew = defaultRew;
% % % % % % % % % % % %     s.act.bdyBumpRew = defaultRew;
% % % % % % % % % % % %     s.act.bdyGoalRew = defaultRew;
% % % % % % % % % % % %     s.act.bdyThreatRew = defaultRew;
% % % % % % % % % % % %     s.act.bdyMoveRew = defaultRew;
% % % % % % % % % % % %     
% % % % % % % % % % % %     s=DefaultSettings(s);
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     s.rl.maxRetr=100;
% % % % % % % % % % % %     
% % % % % % % % % % % %     extraInputs = [0 1];
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     s.fl.perfWhileLearn = 1;
% % % % % % % % % % % %     s.prf.nRep = 10;
% % % % % % % % % % % %     s.prf.skipBatches = 49;
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Run each type of model
% % % % % % % % % % % %     for iM = 1:length(extraInputs)
% % % % % % % % % % % %         
% % % % % % % % % % % %         s.fl.newNet = 1;
% % % % % % % % % % % %         s.fl.newTable = 1;
% % % % % % % % % % % %         
% % % % % % % % % % % %         % Change specific settings of the model
% % % % % % % % % % % %         s.fl.extraInput = extraInputs(iM);
% % % % % % % % % % % %         
% % % % % % % % % % % %         % Run model with new settings (extra Q-inputs)
% % % % % % % % % % % %         [s w storedEps net Qtable]=RunRLfun(s,[],[],[],QPPS);
% % % % % % % % % % % %         
% % % % % % % % % % % %         % Relearn Q with different netsizes
% % % % % % % % % % % %         [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable,QPPS);
% % % % % % % % % % % %         % Store network results Structure
% % % % % % % % % % % %         rtRS(iM).s=s; rtRS(iM).Qtable=Qtable; rtRS(iM).w=w; rtRS(iM).net=net;
% % % % % % % % % % % %         rtRS(iM).perf = perf;
% % % % % % % % % % % %         
% % % % % % % % % % % %         bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\';
% % % % % % % % % % % %         svNm = [bFld 'Fig_Dist_Pos_Dependence_RTtask_BaseRew_' num2str(defaultRew) '_noHist_EasierTask_Only1batch_smallNet_QPPS2_v2.mat'];
% % % % % % % % % % % %         save(svNm,'rtRS','QPPS','-v7.3');
% % % % % % % % % % % %         
% % % % % % % % % % % %     end
% % % % % % % % % % % % end
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % % % [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % s.rl.maxRetr=100;
% % % % % % % % % % % % % % % % s.fl.newNet
% % % % % % % % % % % % % % % % for ii=1:100
% % % % % % % % % % % % % % % %     [net Qtable] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % % % % % % % % % % % % % % %     ii
% % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % save('F:\Projects\DPPS\DefenseAgent\Results\Net_RTtask1_NormNoise_HarderTask_ExtraInput_0BaseRew_SimpleQPPS.mat','net','s','w','storedEps');
% % % % % % % % % % % % % % save('F:\Projects\DPPS\DefenseAgent\Results\Net_RTtask1_NormNoise_HarderTask_ExtraInput_0BaseRew_SimpleQPPS.mat','net','s','w','storedEps');
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % % % % AND from here, it is the non-QPPS version
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % s.fl.extraInput=0;
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % s.fl.newNet = 1;
% % % % % % % % % % % % % % s.fl.newTable = 1;
% % % % % % % % % % % % % % [s w storedEps net Qtable]=RunRLfun(s)
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % s.rl.maxRetr=100;
% % % % % % % % % % % % % % [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % for ii=1:100
% % % % % % % % % % % % % % % %     [net Qtable] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % % % % % % % % % % % % % % %     ii
% % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % save('F:\Projects\DPPS\DefenseAgent\Results\Net_RTtask1_NormNoise_HarderTask_ExtraInput_0BaseRew_SimpleQPPS.mat','net','s','w','storedEps');
% % % % % % % % % % % % % % save('F:\Projects\DPPS\DefenseAgent\Results\Net_RTtask1_NormNoise_HarderTask_ExtraInput_0BaseRew_NO_extra_Q.mat','net','s','w','storedEps');



% % % % % % % % % % % % %% Data with MULTIPLE LIMBS,FOR NEW VERSION OF RT TASK
% % % % % % % % % % % % clear s rS ntRS rtRS ytRS olRS
% % % % % % % % % % % % 
% % % % % % % % % % % % disp('Here 6')
% % % % % % % % % % % % 
% % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % % Base settings
% % % % % % % % % % % % s.gol.alSpR = 1;
% % % % % % % % % % % % s.gol.alSpC = [0 0];
% % % % % % % % % % % % s.gol.randSpr = [2 2];
% % % % % % % % % % % % s.lp.b1Siz = 1e4;
% % % % % % % % % % % % s.lp.minExp = 4e7; % minumum number of actions before training sarts (big here so no training)
% % % % % % % % % % % % s.lp.dispA = 1e6; % only show action count very infreequently
% % % % % % % % % % % % % % % s.wrld.size = [14 15];
% % % % % % % % % % % % s.wrld.size = [10 15]; % $$$ DOING THIS ONE for simplicity
% % % % % % % % % % % % s.wrld.resetType = 'BottomTop_InfLR';
% % % % % % % % % % % % s.rp.maxActions = 4e6 ;
% % % % % % % % % % % % s.act.bdyGoalRew = 2;
% % % % % % % % % % % % s.lmb.startCol=8;
% % % % % % % % % % % % s.fl.trainTable=1; s.fl.trainNet=1;
% % % % % % % % % % % % s.rl.maxRetr=51;
% % % % % % % % % % % % s.fl.bdyMov = 1;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.perfWhileLearn = 1;
% % % % % % % % % % % % s.prf.nRep = 10;
% % % % % % % % % % % % s.prf.skipBatches=50;
% % % % % % % % % % % % 
% % % % % % % % % % % % netSizes={[13 11 9 7 5],[9 9 9 9 9],[5 7 9 11 13]};
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.rtTask=1;
% % % % % % % % % % % % s.actrtRew = 0;
% % % % % % % % % % % % 
% % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % % Run the model
% % % % % % % % % % % % [s w storedEps net Qtable] = RunRLfun(s);
% % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % % Run each type of model
% % % % % % % % % % % % for iM = 1:length(netSizes)
% % % % % % % % % % % %     
% % % % % % % % % % % %     s.fl.newNet = 1;
% % % % % % % % % % % %     s.fl.newTable = 1;
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Change specific settings of the model
% % % % % % % % % % % %     s.lp.netS = netSizes{iM};
% % % % % % % % % % % % 
% % % % % % % % % % % %     % Relearn Q with different netsizes
% % % % % % % % % % % %     [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % % % % % % % % % % %     % Store network results Structure
% % % % % % % % % % % %     rS(iM).s=s; rS(iM).Qtable=Qtable; rS(iM).w=w; rS(iM).net=net;
% % % % % % % % % % % %     rS(iM).perf = perf;
% % % % % % % % % % % %     
% % % % % % % % % % % %     bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\';
% % % % % % % % % % % %     svNm = [bFld 'Fig_Dist_Pos_Dependence_2LIMBS_FORrtTASK_v2.mat'];
% % % % % % % % % % % %     save(svNm,'rS','-v7.3');
% % % % % % % % % % % %     
% % % % % % % % % % % % end


% % % % % % % % % % % % %% Run a RT scenario with and without loaded QPPS, BUT WITH MULTIPLE LIMBS.
% % % % % % % % % % % % %  AND the re-learning of previously learned actions
% % % % % % % % % % % % clear s rS ntRS rtRS ytRS olRS
% % % % % % % % % % % % 
% % % % % % % % % % % % disp('Here 7')
% % % % % % % % % % % % 
% % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % % % First load normal Q-values of multi-limb model
% % % % % % % % % % % % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_2LIMBS_FORrtTASK_v2.mat');
% % % % % % % % % % % % 
% % % % % % % % % % % % iM=1;
% % % % % % % % % % % % 
% % % % % % % % % % % % s = rS(iM).s;
% % % % % % % % % % % % w = rS(iM).w;
% % % % % % % % % % % % net_PPS = rS(iM).net;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.plt.rowLims=[1.5 s.wrld.size(1)-1.5];
% % % % % % % % % % % % s.plt.colLims=[1.5 s.wrld.size(2)-2.5];
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % % Then run the reaction time scenario
% % % % % % % % % % % % % $$$ I'm redoing it witht the proper world size
% % % % % % % % % % % % 
% % % % % % % % % % % % % s.lp.netS = [9 9];
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.newTable = 1;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.rp.maxActions = 2e5 ;
% % % % % % % % % % % % s.lp.minExp = 2e6 ;
% % % % % % % % % % % % 
% % % % % % % % % % % % % s.rp.maxActions = 2e6 ;
% % % % % % % % % % % % % s.lp.minExp = 2e7 ;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.rtTask=1;
% % % % % % % % % % % % s.fl.hist=0;
% % % % % % % % % % % % s.actrtRew = 1;
% % % % % % % % % % % % 
% % % % % % % % % % % % % s.act.Name={'LEFT','STAY','RIGHT','BUTTON'};
% % % % % % % % % % % % % % s.act.Name={'lmbLEFT_bdyLEFT', 'lmbLEFT_bdySTAY', 'lmbLEFT_bdyRIGHT', ...
% % % % % % % % % % % % % %             'lmbSTAY_bdyLEFT', 'lmbSTAY_bdySTAY', 'lmbSTAY_bdyRIGHT', ...
% % % % % % % % % % % % % %             'lmbRIGHT_bdyLEFT','lmbRIGHT_bdySTAY','lmbRIGHT_bdyRIGHT'};
% % % % % % % % % % % % 
% % % % % % % % % % % % % s.act.Name={'lmbLEFT_bdyLEFT', 'lmbLEFT_bdySTAY', 'lmbLEFT_bdyRIGHT', ...
% % % % % % % % % % % % %             'lmbBUTTON_bdyLEFT', 'lmbBUTTON_bdySTAY', 'lmbBUTTON_bdyRIGHT', ...
% % % % % % % % % % % % %             'lmbRIGHT_bdyLEFT','lmbRIGHT_bdySTAY','lmbRIGHT_bdyRIGHT'};
% % % % % % % % % % % % 
% % % % % % % % % % % % % s.act.Name={'lmbLEFTBUTTON_bdyLEFT', 'lmbLEFTBUTTON_bdySTAY', 'lmbLEFTBUTTON_bdyRIGHT', ...
% % % % % % % % % % % % %             'lmbSTAY_bdyLEFT', 'lmbSTAY_bdySTAY', 'lmbSTAY_bdyRIGHT', ...
% % % % % % % % % % % % %             'lmbRIGHT_bdyLEFT','lmbRIGHT_bdySTAY','lmbRIGHT_bdyRIGHT'};
% % % % % % % % % % % % 
% % % % % % % % % % % % s.act.Name={'lmbLEFT_bdyLEFT', 'lmbLEFT_bdySTAY', 'lmbLEFT_bdyRIGHT', ...
% % % % % % % % % % % %             'lmbBUTTON_bdyLEFT', 'lmbBUTTON_bdySTAY', 'lmbBUTTON_bdyRIGHT', ...
% % % % % % % % % % % %             'lmbRIGHT_bdyLEFT','lmbRIGHT_bdySTAY','lmbRIGHT_bdyRIGHT'};
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % % % s.fl.dspm=1;
% % % % % % % % % % % % s.fl.dspm=0;
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % s.act.rtRew=1;
% % % % % % % % % % % % s.rtt.threshold=3;
% % % % % % % % % % % % 
% % % % % % % % % % % % for defaultRew = [0, -0.01]
% % % % % % % % % % % %     
% % % % % % % % % % % % %     clear rS
% % % % % % % % % % % %     
% % % % % % % % % % % %     s.act.BumpRew = defaultRew;
% % % % % % % % % % % %     s.act.GoalRew = defaultRew;
% % % % % % % % % % % %     s.act.ThreatRew = defaultRew;
% % % % % % % % % % % %     s.act.DefaultMoveRew = defaultRew;
% % % % % % % % % % % %     s.act.bdyBumpRew = defaultRew;
% % % % % % % % % % % %     s.act.bdyGoalRew = defaultRew;
% % % % % % % % % % % %     s.act.bdyThreatRew = defaultRew;
% % % % % % % % % % % %     s.act.bdyMoveRew = defaultRew;
% % % % % % % % % % % %     
% % % % % % % % % % % %     s=DefaultSettings(s);
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % % % %     s.rl.maxRetr=100;
% % % % % % % % % % % %     s.rl.maxRetr=10; %
% % % % % % % % % % % % 
% % % % % % % % % % % %     
% % % % % % % % % % % %     fromScratch = [0 1];
% % % % % % % % % % % %     
% % % % % % % % % % % %     
% % % % % % % % % % % %     s.fl.perfWhileLearn = 1;
% % % % % % % % % % % %     s.prf.nRep = 10;
% % % % % % % % % % % %     s.prf.skipBatches = 9;
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Run each type of model
% % % % % % % % % % % %     for iM = 1:length(fromScratch)
% % % % % % % % % % % %         
% % % % % % % % % % % %         
% % % % % % % % % % % %         s.fl.newTable = 1;
% % % % % % % % % % % %         
% % % % % % % % % % % %         % Change specific settings of the model
% % % % % % % % % % % %         s.fl.newNet = fromScratch(iM);
% % % % % % % % % % % %         
% % % % % % % % % % % %         % Run model with new settings (extra Q-inputs)
% % % % % % % % % % % %         
% % % % % % % % % % % % %         RunRLfun(s,net,Qtable)
% % % % % % % % % % % % %         [s w storedEps net Qtable]=RunRLfun(s,[],[],[],QPPS);
% % % % % % % % % % % %         [s w storedEps net Qtable]=RunRLfun(s,net_PPS);
% % % % % % % % % % % %         
% % % % % % % % % % % %         % Relearn Q with different netsizes
% % % % % % % % % % % % %         [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable,QPPS);
% % % % % % % % % % % %         [net Qtable perf] = RelearnNNFun(s,w,storedEps,net_PPS,Qtable);
% % % % % % % % % % % %         % Store network results Structure
% % % % % % % % % % % %         rtRS(iM).s=s; rtRS(iM).Qtable=Qtable; rtRS(iM).w=w; rtRS(iM).net=net;
% % % % % % % % % % % %         rtRS(iM).perf = perf;
% % % % % % % % % % % %         
% % % % % % % % % % % %         bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\';
% % % % % % % % % % % %         svNm = [bFld 'Fig_Dist_Pos_Dependence_RTtask_BaseRew_' num2str(defaultRew) '_noHist_EasierTask_RELEARN_9ACTIONS_lmbButtonReact_v2.mat'];
% % % % % % % % % % % %         save(svNm,'rtRS','-v7.3');
% % % % % % % % % % % %         
% % % % % % % % % % % %     end
% % % % % % % % % % % % end


% % % % % % % % % % % % %% Find performance of different network sizes for velocity and direction dependence
% % % % % % % % % % % % clear s rS ntRS rtRS ytRS olRS
% % % % % % % % % % % % 
% % % % % % % % % % % % disp('Here 8')
% % % % % % % % % % % % 
% % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % % Base settings
% % % % % % % % % % % % s.gol.alSpR = [1 2 3];
% % % % % % % % % % % % s.gol.alSpC = [-2 -1 0 1 2];
% % % % % % % % % % % % s.gol.randSpr = [2 2];
% % % % % % % % % % % % s.lp.b1Siz = 1e4;
% % % % % % % % % % % % s.wrld.size = [14 15];
% % % % % % % % % % % % s.wrld.resetType = 'BottomTop_InfLR';
% % % % % % % % % % % % s.lmb.startCol=8;
% % % % % % % % % % % % s.fl.trainTable=1; s.fl.trainNet=1;
% % % % % % % % % % % % s.rl.maxRetr=51;
% % % % % % % % % % % % % s.rl.maxRetr=1;
% % % % % % % % % % % % s.fl.hist = 1;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.perfWhileLearn = 1;
% % % % % % % % % % % % s.prf.nRep = 10;
% % % % % % % % % % % % s.prf.skipBatches=5;
% % % % % % % % % % % % 
% % % % % % % % % % % % netSizes={12.*ones(1,4),12.*ones(1,5),12.*ones(1,6), ...
% % % % % % % % % % % %     12.*ones(1,7),12.*ones(1,8),12.*ones(1,9), ...
% % % % % % % % % % % %     13.*ones(1,9),14.*ones(1,9)};
% % % % % % % % % % % % 
% % % % % % % % % % % % % Run the model and learn Q (table)
% % % % % % % % % % % % [s w storedEps net Qtable] = RunRLfun(s);
% % % % % % % % % % % % 
% % % % % % % % % % % % % Run each type of model
% % % % % % % % % % % % for iM = 1:length(netSizes)
% % % % % % % % % % % %     
% % % % % % % % % % % %     s.fl.newNet = 1;
% % % % % % % % % % % %     s.fl.newTable = 1;
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Change specific settings of the model
% % % % % % % % % % % %     s.lp.netS = netSizes{iM};
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Relearn Q with different netsizes
% % % % % % % % % % % %     [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % % % % % % % % % % %     % Store network results Structure
% % % % % % % % % % % %     ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
% % % % % % % % % % % %     ntRS(iM).perf = perf;
% % % % % % % % % % % %     
% % % % % % % % % % % %     bFld = 'F:\Projects\DPPS\DefenseAgent\Results\Performance\VelocityDirection\NetSizes\';
% % % % % % % % % % % %     svNm = [bFld 'Mults_of_12_51Batch_V3'];
% % % % % % % % % % % %     save(svNm,'ntRS','-v7.3');
% % % % % % % % % % % %     
% % % % % % % % % % % % end



%% Data for plot $3$ about velocity and direction dependence, and corresponding stats
clear s rS ntRS rtRS ytRS olRS

disp('Here 9')

addpath('F:\Projects\DPPS\DefenseAgent\Scripts\TESTCHECK - from supercomputer')

% -------------------------------------------------------------------------
% Base settings
s.gol.alSpR = [1 2 3];
s.gol.alSpC = [-2 -1 0 1 2];
s.gol.randSpr = [2 2];
s.lp.b1Siz = 1e4;
s.wrld.size = [14 15];
s.wrld.resetType = 'BottomTop_InfLR';
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
% s.rl.maxRetr=51;
s.rl.maxRetr=101;
% s.rl.maxRetr=1;
s.fl.hist = 1;

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=5;

netSizes={[8 9 10 11 12 13 14 15 16], 12.*ones(1,9), [16 15 14 13 12 11 10 9 8]};

% Run the model and learn Q (table)
[s w storedEps net Qtable] = RunRLfun(s);

% Run each type of model
for iM = 1:length(netSizes)
    
    s.fl.newNet = 1;
    s.fl.newTable = 1;
    
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};
    
    % Relearn Q with different netsizes
%     [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
    [net Qtable perf] = RelearnNNFunSC(s,w,storedEps,net,Qtable);
    % Store network results Structure
    rS(iM).s=s; rS(iM).Qtable=Qtable; rS(iM).w=w; rS(iM).net=net;
    rS(iM).perf = perf;
    
    bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\';
    svNm = [bFld 'Fig_Vel_Dir_Dependence_' num2str(s.rl.maxRetr) 'Batches_OLDRUNRELEARN'];
    save(svNm,'rS','-v7.3');
    
end


%%

% % % % % % % % % % % % %% Find performance of different network sizes for tool use
% % % % % % % % % % % % clear s rS ntRS rtRS ytRS olRS
% % % % % % % % % % % % 
% % % % % % % % % % % % disp('Here 10')
% % % % % % % % % % % % 
% % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % % Base settings
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.ToolChange=1;
% % % % % % % % % % % % s.lmb.ToolRows=[3];
% % % % % % % % % % % % s.lmb.ToolProb=0.5;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.hist = 0;
% % % % % % % % % % % % s.gol.alSpR = [1];
% % % % % % % % % % % % s.gol.alSpC =[0 0];
% % % % % % % % % % % % s.gol.randSpr = 0;
% % % % % % % % % % % % s.rp.maxActions = 2e6 ;
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % s.gol.randSpr = [2 2];
% % % % % % % % % % % % s.lp.b1Siz = 1e4;
% % % % % % % % % % % % s.wrld.size = [14 15];
% % % % % % % % % % % % s.wrld.resetType = 'BottomTop_InfLR';
% % % % % % % % % % % % s.lmb.startCol=8;
% % % % % % % % % % % % s.fl.trainTable=1; s.fl.trainNet=1;
% % % % % % % % % % % % s.rl.maxRetr=51;
% % % % % % % % % % % % % % % s.rl.maxRetr=1;
% % % % % % % % % % % % s.fl.hist = 1;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.perfWhileLearn = 1;
% % % % % % % % % % % % s.prf.nRep = 10;
% % % % % % % % % % % % s.prf.skipBatches=5;
% % % % % % % % % % % % 
% % % % % % % % % % % % netSizes={12.*ones(1,3),12.*ones(1,4),12.*ones(1,5), ...
% % % % % % % % % % % %     12.*ones(1,6),12.*ones(1,7),12.*ones(1,8)};
% % % % % % % % % % % % 
% % % % % % % % % % % % % Run the model and learn Q
% % % % % % % % % % % % [s w storedEps net Qtable] = RunRLfun(s);
% % % % % % % % % % % % 
% % % % % % % % % % % % % Run each type of model
% % % % % % % % % % % % for iM = 1:length(netSizes)
% % % % % % % % % % % %     
% % % % % % % % % % % %     s.fl.newNet = 1;
% % % % % % % % % % % %     s.fl.newTable = 1;
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Change specific settings of the model
% % % % % % % % % % % %     s.lp.netS = netSizes{iM};
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Relearn Q with different netsizes
% % % % % % % % % % % %     [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % % % % % % % % % % %     % Store network results Structure
% % % % % % % % % % % %     ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
% % % % % % % % % % % %     ntRS(iM).perf = perf;
% % % % % % % % % % % %     
% % % % % % % % % % % %     bFld = 'F:\Projects\DPPS\DefenseAgent\Results\Performance\ToolUse\NetSizes\';
% % % % % % % % % % % %     svNm = [bFld 'Mults_of_12_51Batch_v3'];
% % % % % % % % % % % %     save(svNm,'ntRS','-v7.3');
% % % % % % % % % % % %     
% % % % % % % % % % % % end


% % % % % % % % % % % % %% Data for plot $4$ about tool use
% % % % % % % % % % % % clear s rS ntRS rtRS ytRS olRS
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.ToolChange=1;
% % % % % % % % % % % % s.lmb.ToolRows=[5];
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.hist = 1;
% % % % % % % % % % % % s.gol.alSpR = [1];
% % % % % % % % % % % % s.gol.alSpC =[0 0];
% % % % % % % % % % % % s.gol.randSpr = [2 2];
% % % % % % % % % % % % % % % s.rp.maxActions = 4e6 ;
% % % % % % % % % % % % s.rp.maxActions = 2e6 ;
% % % % % % % % % % % % s.lp.minExp = 2e7;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.lp.b1Siz = 1e4;
% % % % % % % % % % % % s.wrld.size = [14 15];
% % % % % % % % % % % % s.wrld.resetType = 'BottomTop_InfLR';
% % % % % % % % % % % % s.lmb.startCol=8;
% % % % % % % % % % % % s.fl.trainTable=1; s.fl.trainNet=1;
% % % % % % % % % % % % s.rl.maxRetr=51;
% % % % % % % % % % % % % % % s.rl.maxRetr=2;
% % % % % % % % % % % % s.fl.hist = 0;
% % % % % % % % % % % % 
% % % % % % % % % % % % % % % s.fl.perfWhileLearn = 1;
% % % % % % % % % % % % s.fl.perfWhileLearn = 1;
% % % % % % % % % % % % s.prf.nRep = 10;
% % % % % % % % % % % % s.prf.skipBatches=25;
% % % % % % % % % % % % 
% % % % % % % % % % % % netSizes={[8 10 12 14 16], 12.*ones(1,5), [16 14 12 10 8]};
% % % % % % % % % % % % 
% % % % % % % % % % % % s.lmb.ToolProb=0.5;
% % % % % % % % % % % % % Run the model and learn Q
% % % % % % % % % % % % [s w storedEps ] = RunRLfun(s); % Run up to here to check that ther is at least 1 ep with a tool
% % % % % % % % % % % % % % % [s w storedEps net Qtable] = RunRLfun(s); % Run up to here to check that ther is at least 1 ep with a tool
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % % Run each type of model
% % % % % % % % % % % % for iM = 1:length(netSizes)
% % % % % % % % % % % %     
% % % % % % % % % % % % 
% % % % % % % % % % % %     s.fl.newNet = 1;
% % % % % % % % % % % %     s.fl.newTable = 1;
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Change specific settings of the model
% % % % % % % % % % % %     s.lp.netS = netSizes{iM};
% % % % % % % % % % % % 
% % % % % % % % % % % %     
% % % % % % % % % % % %     % ---------------------------------------------------------------------
% % % % % % % % % % % %     % Create stored epochs which have only non-tool epochs
% % % % % % % % % % % %     tmpS=cell2mat(storedEps.S);
% % % % % % % % % % % %     % no tool stored eps
% % % % % % % % % % % %     noToolInd = find(tmpS(:,4) == 0);
% % % % % % % % % % % %     sEnoTool.S = cell(length(noToolInd)+1,1);
% % % % % % % % % % % %     sEnoTool.prvS = sEnoTool.S;
% % % % % % % % % % % %     for iInd = 1:length(noToolInd)
% % % % % % % % % % % %         sEnoTool.S{iInd} = storedEps.S{noToolInd(iInd)};
% % % % % % % % % % % %         sEnoTool.prvS{iInd} = storedEps.prvS{noToolInd(iInd)};
% % % % % % % % % % % %     end
% % % % % % % % % % % %     sEnoTool.R = storedEps.R(noToolInd);
% % % % % % % % % % % %     sEnoTool.A = storedEps.A(noToolInd);
% % % % % % % % % % % %     % Add one tool epoch, so that the network doesn't ignore ToolUse inputs
% % % % % % % % % % % %     lastToolInd = find(tmpS(:,4) == 1 & abs(storedEps.R)<0.5,1);
% % % % % % % % % % % %     sEnoTool.S{end} = storedEps.S{lastToolInd};
% % % % % % % % % % % %     sEnoTool.prvS{end} = storedEps.prvS{lastToolInd};
% % % % % % % % % % % %     sEnoTool.R = [ sEnoTool.R; storedEps.R(lastToolInd) ];
% % % % % % % % % % % %     sEnoTool.A = [ sEnoTool.A; storedEps.A(lastToolInd) ];
% % % % % % % % % % % % 
% % % % % % % % % % % %     % Learn the no-tool Q values
% % % % % % % % % % % % % % % % %     s.rl.maxRetr=61; % $$$ HERE HERE HERE
% % % % % % % % % % % %     s.fl.newNet=1;
% % % % % % % % % % % %     s.fl.trainTable=0; s.fl.trainNet=1;
% % % % % % % % % % % %     [net Qtable perf] = RelearnNNFun(s,w,sEnoTool); 
% % % % % % % % % % % %     % Store Only Limb network results Structure
% % % % % % % % % % % %     olRS(iM).s=s; olRS(iM).Qtable=Qtable; olRS(iM).w=w; olRS(iM).net=net;
% % % % % % % % % % % %     olRS(iM).perf = perf;
% % % % % % % % % % % % 
% % % % % % % % % % % %     % Then adapt the same network, but with tool-use Q-values
% % % % % % % % % % % % % % % % %     s.rl.maxRetr=201; % $$$ HERE HERE HERE
% % % % % % % % % % % %     s.fl.newNet=0;
% % % % % % % % % % % %     [net Qtable] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % % % % % % % % % % %     % Store Yes Tool network results Structure
% % % % % % % % % % % %     ytRS(iM).s=s; ytRS(iM).Qtable=Qtable; ytRS(iM).w=w; ytRS(iM).net=net;
% % % % % % % % % % % %     ytRS(iM).perf = perf;
% % % % % % % % % % % %     
% % % % % % % % % % % %     bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\ToolUse\';
% % % % % % % % % % % %     svNm = [bFld 'Tool_Pre_post_51_51Batch_ToolPos5_moreRandSpr_NoHist_v2'];
% % % % % % % % % % % %     save(svNm,'olRS','ytRS','-v7.3');
% % % % % % % % % % % %     
% % % % % % % % % % % % end


% % % % % % % % % % % % %% Data for plot $4$ about tool use VERSION 2
% % % % % % % % % % % % clear s rS ntRS rtRS ytRS olRS
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.ToolChange=1;
% % % % % % % % % % % % s.lmb.ToolRows=[5];
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.hist = 1;
% % % % % % % % % % % % s.gol.alSpR = [1];
% % % % % % % % % % % % s.gol.alSpC =[0 0];
% % % % % % % % % % % % s.gol.randSpr = [2 2];
% % % % % % % % % % % % % % % s.rp.maxActions = 4e6 ;
% % % % % % % % % % % % s.rp.maxActions = 2e6 ;
% % % % % % % % % % % % s.lp.minExp = 2e7;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.lp.b1Siz = 1e4;
% % % % % % % % % % % % s.wrld.size = [14 15];
% % % % % % % % % % % % s.wrld.resetType = 'BottomTop_InfLR';
% % % % % % % % % % % % s.lmb.startCol=8;
% % % % % % % % % % % % s.fl.trainTable=1; s.fl.trainNet=1;
% % % % % % % % % % % % s.rl.maxRetr=51;
% % % % % % % % % % % % % % % s.rl.maxRetr=2;
% % % % % % % % % % % % s.fl.hist = 0;
% % % % % % % % % % % % 
% % % % % % % % % % % % % % % s.fl.perfWhileLearn = 1;
% % % % % % % % % % % % s.fl.perfWhileLearn = 1;
% % % % % % % % % % % % s.prf.nRep = 10;
% % % % % % % % % % % % s.prf.skipBatches=25;
% % % % % % % % % % % % 
% % % % % % % % % % % % netSizes={[8 10 12 14 16], 12.*ones(1,5), [16 14 12 10 8]};
% % % % % % % % % % % % 
% % % % % % % % % % % % s.lmb.ToolProb=0.5;
% % % % % % % % % % % % % Run the model and learn Q
% % % % % % % % % % % % [s w storedEps ] = RunRLfun(s); % Run up to here to check that ther is at least 1 ep with a tool
% % % % % % % % % % % % % % % [s w storedEps net Qtable] = RunRLfun(s); % Run up to here to check that ther is at least 1 ep with a tool
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % % Run each type of model
% % % % % % % % % % % % for iM = 1:length(netSizes)
% % % % % % % % % % % %     
% % % % % % % % % % % % 
% % % % % % % % % % % %     s.fl.newNet = 1;
% % % % % % % % % % % %     s.fl.newTable = 1;
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Change specific settings of the model
% % % % % % % % % % % %     s.lp.netS = netSizes{iM};
% % % % % % % % % % % % 
% % % % % % % % % % % %     % $$$$ HERE UPDATE IT SO IT DOES THE TOOL STUFF RIGHT
% % % % % % % % % % % %     % ---------------------------------------------------------------------
% % % % % % % % % % % %     % Create stored epochs which do not receive reward for tool touch
% % % % % % % % % % % %     tmpS=cell2mat(storedEps.S);
% % % % % % % % % % % %     tmpPrvS = cell2mat(storedEps.prvS);
% % % % % % % % % % % %     % no tool stored eps
% % % % % % % % % % % %     sEnoTool.S = storedEps.S;
% % % % % % % % % % % %     sEnoTool.prvS = storedEps.prvmS; 
% % % % % % % % % % % %     sEnoTool.R = storedEps.R;
% % % % % % % % % % % %     % replace rewarded with non-rewarded case
% % % % % % % % % % % %     tmpInd = (sEnoTool.R > 1 & tmpS(:,4) == 1 & tmpPrvS(:, 2) < 9 );
% % % % % % % % % % % %     sEnoTool.R(tmpInd) = storedEps.R(tmpInd) - 2;
% % % % % % % % % % % %     sEnoTool.A = storedEps.A;
% % % % % % % % % % % % 
% % % % % % % % % % % %     % Learn the no-tool Q values
% % % % % % % % % % % % % % % % %     s.rl.maxRetr=61; % $$$ HERE HERE HERE
% % % % % % % % % % % %     s.fl.newNet=1;
% % % % % % % % % % % %     s.fl.trainTable=0; s.fl.trainNet=1;
% % % % % % % % % % % %     [net Qtable perf] = RelearnNNFun(s,w,sEnoTool); 
% % % % % % % % % % % %     % Store Only Limb network results Structure
% % % % % % % % % % % %     olRS(iM).s=s; olRS(iM).Qtable=Qtable; olRS(iM).w=w; olRS(iM).net=net;
% % % % % % % % % % % %     olRS(iM).perf = perf;
% % % % % % % % % % % % 
% % % % % % % % % % % %     % Then adapt the same network, but with tool-use Q-values
% % % % % % % % % % % % % % % % %     s.rl.maxRetr=201; % $$$ HERE HERE HERE
% % % % % % % % % % % %     s.fl.newNet=0;
% % % % % % % % % % % %     [net Qtable] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % % % % % % % % % % %     % Store Yes Tool network results Structure
% % % % % % % % % % % %     ytRS(iM).s=s; ytRS(iM).Qtable=Qtable; ytRS(iM).w=w; ytRS(iM).net=net;
% % % % % % % % % % % %     ytRS(iM).perf = perf;
% % % % % % % % % % % %     
% % % % % % % % % % % %     bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\ToolUse\';
% % % % % % % % % % % %     svNm = [bFld 'Tool_Pre_post_51_51Batch_ToolPos5_moreRandSpr_NoHist_ToolAlwaysPresentButHalfNoReward'];
% % % % % % % % % % % %     save(svNm,'olRS','ytRS','-v7.3');
% % % % % % % % % % % %     
% % % % % % % % % % % % end




% % % % % % % % % % % % %% Find performance of different network sizes for valence effects - pos neg, plus magnitude
% % % % % % % % % % % % clear s rS ntRS rtRS ytRS olRS
% % % % % % % % % % % % 
% % % % % % % % % % % % disp('Here A')
% % % % % % % % % % % % 
% % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % % Base settings
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.thr = 1;
% % % % % % % % % % % % s.act.ThreatRew = -4;
% % % % % % % % % % % % s.act.DefaultMoveRew = -0.1;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.hist = 0;
% % % % % % % % % % % % s.gol.alSpR = [1];
% % % % % % % % % % % % s.gol.alSpC =[0 0];
% % % % % % % % % % % % s.gol.randSpr = [2 2];
% % % % % % % % % % % % s.thr.alSpR = [1];
% % % % % % % % % % % % s.thr.alSpC =[0 0];
% % % % % % % % % % % % s.thr.randSpr = [2 2];
% % % % % % % % % % % % s.rp.maxActions = 2e6 ;
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % s.lp.b1Siz = 1e4;
% % % % % % % % % % % % s.wrld.size = [14 15];
% % % % % % % % % % % % s.wrld.resetType = 'BottomTop_InfLR';
% % % % % % % % % % % % s.lmb.startCol=8;
% % % % % % % % % % % % % % % s.fl.trainTable=1; s.fl.trainNet=1;
% % % % % % % % % % % % s.rl.maxRetr=51; 
% % % % % % % % % % % % % s.rl.maxRetr=1;
% % % % % % % % % % % % s.fl.hist = 0;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.perfWhileLearn = 1;
% % % % % % % % % % % % s.prf.nRep = 10;
% % % % % % % % % % % % s.prf.skipBatches=25;
% % % % % % % % % % % % 
% % % % % % % % % % % % netSizes={9.*ones(1,3),9.*ones(1,4),9.*ones(1,5), ...
% % % % % % % % % % % %     9.*ones(1,6),9.*ones(1,7)};
% % % % % % % % % % % % 
% % % % % % % % % % % % % Run the model and learn Q
% % % % % % % % % % % % [s w storedEps net Qtable] = RunRLfun(s);
% % % % % % % % % % % % 
% % % % % % % % % % % % % Run each type of model
% % % % % % % % % % % % for iM = 1:length(netSizes)
% % % % % % % % % % % %     
% % % % % % % % % % % %     s.fl.newNet = 1;
% % % % % % % % % % % %     s.fl.newTable = 1;
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Change specific settings of the model
% % % % % % % % % % % %     s.lp.netS = netSizes{iM};
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Relearn Q with different netsizes
% % % % % % % % % % % %     [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % % % % % % % % % % %     % Store network results Structure
% % % % % % % % % % % %     ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
% % % % % % % % % % % %     ntRS(iM).perf = perf;
% % % % % % % % % % % %     
% % % % % % % % % % % %     bFld = 'F:\Projects\DPPS\DefenseAgent\Results\Performance\Valence\NetSizes\';
% % % % % % % % % % % %     svNm = [bFld 'Mults_of_9_51Batch_plus2_minus4_rewards_NoHist_DefRew-01_v2'];
% % % % % % % % % % % %     save(svNm,'ntRS','-v7.3');
% % % % % % % % % % % %     
% % % % % % % % % % % % end


% % % % % % % % % % % % %% Find performance of different network sizes for valence effects - magnitude
% % % % % % % % % % % % clear s rS ntRS rtRS ytRS olRS
% % % % % % % % % % % % 
% % % % % % % % % % % % disp('Here B')
% % % % % % % % % % % % 
% % % % % % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % % % % % % Base settings
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.thr = 1;
% % % % % % % % % % % % s.act.ThreatRew = 4;
% % % % % % % % % % % % s.act.DefaultMoveRew = -0.1;
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.hist = 0;
% % % % % % % % % % % % s.gol.alSpR = [1];
% % % % % % % % % % % % s.gol.alSpC =[0 0];
% % % % % % % % % % % % s.gol.randSpr = [2 2];
% % % % % % % % % % % % s.thr.alSpR = [1];
% % % % % % % % % % % % s.thr.alSpC =[0 0];
% % % % % % % % % % % % s.thr.randSpr = [2 2];
% % % % % % % % % % % % s.rp.maxActions = 2e6 ;
% % % % % % % % % % % % 
% % % % % % % % % % % % 
% % % % % % % % % % % % s.lp.b1Siz = 1e4;
% % % % % % % % % % % % s.wrld.size = [14 15];
% % % % % % % % % % % % s.wrld.resetType = 'BottomTop_InfLR';
% % % % % % % % % % % % s.lmb.startCol=8;
% % % % % % % % % % % % % % % s.fl.trainTable=1; s.fl.trainNet=1;
% % % % % % % % % % % % s.rl.maxRetr=51; 
% % % % % % % % % % % % % s.rl.maxRetr=1;
% % % % % % % % % % % % s.fl.hist = 0;
% % % % % % % % % % % % 
% % % % % % % % % % % % s.fl.perfWhileLearn = 1;
% % % % % % % % % % % % s.prf.nRep = 10;
% % % % % % % % % % % % s.prf.skipBatches=25;
% % % % % % % % % % % % 
% % % % % % % % % % % % % netSizes={12.*ones(1,3),12.*ones(1,4),12.*ones(1,5), ...
% % % % % % % % % % % % %     12.*ones(1,6),12.*ones(1,7),12.*ones(1,8)};
% % % % % % % % % % % % 
% % % % % % % % % % % % netSizes={9.*ones(1,3),9.*ones(1,4),9.*ones(1,5), ...
% % % % % % % % % % % %     9.*ones(1,6),9.*ones(1,7),9.*ones(1,8)};
% % % % % % % % % % % % 
% % % % % % % % % % % % % Run the model and learn Q
% % % % % % % % % % % % [s w storedEps net Qtable] = RunRLfun(s);
% % % % % % % % % % % % 
% % % % % % % % % % % % % Run each type of model
% % % % % % % % % % % % for iM = 1:length(netSizes)
% % % % % % % % % % % %     
% % % % % % % % % % % %     s.fl.newNet = 1;
% % % % % % % % % % % %     s.fl.newTable = 1;
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Change specific settings of the model
% % % % % % % % % % % %     s.lp.netS = netSizes{iM};
% % % % % % % % % % % %     
% % % % % % % % % % % %     % Relearn Q with different netsizes
% % % % % % % % % % % %     [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % % % % % % % % % % %     % Store network results Structure
% % % % % % % % % % % %     ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
% % % % % % % % % % % %     ntRS(iM).perf = perf;
% % % % % % % % % % % %     
% % % % % % % % % % % %     bFld = 'F:\Projects\DPPS\DefenseAgent\Results\Performance\Valence\NetSizes\';
% % % % % % % % % % % %     svNm = [bFld 'Mults_of_9_51Batch_plus2_plus4_rewards_NoHist_DefRew-01_v2'];
% % % % % % % % % % % %     save(svNm,'ntRS','-v7.3');
% % % % % % % % % % % %     
% % % % % % % % % % % % end



%% Data for plot $5$ about valence - all positive, with large moving cost
clear s rS ntRS rtRS ytRS olRS

% -------------------------------------------------------------------------
% Base settings

s.fl.thr = 1;
s.act.GoalRew = 1;
s.act.ThreatRew = 10;
s.act.DefaultMoveRew = -0.1;


s.fl.hist = 0;
s.gol.alSpR = [1];
s.gol.alSpC =[0 0];
s.gol.randSpr = [0 0];
s.thr.alSpR = [1];
s.thr.alSpC =[0 0];
s.thr.randSpr = [2 2];
s.rp.maxActions = 2e6 ;


s.lp.b1Siz = 1e4;
s.wrld.size = [14 15];
s.wrld.resetType = 'BottomTop_InfLR';
s.lmb.startCol=8;
% % % s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr=51; 
% s.rl.maxRetr=2;
s.fl.hist = 0;

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=25;

% netSizes={[6 8 10 12 14 16 18], [12 12 12 12 12 12 12], [18 16 14 12 10 8 6]};

netSizes={[6 10 14 18],    [12 12 12 12 ],    [18 14 10 6], ...
          [6 10 14 18 22], [12 12 12 12 12 ], [22 18 14 10 6], ...
          [15 15 15 15], [10 10 10 10 10 10 10]};


% Run the model and learn Q
[s w storedEps net Qtable] = RunRLfun(s);

% Run each type of model
for iM = 1:length(netSizes)
    
    s.fl.newNet = 1;
    s.fl.newTable = 1;
    
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};
    
%     net = rS(iM).net;
    
    % Relearn Q with different netsizes
    [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
    % Store network results Structure
    rS(iM).s=s; rS(iM).Qtable=Qtable; rS(iM).w=w; rS(iM).net=net;
    rS(iM).perf = perf;
    
    bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\';
    svNm = [bFld 'Valence_51Batch_Plus1_Plus10_minus01_movecost_NoHist_2_2RandSpr_v3_MoreNetSizes'];
    save(svNm,'rS','-v7.3');
    
    
end
%%




%% VERSION ORIGINAL)) -- Data for plot $5$ about valence - all POSNEG, with large moving cost
clear s rS ntRS rtRS ytRS olRS

% addpath('F:\Projects\DPPS\DefenseAgent\Scripts\TESTCHECK - from supercomputer')

% -------------------------------------------------------------------------
% Base settings
for iRep = 1:15 % $$$ HERE 2023: ADD IN WHICH UNIT TYPES ARE USED
    
    clear rS

s.fl.thr = 1;
s.act.ThreatRew = -2;
s.act.DefaultMoveRew = -0.1;


s.fl.hist = 0;
s.gol.alSpR = [1];
s.gol.alSpC =[0 0];
% % % s.gol.randSpr = [2 2];
s.gol.randSpr = [0 0];
s.thr.alSpR = [1];
s.thr.alSpC =[0 0];
s.thr.randSpr = [2 2];
% % % s.thr.randSpr = [0 0];
% % % s.rp.maxActions = 2e6 ;
s.rp.maxActions = 4e6 ;


s.lp.minExp = 4e7;
s.lp.b1Siz = 1e4;
s.wrld.size = [14 15];
s.wrld.resetType = 'BottomTop_InfLR';
s.lmb.startCol=8;
% % % s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr=51; 
% % % s.rl.maxRetr=2;
s.fl.hist = 0;

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=50;

% $$$ I'm trying the bigger net sizes now
% netSizes={[6 8 10 12 14 16 18], [12 12 12 12 12 12 12], [18 16 14 12 10 8 6]};
netSizes={[6 10 14 18], [12 12 12 12 ], [18 14 10 6]};


% Run the model and learn Q
s.lp.neurTypes = 'radbas'; %'softmax'; % 'logsig'; %'radbas'; % 'tribas'; %'tansig';  %'hardlim'; %'softmax'; %poslin'; % 'purelin' 
s.lp.TrnFn = 'trainlmL1';
% s.lp.TrnFn = 'trainlm';
s = DefaultSettings(s);
[s w storedEps net Qtable] = RunRLfun(s);

% Run each type of model
for iM = 1:length(netSizes)
    
    s.fl.newNet     = 1;
    s.fl.newTable   = 1;
    s.fl.trainNet   = 1;
    s.act.numA      = 3;
    s.act.Name      = {'LEFT','STAY','RIGHT'};
    s               = DefaultSettings(s);
%     
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};
    
%     net = rS(iM).net;

%     % $$$ HERE HERE 
%     addpath('F:\Projects\DPPS\DefenseAgent\Scripts\TESTCHECK - from supercomputer')
    
    % Relearn Q with different netsizes
%     [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
    [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
    % Store network results Structure
    rS(iM).s=s; rS(iM).Qtable=Qtable; rS(iM).w=w; rS(iM).net=net;
    rS(iM).perf = perf;
    
    bFld = 'Results\ForFigures\Valence\';

% % %     svNm = [bFld 'Valence_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_' s.lp.neurTypes '_B_V' num2str(iRep)];
    svNm = [bFld 'Valence_L1_Regularization_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_' s.lp.neurTypes '_B_V' num2str(iRep)];
    save(svNm,'rS','-v7.3');

    iM
    iRep
    
end


end


%% Data for plot $6$ about successor states
clear s rS ntRS rtRS ytRS olRS

% -------------------------------------------------------------------------
% Base settings
s.gol.alSpR = 1;
s.gol.alSpC = [0 0];
s.gol.randSpr = [0 3];
s.thr.alSpR = 1;
s.thr.alSpC = [0 0];
s.thr.randSpr = [2 2];
s.lp.b1Siz = 1e4;
s.wrld.size = [14 15];
s.wrld.resetType = 'BottomTop_InfLR';
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr=51;

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=25;

s.fl.thr = 1;

netSizes={[12 10 8 6],[9 9 9 9],[6 8 10 12], ...
          [12 10 8 6],[9 9 9 9],[6 8 10 12]};
      
gRews = [2  2  2  0  0  0];
tRews = [0  0  0 -2 -2 -2];


% Run each type of model
for iM=1:length(netSizes)
    
    % -------------------------------------------------------------------------
    % Run the model
    s.act.GoalRew = gRews(iM);
    s.act.ThreatRew = tRews(iM);
    [s w storedEps net Qtable] = RunRLfun(s);
    % -------------------------------------------------------------------------
    
    s.fl.newNet = 1;
    s.fl.newTable = 1;
    
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};
    
    
    [net Qtable] = RelearnNNFun(s,w,storedEps,net,Qtable);
    % Store results Structure
    rS(iM).s=s; rS(iM).w=w; rS(iM).net=net; 
    %rS(iM).Qtable=Qtable;
    if iM == length(netSizes)
        rS(iM).storedEps=storedEps;
    end
    
    bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\SuccessorState\';
    svNm = [bFld 'SuccessorState_51Batch_Plus2_OR_Minus2_movecost_NoHist_v3'];
    save(svNm,'rS','-v7.3');
end

%% Data for response to reviewer: wasp example
% INCLUDES 2 speeds, and limb being bound to body
% $$$ make sure the move costs happen at the right times

clear s rS ntRS rtRS ytRS olRS

% -------------------------------------------------------------------------
% Base settings


s.gol.alSpC = [0 0];
s.gol.randSpr = [0 0];
s.lp.b1Siz = 1e4;
s.lp.minExp = 4e7; % minumum number of actions before training sarts (big here so no training)
s.lp.dispA = 1e6; % only show action count very infreequently
s.wrld.size = [10 15];
% % % s.wrld.size = [10 21]; % $$$ DOING THIS ONE for simplicity
s.wrld.resetType = 'BottomTop_InfLR';
s.rp.maxActions = 4e6 ;
s.act.GoalRew = 2; 
s.act.bdyGoalRew = 0; % Only reward the agent for grabbing a goal with its hands
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr=100;
s.fl.bdyMov = 1;

% ---------------------------
% NEW BITS
s.fl.hist = 1;
s.act.bdyMovesLimb = 1;
s.act.bdyLimbProx  = 1;

s.act.DefaultMoveRew = -1e-2;
s.act.bdyMoveRew     = -5e-2;

s.gol.alSpR = [1 2];
% % % netSizes={[8 9 10 11 12 13 14 15 16], 12.*ones(1,9), [16 15 14 13 12 11 10 9 8]};
% % % netSizes={12.*ones(1,9), 12.*ones(1,12)};
% % % % % % netSizes={12.*ones(1,12)};

% netSizes={15.*ones(1,15), 20.*ones(1,20), 30.*ones(1,30)};

netSizes={15.*ones(1,20)};

s.rl.maxRetr=21;
s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches= 5;


s.fl.dspm = 0;

s.bdy.startCol = 10;
s.lmb.startCol = 10;

% ------------------------------------------
% $$$ WHEN ADDING IN THREAT
s.fl.thr = 1;
s.thr.alSpR = [1 2];
s.thr.alSpC = [0 0];
s.thr.randSpr = [2 2];

s.act.GoalRew = 2; 
s.act.bdyGoalRew = 0;


s.act.ThreatRew = -.5; % Put negative reward that isn't as big as body reward
s.act.bdyThreatRew = -4;
% $$$
% ------------------------------------------

% -------------------------------------------------------------------------
% Run the model
[s w storedEps net Qtable] = RunRLfun(s);
% -------------------------------------------------------------------------
% Run each type of model
for iM = 1:length(netSizes)
    
    s.fl.newNet = 1;
    s.fl.newTable = 1;
    
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};

    % Relearn Q with different netsizes
    [net Qtable perf] = RelearnNNFunSC(s,w,storedEps,net,Qtable);
    % Store network results Structure
    ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
    ntRS(iM).perf = perf;
    
    bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\';
% % %     svNm = [bFld 'TESTING_WASP_fixedlimbmove_increasedMovePunishment_MoreForBody_OnlyThreat_v2.mat'];
    svNm = [bFld 'TESTING_WASP_fixedlimbmove_increasedMovePunishment_MoreForBody_GOALANDTHREAT_BIGGERNETS_MoreLearning.mat'];
    save(svNm,'ntRS','storedEps','-v7.3');
    iM
end


%% $$$ Data for response to reviewer: Reward only after eating:

clear s rS ntRS rtRS ytRS olRS


% -------------------------------------------------------------------------
% Base settings
s.gol.alSpC = [0 0];
s.gol.randSpr = [0 0];
s.lp.b1Siz = 1e4;
s.lp.minExp = 4e7; % minumum number of actions before training sarts (big here so no training)
s.lp.dispA = 1e6; % only show action count very infreequently
s.wrld.size = [14 15];
s.wrld.resetType = 'BottomTop_InfLR';
s.rp.maxActions = 4e6 ;
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr= 100;


netSizes={[16 14 12 10 8 6]};


% ---------------------------
% Unique settings
s.fl.grabAndEat = 1;
s.act.eatRew = 1;
s.act.GoalRew = 0; 

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=10;

s.lp.gamma = 0.8;

% -------------------------------------------------------------------------
% Run the model
[s w storedEps net Qtable] = RunRLfun(s);
% -------------------------------------------------------------------------
% Learn the Q-values
for iM = 1:length(netSizes)
    
    s.fl.newNet = 1;
    s.fl.newTable = 1;
    
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};

    % Relearn Q with different netsizes
    [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
    % Store network results Structure
    ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
    ntRS(iM).perf = perf;
    
    bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\';
    svNm = [bFld 'EatingModel.mat'];
    save(svNm,'ntRS','storedEps','-v7.3');
    iM
end

%% Try a 'defend zone' model?
% $$$ New shield model uses body and limb, but limb contact doesn't offer
% ANY reward, and body contact punishes


% -------------------------------------------------------------------------
% Base settings
s.gol.alSpC = [0 0];
s.gol.randSpr = [0 0];
s.lp.b1Siz = 1e4;
s.lp.minExp = 4e7; % minumum number of actions before training sarts (big here so no training)
s.lp.dispA = 1e6; % only show action count very infreequently
s.wrld.size = [14 15];
s.wrld.resetType = 'BottomTop_InfLR';
s.rp.maxActions = 4e6 ;
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr= 100;



s.fl.bdyMov = 1;



netSizes={[16 14 12 10 8 6]};


% ---------------------------
% Unique settings
s.fl.grabAndEat =  0;
s.fl.defendZone =  0; %1;
s.act.shieldRew =  0%; -1;
s.act.GoalRew   =  0;

% $$$ NEW -----
s.act.bdyGoalRew = -2;
s.fl.bdyMov      =  1;

s.gol.randSpr   = [2 2];

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=10;



s.lp.gamma = 0.8;

% -------------------------------------------------------------------------
% Run the model
[s w storedEps net Qtable] = RunRLfun(s);


% % % ========================================================
% % % THIS WAS AN ATTEMPT BUT IT DOESNT ACTUALLY WORK. STILL IT'S TOO NEAT AND
% % % POTENTIALLY USEFUL TO DELETE
% % % Modify storedEps so that the agent gets punished at particular locations
% % storedEpsDefend = storedEps;
% % storedEpsDefend.R(storedEpsDefend.R > 0) = storedEpsDefend.R(storedEpsDefend.R > 0) - 1;
% % ss = cell2mat(storedEps.S);
% % ps = cell2mat(storedEps.prvS);
% % % % % % This is the normal reward rule
% % % % % toBeRewarded = ( (ps(:,1) - ps(:,3) ==  1 & ps(:,2) == 11 ) & storedEps.A(:) == 1 ) | ...
% % % % %                ( (ps(:,1) - ps(:,3) ==  0 & ps(:,2) == 11 ) & storedEps.A(:) == 2 ) | ...
% % % % %                ( (ps(:,1) - ps(:,3) == -1 & ps(:,2) == 11 ) & storedEps.A(:) == 3 );
% % % And modify stored eps to not have the 'holding something' input
% % storedEpsDefend.R(toBeRewarded) = storedEpsDefend.R(toBeRewarded) + 1;
% % storedEpsDefend.S     = cellfun(@(x) x(1:3), storedEps.S,'UniformOutput',0); 
% % storedEpsDefend.prvS  = cellfun(@(x) x(1:3), storedEps.prvS,'UniformOutput',0); 
% % % ========================================================


% -------------------------------------------------------------------------
% Learn the Q-values
for iM = 1:length(netSizes)
    
    s.fl.newNet = 1;
    s.fl.newTable = 1;
    
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};

    % Relearn Q with different netsizes
    [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
    % Store network results Structure
    ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
    ntRS(iM).perf = perf;
    
    bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\';
    svNm = [bFld 'ShieldModel_BodyMoves.mat'];
    save(svNm,'ntRS','storedEps','-v7.3');
    iM
end


%% Data for response to reviewer: on-policy learning

% -------------------------------------------------------------------------
% Base settings
s = DefaultSettings();
s.gol.alSpR = 1;
s.gol.alSpC = [0 0];
s.gol.randSpr = [2 2];
s.lp.b1Siz = 1e4;
s.wrld.size = [14 15];
s.wrld.resetType = 'BottomTop_InfLR';
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr=20;

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=5;

netSizes={[12 10 8 6],[9 9 9 9],[6 8 10 12]};


% --------------------------------------
% Unique settings
s.lp.epsilon = 0.8;
% Policy: unbalanced-epsilon-greedy
unbPi = @(optimal_action, n_actions, randNum) (randNum > s.lp.epsilon) * optimal_action + (randNum <= s.lp.epsilon) * 3;


allGoalRewards = [-2 2] ;



% Create threat stored epochs
storedEpsThr = storedEps;
storedEpsThr.R(storedEps.R > 0.1) = - storedEpsThr.R(storedEps.R > 0.1);

for iRew = 2:length(allGoalRewards)


s.act.GoalRew = allGoalRewards(iRew);


for iDepth = 1:3
    % Run each type of model

%     if iDepth == 1
%         % Run the model
%         [s w storedEps net Qtable] = RunRLfun(s);
%     end


    for iM=1:length(netSizes)


        if iDepth == 1
            s.fl.newNet = 1;
            s.fl.newTable = 1;
    

            % Make sure the correct stored epochs are used
            if allGoalRewards(iRew) < 0
                storedEpsSARSA{iDepth}      = storedEpsThr;
                storedEpsSARSAUnb{iDepth}   = storedEpsThr;
                storedEpsQ{iDepth}          = storedEpsThr;
            else
                storedEpsSARSA{iDepth}      = storedEps;
                storedEpsSARSAUnb{iDepth}   = storedEps;
                storedEpsQ{iDepth}          = storedEps;
            end

            netSARSA    = net;
            netSARSAUnb = net;
            netQ        = net;
        else
            % Set parameters for nth run, using the network to make decisions
            s.fl.newNet = 0;
            s.fl.newTable = 1;

            s.fl.newNet     = 0;
            s.lp.retr       = 1000000;
            s.rp.maxActions =  100000;
            % Change specific settings of the model
            s.lp.netS = netSizes{iM};

            s.lp.alg = 'SARSA';
            [s w storedEpsSARSA{iDepth} netSARSA Qtable] = RunRLfun(s,rSsarsa(iM,iDepth-1).net,[]);
                        
            sUnb = s;
            sUnb.act.pi = unbPi;
            [sUnb w storedEpsSARSAUnb{iDepth} netSARSAUnb Qtable] = RunRLfun(sUnb,rSsarsaUnb(iM,iDepth-1).net,[]);



            s.lp.alg = 'Q';
            [s w storedEpsQ{iDepth} netQ Qtable] = RunRLfun(s,rS(iM,iDepth-1).net,[]);
            
            netSARSA    = rSsarsa(iM,iDepth-1).net;
            netSARSAUnb = rSsarsaUnb(iM,iDepth-1).net;
            netQ        = rS(iM,iDepth-1).net;
        end

        % Change specific settings of the model
        s.lp.netS = netSizes{iM};

        s.lp.alg = 'SARSA';
        [net Qtable perf] = RelearnNNFun(s, w, storedEpsSARSA{iDepth}, ...
            netSARSA, Qtable);
        % Store results Structure: SARSA
        rSsarsa(iM,iDepth).s=s; rSsarsa(iM,iDepth).Qtable=Qtable;
        rSsarsa(iM,iDepth).w=w; rSsarsa(iM,iDepth).net=net;
        rSsarsa(iM,iDepth).perf=perf;

        
        sUnb = s;
        sUnb.act.pi = unbPi;
        [net Qtable perf] = RelearnNNFun(sUnb, w, storedEpsSARSAUnb{iDepth}, ...
            netSARSAUnb, Qtable);
        % Store results Structure: SARSA but with UNBALANCED epsilon policy
        rSsarsaUnb(iM,iDepth).s=s; rSsarsaUnb(iM,iDepth).Qtable=Qtable;
        rSsarsaUnb(iM,iDepth).w=w; rSsarsaUnb(iM,iDepth).net=net;
        rSsarsaUnb(iM,iDepth).perf=perf;


        s.lp.alg = 'Q';
        [net Qtable perf] = RelearnNNFun(s, w, storedEpsQ{iDepth}, ...
            netQ, Qtable);
        % Store results Structure: Q-learning
        rS(iM,iDepth).s=s; rS(iM,iDepth).Qtable=Qtable;
        rS(iM,iDepth).w=w; rS(iM,iDepth).net=net;
        rS(iM,iDepth).perf=perf;


        iM
    end

    if allGoalRewards(iRew) < 0
        rSsarsathr      = rSsarsa; 
        rSsarsaUnbthr   = rSsarsaUnb;
        rSthr           = rS;
    end
    
    svNm='F:\Projects\DPPS\DefenseAgent\Results\ForFigures\SARSA_WithUnbalanced_AndThreat.mat';
    save(svNm,'rSsarsa','rSsarsaUnb','rS','rSsarsathr','rSsarsaUnbthr','rSthr','-v7.3')

    iDepth
end

end





%% Data for response to reviewers: tool use
clear s rS ntRS rtRS ytRS olRS

% Load previous tool results
load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\ToolUse\Tool_Pre_post_51_51Batch_ToolPos5_moreRandSpr_NoHist_ToolAlwaysPresentButHalfNoReward.mat')


s.fl.ToolChange=1;
s.lmb.ToolRows=[5];

s.fl.hist = 1;
s.gol.alSpR = [1];
s.gol.alSpC =[0 0];
s.gol.randSpr = [2 2];
% % % s.rp.maxActions = 4e6 ;
s.rp.maxActions = 2e6 ;
s.lp.minExp = 2e7;

s.lp.b1Siz = 1e4;
s.wrld.size = [14 15];
s.wrld.resetType = 'BottomTop_InfLR';
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr=1;
% % % s.rl.maxRetr=2;
s.fl.hist = 0;

% % % s.fl.perfWhileLearn = 1;
s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches = 1;

netSizes={[8 10 12 14 16], 12.*ones(1,5), [16 14 12 10 8]};

s.lmb.ToolProb=0.5;
% Run the model and learn Q
[s w storedEps ] = RunRLfun(s); % Run up to here to check that ther is at least 1 ep with a tool
% % % [s w storedEps net Qtable] = RunRLfun(s); % Run up to here to check that ther is at least 1 ep with a tool

% Start from the 'only limb' model, and learn tool use, storing the network
% after every single batch
for iM = 1:length(netSizes)
    
    % Change specific settings of the model
    s.lp.netS = netSizes{iM};

    % Set the first network to be equal to the no-tool network
    incrToolRS(iM,1) = olRS(iM);

    for iBatch = 7:51
        % Then adapt the no-tool network, but with tool-use Q-values
        s.fl.newNet=0;
        [net Qtable perf] = RelearnNNFun(s,w,storedEps,incrToolRS(iM,iBatch - 1).net,Qtable);
        % Store Yes Tool network results Structure
        incrToolRS(iM,iBatch).s=s; incrToolRS(iM,iBatch).Qtable=Qtable; incrToolRS(iM,iBatch).w=w; incrToolRS(iM,iBatch).net=net;
        incrToolRS(iM,iBatch).perf = perf;

        bFld = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\';
        svNm = [bFld 'Tool_ReviewerResponse_Pre_post_51_51Batch_ToolPos5_moreRandSpr_NoHist_V2'];
        save(svNm,'incrToolRS','-v7.3');
    end
    
end



%% Make videos showing what's happening during training

load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_SmallNetRelearn_B_V15.mat')
% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_v5_randRC.mat');
% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_2LIMBS_v2.mat'); rS = ntRS;
% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\Fig_Vel_Dir_Dependence_101Batches_OLDRUNRELEARN.mat')


iM = 1;

net                 = rS(iM).net;
s                   = DefaultSettings(rS(iM).s);
s.fl.newNet         = 0;
s.fl.dspm           = 1;
s.rp.maxActions     = 100;
s.fl.vid            = 1;
s.plt.fancyWorld    = 1;
s.plt.vidFR         = 7.5;
% s.plt.vidFileName   = 'F:\Projects\DPPS\DefenseAgent\Documentation\Figures\Finalising_All\TEST.avi';
s.plt.vidFileName   = 'F:\Projects\DPPS\DefenseAgent\Documentation\Figures\Finalising_All\LearnVidLate1d4.avi';

s.plt.grid = 1;

[s2, w2, storedEps, ~, ~]    = RunRLfun(s,net);


% rS(iM).perf.rewPerAct(end,:)

sum(storedEps.R)./numel(storedEps.R)


%% Also show pretraining

iM = 1;

net                 = rS(iM).net;
s                   = DefaultSettings(rS(iM).s);
s.fl.newNet         = 1;
s.fl.dspm           = 1;
s.rp.maxActions     = 100;
s.fl.vid            = 1;
s.plt.fancyWorld    = 1;
s.plt.vidFR         = 7.5;
% s.plt.vidFileName   = 'F:\Projects\DPPS\DefenseAgent\Documentation\Figures\Finalising_All\TEST2_grid_images_redbluefigsAndBackgrounds.avi';
s.plt.vidFileName   = 'F:\Projects\DPPS\DefenseAgent\Documentation\Figures\Finalising_All\LearnVidNoTrain1d4.avi';

s.plt.grid = 1;

[s2, w2, storedEps, ~, ~]    = RunRLfun(s);


%% Make video showing the limb bound to body


iM = 1;

net                  = rS(iM).net;
s2                   = DefaultSettings(rS(iM).s);
s2.fl.newNet         = 0;
s2.fl.dspm           = 1;
s2.rp.maxActions     = 100;
s2.fl.vid            = 1;
s2.plt.fancyWorld    = 1;
s2.plt.vidFR         = 5; %7.5;
s2.plt.vidFileName   = 'F:\Projects\DPPS\DefenseAgent\Documentation\Figures\Finalising_All\TESTWASP_example_limbWithBody7.avi';

s2.plt.grid = 1;


s2.bdy.startCol = 11;
s2.lmb.startCol = 11;

[s3, w3, storedEps, ~, ~]    = RunRLfun(s2,net);


% rS(iM).perf.rewPerAct(end,:)

sum(storedEps.R)./numel(storedEps.R)

%%

% % % %% Find performance of different network sizes for FULL EVERYTHING MODEL
% % % clear s rS ntRS rtRS ytRS olRS
% % % 
% % % % -------------------------------------------------------------------------
% % % % Base settings
% % % 
% % % % BODY
% % % s.fl.bdyMov = 1;
% % % s.act.bdyGoalRew = 2;
% % % s.act.bdyThreatRew = -4;
% % % 
% % % % TOOL
% % % s.fl.ToolChange = 1;
% % % s.lmb.ToolRows = [5];
% % % s.lmb.ToolProb = 0.5;
% % % 
% % % % THREAT
% % % s.fl.thr = 1;
% % % s.act.ThreatRew = -2;
% % % s.act.DefaultMoveRew = -0.1;
% % % 
% % % % KINEMATICS
% % % s.fl.hist = 0;
% % % s.gol.alSpR = [1 2 3];
% % % s.gol.alSpC = [-2 -1 0 1 2];
% % % s.gol.randSpr = [2 2];
% % % s.thr.alSpR = [1 2 3];
% % % s.thr.alSpC = [-2 -1 0 1 2];
% % % s.thr.randSpr = [0 3];
% % % 
% % % 
% % % % LEARNING PARAMS
% % % s.rp.maxActions = 4e6 ;
% % % s.lp.minExp = 4e7; % minumum number of actions before training sarts (big here so no training)
% % % s.lp.sWs=4e7; % minumum number of actions before saving sarts (big here so no training)
% % % s.lp.dispA = 1e6; % only show action count very infreequently
% % % % s.lp.b1Siz = 1e4;
% % % s.lp.bSiz = 14e4;
% % % 
% % % % WORLD
% % % s.wrld.size = [14 15];
% % % s.wrld.resetType = 'BottomTop_InfLR';
% % % s.lmb.startCol=8;
% % % 
% % % 
% % % % RELEARNING 
% % % % % % s.fl.trainTable=1; s.fl.trainNet=1;
% % % % % % s.rl.maxRetr=51; 
% % % s.rl.maxRetr=1;
% % % s.fl.hist = 1;
% % % 
% % % % % % s.fl.perfWhileLearn = 1;
% % % s.fl.perfWhileLearn = 1;
% % % s.prf.nRep = 10;
% % % s.prf.skipBatches=25;
% % % 
% % % % netSizes={  12.*ones(1,4),12.*ones(1,5), 12.*ones(1,6),...
% % % %             12.*ones(1,7),12.*ones(1,8),12.*ones(1,8)};
% % % 
% % % % netSizes={  12.*ones(1,7),12.*ones(1,9), 12.*ones(1,11),...
% % % %             20.*ones(1,5),20.*ones(1,7),20.*ones(1,9)};
% % % 
% % % % netSizes={  12.*ones(1,14), ...
% % % %             20.*ones(1,10)};
% % %         
% % % netSizes={  30.*ones(1,12)};
% % % 
% % % % netSizes={  [] };
% % % 
% % % 
% % % % Run the model and learn Q
% % % [s w storedEps net Qtable] = RunRLfun(s);

% $$$ DONE UNTIL HERE


%%

s.lp.ShwW = 1;
s.lp.mxF = 2;

s.rl.maxRetr=1; 
s.prf.skipBatches=1;
totalBatch = 1;

s.fl.newNet = 1;

s.lp.b1Siz = 14e4;
s.lp.bSiz = 14e4;

for iRep = 11:20
    
    % Run each type of model
    for iM = 1:length(netSizes)
        
        if iRep >1
            s.fl.newNet = 0;
        end
        s.fl.newTable = 0;
        
        % Change specific settings of the model
        s.lp.netS = netSizes{iM};
        
        % Relearn Q with different netsizes
        %     [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
        [net Qtable perf] = RelearnNNFun(s,w,storedEps,net);
        % Store network results Structure
        ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
        ntRS(iM).perf = perf;
        
        totalBatch = totalBatch + s.rl.maxRetr;
        
        bFld = 'F:\Projects\DPPS\DefenseAgent\Results\Performance\FullModel\NetSizes\';
        svNm = [bFld 'Body_30_' num2str(totalBatch) 'Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE_v2'];
        save(svNm,'ntRS','-v7.3');
        
        
    end
    
end

%% Make data for new learning algorithm test


%% Find performance of different network sizes for MULTIPLE LIMBS and proximity position dependence
clear s rS ntRS rtRS ytRS olRS

% -------------------------------------------------------------------------
% Base settings
s.gol.alSpR = 1;
s.gol.alSpC = [0 0];
s.gol.randSpr = [2 2];
s.lp.b1Siz = 1e4;
s.lp.minExp = 4e7; % minumum number of actions before training sarts (big here so no training)
s.lp.dispA = 1e6; % only show action count very infreequently


s.wrld.size = [14 15];
% s.wrld.size = [10 15]; % $$$ DOING THIS ONE for simplicity
s.wrld.resetType = 'BottomTop_InfLR';
s.rp.maxActions = 4e6 ;
s.lmb.startCol=8;
s.fl.trainTable=1; s.fl.trainNet=1;
s.rl.maxRetr=101;

s.fl.thr = 1;


% NEW SETTINGS SPECIFIC TO DEEP LEARNING
s.lp.n1stTr = 1000;
s.lp.maxRetrainSteps = 10;

s.fl.perfWhileLearn = 1;
s.prf.nRep = 10;
s.prf.skipBatches=50;

s.lp.b1Siz = 400000;
s.lp.bSiz = 200000;

% netSizes={[9],[9 9],[9 9 9],[9 9 9 9],[9 9 9 9 9],[9 9 9 9 9], ...
%     9*ones(1,6),9*ones(1,7),9*ones(1,8),9*ones(1,9),9*ones(1,10),9*ones(1,11)};

% -------------------------------------------------------------------------
% Run the model
[s w storedEps net Qtable] = RunRLfun(s);
% -------------------------------------------------------------------------

s.fl.newNet = 1;
s.fl.deepNet = 1;

tic
[net Qtable perf] = RelearnNNFun(s,w,storedEps);
toc
% % % % Run each type of model
% % % for iM = 1:length(netSizes)
% % %     
% % %     s.fl.newNet = 1;
% % %     s.fl.newTable = 1;
% % %     
% % %     % Change specific settings of the model
% % %     s.lp.netS = netSizes{iM};
% % % 
% % %     % Relearn Q with different netsizes
% % %     [net Qtable perf] = RelearnNNFun(s,w,storedEps,net,Qtable);
% % %     % Store network results Structure
% % %     ntRS(iM).s=s; ntRS(iM).Qtable=Qtable; ntRS(iM).w=w; ntRS(iM).net=net;
% % %     ntRS(iM).perf = perf;
% % %     
% % %     bFld = 'F:\Projects\DPPS\DefenseAgent\Results\Performance\ProximityPosition\NetSizes\';
% % %     svNm = [bFld '2LIMBS_Mults_of_9_V3'];
% % %     save(svNm,'ntRS','-v7.3');
% % %     
% % % end

%% 

s.rl.maxRetr=1;


% NEW SETTINGS SPECIFIC TO DEEP LEARNING
s.lp.n1stTr = 2000;
s.lp.b1Siz = 4e6;
s.lp.bSiz = 4e6;

s.lp.deepMiniBatchSize = 5e5;

% s.fl.newNet = 0;
% [net Qtable perf] = RelearnNNFun(s,w,storedEps,net);
s.fl.newNet = 1;
[net Qtable perf] = RelearnNNFun(s,w,storedEps);

%% $$$ TEMP - plot
% Settings for plot
sFP=s;
sFP.plt.lmbCol = 3:s.wrld.size(2)-2;
sFP.plt.ON = 1;
sFP.plt.rowLims = [1.5 s.wrld.size(1)-0.5];
sFP.plt.sequentialLimbCols=0;
sFP.plt.pltType = 'Imagesc';

% 
% sFP.plt.meanLimbCols = 1;
% sFP.plt.lmbCol=3:12;

sFP.plt.meanLimbCols = 0;
sFP.plt.lmbCol=5;

% sFP.plt.bdyCol = 9;

sFP.plt.plAct=3;

sFP.plt.rowLag = 1;
sFP.plt.colLag = 0;

sFP.plt.stimRow=[3:size(w.world2D,1)-1];
sFP.plt.stimCol=[2:size(w.world2D,2)-1];
%     sFP.plt.pltType = 'Binned';

sFP = DefaultSettings(sFP);

% 	set(0, 'currentfigure', f);
subplot(1,2,1)
tmp = s.prf.skipBatches;
plot(perf.rewPerAct(1:tmp:end,1),'LineWidth',2);
hold on;

%     f2 = figure,
subplot(1,2,2)
[Q,allNeurAct] = CalcNetOutput(sFP,w,net);
DisplActValsFun(sFP,w,Q);


