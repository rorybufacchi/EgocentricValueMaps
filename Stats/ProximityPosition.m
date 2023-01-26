%% Assess the effect of network size on performance

% % % $$$ NOTE I need to change the generating script to not make new copies of
% % % this
% % load('F:\Projects\DPPS\DefenseAgent\Results\Performance\ProximityPosition\NetSizes\Mults_of_9_v2.mat');

% % load('F:\Projects\DPPS\DefenseAgent\Results\Performance\ProximityPosition\NetSizes\2LIMBS_Mults_of_9.mat')

% % % % % % load('F:\Projects\DPPS\DefenseAgent\Results\Performance\ProximityPosition\NetSizes\2LIMBS_Mults_of_9_V4.mat')

addpath('F:\Projects\DPPS\DefenseAgent\Scripts\rlsimplepps')
addpath('F:\Projects\DPPS\DefenseAgent\Scripts')
addpath(genpath('F:\Programs\Matlab\Utilities\'))
addpath('F:\Programs\Matlab\Utilities\plotting\colormaps\')

% % % % % % % f = figure,
% % % % % % for iM = 1:length(ntRS)
% % % % % %     
% % % % % %     s = ntRS(iM).s;
% % % % % %     w = ntRS(iM).w;
% % % % % %     net = ntRS(iM).net;
% % % % % %     
% % % % % %     % Settings for plot
% % % % % %     sFP=s;
% % % % % %     sFP.plt.lmbCol = 3:s.wrld.size(2)-2;
% % % % % %     sFP.plt.ON = 1;
% % % % % %     sFP.plt.rowLims = [1.5 s.wrld.size(1)-0.5];
% % % % % %     sFP.plt.sequentialLimbCols=0;
% % % % % %     
% % % % % %     sFP.plt.pltType = 'Imagesc';
% % % % % %     
% % % % % %     % %     sFP.plt.meanLimbCols = 1;
% % % % % %     % %     sFP.plt.lmbCol=3:12;
% % % % % %     
% % % % % %     sFP.plt.meanLimbCols = 0;
% % % % % %     sFP.plt.lmbCol=11;
% % % % % %     
% % % % % %     sFP.plt.bdyCol=3;
% % % % % %     
% % % % % %     sFP.plt.plAct=2;
% % % % % %     
% % % % % % %     sFP.plt.colLag = -1;
% % % % % %     
% % % % % %     sFP.plt.stimRow=[3:size(w.world2D,1)-1];
% % % % % %     sFP.plt.stimCol=[2:size(w.world2D,2)-1];
% % % % % % %     sFP.plt.pltType = 'Binned';
% % % % % %     
% % % % % %     sFP = DefaultSettings(sFP);
% % % % % %     
% % % % % %     
% % % % % %     % 	set(0, 'currentfigure', f);
% % % % % %     subplot(1 + ceil(length(ntRS)/2),2,1)
% % % % % %     tmp = s.prf.skipBatches;
% % % % % %     plot(ntRS(iM).perf.rewPerAct(1:tmp:end,1),'LineWidth',2);
% % % % % %     hold on;
% % % % % %     
% % % % % %     %     f2 = figure,
% % % % % %     subplot(1 + ceil(length(ntRS)/2),2,1+iM)
% % % % % %     [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
% % % % % %     DisplActValsFun(sFP,w,Q);
% % % % % %     
% % % % % % end
% % % % % % subplot(1 + ceil(length(ntRS)/2),2,1)
% % % % % % % set(0, 'currentfigure', f);
% % % % % % legend('1','2','3','4','5','6','7','8','9','10','11');
% % % % % % % hold off



%% Work with the given network sizes

% Base model with 1 limb
% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_v4.mat');
load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_v5_randRC.mat');
addpath(genpath('F:\Projects\DPPS\DefenseAgent\Scripts\rlsimplepps'));
tmprS = rS;
if ~isfield(tmprS,'perf')
    for iM=1:length(tmprS)
        tmprS(iM).perf = [];
    end
end
rS = tmprS;


% Base model with multiple limbs
load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_2LIMBS_v2.mat');
if ~isfield(ntRS,'perf')
    for iM=1:length(tmprS)
        ntRS(iM).perf = [];
    end
end
tmprS(end+1:end+length(ntRS)) = ntRS;
rS = tmprS;

% -------------------------------------------------------------------------
% extra model which includes all other stuff ($$$ MAYBE except multple limbs)
load('F:\Projects\DPPS\DefenseAgent\Results\Performance\FullModel\NetSizes\Body_50_130Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE_FewerSpeeds_SmallWorld_V3.mat');
% rS(length(rS)+1).s = s;
% rS(end).w = w;
% rS(end).net = net;
rS(length(rS)+1).s  = ntRS.s;
rS(end).w           = ntRS.w;
rS(end).net         = ntRS.net;
rS(end).perf        = ntRS.perf;

clear tmprS


% % % extra model which includes multiple limbs
% % load('F:\Projects\DPPS\DefenseAgent\Results\Net_Goal_3_2_Rew_BodyAndHand_BodyMoves_InfLR.mat')
% % rS(length(rS)+1).s = s;
% % rS(end).w = w;
% % rS(end).net = net;

%% Also load up the 'reaction time' task

% $$$ THIS ONE NEEDS TO BE REDONE

% % load('F:\Projects\DPPS\DefenseAgent\Results\Net_RTtask1_NormNoise_HarderTask_ExtraInput_0BaseRew_NO_extra_Q.mat')
% % rtRS(1).s = s;
% % rtRS(1).w = w;
% % rtRS(1).net = net;
% % rtRS(1).Qtable = [];
% % rtRS(1).perf = [];
% % load('F:\Projects\DPPS\DefenseAgent\Results\Net_RTtask1_NormNoise_HarderTask_ExtraInput_0BaseRew_SimpleQPPS.mat')
% % rtRS(2).s = s;
% % rtRS(2).w = w;
% % rtRS(2).net = net;
% % rtRS(2).Qtable = [];
% % rtRS(2).perf = [];

% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_RTtask_BaseRew_-0.01_noHist_EasierTask_RELEARN_9ACTIONS_lmbButtonReact.mat')
load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_RTtask_BaseRew_-0.01_noHist_EasierTask_Only1batch_smallNet_QPPS2_v2.mat')
% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_RTtask.mat');
%% -------------------------------------------------------------------------
% Loop through current Model, to find limb-specific activations

clear allNeurAct_PlusBody Q_PlusBody
for cM = 1:length(rS)
    
    s = rS(cM).s;
    w = rS(cM).w;
    net = rS(cM).net;
    Qtable = rS(cM).Qtable;
    
    s = DefaultSettings(s);
    s.plt.plotThrFl = 0;
    
    % settings for plot
    sFP = s;
    %     sFP.plt.lmbCol = 2:s.wrld.size(2)-1;
    sFP.plt.lmbCol              = 3:s.wrld.size(2)-2;
    sFP.plt.ON                  = 0;
    sFP.plt.sequentialLimbCols  = 0;
    sFP.plt.stimRow             = [3:size(w.world2D,1)-1];
    sFP.plt.stimCol             = [2:size(w.world2D,2)-1];
    sFP.plt.pltType             = 'Binned';
    
    if s.fl.bdyMov == 1
        clear Q_PlusBody allNeurAct_PlusBody
        for iBdyCol = 1:size(w.world2D,2)
            s.plt.bdyCol = iBdyCol;
            
            [Q_tmp, allNeurAct_tmp] = CalcNetOutput(s,w,net);
            
            %             if cM>1
            %                 [Q_PlusBody]  =  MatExpand(Q_PlusBody, Q_tmp);
            %                 [allNeurAct_PlusBody]  =  MatExpand(allNeurAct_PlusBody, allNeurAct_tmp);
            %             end
            
            Q_PlusBody(:,:,:,:,:,iBdyCol) = Q_tmp;
            allNeurAct_PlusBody(:,:,:,:,:,:,iBdyCol) = allNeurAct_tmp;
        end
        
        allNeurAct_MeanBody = nanmean(allNeurAct_PlusBody(:,:,:,:,:,:,2:size(w.world2D,2)-1),7);
        allNeurAct_MeanLimb = permute(nanmean( ...
            allNeurAct_PlusBody(:,:,:,2:size(w.world2D,2)-1,:,:,:),4),[1 2 3 7 5 6 4]);
        
        Q_MeanBody = nanmean(Q_PlusBody(:,:,:,:,:,2:size(w.world2D,2)-1),6);
        Q_MeanLimb = permute(nanmean(Q_PlusBody(:,2:size(w.world2D,2)-1,:,:,:,:),2),[1 6 3 4 5 2]);
        
        
        % % %         % $$$ Calculate the limb and body activations for when the limb and
        % % %         % body are maximally separated.
        % % %         % OK SO, turns out this isn'r actually necessary, and the problem
        % % %         % is htat it removes the cDAbs correlation, because th limb is
        % % %         % always in a specific position relative to the body. SO, I now
        % % %         % feel more comfortable just going with the straight up averaging
        % % %         for bdyC = 1:s.wrld.size(2)
        % % %
        % % %             d = floor((s.wrld.size(2)-2)/2);
        % % %             tmp = bdyC + [d, - d];
        % % %             tmp = tmp + [1:ceil(d/2)]' - floor (d/4); % Use a few limb columns around the maximal distance as well
        % % %             lmbC = unique(tmp(  tmp > 1 & tmp < s.wrld.size(2)-1 ));
        % % % %             lmbC=2:14;
        % % %
        % % %             % $$$ NOW I JUST need to average it over the correct limbs
        % % %
        % % %             % Pick the limb position maximally distant
        % % %             % The -1 is because the body columns don't run over
        % % %             % everything.. BUt then neither do the limb cols?..
        % % %             allNeurAct_MeanBody(:,:,:,bdyC,:,:) = nanmean(allNeurAct_PlusBody(:,:,:,bdyC,:,:,lmbC),7); % maybe last dim -1
        % % %             allNeurAct_MeanLimb(:,:,:,bdyC,:,:) = nanmean(allNeurAct_PlusBody(:,:,:,lmbC,:,:,bdyC),4);
        % % %
        % % %             Q_MeanBody(:,bdyC,:,:,:) = nanmean(Q_PlusBody(:,bdyC,:,:,:,lmbC),6);
        % % %             Q_MeanLimb(:,bdyC,:,:,:) = nanmean(Q_PlusBody(:,lmbC,:,:,:,bdyC),2);
        % % %         end
        
        
        
        allNeurAct = allNeurAct_MeanBody;
        Q = Q_MeanBody;
        
        % This returns a correlation value for each (neuron, layer, limbcol)0.05;
        [rS(cM).bdy.rDistN, rS(cM).bdy.pDistN, rS(cM).bdy.hProxN, aD, rD, cD, ~, ~, rS(cM).bdy.rUDistN, rS(cM).bdy.rLDistN] = ...
            CalcDistCorr(sFP,w,permute(allNeurAct_MeanLimb,[3 4 1 2 5 6]));
        % and this for each action q-value
        [rS(cM).bdy.rDistQ, rS(cM).bdy.pDistQ, rS(cM).bdy.hProxQ, aDQ, rDQ, cDQ, ~, ~, rS(cM).bdy.rUDistQ, rS(cM).bdy.rLDistQ] = ...
            CalcDistCorr(sFP,w,Q_MeanLimb);
        
        % This returns a correlation value for each (neuron, layer, limbcol)0.05;
        [rS(cM).lmb.rDistN, rS(cM).lmb.pDistN, rS(cM).lmb.hProxN, aD, rD, cD, ~, ~, rS(cM).lmb.rUDistN, rS(cM).lmb.rLDistN] = ...
            CalcDistCorr(sFP,w,permute(allNeurAct_MeanBody,[3 4 1 2 5 6]));
        % and this for each action q-value
        [rS(cM).lmb.rDistQ, rS(cM).lmb.pDistQ, rS(cM).lmb.hProxQ, aDQ, rDQ, cDQ, ~, ~, rS(cM).lmb.rUDistQ, rS(cM).lmb.rLDistQ] = ...
            CalcDistCorr(sFP,w,Q_MeanBody);
    else
        [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
    end
    
    % This returns a correlation value for each (neuron, layer, limbcol)0.05;
    [rS(cM).rDistN, rS(cM).pDistN, rS(cM).hProxN, aD, rD, cD, ~, ~, rS(cM).rUDistN, rS(cM).rLDistN] = ...
        CalcDistCorr(sFP,w,permute(allNeurAct,[3 4 1 2 5 6]));
    % and this for each action q-value
    [rS(cM).rDistQ, rS(cM).pDistQ, rS(cM).hProxQ, aDQ, rDQ, cDQ, ~, ~, rS(cM).rUDistQ, rS(cM).rLDistQ] = ...
        CalcDistCorr(sFP,w,Q);
    
    
    
    % Return neural activation for penultimate network
    if cM == 4
        QforPlot_MeanBody = Q_MeanBody;
        QforPlot_MeanLimb = Q_MeanLimb;
        neurActForPlot_MeanBody = permute(allNeurAct_MeanBody,[3 4 1 2 5 6]);
        neurActForPlot_MeanLimb = permute(allNeurAct_MeanLimb,[3 4 1 2 5 6]);
    end
    
    
    % Return neural activation for first only limb network
    if cM < 4
        neurActForPlot{cM} = permute(allNeurAct,[3 4 1 2 5 6]);
    else
        neurActForPlot{cM} = allNeurAct_PlusBody;
    end
    
end

%% Display proximity-dependence stats

% For all the models
[rMat] = DisplayProxStats(rS);

% For only models with multiple limbs - split the appropriate structure up
iM2 = 1; iM3 = 1;
for iM=1:length(rS)
    % load any possible extra settings
    s = rS(iM).s;
    s = DefaultSettings(s);
    
    
    if s.fl.bdyMov == 1
        
        % -----------------------------------------------------------------
        % Make a new substructure for every limb, rSML.
        % AND also make a structure which includes all models, but each limb as
        % a separate model, rSAll
        sFields = fields(rS(iM).lmb);
        for iF = 1:length(sFields)
            rSML(iM2).(sFields{iF})        = rS(iM).lmb.(sFields{iF});
        end
        rSML(iM2).s      = rS(iM).s;
        rSML(iM2).Qtable = rS(iM).Qtable;
        rSML(iM2).w      = rS(iM).w;
        rSML(iM2).net    = rS(iM).net;
        
        sFields = fields(rSML(iM2));
        for iF = 1:length(sFields)
            rSAll(iM3).(sFields{iF})        = rSML(iM2).(sFields{iF});
        end
        
        iM2 = iM2 + 1;
        iM3 = iM3 + 1;
        
        sFields = fields(rS(iM).lmb);
        for iF = 1:length(sFields)
            rSML(iM2).(sFields{iF})        = rS(iM).bdy.(sFields{iF});
        end
        rSML(iM2).s      = rS(iM).s;
        rSML(iM2).Qtable = rS(iM).Qtable;
        rSML(iM2).w      = rS(iM).w;
        rSML(iM2).net    = rS(iM).net;
        
        sFields = fields(rSML(iM2));
        for iF = 1:length(sFields)
            rSAll(iM3).(sFields{iF})        = rSML(iM2).(sFields{iF});
        end
        
        iM2 = iM2 + 1;
        iM3 = iM3 + 1;
        
    else
        rSAll(iM3) = rS(iM);
        iM3 = iM3 + 1;
    end
end

% So multiple limb results here
disp('-----------')
disp('Multiple limb results')
[rMat_ML] = DisplayProxStats(rSML);

disp('-----------')
disp('All results')
[rMat_All] = DisplayProxStats(rSAll);
% legend('limb','body')

disp('-----------')
disp('Q correlation with proximity')
disp(rMat_All.pDistQ)

disp([''])
disp(['Proportion of neurons that correlate with proximity'])
disp(rMat_All.propCorrNeur)
disp(['Av: ' num2str(nanmean(rMat_All.propCorrNeur)) '+-' num2str(nanstd(rMat_All.propCorrNeur)) ])

% also find proportion in last layer
for iM = 1:length(rSAll)
    llCorr(iM) = rMat_All.propCorrNeurPerLay(length(rSAll(iM).s.lp.netS),iM);
end
disp(['Last layer: ' num2str(nanmean(llCorr)) '+-' num2str(nanstd(llCorr)) ])

%% For RT task

% number of repetitions of fake reaction time task.
nRep = 50;

clear Q

for iM=1:length(rtRS)+1
    
    % For iM == 3, use QPPS but full of the median QPPS outpu
    if iM < 3
        cM = iM;
        QPPS2 = QPPS;
    else 
        cM = 2;
        QPPS2(:,:,:,:,:) = 0;nanmedian(QPPS(:));
    end
    
    
    sFP = DefaultSettings(rtRS(cM).s); 
    w = rtRS(cM).w;
    net = rtRS(cM).net;
    Qtable = rtRS(cM).Qtable;
    
    % settings for plot
    sFP.plt.lmbCol = 3:sFP.wrld.size(2)-2;
    sFP.plt.ON = 0;
    sFP.plt.sequentialLimbCols = 0;
    sFP.plt.stimRow = [3:size(w.world2D,1)-1];
    sFP.plt.stimCol = [2:size(w.world2D,2)-1];
    sFP.plt.pltType = 'Binned';
    
    sFP.plt.rtTouchVal = sFP.rtt.threshold;
    
    for iRep = 1:nRep
        sFP.plt.rtTouchVal = sFP.rtt.threshold+sFP.rtt.NoiseFun(sFP.rtt.mu,sFP.rtt.std);
        touchVals(iRep) = sFP.plt.rtTouchVal;
        [Q(:,:,:,:,:,iRep),allNeurAct] = CalcNetOutput(sFP,w,net,QPPS2);
    end
    
    
    qDiff = squeeze(Q(:,:,:,:,4,:)-Q(:,:,:,:,2,:)) ;
    buttonPressed = qDiff > 0;
    rtRS(iM).avBP = nanmean(buttonPressed,5) ; % Button press average
    rtRS(iM).Q = Q;
    
    % [rtRS(iM).rDistQ, rtRS(iM).pDistQ, rtRS(iM).hProxQ, aDQ, rDQ, cDQ] = ...
    %          CalcDistCorr(sFP,w,buttonPressed);
    
    % $$$ THIS WORKED!! (with 100 reps ) --> decide what to do next (obviously
    % write it in the methods, and the p-value in the main text, and make a
    % figure. BUT THEN what?
    [rtRS(iM).rDistBP, rtRS(iM).pDistBP, rtRS(iM).hProxBP, aDQ, rDQ, cDQ] = ...
        CalcDistCorr(sFP,w,rtRS(iM).avBP);
    
    [rtRS(iM).rDistQ, rtRS(iM).pDistQ, rtRS(iM).hProxQ, aDQ, rDQ, cDQ] = ...
        CalcDistCorr(sFP,w,nanmean(rtRS(iM).Q,6));
    
    %      % This returns a correlation value for each (neuron, layer, limbcol)0.05;
    %      [rtRS(iM).rDistN, rtRS(iM).pDistN, rtRS(iM).hProxN, aD, rD, cD] = ...
    %          CalcDistCorr(sFP,w,permute(allNeurAct,[3 4 1 2 5 6]));
    %      % and this for each action q-value
    %      [rtRS(iM).rDistQ, rtRS(iM).pDistQ, rtRS(iM).hProxQ, aDQ, rDQ, cDQ] = ...
    %          CalcDistCorr(sFP,w,Q);
end

[rMat_RS] = DisplayProxStats(rtRS); % $$$ NEXT CALCULATE THIS

%% $$$ I NEED TO FIGURE OUT HOW THE BODY/LIMB THING COMES INTO EFFECT HERE
% 
% iM = 2;
% 
% rtRS(iM).s.fl.extraInput
% 
% sFP = DefaultSettings(sFP);
% 
% sFP.plt.pltType = 'Imagesc';
% % sFP.plt.pltType = 'Binned';
% sFP.plt.ON = 1;
% sFP.plt.meanLimbCols = 1;
% 
% sFP.plt.bdyCol = 2;
% 
% 
% sFP.plt.stimRow = [3:size(w.world2D,1)-1];
% 
% sFP.plt.stimCol = [2:size(w.world2D,2)-1];
% sFP.plt.lmbCol = [2:size(w.world2D,2)-1];
% % sFP.plt.stimCol = 8;
% % sFP.plt.lmbCol = 8;
% % sFP.plt.stimCol = [6:10];
% % sFP.plt.lmbCol = 8;
% 
% sFP.plt.fS.pltType = 'bar';
% % sFP.plt.fS.pltType = 'shadeplot';
% sFP.plt.fS.nBins = 5;
% 
% sFP.plt.varargs = {'LineWidth',2};
% 
% sFP.plt.distanceType = 'Absolute';
% % sFP.plt.distanceType = 'Column';
% % sFP.plt.distanceType = 'AbsColumn';
% % sFP.plt.distanceType = 'AbsRow';
% % sFP.plt.distanceType = 'Row';
% 
% figure,
% sFP.plt.plAct = 1;
% DisplActValsFun(sFP,w, rtRS(iM).avBP);
% title('ButtonPressNess')
% 
% colormap(whitetocol(100,[0 0 0.7])); 
% 
% figure,
% sFP.plt.plAct = 3;
% DisplActValsFun(sFP,w, nanmean(rtRS(iM).Q,6));
% 
% colormap(whitetocol(100,[0 0 0.7])); 
% 
% 
% rtRS(iM).pDistQ.max
% rtRS(iM).rDistQ.maxP

% % % 
% % % %% For RT task VERSION 2
% % % 
% % % % number of repetitions of fake reaction time task.
% % % nRep = 50;
% % % 
% % % clear Q
% % % 
% % % for iM=1:length(rtRS)
% % %     
% % %     
% % %     sFP = DefaultSettings(rtRS(iM).s); 
% % %     w = rtRS(iM).w;
% % %     net = rtRS(iM).net;
% % %     Qtable = rtRS(iM).Qtable;
% % %     
% % %     % settings for plot
% % %     sFP.plt.lmbCol = 3:sFP.wrld.size(2)-2;
% % %     sFP.plt.ON = 0;
% % %     sFP.plt.sequentialLimbCols = 0;
% % %     sFP.plt.stimRow = [3:size(w.world2D,1)-1];
% % %     sFP.plt.stimCol = [2:size(w.world2D,2)-1];
% % %     sFP.plt.pltType = 'Binned';
% % %     
% % %     sFP.plt.rtTouchVal = sFP.rtt.threshold;
% % %     
% % %     for iRep = 1:nRep
% % %         sFP.plt.rtTouchVal = sFP.rtt.threshold+sFP.rtt.NoiseFun(sFP.rtt.mu,sFP.rtt.std);
% % %         touchVals(iRep) = sFP.plt.rtTouchVal;
% % %         [Q(:,:,:,:,:,iRep),allNeurAct] = CalcNetOutput(sFP,w,net);
% % %     end
% % %     
% % %     
% % %     qDiff = squeeze(Q(:,:,:,:,4,:)-Q(:,:,:,:,2,:)) ;
% % %     buttonPressed = qDiff > 0;
% % %     rtRS(iM).avBP = nanmean(buttonPressed,5) ; % Button press average
% % %     rtRS(iM).Q = Q;
% % %     
% % %     % [rtRS(iM).rDistQ, rtRS(iM).pDistQ, rtRS(iM).hProxQ, aDQ, rDQ, cDQ] = ...
% % %     %          CalcDistCorr(sFP,w,buttonPressed);
% % %     
% % %     % $$$ THIS WORKED!! (with 100 reps ) --> decide what to do next (obviously
% % %     % write it in the methods, and the p-value in the main text, and make a
% % %     % figure. BUT THEN what?
% % %     [rtRS(iM).rDistBP, rtRS(iM).pDistBP, rtRS(iM).hProxBP, aDQ, rDQ, cDQ] = ...
% % %         CalcDistCorr(sFP,w,rtRS(iM).avBP);
% % %     
% % %     [rtRS(iM).rDistQ, rtRS(iM).pDistQ, rtRS(iM).hProxQ, aDQ, rDQ, cDQ] = ...
% % %         CalcDistCorr(sFP,w,nanmean(rtRS(iM).Q,6));
% % %     
% % %     %      % This returns a correlation value for each (neuron, layer, limbcol)0.05;
% % %     %      [rtRS(iM).rDistN, rtRS(iM).pDistN, rtRS(iM).hProxN, aD, rD, cD] = ...
% % %     %          CalcDistCorr(sFP,w,permute(allNeurAct,[3 4 1 2 5 6]));
% % %     %      % and this for each action q-value
% % %     %      [rtRS(iM).rDistQ, rtRS(iM).pDistQ, rtRS(iM).hProxQ, aDQ, rDQ, cDQ] = ...
% % %     %          CalcDistCorr(sFP,w,Q);
% % % end
% % % 
% % % 
% % % %%
% % % 
% % % % $$$ figure out why this isn't p[lotting
% % % 
% % % % DisplActValsFun(sFP,w, nanmean(Q(:,:,:,:,:,38),6));
% % % 
% % % iM = 1;
% % % 
% % % % sFP.plt.pltType = 'Binned';
% % % sFP.plt.pltType = 'Imagesc';
% % % 
% % % % sFP.plt.fS.pltType = 'shadeplot'
% % % 
% % % sFP.plt.plAct = 1;
% % % sFP.plt.OneActFl = 1;
% % % sFP.plt.ON = 1;
% % % sFP.plt.meanLimbCols = 1;
% % % sFP.plt.lmbCol = 2:14;
% % % % sFP.plt.lmbCol = 8;
% % % sFP.plt.bdyCol = 10;
% % % sFP.plt.plotThrFl = 0;
% % % 
% % % 
% % % DisplActValsFun(sFP,w, rtRS(iM).avBP);
% % % 
% % % colormap(whitetocol(100,[0 0 0.7]));
% % % 
% % % caxis([0 1])

%%
%  figure,plot(Q(:),'x')





% $$$ HOW AM I going to show that RESPONSE RATE correlates? Need to run
% assess the network with tactile input, and see what the chosen action
% is. THEN, if the probability of choosing button press is affected by
% distance, we're good!

%  disp(['Detection correlation with


% $$$ SO TRY: use the cD etc returned by the Calcdistcorr function, to
% correlate all repetitions of Q simultaneously with cD.
% $$$ OR: run more repetitions of Q - trying that now
% $$$ NOTE: the buttonPressed on its own (not average) wasn't working very
% well because in some situations, the button isn't pressed at all..

%% Check performance of networks

% % % for iM = 1:length(rS)
% % %
% % %     s = rS(iM).s;
% % %     w = rS(iM).w;
% % %     net = rS(iM).net;
% % %     Qtable = rS(iM).Qtable;
% % %
% % %     % SET environment to simple default
% % %     s.gol.alSpR = 1;
% % %     s.gol.alSpC = [0 0];
% % %     s.gol.randSpr = [0 3];
% % %     s.lp.b1Siz = 1e4;
% % %     s.wrld.resetType = 'BottomTop_InfLR';
% % %     s.lmb.startCol=8;
% % %     s.act.GoalRew = 2;
% % %     s.act.bdyGoalRew = 2;
% % %
% % %     s = DefaultSettings(s);
% % %
% % %     s.prf.nRep = 10;
% % %
% % %     [rewPerAct(iM,:), rewPerActSD(iM,:), rewSums(:,:,iM), sE{iM}] = CalcNetPerf(s,net);
% % %
% % % end
% % % 
% % % 
% % % %% $$$ This next bit could be converted into plots if neccessary
% % % 
% % % 
% % % cM = 14;
% % % clear allNeurAct Q
% % % 
% % % rC = rSML;
% % % 
% % % s = rC(cM).s;
% % % w = rC(cM).w;
% % % net = rC(cM).net;
% % % Qtable = rC(cM).Qtable;
% % % 
% % % s = DefaultSettings(s);
% % % 
% % % sFP = s;
% % % sFP.plt.ON = 1;
% % % 
% % % sFP.plt.plAct = 5;
% % % 
% % % 
% % % 
% % % 
% % % % figure,
% % % sFP.plt.ON = 1;
% % % sFP.plt.rowLims = [1.5 s.wrld.size(1)-0.5];
% % % sFP.plt.sequentialLimbCols = 1;
% % % sFP.plt.meanLimbCols = 0;
% % % sFP.plt.lmbCol=12;
% % % sFP.plt.bdyCol=5;
% % % % sFP.plt.lmbCol = 2:s.wrld.size(2)-1;
% % % [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
% % % DisplActValsFun(sFP,w,Q);
% % % % GridOverImage(fS,gca);
% % % colormap(whitetocol(100,[0 0 0.7]));
% % % 
% % % 
% % % %%
% % % 
% % % % $$$ I NEED TO UNDERSTAND WHAT'S GOING ON WITH THE MULTIPLE LIMB COLUMNS NEXT
% % % % $$$ I NEED TO UNDERSTAND WHAT'S GOING ON WITH THE MULTIPLE LIMB COLUMNS NEXT
% % % 
% % % % figure,
% % % sFP.plt.ON = 1;
% % % sFP.plt.pltType = 'Imagesc';
% % % sFP.plt.sequentialLimbCols=1;
% % % sFP.plt.meanLimbCols=0;
% % % sFP.plt.bdyCol = 8;
% % % sFP.plt.lmbCol = 11; % $$$ MAYBE later make it work with average limb cols better
% % % % THe way to do that would be to go into the % TEST ----- section in
% % % % NetAnalysis, and find out what else needs to be changed to make it work
% % % [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
% % % NetAnalysis(sFP,w,allNeurAct);
% % % colormap(whitetocol(100,[0 0 0.7]))





% % % % BITS:
% % % for iL = 1:length(s.lp.netS)
% % %     for iN = 1:s.lpnetS(iL)
% % %         maxP
% % %     end
% % % end
