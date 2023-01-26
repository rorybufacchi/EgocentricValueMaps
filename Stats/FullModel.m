%% Assess the effect of network size on performance
addpath(genpath('D:\Old_D\Rory\Agent\Scripts'))

% load('D:\Old_D\DPPS\DefenseAgent\Results\Performance\ToolUse\NetSizes\Mults_of_12_51Batch_V2');

% load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ToolUSe\Tool_Pre_post_61_201Batch_v2');

% load('D:\Old_D\DPPS\DefenseAgent\Results\Performance\FullModel\NetSizes\Body_12_20_1Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01.mat')



% clear ntRS
% load('D:\Old_D\Rory\Agent\Results\Performance\FullModel\NetSizes\Body_pyramid1_2Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE.mat')

% % % load('D:\Old_D\Rory\Agent\Results\Performance\FullModel\NetSizes\Body_30_2Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE.mat')
% % % tmp = ntRS;
% % % 
% % % clear ntRS
% % % load('D:\Old_D\Rory\Agent\Results\Performance\FullModel\NetSizes\Body_30_30Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE.mat')
% % % 
% % % tmp(end+1)= ntRS;
% % % tmp(1).perf.rewPerAct(2,:) = tmp(2).perf.rewPerAct(1,:);
% % % 
% % % clear ntRS
% % % % $$$ THIS IS RELEARNED at least once
% % % load('D:\Old_D\Rory\Agent\Results\Performance\FullModel\NetSizes\Body_30_50Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE.mat')
% % % 
% % % tmp(end+1)= ntRS;
% % % tmp(1).perf.rewPerAct(3,:) = tmp(3).perf.rewPerAct(1,:);
% % % 
% % % clear ntRS
% % % % $$$ THIS IS RELEARNED at least once
% % % load('D:\Old_D\Rory\Agent\Results\Performance\FullModel\NetSizes\Body_30_61Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE.mat')
% % % 
% % % tmp(end+1)= ntRS;
% % % tmp(1).perf.rewPerAct(4,:) = tmp(4).perf.rewPerAct(1,:);
% % % 
% % % ntRS = tmp;



% clear ntRS
% load('D:\Old_D\Rory\Agent\Results\Performance\FullModel\NetSizes\Body_pyramid1_2Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE.mat')
% 
% tmp(end+1)= ntRS;
% tmp(1).perf.rewPerAct(2,:) = tmp(2).perf.rewPerAct(1,:);


% f = figure,
for iM = 1:1 %length(ntRS)
    
    s = ntRS(iM).s;
    w = ntRS(iM).w;
    net = ntRS(iM).net;
    
    % Settings for plot
    sFP=s;
    sFP.plt.lmbCol = 3:s.wrld.size(2)-2;
    sFP.plt.ON = 1;
    sFP.plt.rowLims = [1.5 s.wrld.size(1)-0.5];
    sFP.plt.sequentialLimbCols=0;   
    sFP.plt.pltType = 'Imagesc';
    
    sFP.plt.OneActFl = 1;
    sFP.plt.intrpFact = 1;
    
    sFP.plt.plotThrFl = 1;
    
    sFP.plt.ToolPresent = 0;
    
%     sFP.plt.meanLimbCols = 1;
%     sFP.plt.lmbCol=3:11;
    
    sFP.plt.meanLimbCols = 0;
    sFP.plt.lmbCol = 6;
    
    sFP.plt.bdyCol = 8;
    
    sFP.plt.plAct=4;
    
    sFP.plt.rowLag = 1;
    sFP.plt.colLag = 0;
    
    sFP.plt.stimRow=[3:size(w.world2D,1)-1];
    sFP.plt.stimCol=[2:size(w.world2D,2)-1];
    %     sFP.plt.pltType = 'Binned';
    
    sFP = DefaultSettings(sFP);
    
    
    sFP.plt.meanOSVfl = 0;
% % %     
% % %     % 	set(0, 'currentfigure', f);
% % %     subplot(1 + ceil(length(ntRS)/2),2,1)
% % %     tmp = s.prf.skipBatches;
% % %     plot(ntRS(iM).perf.rewPerAct(1:tmp:end,1),'LineWidth',2);
% % %     hold on;
    
    %     f2 = figure,
%     subplot(1 + ceil(length(ntRS)/2),2,1+iM)
    [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
    DisplActValsFun(sFP,w,Q);
    
%     caxis([-0.1 0.1])
    
end
% % % subplot(1 + ceil(length(ntRS)/2),2,1)
% % % % set(0, 'currentfigure', f);
% % % legend('1','2','3','4','5','6','7','8','9','10','11');
% % % % hold off



%% -------------------------------------------------------------------------
% Loop through current Model, to find limb-specific activations

clear allNeurAct_PlusBody Q_PlusBody
for iM = 1:length(rS)
    
    s = rS(iM).s;
    w = rS(iM).w;
    net = rS(iM).net;
    Qtable = rS(iM).Qtable;
    
    s = DefaultSettings(s);
    
    % settings for plot
    sFP = s;
    %     sFP.plt.lmbCol = 2:s.wrld.size(2)-1;
    sFP.plt.lmbCol = 3:s.wrld.size(2)-2;
    sFP.plt.ON = 0;
    sFP.plt.sequentialLimbCols = 0;
    sFP.plt.stimRow = [3:size(w.world2D,1)-1];
    sFP.plt.stimCol = [2:size(w.world2D,2)-1];
    sFP.plt.pltType = 'Binned';
    
    % =====================================================================
    sFP.plt.plotThrFl = 0;
    sFP.plt.ToolPresent = 1;
    sFP.plt.DistFromTool = 1;
    % =====================================================================
    
    
    if s.fl.bdyMov == 1
        clear Q_PlusBody allNeurAct_PlusBody
        for iBdyCol = 1:size(w.world2D,2)
            s.plt.bdyCol = iBdyCol;
            
            [Q_tmp, allNeurAct_tmp] = CalcNetOutput(s,w,net);
    
            Q_PlusBody(:,:,:,:,:,iBdyCol) = Q_tmp;
            allNeurAct_PlusBody(:,:,:,:,:,:,iBdyCol) = allNeurAct_tmp;
        end
        
        allNeurAct_MeanBody = nanmean(allNeurAct_PlusBody(:,:,:,:,:,:,2:size(w.world2D,2)-1),7);
        allNeurAct_MeanLimb = permute(nanmean( ...
            allNeurAct_PlusBody(:,:,:,2:size(w.world2D,2)-1,:,:,:),4),[1 2 3 7 5 6 4]);
        
        Q_MeanBody = nanmean(Q_PlusBody(:,:,:,:,:,2:size(w.world2D,2)-1),6);
        Q_MeanLimb = permute(nanmean(Q_PlusBody(:,2:size(w.world2D,2)-1,:,:,:,:),2),[1 6 3 4 5 2]);
        

        
        
        allNeurAct = allNeurAct_MeanBody;
        Q = Q_MeanBody;
        
        % This returns a correlation value for each (neuron, layer, limbcol)0.05;
        [rS(iM).bdy.rDistN, rS(iM).bdy.pDistN, rS(iM).bdy.hProxN, aD, rD, cD] = ...
            CalcDistCorr(sFP,w,permute(allNeurAct_MeanLimb,[3 4 1 2 5 6]));
        % and this for each action q-value
        [rS(iM).bdy.rDistQ, rS(iM).bdy.pDistQ, rS(iM).bdy.hProxQ, aDQ, rDQ, cDQ] = ...
            CalcDistCorr(sFP,w,Q_MeanLimb);
        
        % This returns a correlation value for each (neuron, layer, limbcol)0.05;
        [rS(iM).lmb.rDistN, rS(iM).lmb.pDistN, rS(iM).lmb.hProxN, aD, rD, cD] = ...
            CalcDistCorr(sFP,w,permute(allNeurAct_MeanBody,[3 4 1 2 5 6]));
        % and this for each action q-value
        [rS(iM).lmb.rDistQ, rS(iM).lmb.pDistQ, rS(iM).lmb.hProxQ, aDQ, rDQ, cDQ] = ...
            CalcDistCorr(sFP,w,Q_MeanBody);
    else
        [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
    end
    
    % This returns a correlation value for each (neuron, layer, limbcol)0.05;
    [rS(iM).rDistN, rS(iM).pDistN, rS(iM).hProxN, aD, rD, cD] = ...
        CalcDistCorr(sFP,w,permute(allNeurAct,[3 4 1 2 5 6]));
    % and this for each action q-value
    [rS(iM).rDistQ, rS(iM).pDistQ, rS(iM).hProxQ, aDQ, rDQ, cDQ] = ...
        CalcDistCorr(sFP,w,Q);
    
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

[rMat_ML] = DisplayProxStats(rSML);

[rMat_All] = DisplayProxStats(rSAll);
% legend('limb','body')

disp('Q correlation with proximity')
disp(rMat_All.pDistQ)

disp([''])
disp(['Proportion of neurons that correlate with proximity'])
disp(rMat_All.propCorrNeur)
disp(['Av: ' num2str(nanmean(rMat_All.propCorrNeur)) '+-' num2str(nanstd(rMat_All.propCorrNeur)) ])

% also find proportion in last layer
for iM = 1:length(rSML)
    llCorr(iM) = rMat_All.propCorrNeurPerLay(length(rSML(iM).s.lp.netS),iM);
end
disp(['Last layer: ' num2str(nanmean(llCorr)) '+-' num2str(nanstd(llCorr)) ])
