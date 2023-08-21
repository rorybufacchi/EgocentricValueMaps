%% Effects of policy

% % % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\SARSA_WithUnbalanced_AndThreat.mat')
% % % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\SARSA_WithUnbalanced.mat')

f.Policy.f = figure('Position',[20 20 1800 900]);

sgtitle('POSITIVE VALENCE STIMULI                                                              NEGATIVE VALENCE STIMULI')


fS.lineYLim = [-0.01 0.55];
fS.lineXLim = [ 0   0.7];
fS.xPoints  = [5 5 ; 10 10]';
fS.cAxes    = [-1.8 1.8];

fS.gridXstart   = -5.5;
fS.gridXstep    =  1;
fS.gridYstart   =  2.5;
fS.gridYstep    =  1;

iPl = 0; iAx = 0;
for iD = 1:size(rS,2)

    % =====================================================================
    % Positive valence stimuli

    % -----------------------------------
    % Q-learning
    iPl = iPl + 1;
    iAx = iAx + 1;
    f.Policy.ax{iAx} = subplot(3,6,iPl);
    for iM = 1:3
        [tmpQ(:,:,:,:,:,iM),allNeurAct] = CalcNetOutput(rS(iM,iD).s,rS(iM,iD).w,rS(iM,iD).net);
    end
    Q = nanmean(tmpQ,6);
    [ax, plQhoriz(:,1,1,iD), plQvert(:,1,1,iD)] = plot2D(rS,fS,Q,iM,iD);
    axes(ax{1});
    caxis(fS.cAxes); %caxis([0 2])
    title(['Q learning, policy update nr:' num2str(iD)]);

    % -----------------------------------
    % SARSA
    iPl = iPl+1;
    iAx = iAx + 1;
    f.Policy.ax{iAx} = subplot(3,6,iPl);
    for iM = 1:3
        [tmpQ(:,:,:,:,:,iM),allNeurAct] = CalcNetOutput(rSsarsa(iM,iD).s,rSsarsa(iM,iD).w,rSsarsa(iM,iD).net);
    end
    Q = nanmean(tmpQ,6);
    [ax, plQhoriz(:,1,2,iD), plQvert(:,1,2,iD)] = plot2D(rSsarsa,fS,Q,iM,iD);
    axes(ax{1});
    caxis(fS.cAxes); %caxis([0 2])
    title(['SARSA, policy update nr:' num2str(iD)]);
   
    % -----------------------------------
    % SARSA with biased-direction epsilon-greedy policy
    iPl = iPl+1;
    iAx = iAx + 1;
    f.Policy.ax{iAx} = subplot(3,6,iPl);
    for iM = 1:3
        [tmpQ(:,:,:,:,:,iM),allNeurAct] = CalcNetOutput(rSsarsaUnb(iM,iD).s,rSsarsaUnb(iM,iD).w,rSsarsaUnb(iM,iD).net);
    end
    Q = nanmean(tmpQ,6);
    [ax, plQhoriz(:,1,3,iD), plQvert(:,1,3,iD)] = plot2D(rSsarsaUnb,fS,Q,iM,iD);
    axes(ax{1});
    caxis(fS.cAxes); %caxis([0 2])
    title(['UNbalanced SARSA, policy update nr:' num2str(iD)]);


    % =====================================================================
    % Negative valence stimuli


    % -----------------------------------
    % Q-learning
    iPl = iPl + 1;
    iAx = iAx + 1;
    f.Policy.ax{iAx} = subplot(3,6,iPl);
    for iM = 1:3
        [tmpQ(:,:,:,:,:,iM),allNeurAct] = CalcNetOutput(rSthr(iM,iD).s,rSthr(iM,iD).w,rSthr(iM,iD).net);
    end
    Q = nanmean(tmpQ,6);
    [ax, plQhoriz(:,2,1,iD), plQvert(:,2,1,iD)] = plot2D(rSthr,fS,Q,iM,iD);
    axes(ax{3}); xlim(-flip(fS.lineXLim));
    axes(ax{2}); ylim(-flip(fS.lineYLim));
    axes(ax{1}); 
    caxis(fS.cAxes); %caxis([-2 0])
    title(['Q learning, policy update nr:' num2str(iD)]);


    % -----------------------------------
    % SARSA
    iPl = iPl+1;
    iAx = iAx + 1;
    f.Policy.ax{iAx} = subplot(3,6,iPl);
    for iM = 1:3
        [tmpQ(:,:,:,:,:,iM),allNeurAct] = CalcNetOutput(rSsarsathr(iM,iD).s,rSsarsathr(iM,iD).w,rSsarsathr(iM,iD).net);
    end
    Q = nanmean(tmpQ,6);
    [ax, plQhoriz(:,2,2,iD), plQvert(:,2,2,iD)] = plot2D(rSsarsathr,fS,Q,iM,iD);
    axes(ax{3}); xlim(-flip(fS.lineXLim));
    axes(ax{2}); ylim(-flip(fS.lineYLim));
    axes(ax{1}); 
    caxis(fS.cAxes); %caxis([-2 0])
    title(['SARSA, policy update nr:' num2str(iD)]);


   
    % -----------------------------------
    % SARSA with biased-direction epsilon-greedy policy
    iPl = iPl+1;
    iAx = iAx + 1;
    f.Policy.ax{iAx} = subplot(3,6,iPl);
    for iM = 1:3
        [tmpQ(:,:,:,:,:,iM),allNeurAct] = CalcNetOutput(rSsarsaUnbthr(iM,iD).s,rSsarsaUnbthr(iM,iD).w,rSsarsaUnbthr(iM,iD).net);
    end
    Q = nanmean(tmpQ,6);
    [ax, plQhoriz(:,2,3,iD), plQvert(:,2,3,iD)] = plot2D(rSsarsaUnbthr,fS,Q,iM,iD);
    axes(ax{3}); xlim(-flip(fS.lineXLim));
    axes(ax{2}); ylim(-flip(fS.lineYLim));
    axes(ax{1}); 
    caxis(fS.cAxes); %caxis([-2 0])
    title(['UNbalanced SARSA, policy update nr:' num2str(iD)]);



%     colormap(flip(redbluecmapRory))
%     colormap(parula(9))
%     colormap(hsv)
    colormap(flip(redbluecmapRory(11,5)))
end


% Plot the horizontal and vertical means in the same plots, to show effects

% pos, stimtype, algorithm, depth

vertPos = (size(plQvert, 1):-1:1) - 2;
horzPos = (size(plQhoriz,1):-1:1) - 7;

cDepths = 2;

algNames = {'Q-learning (optimisitc)', 'SARSA', 'left-biased SARSA'};
trainedOn = {'Initial random movement','1st policy-dependent movement','2nd policy-dependent movement'}

lineSettings = {'LineWidth',2};

f.PolicyLinesOptimism.f = figure('Position',[20 50 1200 1200]);

sgtitle('Q-learning gives more optimistic fields')

subplot(2,2,1)
plot(vertPos,  squeeze(  nanmean(plQvert(:,1,:,cDepths),4 )  ),lineSettings{:}); % Average over depth
title('Vertical distance; Positive stim: optimism results in bigger fields')
xlabel('Stim Dist')
ylabel('Q-value')
xlim([1 10]);
% ylim([-0.6 0.6]);
ylim([0 0.6]);
legend(algNames{:});
box off; hold on; grid on
plot([1 10],[0 0],'-k')

subplot(2,2,3)
plot(vertPos,  squeeze(  nanmean(plQvert(:,2,:,cDepths),4 )  ),lineSettings{:}); % Average over depth
title('Vertical distance; Negative stim: optimism results in smaller fields')
xlabel('Stim Dist')
ylabel('Q-value')
xlim([1 10]);
% ylim([-0.6 0.6]);
ylim([-0.6 0]);
box off; hold on; grid on
plot([1 10],[0 0],'-k')

subplot(2,2,2)
plot(horzPos,  squeeze(  nanmean(plQhoriz(:,1,:,cDepths),4 )  ),lineSettings{:}); % Average over depth
title('Horizontal distance; Positive stim: optimism results in bigger fields')
xlabel('Stim Dist')
ylabel('Q-value')
% ylim([-0.6 0.6]);
ylim([0 0.6]);
box off; hold on; grid on
plot([-6 6],[0 0],'-k')

subplot(2,2,4)
plot(horzPos,  squeeze(  nanmean(plQhoriz(:,2,:,cDepths),4 )  ),lineSettings{:}); % Average over depth
title('Horizontal distance; Negative stim: optimism results in smaller fields')
xlabel('Stim Dist')
ylabel('Q-value')
% ylim([-0.6 0.6]);
ylim([-0.6 0]);
box off; hold on; grid on
plot([-6 6],[0 0],'-k')


% LEARNING EFFECTS

f.PolicyLinesLearning.f = figure('Position',[20 50 1200 600]);


% sgtitle('Effects of training and experience: POSITIVE VALENCE STIMULI                                                              NEGATIVE VALENCE STIMULI')
sgtitle('Effects of training and experience')

for iP = 1:length(algNames)

%     subplot(2,6,iP)
    subplot(2,6,iP)
%     subplot(3,4,(iP-1) .*4 + 1)
    plot(vertPos,  squeeze(  plQvert(:,1,iP,:)  ),lineSettings{:}); 
    title( ['Vertical, ' algNames{iP} ])
    xlabel('Stim Dist')
    ylabel('Q-value: GOAL')
    box off
    grid on
%     ylim([-.6 .6])
    ylim([0 .6])
    xlim([1 10]);

% % %     axx = gca;
% % %     axxx = axes('Position',axx.Position .* [1 1 1 0.2] + [0, -axx.Position(4).*0.2, 0, 0])  
% % %     bar(squeeze(  nanmean(plQvert(:,:,iP,:),1)))
% % %     title('mean Q')


%     subplot(2,6,iP + 3)
    subplot(2,6,iP + 6)
%     subplot(3,4,(iP-1) .*4 + 1)
    plot(vertPos,  squeeze(  plQvert(:,2,iP,:)  ),lineSettings{:}); 
    title( ['Vertical, ' algNames{iP} ])
    xlabel('Stim Dist')
    ylabel('Q-value: THREAT')
    box off
    grid on
%     ylim([-.6 .6])
    ylim([-.6 0])
    xlim([1 10]);


%     subplot(2,6,iP + 6)
    subplot(2,6,iP + 3)
    plot(horzPos,  squeeze(  plQhoriz(:,1,iP,:)  ),lineSettings{:}); 
    title( ['Horizontal, ' algNames{iP} ])
    xlabel('Stim Dist')
    ylabel('Q-value: GOAL')
    box off
    grid on
%     ylim([-.6 .6])
    ylim([0 .6])

    subplot(2,6,iP + 9)
    plot(horzPos,  squeeze(  plQhoriz(:,2,iP,:)  ),lineSettings{:}); 
    title( ['Vertical, ' algNames{iP} ])
    xlabel('Stim Dist')
    ylabel('Q-value: THREAT')
    box off
    grid on
%     ylim([-.6 .6])
    ylim([-.6 0])


end
legend(trainedOn{:});


f.PolicyBarsLearning.f = figure('Position',[20 50 1200 600]);

% pos, stimtype, algorithm, depth
avQval = squeeze(nanmean(plQvert) + nanmean(plQhoriz));

for iP = 1:length(algNames)

    subplot(1,3,iP)
    bar(squeeze(avQval(:,iP,:)))
    set(gca,'XTickLabel',{'Goal','Threat'})
    ylabel('Average Q')
    title(algNames{iP})
    ylim([-0.5 0.5])
    grid on

end

legend(trainedOn)


%% Wasp scenario

% $$$ MAKE THIS into a loop and plot the best action for all 3
% 'environments' (columns) and hand-body distances (rows)

load('The_BEES_plus_AWAYoption_v2.mat')
% % % load('The_BEES_plus_AWAYoption_handCanMoveFaster.mat')

nDynam = 2;
allHP = -(0:3:6); 
nHandPos = numel(allHP)
iPl = 0;
% cAx = [-max(newQ2(:)) -1];
cAx = [-6 -1];

f.TwoDimPlot.f        = figure('Position',[20 20 1200 1200]);

for iHP = 1:nHandPos

    cHP = allHP(iHP);
    for iDynam = 1:nDynam

        iPl = iPl + 1

        newQ = allQ.qVals{iDynam};
        sTmp = sBdy;
        iAct = 1:size(newQ,1);

        [newQ2 optAct] = max(newQ(iAct,:,:,:),[],1); % Max value
        % newQ2 = mean(newQ(iAct,:,:,:),1); % avg value

        newQ2 = optAct;


        % % % f.TwoDimPlot.ax{1}    = axes('Position',[.05 .1 .4 .8]);
        % % % f.TwoDimPlot.ax{1}    = axes('Position',[.05 .05 .9 .9]);
        f.TwoDimPlot.ax{iPl}    = subplot(nHandPos,nDynam,iPl);


        zPos = sTmp.clc.nearPos(3) - cHP; % POSITION OF SLICE
        imagesc(squeeze(-newQ2(1,:,:,zPos )) );
%         GridOverImage(fS,f.TwoDimPlot.ax{1});

        % Show where the hand and the body are
        hold on
        plSymb = {'xk','ok'};
        inSlice = find(sBdy.clc.startSZ == zPos);
        for cVol = inSlice
            for iSplitInd = 1:length(sBdy.clc.rewSplitInd)
                if cVol >= sBdy.clc.rewSplitInd(iSplitInd) & cVol < sBdy.clc.rewSplitInd(iSplitInd+1)
                    plot(sBdy.clc.startSC(cVol),sBdy.clc.startSR(cVol),plSymb{iSplitInd},'LineWidth',2)
                end
            end
        end
        hold off;

        % caxis auto
        caxis(cAx);


% % %         xlim([20  41]);
% % %         ylim([40  61]);

% % %         xlim([23  38]);
% % %         ylim([43  58]);

        xlim([12 49]);
        ylim([22 60]);

    end
end


f.TwoDimPlot.ax{iPl+1} = axes('Position', [0.9 0.05 0.1 0.9])
colorbar;
caxis(cAx);
f.TwoDimPlot.ax{iPl+1}.Visible = 'off';

% ------------------------------------------
% Create a custom colormap with 'stay still' as white
custom_colormap = jet;
% Replace the color for the specific value with white
custom_colormap(end, :) = [1, 1, 1];
colormap(custom_colormap);
% ------------------------------------------

%% Effects of neuron type

iM = 1;

neurTypes   = {'tansig','logsig','softmax','poslin','purelin','tribas','radbas'};
% % % neurTypes   = {'softmax','logsig','tribas','radbas','poslin','purelin','tansig'};
regTypes    = {'Valence_51','L1'};
regNames    = {'No','L1'};



for iTyp = 1:length(neurTypes);
for iReg = 1:length(regTypes);


bFold = 'F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\';
cFiles = dir([bFold '*' regTypes{iReg} '*_' neurTypes{iTyp} '_B_V*.mat'])


for iRun = 1:length(cFiles)

load([bFold cFiles(iRun).name]);

rSall(1:numel(rS),iRun,iTyp,iReg) = rS(:);

% % % end

% F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_L1_Regularization_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_softmax_B_V4.mat'


% Loop through all results and calculate structure metrics
% % % for iRun = 1:size(rSall,2)
    for iM = 1:size(rSall,1)
        try
    rSall(iM,iRun,iTyp,iReg).s.plt.nPm = 100 ; % change to 100 later

    [rSall(iM,iRun,iTyp,iReg),allNetAR(iM,iRun,iTyp,iReg),allSumBinNPFdist(iM,iRun,iTyp,iReg),allSumBinNPFpmdist(:,iM,iRun,iTyp,iReg),allBinNPF(:,:,iM,iRun,iTyp,iReg),allBinNPFpm(:,:,:,iM,iRun,iTyp,iReg)] = AssessStructure(rSall(iM,iRun,iTyp,iReg));


% % %         [rSallTMP,allNetARTMP,allSumBinNPFdistTMP,allSumBinNPFpmdistTMP,allBinNPFTMP,allBinNPFpmTMP] = AssessStructure(rSall(iM,iRun,iTyp,iReg));

        catch
            if iTyp > 1
                warning('whoops:')
                iM
                iRun
                iTyp
                iReg
            end
        end
% % %     [TMPrSall,TMPallNetAR,TMPallSumBinNPFdist,TMPallSumBinNPFpmdist,TMPallBinNPF,TMPallBinNPFpm] = ...
% % %         AssessStructure(rSall(iRun));
    end

end

end
end

% Save analysed neuron types
save('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\NeurTypesAnalysed_V3.mat','-v7.3')

%% $$$ Plot some kind of order metrics or something

% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\NeurTypesAnalysed_V2.mat')

% binEdges = -4:.5:10;
binEdges    = -7:.5:7;
cOrder      = {[0 0 .7], [.7 0 0], [.8 .6 0]};

f.NeurTypes.f = figure('Position',[20 -20 600 9000]);

iPl = 1;

% % (iM,iRun,iTyp,iReg)

for cTyp = 1:7
for cReg = 1:2
subplot(numel(neurTypes),numel(regNames),iPl)

% % % for cReg = 1:2
% % % for cTyp = 1:7
% % % subplot(numel(regNames),numel(neurTypes),iPl)

iPl = iPl + 1;



cM   = 1:3;

rsTmp = rSall(cM,:,cTyp,cReg); rsTmp = rsTmp(:);
netARTmp = allNetAR(cM,:,cTyp,cReg); netARTmp = netARTmp(:);
binNPFtmp = allBinNPF(:,:,cM,:,cTyp,cReg); binNPFtmp = binNPFtmp(:,:,:);
binNPFpmtmp = allBinNPFpm(:,:,:,cM,:,cTyp,cReg); binNPFpmtmp = binNPFpmtmp(:,:,:,:);

inclRuns = arrayfun(@(iRun) ~isempty(netARTmp(iRun).G) , 1:size(netARTmp,1)) ;


% % In case not everything has finished calculating, only include a few
% % % inclRuns = min( [ find(arrayfun(@(iRun) ~isempty(rSall(1,iRun,cTyp,cReg).w) , 1:size(rsTmp,2)),1,"last")  ... 
% % %                   find(arrayfun(@(iRun) ~isempty(rSall(3,iRun,cTyp,cReg).w) , 1:size(rsTmp,2)),1,"last") ]) ;

rsTmp       = rsTmp(inclRuns);
netARTmp    = netARTmp(inclRuns);
binNPFtmp   = binNPFtmp(:,:,inclRuns);
binNPFpmtmp = binNPFpmtmp(:,:,:,inclRuns);




if sum(inclRuns) > 0

tmpP = arrayfun(@(x) netARTmp(x).difdifP , [1:numel(netARTmp)] );
tmpR = arrayfun(@(x) netARTmp(x).difdifR , [1:numel(netARTmp)] ); 

[dmy corrP] = fdr(tmpP);

disp('mean rho')
nanmean(tmpR);
disp('std rho')
nanstd(tmpR);


disp('--------------')
sprintf('%i out of %i networks show structure using this metric', sum(corrP < .05), numel(corrP))

% % % figure,
% % % histogram(tmpP,20)


structMetric    = squeeze(sum(abs(binNPFtmp(:,:,:)),[1 2]));
structMetricPm  = squeeze(sum(abs(binNPFpmtmp(:,:,:,:)),[1 2]));

clear numSmaller
for iM = 1:size(binNPFtmp,3)
    numSmaller(iM) = sum(structMetricPm(:,iM) > structMetric(iM))
end


structureOverChance = (structMetric - structMetricPm'   );

% % % figure,histogram(structureOverChance(:))


[hh pp] = ttest(structureOverChance(:))

[hhh ppp ci statsss] =ttest2(structMetric,structMetricPm(:))
disp('struct metric mean')
nanmean(structMetric)

disp('permuted struct metric mean')
nanmean(structMetricPm(:))


% % % figure,histogram(structMetric,10,'Normalization','probability'); hold on
% % % histogram(structMetricPm(:),10,'Normalization','probability');

% for large network
allPerf = arrayfun( @(im) sum(rsTmp(im).perf.rewPerAct(end,:)) , [1:numel(rsTmp)] );
allNetW = arrayfun( @(im) rsTmp(im).s.lp.netS(end) , [1:numel(rsTmp)] );



netShapes = unique(allNetW);
nNS = length(netShapes);




for iNS = 1:nNS
    
    inclM = allNetW == netShapes(iNS);
    
    tmpPM = structureOverChance(inclM,:);
    h = histogram(tmpPM(:),'Normalization','probability','BinEdges',binEdges,...
        'FaceColor',cOrder{iNS},'FaceAlpha',0.3, ...
        'EdgeColor','k','EdgeAlpha', 0.2 ); hold on
    tmpPM2(1:numel(tmpPM(:)),iNS) = tmpPM(:);    
% % %     xlabel('Increase in network order over chance')
    hold on
    splineH{iNS} = PlotHistSpline(h,'LineWidth',2,'Color',cOrder{iNS});
    xlim([binEdges(1) binEdges(end)])
end

yLims = ylim;
plot([0 0],yLims,'-.k');
ylim([0 yLims(2)]);


end

title([neurTypes{cTyp} ', ' regNames{cReg} ' Reg'])


end
end

legend([splineH{:}],'Narrowing networks', 'constant width networks', 'Widening networks')


%% Plot eating version

iM = 1;

f.Eating.f = figure('Position',[50 200 1600 400]);

eatRS(iM).s.plt.meanLimbCols  = 1
eatRS(iM).s.plt.holdingRew    = 0;
eatRS(iM).s.plt.lmbCol        = 10;
eatRS(iM).s.plt.plAct         = 2;

eatRS(iM).s.plt.colLims = [1.5 13.5];
eatRS(iM).s.plt.rowLims = [1.5 13.5];
fS.gridXstart = -5.5;
fS.gridXstep = 1;
fS.gridYstart = 2.5;
fS.gridYstep = 1;


subplot(2,9,[1 2 10 11])
eatRS(iM).s.plt.lmbCol  = 2:14;
[Q,allNeurAct] = CalcNetOutput(eatRS(iM).s,eatRS(iM).w,eatRS(iM).net);
DisplActValsFun(eatRS(iM).s,eatRS(iM).w,-Q);
colorbar off
GridOverImage(fS,gca)
caxis([-1 1]);
title('Positive stimulus, average Q-field')

mouthOffset = [-4:2:4]; % [-2:2] ; 
for iO = 1:length(mouthOffset) % loop through offsets

    % current offset
    cLimbCol = 8 - mouthOffset(iO);
    subplot(2,9,2 + iO)
    eatRS(iM).s.plt.lmbCol  = cLimbCol;
    [Q,allNeurAct] = CalcNetOutput(eatRS(iM).s,eatRS(iM).w,eatRS(iM).net);
    DisplActValsFun(eatRS(iM).s,eatRS(iM).w,-Q);
    colorbar off
    GridOverImage(fS,gca)
    caxis([-1 1]);
    title(['REWARD: Mouth at ' num2str(mouthOffset(iO)) ])
end



colormap(redbluecmapRory)

% Shield version --> stationary dangerzone


iM = 1;

% % % f.Eating.f = figure('Position',[50 200 1600 400]);y

shieldRS(iM).s.plt.meanLimbCols  = 1
shieldRS(iM).s.plt.holdingRew    = 0;
shieldRS(iM).s.plt.lmbCol        = 10;
shieldRS(iM).s.plt.plAct         = 2;

shieldRS(iM).s.plt.colLims = [1.5 13.5];
shieldRS(iM).s.plt.rowLims = [1.5 13.5];
fS.gridXstart = -5.5;
fS.gridXstep = 1;
fS.gridYstart = 2.5;
fS.gridYstep = 1;


subplot(2,9,[8 9 17 18])
shieldRS(iM).s.plt.lmbCol  = 2:14;
[Q,allNeurAct] = CalcNetOutput(shieldRS(iM).s,shieldRS(iM).w,shieldRS(iM).net);
DisplActValsFun(shieldRS(iM).s,shieldRS(iM).w,-Q);
colorbar off
GridOverImage(fS,gca)
caxis([-.1 .1]);
title('Negative stimulus, only harms mouth, average Q-field')

mouthOffset = [-4:2:4]; % [-2:2] ; 
for iO = 1:length(mouthOffset) % loop through offsets

    % current offset
    cLimbCol = 8 - mouthOffset(iO);
    subplot(2,9,11 + iO)
    shieldRS(iM).s.plt.lmbCol  = cLimbCol;
    [Q,allNeurAct] = CalcNetOutput(shieldRS(iM).s,shieldRS(iM).w,shieldRS(iM).net);
    DisplActValsFun(shieldRS(iM).s,shieldRS(iM).w,-Q);
    colorbar off
    GridOverImage(fS,gca)
    caxis([-.4 .4]);
    title(['DANGER: Mouth at ' num2str(mouthOffset(iO)) ])
end



colormap(redbluecmapRory)

%% $$$ Shield version --> moving body
% $$$ I need to see whether the body-moving shield version works.

M = 1;

f.Eating.f = figure('Position',[50 200 1600 400]);

ntRS(iM).s.plt.meanLimbCols  = 1;
ntRS(iM).s.plt.holdingRew    = 0;
ntRS(iM).s.plt.lmbCol        = 10;
ntRS(iM).s.plt.plAct         = 3;

ntRS(iM).s.plt.colLims = [1.5 13.5];
ntRS(iM).s.plt.rowLims = [1.5 13.5];
fS.gridXstart = -5.5;
fS.gridXstep = 1;
fS.gridYstart = 2.5;
fS.gridYstep = 1;


subplot(2,9,[1 2 10 11])
ntRS(iM).s.plt.lmbCol  = 2:14;

bdyCols = 2:14;
clear qTmp
for iC = 1:length(bdyCols)

    ntRS(iM).s.plt.bdyCol = bdyCols(iC);
    [qTmp(:,:,:,:,:,iC),allNeurAct] = CalcNetOutput(ntRS(iM).s,ntRS(iM).w,ntRS(iM).net);
end
Q = nanmean(qTmp,6);

% $$$ IS THE DIFFERENCE that the limb's position changes the body's
% receptive field?...
ntRS(iM).s.plt.plAct         = 5;
DisplActValsFun(ntRS(iM).s,ntRS(iM).w,Q);
% % % DisplActValsFun(ntRS(iM).s,ntRS(iM).w,sum(Q,5)./3);

% % % colorbar off
GridOverImage(fS,gca)
% caxis([-2 2]);
% % % caxis([-0.01 0.01])
title('Positive stimulus, average Q-field')

mouthOffset = [-4:2:4]; % [-2:2] ; 

%%
ntRS(iM).s.plt.meanLimbCols  = 0;
for iO = 1:length(mouthOffset) % loop through offsets

    % current offset
    cBdyCol = 8 - mouthOffset(iO);
    ntRS(iM).s.plt.lmbCol  = 8;
    ntRS(iM).s.plt.bdyCol = cBdyCol;

    subplot(2,9,2 + iO)
    [Q,allNeurAct] = CalcNetOutput(ntRS(iM).s,ntRS(iM).w,ntRS(iM).net);
    DisplActValsFun(ntRS(iM).s,ntRS(iM).w,Q);
% % %     DisplActValsFun(ntRS(iM).s,ntRS(iM).w,sum(Q,5)./3);
    colorbar off
    GridOverImage(fS,gca)
    caxis([-2 2]);
%     caxis([-0.01 0.01])
    title(['Mouth at ' num2str(mouthOffset(iO)) ])
end



colormap(redbluecmapRory)



%% $$$ Incremental tool improvements

iM = 1;

% plotReps = 1:length(incrToolRS);
plotReps = 1:51;


% toolPerfs = cell2mat(arrayfun(@(iR) incrToolRS(iM,iR).perf.rewPerAct(end,1), plotReps, 'UniformOutput', false)');

% toolPerfs = cell2mat(arrayfun(@(iR) [incrToolRS(1,iR).perf.rewPerAct(end,1) ... 
%                                      incrToolRS(2,iR).perf.rewPerAct(end,1) ... 
%                                      incrToolRS(3,iR).perf.rewPerAct(end,1)], ...
%                                      plotReps, 'UniformOutput', false)');

toolPerfs = cell2mat(arrayfun(@(iR) [incrToolRS(1,iR).perf.rewPerAct(end,1) ... 
                                     incrToolRS(2,iR).perf.rewPerAct(end,1) ... 
                                     incrToolRS(3,iR).perf.rewPerAct(end,1)], ...
                                     plotReps, 'UniformOutput', false)');

figure,plot(toolPerfs,'LineWidth',2)      
ylabel('Performance (score/timestep)')
xlabel('Batch number')

legend('Widening','Constant','Narrowing')

%% $$$ Find the number of peaks for each training step

% $$$ Something is going wrong with iD beyond 23?...
tic
for iM = 1:size(incrToolRS,1)
    for iD = 1:size(incrToolRS,2)

        s = incrToolRS(iM,iD).s;
        w = incrToolRS(iM,iD).w;

        Qtable = incrToolRS(iM,iD).Qtable;

        s = DefaultSettings(s);
        s.plt.plotThrFl = 0;

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
        % Calculate activity WRT limb location
        sFP.plt.DistFromTool = 0;
        sFP.plt.stimRow = [size(w.world2D,1)-6:size(w.world2D,1)-3]; 
        % ---------------------------------------------------------------------
        net = incrToolRS(iM,iD).net;
        sFP.plt.ToolPresent = 1;
        [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
        [incrToolRS(iM,iD).rDistN, incrToolRS(iM,iD).pDistN, incrToolRS(iM,iD).hProxN, aD, rD, cD]    = CalcDistCorr(sFP,w,permute(allNeurAct,[3 4 1 2 5 6]));
        [incrToolRS(iM,iD).rDistQ, incrToolRS(iM,iD).pDistQ, incrToolRS(iM,iD).hProxQ, aDQ, rDQ, cDQ] = CalcDistCorr(sFP,w,Q);
        [incrToolRS(iM,iD).npks,   incrToolRS(iM,iD).pks,    incrToolRS(iM,iD).pkLocs]                = RowPeakFind(sFP,Q);
        Qall(:,:,:,:,:,iM,iD,1) = Q;


        for iN = 1:size(allNeurAct,6)
            for iL = 1:size(allNeurAct,5)
                npksN(:,iL,iN,iM,iD) = RowPeakFind(sFP, permute(allNeurAct,[3 4 1 2 5 6]));
            end
        end

        % Also calculate values for no tool
        sFP.plt.ToolPresent = 0;
        [QnoTool,allNeurAct] = CalcNetOutput(sFP,w,net);
        QallNoTool(:,:,:,:,:,iM,iD,1) = QnoTool;

        for iN = 1:size(allNeurAct,6)
            for iL = 1:size(allNeurAct,5)
                npksNnoTool(:,iL,iN,iM,iD) = RowPeakFind(sFP, permute(allNeurAct,[3 4 1 2 5 6]));
            end
        end
        
    end
end
toc



% % % % Calculate number of neuron peaks with particular number of peaks in each layer
% % % tic
% % % 
% % %     for iD = 1:size(nAall,8)
% % %     for iM = 1:size(nAall,7)
% % %         for iN = 1:size(nAall,6)
% % %         for iL = 1:size(nAall,5)
% % %             tmpNA = nAall(:,:,:,:,iL,iN,iM,iD);
% % %             if min(isnan(tmpNA(:))) == 1
% % %                 npksN(:,iL,iN,iM,iD) = NaN;
% % %             else
% % % npksN(:,iL,iN,iM,iD)  = RowPeakFind(sFP,nAall(:,:,:,:,iL,iN,iM,iD));
% % %             end
% % %         end
% % %         end
% % %     end
% % %     end
% % % 
% % % toc

save('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Tool_ReviewerResponse_Pre_post_51_51Batch_ToolPos5_moreRandSpr_NoHist_PEAKSANALYSED.mat','-v7.3')

%% Plot tool results


fS.lineYLim = [-0.01 0.55];
fS.lineXLim = [ 0   0.7];
fS.xPoints  = [5 5 ; 10 10]';
fS.cAxes    = [-1.8 1.8];

fS.gridXstart   = -5.5;
fS.gridXstep    =  1;
fS.gridYstart   =  2.5;
fS.gridYstep    =  1;




f.Tool.f = figure('Position',[20 20 1200 600]);

% plot performance
yyaxis right

mnPrf = mean(toolPerfs,2);
sdPrf = std(toolPerfs,0,2);
plot(mnPrf','r-o','LineWidth',2); hold on
opts.c = [.7 0 0];
ShadedPlot((1:51),mnPrf,mnPrf - sdPrf, mnPrf + sdPrf,opts);

ylim([0 .3])

ylabel('Reward per timestep')

% iM = 1; iD = 23;
clear allPks
for iM = 1:size(incrToolRS,1)
    for iD  = 1:size(incrToolRS,2)
        allPks(iD,iM,:)     = incrToolRS(iM,iD).npks;
        mnPks(iM,iD)        = nanmean(incrToolRS(iM,iD).npks);
        sdPks(iM,iD)        = nanstd(incrToolRS(iM,iD).npks);
        prcTilePks(:,iM,iD)   = prctile(incrToolRS(iM,iD).npks,[5 95]);
    end 
end
allPks = allPks(:,:);

% Average across model architectures
mnPks = nanmean(mnPks);
sdPks = sum(sqrt(sdPks .^2 ));


% subplot(4,3,[1:9]);

yyaxis left

[N C] = hist3( [repmat(  (1:size(allPks,1))', [size(allPks,2) 1]), allPks(:) ],[51 3]);

% imagesc(N' ./ max(N(:))); axis xy; hold on

f.Tool.ax{1} = gca;

plot(mnPks','b-o','LineWidth',2); hold on
opts.c = [0 0 .7];
opts.PlotMean = 0;
ShadedPlot((1:51),mnPks,mnPks - sdPks,mnPks + sdPks,opts);


% ylim([0 3.5])
ylim([0.8 3])
xlabel('Training batch since exposure to tool')
ylabel('Number of peaks in receptive field')


title('A receptive field grows around the tooltip as the agent learns to use the tool')

xlim([0 21])

box off



% Plot receptive fields
sFP.plt.pltType         = 'Imagesc';
sFP.plt.ON              = 1;
sFP.plt.meanLimbCols    = 1;
sFP.plt.rowLims         = [1.5 13.5];


% Plot the tool-fields as the agent learns
plBatchs = [1 2 3 5 10 20];
% Space out the sub-plots equally along the x-axis
% subPlotxPos = f.Tool.ax{1}.Position(3) ./ (numel(plBatchs) + 2);  
subPlotxPos = f.Tool.ax{1}.Position(3) ./ (numel(plBatchs) + 1 );  
cAx = 2;
for iPl = 1:length(plBatchs)

%     subplot(1,length(plBatchs),iPl);
%     subplot(4,3,9 + iPl);

    % Make an inset plot
    f.Tool.ax{cAx} = axes('Position', ... 
        f.Tool.ax{1}.Position .* [1, 1, 0, 0.2] + [ iPl .* subPlotxPos - subPlotxPos./2 ,  f.Tool.ax{1}.Position(4) .* 0.75 , subPlotxPos .* 0.8,  0]  );
    
    DisplActValsFun(sFP,incrToolRS(end,end).w,mean(Qall(:,:,:,:,:,:,plBatchs(iPl)),6));
%     DisplActValsFun(sFP,incrToolRS(end,end).w,mean(QallNoTool(:,:,:,:,:,:,plBatchs(iPl)),6));
    
    title(['Batch ' num2str(plBatchs(iPl))]);

    set(f.Tool.ax{cAx},'xTick',[])
    set(f.Tool.ax{cAx},'yTick',[])

    colorbar off
    caxis([0 2])

    GridOverImage(fS,gca);


    DrawLineBetweenAxes(f.Tool.ax{1},plBatchs(iPl), mnPks(plBatchs(iPl)), ...
                        f.Tool.ax{cAx}, 0, 1.5);



    cAx = cAx + 1;
end

colormap(whitetocol(100,[0 0 0.7 ],[1]));



% plot(allPks(:,:),'.','LineWidth',0.5,'Color',[.8 .8 .8])


% % % allC = {'b','r','y'};
% % % opts.PlotMean = 0;
% % % 
% % % 
% % % for iM =1:1
% % %     opts.c = allC{iM}
% % %     ShadedPlot(1:51,mnPks(iM,:),mnPks(iM,:) - sdPks(iM,:),mnPks(iM,:) + sdPks(iM,:),opts)
% % % end
% % % 
% % % opts.c = [.7 .7 .7];
% % % ShadedPlot(repmat(1:51,[3 1]),mnPks,mnPks - sdPks,mnPks + sdPks,opts)
% % % 
% % % xlabel('batch number')
% % % ylabel('number of peaks')
% % % 
% % % figure, plot(toolPerfs)

% $$$ Neurons below

% % % % (:,iL,iN,iM,iD) 
% % % cBatch = 1;
% % % figure,plot(squeeze(nanmean(npksN(:,:,:,:,:),[1 3 4])))
% % % legend(1:6)


%% $$$ Create subspace analysis for separate networks

% $$$ Do PCA?
% Rows of X correspond to observations and columns correspond to variables.

switch s.nta.comparison
    case 'Valence'
        s2=DefaultSettings(s);
        s2.plt.ON = 0;
        s2.plt.lmbCol=2:s.wrld.size(2)-1;
        s2.plt.plotThrFl=0;
        [Q,allNeurAct] = CalcNetOutput(s2,w,net); [rGl pGl covGl] = NetAnalysis(s2,w,allNeurAct); % Goal correlation
        % stimrow, stimcol, limbrow, limbcol
        glAct=squeeze(nanmean(allNeurAct(:,:,:,s2.plt.lmbCol,s2.plt.startLayer:s2.plt.stopLayer,:),[3]));
        s2.plt.plotThrFl=1;
        [Q,allNeurAct] = CalcNetOutput(s2,w,net); [rThr pThr covThr] = NetAnalysis(s2,w,allNeurAct); % Threat correlation
        thrAct=squeeze(nanmean(allNeurAct(:,:,:,s2.plt.lmbCol,s2.plt.startLayer:s2.plt.stopLayer,:),[3]));

end


for thrR = 1:s.wrld.size(1)
    for thrC = 1:s.wrld.size(2)


        s2=DefaultSettings(s);
        s2.plt.ON = 0;
        s2.plt.lmbCol=2:s.wrld.size(2)-1;
        s2.plt.plotThrFl=0;

        s2.plt.otherStateVars = [thrR thrC];

        [tmpQ,tmpNeurAct] = CalcNetOutput(s2,w,net); 
        % stimrow, stimcol, limbrow, limbcol
        glByThrAct(:,:,:,:,:,thrR,thrC) = squeeze(nanmean(tmpNeurAct(:,:,:,s2.plt.lmbCol,s2.plt.startLayer:s2.plt.stopLayer,:),3));


    end
end



%% $$$ TEMP BREAK --> just on goal activations

% % % glActFlat = permute(glAct,[4 5 1 2 3]);
% % % glActFlat = glActFlat(:,:,:);
% % % iN = 1;
% % % iL = 1;
% % % pca(squeeze(glActFlat(iN,iL,:)) )

% Flatten neurons and layers
glActFlat = glAct(:,:,:,:); 
% Flatten world configurations
glActFlat = permute(glActFlat,[4 1 2 3]);
glActFlat = glActFlat(:,:);

% Remove NaNs
glActFlat = glActFlat(~isnan(glActFlat(:,1)),:);

[coeff,score,~,~,explained] = pca(glActFlat');
% plot(score(:,1), score(:,2), 'o'); % plotting first two principal components

plot3(score(:,1), score(:,2), score(:,3), 'o'); % plotting first two principal components

plot3(score(:,1), score(:,2), score(:,3), 'o'); % plotting first two principal components


%%

% Create flattened activations for PCA
glByThrActFlat = permute(glByThrAct,[1 2 3 6 7 4 5]);
unwrapSize     = size(glByThrActFlat);
glByThrActFlat = permute(glByThrActFlat(:,:,:,:,:,:),[6 1 2 3 4 5]);
glByThrActFlat = glByThrActFlat(:,:);


% Create 'explanatory variables'
thrRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[1 2 3 7]) ] );  
thrRow = permute(thrRow, [2 3 4 1 5]);
thrCol = repmat( (1:s.wrld.size(2))' , [1, size(glByThrAct,[1 2 3 6]) ] );  
thrCol = permute(thrCol, [2 3 4 5 1]);

golRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[2 3 6 7]) ] );  
golRow = permute(golRow, [1 2 3 4 5]);
golCol = repmat( (1:s.wrld.size(2))' , [1, size(glByThrAct,[1 3 6 7]) ] );  
golCol = permute(golCol, [2 1 3 4 5]);

lmbCol = repmat( (2:s.wrld.size(2)-1)' , [1, size(glByThrAct,[1 2 6 7]) ] );  
lmbCol = permute(lmbCol, [2 3 1 4 5]);

[coeff,score,~,~,explained] = pca(glByThrActFlat( ~isnan(glByThrActFlat(:,1)),:)');

golDist = sqrt((golCol - lmbCol).^2 + (golRow - w.lmb.row).^2);
thrDist = sqrt((thrCol - lmbCol).^2 + (thrRow - w.lmb.row).^2);

minStimDist = min(cat(6,golDist,thrDist),[],6);
maxStimDist = max(cat(6,golDist,thrDist),[],6);
avStimDist  = mean(cat(6,golDist,thrDist),6);

% % % %% Do tSNE
% % % % % % figure,
% % % % % % sneY = tsne(glByThrActFlat( ~isnan(glByThrActFlat(:,1)),:)');
% % % % % % scatter(Y(:,1), Y(:,2));
% % % % % % sneY3 = [sneY, ones([size(sneY,1) 1])];
% % % 
% % % %% Do diffusion mapping
% % % 
% % % tmpGlByThrActFlat = glByThrAct(2:2:end-1,2:2:end,1:2:end,:,:,2:2:end-1,2:2:end);
% % % tmpGlByThrActFlat = permute(tmpGlByThrActFlat,[1 2 3 6 7 4 5]);
% % % unwrapSizeDownSample = size(tmpGlByThrActFlat);
% % % tmpGlByThrActFlat = tmpGlByThrActFlat(:,:,:,:,:,:);
% % % tmpGlByThrActFlat = permute(tmpGlByThrActFlat,[6 1 2 3 4 5]);
% % % tmpGlByThrActFlat = tmpGlByThrActFlat(:,:);
% % % 
% % % neuralData = tmpGlByThrActFlat( ~isnan(tmpGlByThrActFlat(:,1)),:)';
% % % % neuralData = neuralData(1:nStates,:);
% % % 
% % % tic
% % % 
% % % % Compute distance matrix
% % % distanceMatrix = pdist2(neuralData, neuralData);
% % % % Create Kernel matrix
% % % sigma = 1.0; % Scale parameter for Gaussian kernel
% % % kernelMatrix = exp(-distanceMatrix.^2 / (2 * sigma^2));
% % % % Create normalized diffusion matrix
% % % rowSums = sum(kernelMatrix, 2);
% % % diffusionMatrix = diag(1 ./ sqrt(rowSums)) * kernelMatrix * diag(1 ./ sqrt(rowSums));
% % % % Compute eigenvectors and eigenvalues
% % % [eigVectors, eigValues] = eig(diffusionMatrix);
% % % % Select eigenvectors for embedding
% % % k = 3; % Number of dimensions for the embedding
% % % embedding = eigVectors(:, 2:(k+1));
% % % % Plot embedded data
% % % scatter3(embedding(:,1), embedding(:,2), embedding(:,3)); % For a 3D plot
% % % 
% % % 
% % % toc



%% Plot the pcs

figure, 

% 
% % % Downsample
% % tmpScore = reshape(score', [size(score,2) unwrapSize(1:5)] ); % score'
% % tmpScore = tmpScore(:,2:2:end-1,2:2:end,1:2:end,2:2:end-1,2:2:end);
% tmpScore = reshape(embedding', [size(embedding,2) unwrapSizeDownSample(1:5)] ); 
% tmpCl    = thrDist(2:2:end-1,2:2:end,1:2:end,2:2:end-1,2:2:end); 
% tmpSz    = thrDist(2:2:end-1,2:2:end,1:2:end,2:2:end-1,2:2:end);

% Downsample only a little
% tmpScore = reshape(sneY3', [size(sneY3,2) unwrapSize(1:5)] ); 
tmpScore = reshape(score', [size(score,2) unwrapSize(1:5)] ); 
tmpScore = tmpScore(:,:,2:2:end,1:2:end,:,2:2:end);
tmpCl    = thrRow(:,2:2:end,1:2:end,:,2:2:end); 
tmpSz    = golRow(:,2:2:end,1:2:end,:,2:2:end);



tmpScore = tmpScore(:,:)';
tmpCl    = tmpCl(:);
tmpSz    = tmpSz(:);

tmpSz = 10 .* (tmpSz - min(tmpSz)) ./ max(tmpSz) + 1;

value_min = min(tmpCl);
value_max = max(tmpCl);

% Step 2: Choose a colormap
colormap_name = 'jet'; % You can use any built-in colormap or create a custom one

% Step 3: Apply the colormap to your data
nrmCl = (tmpCl - value_min) / (value_max - value_min); % Normalize data to [0, 1]
% color_map = colormap(whitetocol(100,[0 0 0.7])); % Apply the colormap to the normalized data
color_map = colormap(redbluecmapRory); % Apply the colormap to the normalized data
colors = interp1(linspace(0, 1, size(color_map, 1)), color_map, nrmCl); % Interpolate colors



scatter3(tmpScore(:,1), tmpScore(:,2), tmpScore(:,3), tmpSz , colors , 'filled'); 
% scatter3(tmpScore(:,4), tmpScore(:,5), tmpScore(:,6), tmpSz , colors , 'filled'); 


%% Compute lines of best fit for each variable of interest

dS =  1; % Start dimension 
dE =  3; % End dimension [max 84]

[b1 b1Int] = regress(tmpCl, [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);
[b2 b2Int] = regress(tmpSz, [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);

figure,histogram2(tmpCl,(b1(end) + sum(b1(1:end-1).*tmpScore(:,dS:dE)')'),'FaceColor','Flat')
figure,histogram2(tmpSz,(b2(end) + sum(b2(1:end-1).*tmpScore(:,dS:dE)')'),'FaceColor','Flat')

% figure,plot(tmpCl,(b1(end) + sum(b1(1:end-1).*tmpScore(:,dS:dE)')'),'.')
% figure,plot(tmpSz,(b2(end) + sum(b2(1:end-1).*tmpScore(:,dS:dE)')'),'.')


dotProd     = dot(b1(1:end-1), b2(1:end-1)) ;
dotProdAbs  = abs(dotProd) ;
cosTheta    = dotProd ./ (norm(b1(1:end-1)) * norm(b2(1:end-1)));
cosThetaAbs = abs(cosThetaAbs);
theta       = acos(cosTheta);


%% Compute similarity of encoding direction as a function of dimenstion depth

nPCs = 3;

assessDims = 1:size(tmpScore,2) - (nPCs-1);

for iD = assessDims

    dS =  iD; % Start dimension
    dE =  iD + nPCs -1; % End dimension [max 84]

    [b1 b1Int] = regress(tmpCl, [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);
    [b2 b2Int] = regress(tmpSz, [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);

    dotProd(iD)     = dot(b1(1:end-1), b2(1:end-1));
    dotProdAbs(iD)  = abs(dotProd(iD));
    cosTheta(iD)    = dotProd(iD) ./ (norm(b1(1:end-1)) * norm(b2(1:end-1)));
    cosThetaAbs(iD) = abs(cosTheta(iD));
    theta(iD)       = acos(cosTheta(iD));
    thetaAbs(iD)    = abs(theta(iD) - pi./2);

end

% figure,
% plot(assessDims , similarityMetric,'x','LineWidth',2); hold on
% plot(assessDims , dotProdAbs,'-o','LineWidth',2);
% plot(assessDims , cosTheta,'x','LineWidth',2);
% plot(assessDims , cosThetaAbs,'-x','LineWidth',2);
% plot(assessDims , theta,'-x','LineWidth',2);
% plot(assessDims , thetaAbs,'-o','LineWidth',2);
% legend('DotProd', 'DotProdAbs', 'CosTheta', 'CosThetaAbs','Theta','ThetaAbs')
% xlabel('dimension of PCA')
% ylabel('similarity')


figure('Position',[50 50 600 900]),
subplot(3,1,1)
plot(assessDims , theta,'-x','LineWidth',2); hold on
plot(assessDims , ones(size(assessDims)) .* pi./2,'-.k');
title('theta')
xlabel('dimension of PCA')
ylabel('angle between vectors (pi/2 is orthogonal)')


subplot(3,1,2)
plot(assessDims , cosThetaAbs,'-x','LineWidth',2); hold on
plot(assessDims , zeros(size(assessDims)) ,'-.k');
title('costhetaAbs')
xlabel('dimension of PCA')
ylabel('overlap of vectors (0 is orthogonal)')

subplot(3,1,3)
plot(assessDims , thetaAbs,'-x','LineWidth',2); hold on
plot(assessDims , ones(size(assessDims)) .* pi./2,'-.k');
title('theta abs diff from pi/2')
xlabel('dimension of PCA')
ylabel('abs angle diff between vectors (0 is orthogonal)')


%% $$$ NEXT change this to loop over all models, and extract similar lines 
%      for all --> see if there is consistent, model-independent pattern


neurTypes   = {'tansig','logsig','softmax','poslin','purelin','tribas','radbas'};
% % % neurTypes   = {'softmax','logsig','tribas','radbas','poslin','purelin','tansig'};
regTypes    = {'Valence_51','L1'};
regNames    = {'No','L1'};

sAll.plt.ON = 0;



for iTyp = 7:length(neurTypes)
for iReg = 1:length(regTypes)


bFold = 'Results\ForFigures\Valence\';
cFiles = dir([bFold '*' regTypes{iReg} '*_' neurTypes{iTyp} '_B_V*.mat'])

clear glByThrAct
for iRun = 1:length(cFiles)

    load([bFold cFiles(iRun).name]);

    rSall(1:numel(rS),iRun,iTyp,iReg) = rS(:);

    for iM = 1:size(rSall,1)

        s   = rS(iM).s;
        net = rS(iM).net;
        w   = rS(iM).w;

        % Extract data for PCA
        for thrR = 1:s.wrld.size(1)
            for thrC = 1:s.wrld.size(2)
                s2=DefaultSettings(s);
                s2.plt.ON = 0;
                s2.plt.lmbCol=2:s.wrld.size(2)-1;
                s2.plt.plotThrFl=0;

                s2.plt.otherStateVars = [thrR thrC];

                [tmpQ,tmpNeurAct] = CalcNetOutput(s2,w,net);

                % Initialise activations
                if iM == 1 && thrR == 1 && thrC == 1
                    glByThrAct = nan( size(tmpNeurAct,[1 2 2 5 6 3 4]) - [0 0 2 0 0 0 0] );
                end

                % stimrow, stimcol, limbrow, limbcol
                glByThrAct(:,:,:,:,1:size(tmpNeurAct,6),thrR,thrC) = ...
                    squeeze(nanmean(tmpNeurAct(:,:,:,s2.plt.lmbCol,s2.plt.startLayer:s2.plt.stopLayer,:),3));
            end
        end

         save(['C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\DimensionReduction\NetworkActivations\' ...
          'FullNetActivity_NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(iRun) '_NetArch_' num2str(iM) '.mat'], ...
          'glByThrAct','-v7.3')
  

    end

disp('done:')
iTyp
iReg
iRun

end

save('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\DimensionReduction\FullWorkSpace.mat','-v7.3')

end
end

%% $$$ Plot out the PCA results for all the models and neuron types


clear dotProd dotProdAbs cosTheta cosThetaAbs theta thetaAbs crossProd crossProdSize crossProdRelMax

clear b1 b1ScaleDimSd b1ScaleDimSdExpVar b2 b2ScaleDimSd b2ScaleDimSdExpVar

for iTyp = 1:length(neurTypes);
for iReg = 1:length(regTypes);


bFold = 'Results\ForFigures\Valence\';
cFiles = dir([bFold '*' regTypes{iReg} '*_' neurTypes{iTyp} '_B_V*.mat']);


for iRun = 1:length(cFiles)

for iM = 1:size(rSall,1)

load(['C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\DimensionReduction\NetworkActivations\' ...
          'FullNetActivity_NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(iRun) '_NetArch_' num2str(iM) ' .mat']);

% $$$ HERE figure out how to deal with the bloody model architectures --> I
% think just store them separately?

  % Create flattened activations for PCA
    glByThrActFlat = permute(glByThrAct,[1 2 3 6 7 4 5]);
    unwrapSize     = size(glByThrActFlat);
    glByThrActFlat = permute(glByThrActFlat(:,:,:,:,:,:),[6 1 2 3 4 5]);
    glByThrActFlat = glByThrActFlat(:,:);
    
    % Create 'explanatory variables'
    thrRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[1 2 3 7]) ] );  
    thrRow = permute(thrRow, [2 3 4 1 5]);
    thrCol = repmat( (1:s.wrld.size(2))' , [1, size(glByThrAct,[1 2 3 6]) ] );  
    thrCol = permute(thrCol, [2 3 4 5 1]);
    
    golRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[2 3 6 7]) ] );  
    golRow = permute(golRow, [1 2 3 4 5]);
    golCol = repmat( (1:s.wrld.size(2))' , [1, size(glByThrAct,[1 3 6 7]) ] );  
    golCol = permute(golCol, [2 1 3 4 5]);
    
    lmbCol = repmat( (2:s.wrld.size(2)-1)' , [1, size(glByThrAct,[1 2 6 7]) ] );  
    lmbCol = permute(lmbCol, [2 3 1 4 5]);

    golDist = sqrt((golCol - lmbCol).^2 + (golRow - w.lmb.row).^2);
    thrDist = sqrt((thrCol - lmbCol).^2 + (thrRow - w.lmb.row).^2);
    
    minStimDist = min(cat(6,golDist,thrDist),[],6);
    maxStimDist = max(cat(6,golDist,thrDist),[],6);
    avStimDist  = mean(cat(6,golDist,thrDist),6);
    
    % Perform PCA
    [coeff,score,~,~,explained] = pca(glByThrActFlat( ~isnan(glByThrActFlat(:,1)),:)');

    

    % Compute similarity of encoding direction as a function of dimenstion depth
    allVars     = {'thrRow','thrCol','lmbCol','lmbCol','golDist',...
               'minStimDist','minStimDist'};
    allVarsComp = {'golRow','golCol','thrCol','golCol','thrDist',...
               'maxStimDist','avStimDist'};
    
    % Downsample only a little
    tmpScore = reshape(score', [size(score,2) unwrapSize(1:5)] );
    tmpScore = tmpScore(:,:,2:2:end,1:2:end,:,2:2:end);
    tmpScore = tmpScore(:,:)';

    % Set up variables for computing similarity of encoding direction
    nPCs = 3;
    assessDims = 1:size(tmpScore,2) - (nPCs-1);
    for iV = 1:numel(allVars)

        eval(['tmpCl    = ' allVars{iV} '(:,2:2:end,1:2:end,:,2:2:end);']);
        eval(['tmpSz    = ' allVarsComp{iV} '(:,2:2:end,1:2:end,:,2:2:end);']);

        tmpCl    = tmpCl(:);
        tmpSz    = tmpSz(:);
        for iD = assessDims
            dS =  iD; % Start dimension
            dE =  iD + nPCs - 1; % End dimension [max 84]
            
            % Find vector that best explains the variable of interest A
            [b1(:,iD,iV,iTyp,iReg,iRun,iM), b1Int, ~, ~, stats] = regress(tmpCl, [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);
            % Rescale vector by the scale of the dimensions
            b1ScaleDimSd(:,iD,iV,iTyp,iReg,iRun,iM) = b1(:,iD,iV,iTyp,iReg,iRun,iM) .* [std(tmpScore(:,dS:dE)) 1]';
            % Rescale vector by the explained variance
            b1ScaleDimSdExpVar(:,iD,iV,iTyp,iReg,iRun,iM) = b1ScaleDimSd(:,iD,iV,iTyp,iReg,iRun,iM) .* sqrt(stats(1));
            % Find vector that best explains the variable of interest B
            [b2(:,iD,iV,iTyp,iReg,iRun,iM), b2Int, ~, ~, stats] = regress(tmpSz, [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);
            % Rescale vector by the scale of the dimensions
            b2ScaleDimSd(:,iD,iV,iTyp,iReg,iRun,iM) = b2(:,iD,iV,iTyp,iReg,iRun,iM) .* [std(tmpScore(:,dS:dE)) 1]';
            % Rescale vector by the explained variance
            b2(:,iD,iV,iTyp,iReg,iRun,iM) = b2ScaleDimSdExpVar(:,iD,iV,iTyp,iReg,iRun,iM) .* sqrt(stats(1));

% % %             dotProd(iD,iV,iTyp,iReg,iRun,iM)     = dot(b1(1:end-1), b2(1:end-1));
% % %             dotProdAbs(iD,iV,iTyp,iReg,iRun,iM)  = abs(dotProd(iD,iV,iTyp,iReg,iRun,iM));
% % %             cosTheta(iD,iV,iTyp,iReg,iRun,iM)    = dotProd(iD,iV,iTyp,iReg,iRun,iM) ./ (norm(b1(1:end-1)) * norm(b2(1:end-1)));
% % %             cosThetaAbs(iD,iV,iTyp,iReg,iRun,iM) = abs(cosTheta(iD,iV,iTyp,iReg,iRun,iM));
% % %             theta(iD,iV,iTyp,iReg,iRun,iM)       = acos(cosTheta(iD,iV,iTyp,iReg,iRun,iM));
% % %             thetaAbs(iD,iV,iTyp,iReg,iRun,iM)    = abs(theta(iD,iV,iTyp,iReg,iRun,iM) - pi./2);
% % % 
% % %             crossProd(:,iD,iV,iTyp,iReg,iRun,iM)     = cross(b1(1:end-1), b2(1:end-1));
% % %             crossProdSize(iD,iV,iTyp,iReg,iRun,iM)   = norm(crossProd(:,iD,iV,iTyp,iReg,iRun,iM));
% % %             crossProdRelMax(iD,iV,iTyp,iReg,iRun,iM) = crossProdSize(iD,iV,iTyp,iReg,iRun,iM) ./ ...
% % %                                                     (norm(b1(1:end-1)) .* norm(b2(1:end-1)));
        end


        % Plot the PCs
        if sAll.plt.ON == 1 && iV == 1
            tmpSz = 10 .* (tmpSz - min(tmpSz)) ./ max(tmpSz) + 1;

            value_min = min(tmpCl);
            value_max = max(tmpCl);

            % Step 2: Choose a colormap
            colormap_name = 'jet'; % You can use any built-in colormap or create a custom one

            % Step 3: Apply the colormap to your data
            nrmCl = (tmpCl - value_min) / (value_max - value_min); % Normalize data to [0, 1]
            color_map = colormap(redbluecmapRory); % Apply the colormap to the normalized data
            colors = interp1(linspace(0, 1, size(color_map, 1)), color_map, nrmCl); % Interpolate colors

            figure
            scatter3(tmpScore(:,1), tmpScore(:,2), tmpScore(:,3), tmpSz , colors , 'filled');
            disp('testing til here')

            % Save figure
            saveas(gcf,['C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\' ...
                'Results\ForFigures\DimensionReduction\PCAs\NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(iRun) '_NetArch_' num2str(iM) ] , 'fig')

            close all
        end

    end


    

disp('done:')
iTyp
iReg
iRun

end
end
end
end

save(['C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\DimensionReduction\Similarities\' ...
          'SimilarityScores.mat'], ...
          'b1','b1ScaleDimSd','b1ScaleDimSdExpVar','b2','b2ScaleDimSd','b2ScaleDimSdExpVar','-v7.3');


    %% $$$ PLOT

    % $$$ Here I'm re-storing the logsig run 1. Then I should rerun the PCA business, including all the different network architectures

    % $$$ HERE - I just tried rerunning everything :D

    % $$$ should als redo for each model architecture
% $$$ And should figure out why b1 and b2 are still getting so huge --> maybe do median split thing
% $$$ AND I should store b1 and b2, AS WELL AS OTHER VECTORS
% $$$ Then, just below here, I can recompute the dot-products etc etc
    
    iV   = 1; % comparison Variables
    iTyp = 7; % Neuron type
    iReg = 2; % Regression type
    

figure('Position',[2000 50 600 900]),
subplot(4,1,1)
plot(assessDims , squeeze(theta(:,iV,iTyp,iReg,:))); hold on
plot(assessDims , squeeze(mean(theta(:,iV,iTyp,iReg,:),5)) ,'-k','LineWidth',2); 
plot(assessDims , ones(size(assessDims)) .* pi./2,'-.k');
title('theta')
xlabel('dimension of PCA')
ylabel('angle between vectors (pi/2 is orthogonal)')


subplot(4,1,2)
plot(assessDims , squeeze(cosThetaAbs(:,iV,iTyp,iReg,:))); hold on
plot(assessDims , squeeze(mean(cosThetaAbs(:,iV,iTyp,iReg,:),5)),'-xk','LineWidth',2); 
plot(assessDims , zeros(size(assessDims)) ,'-.k');
title('costhetaAbs')
xlabel('dimension of PCA')
ylabel('overlap of vectors (0 is orthogonal)')

subplot(4,1,3)
plot(assessDims , squeeze(thetaAbs(:,iV,iTyp,iReg,:))); hold on
plot(assessDims , squeeze(mean(thetaAbs(:,iV,iTyp,iReg,:),5)),'-xk','LineWidth',2); 
plot(assessDims , ones(size(assessDims)) .* pi./2,'-.k');
title('theta abs diff from pi/2')
xlabel('dimension of PCA')
ylabel('abs angle diff between vectors (0 is orthogonal)')

subplot(4,1,4)
plot(assessDims , squeeze(dotProd(:,iV,iTyp,iReg,:))); hold on
plot(assessDims , squeeze(mean(dotProd(:,iV,iTyp,iReg,:),5)),'-xk','LineWidth',2); 
plot(assessDims , zeros(size(assessDims)) ,'-.k');
title('dot product')
xlabel('dimension of PCA')
ylabel('dot product between vectors (0 is orthogonal)')
ylim([-1 1])


figure('Position',[2000 50 600 900]),
subplot(4,1,1)
plot(assessDims , squeeze(crossProdRelMax(:,iV,iTyp,iReg,:))); hold on
plot(assessDims , squeeze(mean(crossProdRelMax(:,iV,iTyp,iReg,:),5)) ,'-k','LineWidth',2); 
plot(assessDims , ones(size(assessDims)) .* 1,'-.k');
title('xprod relative to maximum')
xlabel('dimension of PCA')
ylabel('xProd')

subplot(4,1,2)
plot(assessDims , squeeze(crossProdSize(:,iV,iTyp,iReg,:))); hold on
plot(assessDims , squeeze(mean(crossProdSize(:,iV,iTyp,iReg,:),5)) ,'-k','LineWidth',2); 
plot(assessDims , ones(size(assessDims)) .* 0 ,'-.k');
title('xprod size')
xlabel('dimension of PCA')
ylabel('xProd')


% $$$ NOW scaling the vectors by the variance to get more indicative
% vectors of the effect size

% $$$ For some reason it's still huge --> is B1 wrong or shold I do center
% of mass of top and bottom half?


%% $$$ QUESTION: DO I need to do it on the goal x threat data? 
%      --> I guess so, because i want to find separate subspaces for each,
%      right? Like have a 'goal close' and a 'threat close' type situation.
%      So then I need to
%      A) Figure out how to 'reconstruct' which element of each of the
%      score rows corresponds to exactly which element of the activity
%      matrix, and
%      B) Creat allNeurAct in which there is a dimension for goals AND for
%      threats

%% Functions



function [ax, plQhoriz, plQvert] = plot2D(rS,fS,Q,iM,iD);

    rS(iM,iD).s.plt.lmbCol = 2:14;
    rS(iM,iD).s.plt.meanLimbCols = 1;
    rS(iM,iD).s.plt.colLims = [.5 14.5];
    rS(iM,iD).s.plt.rowLims = [.5 14.5];
    DisplActValsFun(rS(iM,iD).s,rS(iM,iD).w,Q);
    colorbar off
    GridOverImage(fS,gca);
    xlim([-4.5 4.5]); ylim([1.5 13.5]);
    axis off

    ax{1} = gca;

    % Plot average OR line through Q-values

    ax{2} = axes('Position',ax{1}.Position .* [1 1 1 0] + [0 -0.05 0 0.05]  );
    plQhoriz = squeeze(nanmean(Q(1,8,2:11,2:end-1,:),[3 5])); % lmb R, C, stim R, C, A
% % %     plQ = squeeze(nanmean(Q(1,8,9,2:end-1,:),[3 5])); % lmb R, C, stim R, C, A
    plot((1:size(Q,4))' - 0.5 , [0; plQhoriz; 0 ],'LineWidth',2); hold on
% % %     plot(fS.xPoints,[-2 2],'-.k');
    xLims = xlim;
% % %     plot(xLims,[0.25 0.25 ; -0.25 -0.25]','-.k');
    xlim(xLims);
    hold on
% % %     plot([7.5 7.5],[-2 2],'-.k');
    ylim(fS.lineYLim);
    grid on
%     ax{2}.Visible = 'off';
    ax{2}.XAxis.Visible = 'off';
    ax{2}.YAxis.Visible = 'off';
    


    ax{3} = axes('Position',ax{1}.Position .* [1 1 0 1] + [-0.03 0 0.03 0]  );
    plQvert = squeeze(nanmean(Q(1,8,2:end-1,2:end-1,:),[4 5])); % lmb R, C, stim R, C, A
% % %     plQ = squeeze(nanmean(Q(1,8,2:end-1,8,:),[4 5])); % lmb R, C, stim R, C, A
    plot( plQvert, (numel(plQvert):-1:1)' - 0.5 , 'LineWidth',2); hold on
    xlim(fS.lineXLim);
% % %     plot([-2 2],[4 4],'-.k');
    grid on
%     ax{3}.Visible = 'off';
    ax{3}.XAxis.Visible = 'off';
    ax{3}.YAxis.Visible = 'off';
    

end


function [] = DrawLineBetweenAxes(ax_main,x_main,y_main, ax_sub,x_sub,y_sub)

% Assuming ax_main and ax_sub are your main and subplot axes, and pt_main and pt_sub are the point coordinates in the respective axes
pt_main = [x_main, y_main]; % Replace with actual main axes point
pt_sub = [x_sub, y_sub]; % Replace with actual subplot point

% % % % Convert axes coordinates to normalized figure coordinates
% % % pt_main_fig = ax_main.Position(1:2) + pt_main .* ax_main.Position(3:4);
% % % pt_sub_fig = ax_sub.Position(1:2) + pt_sub .* ax_sub.Position(3:4);


% Normalise units
pt_main_normalized = [(pt_main(1) - ax_main.XLim(1)) / diff(ax_main.XLim), (pt_main(2) - ax_main.YLim(1)) / diff(ax_main.YLim)];
pt_sub_normalized = [(pt_sub(1) - ax_sub.XLim(1)) / diff(ax_sub.XLim), (pt_sub(2) - ax_sub.YLim(1)) / diff(ax_sub.YLim)];


pt_main_fig = ax_main.Position(1:2) + pt_main_normalized .* ax_main.Position(3:4);
pt_sub_fig = ax_sub.Position(1:2) + pt_sub_normalized .* ax_sub.Position(3:4);


% Draw the annotation line
annotation('line', [pt_main_fig(1), pt_sub_fig(1)], [pt_main_fig(2), pt_sub_fig(2)] , 'LineStyle', '--' );

end


