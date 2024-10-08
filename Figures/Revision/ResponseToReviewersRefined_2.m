%% Effects of policy

addpath(genpath('Scripts\EgocentricValueMaps'))

load('Figures\SARSA_WithUnbalanced_AndThreat.mat')
% load('Results\ForFigures\SARSA_WithUnbalanced_AndThreat.mat')
% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\SARSA_WithUnbalanced_AndThreat.mat')
% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\SARSA_WithUnbalanced.mat')

f.Policy.f = figure('Position',[20 20 1800 900]);

sgtitle('POSITIVE VALENCE STIMULI                                                              NEGATIVE VALENCE STIMULI')


fS.lineYLim = [-0.01 0.55];
fS.lineXLim = [ 0   0.7];
fS.xPoints  = [5 5 ; 10 10]';
% fS.cAxes    = [-1.8 1.8];
fS.cAxes    = [-2 2];

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
    Q = nanmean(tmpQ,[5 6]);
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
    Q = nanmean(tmpQ,[5 6]);
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
    Q = nanmean(tmpQ,[5 6]);
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
    Q = nanmean(tmpQ,[5 6]);
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
    Q = nanmean(tmpQ,[5 6]);
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
    Q = nanmean(tmpQ,[5 6]);
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


f.Colourbars.f = figure('Position',[20 50 400 600]);
colorbar;
caxis(fS.cAxes);
colormap(flip(redbluecmapRory(11,5)))

% % Savefigure
% allFields = fields(f);
% for iF = 1:length(allFields)
%     cF = allFields{iF};
%   
%     set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
%     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Revision\LearningAlgorithm\' cF '.eps'] , 'epsc')
%     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Revision\LearningAlgorithm\' cF '.pdf'] , 'pdf')
% 
% end


%% Wasp scenario

% $$$ MAKE THIS into a loop and plot the best action for all 3
% 'environments' (columns) and hand-body distances (rows)

load('Data\The_BEES_plus_AWAYoption_v2.mat')
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


%% For wasp, also plot value of each action individually


% $$$ Maybe also do gridoverimage?

% cAx = [-1 1];


% % %         xlim([20  41]);
% % %         ylim([40  61]);

% % %         xlim([23  38]);
% % %         ylim([43  58]);

xLims = [12 49];
yLims = [22 59];

% % % xLims = [20 41];
% % % yLims = [41 60];


for iAct = 1:6
figN = ['TwoDimPlotQAct_TiesPlottedCorrectly' num2str(iAct)];
% % % figN = ['TwoDimPlotQmean'];

% f.(figN).f        = figure('Position',[20 20 1200 1200]);
f.(figN).f        = figure('Position',[20 20 750 1200]);

iPl = 0;

for iHP = 1:nHandPos

    cHP = allHP(iHP);
    for iDynam = 1:nDynam

        if ~ (iDynam == 1 & iAct == 6)

        iPl = iPl + 1;

        newQ = allQ.qVals{iDynam};
        sTmp = sBdy;
        

% % %         newQ2 = newQ(iAct,:,:,:); % Max value

% % %         % Mean action value
% % %         newQ2 = mean(-newQ,1); % Max value


% % %         % Differential action value compared to the mean or max of the others
% % %         blAct = 1:size(newQ,1);
% % %         blAct(blAct == iAct) = [];
% % % %         newQ2 = newQ(iAct,:,:,:) - mean(newQ(blAct,:,:,:),1);
% % %         newQ2 = newQ(iAct,:,:,:) - max(newQ(blAct,:,:,:),[],1);


        % ONLY POSITIVE part of the Differential action value compared to 
        % the mean or max of the others
        blAct                   = 1:size(newQ,1);
        blAct(blAct == iAct)    = [];
        newQtmp = newQ;
% % %         % For moving, give all actions some tiny penalty
% % %         newQtmp(2:end,:,:,:) = newQtmp(2:end,:,:,:) - 0.001;
        newQ2                   = newQtmp(iAct,:,:,:) - max(newQtmp(blAct,:,:,:),[],1);        
% % %         newQ2(newQ2 == 0)       = 0.00000001; % If there is a tie, set it to a tiny value

% % %         newQ2(newQ2 < 0)        = 0;
% % %         % Set tot NaN for graying out later
% % %         newQ2(newQ2 == 0 )       = NaN;

        
        % Set tot NaN for graying out later
        newQ2(newQ2 < 0)                = NaN;
        newQ2(newQ2 == 0 & newQtmp(1,:,:,:)==0)  = NaN;
        


        f.(figN).ax{iPl}    = subplot(nHandPos,nDynam,iPl);

        zPos = sTmp.clc.nearPos(3) - cHP; % POSITION OF SLICE
% % %         if iAct == 1
            % Set an arbitrary point very close to 0, so that the colourmap
            % still works
            newQ2(1,1,1,zPos) = 0.00000001;
% % %         end

        imagesc(squeeze(newQ2(1,:,:,zPos )) ); hold on
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

%         caxis auto
%         colorbar
        cAXES = caxis;
%         caxis(max(abs(cAXES)) .* [-1 1]);
%         caxis([0.0001 max(cAXES)]);

% % %         title(['Act' num2str(iAct) '.col' num2str(cAXES)])

        xlim(xLims);
        ylim(yLims);

        % PLOT THE mask for when the action is not taken
        maskDat = isnan(squeeze(newQ2(1,:,:,zPos )));
        ax2 = axes('Position', f.(figN).ax{iPl}.Position);
        h = imagesc(maskDat);
        % Adjust the AlphaData property of the image object
        set(h, 'AlphaData', 0.2 .* maskDat);
        colormap(ax2,whitetocol(2,[0 0 0]));

        xlim(xLims);
        ylim(yLims);

        axis off;
        linkaxes([f.(figN).ax{iPl}, ax2]); % link the two axes for zooming, panning, etc.

        % Hide the top axes box but retain the ability to interact with it
        set(ax2, 'Color', 'none', 'Box', 'off');

        colormap(f.(figN).ax{iPl},whitetocol(256,[0 0.7 0]));
% % %         colormap(f.(figN).ax{iPl},whitetocol(256,[0 0 0]));
        axis off


        end

    end

% % %     sgtitle(['ACTION' num2str(iAct)]);
end

% colormap(redbluecmapRory)

% % % % Add grey for Nan
% % % if any(isnan(newQ2(:)))
% % %     colormap([ [.95 .95 .95] ;  whitetocol(256000,[0 0.7 0]) ]);
% % % else
% % %     colormap(whitetocol(256,[0 0.7 0]));
% % % end

f.(figN).ax{iPl+1} = axes('Position', [0.9 0.05 0.1 0.9]);
colorbar;
% caxis(cAx);
f.(figN).ax{iPl+1}.Visible = 'off';

end


%% Save figures for wasp stuff
allFields = fields(f);
for iF = 1:length(allFields)5
    
    cF = allFields{iF};
  
    set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['Documentation\Figures\Revision\Wasp\' cF '.eps'] , 'epsc')
    saveas(f.(cF).f,['Documentation\Figures\Revision\Wasp' cF '.pdf'] , 'pdf')
% % %     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Revision\NeurTypes\CartDist' cF '.eps'] , 'epsc')
% % %     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Revision\NeurTypes\CartDist' cF '.pdf'] , 'pdf')

end



%% $$$ Plot some kind of order metrics or something --> $$$ !!! POSSIBLY REMOVE

% % % load('Results\ForFigures\NeurTypesAnalysed_V5c.mat')

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

% $$$ HERE CAN ALSO PLOT AGAINST PERFORMANCE

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

% % % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\EatingModel.mat')
% % % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\EatingModel_WithRandomMove.mat')


iM = 1;

f.Eating.f = figure('Position',[50 200 1600 400]);

eatRS = ntRS;

eatRS(iM).s.plt.meanLimbCols  = 1
eatRS(iM).s.plt.holdingRew    = 0;
eatRS(iM).s.plt.lmbCol        = 10;
eatRS(iM).s.plt.plAct         = 2;

fS.lineYLim = [-0.2 0];
fS.lineXLim = [-0.2 0];
fS.xPoints  = [5 5 ; 10 10]';
% fS.cAxes    = [-1.8 1.8];
% fS.cAxes    = [-2 2];
fS.cAxes    = [-.4 .4];


eatRS(iM).s.plt.colLims = [1.5 13.5];
eatRS(iM).s.plt.rowLims = [1.5 13.5];
fS.gridXstart = -5.5;
fS.gridXstep = 1;
fS.gridYstart = 2.5;
fS.gridYstep = 1;



subplot(2,6,1)
eatRS(iM).s.plt.lmbCol  = 2:14;
[Q,allNeurAct] = CalcNetOutput(eatRS(iM).s,eatRS(iM).w,eatRS(iM).net);
DisplActValsFun(eatRS(iM).s,eatRS(iM).w,-Q);
colorbar off
GridOverImage(fS,gca)
caxis(fS.cAxes);
title('Positive stimulus, average Q-field')

MakeSidePlots(gca,Q,fS,eatRS(iM))
mouthOffset = [-4:2:4]; % [-2:2] ; 
for iO = 1:length(mouthOffset) % loop through offsets

    % current offset
    cLimbCol = 8 - mouthOffset(iO);
    subplot(2,6,1 + iO);
    eatRS(iM).s.plt.lmbCol  = cLimbCol;
    

    
    [tmpQ,allNeurAct] = CalcNetOutput(eatRS(iM).s,eatRS(iM).w,eatRS(iM).net);
    Q = nanmean(tmpQ,[5 6]);

    eatRS(iM).s.plt.plAct = 1;
    DisplActValsFun(eatRS(iM).s,eatRS(iM).w,-Q);


    colorbar off
    GridOverImage(fS,gca)
    caxis(fS.cAxes);
% % %     title(['REWARD: Mouth at ' num2str(mouthOffset(iO)) ])

    MakeSidePlots(gca,Q,fS,eatRS(iM))
end


colormap(redbluecmapRory)

% % % 
% % % % Shield version --> stationary dangerzone
% % % 
% % % 
% % % iM = 1;
% % % 
% % % % % % f.Eating.f = figure('Position',[50 200 1600 400]);y
% % % 
% % % shieldRS(iM).s.plt.meanLimbCols  = 1
% % % shieldRS(iM).s.plt.holdingRew    = 0;
% % % shieldRS(iM).s.plt.lmbCol        = 10;
% % % shieldRS(iM).s.plt.plAct         = 2;
% % % 
% % % shieldRS(iM).s.plt.colLims = [1.5 13.5];
% % % shieldRS(iM).s.plt.rowLims = [1.5 13.5];
% % % fS.gridXstart = -5.5;
% % % fS.gridXstep = 1;
% % % fS.gridYstart = 2.5;
% % % fS.gridYstep = 1;
% % % 
% % % 
% % % subplot(2,9,[8 9 17 18])
% % % shieldRS(iM).s.plt.lmbCol  = 2:14;
% % % [Q,allNeurAct] = CalcNetOutput(shieldRS(iM).s,shieldRS(iM).w,shieldRS(iM).net);
% % % DisplActValsFun(shieldRS(iM).s,shieldRS(iM).w,-Q);
% % % colorbar off
% % % GridOverImage(fS,gca)
% % % caxis([-.1 .1]);
% % % title('Negative stimulus, only harms mouth, average Q-field')
% % % 
% % % mouthOffset = [-4:2:4]; % [-2:2] ; 
% % % for iO = 1:length(mouthOffset) % loop through offsets
% % % 
% % %     % current offset
% % %     cLimbCol = 8 - mouthOffset(iO);
% % %     subplot(2,9,11 + iO)
% % %     shieldRS(iM).s.plt.lmbCol  = cLimbCol;
% % %     [Q,allNeurAct] = CalcNetOutput(shieldRS(iM).s,shieldRS(iM).w,shieldRS(iM).net);
% % %     DisplActValsFun(shieldRS(iM).s,shieldRS(iM).w,-Q);
% % %     colorbar off
% % %     GridOverImage(fS,gca)
% % %     caxis([-.4 .4]);
% % %     title(['DANGER: Mouth at ' num2str(mouthOffset(iO)) ])
% % % end
% % % 
% % % 
% % % 
% % % colormap(redbluecmapRory)



% Save figures for eating
allFields = fields(f);
for iF = 1:length(allFields)
    cF = allFields{iF};
  
    set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Revision\Eating\' cF '.eps'] , 'epsc')
    saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Revision\Eating\' cF '.pdf'] , 'pdf')
% % %     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Revision\NeurTypes\CartDist' cF '.eps'] , 'epsc')
% % %     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Revision\NeurTypes\CartDist' cF '.pdf'] , 'pdf')

end
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
                npksN(:,iL,iN,iM,iD) = RowPeakFind(sFP, permute(allNeurAct(:,:,:,:,iL,iN),[3 4 1 2 5 6]));
            end
        end

        % Also calculate values for no tool
        sFP.plt.ToolPresent = 0;
        [QnoTool,allNeurAct] = CalcNetOutput(sFP,w,net);
        QallNoTool(:,:,:,:,:,iM,iD,1) = QnoTool;

        for iN = 1:size(allNeurAct,6)
            for iL = 1:size(allNeurAct,5)
                npksNnoTool(:,iL,iN,iM,iD) = RowPeakFind(sFP, permute(allNeurAct(:,:,:,:,iL,iN),[3 4 1 2 5 6]));
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

save('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Tool_ReviewerResponse_Pre_post_51_51Batch_ToolPos5_moreRandSpr_NoHist_PEAKSANALYSED_V2.mat','-v7.3')


%% Calculate Incremental tool improvements


% % % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Tool_ReviewerResponse_Pre_post_51_51Batch_ToolPos5_moreRandSpr_NoHist_PEAKSANALYSED.mat')

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



%% Plot tool results


fS.lineYLim = [-0.55 0.01 ];
fS.lineXLim = [-0.8  0];
fS.xPoints  = [5 5 ; 10 10]';
fS.cAxes    = [-1.8 1.8];

fS.gridXstart   = -5.5;
fS.gridXstep    =  1;
fS.gridYstart   =  2.5;
fS.gridYstep    =  1;


% % % fS2 = fS;
% % % fS2.lineXLim = [-0.8 0];
% % % fS2.lineYLim = [-0.8 0];




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
allPks2 = allPks;
allPks  = allPks(:,:);

% Average across model architectures
mnPks2 = mnPks;
mnPks  = nanmean(mnPks);
sdPks = sdPks ./ sqrt(size(allPks2,3)); % convert to standard error
sdPks  = sum(sqrt(sdPks .^2 ));


% subplot(4,3,[1:9]);

yyaxis left

[N C] = hist3( [repmat(  (1:size(allPks,1))', [size(allPks,2) 1]), allPks(:) ],[51 3]);

% imagesc(N' ./ max(N(:))); axis xy; hold on

f.Tool.ax{1} = gca;

plot(mnPks','b-o','LineWidth',2); hold 
plot(mnPks2');
opts.c = [0 0 .7];
opts.PlotMean = 0;
ShadedPlot((1:51),mnPks,mnPks - sdPks,mnPks + sdPks,opts);
plot(mnPks','b-o','LineWidth',2); hold 


% ylim([0 3.5])
ylim([0.8 3])
xlabel('Training batch since exposure to tool')
ylabel('Number of peaks in receptive field')


title('A receptive field grows around the tooltip as the agent learns to use the tool')

xlim([0 51])

box off



% Plot receptive fields
sFP.plt.pltType         = 'Imagesc';
sFP.plt.ON              = 1;
sFP.plt.meanLimbCols    = 1;
sFP.plt.rowLims         = [1.5 13.5];


% Plot the tool-fields as the agent learns
plBatchs = [1 2 3 7 51];
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


    % $$$ Make side plots
    MakeSidePlots(f.Tool.ax{cAx},mean(Qall(:,:,:,:,:,:,plBatchs(iPl)),6),fS,incrToolRS(1));

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


% Savefigure
allFields = fields(f);
for iF = 1:length(allFields)
    cF = allFields{iF};
  
    set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['Results\ForFigures\Revision\ToolUse\' cF '.eps'] , 'epsc')
    saveas(f.(cF).f,['Results\ForFigures\Revision\ToolUse\' cF '.pdf'] , 'pdf')

end



%% Do tool stats



% Calculate number of peaks without a tool --------------------------------
tic
clear toolPks noToolPks
for iM = 1:size(incrToolRS,1)
    for iD = 1:size(incrToolRS,2)
        toolPks(iD,iM,:)     = incrToolRS(iM,iD).npks;
        sFP.plt.ToolPresent = 0;
        [noToolPks(iD,iM,:), dmy1, dmy2]  = ...
            RowPeakFind(sFP,QallNoTool(:,:,:,:,:,iM,iD,1));
    end
end
toc

% Display stats for number of tool fields 
disp('After training, tool vs no tool')
[pTl hTl stats] = signrank(toolPks(end,:),noToolPks(end,:),'method','approximate')
effSizR1 = stats.zval ./ sqrt(numel(toolPks(end,:).*2))


% Effect of training ------------------------------------------------------
batchNum = repmat((1:size(toolPks,1))',[1 size(toolPks,2) size(toolPks,3)]);
tmpX = permute(batchNum,[2 1 3]); tmpX = tmpX(:,:);
tmpY = permute(toolPks,[2 1 3]);  tmpY = tmpY(:,:);
[rhorho pp] = corr(tmpX',tmpY');

pp = pp(1,:);
rhorho = rhorho(1,:);

[dmy1 dmy2 pp] = fdr(pp);

disp('Minimum correlation, p:')
[maxP maxInd] = max(pp(:))
disp('Minimum correlation, rho:')
rhorho(maxInd)

disp('mean rho')
nanmean(rhorho)
disp('std rho')
nanstd(rhorho)

disp('mean peaks with tool:')
mean(toolPks(end,:))
disp('mean sd with tool:')
std(toolPks(end,:))

disp('mean peaks without tool:')
mean(noToolPks(end,:))
disp('mean sd without tool:')
std(noToolPks(end,:))


%% Convert to format that can be used by old neuron analysis


olRS = incrToolRS(:,1);
ytRS = incrToolRS(:,end);

[rMatOL] = DisplayProxStats(olRS);
[rMatYT] = DisplayProxStats(ytRS);

nPksNLikeOld = cat(5,npksNnoTool(:,:,:,:,1),npksN(:,:,:,:,1),npksNnoTool(:,:,:,:,end),npksN(:,:,:,:,end));

rMatOL.npksN      = nPksNLikeOld;
rMatYT.npksN      = nPksNLikeOld; 


rMat.nPeaksPerNeur = squeeze(nanmean(nPksNLikeOld(:,:,:,:,:),[1]));
rMat.peaksPerLay   = squeeze(nanmean(nPksNLikeOld(:,:,:,:,:),[1 2]));




%%

% $$$ USE THESE BITS TO MAKE IT WORK#

% $$$ TO DO THIS, I need to shape the new data back into the old format
% --> It's incrToolRS that I need to change into olRS or whatever it was







npksNnoTool

QallNoTool
%% $$$ Create subspace analysis for separate networks

% % % % Rows of X correspond to observations and columns correspond to variables.
% % % 
% % % switch s.nta.comparison
% % %     case 'Valence'
% % %         s2=DefaultSettings(s);
% % %         s2.plt.ON = 0;
% % %         s2.plt.lmbCol=2:s.wrld.size(2)-1;
% % %         s2.plt.plotThrFl=0;
% % %         [Q,allNeurAct] = CalcNetOutput(s2,w,net); [rGl pGl covGl] = NetAnalysis(s2,w,allNeurAct); % Goal correlation
% % %         % stimrow, stimcol, limbrow, limbcol
% % %         glAct=squeeze(nanmean(allNeurAct(:,:,:,s2.plt.lmbCol,s2.plt.startLayer:s2.plt.stopLayer,:),[3]));
% % %         s2.plt.plotThrFl=1;
% % %         [Q,allNeurAct] = CalcNetOutput(s2,w,net); [rThr pThr covThr] = NetAnalysis(s2,w,allNeurAct); % Threat correlation
% % %         thrAct=squeeze(nanmean(allNeurAct(:,:,:,s2.plt.lmbCol,s2.plt.startLayer:s2.plt.stopLayer,:),[3]));
% % % 
% % % end
% % % 
% % % 
% % % for thrR = 1:s.wrld.size(1)
% % %     for thrC = 1:s.wrld.size(2)
% % % 
% % % 
% % %         s2=DefaultSettings(s);
% % %         s2.plt.ON = 0;
% % %         s2.plt.lmbCol=2:s.wrld.size(2)-1;
% % %         s2.plt.plotThrFl=0;
% % % 
% % %         s2.plt.otherStateVars = [thrR thrC];
% % % 
% % %         [tmpQ,tmpNeurAct] = CalcNetOutput(s2,w,net); 
% % %         % stimrow, stimcol, limbrow, limbcol
% % %         glByThrAct(:,:,:,:,:,thrR,thrC) = squeeze(nanmean(tmpNeurAct(:,:,:,s2.plt.lmbCol,s2.plt.startLayer:s2.plt.stopLayer,:),3));
% % % 
% % % 
% % %     end
% % % end



%% $$$ TEMP BREAK --> just on goal activations

% % % % % % glActFlat = permute(glAct,[4 5 1 2 3]);
% % % % % % glActFlat = glActFlat(:,:,:);
% % % % % % iN = 1;
% % % % % % iL = 1;
% % % % % % pca(squeeze(glActFlat(iN,iL,:)) )
% % % 
% % % % Flatten neurons and layers
% % % glActFlat = glAct(:,:,:,:); 
% % % % Flatten world configurations
% % % glActFlat = permute(glActFlat,[4 1 2 3]);
% % % glActFlat = glActFlat(:,:);
% % % 
% % % % Remove NaNs
% % % glActFlat = glActFlat(~isnan(glActFlat(:,1)),:);
% % % 
% % % [coeff,score,~,~,explained] = pca(glActFlat');
% % % % plot(score(:,1), score(:,2), 'o'); % plotting first two principal components
% % % 
% % % plot3(score(:,1), score(:,2), score(:,3), 'o'); % plotting first two principal components
% % % 
% % % plot3(score(:,1), score(:,2), score(:,3), 'o'); % plotting first two principal components


%%

% % % % Create flattened activations for PCA
% % % glByThrActFlat = permute(glByThrAct,[1 2 3 6 7 4 5]);
% % % unwrapSize     = size(glByThrActFlat);
% % % glByThrActFlat = permute(glByThrActFlat(:,:,:,:,:,:),[6 1 2 3 4 5]);
% % % glByThrActFlat = glByThrActFlat(:,:);
% % % 
% % % 
% % % % Create 'explanatory variables'
% % % thrRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[1 2 3 7]) ] );  
% % % thrRow = permute(thrRow, [2 3 4 1 5]);
% % % thrCol = repmat( (1:s.wrld.size(2))' , [1, size(glByThrAct,[1 2 3 6]) ] );  
% % % thrCol = permute(thrCol, [2 3 4 5 1]);
% % % 
% % % golRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[2 3 6 7]) ] );  
% % % golRow = permute(golRow, [1 2 3 4 5]);
% % % golCol = repmat( (1:s.wrld.size(2))' , [1, size(glByThrAct,[1 3 6 7]) ] );  
% % % golCol = permute(golCol, [2 1 3 4 5]);
% % % 
% % % lmbCol = repmat( (2:s.wrld.size(2)-1)' , [1, size(glByThrAct,[1 2 6 7]) ] );  
% % % lmbCol = permute(lmbCol, [2 3 1 4 5]);
% % % 
% % % [coeff,score,~,~,explained] = pca(glByThrActFlat( ~isnan(glByThrActFlat(:,1)),:)');
% % % 
% % % golDist = sqrt((golCol - lmbCol).^2 + (golRow - w.lmb.row).^2);
% % % thrDist = sqrt((thrCol - lmbCol).^2 + (thrRow - w.lmb.row).^2);
% % % 
% % % minStimDist = min(cat(6,golDist,thrDist),[],6);
% % % maxStimDist = max(cat(6,golDist,thrDist),[],6);
% % % avStimDist  = mean(cat(6,golDist,thrDist),6);

%% Do tSNE

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

% % % figure, 
% % % 
% % % % 
% % % % % % Downsample
% % % % % tmpScore = reshape(score', [size(score,2) unwrapSize(1:5)] ); % score'
% % % % % tmpScore = tmpScore(:,2:2:end-1,2:2:end,1:2:end,2:2:end-1,2:2:end);
% % % % tmpScore = reshape(embedding', [size(embedding,2) unwrapSizeDownSample(1:5)] ); 
% % % % tmpCl    = thrDist(2:2:end-1,2:2:end,1:2:end,2:2:end-1,2:2:end); 
% % % % tmpSz    = thrDist(2:2:end-1,2:2:end,1:2:end,2:2:end-1,2:2:end);
% % % 
% % % % Downsample only a little
% % % % tmpScore = reshape(sneY3', [size(sneY3,2) unwrapSize(1:5)] ); 
% % % tmpScore = reshape(score', [size(score,2) unwrapSize(1:5)] ); 
% % % tmpScore = tmpScore(:,:,2:2:end,1:2:end,:,2:2:end);
% % % tmpCl    = thrRow(:,2:2:end,1:2:end,:,2:2:end); 
% % % tmpSz    = golRow(:,2:2:end,1:2:end,:,2:2:end);
% % % 
% % % 
% % % 
% % % tmpScore = tmpScore(:,:)';
% % % tmpCl    = tmpCl(:);
% % % tmpSz    = tmpSz(:);
% % % 
% % % tmpSz = 10 .* (tmpSz - min(tmpSz)) ./ max(tmpSz) + 1;
% % % 
% % % value_min = min(tmpCl);
% % % value_max = max(tmpCl);
% % % 
% % % % Step 2: Choose a colormap
% % % colormap_name = 'jet'; % You can use any built-in colormap or create a custom one
% % % 
% % % % Step 3: Apply the colormap to your data
% % % nrmCl = (tmpCl - value_min) / (value_max - value_min); % Normalize data to [0, 1]
% % % % color_map = colormap(whitetocol(100,[0 0 0.7])); % Apply the colormap to the normalized data
% % % color_map = colormap(redbluecmapRory); % Apply the colormap to the normalized data
% % % colors = interp1(linspace(0, 1, size(color_map, 1)), color_map, nrmCl); % Interpolate colors
% % % 
% % % 
% % % 
% % % scatter3(tmpScore(:,1), tmpScore(:,2), tmpScore(:,3), tmpSz , colors , 'filled'); 
% % % % scatter3(tmpScore(:,4), tmpScore(:,5), tmpScore(:,6), tmpSz , colors , 'filled'); 


%% Compute lines of best fit for each variable of interest
% % % 
% % % dS =  1; % Start dimension 
% % % dE =  3; % End dimension [max 84]
% % % 
% % % [b1 b1Int] = regress(tmpCl, [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);
% % % [b2 b2Int] = regress(tmpSz, [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);
% % % 
% % % figure,histogram2(tmpCl,(b1(end) + sum(b1(1:end-1).*tmpScore(:,dS:dE)')'),'FaceColor','Flat')
% % % figure,histogram2(tmpSz,(b2(end) + sum(b2(1:end-1).*tmpScore(:,dS:dE)')'),'FaceColor','Flat')
% % % 
% % % % figure,plot(tmpCl,(b1(end) + sum(b1(1:end-1).*tmpScore(:,dS:dE)')'),'.')
% % % % figure,plot(tmpSz,(b2(end) + sum(b2(1:end-1).*tmpScore(:,dS:dE)')'),'.')
% % % 
% % % 
% % % dotProd     = dot(b1(1:end-1), b2(1:end-1)) ;
% % % dotProdAbs  = abs(dotProd) ;
% % % cosTheta    = dotProd ./ (norm(b1(1:end-1)) * norm(b2(1:end-1)));
% % % cosThetaAbs = abs(cosThetaAbs);
% % % theta       = acos(cosTheta);


%% Compute similarity of encoding direction as a function of dimension depth
% % % 
% % % nPCs = 3;
% % % 
% % % assessDims = 1:size(tmpScore,2) - (nPCs-1);
% % % 
% % % for iD = assessDims
% % % 
% % %     dS =  iD; % Start dimension
% % %     dE =  iD + nPCs -1; % End dimension [max 84]
% % % 
% % %     [b1 b1Int] = regress(tmpCl, [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);
% % %     [b2 b2Int] = regress(tmpSz, [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);
% % % 
% % %     dotProd(iD)     = dot(b1(1:end-1), b2(1:end-1));
% % %     dotProdAbs(iD)  = abs(dotProd(iD));
% % %     cosTheta(iD)    = dotProd(iD) ./ (norm(b1(1:end-1)) * norm(b2(1:end-1)));
% % %     cosThetaAbs(iD) = abs(cosTheta(iD));
% % %     theta(iD)       = acos(cosTheta(iD));
% % %     thetaAbs(iD)    = abs(theta(iD) - pi./2);
% % % 
% % % end
% % % 
% % % % figure,
% % % % plot(assessDims , similarityMetric,'x','LineWidth',2); hold on
% % % % plot(assessDims , dotProdAbs,'-o','LineWidth',2);
% % % % plot(assessDims , cosTheta,'x','LineWidth',2);
% % % % plot(assessDims , cosThetaAbs,'-x','LineWidth',2);
% % % % plot(assessDims , theta,'-x','LineWidth',2);
% % % % plot(assessDims , thetaAbs,'-o','LineWidth',2);
% % % % legend('DotProd', 'DotProdAbs', 'CosTheta', 'CosThetaAbs','Theta','ThetaAbs')
% % % % xlabel('dimension of PCA')
% % % % ylabel('similarity')
% % % 
% % % 
% % % figure('Position',[50 50 600 900]),
% % % subplot(3,1,1)
% % % plot(assessDims , theta,'-x','LineWidth',2); hold on
% % % plot(assessDims , ones(size(assessDims)) .* pi./2,'-.k');
% % % title('theta')
% % % xlabel('dimension of PCA')
% % % ylabel('angle between vectors (pi/2 is orthogonal)')
% % % 
% % % 
% % % subplot(3,1,2)
% % % plot(assessDims , cosThetaAbs,'-x','LineWidth',2); hold on
% % % plot(assessDims , zeros(size(assessDims)) ,'-.k');
% % % title('costhetaAbs')
% % % xlabel('dimension of PCA')
% % % ylabel('overlap of vectors (0 is orthogonal)')
% % % 
% % % subplot(3,1,3)
% % % plot(assessDims , thetaAbs,'-x','LineWidth',2); hold on
% % % plot(assessDims , ones(size(assessDims)) .* pi./2,'-.k');
% % % title('theta abs diff from pi/2')
% % % xlabel('dimension of PCA')
% % % ylabel('abs angle diff between vectors (0 is orthogonal)')


%% $$$ NEXT change this to loop over all models, and extract similar lines 
%      for all --> see if there is consistent, model-independent pattern


% % % % FOR SMALL NETWORKS
% % % neurTypes   = {'tansig','logsig','softmax','poslin','purelin','tribas','radbas'};
% % % % % % neurTypes   = {'softmax','logsig','tribas','radbas','poslin','purelin','tansig'};
% % % regTypes    = {'Valence_51','L1'};
% % % regNames    = {'No','L1'};
% % % saveName = 'Results\ForFigures\DimensionReduction\FullWorkSpace.mat'

% FOR BIG NETWORKS
neurTypes   = {'SuperCompRelearn'};
regTypes    = {'Valence_51'};
regNames    = {'No'};
saveName = 'Results\ForFigures\DimensionReduction\FullWorkSpace_BigNets.mat';

sAll.plt.ON = 0;



for iTyp = 1:length(neurTypes)
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
                s2.plt.startLayer = 1;
                s2.plt.stopLayer  = numel(s2.lp.netS);

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

         save(['Results\ForFigures\DimensionReduction\NetworkActivations\' ...
          'FullNetActivity_NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(iRun) '_NetArch_' num2str(iM) '.mat'], ...
          'glByThrAct','-v7.3')
  

    end

disp('done:')
iTyp
iReg
iRun

end

save(saveName,'-v7.3')

end
end

%% $$$ Plot out the PCA results for all the models and neuron types

% % % load('Results\ForFigures\DimensionReduction\FullWorkSpace.mat')

clear dotProd dotProdAbs cosTheta cosThetaAbs theta thetaAbs crossProd crossProdSize crossProdRelMax

clear b1 b1ScaleDimSd b1ScaleDimSdExpVar b2 b2ScaleDimSd b2ScaleDimSdExpVar
tic
for iTyp = 1:length(neurTypes)


for iReg = 1:length(regTypes)


bFold = 'Results\ForFigures\Valence\';
cFiles = dir([bFold '*' regTypes{iReg} '*_' neurTypes{iTyp} '_B_V*.mat']);


for iRun = 1:length(cFiles)

for iM = 1:size(rSall,1)

load(['Results\ForFigures\DimensionReduction\NetworkActivations\' ...
          'FullNetActivity_NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(iRun) '_NetArch_' num2str(iM) '.mat']);

% % % load(['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\DimensionReduction\NetworkActivations\' ...
% % %           'FullNetActivity_NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(iRun) '_NetArch_' num2str(iM) '.mat']);

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
    % Replace with zero if the PCA doesn't work
    if size(score,2) == 0
        score = zeros([size(glByThrActFlat,2) 3]);
    end

    

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
            b2ScaleDimSdExpVar(:,iD,iV,iTyp,iReg,iRun,iM) = b2ScaleDimSd(:,iD,iV,iTyp,iReg,iRun,iM) .* sqrt(stats(1));

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
% % %             saveas(gcf,['C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\' ...
% % %                 'Results\ForFigures\DimensionReduction\PCAs\NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(iRun) '_NetArch_' num2str(iM) ] , 'fig')
            saveas(gcf,['Results\ForFigures\DimensionReduction\PCAs\NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(iRun) '_NetArch_' num2str(iM) ] , 'fig')

            close all
        end
        toc

    end


    

disp('done:')
iTyp
iReg
iRun

end
end
end
end

% save(['C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\DimensionReduction\Similarities\' ...
%           'SimilarityScores.mat'], ...
%           'b1','b1ScaleDimSd','b1ScaleDimSdExpVar','b2','b2ScaleDimSd','b2ScaleDimSdExpVar','-v7.3');


% % % save(['Results\ForFigures\DimensionReduction\Similarities\' ...
% % %           'SimilarityScores_BigNetworks.mat'], ...
% % %           'b1','b1ScaleDimSd','b1ScaleDimSdExpVar','b2','b2ScaleDimSd','b2ScaleDimSdExpVar','-v7.3');



%% Convert to useable quantities


% load('Results\ForFigures\DimensionReduction\Similarities\SimilarityScores_V2.mat');
% load('Results\ForFigures\DimensionReduction\Similarities\SimilarityScores_BigNetworks.mat');

TmpNrm  = @(x) squeeze(sqrt(sum(x.^2,1)));

b1Tmp = b1;
b2Tmp = b2;

% b1Tmp = b1ScaleDimSd;
% b2Tmp = b2ScaleDimSd;

% b1Tmp = b1ScaleDimSdExpVar;
% b2Tmp = b2ScaleDimSdExpVar;

dotProd     = squeeze(dot(b1Tmp(1:end-1,:,:,:,:,:,:), b2Tmp(1:end-1,:,:,:,:,:,:)));
dotProdAbs  = abs(dotProd);
cosTheta    = dotProd ./ (TmpNrm(b1Tmp(1:end-1,:,:,:,:,:,:)) .* TmpNrm(b2Tmp(1:end-1,:,:,:,:,:,:)));
cosThetaAbs = abs(cosTheta);
theta       = acos(cosTheta);
thetaAbs    = abs(theta - pi./2);

crossProd       = cross(b1Tmp(1:end-1,:,:,:,:,:,:), b2Tmp(1:end-1,:,:,:,:,:,:));
crossProdSize   = TmpNrm(crossProd);
crossProdRelMax = crossProdSize ./ ...
                   (TmpNrm(b1Tmp(1:end-1,:,:,:,:,:,:)) .* TmpNrm(b2Tmp(1:end-1,:,:,:,:,:,:)));
    

%% $$$ Create permuted b2 --> permute ACROSS RUNS

% % % % (vect coords,iDim,iVar,iTyp,iReg,iRun,iM)
% % % 
% % % tic
% % % 
% % % nPrm = 1000;
% % % 
% % % b2TmpPrm = repmat(b2Tmp,[1,1,1,1,1,1,1, nPrm]);
% % % 
% % % 
% % % % Define the target array size without the last dimension
% % % aS = size(b2TmpPrm);
% % % aS = aS(1:end-1);
% % % 
% % % % Use a nested loop to fill in the last dimension
% % % for iP = 1:nPrm
% % % for i = 1:aS(1), for j = 1:aS(2), for k = 1:aS(3), for l = 1:aS(4), for m = 1:aS(5), for n = 1:aS(7), 
% % %     b2TmpPrm(i, j, k, l, m, :, n, iP) = b2Tmp(i, j, k, l, m, randperm(aS(6)), n);
% % % end, end, end, end, end, end
% % % end
% % % 
% % % toc
% % % 
% % % tic
% % % for iP = 1:nPrm
% % %     crossProdSizePerm(:,:,:,:,:,:,iP) = TmpNrm(cross(b1Tmp(1:end-1,:,:,:,:,:,:), b2TmpPrm(1:end-1,:,:,:,:,:,:,iP)));
% % %     dotProdAbsPerm(:,:,:,:,:,:,iP)    = abs(squeeze(dot(b1Tmp(1:end-1,:,:,:,:,:,:), b2TmpPrm(1:end-1,:,:,:,:,:,:,iP)))  );
% % %     thetaPerm(:,:,:,:,:,:,iP)         = acos(  squeeze(dot(b1Tmp(1:end-1,:,:,:,:,:,:), b2TmpPrm(1:end-1,:,:,:,:,:,:,iP))) ./ ...
% % %         (TmpNrm(b1Tmp(1:end-1,:,:,:,:,:,:)) .* TmpNrm(b2TmpPrm(1:end-1,:,:,:,:,:,:,iP)))  );
% % % end
% % % toc

%% Scatter plots 
% % % load('C:\Users\Rory Bufacchi\OneDrive\Projects\DPPS\DefenseAgent\Results\ForFigures\DimensionReduction\FullWorkSpace_BigNets.mat')

% % % iMM = 1;
% % % for iRun = 1:size(rSall,2)
% % %     for iM = 1:size(rSall,1)
% % %         tmpPerf(iMM) = sum(rSall(iM,iRun).perf.rewPerAct(end,:));
% % %         iMM = iMM + 1;
% % %     end
% % % end


    %% Line plots 

% % %     % Reshape in case only 1 regularization and type is present (i.e. if
% % %     % only looking at the big networks)
% % %     if size(crossProdSize,5) <= 1
% % %         crossProdSize       = permute(crossProdSize,[1 2 5 6 3 4 7]);
% % %         dotProdAbs          = permute(dotProdAbs,[1 2 5 6 3 4 7]);
% % %         crossProdSizePerm   = permute(crossProdSizePerm,[1 2 5 6 3 4 7 8]);
% % %         dotProdAbsPerm      = permute(dotProdAbsPerm,[1 2 5 6 3 4 7 8]);
% % %     end
% % % 
% % % % % %     % For small networks
% % % % % %     iV   = 5; %1; % comparison Variables % $$$ COMAPRISON 6, MIN VS MAX STIMDIST is actually CRAZY ALIGNED for early PCs (also 7, min vs avstimdist)
% % % % % %     iTyp = 1; % Neuron type
% % % % % %     iReg = 2; % Regularisation type
% % % 
% % %     % For big networks
% % %     iV   = 1; %1; % comparison Variables % $$$ COMAPRISON 6, MIN VS MAX STIMDIST is actually CRAZY ALIGNED for early PCs (also 7, min vs avstimdist)
% % %     iTyp = 1; % Neuron type
% % %     iReg = 1; % Regularisation type
% % % 
% % % 
% % % 
% % %     % Compute similarity of encoding direction as a function of dimenstion depth
% % %     allVars     = {'thrRow','thrCol','lmbCol','lmbCol','golDist',...
% % %                'minStimDist','minStimDist'};
% % %     allVarsComp = {'golRow','golCol','thrCol','golCol','thrDist',...
% % %                'maxStimDist','avStimDist'};
% % % 
% % % 
% % % allVars{iV}
% % % allVarsComp{iV}
% % % neurTypes{iTyp}
% % % 
% % %     assessDims = 1:size(b1,2);
% % % 
% % % 
% % % 
% % % % figure('Position',[2000 50 600 900]),
% % % figure('Position',[50 50 600 900]),
% % % 
% % % % -------------
% % % subplot(4,1,1)
% % % plot(assessDims , squeeze(theta(:,iV,iTyp,iReg,:))); hold on
% % % plot(assessDims , squeeze(nanmean(theta(:,iV,iTyp,iReg,:),5)) ,'-k','LineWidth',2); 
% % % plot(assessDims , ones(size(assessDims)) .* pi./2,'-.k');
% % % title('theta')
% % % xlabel('dimension of PCA')
% % % ylabel('angle between vectors (pi/2 is orthogonal)')
% % % 
% % % % -------------
% % % subplot(4,1,2)
% % % plot(assessDims , squeeze(cosThetaAbs(:,iV,iTyp,iReg,:))); hold on
% % % plot(assessDims , squeeze(nanmean(cosThetaAbs(:,iV,iTyp,iReg,:),5)),'-xk','LineWidth',2); 
% % % plot(assessDims , zeros(size(assessDims)) ,'-.k');
% % % title('costhetaAbs')
% % % xlabel('dimension of PCA')
% % % ylabel('overlap of vectors (0 is orthogonal)')
% % % 
% % % % -------------
% % % subplot(4,1,3)
% % % plot(assessDims , squeeze(thetaAbs(:,iV,iTyp,iReg,:))); hold on
% % % plot(assessDims , squeeze(nanmean(thetaAbs(:,iV,iTyp,iReg,:),5)),'-xk','LineWidth',2); 
% % % plot(assessDims , ones(size(assessDims)) .* pi./2,'-.k');
% % % title('theta abs diff from pi/2')
% % % xlabel('dimension of PCA')
% % % ylabel('abs angle diff between vectors (0 is orthogonal)')
% % % 
% % % % -------------
% % % subplot(4,1,4)
% % % plot(assessDims , squeeze(dotProdAbs(:,iV,iTyp,iReg,:))); hold on
% % % plot(assessDims , squeeze(nanmean(dotProdAbs(:,iV,iTyp,iReg,:),5)),'-xk','LineWidth',2); 
% % % plot(assessDims , zeros(size(assessDims)) ,'-.k');
% % % 
% % % tmpD = nanmean(dotProdAbsPerm(:,iV,iTyp,iReg,:,:,:,:),7); % $$$
% % % plot(assessDims , squeeze(nanmedian(tmpD(:,1,1,1,:),5)  ) ,'-b','LineWidth',3); hold on
% % % 
% % % title('dot product abs')
% % % xlabel('dimension of PCA')
% % % ylabel('dot product between vectors (0 is orthogonal)')
% % % ylim([-1 1])
% % % 
% % % 
% % % figure('Position',[800 50 600 900]),
% % % % -------------
% % % subplot(4,1,1)
% % % plot(assessDims , squeeze(crossProdRelMax(:,iV,iTyp,iReg,:))); hold on
% % % plot(assessDims , squeeze(nanmean(crossProdRelMax(:,iV,iTyp,iReg,:),5)) ,'-k','LineWidth',2); 
% % % plot(assessDims , ones(size(assessDims)) .* 1,'-.k');
% % % title('xprod relative to maximum')
% % % xlabel('dimension of PCA')
% % % ylabel('xProd')
% % % 
% % % % -------------
% % % subplot(4,1,2)
% % % plot(assessDims , squeeze(crossProdSize(:,iV,iTyp,iReg,:))); hold on
% % % plot(assessDims , squeeze(nanmedian(crossProdSize(:,iV,iTyp,iReg,:),5)) ,'-k','LineWidth',2); 
% % % plot(assessDims , ones(size(assessDims)) .* 0 ,'-.k');
% % % 
% % % tmpD = nanmean(crossProdSizePerm(:,iV,iTyp,iReg,:,:,:,:),7); % $$$
% % % plot(assessDims , squeeze(nanmedian(tmpD(:,1,1,1,:),5)  ) ,'-b','LineWidth',3); hold on
% % % 
% % % title('xprod size')
% % % xlabel('dimension of PCA')
% % % ylabel('xProd')
% % % 
% % % % -------------
% % % subplot(4,1,3)
% % % tmpD = abs(squeeze(crossProdSize(:,iV,iTyp,iReg,:))) ./ (  abs(squeeze(crossProdSize(:,iV,iTyp,iReg,:)))  +  abs(squeeze(dotProd(:,iV,iTyp,iReg,:)))  );
% % % plot(assessDims ,tmpD); hold on
% % % plot(assessDims , nanmean(tmpD,2) ,'-k','LineWidth',2); 
% % % plot(assessDims , ones(size(assessDims)) .* 0 ,'-.k');
% % % title('abs(xprod) vs abs(xprod) + abs(dotprod)')
% % % xlabel('dimension of PCA')
% % % ylabel('xProd')



%% Plots for each neuron type


% binEdges = -10:2.5:10;
% binEdges = -5:5;



cOrder      = {[0 0 .7], [.7 0 0], [.8 .6 0]};


% % % % FOR SMALL NETWORKS
% % % neurTypes   = {'tansig','logsig','softmax','poslin','purelin','tribas','radbas'};
% % % regTypes    = {'Valence_51','L1'};
% % % regNames    = {'No','L1'};
% % % az = -74.3145;
% % % el = 43.0968;
% % % binEdges = -15:2.5:15;
% % % f.NeurTypesPCA.f         = figure('Position',[20 -20 900 900]);
% % % % f.NeurTypesPCAAngles.f   = figure('Position',[20 -20 900 900]);
% % % f.NeurTypesPCAExamplesLeast.f = figure('Position',[20 -20 900 900]);
% % % f.NeurTypesPCAExamplesMost.f  = figure('Position',[20 -20 900 900]);
% % % f.NeurTypesPCAExamplesLeastSubSection.f = figure('Position',[20 -20 1600 600]);

% FOR BIG NETWORKS
neurTypes   = {'SuperCompRelearn'};
regTypes    = {'Valence_51'};
regNames    = {'No'};
az = 39.0563;
el = 52.8667;
binEdges = -25:2.5:7.5;
f.NeurTypesPCA.f         = figure('Position',[20 -20 450 900]);
% f.NeurTypesPCAAngles.f   = figure('Position',[20 -20 900 900]);
f.NeurTypesPCAExamplesLeast.f = figure('Position',[20 -20 600 900]);
f.NeurTypesPCAExamplesMost.f  = figure('Position',[20 -20 600 900]);
f.NeurTypesPCAExamplesLeastSubSection.f = figure('Position',[20 -20 1600 600]);


s.wrld.size = [14 15];

cM   = 1:3;

iV    = 1; % 1 % comparison Variables % $$$ COMAPRISON 6, MIN VS MAX STIMDIST is actually CRAZY ALIGNED for early PCs (also 7, min vs avstimdist)

iPl  = 1;
iPl2 = 1;
iPl3 = 1;

for iTyp = 1:length(neurTypes) % Neuron type
    for iReg =  1:length(regTypes) % Regularisation type



        figure(f.NeurTypesPCA.f);
% % %         subplot(numel(neurTypes),numel(regNames),iPl)
% % % 
% % %         for iM = cM
% % %             % $$$ Here think about how to properly plot the effects? Maybe just two
% % %             % distributions? That could be neat, one blue, one red?
% % %             %         tmpPM = nanmean(crossProdSizePerm(1,iV,iTyp,iReg,:,iM,:),7);
% % %             %         tmpPM = crossProdSizePerm(1,iV,iTyp,iReg,:,iM,:);
% % %             %         tmpPM = crossProdSize(1,iV,iTyp,iReg,:,iM);
% % % 
% % % % % % % % %                         tmpPM = crossProdSize(1,iV,iTyp,iReg,:,iM) - crossProdSizePerm(1,iV,iTyp,iReg,:,iM,:);
% % % % % %             tmpPM = dotProdAbs(1,iV,iTyp,iReg,:,iM) - nanmean(dotProdAbsPerm(1,iV,iTyp,iReg,:,iM,:),7);
% % %             tmpPM = dotProdAbs(1,iV,iTyp,iReg,:,iM) - dotProdAbsPerm(1,iV,iTyp,iReg,:,iM,:);
% % % 
% % %             tmpPM = tmpPM(:);
% % % 
% % %             h = histogram(tmpPM(:),'Normalization','probability','BinEdges',binEdges,...
% % %                 'FaceColor',cOrder{iM},'FaceAlpha',0.3, ...
% % %                 'EdgeColor','k','EdgeAlpha', 0.2 ); hold on
% % %             % % %                     h = histogram(tmpPM(:),'Normalization','probability',...
% % %             % % %                         'FaceColor',cOrder{iM},'FaceAlpha',0.3, ...
% % %             % % %                         'EdgeColor','k','EdgeAlpha', 0.2 ); hold on
% % %             % % %                             tmpPM2(1:numel(tmpPM(:)),iM) = tmpPM(:);
% % %             hold on
% % %             splineH{iM} = PlotHistSpline(h,'LineWidth',2,'Color',cOrder{iM});
% % %             % % %         xlim([binEdges(1) binEdges(end)])
% % %             xlim([binEdges(1) -binEdges(1)]);
% % %             ylim([0 0.7]);
% % %             plot([0 0],[0 0.7],'-.k')
% % %         end

        % % %         figure(f.NeurTypesPCAAngles.f);
        % % %         subplot(numel(neurTypes),numel(regNames),iPl)
        % % %
        % % %         for iM = cM
        % % %             tmpTPM = real(squeeze(thetaPerm(1,iV,iTyp,iReg,:,iM,:))) ;
        % % %             polarhistogram(tmpTPM(:),12,'Normalization','probability'); hold on
        % % %             tmpT   = squeeze(theta(1,iV,iTyp,iReg,:,iM));
        % % %             polarhistogram(tmpT,12,'Normalization','probability'); hold on
        % % %         end



        iPl = iPl + 1;

% % %         tmpPMforChoosePlot = squeeze(crossProdSize(1,iV,iTyp,iReg,:,:) - mean(crossProdSizePerm(1,iV,iTyp,iReg,:,:,:),7));
        tmpPMforChoosePlot = squeeze(dotProdAbs(1,iV,iTyp,iReg,:,:) - mean(dotProdAbsPerm(1,iV,iTyp,iReg,:,:,:),7));

        for iExampFig = 1:2
            if iExampFig == 1
                figure(f.NeurTypesPCAExamplesLeast.f);
                tmpMin = nanmin(squeeze(tmpPMforChoosePlot),[],'all');
            elseif iExampFig == 2
                figure(f.NeurTypesPCAExamplesMost.f);
                tmpMin = nanmax(squeeze(tmpPMforChoosePlot),[],'all');
            end

            %         tmpPMforChoosePlot = squeeze(dotProdAbs(1,iV,iTyp,iReg,:,:) - mean(dotProdAbsPerm(1,iV,iTyp,iReg,:,:,:),7));
            %         tmpMin = nanmin(squeeze(tmpPMforChoosePlot),[],'all');
            %         tmpMin = nanmedian(squeeze(tmpPMforChoosePlot),'all');
            [cRun cM] = find(tmpMin == tmpPMforChoosePlot,1);

%             cRun = 3; cM = 2;

            if isempty(cRun)
                cRun = 1;
                cM   = 1;
            end

% % %             load(['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\DimensionReduction\NetworkActivations\' ...
% % %                 'FullNetActivity_NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(cRun) '_NetArch_' num2str(cM) '.mat']);

            load(['Results\ForFigures\DimensionReduction\NetworkActivations\' ...
                'FullNetActivity_NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(cRun) '_NetArch_' num2str(cM) '.mat']);


            % Create flattened activations for PCA
            glByThrActFlat = permute(glByThrAct,[1 2 3 6 7 4 5]);
            unwrapSize     = size(glByThrActFlat);
            glByThrActFlat = permute(glByThrActFlat(:,:,:,:,:,:),[6 1 2 3 4 5]);
            glByThrActFlat = glByThrActFlat(:,:);

            % Create 'explanatory variables'
            thrRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[1 2 3 7]) ] );
            thrRow = permute(thrRow, [2 3 4 1 5]);

            golRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[2 3 6 7]) ] );
            golRow = permute(golRow, [1 2 3 4 5]);


            % Perform PCA
            [coeff,score,~,~,explained] = pca(glByThrActFlat( ~isnan(glByThrActFlat(:,1)),:)');
            % Replace with zero if the PCA doesn't work
            if size(score,2) == 0
                score = zeros([size(glByThrActFlat,2) 3]);
            end


            % Downsample only a little
            tmpScore = reshape(score', [size(score,2) unwrapSize(1:5)] );
            tmpScore = tmpScore(:,:,2:2:end,1:2:end,:,2:2:end);
            tmpScore = tmpScore(:,:)';

            for iVar = 1:2

                if iExampFig == 1
                    subplot(numel(neurTypes), numel(regNames).*2, iPl2);
                    iPl2 = iPl2 + 1;
                elseif iExampFig == 2
                    subplot(numel(neurTypes), numel(regNames).*2, iPl3);
                    iPl3 = iPl3 + 1;
                end

                if iVar == 1
                    tmpCl = golRow(:,2:2:end,1:2:end,:,2:2:end);
                elseif iVar == 2
                    tmpCl = thrRow(:,2:2:end,1:2:end,:,2:2:end);
                end
                tmpCl = tmpCl(:);
                tmpSz = 5;

                value_min = min(tmpCl);
                value_max = max(tmpCl);

                colormap_name = 'jet'; % You can use any built-in colormap or create a custom one

                nrmCl = (tmpCl - value_min) / (value_max - value_min); % Normalize data to [0, 1]
% % %                 color_map = colormap(redbluecmapRory); % Apply the colormap to the normalized data
                
                if iVar == 1
%                     color_map = colormap(coltocol(100,[0 0.7 0],[0 0 0.7])); % Apply the colormap to the normalized data 
                    color_map = colormap(coltocol(100,[0.2 0.2 0.2],[0 0 0.7])); % Apply the colormap to the normalized data 
                elseif iVar == 2
                    color_map = colormap(coltocol(100,[0.2 0.2 0.2],[0.7 0 0])); % Apply the colormap to the normalized data 
                end

                % $$$ HERE HERE
%                 color_map = colormap(coltocol(100,[0 1 0], [0 0 .7 ],[1]));
%                 color_map = colormap(redbluecmapRory); 
%                 color_map = color_map(:,[1 3 2]);
%                 color_map = color_map(:,[2 3 1]);

                
                colors = interp1(linspace(0, 1, size(color_map, 1)), color_map, nrmCl); % Interpolate colors

                scatter3(tmpScore(:,1), tmpScore(:,2), tmpScore(:,3), tmpSz , colors , 'filled');
                %             scatter(tmpScore(:,1), tmpScore(:,2), tmpSz , colors , 'filled');

                if iVar == 1
                    title('Goal row')
                elseif iVar == 2
                    title('Threat row')
                end

                
                % Make clearer mini-example
                if iExampFig == 1 & iTyp == 1 & iReg == 1
                    figure(f.NeurTypesPCAExamplesLeastSubSection.f);

                    if iVar == 1
                        subplot(1,2,1);
                    elseif iVar == 2
                        subplot(1,2,2);
                    end

                    scatter3(tmpScore(:,1), tmpScore(:,2), tmpScore(:,3), tmpSz , colors , 'filled');

                    view([az el]);
% % %                     xlim([-3 3.8]);
% % %                     ylim([-4.5, 2]);

                    if iVar == 1
                        title('Goal row')
                    elseif iVar == 2
                        title('Threat row')
                    end

                    cBarLims(iVar,1) = min(nrmCl(:));
                    cBarLims(iVar,2) = max(nrmCl(:));

                    % % %                     figure(f.NeurTypesPCAExamplesColourBars.f)
                    % % %                     subplot(1,2,iVar)
                    colormap(gca, color_map);
                    caxis(cBarLims(iVar,:));
                    colorbar;

                    % Move back to the other figure
                    if iExampFig == 1
                        figure(f.NeurTypesPCAExamplesLeast.f);
                        tmpMin = nanmin(squeeze(tmpPMforChoosePlot),[],'all');
                    elseif iExampFig == 2
                        figure(f.NeurTypesPCAExamplesMost.f);
                        tmpMin = nanmax(squeeze(tmpPMforChoosePlot),[],'all');
                    end
                end
                
            end

            if iExampFig == 1
                sgtitle('Least orthogonal examples');
            elseif iExampFig == 2
                sgtitle('most orthogonal examples');
            end

        end


    end
end


%% Collect data for performance vs PCA orthogonality

% % % load('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\NeurTypesAnalysed_V5c.mat')
% % % load('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\NeurTypesAnalysed_BIGNETS_temp.mat')

% FOR SMALL NETWORKS
neurTypes   = {'tansig','logsig','softmax','poslin','purelin','tribas','radbas'};
regTypes    = {'Valence_51','L1'};
regNames    = {'No','L1'};

% % % 
% % % % FOR BIG NETWORKS
% % % neurTypes   = {'SuperCompRelearn'};
% % % regTypes    = {'Valence_51'};
% % % regNames    = {'No'};

tic
for iTyp = 1:length(neurTypes) % Neuron type
    for iReg = 1:numel(regTypes) % Regularisation type



        % Prepare base files for loading
        bFold = 'Results\ForFigures\Valence\';
        cFiles = dir([bFold '*' regTypes{iReg} '*_' neurTypes{iTyp} '_B_V*.mat']);
        
        for iRun = 1:length(cFiles)
            load([bFold cFiles(iRun).name]);
            for iM = 1:numel(rS)
                rewPerAct(iTyp,iReg,iRun,iM) = sum(rS(iM).perf.rewPerAct(end,:));
            end
        end
    end
    iTyp
end
toc


%% Effects of neuron type
tic
iM = 1;

% FOR SMALL NETWORKS
neurTypes   = {'tansig','logsig','softmax','poslin','purelin','tribas','radbas'};
% % % neurTypes   = {'softmax','logsig','tribas','radbas','poslin','purelin','tansig'};
regTypes    = {'Valence_51','L1'};
regNames    = {'No','L1'};

% 
% % FOR BIG NETWORKS
% neurTypes   = {'SuperCompRelearn'};
% regTypes    = {'Valence_51'};
% regNames    = {'No'};


for iTyp = 1:length(neurTypes);
for iReg = 1:length(regTypes);


bFold = 'Results\ForFigures\Valence\';
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
    rSall(iM,iRun,iTyp,iReg).s.plt.nPm = 1 ; % $$$ change to 100 if necessary later

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

% % % % Save analysed neuron types
% % % save('Results\ForFigures\NeurTypesAnalysed_V5.mat','-v7.3')
% % % save('Results\ForFigures\NeurTypesAnalysed_BIGNETS_temp.mat','-v7.3')
toc

%% Calculate necessary variables for performing the LMEs

tic
clear layNum layNumDiff nodeDists stimPrefSim
for iTyp = 1:length(neurTypes);
for iReg = 1:length(regTypes);
for iM = 1:size(allNetAR,1)
    for iRun = 1:size(allNetAR,2)

    netAR   = allNetAR(iM,iRun,iTyp,iReg);
    rS      = rSall(iM,iRun,iTyp,iReg);

    % Find performance of network
    allPerf(iM,iRun,iTyp,iReg) = sum(rS(1).perf.rewPerAct(end,:));
    allNetW(iM,iRun,iTyp,iReg) = rS(1).s.lp.netS(end);

    % Make an indicator of layer depth
    tmp = netAR.A_B_rat';
    layNumTmp = [1:size(netAR.A_B_rat,1)] .*    ones(size(tmp));
    layNumTmp = layNumTmp(:);
    layNumTmp(isnan(tmp(:))) = [];
    layNum(:,iM,iRun,iTyp,iReg) = layNumTmp;
    layNumDiff(:,:,iM,iRun,iTyp,iReg)  = abs(bsxfun(@minus, layNumTmp', layNumTmp));

    A_B_rat = NanRemFlatten(netAR.A_B_rat'); classType = 'rat';
    G = netAR.G;

    p = plot(G,'Layout','force','WeightEffect','inverse','UseGravity','on');

    % First make a distance matrix
    for iNode=1:size(G.Nodes,1)
        nodeDists(:,iNode,iM,iRun,iTyp,iReg)    = sqrt(sum(([p.XData(iNode),p.YData(iNode)] - [p.XData(:),p.YData(:)])'.^2));
    end
% % %     % Also permute the distance matrix
% % %     for iPm = 1:100
% % %         permutedNodes = randperm(numel(p.XData));
% % %         for iNode=1:size(G.Nodes,1)
% % %             nodeDistsPm(:,iNode,iM,iRun,iTyp,iReg,iPm)    = sqrt(sum(([p.XData(permutedNodes(iNode)),p.YData(permutedNodes(iNode))] - [p.XData(:),p.YData(:)])'.^2));
% % %         end
% % %     end

    % Then make a stimulus preference matrix
    stimPrefSimTmp = 1 - abs(A_B_rat - A_B_rat');
    stimPrefSim(:,:,iM,iRun,iTyp,iReg)  = stimPrefSimTmp;

    % Calculate the correlation between the two
    nodeDistsTmp = nodeDists(:,:,iM,iRun,iTyp,iReg);
    [rho(iM,iRun,iTyp,iReg) pval(iM,iRun,iTyp,iReg)] = corr(stimPrefSimTmp(nodeDistsTmp ~= 0),nodeDistsTmp(nodeDistsTmp ~= 0));
%     [rho(iM) pval(iM)] = corr(stimPrefSimTmp(:),nodeDistsTmp(:));
    % $$$ I need to account for stimulus proximity being 0

    end
end
end
end
toc

% pvalTmp = pval(:,1:15,:,:);
pvalTmp = pval(:,1:15,1,1);
[pValThr, pValcor, pValAdj] = fdr(pvalTmp(:));


fprintf('%i out of %i networks show structure using this metric\n', sum(pValAdj(:) < 0.05), numel(pValAdj));

disp('mean rho')
nanmean(rho(:))
disp('std rho')
nanstd(rho(:))

%% Perform LME for all models independently, accounting for layer differences

% % % load('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\NeurTypesAnalysed_V5c.mat')
% % % load('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\NeurTypesAnalysed_BIGNETS_temp.mat')


overallModelNum    = repmat(  permute( reshape( 1:numel(allNetAR) , size(allNetAR) ), [5 6 1 2 3 4]) , [size(stimPrefSim,[1 2]) 1 1 1 1]);
layNumOrg   = repmat(permute(layNum,[1 6 2 3 4 5]) , [1 size(layNum,1) 1 1]);

allPerfRep  = repmat(permute(allPerf,[5 6 1 2 3 4 ]) , [size(layNumOrg,[1 2]) 1 1 1 1]);
allNetWRep  = repmat(permute(allNetW,[5 6 1 2 3 4 ]) , [size(layNumOrg,[1 2]) 1 1 1 1]);

allegedNoStructure      = pValAdj > 0.05;
allegedNoStructureRep   = repmat(permute(allegedNoStructure,[5 6 1 2 3 4 ]) , [size(layNumOrg,[1 2]) 1 1 1 1]);

% iM,iRun,iTyp,iReg
allMod = arrayfun(@(i) i .* ones(size(layNumOrg(:,:,1,:,:,:)) ) , 1:size(layNumOrg,3) ,'UniformOutput', false );
allMod = cat(3,allMod{:});
allRun = arrayfun(@(i) i .* ones(size(layNumOrg(:,:,:,1,:,:)) ) , 1:size(layNumOrg,4) ,'UniformOutput', false );
allRun = cat(4,allRun{:});
allTyp = arrayfun(@(i) i .* ones(size(layNumOrg(:,:,:,:,1,:)) ) , 1:size(layNumOrg,5) ,'UniformOutput', false);
allTyp = cat(5,allTyp{:});
allReg = arrayfun(@(i) i .* ones(size(layNumOrg(:,:,:,:,:,1)) ) , 1:size(layNumOrg,6) ,'UniformOutput', false );
allReg = cat(6,allReg{:});


tbl                             = table(stimPrefSim(:));
tbl.Properties.VariableNames    = {'StimulusPreferenceSimilarity'};
% % % tbl.netOutputLayerWidth         = allNetW'; 
tbl.nodeDists                   = nodeDists(:);

% tbl.layNumDiff                  = abs(layNumDiff(:));
% tbl.layNumDiff                  = layNumDiff(:);
tbl.layNumDiff                  = categorical(abs(layNumDiff(:)));

tbl.layNumOrg                   = layNumOrg(:);
tbl.modelNum                    = overallModelNum(:);

tbl.allPerf                     = allPerfRep(:);
tbl.allNetW                     = categorical(allNetWRep(:));

tbl.allegedNoStructure          = allegedNoStructureRep(:);

tbl.allMod                      = allMod(:);
tbl.allRun                      = allRun(:);
tbl.allTyp                      = allTyp(:);
tbl.allReg                      = allReg(:);


tmpTbl = tbl(tbl.nodeDists ~= 0, :);
% tmpTbl = tbl(tbl.nodeDists ~= 0 & tbl.layNumDiff == 0,:);
% tmpTbl = tbl(tbl.nodeDists ~= 0 & tbl.layNumDiff == 0 & tbl.layNumOrg  == 4,:);
% tmpTbl = tbl(tbl.nodeDists ~= 0 & tbl.layNumOrg  == 4,:);

% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists*layNumDiff + (1|modelNum)')
% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists + layNumDiff + (1|modelNum)')
% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists + (1|modelNum)')
% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists + allNetW + allPerf + (1|modelNum)')

% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists*allTyp*allReg + (1|modelNum)')

% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists + nodeDists:allNetW  + (1|modelNum)')

% anova(lme)

% Run LME for each network
tic
for iTyp = 1:length(neurTypes);
for iReg = 1:length(regTypes);
for iM = 1:size(allNetAR,1)
    for iRun = 1:size(allNetAR,2)
    tmpTbl          = tbl(tbl.nodeDists ~= 0 & ...
                          tbl.allMod == iM   & ...
                          tbl.allRun == iRun   & ...
                          tbl.allReg == iReg   & ...
                          tbl.allTyp == iTyp   ,:);

tmpTbl.allMod = categorical(tmpTbl.allMod);
tmpTbl.allRun = categorical(tmpTbl.allRun);
tmpTbl.allTyp = categorical(tmpTbl.allTyp);
tmpTbl.allReg = categorical(tmpTbl.allReg);


%     tmpLme          = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists*layNumDiff + (1|modelNum)');
    tmpLme          = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ layNumDiff:nodeDists + nodeDists + (1|modelNum)');


% % %     tmpTbl          = tbl(tbl.nodeDists ~= 0 & tbl.modelNum == iM & tbl.layNumDiff == 0,:);
% % %     tmpLme          = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists + (1|modelNum)');
    tmpAnova                = anova(tmpLme);
    estimateDistEff(iM,iRun,iTyp,iReg)     = tmpLme.Coefficients{2,2};
    estimateStandErr(iM,iRun,iTyp,iReg)    = tmpLme.Coefficients{2,3};
    tStatsStructure(iM,iRun,iTyp,iReg)     = tmpLme.Coefficients{2,4};
    pNodeDists(iM,iRun,iTyp,iReg)          = tmpLme.Coefficients{2,6};

% % %     % Also do it for a permuted network, just to check the null distribution isn't way away from 0
% % %     for iPm = 1:100
% % %          nodeDistsTmp     = nodeDistsPm(:,:,:,iPm);
% % %          nodeDistsTmp     = nodeDistsTmp(:);
% % %          tmpTbl.nodeDists = nodeDistsTmp(tbl.nodeDists ~= 0 & tbl.modelNum == iM);
% % %          tmpLme           = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists*layNumDiff + (1|modelNum)');
% % %          
% % %          pNodeDistsPm(iM,iPm)       = tmpLme.Coefficients{2,6};
% % %          estimateDistEffPm(iM,iPm)  = tmpLme.Coefficients{2,2};
% % %     end
    end
end
end
end
toc

pTmp        = pNodeDists(:,1:15);
estTmp      = estimateDistEff(:,1:15);
estErrTmp   = estimateStandErr(:,1:15);
tStatTmp    = tStatsStructure(:,1:15); 

[pValThr, pValcor, pValLMEAdj] = fdr(pTmp(:));


fprintf('%i out of %i networks show structure using this metric\n', sum(pValLMEAdj(:) < 0.05), numel(pValLMEAdj));

dimsToSum = [1 2];

disp('Weighted mean estimate of node distance')
meanEffs = squeeze(sum(  (estTmp ./ estErrTmp.^2) , dimsToSum) ./ sum(1 ./ estErrTmp.^2  , dimsToSum ))
disp('standard error of effect of node distance')
sdEffs   = squeeze(sqrt(1 ./ sum(1 ./ estErrTmp.^2  , dimsToSum )))


disp('Weighted mean tstat of node distance')
meanEffs = mean(tStatTmp(:))
disp('standard deviation of tstat of node distance')
sdEffs   = std(tStatTmp(:))

% % % save('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\NeurTypesAnalysed_V5.mat','-v7.3')

%% $$$ Use the projection method to perform parametric stats for the PCA?

% % % % First project onto the goal and threat axes
% % % 
% % % % Select the particalar comparisons and PCA dimensions of interest
% % % iV = 1;
% % % iD = 1;
% % % 
% % % for iTyp = 1:length(neurTypes)
% % % 
% % % 
% % % for iReg = 1:length(regTypes)
% % % 
% % % 
% % % bFold = 'Results\ForFigures\Valence\';
% % % cFiles = dir([bFold '*' regTypes{iReg} '*_' neurTypes{iTyp} '_B_V*.mat']);
% % % 
% % % 
% % % for iRun = 1:length(cFiles)
% % % 
% % % for iM = 1:size(rSall,1)
% % % 
% % % load(['Results\ForFigures\DimensionReduction\NetworkActivations\' ...
% % %           'FullNetActivity_NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(iRun) '_NetArch_' num2str(iM) '.mat']);
% % % 
% % % % % % load(['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\DimensionReduction\NetworkActivations\' ...
% % % % % %           'FullNetActivity_NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(iRun) '_NetArch_' num2str(iM) '.mat']);
% % % 
% % % % $$$ HERE figure out how to deal with the bloody model architectures --> I
% % % % think just store them separately?
% % % 
% % %   % Create flattened activations for PCA
% % %     glByThrActFlat = permute(glByThrAct,[1 2 3 6 7 4 5]);
% % %     unwrapSize     = size(glByThrActFlat);
% % %     glByThrActFlat = permute(glByThrActFlat(:,:,:,:,:,:),[6 1 2 3 4 5]);
% % %     glByThrActFlat = glByThrActFlat(:,:);
% % % 
% % %     
% % % % % %     % ============================I THINK THIS CAN GO =====================
% % % % % %     % Create 'explanatory variables'
% % % % % %     thrRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[1 2 3 7]) ] );  
% % % % % %     thrRow = permute(thrRow, [2 3 4 1 5]);
% % % % % %     
% % % % % %     golRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[2 3 6 7]) ] );  
% % % % % %     golRow = permute(golRow, [1 2 3 4 5]);   
% % % % % %     % ============================I THINK THIS CAN GO =====================
% % % 
% % % 
% % %     % Perform PCA
% % %     [coeff,score,~,~,explained] = pca(glByThrActFlat( ~isnan(glByThrActFlat(:,1)),:)');
% % %     % Replace with zero if the PCA doesn't work
% % %     if size(score,2) == 0
% % %         score = zeros([size(glByThrActFlat,2) 3]);
% % %     end
% % % 
% % % 
% % % 
% % %     
% % %     % Downsample more here than when creating the basis vectors, because
% % %     % otherwise we end up with a stupidly large array of vectors between
% % %     % states
% % %     tmpScore = reshape(score', [size(score,2) unwrapSize(1:5)] );
% % %     tmpScore = tmpScore(:,2:3:end,2:3:end,1:3:end,2:3:end,2:3:end);
% % %     tmpScore = tmpScore(:,:)';
% % % 
% % % 
% % % 
% % %     % Extract attributes of interest for given state
% % %     allVars     = {'thrRow','thrCol','lmbCol','lmbCol','golDist',...
% % %                'minStimDist','minStimDist'};
% % %     allVarsComp = {'golRow','golCol','thrCol','golCol','thrDist',...
% % %                'maxStimDist','avStimDist'};
% % %     eval(['tmpCl    = ' allVars{iV} '(2:3:end,2:3:end,1:3:end,2:3:end,2:3:end);']);
% % %     eval(['tmpSz    = ' allVarsComp{iV} '(2:3:end,2:3:end,1:3:end,2:3:end,2:3:end);']);
% % %     tmpCl    = tmpCl(:);
% % %     tmpSz    = tmpSz(:);
% % % 
% % % 
% % %     % Create a vector between each pair of points in the 1st 3 dimensions
% % %     % of the score
% % %     dS = 1;
% % %     dE = 3;
% % %     pointToPointVecs = tmpScore(:,dS:dE)' - permute(tmpScore(:,dS:dE)',[1 3 2]);
% % %     % Create a vector between each pair of points in goal-threat (fully
% % %     % orthogonal) space
% % %     goalThreatSpace         = [tmpCl, tmpSz];
% % %     pointToPointVecsGolThr  = goalThreatSpace' - permute(goalThreatSpace',[1 3 2]);
% % % 
% % %     tic
% % %     % Remove diagonal (because distance is 0), and diagonally symmetric elements
% % %     for iR = 1:size(pointToPointVecsGolThr,2)
% % %         pointToPointVecs(:,iR:end,iR) = NaN;
% % %         pointToPointVecsGolThr(:,iR:end,iR) = NaN;
% % %     end
% % %     pointToPointVecs        = pointToPointVecs(:,:);
% % %     pointToPointVecsGolThr  = pointToPointVecsGolThr(:,:);
% % %     pointToPointVecs(:,isnan(pointToPointVecs(1,:))) = [];
% % %     pointToPointVecsGolThr(:,isnan(pointToPointVecsGolThr(1,:))) = [];
% % %     toc
% % % 
% % %     % $$$ Calculate dot product between each pair of vectors? Will that
% % %     % make it huge? I think that will make it too huge...
% % % 
% % % 
% % % %     figure, histogram2(sqrt(sum(pointToPointVecsGolThr.^2)),sqrt(sum(pointToPointVecs.^2)),'Facecolor','flat')
% % % 
% % %     figure, histogram2(pointToPointVecsGolThr(1,:),projOntoGoal,'Facecolor','flat')
% % % 
% % %     figure, histogram2(pointToPointVecsGolThr(2,:),projOntoThreat,'Facecolor','flat')
% % % 
% % % 
% % % 
% % % 
% % %     % Remove zero distance
% % %     pointToPointVecsGolThr(:,sum(pointToPointVecs,1) == 0) = [];
% % %     pointToPointVecs(:,sum(pointToPointVecs,1) == 0) = [];
% % % 
% % % 
% % % 
% % % 
% % %     % Normalise
% % %     pointToPointVecs = pointToPointVecs ./ sqrt(sum(pointToPointVecs.^2 ,1));
% % %     
% % % 
% % %     % Project that vector onto the goal and threat-preferring axes
% % %     % i.e. dot product
% % %     tmpGlVec        = b1ScaleDimSd(1:3,iD,iV,iTyp,iReg,iRun,iM) ./ sqrt(sum(b1ScaleDimSd(1:3,iD,iV,iTyp,iReg,iRun,iM).^2,1));
% % %     tmpThrVec       = b2ScaleDimSd(1:3,iD,iV,iTyp,iReg,iRun,iM) ./ sqrt(sum(b2ScaleDimSd(1:3,iD,iV,iTyp,iReg,iRun,iM).^2,1));
% % %     projOntoGoal    = sum(tmpGlVec .* pointToPointVecs , 1);
% % % %     projOntoGoal    = unique(abs(projOntoGoal),'stable');
% % %     projOntoThreat  = sum(tmpThrVec .* pointToPointVecs , 1);
% % % %     projOntoThreat  = unique(abs(projOntoThreat),'stable');
% % % 
% % %     [rrr(iM,iRun,iTyp,iReg) ppp(iM,iRun,iTyp,iReg)] = corr(projOntoGoal',projOntoThreat');
% % % 
% % % end
% % % end
% % % end
% % % end
% % % 
% % % % $$$ Tomorrow, double check that the dot product isn't biased to being
% % % % smaller than 0.5 --> it is. 
% % % 
% % % % $$$ WAIT A MINUTE, did I just create a random cloud of points and then use that to show that two vectors are more or less orthogonal?...
% % % 
% % % % $$$ What do I do about it?
% % % 
% % % % $$$ --> let's see what the shitty iTyps show. Cause there shouldn't be
% % % % much structure there



%% $$$ Use PROPER permutation method And/or bootstrapping to perform stats for the PCA?

load(['Results\ForFigures\DimensionReduction\NetworkActivations\' ...
          'FullNetActivity_NeurType_' neurTypes{1} '_RegType_' regTypes{iReg} '_Run_' num2str(1) '_NetArch_' num2str(1) '.mat']);

% Select the particalar comparisons and PCA dimensions of interest
iV = 1;
iD = 1;

% Create 'explanatory variables'
thrRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[1 2 3 7]) ] );
thrRow = permute(thrRow, [2 3 4 1 5]);

golRow = repmat( (1:s.wrld.size(1))' , [1, size(glByThrAct,[2 3 6 7]) ] );
golRow = permute(golRow, [1 2 3 4 5]);

tic
for iTyp = 1:length(neurTypes)


for iReg = 1:length(regTypes)


% % % bFold = 'Results\ForFigures\Valence\';
bFold = 'Results\ForFigures\Valence\';
cFiles = dir([bFold '*' regTypes{iReg} '*_' neurTypes{iTyp} '_B_V*.mat']);


for iRun = 1:length(cFiles)

for iM = 1:size(rSall,1)

load(['Results\ForFigures\DimensionReduction\NetworkActivations\' ...
          'FullNetActivity_NeurType_' neurTypes{iTyp} '_RegType_' regTypes{iReg} '_Run_' num2str(iRun) '_NetArch_' num2str(iM) '.mat']);


  % Create flattened activations for PCA
    glByThrActFlat = permute(glByThrAct,[1 2 3 6 7 4 5]);
    unwrapSize     = size(glByThrActFlat);
    glByThrActFlat = permute(glByThrActFlat(:,:,:,:,:,:),[6 1 2 3 4 5]);
    glByThrActFlat = glByThrActFlat(:,:);

    

    % Perform PCA
    [coeff,score,~,~,explained] = pca(glByThrActFlat( ~isnan(glByThrActFlat(:,1)),:)');
    % Replace with zero if the PCA doesn't work
    if size(score,2) == 0
        score = zeros([size(glByThrActFlat,2) 3]);
    end


    % Downsample more here than when creating the basis vectors, because
    % otherwise we end up with a stupidly large array of vectors between
    % states
    tmpScore = reshape(score', [size(score,2) unwrapSize(1:5)] );
    tmpScore = tmpScore(:,:,2:2:end,1:2:end,:,2:2:end);
    tmpScore = tmpScore(:,:)';



    % Extract attributes of interest for given state
    allVars     = {'thrRow','thrCol','lmbCol','lmbCol','golDist',...
               'minStimDist','minStimDist'};
    allVarsComp = {'golRow','golCol','thrCol','golCol','thrDist',...
               'maxStimDist','avStimDist'};
    eval(['tmpCl    = ' allVars{iV} '(:,2:2:end,1:2:end,:,2:2:end);']);
    eval(['tmpSz    = ' allVarsComp{iV} '(:,2:2:end,1:2:end,:,2:2:end);']);
    tmpCl    = tmpCl(:);
    tmpSz    = tmpSz(:);

    
    nPm = 1000;
    dS = 1; dE = 3;
    clear dotProdPm

    % Perform permutations
    for iPm = 1:nPm
        permOrd = randperm(size(tmpCl,1));

        % Find vector that best explains the variable of interest A
        [b1Tmp, ~, ~, ~, stats] = regress(tmpCl(permOrd), [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);

        % Find vector that best explains the variable of interest B        
        [b2Tmp, ~, ~, ~, stats] = regress(tmpSz(permOrd), [tmpScore(:,dS:dE), ones(size(tmpScore, 1), 1)]);

        % Calculate dot and (non-normalised) cross-products
        dotProdPm(iPm)  = sum(b1Tmp(1:3)./norm(b1Tmp(1:3)) .* b2Tmp(1:3)./norm(b2Tmp(1:3)));
        crossProdPm(iPm) = norm(cross(b1Tmp(1:3) , b2Tmp(1:3)));
    end

    % Perform bootstrapping
    for iPm = 1:nPm
        permOrd = randi(size(tmpCl,1),size(tmpCl));

        % Find vector that best explains the variable of interest A
        [b1Tmp, ~, ~, ~, stats] = regress(tmpCl(permOrd), [tmpScore(permOrd,dS:dE), ones(size(tmpScore, 1), 1)]);
              % Find vector that best explains the variable of interest B
        [b2Tmp, ~, ~, ~, stats] = regress(tmpSz(permOrd), [tmpScore(permOrd,dS:dE), ones(size(tmpScore, 1), 1)]);

        % Calculate dot and (non-normalised) cross-products
        dotProdBs(iPm)  = sum(b1Tmp(1:3)./norm(b1Tmp(1:3)) .* b2Tmp(1:3)./norm(b2Tmp(1:3)));
        crossProdBs(iPm) = norm(cross(b1Tmp(1:3) , b2Tmp(1:3)));
    end
    
    
toc



    
    % Calculate dot and (non-normalised) cross-products
    dotProdBase = sum( b1(1:3,iD,iV,iTyp,iReg,iRun,iM)./norm(b1(1:3,iD,iV,iTyp,iReg,iRun,iM)) .* ...
        b2(1:3,iD,iV,iTyp,iReg,iRun,iM)./norm(b2(1:3,iD,iV,iTyp,iReg,iRun,iM)) );
    
    crossProdBase = norm(cross( b1(1:3,iD,iV,iTyp,iReg,iRun,iM) , ...
        b2(1:3,iD,iV,iTyp,iReg,iRun,iM) ));
    
    alldotProd(iM,iRun,iTyp,iReg,iD,iV)   = dotProdBase;
    allcrossProd(iM,iRun,iTyp,iReg,iD,iV) = crossProdBase;
    
    pDotProd(iM,iRun,iTyp,iReg,iD,iV)    = sum( abs(dotProdBase)   - abs(dotProdPm) > 0) ./ nPm;
    pCrossProd(iM,iRun,iTyp,iReg,iD,iV)  = sum( abs(crossProdBase) - abs(crossProdPm) < 0) ./ nPm;
    
    alldotProdPm(iM,iRun,iTyp,iReg,iD,iV,:)   = dotProdPm;
    allcrossProdPm(iM,iRun,iTyp,iReg,iD,iV,:) = crossProdPm;
    
    alldotProdBs(iM,iRun,iTyp,iReg,iD,iV,:)   = dotProdBs;
    allcrossProdBs(iM,iRun,iTyp,iReg,iD,iV,:) = crossProdBs;

        
toc
end
end
end
end
toc

% save('C:\Users\Rory Bufacchi\Documents\Projects\DefenseAgent\Results\ForFigures\NeurTypesAnalysed_V5.mat','-v7.3')


%% Finalise pca stats and plot scatter plus histograms

% $$$ HERE
% $$$ Next plot PCA and structure stuff for all the little networks here

% First demonstrate PCA consistency


% % % pValBootstrapPCA = sum(abs(alldotProdBs) > abs(cos(deg2rad(75))),7) ./ nPm;
iD = 1; iV = 1;

for iTyp = 1:length(neurTypes)
for iReg = 1:length(regTypes)
for iRun = 1:length(cFiles)
for iM = 1:size(rSall,1)

    try
    [pTmp,hTmp,stats] = signrank( ...
        abs(squeeze(alldotProdBs(iM,iRun,iTyp,iReg,iD,iV,:))),  ...
    abs(cos(deg2rad(75))) ...
    ,'method','approximate','tail','left');
    pValBootstrapPCA(iM,iRun,iTyp,iReg,iD,iV) = pTmp;
    signRankPCA(iM,iRun,iTyp,iReg,iD,iV) = stats.signedrank;
    zPCA(iM,iRun,iTyp,iReg,iD,iV) = stats.zval;
    catch
        pValBootstrapPCA(iM,iRun,iTyp,iReg,iD,iV) = 1;
        signRankPCA(iM,iRun,iTyp,iReg,iD,iV) = 0;
        zPCA(iM,iRun,iTyp,iReg,iD,iV) = 0;
    end

end
end
end
end

disp('effect size')
effSizePCA    = zPCA ./ sqrt(nPm)
effSizePCAdot = abs(alldotProd) -  abs(cos(deg2rad(75))) ./ ...
    std(...
    abs(alldotProdBs) - abs(cos(deg2rad(75))) ...
    ,0,7)
effSizePCAcross = abs(allcrossProd) ./ ...
    std(abs(allcrossProdBs),0,7)

disp('mean Z')
nanmean(zPCA(:),7)

tmpP = pValBootstrapPCA(:,1:15);
tmpCos = alldotProd(:);
tmpSin = sin(acos(alldotProd(:)));

[pValThr, pValcor, pValPCAAdj] = fdr(tmpP(:));

fprintf('%i out of %i networks show structure using this metric\n', sum(pValPCAAdj(:) < 0.05), numel(pValPCAAdj));

disp('mean rho')
rad2deg(atan2(mean(tmpSin(~isnan(tmpSin))) , mean(tmpCos(~isnan(tmpSin)))))
disp('std rho')
rad2deg(atan2(std(tmpSin(~isnan(tmpSin))) , std(tmpCos(~isnan(tmpSin)))))


tStatsStructureTmp = tStatsStructure(:,1:15,:,:,:);
allPerfTmp         = allPerf(:,1:15,:,:,:);
% % % tStatsStructureTmp = tStatsStructure(:,1:15,1,1,1);
% % % allPerfTmp         = allPerf(:,1:15,1,1,1);
disp('correlation between structure t-stat and network performance')
[rhorho pvalpval] = corr(tStatsStructureTmp(:),allPerfTmp(:))



% % % metDefs = {'estimateDistEff','alldotProd','estimateDistEff'};
% % % yDefs   = {'allPerf'        ,'allPerf'   ,'alldotProd'};

absalldotProd = abs(alldotProd);
% metDefs = {'estimateDistEff','alldotProd'};
% yDefs   = {'allPerf'        ,'allPerf'   };
% metDefs = {'tStatsStructure','absalldotProd'};
% yDefs   = {'allPerf'        ,'allPerf'   };
metDefs = {'tStatsStructure'};
yDefs   = {'allPerf'        };
% yDefs   = {'rhoNewTasks'        };





%% Metric count

% f.StructureVsReward.f        = figure('Position',[50,200,1200,500]); 


rowDefs = neurTypes;
for iReg =1:length(regTypes)
    f.(['StructureVsReward' regTypes{iReg}]).f        = figure('Position',[50,20,300,1800]); 
for iRow = 1:numel(rowDefs) 
for iMet = 1:length(metDefs)

    % Now scatter plot with side histograms
    eval([' tmpX = ' metDefs{iMet} '(:,1:15,' num2str(iRow) ',' num2str(iReg) ');']);
    eval([' tmpY = ' yDefs{iMet}   '(:,1:15,' num2str(iRow) ',' num2str(iReg) ');']);
    tmpX = tmpX(:);
    tmpY = tmpY(:);


    axS.xOffset = 0.05;
    axS.yOffset = 0.05;

    axS.xSpace  = 0.4./length(metDefs);
    axS.ySpace  = 0.4./length(rowDefs);


    ax{1} = axes('Position',[axS.xOffset + ...
        axS.xSpace  + (iMet-1) .* ((1 - 2.*axS.xOffset) ./ length(metDefs)) , ...
        axS.yOffset + (iRow-1) .* ((1 - 2.*axS.yOffset) ./ length(rowDefs)), ...
        (  ((1 - 2.*axS.xOffset) ./ length(metDefs)) - axS.xSpace ) , ...
        (  ((1 - 2.*axS.yOffset) ./ length(rowDefs)) - axS.ySpace ) ]);

%     plot(tmpX,tmpY,'o')
    plot(tmpX,tmpY,'.')
    set(ax{1},'Visible','off')
    yLims = [-.1 .1];
    ylim(yLims)

    xLims = [-20 20];
%     xLims = xlim;
%     yLims = ylim;
    if iMet == 1
        xlim(max(abs(xLims)) .* [-1 1] );
        %     ylim(max(abs(yLims)) .* [-1 1] );
        xLims = xlim;
        %     yLims = ylim;
    end
    lsline
%     % Put correlation in the title
%     title(CorrelationTitle(tmpX,tmpY));


    

    ax{2} = axes('Position',ax{1}.Position .* [1 1 1 0] + [0 -axS.ySpace.*.9 0 axS.ySpace.*.9] + [0 axS.ySpace+ax{1}.Position(4) 0 0] );
    h = histogram(tmpX,linspace(xLims(1),xLims(2),11)); hold on
    % set(ax{2},'ydir','r')
    box off
    xlim(xLims);
    [splineHandle] = PlotHistSpline(h,'LineWidth',2,'Color','b')
%     title(metDefs{iMet});
    % Put correlation in the title
    title(CorrelationTitle(tmpX,tmpY));

    ax{3} = axes('Position',ax{1}.Position .* [1 1 0 1] + [-axS.xSpace.*.9  0 axS.xSpace.*.9  0]  );
    h = histogram(tmpY,linspace(yLims(1),yLims(2),11),'orientation','horizontal'); hold on
    [splineHandle] = PlotHistSpline(h,'LineWidth',2,'Color','b')
    set(ax{3},'xdir','r')
    box off
    ylim(yLims);
    title([yDefs{iMet} ' ' rowDefs{iRow}]);

end
end
sgtitle(regTypes{iReg})


f.(['PCAangles' regTypes{iReg}]).f        = figure('Position',[50,20,400,1800]); 
for iRow = 1:numel(rowDefs) 

% Also plot angle distributions
allAnglePm = acos(alldotProdPm);
allAngleBs = acos(alldotProdBs);


subplot(numel(rowDefs) ,2,(iRow -1) .* 2 + 1)
allAngToPl = squeeze(allAngleBs(:,1:15,iRow,iReg,1,1,:));
h = polarhistogram(real(allAngToPl),13,'Normalization','probability'); hold on
[splineHandle] = PlotHistSpline(h,'LineWidth',2,'Color','b')
rLims = rlim;
pAx = gca;
pAx.ThetaLim = [0 180];
% title('angle between goal and threat')
title('pValues?')

subplot(numel(rowDefs) ,2,(iRow -1) .* 2 + 2)
allAngToPl = squeeze(allAnglePm(:,1:15,iRow,iReg,1,1,:));
h = polarhistogram(real(allAngToPl),13,'Normalization','probability'); hold on
[splineHandle] = PlotHistSpline(h,'LineWidth',2,'Color','b')
title('permuted angle between goal and threat')
rlim(rLims);
pAx = gca;
pAx.ThetaLim = [0 180];


end
sgtitle(regTypes{iReg})

end

%%

cReg = 1:2;

% Plot performance vs orthogonality and structure
metDefs = {'tStatsStructure'};
yDefs   = {'allPerf'        };

% $$$ FOr structure, check out line 128

% figure('Position',[50,500,1800,400])
f.PCAvsRew.f            = figure('Position',[50,200,600,600]);
f.StructvsRew.f         = figure('Position',[50,200,600,600]);
f.PCAvsStructvsRew.f    = figure('Position',[50,200,600,600]);
f.StructvsReconsstruct.f= figure('Position',[50,200,600,600]);
f.PCAvsReconsstruct.f   = figure('Position',[50,200,600,600]);

cM = 1:3;
cReg = 1:2;

iPl = 1;

figure(f.PCAvsRew.f)
% PCA ORTHOGONALITY -------------------------------------------------------
tmpX = absalldotProd(cM,:,:,cReg);
tmpY = allPerf(cM,:,:,cReg);
scatter(tmpX(~isinf(tmpX(:)) & ~isnan(tmpX(:))), tmpY(~isinf(tmpX(:)) & ~isnan(tmpX(:))),'filled'); hold on
lsline
tmpX = permute(absalldotProd(cM,:,:,cReg),[3 1 2 4]); tmpX = tmpX(:,:);
tmpY = permute(allPerf(cM,:,:,cReg),[3 1 2 4]); tmpY = tmpY(:,:);
scatter(tmpX',tmpY','filled'); hold on
xlabel('AbsdDot product (paralellity)')
ylabel('Performance(reward/timestep)')
% [xQuery, yAvg] = SlidingWindowAverage(tmpX(:), tmpY(:), linspace(min(tmpX(:)),max(tmpX(:)),100), min(tmpX(:)-))
[xQuery, yAvg] = SlidingWindowAverage(tmpX(:), tmpY(:));
[xQuery, yStd] = SlidingWindowFunction(tmpX(:), tmpY(:),@(x)std(x)./sqrt(numel(x)) );
shOpts.PlotMean = 0; shOpts.c = [0.7 0.7 0.7];
plot(xQuery,yAvg,'k','LineWidth',2); hold on
ShadedPlot(xQuery,yAvg, yAvg - yStd, yAvg + yStd, shOpts)



% NETWORK STRUCTURE -------------------------------------------------------
figure(f.StructvsRew.f)
tmpX =  tStatsStructure(cM,:,:,cReg);
tmpY = allPerf(cM,:,:,cReg);
scatter(tmpX(:), tmpY(:),'filled'); hold on
lsline
tmpX = permute(tStatsStructure(cM,:,:,cReg),[3 1 2 4]); tmpX = tmpX(:,:);
tmpY = permute(allPerf(cM,:,:,cReg),[3 1 2 4]); tmpY = tmpY(:,:);
scatter(tmpX',tmpY','filled'); hold on
% ylim([-0.1 0.1])
% xlim([-20 20])
[xQuery, yAvg] = SlidingWindowAverage(tmpX(:), tmpY(:));
[xQuery, yStd] = SlidingWindowFunction(tmpX(:), tmpY(:),@(x)std(x)./sqrt(numel(x)) );
plot(xQuery,yAvg,'k','LineWidth',2); hold on
ShadedPlot(xQuery,yAvg, yAvg - yStd, yAvg + yStd, shOpts)
xlabel('Network structure (t-stat)')
ylabel('Performance(reward/timestep)')



figure(f.StructvsReconsstruct.f)
% PCA VS RECONSIRUCTION ---------------------------------------------
tmpX = tStatsStructureAll (cM,:,:,cReg,:);
tmpY = squeeze(rhoAll(2,cM,:,:,cReg,:));
% % % h = histogram2(tmpX(~isinf(tmpX(:)) & ~isnan(tmpX(:))), tmpY(~isinf(tmpX(:)) & ~isnan(tmpX(:))),'FaceColor','Flat');
% % % densCol = h.Values ./ sum(h.Values,2);
% % % [XX YY] = meshgrid(h.XBinEdges(1:end-1) + diff(h.XBinEdges(1:2)), h.YBinEdges(1:end-1) + diff(h.YBinEdges(1:2)) )
% % % imagesc(h.XBinEdges(1:end-1) + diff(h.XBinEdges(1:2)) ./ 2, ...
% % %         h.YBinEdges(1:end-1) + diff(h.YBinEdges(1:2)) ./ 2,densCol'); axis xy
% % % colormap(whitetocol(256,[0 0 0]));
% % % hold on
scatter(tmpX(~isinf(tmpX(:)) & ~isnan(tmpX(:))), tmpY(~isinf(tmpX(:)) & ~isnan(tmpX(:))),'.'); hold on
lsline

tmpX = permute(tStatsStructureAll(cM,:,:,cReg,:),[3 1 2 4 5]); tmpX = tmpX(:,:);
tmpY = permute(squeeze(rhoAll(2,cM,:,:,cReg,:)),[3 1 2 4 5]); tmpY = tmpY(:,:);
% % % scatter(tmpX',tmpY','.'); hold on
scatter(tmpX',tmpY',10,'Filled'); hold on
xlabel('Network structure (t-stat)')
ylabel('Reconstruction quality (rho)')
% [xQuery, yAvg] = SlidingWindowAverage(tmpX(:), tmpY(:), linspace(min(tmpX(:)),max(tmpX(:)),100), min(tmpX(:)-))
[xQuery, yAvg] = SlidingWindowAverage(tmpX(:), tmpY(:));
[xQuery, yStd] = SlidingWindowFunction(tmpX(:), tmpY(:),@(x)std(x)./sqrt(numel(x)) );
shOpts.PlotMean = 0; shOpts.c = [0.7 0.7 0.7];
plot(xQuery,yAvg,'k','LineWidth',2); hold on
ShadedPlot(xQuery,yAvg, yAvg - yStd, yAvg + yStd, shOpts)
ylim([0 1]);
xlim([-30 10]);




figure(f.PCAvsReconsstruct.f)
% PCA VS RECONSIRUCTION ---------------------------------------------
tmpX = absalldotProdAll(cM,:,:,cReg,:);
tmpY = squeeze(rhoAll(2,cM,:,:,cReg,:));
% % % h = histogram2(tmpX(~isinf(tmpX(:)) & ~isnan(tmpX(:))), tmpY(~isinf(tmpX(:)) & ~isnan(tmpX(:))),'FaceColor','Flat');
% % % densCol = h.Values ; %./ sum(h.Values,2);
% % % [XX YY] = meshgrid(h.XBinEdges(1:end-1) + diff(h.XBinEdges(1:2)), h.YBinEdges(1:end-1) + diff(h.YBinEdges(1:2)) )
% % % imagesc(h.XBinEdges(1:end-1) + diff(h.XBinEdges(1:2)) ./ 2, ...
% % %         h.YBinEdges(1:end-1) + diff(h.YBinEdges(1:2)) ./ 2,densCol'); axis xy
% % % colormap(whitetocol(256,[0 0 0]));
% % % hold on
scatter(tmpX(~isinf(tmpX(:)) & ~isnan(tmpX(:))), tmpY(~isinf(tmpX(:)) & ~isnan(tmpX(:))),'.'); hold on
lsline

tmpX = permute(absalldotProdAll(cM,:,:,cReg,:),[3 1 2 4 5]); tmpX = tmpX(:,:);
tmpY = permute(squeeze(rhoAll(2,cM,:,:,cReg,:)),[3 1 2 4 5]); tmpY = tmpY(:,:);
% % % scatter(tmpX',tmpY','.'); hold on
scatter(tmpX',tmpY',10,'Filled'); hold on
xlabel('AbsdDot product (paralellity)')
ylabel('Reconstruction quality (rho)')
% [xQuery, yAvg] = SlidingWindowAverage(tmpX(:), tmpY(:), linspace(min(tmpX(:)),max(tmpX(:)),100), min(tmpX(:)-))
[xQuery, yAvg] = SlidingWindowAverage(tmpX(:), tmpY(:));
[xQuery, yStd] = SlidingWindowFunction(tmpX(:), tmpY(:),@(x)std(x)./sqrt(numel(x)) );
shOpts.PlotMean = 0; shOpts.c = [0.7 0.7 0.7];
plot(xQuery,yAvg,'k','LineWidth',2); hold on
ShadedPlot(xQuery,yAvg, yAvg - yStd, yAvg + yStd, shOpts)




%%

cReg = 1:2;

for iTyp = 1:length(neurTypes)

%     subplot(1,length(neurTypes),iPl)
%     iPl = iPl + 1;

figure(f.PCAvsRew.f)
% PCA ORTHOGONALITY -------------------------------------------------------
tmpX = log(squeeze(dotProdAbs(1,1,iTyp,cReg,:,cM)));
% tmpX = squeeze(crossProdSize(1,1,iTyp,cReg,:,cM));
% tmpX = squeeze(crossProdRelMax(1,1,iTyp,cReg,:,cM));

tmpY = rewPerAct(iTyp,cReg,:,cM);

scatter(tmpX(:), tmpY(:),'filled'); hold on

% lsline
ylim([-0.1 0.1])
% xlim([0 10])
% xlim([-6 2])
xlabel('Dot product (paralellity)')
ylabel('Reward per action')
sgtitle(['CM ' num2str(cM) ])


figure(f.StructvsRew.f)
% NETWORK STRUCTURE -------------------------------------------------------
% dotprob: (iD,iV,iTyp,iReg,iRun,iM)
% allbinNPF: (dim1, dim2, iM,iRun,iTyp,iReg)
tmpX = permute(squeeze(sum(abs(allBinNPF(:,:,cM,:,iTyp,cReg)),[1 2])), [3 2 1]);
% % % % Subtract the permuted average
% % % if numel(cM) > 1
% % %     tmpX = tmpX - permute(squeeze(mean(sum(abs(allBinNPFpm(:,:,:,:,:,iTyp,:)),[1 2]),3)), [3 2 1]);
% % % else
% % %     tmpX = tmpX - permute(squeeze(mean(sum(abs(allBinNPFpm(:,:,:,cM,:,iTyp,:)),[1 2]),[3 5])), [2 1 3]);
% % % end

tmpY = rewPerAct(iTyp,cReg,:,cM);

scatter(tmpX(:), tmpY(:),'filled'); hold on
% lsline
ylim([-0.1 0.1])
xlim([0 20])
xlabel('Network structure')
ylabel('Reward per action')
legend(['all','linear fit','smoothed fit',neurTypes])


figure(f.PCAvsStructvsRew.f)
% BOTH -------------------------------------------------------
tmpX = log(squeeze(dotProdAbs(1,1,iTyp,cReg,:,cM)));
tmpY = permute(squeeze(sum(abs(allBinNPF(:,:,cM,:,iTyp,cReg)),[1 2])), [3 2 1]);
tmpZ = rewPerAct(iTyp,cReg,:,cM);

scatter3(tmpX(:), tmpY(:), tmpZ(:), 'filled'); hold on
% lsline
xlabel('Dot product (paralellity)')
ylabel('Network structure')
zlabel('Reward per action')


end




%  $$$ dotProdAbs(1,1,iTyp,iReg,:,iM)

%% Plot icons for neuron types

% Input strength range
x = -5:0.1:5;

% Initialize a figure
f.NeurTypesSketch.f         = figure('Position',[20 -20 900 900]);

for i = 1:length(neurTypes)
    % Compute neuron output for given input strength
    switch neurTypes{i}
        case 'tansig'
            y = tansig(x);
        case 'logsig'
            y = logsig(x);
        case 'softmax'
            % Softmax is a bit unique since it's typically used for multi-dimensional input
            % For our example, let's consider a 2-D input where the second dimension is just negative of the first.
            y = softmax([x; -x]);
            y = y(1, :); % consider one dimension for display
        case 'poslin'
            y = poslin(x);
        case 'purelin'
            y = purelin(x);
        case 'tribas'
            y = tribas(x);
        case 'radbas'
            y = radbas(x);
        otherwise
            y = zeros(size(x));
    end
    
    % Create a subplot and plot the neuron's response
    subplot(4, 2, i);
    plot(x, y, 'LineWidth', 1.5);
    title(neurTypes{i});
    xlabel('Input Strength');
    ylabel('Output');
    grid on;
end

% Adjust layout
tight_layout = true;
if tight_layout
    sgtitle('Response Functions of Various Neuron Types');
end



%% Save figures for neuron types
allFields = fields(f);
for iF = 1:length(allFields)
    cF = allFields{iF};
  
    set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['Results\ForFigures\Revision\NeurTypes\' cF '.eps'] , 'epsc')
    saveas(f.(cF).f,['Results\ForFigures\Revision\NeurTypes\' cF '.pdf'] , 'pdf')
% % %     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Revision\NeurTypes\CartDist' cF '.eps'] , 'epsc')
% % %     saveas(f.(cF).f,['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Revision\NeurTypes\CartDist' cF '.pdf'] , 'pdf')

end


%% Functions



function [ax, plQhoriz, plQvert] = plot2D(rS,fS,Q,iM,iD);

    
    rS(iM,iD).s.plt.contours = 1;
    rS(iM,iD).s.plt.contourVals = [-0.25, 0.25];
    
    rS(iM,iD).s.plt.plAct = 1;

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


function [] = MakeSidePlots(baseAx,Q,fS,rS)

    Q = - Q;

    ax{1} = baseAx;

    % Plot average OR line through Q-values
    ax{2} = axes('Position',ax{1}.Position .* [1 1 1 0] + [0 -0.1 0 0.1]  );
% % %     plQhoriz = flip(squeeze(nanmean(Q(1,rS.s.plt.lmbCol,2:11,2:end-1,:),[2 3 5]))); % lmb R, C, stim R, C, A
    plQhoriz = flip(squeeze(nanmean(Q(1,rS.s.plt.lmbCol,:,3:end-2,:),[2 3 5]))); % lmb R, C, stim R, C, A
    plot((1:numel(plQhoriz))' - 0.5  , [plQhoriz],'-','LineWidth',2); hold on
    xLims = xlim;
    xlim([0 numel(plQhoriz) ]);
    hold on
    ylim(fS.lineYLim);
    grid on
    ax{2}.XAxis.Visible = 'off';
    ax{2}.YAxis.Visible = 'off';
    ax{3} = axes('Position',ax{1}.Position .* [1 1 0 1] + [-0.03 0 0.03 0]  );
% % %     plQvert = squeeze(nanmean(Q(1,rS.s.plt.lmbCol,2:end-1,2:end-1,:),[2 4 5])); % lmb R, C, stim R, C, A
    plQvert = squeeze(nanmean(Q(1,rS.s.plt.lmbCol,2:end-1,:,:),[2 4 5])); % lmb R, C, stim R, C, A
    plot( plQvert, (numel(plQvert):-1:1)' - 0.5 , '-','LineWidth',2); hold on
    xlim(fS.lineXLim);
    ylim([0 12])
    grid on
    ax{3}.XAxis.Visible = 'off';
    ax{3}.YAxis.Visible = 'off';

end


