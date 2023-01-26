%% Run ToolUse or load this:

load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ToolUse\WorkSpace_ForFigures.mat');

addpath('D:\Old_D\DPPS\DefenseAgent\Scripts\rlsimplepps')
addpath('D:\Old_D\DPPS\DefenseAgent\Scripts')
addpath('D:\Old_D\Programs\Matlab\Utilities\')
addpath('D:\Old_D\Programs\Matlab\Utilities\plotting\')
addpath('D:\Old_D\Programs\Matlab\Utilities\plotting\colormaps\')

%% Plot heatmaps for low valence, high valence and for negative positive valence?

f.ToolMap.f = figure('Position',[20 20 800 800])

cM = 1;

sFP = DefaultSettings(olRS(cM).s);
net = olRS(cM).net;
w   = olRS(cM).w;

sFP.plt.lmbCol         = 2:14;
sFP.plt.stimRow        = [1:size(w.world2D,1)];
sFP.plt.stimCol        = [2:size(w.world2D,1)];
sFP.plt.rowLims        = [1.5 13.5];
sFP.plt.meanLimbCols   = 1;
sFP.plt.plAct          = 2;



fS.gridXstart = -4.5;
fS.gridXstep = 1;
fS.gridYstart = 3.5;
fS.gridYstep = 1;


allRowLags = olRS(cM).s.gol.alSpR;

condN = {'No Tool, Pre Training' , 'Tool, Pre Training', ...
         'No Tool, Post Training', 'Tool, post Training'};

for iCond = 1:4

        ax{iCond} = subplot(2,2,iCond);
                
%         sFP.plt.intrpFact = .25;
        sFP.plt.intrpFact = 1;
        
        
        % $$$ HERE, need to change this to take q values from the yes and
        % no tool section - so gotta create those first of course :)
        
        % Plot average over all models
        %                   Qall(:,:,:,:,:,iCL,iRL,cM)
        qTemp = nanmean(Qall(:,:,:,:,:,:,iCond),[6]);
        
        DisplActValsFun(sFP,w,-qTemp)
        GridOverImage(fS,ax{iCond});
        
        box off
        
        caxis([-2 2])
%         caxis([0 2])
        
        colorbar off
        
        title([condN{iCond}])

end

colormap(redbluecmapRory(5,5))
% colormap(whitetocol(100,[0 0 0.7]));
% $$$ NOW THINK WHAT TO DO WITH THE INDIVIDUAL NEURONS?
% $$$ MAYBE proportion of neurons with multiple peaks?

% Bar plot showing number of peaks?
for iCond = 1:4
    
    axb = axes('Position',[ [ax{iCond}.Position(1:2) + [ax{iCond}.Position(3).*.1 ax{iCond}.Position(3).*.5] ] [ ax{iCond}.Position(3:4).*[0.15 .4]] ] ,'Box','off');

    plot([0 2; 0 2]', [1 1 ; 2 2]' ,'-.k','LineWidth',0.5); hold on
    
    tmpM  = nanmean(npksVec(:,iCond));
    tmpSD = nanstd(npksVec(:,iCond));
    bar(tmpM);
    hold on
    plot([1 1],tmpM + [ - tmpSD, tmpSD],'k','LineWidth',2);
        
    ylim([0 3])
    xlim([0 2])
    
    title('n peaks')
    
    set(axb,'xTickLabels',[])
    yticks([0 1 2])
    
end

sgtitle('Value fields expand dynamically with tool use, but only after training')

% % % %% Plot 'extent' barplots for low/high valence. AND POSSIBLY for neg pos valence
% % % 
% % % figure,
% % % 
% % % tmpM    = 11 - nanmean(allOverThresh(:,:),2);
% % % tmpSD   = nanstd(allOverThresh(:,:),[],2);
% % % plot([0 6; 0 6; 0 6]', [1 1 ; 2 2; 3 3]' ,'-.k','LineWidth',0.5); hold on
% % % bar( [2 4], tmpM ); hold on
% % % plot([2 2 ; 4 4]',[tmpM + [ - tmpSD'; tmpSD']']','k','LineWidth',2);
% % % 
% % % xlim([.5 5.5])
% % % 
% % % set(gca,'xTick',[2 4])
% % % 
% % % xlabel('Stim valence')
% % % ylabel('Field size')


% % % %% Plot number of peaks as a function of layer depth
% % % 
% % % 
% % % figure,
% % % 
% % % 
% % % for iCond = 1:4
% % % 
% % %     ax{iCond} = subplot(2,2,iCond);
% % %     
% % %     plot([0 6; 0 6]', [1 1 ; 2 2]' ,'-.k','LineWidth',0.5); hold on
% % %                 
% % %     % npksN(rows,iL,iN,iM,iCond) 
% % %     tmpM    = nanmean(npksN(:,:,:,:,iCond),[1,3,4]);
% % %     tmpD    = nanmean(npksN(:,:,:,:,iCond),[1]);
% % %     tmpSD   = nanstd(tmpD(:,:,:,:),[],[3,4]);
% % %     bar(tmpM) % one value for each of the 5 layers
% % %     
% % %     
% % %     arrayfun(@(n) plot([n n],tmpM(n) + [ - tmpSD(n), tmpSD(n)],'k','LineWidth',2) ,1:length(tmpM));
% % %     
% % %     title(['mean #peaks in neurons when ' condN{iCond}])
% % %     
% % %     ylim([0 3])
% % %     xlim([0 6])
% % %     
% % %     xlabel('layer')
% % %     ylabel('#peaks')
% % % 
% % % end
% % % % colormap(whitetocol(100,[0 0 0.7]));
% % % colormap(redbluecmapRory(5,5))
% % % 
% % % 
% % % 
% % % %% Plot proportion of neurons that increase their number of peaks as a function of layer depth
% % % 
% % % 
% % % figure,
% % % 
% % % subplot(1,2,1)
% % % fS.nBins = 3;
% % % fS.pltType = 'bar';
% % % plot([0 1]', [.5 .5]' ,'-.k','LineWidth',0.5); hold on
% % % BinPlot(rMatYT.relLayDep,rMat.propOfCorrToolRespPreTrain,fS,'LineWidth',2)
% % % xlabel('binned layer depth')
% % % ylabel('binned proportion of field neurons')
% % % ylim([0 1])
% % % xlim([0 1])
% % % title('PRE-TRAINING')
% % % 
% % % 
% % % subplot(1,2,2)
% % % plot([0 1]', [.5 .5]' ,'-.k','LineWidth',0.5); hold on
% % % BinPlot(rMatYT.relLayDep,rMat.propOfCorrToolRespPostTrain,fS,'LineWidth',2)
% % % xlabel('binned layer depth')
% % % ylabel('binned proportion of field neurons')
% % % ylim([0 1])
% % % xlim([0 1])
% % % title('POST-TRAINING')
% % % 
% % % 
% % % sgtitle('proportion of field neurons which expand field in response to tool')


%% Plot proportion of neurons that increase their number of peaks after training as a function of layer depth


f.ToolNeurs.f = figure('Position',[20 20 400 800]);

clear propNeurs
nSplits = size(rMat.propOfCorrToolRespPreTrain,1);

for iC = 1:1
for iSpl = 1:size(rMat.propOfCorrToolRespPreTrain,1)
    
    tmpD = rMat.propOfCorrToolRespPrePostTrain;
    
    for iM = 1:size(tmpD,2)
             
        propNeurs(iM,iSpl) = nanmean(tmpD(iSpl,iM),[2 1]);
        if propNeurs(iM,iSpl) == 0
            propNeurs(iM,iSpl) = .001
        end
    end
    
    subplot(nSplits,1,iSpl )
%     p = pie([nanmean(propNeurs(:,iSpl),1) ; 1 - nanmean(propNeurs(:,iSpl),1)],{'',''});
    p = pie([nanmean(propNeurs(:,iSpl),1) ; 1 - nanmean(propNeurs(:,iSpl),1)]);

    set(p,'EdgeColor','k')
    
    if iSpl == 1
        title('Proportion of Neurons with receptive fields that expand upon tool use');
    end
    
    
caxis([1 2]);
end
end



colormap(coltowhite(100,[0 0.2 .9]));





% %% Plot proportion of neurons that increase their number of peaks as a function of layer depth
% 
% clear propNeurs
% nSplits = size(rMat.propOfCorrToolRespPreTrain,1);
% figure('Position',[20 20 600 800]);
% for iC = 1:2
% for iSpl = 1:size(rMat.propOfCorrToolRespPreTrain,1)
%     
%     if iC == 1
%         tmpD = rMat.propOfCorrToolRespPreTrain;
%     else
%         tmpD = rMat.propOfCorrToolRespPostTrain;
%     end
%     
%     for iM = 1:size(tmpD,2)
%              
%         propNeurs(iM,iSpl) = nanmean(tmpD(iSpl,iM),[2 1]);
%         if propNeurs(iM,iSpl) == 0
%             propNeurs(iM,iSpl) = .001
%         end
%     end
%     
%     subplot(nSplits,2,iSpl.*2 - 2 + iC )
%     p = pie([nanmean(propNeurs(:,iSpl),1) ; 1 - nanmean(propNeurs(:,iSpl),1)],{'',''});
% 
%     set(p,'EdgeColor','k')
%     
%     if iC == 1 & iSpl == 1
%         title('Before Tool Training')
%     elseif iC == 2 & iSpl == 1
%         title('After Tool Training')
%     end
%     
% %     title([layNames{iSpl} ' layers'])
%     
% caxis([1 2]);
% end
% end
% 
% sgtitle('Proportion of Neurons with receptive fields that expand upon tool use')
% 
% colormap(coltowhite(100,[0 0.2 .9]));


%% Save figures


allFields = fields(f);
for iF = 1:length(allFields)
    cF = allFields{iF};
    
    set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ToolUse\BitsAndPieces\FromMatlab\' cF '_V3.eps'] , 'epsc')
    saveas(f.(cF).f,['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ToolUse\BitsAndPieces\FromMatlab\' cF '_V3.pdf'] , 'pdf')
end

