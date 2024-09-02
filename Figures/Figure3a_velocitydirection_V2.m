%% MAKE SURE TO RUN VELOCITYDIRECTION FIRST

% OR just load this
% load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\BitsAndPieces\VelDirWorkSpace.mat');
% load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\BitsAndPieces\VelDirWorkSpace_v2.mat');

%% Plot heatmaps

f.Speed.f = figure('Position',[20 20 1200 400]);
% subplot(2,2,1)

cM = 1;

sFP = DefaultSettings(rS(cM).s);
net = rS(cM).net;
w   = rS(cM).w;

sFP.plt.lmbCol         = 2:14;
sFP.plt.stimRow        = [1:size(w.world2D,1)];
sFP.plt.stimCol        = [2:size(w.world2D,1)];
sFP.plt.rowLims        = [1.5 13.5];
sFP.plt.meanLimbCols   = 1;
sFP.plt.plAct          = 2;


sFP.plt.rowLag = 1;
sFP.plt.colLag = 0;

sFP.plt.intrpFact = 1;

fS.gridXstart = -4.5;
fS.gridXstep = 1;
fS.gridYstart = 3.5;
fS.gridYstep = 1;

cDat = [0 0 1 ; .5 0 .5 ; .7 0 0];


allRowLags = rS(cM).s.gol.alSpR;
allColLags = rS(cM).s.gol.alSpC;

% -------------------------------------------------------------------------
% Plot velocity effect
plotRowLags = {[1] , [2], [3]};
plotColLags = {[-2:2] , [-2:2], [-2:2]};
nPlts       = length(plotRowLags);
plNames     = {'Slow','Medium','Fast'};
for iPl = 1:nPlts
    
    vAx{iPl} = subplot(1,nPlts,iPl);
    %         subplot(nPlts,1,iPl);
    
    % Plot average over all models
    qTemp = nanmean(QallOrig(:,:,:,:,:,...
        find(ismember(allColLags,plotColLags{iPl})),...
        find(ismember(allRowLags,plotRowLags{iPl})),  1: end-1 ),[6 7 8]);
    
    for iAct = 1:length(qTemp)
        qTemp(:,:,:,:,iAct) = nanmean(qTemp,5);
    end
    
    DisplActValsFun(sFP,w,qTemp); hold on
    caxis([0 2])
    colorbar off
    GridOverImage(fS,vAx{iPl});
    xLims = xlim;
    yLims = ylim;
    plot([0 0],yLims,'-.k')
    plot(xLims,[10.5 10.5],'-.k')
    
    vAx{iPl}.YAxis.Visible = 'off';
    vAx{iPl}.Box = 'off';
    
    title(plNames{iPl})
    
% % %         % % IN SEPARATE INSET FOR EACH VELOCITY
% % %         % 1 Dimensional extent plot
% % %         axb = axes('Position',[ [vAx{iPl}.Position(1:2) + [vAx{iPl}.Position(3).*.1 vAx{iPl}.Position(4).*.5 ] ] [ vAx{iPl}.Position(3:4).*[0.3 .4]] ] ,'Box','off');
% % %         tmpD = allOverThresh(:,:);
% % %         bar(12 - nanmean(tmpD(iPl,:),2)); hold on
% % %         plot([1 1]', 12 - [nanmean(tmpD(iPl,:),2) - nanstd(tmpD(iPl,:),[],2)./2  ...
% % %              nanmean(tmpD(iPl,:),2) + nanstd(tmpD(iPl,:),[],2)./2  ]',...
% % %              'k','LineWidth',2);
% % %         ylim([0 5])
% % %         axb.XAxis.Visible = 'off';
% % %         axb.Box = 'off';
% % %         title('Field Size')
    
    % 1 Dimensional extent plot
    if iPl == nPlts
        axb = axes('Position',[ [vAx{1}.Position(1:2) + [vAx{1}.Position(3).*.1 vAx{1}.Position(4).*.5 ] ] [ vAx{1}.Position(3:4).*[0.3 .4]] ] ,'Box','off');
        tmpD = allOverThresh(:,:);
        b = bar(12 - nanmean(tmpD(:,:),2),'FaceColor','flat'); hold on
        b.CData = cDat;
        for iiPl = 1:nPlts
        plot([iiPl iiPl]', 12 - [nanmean(tmpD(iiPl,:),2) - nanstd(tmpD(iiPl,:),[],2)./2  ...
            nanmean(tmpD(iiPl,:),2) + nanstd(tmpD(iiPl,:),[],2)./2  ]',...
            'k','LineWidth',2);

% % %         plot(iiPl .* ones([size(tmpD,2) 1]) + rand([size(tmpD,2) 1]).*0.4, 12 - tmpD(iiPl,:),'.');

        end
        boxplot(12-tmpD');
        plotRowAveragesWithSpread(12 - tmpD);
        ylim([0 8])
        axb.XAxis.Visible = 'off';
        axb.Box = 'off';
        title('Field Size')
        xlim([.5 3.5])
    end
    
   
end

for iPl = 1:nPlts
    colormap(vAx{iPl},whitetocol(100,cDat(iPl,:)));
end

sgtitle('Value fields expand in response to faster stimuli')

% -------------------------------------------------------------------------
% Plot direction effect

plotRowLags = {[1:3] , [1:3] , [1:3]};
plotColLags = { 2 0 -2 };
nPlts       = length(plotRowLags);
plNames     = {'Left','Center','Right'};
f.Dir.f = figure('Position',[20 20 1200 400]);
for iPl = 1:nPlts

%         subplot(nPlts,1,iPl);
        dAx{iPl} = subplot(1,nPlts,iPl);
        
        % Plot average over all models
        qTemp = nanmean(QallOrig(:,:,:,:,:,...
            find(ismember(allColLags,plotColLags{iPl})),...
            find(ismember(allRowLags,plotRowLags{iPl})),  1: end-1 ),[6 7 8]);
        
        for iAct = 1:length(qTemp)
            qTemp(:,:,:,:,iAct) = nanmean(qTemp,5);
        end
        
        DisplActValsFun(sFP,w,qTemp); hold on
        yLims = ylim;
        plot([fig0 0],yLims,'-.k')
        caxis([0 2])
        colorbar off
        GridOverImage(fS,vAx{iPl});
        xLims = xlim;
        yLims = ylim;
        plot([0 0],yLims,'-.k')
        plot(xLims,[10.5 10.5],'-.k')
        
        dAx{iPl}.YAxis.Visible = 'off';
        dAx{iPl}.Box = 'off';
        
        title(plNames{iPl})
        
        % 1 Dimensional extent plot
        cCL = find(ismember(allColLags,plotColLags{iPl}));
        axb = axes('Position',[ [vAx{iPl}.Position(1:2) + [vAx{iPl}.Position(3).*.25 vAx{iPl}.Position(4).*.7 ] ] [ vAx{iPl}.Position(3:4).*[0.5 .2]] ] ,'Box','off');
        tmpD = -alloverThreshC(:,:);
        b = bar(- nanmean(tmpD(cCL,:),2),'FaceColor',cDat(1,:)); hold on
        plot([1 1]', - [nanmean(tmpD(cCL,:),2) - nanstd(tmpD(cCL,:),[],2)./2  ...
             nanmean(tmpD(cCL,:),2) + nanstd(tmpD(cCL,:),[],2)./2  ]',...
             'k','LineWidth',2); 
        ylim([-2 2])
        axb.XAxis.Visible = 'off';
        axb.Box = 'off';
        view([90 -90])
        title('Field Extent')
        

end
% colormap(redbluecmapRory(20,20))
colormap(whitetocol(100,cDat(2,:)));

sgtitle('Value field shape depend on stimulus direction')







%% Show proportion of neurons that care about (directional) velocity


f.Prop.f = figure('Position',[20 20 400 800]);


fS.nDepthBins = 3;
fS.npropBins = 3;

% layNames = {'early','mid','late'};
% layNames = {'early','early-mid','mid-late','late'};
layNames = {'early','early-mid','mid-late','late','','','','',''};

clear splGroups;
% layNames = {'all'};

% create depth groups so that I can make pie charts
nSplits = length(layNames);
clear splGroups
for iM  = 1:4
    nLay = length(rS(iM).s.lp.netS);
    grSize = round(nLay./nSplits);
    cInd =1;
    for iSpl = 1:nSplits-1
        splGroups{iM,iSpl} = cInd : iSpl .* grSize;
        cInd = iSpl.* grSize + 1;
    end
    splGroups{iM,iSpl + 1} = cInd : nLay;
end


clear propNeurs
for iSpl = 1:nSplits
    
    for iM = 1:length(rS)
             
        cL = splGroups{iM,iSpl};
        
%         propNeurs(iM,iSpl) = nanmean(rMat.propOfCorrVELNeur(cL,iM),[2 1])
        propNeurs(iM,iSpl) = nanmean(rMat.propCorrAndVELNeur(cL,iM),[2 1])
    end
   
    subplot(nSplits,1,iSpl)
%     p = pie([nanmean(propNeurs(:,iSpl),1) ; 1 - nanmean(propNeurs(:,iSpl),1)],{'',''})
    p = pie([nanmean(propNeurs(:,iSpl),1) ; 1 - nanmean(propNeurs(:,iSpl),1)])

    set(p,'EdgeColor','k')
    
%     title([layNames{iSpl} ' layers'])
    
end

sgtitle('Proportion of field neurons with field size depending on velocity (P<0.05)')

legend('Velocity Tuning')

colormap(coltowhite(100,[0 0.2 .9]));

%% Make a venn diagram of neuron sensitivity

f.OverallProp.f = figure('Position',[20 20 400 400]);

addpath('D:\Old_D\Programs\Matlab\Utilities\plotting\venn')


disp(' ')
disp(' -------------------- ')
disp('Proportion of velocity neurons overall')
rMat.overallPropNeur = ...
                 sum(rMat.pDistN<0.05 ,[1 2]) ./ sum(~isnan(rMat.pDistN),[1 2]);
nanmean(rMat.overallPropNeur)
nanstd(rMat.overallPropNeur)




disp(' ')
disp(' -------------------- ')
disp('Proportion of velocity neurons overall')
rMat.overallVELneurOfCorrProp = ...
                 sum(rMat.pDistN<0.05 & ...
                     pCLN<0.05 & ...
                     pRLN<0.05,[1 2]) ./ sum(rMat.pDistN<0.05,[1 2]);
nanmean(rMat.overallVELneurOfCorrProp)
nanstd(rMat.overallVELneurOfCorrProp)

disp(' ')
disp(' -------------------- ')
disp('Proportion of row lag neurons overall')
rMat.overallRLneurOfCorrProp = ...
                 sum(rMat.pDistN<0.05 & ...
                     pRLN<0.05 ,[1 2]) ./ sum(rMat.pDistN<0.05,[1 2]);
nanmean(rMat.overallRLneurOfCorrProp)
nanstd(rMat.overallRLneurOfCorrProp)

disp(' ')
disp(' -------------------- ')
disp('Proportion of column lag neurons overall')
rMat.overallCLneurOfCorrProp = ...
                 sum(rMat.pDistN<0.05 & ...
                     pCLN<0.05 ,[1 2]) ./ sum(rMat.pDistN<0.05,[1 2]);
nanmean(rMat.overallCLneurOfCorrProp)
nanstd(rMat.overallCLneurOfCorrProp)




disp(' ')
disp(' -------------------- ')
disp('Proportion of velocity neurons overall')
rMat.overallVELneurAndCorrProp = ...
                 sum(rMat.pDistN<0.05 & ...
                     pCLN<0.05 & ...
                     pRLN<0.05,[1 2]) ./ sum(~isnan(rMat.pDistN),[1 2]);
nanmean(rMat.overallVELneurAndCorrProp)
nanstd(rMat.overallVELneurAndCorrProp)

disp(' ')
disp(' -------------------- ')
disp('Proportion of row lag neurons overall')
rMat.overallRLneurAndCorrProp = ...
                 sum(rMat.pDistN<0.05 & ...
                     pRLN<0.05 ,[1 2]) ./ sum(~isnan(rMat.pDistN),[1 2]);
nanmean(rMat.overallRLneurAndCorrProp)
nanstd(rMat.overallRLneurAndCorrProp)

disp(' ')
disp(' -------------------- ')
disp('Proportion of column lag neurons overall')
rMat.overallCLneurAndCorrProp = ...
                 sum(rMat.pDistN<0.05 & ...
                     pCLN<0.05 ,[1 2]) ./ sum(~isnan(rMat.pDistN),[1 2]);
nanmean(rMat.overallCLneurAndCorrProp)
nanstd(rMat.overallCLneurAndCorrProp)



disp(' ')
disp(' -------------------- ')
disp('Proportion of velocity neurons overall')
rMat.overallVELneurProp = ...
                 sum(pCLN<0.05 & ...
                     pRLN<0.05,[1 2]) ./ sum(~isnan(rMat.pDistN),[1 2]);
nanmean(rMat.overallVELneurProp)
nanstd(rMat.overallVELneurProp)

disp(' ')
disp(' -------------------- ')
disp('Proportion of row lag neurons overall')
rMat.overallRLneurProp = ...
                 sum(pRLN<0.05 ,[1 2]) ./ sum(~isnan(rMat.pDistN),[1 2]);
nanmean(rMat.overallRLneurProp)
nanstd(rMat.overallRLneurProp)

disp(' ')
disp(' -------------------- ')
disp('Proportion of column lag neurons overall')
rMat.overallCLneurProp = ...
                 sum(pCLN<0.05 ,[1 2]) ./ sum(~isnan(rMat.pDistN),[1 2]);
nanmean(rMat.overallCLneurProp)
nanstd(rMat.overallCLneurProp)




[H S] = venn([nanmean(rMat.overallPropNeur), nanmean(rMat.overallCLneurProp), nanmean(rMat.overallRLneurProp)], ...
    [nanmean(rMat.overallCLneurAndCorrProp) nanmean(rMat.overallRLneurAndCorrProp) nanmean(rMat.overallVELneurProp) nanmean(rMat.overallVELneurAndCorrProp)])

% tmpA = 1 + sum(S.ZoneArea) - sum(S.ZonePop)
% tmpA = sum(S.CirclePop) - sum(S.IntersectPop(1:3)) + 2.*S.IntersectPop(end)
tmpA = sum(S.CircleArea) - sum(S.IntersectArea(1:3)) + 2.*S.IntersectArea(end);
% tmpA = sum(S.CirclePop) - sum(S.IntersectPop) 

hold on


t = linspace(0, 2*pi);
r = sqrt(tmpA./pi);
x = r*cos(t) + .138;
y = r*sin(t)+ .045;
patch(x, y, 'g','FaceAlpha',.1)
axis equal


cAx = gca; 
cAx.Visible = 'off';

legend('Body-part centred', 'velocity-sensitive','direction-sensitive','all')


% 
% 
% fS.nBins = 4;
% fS.pltType = 'bar';
% figure,
% BinPlot(rMat.relLayDep,rMat.propOfCorrRLExpNeur,fS,'LineWidth',2)
% xlabel('binned layer depth')
% ylabel('binned proportion of field neurons')
% 
% title('Proportion of field neurons with field size depending on velocity (P<0.05)')





%% Save figures

allFields = fields(f);
for iF = 1:length(allFields)
    cF = allFields{iF};
       
    set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['Results\ForFigures\VelocityDirection\BitsAndPieces\FromMatlab\' cF '_V3.eps'] , 'epsc')
    saveas(f.(cF).f,['Results\ForFigures\VelocityDirection\BitsAndPieces\FromMatlab\' cF '_V3.pdf'] , 'pdf')
end

% % % 
% % % 
% % % %% Boxplot of velocity effect
% % % 
% % % figure,
% % % tmpD = allOverThresh(:,:);
% % % bar(12 - nanmean(tmpD,2)); hold on
% % % plot([1 1 ; 2 2 ; 3 3]', 12 - [nanmean(tmpD,2) - nanstd(tmpD,[],2)./2  ...
% % %              nanmean(tmpD,2) + nanstd(tmpD,[],2)./2  ]',...
% % %              'k','LineWidth',2);
% % % set(gca,'xTickLabels',{'slow','medium','fast'})
% % % xlabel('Incoming stimulus velocity')
% % % ylabel('Extent of limb-centric fields (pixels)');
% % % title('Value fields extend when stimulus velocity increases')
% % % 
% % % %% Boxplot of direction effect
% % % 
% % % 
% % % figure,
% % % tmpD = alloverThreshC(:,:);
% % % bar( - nanmean(tmpD,2)); hold on
% % % plot([1 1 ; 2 2 ; 3 3 ; 4 4 ; 5 5]', - [nanmean(tmpD,2) - nanstd(tmpD,[],2)./2  ...
% % %              nanmean(tmpD,2) + nanstd(tmpD,[],2)./2  ]',...
% % %              'k','LineWidth',2); 
% % % set(gca,'xTickLabels',{'Fast Left','Slow Left','Front','Slow Right','Fast Right'})
% % % xlabel('Incoming stimulus Direction')
% % % ylabel('Displacement of limb-centric fields (pixels)');
% % % title('Value fields expand in the direction of the oncoming stimulus')
% % % 

% % % 
% % % %% Plot of proportion of neurons that expand their reciptive fields
% % % 
% % % fS.nDepthBins = 3;
% % % fS.npropBins = 3;
% % % figure,
% % % 
% % % % layNames = {'early','mid','late'};
% % % layNames = {'early','early-mid','mid-late','late'};
% % % 
% % % clear splGroups;
% % % % layNames = {'all'};
% % % 
% % % % create depth groups so that I can make pie charts
% % % nSplits = length(layNames);
% % % clear splGroups
% % % for iM  = 1:4
% % %     nLay = length(rS(iM).s.lp.netS);
% % %     grSize = round(nLay./nSplits);
% % %     cInd =1;
% % %     for iSpl = 1:nSplits-1
% % %         splGroups{iM,iSpl} = cInd : iSpl .* grSize;
% % %         cInd = iSpl.* grSize + 1;
% % %     end
% % %     splGroups{iM,iSpl + 1} = cInd : nLay;
% % % end
% % % 
% % % 
% % % clear propNeurs
% % % for iSpl = 1:nSplits
% % %     
% % %     for iM = 1:length(rS)
% % %              
% % %         cL = splGroups{iM,iSpl};
% % %         
% % %         propNeurs(iM,iSpl) = nanmean(rMat.propOfCorrVELNeur(cL,iM),[2 1])
% % % %         propNeurs(iM,iSpl) = nanmean(rMat.propCorrNeur(cL,iM),[2 1])
% % % %         propNeurs(iM,iSpl) = nanmean(rMat.propOfCorrRLExpNeur(cL,iM),[2 1]);
% % % % % %         propNeurs(iM,iSpl) = nanmean(rMat.propOfCorrRLNeur(cL,iM),[2 1])
% % % %         propNeurs(iM,iSpl) = nanmean(rMat.propCorrAndRLNeur(cL,iM),[2 1])
% % % %         propNeurs(iM,iSpl) = nanmean(rMat.propCorrAndRLExpNeur(cL,iM),[2 1])
% % %         
% % %     end
% % %    
% % %     subplot(nSplits,1,iSpl)
% % %     p = pie([nanmean(propNeurs(:,iSpl),1) ; 1 - nanmean(propNeurs(:,iSpl),1)])
% % % %     pie([nanmean(propNeurs(:,iSpl),1) ])
% % % 
% % %     set(p,'EdgeColor','none')
% % %     
% % %     title([layNames{iSpl} ' layers'])
% % %     
% % % end
% % % 
% % % sgtitle('Proportion of field neurons with field size depending on velocity (P<0.05)')
% % % 
% % % legend('speed tuning','no speed tuning')
% % % 
% % % colormap(coltowhite(100,[0 0.2 .9]));
% % % 
% % % 
% % % 
% % % fS.nBins = 4;
% % % fS.pltType = 'bar';
% % % figure,
% % % BinPlot(rMat.relLayDep,rMat.propOfCorrRLExpNeur,fS,'LineWidth',2)
% % % xlabel('binned layer depth')
% % % ylabel('binned proportion of field neurons')
% % % 
% % % title('Proportion of field neurons with field size depending on velocity (P<0.05)')
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % %% Plot of proportion of neurons that shift their reciptive fields
% % % 
% % % fS.nDepthBins = 3;
% % % fS.npropBins = 3;
% % % % figure,
% % % 
% % % % layNames = {'early','mid','late'};
% % % layNames = {'early','early-mid','mid-late','late'};
% % % 
% % % % $$$ INSERT THIS:
% % % figure,plot(rMat.relLayDep,rMat.propCorrAndRLNeur)
% % % 
% % % % create depth groups so that I can make pie charts
% % % nSplits = length(layNames);
% % % clear splGroups
% % % for iM  = 1:4
% % %     nLay = length(rS(iM).s.lp.netS);
% % %     grSize = round(nLay./nSplits);
% % %     cInd =1;
% % %     for iSpl = 1:nSplits-1
% % %         splGroups{iM,iSpl} = cInd : iSpl .* grSize;
% % %         cInd = iSpl.* grSize + 1;
% % %     end
% % %     splGroups{iM,iSpl + 1} = cInd : nLay;
% % % end
% % % 
% % % 
% % % clear propNeurs
% % % for iSpl = 1:nSplits
% % %     
% % %     for iM = 1:length(rS)
% % %              
% % %         cL = splGroups{iM,iSpl};
% % % %         propNeurs(iM,iSpl) = nanmean(rMat.propCorrNeur(cL,iM),[2 1])
% % % %         propNeurs(iM,iSpl) = nanmean(rMat.propCorrAndCLNeur(cL,iM),[2 1])
% % %         propNeurs(iM,iSpl) = nanmean(rMat.propOfCorrCLNeur(cL,iM),[2 1])
% % %         
% % %     end
% % %    
% % %     subplot(nSplits,1,iSpl)
% % %     p = pie([nanmean(propNeurs(:,iSpl),1) ; 1 - nanmean(propNeurs(:,iSpl),1)])
% % % %     pie([nanmean(propNeurs(:,iSpl),1) ])
% % % 
% % %     set(p,'EdgeColor','none')
% % %     
% % %     title([layNames{iSpl} ' layers'])
% % %     
% % % end
% % % 
% % % sgtitle('Proportion of field neurons with field that shift with stimulus direction (P<0.05)')
% % % 
% % % legend('direction tuning','no direction tuning')
% % % 
% % % colormap(coltowhite(100,[0 0.2 .9]));
% % % 
% % % 
% % % 
% % % fS.nBins = 4;
% % % fS.pltType = 'bar';
% % % figure,
% % % BinPlot(rMat.relLayDep,rMat.propOfCorrCLNeur,fS,'LineWidth',2)
% % % xlabel('binned layer depth')
% % % ylabel('binned proportion of field neurons')
% % % 
% % % title('Proportion of field neurons with field that shift with stimulus direction (P<0.05)')


%% $$$ SOME NEURAL MEASURES



% %% Lineplot of velocity effect OLD
% 
% % plotRowLags = [1 3];
% % plotColLags = [2 0 -2];
% plotRowLags = [1 2 3];
% plotColLags = [2 1 0 -1 -2];
% nRL = length(plotRowLags);
% nCL = length(plotColLags);
% 
% 
% allVel = sqrt(plotRowLags.^2 + (plotColLags').^2);
% for iRL = 1:nRL
%     for iCL = 1:nCL
%          cVel = allVel(iRL,iCL);
%     end
% end
% 
% %% $$$ NOT WORKING AS INTENDED --> maybe just use the distance at which Q 
% %      is above some kind of baseline value as the 'extent' instead?
% 
% sFP.plt.fitSigmoid = 1
% sFP.plt.stimRow    = 1:12;
% sFP.plt.stimCol    = 8;
% sFP.plt.lmbCol    = 8;
% 
% for iRL = 1:3
% 
% [rDistQ2, pDistQ2, ...
%     hProxQ2, aDQ2, rDQ2, cDQ2, ...
%     sigmQ2, sigMidQ2] = ...
%     CalcDistCorr(sFP,w, nanmean(QallOrig(:,:,:,:,:,...
%             plotColLags(iCL) == allColLags,...
%             plotRowLags(iRL) == allRowLags,  :),8));
%         
%         sigMidQ2.aD
%         
% end
% 



