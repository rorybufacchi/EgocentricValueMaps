addpath(genpath('Scripts\EgocentricValueMaps'))
load('Results\ForFigures\SuccessorState\NetworkAnalysisWorkSpace.mat')

%% Actions chosen by network

fS.gridXstart = -4.5;
fS.gridXstep = 1;
fS.gridYstart = 3.5;
fS.gridYstep = 1;

f.ChAct.f = figure('Position',[20 20 800 400]);

f.ChAct.ax{1} = subplot(1,2,1);
imagesc(squeeze(nanmean(moveType(:,:,:,:,1),[3 4])));
GridOverImage(fS,f.ChAct.ax{1});
caxis([-1 1])
colormap(redbluecmapRory(20,20))
title('GOALS: Limb moves towards stimuli')

box off


f.ChAct.ax{2} = subplot(1,2,2);
imagesc(squeeze(nanmean(moveType(:,:,:,:,2),[3 4]))); 
GridOverImage(fS,f.ChAct.ax{2});
caxis([-1 1])
colormap(redbluecmapRory(20,20))
title('THREATS: Limb moves away from stimuli')

box off

tmpP = f.ChAct.ax{2}.Position;
f.ChAct.ax{3} = axes('Position',[[ tmpP(1:2) + [tmpP(3).*1 tmpP(4).*0 ] ] [ tmpP(3:4).*[0.3 1]] ] ,'Box','off')
f.ChAct.ax{3}.Visible = 'off'
f.ChAct.cb = colorbar('Color','k')
f.ChAct.cb.Position = f.ChAct.cb.Position + [0 0 0.025 0];
f.ChAct.cb.Ticks = [0 1];
f.ChAct.cb.TickLabels = {'Towards','Away'};


%% Separability of network

for iM = [6 9] % 53 FOR BIG NETWORK, 53 is nice
    
    iM
    
    s=rS(iM).s;
    w=rS(iM).w;
    net=rS(iM).net;
    
    sFP = s;
    
    %     sFP.nta.comparison = 'BodyPart';
    sFP.nta.comparison = 'Valence';
    
    % 'Absolute' 'Row' 'Column' 'AbsRow' 'AbsColumn'
    %     sFP.plt.distanceType = 'AbsColumn';
    sFP.plt.distanceType = 'Absolute';
    %     sFP.plt.distanceType = 'Row';
    %     sFP.plt.distanceType = 'AbsRow';
    
    sFP.plt.startLayer = 1;
    sFP.plt.stopLayer  = length(s.lp.netS) ;
    
    sFP.plt.ON = 0;
    
    [netAR,s] = PlotConnections(net,sFP,w);
    
    
    % ===========================================================
    % Try to do a proximity 'colour' thing
    
    %     close all
    
    %     f.NetSep.f = figure('Position',[20 20 800 800]),
    f1 = figure,
%     p = plot(netAR.G,'Layout','force','WeightEffect','inverse','UseGravity',...
%         'on','ArrowSize',0,'NodeLabel',[],'EdgeColor',[.8 .8 .8],'NodeColor','none');
    p = plot(netAR.G,'Layout','force','WeightEffect','inverse','UseGravity',...
        'on','ArrowSize',0,'EdgeColor',[.8 .8 .8],'NodeColor','none');
    hold on
    
    % Make an indicator of layer depth
    tmp = netAR.A_B_rat';
    layNum = [1:size(netAR.A_B_rat,1)] .*    ones(size(tmp));
    layNum = layNum(:);
    layNum(isnan(tmp(:))) = [];
    
    A_B_rat = NanRemFlatten(netAR.A_B_rat'); classType = 'rat';
    %     A_B_rat = NanRemFlatten(netAR.A_B_diff'); classType = 'diff';
    %     A_B_rat = NanRemFlatten(netAR.AorB'); classType = 'bool';
    G = netAR.G;
    
    
    % First sort nodes by distance to all other nodes (i.e. make a distance
    % matrix)
    clear nodeDists closestDists closestNodes
    scatter(p.XData,p.YData,200,[zeros(size(A_B_rat)) 1-A_B_rat A_B_rat], ...
        'o','Filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',.9,'MarkerEdgeColor','k','LineWidth',2)
    for iNode=1:size(G.Nodes,1)
        nodeDists(:,iNode)    = sqrt(sum(([p.XData(iNode),p.YData(iNode)] - [p.XData(:),p.YData(:)])'.^2));
        [closestDists(:,iNode) closestNodes(:,iNode)] = sort(nodeDists(:,iNode));
    end
    
    
    tmpX = p.XData;
    tmpY = p.YData;
    
end

p.LineWidth = .1;

tmpXold = tmpX;
tmpYold = tmpY;

% axis square
SquareSymAxes
axis square


% theta = deg2rad(55);
theta = deg2rad(0);
% theta = deg2rad(-30);
% theta = deg2rad(-45);
% theta = deg2rad(-60);
% theta = deg2rad(-130);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

rotXY = R*[tmpXold;tmpYold];

tmpX = rotXY(1,:);
tmpY = rotXY(2,:);

rotG = G;

% Make it into a 2d plot and convolve with a 2D kernel


% For non-rotated x and y
nBinsX = 100;
nBinsY = 100;
sgm = 5;
sz  = 15;


myGausFilt = images.internal.createGaussianKernel([sgm sgm], [sz sz]);

binXw = (max(tmpX) - min(tmpX))./nBinsX;
binYw = (max(tmpY) - min(tmpY)) ./nBinsY;

% Binned neural preferences
% binNP  = zeros([nBinsX nBinsY]);
binNP  = nan([nBinsX nBinsY]);
% binNPF = zeros([nBinsX nBinsY]);
binNPF = nan([nBinsX nBinsY]);
binCsx = linspace(min(tmpX),max(tmpX),nBinsX);
binCsy = linspace(min(tmpY),max(tmpY),nBinsY);

for iBx = 1:nBinsX
    inX =   tmpX >= binCsx(iBx) - binXw ./2 & ...
        tmpX <= binCsx(iBx) + binXw ./2 ;
    for iBy = 1:nBinsY
        inY =   tmpY >= binCsy(iBy) - binYw ./2 & ...
            tmpY <= binCsy(iBy) + binYw ./2 ;
        
        binNP(iBx,iBy) = nansum(A_B_rat(inX & inY) - .5);
    end
end


%  Pad binNP so that it can be convolved nicely
tmpPad = nan(size(binNP) + [1 1] .*2 .* ceil(sz/2));
tmpPad(ceil(sz/2) + 1 : end - ceil(sz/2) , ceil(sz/2) + 1 : end - ceil(sz/2) ) = binNP;


% Make my own Gaussian filter with hookers and cocaine
basePoint = ceil(sz/2) + 1;


for iBx = 1:nBinsX
    for iBy = 1:nBinsY
        
        binNPF(iBx,iBy) = nansum(tmpPad(basePoint + [iBx - floor(sz/2) : iBx + floor(sz/2)] , ...
            basePoint + [iBy - floor(sz/2) : iBy + floor(sz/2)]  ) .* ...
            myGausFilt, 'all');
        
    end
end

% Find classification distance to neighbours
for iBx = 2:nBinsX-1
    for iBy = 2:nBinsY-1
        tmpNeighbours = binNPF(iBx + [-1 1],iBy + [-1 1]);
        binNPFdist(iBx,iBy) = sum((binNPF(iBx,iBy) - tmpNeighbours(:)).^2);
    end
end
sumBinNPFdist = sum(binNPFdist(:))

f.NetSep.f = figure('Position',[20 20 800 800]),
ImagescInvisNan(1,binCsx,binCsy,binNPF')
axis xy
colormap(redbluecmapRory(5,6))
SymColAxis;
hold on

p = plot(netAR.G,'Layout','force','WeightEffect','inverse','UseGravity',...
    'on','ArrowSize',0,'NodeLabel',[],'EdgeColor',[.5 .5 .5],'NodeColor','none');

p.LineWidth = .1;


colormap jet
hold on

% scatter(tmpX ,tmpY ,150,[A_B_rat zeros(size(A_B_rat)) 1-A_B_rat], ...
%     'o','Filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',.8,'MarkerEdgeColor','k',...
%     'LineWidth',2)
tmpC = [A_B_rat min([A_B_rat , 1 - A_B_rat],[],2) 1-A_B_rat];
tmpC = tmpC ./ max(tmpC,[],2).*.9;
scatter(tmpX ,tmpY ,250,tmpC, ...
    'o','Filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',.8,'MarkerEdgeColor','k',...
    'LineWidth',2)

colormap(redbluecmapRory(100,10))

SquareSymAxes
axis square

box off

title('Networks separate into functional sub-networks')

f.NetSep.ax{1} = gca;
f.NetSep.ax{1}.Visible = 'off';

tmpP = f.NetSep.ax{1}.Position;
f.NetSep.ax{2} = axes('Position',[[ tmpP(1:2) + [tmpP(3).*.85 tmpP(4).*0 ] ] [ tmpP(3:4).*[0.2 1]] ] ,'Box','off')
f.NetSep.ax{2}.Visible = 'off'
f.NetSep.cb = colorbar('Color','k')
f.NetSep.cb.Position = f.NetSep.cb.Position + [0 0 0.025 0];
f.NetSep.cb.Ticks = [0 1];
f.NetSep.cb.TickLabels = {'Goal preferring','Threat preferring'};



%% General separability



netShapes = unique(allNetW);
nNS = length(netShapes)
f.NetSepHist.f = figure('Position',[20 20 800 400]),
for iNS = 1:nNS
    
    inclM = allNetW == netShapes(iNS);
    
    tmpPM = structureOverChance(inclM,:);
    %     histogram(tmpPM(:),linspace(-4,8,7),'Normalization','probability'); hold on
    histogram(tmpPM(:),'Normalization','probability'); hold on
    tmpPM2(:,iNS) = tmpPM(:);
    title('Functional sub-networks naturally arise, and are more likely to arise when information is condensed')
    xlabel('Increase in network order over chance')
    xlim([-4 8])
    
end

box off

legend('Narrowing networks', 'constant width networks', 'Widening networks')

%% Save figures


allFields = fields(f);
for iF = 1:length(allFields)
    cF = allFields{iF};
    
    set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['Results\ForFigures\SuccessorState\BitsAndPieces\FromMatlab\' cF '.eps'] , 'epsc')
    saveas(f.(cF).f,['Results\ForFigures\SuccessorState\BitsAndPieces\FromMatlab\' cF '.pdf'] , 'pdf')
end


%% Flip things

allNeurActThr2  = allNeurActThr;
allNeurActGl2   = allNeurActGl;
allNeurActGl    = allNeurActThr2;
allNeurActGlThr = allNeurActGl2;

QGl     = QThr2;
QThr    = QGl2;


%% Original and reconstructed QThr values (RUN first part of SUCCESSORSTATE FIRST, or load file below)
% load('Results\ForFigures\SuccessorState\ReconstructingWorkSpace_V2.mat')
load('Results\ForFigures\SuccessorState\ReconstructingWorkSpace_V3.mat')

fS.gridXstart = -4.5;
fS.gridXstep = 1;
fS.gridYstart = 3.5;
fS.gridYstep = 1;

Psii=[];

for iTL = 1:6 % time lag $$$ --> making stats for all 7 properly now [for impact prediction] $$$
cTL = timeLagsToTest(iTL);

% cAxes = [0 0.45 ; -4 4; -2 2];
cAxes = [0 0.5 ; -2 2; -2 2];
xAxes = [0 0.5 ; -.5 2; -2 .5];
yAxes = [0 0.5 ; -.5 2; -2 .5];
convN = {'HitProbability', 'ThreatValue','GoalValue'};
maxCols = [.5 0 .5 ; 1 0 0 ; 0 0 1];
xLims = [4.5 13.5];
yLims = [3.5 13.5];

for iM = 1
    
% % %     % Use only PPS neurons
% % %     QGl = permute(allNeurActGl(:,:,:,:,:,:,iM),[3 4 1 2 5 6]);
% % %     for iL = 1:size(rMat.pDistN(:,:,3),1)
% % %         inclN = rMat.pDistN(iL,:,iM) < .05;
% % % %         QGl(:,:,:,:,iL,~inclN) = NaN;
% % %     end
% % %     QGl(:,:,:,:,1:end-floor(size(QGl,5)/2),:)=[]; % Remove early layers
% % %     QGl = QGl(:,:,:,:,:);
    
    % Use only action values (i.e. outputs)
    QGl = QGl2(:,:,:,:,:,iM);
    
    
% Remove Q goal values for which there are NaNs
nanQ = squeeze(max(isnan(QGl),[],[1 2 3 4]));
QGl(:,:,:,:,nanQ) = [];
    
    
for iConversion = 1:3
    
    f2.(convN{iConversion}).f = figure('Position',[20 -20 800/3 1600]);
    
    clear Psii
    
    for iLimbCol=1:length(lmbColsToTest)
        
        % $$$ FIsrst calculate the hitprob...
        % cLR=1; cLC=5; cGlR=1:s.wrld.size(1); cGlC=5;
        cLR=1; cLC=lmbColsToTest(iLimbCol); cGlR=12; cGlC=5;
        
        % FIrst find a particular state cS
        % % S=[w.lmb.col w.goal.row w.goal.col w.thr.row w.thr.col]
        % cS=tmpS(ismember(tmpS(:,2),cGlR) & ismember(tmpS(:,3),cGlC) & ismember(tmpS(:,1),cLC)   ,:);
        
        % Psi=cat(3,squeeze(QGl(cLR,cLC,:,:,:)), squeeze(QThr(cLR,cLC,:,:,:)));
        Psi=squeeze(QGl(cLR,cLC,:,:,:));
        
        Psii(:,:,:,iLimbCol)=Psi;
        
        if iConversion == 2
            hpmaTmpp=squeeze(QThr2(cLR,cLC,:,:,2,iM));
            hpmaFl2(:,:,iLimbCol)=-hpmaTmpp;
        elseif iConversion == 3
            hpmaTmpp=squeeze(QGl2(cLR,cLC,:,:,2,iM));
            hpmaFl2(:,:,iLimbCol)=-hpmaTmpp;
        end
        
    end
    

    % Remove edge locations, so that non-useful numbers don't go into the
    % correlation
% % %     hpmaFl2Tmp = hpmaFl2;
% % %     PsiiTmp = Psii;
    hpmaFl2Tmp = hpmaFl2(1:12,2:end-1,:);
    PsiiTmp = Psii(1:12,2:end-1,:,:);
    if iConversion == 1
% % %         hpmaFl2 = hitProbMatAll(:,:,:,iTL);
        hpmaFl2Tmp = hitProbMatAll(1:12,2:end-1,:,iTL);
        hpmaFl = hpmaFl2Tmp(:);
    elseif iConversion == 2 | iConversion == 3
        hpmaFl=hpmaFl2Tmp(:);
    end
    
    
    PsiiFl=permute(PsiiTmp,[3 1 2 4]); PsiiFl=PsiiFl(:,:);
    
    corrFacts=(PsiiFl'\hpmaFl);
    
    [rho p]=corr(hpmaFl,PsiiFl'*corrFacts);
    rhoAll(iConversion,iM,iTL) = rho;
    pAll(iConversion,iM,iTL)   = p;
    f2.(convN{iConversion}).ax{3} = subplot(3,1,3);
    plot(hpmaFl,PsiiFl'*corrFacts,'.','MarkerEdgeColor',[.5 .5 .5]); hold on
    title([convN{iConversion} ' reconstruction. rho=' num2str(rho) '. p=' num2str(p)])
    if iConversion == 1
        xlabel(['real ' num2str(cTL +1) '-step hit probabilities'])
        ylabel(['reconstructed ' num2str(cTL +1) '-step ' convN{iConversion} 's']);
    elseif iConversion == 2 | iConversion == 3
        xlabel(['real ' convN{iConversion} 's'])
        ylabel(['reconstructed ' convN{iConversion} 's']);
    end
    SquareAxes ; axis square; box off
    xlim(xAxes(iConversion,:));
    ylim(yAxes(iConversion,:));
    xx = xlim;
    plot(xx,xx,'-.k')
    
    hitProbEstimate=squeeze(sum(Psii.*permute(corrFacts,[2 3 1]),3));
    
    f2.(convN{iConversion}).ax{2} = subplot(3,1,2)
    tmpHC=8;
    imagesc(hitProbEstimate(:,:,tmpHC));
    GridOverImage(fS,f2.(convN{iConversion}).ax{2});
    xlim(xLims);
    ylim(yLims);
    caxis(cAxes(iConversion,:))
    title([convN{iConversion} ' Reconstruction at limb col = ' num2str(tmpHC)])
    axis square; box off
    
    f2.(convN{iConversion}).ax{1} = subplot(3,1,1);
    imagesc(hpmaFl2(:,:,tmpHC));    
    GridOverImage(fS,f2.(convN{iConversion}).ax{1});
    xlim(xLims);
    ylim(yLims);
    caxis(cAxes(iConversion,:))
    title(['Real ' convN{iConversion} ' at limb col = ' num2str(tmpHC)])
    axis square; box off
    
    if iConversion == 1
        colormap(whitetocol(100,maxCols(iConversion,:)));
    elseif iConversion == 2  | iConversion == 3
        colormap(redbluecmapRory(20,20))
    end
    
    tmpP = f2.(convN{iConversion}).ax{1}.Position;

    f2.(convN{iConversion}).ax{3} = axes('Position',[[ tmpP(1:2) + [tmpP(3).*0 tmpP(4).*1.1 ] ] [ tmpP(3:4).*[1 1.5]] ] ,'Box','off')
    f2.(convN{iConversion}).ax{3}.Visible = 'off'
    f2.(convN{iConversion}).cb = colorbar('south','Color','k')
%     f2.(convN{iConversion}).cb.Position = f2.(convN{iConversion}).cb.Position + [0 0 0.025 0];
    caxis([cAxes(iConversion,:)]);
    f2.(convN{iConversion}).cb.Ticks = cAxes(iConversion,:);
    
end
end
end

%% Do fdr stats

rhoAll
pAll

for iConv = 1:3
    
    [dmy pAllCorr] = fdr(pAll(iConv,:));
    
    [pMax pMaxInd] = max(pAllCorr(:));
    
    disp(['Pos valence can reconstruct ' convN{iConv} ' : ' ])
    if min(pAllCorr == 0) == 1
        minRho = min(rhoAll(iConv,:));
    else
        minRho = rhoAll(iConv,pMaxInd);
    end
    disp(['P <= ' num2str(pMax) '. |rho| <= |' num2str(minRho) '|' ])


    disp(['mean rho = ' num2str( nanmean(rhoAll(iConv,:,:),'all'))  ])
    disp(['std  rho = ' num2str( squeeze(nanstd(rhoAll(iConv,:,:))))  ])
    
end

%% Save figures


allFields = fields(f2);
for iF = 1:length(allFields)
    cF = allFields{iF};
    
    set(f2.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f2.(cF).f,['Results\ForFigures\SuccessorState\BitsAndPieces\FromMatlab\' cF '_V2.eps'] , 'epsc')
    saveas(f2.(cF).f,['Results\ForFigures\SuccessorState\BitsAndPieces\FromMatlab\' cF '_V2.pdf'] , 'pdf')
end

