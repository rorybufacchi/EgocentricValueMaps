
%% Calculate a mean Q value - replaces all Q values below

% % % load('D:\Old_D\DPPS\DefenseAgent\Results\Net_GoalAndThreat_MultSpeeds_Goal2_Thr-4_Stay-0_1_400Retrains')
thrNet  = net;
thrS    = s;
thrW    = w;

% $$$ SHOULD I TAKE A NETWORK WHICH HAS ONLY SEEN GOALS? --> I thinkk so,
% to reweight into the threats --> I have that anyway, so it's ok, I think

% s = DefaultSettings(s);

% $$$ START AGAIN FROM HERE!
rowLags=[3 2 1];
colLags=[-1 0 1];
iRn=1; % Run number
clear QGltmp QThrtmp allNeurActGltmp allNeurActThrtmp
for iRL=1:length(rowLags)
    for iCL=1:length(colLags)

% % %         load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_v5_randRC.mat');
% % % % % %         load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\Fig_Vel_Dir_Dependence_101Batches_v2');
        net = rS(3).net;
        w   = rS(3).w;
        s   = DefaultSettings(rS(3).s);
        s.plt.rowLag=rowLags(iRL);
        s.plt.colLag=colLags(iCL);
        s.plt.plotThrFl=0;
        [QGltmp(:,:,:,:,:,iRn), tmptmpNA] = CalcNetOutput(s,w,net);
        allNeurActGltmp(:,:,:,:,1:size(tmptmpNA,5),1:size(tmptmpNA,6),iRn) = tmptmpNA;
        allNeurActGltmp(:,:,:,:,size(tmptmpNA,5)+1:end,size(tmptmpNA,6)+1:end,iRn) = NaN;

% % %         load(['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_SuperCompRelearn_B_V2.mat']);
        net = thrNet;
        w   = thrW;
        s   = DefaultSettings(thrS);
        s.plt.rowLag=rowLags(iRL);
        s.plt.colLag=colLags(iCL);
        s.plt.plotThrFl=1;
        [QThrtmp(:,:,:,:,:,iRn),allNeurActThrtmp(:,:,:,:,:,:,iRn)] = CalcNetOutput(s,w,net);        
        
        iRn=iRn+1;
    end
end
QGl=nanmean(QGltmp,6); allNeurActGl=nanmean(allNeurActGltmp,7);
QThr=nanmean(QThrtmp,6); allNeurActThr=nanmean(allNeurActThrtmp,7);

%% Calculate a mean Q value, but with the specifically made data

% load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\SuccessorState\SuccessorState_51Batch_Plus2_OR_Minus2_movecost_NoHist_v3.mat')

% $$$ SHOULD I TAKE A NETWORK WHICH HAS ONLY SEEN GOALS? --> I thinkk so,
% to reweight into the threats --> I have that anyway, so it's ok, I think

% s = DefaultSettings(s);


iRn=1; % Run number
clear QGltmp QThrtmp allNeurActGltmp allNeurActThrtmp


% % %         load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\Fig_Dist_Pos_Dependence_v5_randRC.mat');
% % % % % %         load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\Fig_Vel_Dir_Dependence_101Batches_v2');
net = rS(iRn).net;
w   = rS(iRn).w;
s   = DefaultSettings(rS(iRn).s);

s.plt.plotThrFl=0;
[QGltmp(:,:,:,:,:,iRn), tmptmpNA] = CalcNetOutput(s,w,net);
allNeurActGltmp(:,:,:,:,1:size(tmptmpNA,5),1:size(tmptmpNA,6),iRn) = tmptmpNA;
allNeurActGltmp(:,:,:,:,size(tmptmpNA,5)+1:end,size(tmptmpNA,6)+1:end,iRn) = NaN;

% % %         load(['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_SuperCompRelearn_B_V2.mat']);
net = rS(iRn+3).net;
w   = rS(iRn+3).w;
s   = DefaultSettings(rS(iRn+3).s);
s.plt.plotThrFl=1;
[QThrtmp(:,:,:,:,:,iRn),allNeurActThrtmp(:,:,:,:,:,:,iRn)] = CalcNetOutput(s,w,net);        
        


QGl=nanmean(QGltmp,6); allNeurActGl=nanmean(allNeurActGltmp,7);
QThr=nanmean(QThrtmp,6); allNeurActThr=nanmean(allNeurActThrtmp,7);

%%
% % % %% HERE DO THE reweighting thing=
% % % 
% % % % netGT=net;sGT=s;wGT=w;storedEpsGT=storedEps;
% % % 
% % % clear tbl
% % % 
% % % s=DefaultSettings(s);
% % % 
% % % % s.plt.plotThrFl=0;
% % % % [QGl,allNeurActGl] = CalcNetOutput(s,w,net);
% % % % 
% % % % s.plt.plotThrFl=1;
% % % % [QThr,allNeurActThr] = CalcNetOutput(s,w,net);
% % % 
% % % Psi=permute(QGl,[5 1 2 3 4]);
% % % % Psi=permute(squeeze(allNeurActGl(:,:,:,:,end,1:s.lp.netS(end))),[5 3 4 1 2]);
% % % PsiFlat=Psi(:,:);
% % % 
% % % QGlPerm=permute(QGl,[5 1 2 3 4]);
% % % QGlFlat=QGlPerm(:,:);
% % % 
% % % QThrPerm=permute(QThr,[5 1 2 3 4]);
% % % QThrFlat=QThrPerm(:,:);
% % % 
% % % tbl=table(QThrFlat(1,:)','VariableNames',{'Q1'});
% % % tbl.Q2=QThrFlat(2,:)'; tbl.Q3=QThrFlat(3,:)';
% % % 
% % % tbl.Psi1=PsiFlat(1,:)';
% % % tbl.Psi2=PsiFlat(2,:)';
% % % tbl.Psi3=PsiFlat(3,:)';
% % % % tbl.Psi4=PsiFlat(4,:)';
% % % % tbl.Psi5=PsiFlat(5,:)';
% % % 
% % % tbl.QGl1=QGlFlat(1,:)';
% % % tbl.QGl2=QGlFlat(2,:)';
% % % tbl.QGl3=QGlFlat(3,:)';
% % % 
% % % % glme1= fitglme(tbl,' Q1 ~ -1 + Psi1 + Psi2 + Psi3 + Psi4 + Psi5 ' );
% % % % glme2= fitglme(tbl,' Q2 ~ -1 + Psi1 + Psi2 + Psi3 + Psi4 + Psi5 ' );
% % % % glme3= fitglme(tbl,' Q3 ~ -1 + Psi1 + Psi2 + Psi3 + Psi4 + Psi5 ' );
% % % 
% % % glme1= fitglme(tbl,' Q1 ~ -1 + QGl1 + QGl2 + QGl3 ' );
% % % glme2= fitglme(tbl,' Q2 ~ -1 + QGl1 + QGl2 + QGl3 ' );
% % % glme3= fitglme(tbl,' Q3 ~ -1 + QGl1 + QGl2 + QGl3 ' );
% % % 
% % % w1=glme1.Coefficients.Estimate;
% % % w2=glme2.Coefficients.Estimate;
% % % w3=glme3.Coefficients.Estimate;
% % % 
% % % Q1=squeeze(QThrPerm(1,:,:,:,:));
% % % Q2=squeeze(QThrPerm(2,:,:,:,:));
% % % Q3=squeeze(QThrPerm(3,:,:,:,:));
% % % 
% % % reconstructQThr=cat(5,glme1.predict,glme2.predict,glme3.predict);
% % % reconstructQThr=reshape(reconstructQThr,[size(QGl)]);
% % % 
% % % QThrV2=cat(5,Q1,Q2,Q3);



% % % %% Plot the reconstructed Q-values
% % % figure,
% % % for iQ=1:3
% % % subplot(1,3,iQ)
% % % tmpRQ=reconstructQThr(:,:,:,:,iQ);
% % % QQ=squeeze(QThrPerm(iQ,:,:,:,:));
% % % plot(QQ(:),tmpRQ(:),'.'); 
% % % % [p h]=signrank(QQ(:),tmpRQ(:));
% % % % [h p]=ttest(QQ(:),tmpRQ(:));
% % % [r p]=corrcoef(QQ(:),tmpRQ(:));
% % % title(['Q' num2str(iQ) '. p=' num2str(p(2))]);
% % % hold on;
% % % xLIM=xlim; yLIM=ylim; alLim=[xLIM;yLIM];
% % % unitLine=[max(alLim(:,1)) min(alLim(:,2))];
% % % plot(unitLine,unitLine) 
% % % xlabel('NN Q'); ylabel('Reconstructed Q')
% % % % hist3([tbl.Q1, glme1.predict],[50 50])
% % % legend('data','x==y')
% % % end
% % % 
% % % s.plt.rowLims=[0.5 14.5];
% % % 
% % % % s.plt.intrpFact=0.1;
% % % s.plt.intrpFact=1;
% % % 
% % % s.plt.plAct=2;
% % % s.plt.lmbCol=1:15;
% % % s.plt.meanLimbCols=1;
% % % s.plt.OneActFl=0;
% % % 
% % % s.plt.OneActFl=1;
% % % s.plt.pltType = 'Imagesc';
% % % 
% % % 
% % % s = DefaultSettings(s)
% % %     
% % % figure,
% % % subplot(3,1,1)
% % % DisplActValsFun(s,w,QThr); title('Base Q Thr')
% % % caxis([-4 4])
% % % subplot(3,1,2)
% % % DisplActValsFun(s,w,reconstructQThr); title('Reconstructed Q Thr')
% % % caxis([-4 4])
% % % subplot(3,1,3)
% % % DisplActValsFun(s,w,QGl); title('Base Q Gl')
% % % caxis([-4 4])
% % % 
% % % colormap(redbluecmapRory(20,20))

% % % %% Calculate p(st+1|st,a) from the Q-matrix (over all states) and Psi
% % % 
% % % 
% % % 
% % % % $$$ I need to make extended Qs, which include the speeds maybe?
% % % % $$$ OR just use the single speed Qs
% % % 
% % % 
% % % % QQ=[5 3 4 ; 2 4 4 ; 1 3 6 ; 9 8 4]'
% % % % Psi=[2 4 6 ]'
% % % 
% % % % current limb and goal columns and rows
% % % cLR=5; cLC=5; cGlR=5; cGlC=5;
% % % 
% % % % QQ=QGlFlat;
% % % % Psi=squeeze(QGl(cLR,cLC,cGlR,cGlC,:));
% % % 
% % % % $$$$$ MAYBE try only using a small subsection of all the possible Q
% % % % configurations?
% % % QQ=[QThrFlat ; QGlFlat];
% % % Psi=[squeeze(QGl(cLR,cLC,cGlR,cGlC,:)) ; squeeze(QThr(cLR,cLC,cGlR,cGlC,:))]
% % % 
% % % Psii=[];
% % % QQQ=[];
% % % for cL = 9:-1:7
% % % Q2=permute(squeeze(allNeurActThr(:,:,:,:,cL,1:s.lp.netS(cL))),[5 3 4 1 2]);
% % % QFlat2=Q2(:,:);
% % % QQQ=[QQQ; QFlat2];
% % % Psi2=squeeze(allNeurActGl(cGlR,cGlC,cLR,cLC,cL,1:s.lp.netS(cL)));
% % % Psii=[Psii; Psi2];
% % % end
% % % Psi=Psii;
% % % QQ=QQQ;
% % % 
% % % %%% QQp=Psi
% % % % p=QQ\Psi;
% % % 
% % % [prob,resnorm,residual,exitflag,output,lambda]=lsqnonneg(QQ,Psi);
% % % 
% % % % clear opts
% % % % % opts.LT=true;
% % % % opts.RECT=true;
% % % % [prob]=linsolve(QQ,Psi,opts);
% % % 
% % % % QQ*p
% % % 
% % % pResh=reshape(prob,[size(QGl(:,:,:,:,1))]);
% % % 
% % % idxs=find(pResh~=0);
% % % 
% % % sizes=size(QGl);
% % % [I,J,K,L] = ind2sub(sizes(1:end-1),idxs)
% % % 
% % % unique(pResh)
% % % 
% % % figure,
% % % imagesc(squeeze(nanmean(pResh(1,:,:,:),[1 2])))
% % % % imagesc(squeeze(nanmean(pResh(1,:,:,:),[1 2]))>0)
% % % 

%%
QGl2    = QGl;
QThr2   = QThr;

% % % %%
% % % QGl    = QGl2;
% % % QThr   = QThr2;

%% Calculate actual hit probbilities

lmbColsToTest=2:s.wrld.size(2)-1;
timeLagsToTest = [0 1 2 3 4 5];
maxTL = max(timeLagsToTest);

tmpS=cell2mat(storedEps.S);

hitProbMatAll=[];

for iTL = 1:length(timeLagsToTest)
    cTL = timeLagsToTest(iTL);
    for iLimbCol=1:length(lmbColsToTest)
        
        % $$$ FIsrst calculate the hitprob...
        % cLR=1; cLC=5; cGlR=1:s.wrld.size(1); cGlC=5;
        cLR=1; cLC=lmbColsToTest(iLimbCol); cGlR=12; cGlC=5;

% There was a reward, and the hand was in a particular position
% % % hitS=tmpS(find(storedEps.R>1 & ismember(tmpS(:,1),cLC)  )-1,:);
hitBool = zeros([size(storedEps.R,1) iTL]);
missBool = hitBool;
ccTL = 0:cTL;
for iiTL = 1:length(ccTL)
    hitBool(1 : end - (iiTL-1) , iiTL)  = storedEps.R(iiTL:end ) > 1 ;
    missBool(1 : end - (iiTL-1) , iiTL) = storedEps.R(iiTL:end ) < 1 ;
end
hitBool     = max(hitBool,[],2);
missBool    = min(missBool,[],2); % $$$ MAYBE CHANGE THIS BACK TO MAX?
hitS=tmpS(find( hitBool & ismember(tmpS(:,1),cLC)  )-1,:);
% There NOT was a reward, and the hand was in a particular position
% tmpTimePs=find(abs(storedEps.R)<1 & ismember(tmpS(:,1),cLC)  )-1;
tmpTimePs=find( missBool & ismember(tmpS(:,1),cLC)  )-1;
tmpTimePs(tmpTimePs<=1)=[];
missS=tmpS(tmpTimePs,:);

xHitTmp=unique(hitS(:,7));
yHitTmp=unique(hitS(:,8));
xDimHit=numel(unique(xHitTmp));
yDimHit=numel(unique(yHitTmp));

xMissTmp=unique(missS(:,7));
yMissTmp=unique(missS(:,8));
xDimMiss=numel(unique(xMissTmp));
yDimMiss=numel(unique(yMissTmp));

% hist3([cS(:,2),cS(:,3)],[3 5])
% hist3([missS(:,7),missS(:,8)],[xDimMiss,yDimMiss])
[nMiss c]=hist3([missS(:,7),missS(:,8)],[xDimMiss,yDimMiss]);
[nHit c]=hist3([hitS(:,7),hitS(:,8)],[xDimHit,yDimHit]);

% $$$ THIS IS WRONG when there is acutally infiniteLR - then It is 1 smaller..
missMat=zeros(size(nMiss)+2);
missMat(xMissTmp(1):xMissTmp(end),yMissTmp(1):yMissTmp(end)) = nMiss;
hitMat=zeros(size(nMiss)+2);
hitMat(xHitTmp,yHitTmp) = nHit;

hitProbMat=hitMat./(hitMat+missMat);
hitProbMat(isnan(hitProbMat))=0;

% % % figure,
% % subplot(3,1,1);
% % imagesc(missMat); colorbar; title('missMat')
% % subplot(3,1,2);
% % imagesc(hitMat); colorbar; title('hitMat')
% % subplot(3,1,3);
% % imagesc(hitProbMat); colorbar; title('hitProbMat')

% % % hitProbMatAll=[hitProbMatAll; hitProbMat(:)];
% % % 
% % % Psi=cat(3,squeeze(QGl(cLR,cLC,:,:,:)), squeeze(QThr(cLR,cLC,:,:,:)));
% % % Psi=permute(Psi,[3 1 2]); Psi=Psi(:,:);
% % % 
% % % Psii=[Psii, Psi];

hitProbMatAll(:,:,iLimbCol,iTL) = hitProbMat;


end

end

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
    sFP.plt.lmbCol = 3:s.wrld.size(2)-2;
    sFP.plt.ON = 0;
    sFP.plt.sequentialLimbCols = 0;
    sFP.plt.stimRow = [3:size(w.world2D,1)-1];
    sFP.plt.stimCol = [2:size(w.world2D,2)-1];
    sFP.plt.pltType = 'Binned';
 
    [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
    
    % This returns a correlation value for each (neuron, layer, limbcol)0.05;
    [rS(cM).rDistN, rS(cM).pDistN, rS(cM).hProxN, aD, rD, cD] = ...
        CalcDistCorr(sFP,w,permute(allNeurAct,[3 4 1 2 5 6]));
    % and this for each action q-value
    [rS(cM).rDistQ, rS(cM).pDistQ, rS(cM).hProxQ, aDQ, rDQ, cDQ] = ...
        CalcDistCorr(sFP,w,Q);
    
end

%% Try to get proximity correlations

[rMat] = DisplayProxStats(rS);

% % % %% Redo the linear regression of p from Psi, but using ALL POSSIBLE psi instead
% % % 
% % % % $$$ !!! the system was never trained to explicitly encode hit-probability,
% % % % yet we can almost perfectly reconstruct the hit probability over a
% % % % timescale which the system doesn't even have access to
% % % 
% % % 
% % % % $$$ MAYBE NEXT SHOW AMOUNT OF CORRELATION WITH HIT PROBABILITY PER LAYER?
% % % % $$$ AND maybe also check/think if I should use a different model? i.e. the already-in-use models for the other valence section?
% % % 
% % % % $$$ MAYBE NEXT SHOW AMOUNT OF CORRELATION WITH HIT PROBABILITY PER LAYER?
% % % 
% % % % Use neural activation if wanted
% % % 
% % % % $$$ Use only PPS neurons
% % % QGl = permute(allNeurActGl(:,:,:,:,:,:),[3 4 1 2 5 6]);
% % % QThr = [];
% % % for iL = 1:size(rMat.pDistN(:,:,3),1)
% % %     inclN = rMat.pDistN(iL,:,3) < .05;
% % %     QGl(:,:,:,:,iL,~inclN) = NaN;
% % % end
% % % QGl(:,:,:,:,1:end-3,:)=[]; % Remove early layers
% % % QGl = QGl(:,:,:,:,:);
% % % QThr = nan(size(QGl));
% % % 
% % % % % % % QGl = permute(allNeurActGl(:,:,:,:,:,:),[3 4 1 2 5 6]); QGl = QGl(:,:,:,:,:);
% % % % % % % QThr = permute(allNeurActThr(:,:,:,:,:,:),[3 4 1 2 5 6]); QThr = QThr(:,:,:,:,:);
% % % % % % lastLay = find(~isnan(squeeze(QGl(2,2,2,2,:,1))),1,'last')
% % % % % % QGl = permute(allNeurActGl(:,:,:,:,lastLay,:),[3 4 1 2 5 6]); QGl = QGl(:,:,:,:,:);
% % % % % % lastLay = find(~isnan(squeeze(QThr(2,2,2,2,:,1))),1,'last')
% % % % % % QThr = permute(allNeurActThr(:,:,:,:,lastLay,:),[3 4 1 2 5 6]); QThr = QThr(:,:,:,:,:);
% % % % % QGl = permute(allNeurActGl(:,:,:,:,1:9,1:12),[3 4 1 2 5 6]); QGl = QGl(:,:,:,:,:);
% % % % % QThr = permute(allNeurActThr(:,:,:,:,1:9,1:12),[3 4 1 2 5 6]); QThr = QThr(:,:,:,:,:);
% % % % QGl = permute(allNeurActGl(:,:,:,:,4,5),[3 4 1 2 5 6]); QGl = QGl(:,:,:,:,:);
% % % % QThr = permute(allNeurActThr(:,:,:,:,4,5),[3 4 1 2 5 6]); QThr = QThr(:,:,:,:,:);
% % % 
% % % % tmpS=cell2mat(storedEps.S);
% % % 
% % % % Remove Q goal values for which there are NaNs
% % % nanQ = squeeze(max(isnan(QGl),[],[1 2 3 4]));
% % % QGl(:,:,:,:,nanQ) = [];
% % % nanQ = squeeze(max(isnan(QThr),[],[1 2 3 4]));
% % % QThr(:,:,:,:,nanQ) = [];
% % % 
% % % Psii=[];
% % % 
% % % for iTL = 1:length(timeLagsToTest)
% % %     cTL = timeLagsToTest(iTL);
% % % for iLimbCol=1:length(lmbColsToTest)
% % % 
% % % % $$$ FIsrst calculate the hitprob...
% % % % cLR=1; cLC=5; cGlR=1:s.wrld.size(1); cGlC=5;
% % % cLR=1; cLC=lmbColsToTest(iLimbCol); cGlR=12; cGlC=5;
% % % 
% % % % FIrst find a particular state cS
% % % % % S=[w.lmb.col w.goal.row w.goal.col w.thr.row w.thr.col]
% % % % cS=tmpS(ismember(tmpS(:,2),cGlR) & ismember(tmpS(:,3),cGlC) & ismember(tmpS(:,1),cLC)   ,:);
% % % 
% % % Psi=cat(3,squeeze(QGl(cLR,cLC,:,:,:)), squeeze(QThr(cLR,cLC,:,:,:)));
% % % % Psi=squeeze(QThr(cLR,cLC,:,:,1));
% % % 
% % % Psii(:,:,:,iLimbCol)=Psi;
% % % 
% % % 
% % % end
% % % 
% % % hpmaFl = hitProbMatAll(:,:,:,iTL);
% % % hpmaFl = hpmaFl(:);
% % % PsiiFl=permute(Psii,[3 1 2 4]); PsiiFl=PsiiFl(:,:);
% % % 
% % % corrFacts=(PsiiFl'\hpmaFl); 
% % % 
% % % [rho p]=corr(hpmaFl,PsiiFl'*corrFacts);
% % % figure,
% % % subplot(3,1,1)
% % % plot(hpmaFl,PsiiFl'*corrFacts,'.');
% % % title(['hit prob reconstruction. rho=' num2str(rho) '. p=' num2str(p)])
% % % xlabel(['real ' num2str(cTL +1) '-step hit probabilities'])
% % % ylabel(['reconstructed ' num2str(cTL +1) '-step hit probabilities']);
% % % 
% % % SquareAxes ; axis square
% % % 
% % % hitProbEstimate=squeeze(sum(Psii.*permute(corrFacts,[2 3 1]),3));
% % % 
% % % subplot(3,1,2)
% % % tmpHC=7;
% % % imagesc(hitProbEstimate(:,:,tmpHC)); colorbar
% % % title(['Hit probability estimate at limb col = ' num2str(tmpHC)])
% % % 
% % % subplot(3,1,3)
% % % imagesc(hitProbMatAll(:,:,tmpHC,iTL)); colorbar
% % % title(['REAL hit probability at limb col = ' num2str(tmpHC)])
% % % 
% % % 
% % % 
% % % end
% % % 
% % % 
% % % 
% % % 
% % % %% PROPERLY calculate Q threat from Q Goal
% % % 
% % % % % Using QGl to predict QThr
% % % % Psii=[];
% % % % for cL=9:-1:1
% % % % Psi2=permute(squeeze(allNeurActGl(:,:,:,:,cL,1:s.lp.netS(cL))),[5 3 4 1 2]);
% % % % PsiFlat2=Psi2(:,:);
% % % % Psii=[Psii, PsiFlat2'];
% % % % end
% % % % QQ=QThrFlat';
% % % % Psii=Psii(:,8:10);
% % % 
% % % % Using QGl to predict QThr
% % % QQ=QThrFlat';
% % % Psii=QGlFlat';
% % % 
% % % % % Using QThr to predict QGl
% % % % QQ=QGlFlat';
% % % % Psii=QThrFlat';
% % % 
% % % 
% % % % Psii=PsiFlat';
% % % % Psii=[PsiFlat', PsiFlat2'];
% % % 
% % % % % %%% QT p= QG
% % % % % p=QT\QP;
% % % % % figure,plot(QP,QT*p,'.')
% % % 
% % % %%% QP p= QT
% % % % corrFacts=QP\QT;
% % % % figure,plot(QT,QP*corrFacts,'.'); hold on
% % % 
% % % corrFacts=(Psii\QQ); 
% % % % $$$ TO BE REMOVED!
% % % corrFacts=nanmean(corrFacts,2);
% % % 
% % % figure,plot(QQ,Psii*corrFacts,'.'); hold on
% % % 
% % % xLIM=xlim; yLIM=ylim; alLim=[xLIM;yLIM];
% % % unitLine=[max(alLim(:,1)) min(alLim(:,2))];
% % % plot(unitLine,unitLine) 
% % % 
% % % reconstructQThr=Psii*corrFacts;
% % % 
% % % 
% % % if size(reconstructQThr,5)==1
% % %     reconstructQThr=reshape(reconstructQThr,[size(QGl(:,:,:,:,1))]);
% % %     reconstructQThr=repmat(reconstructQThr,[1 1 1 1 3]);
% % % else
% % %     reconstructQThr=reshape(reconstructQThr,[size(QGl)]);
% % % end
