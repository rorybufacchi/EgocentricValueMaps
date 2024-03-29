%% Assess the effect of network size on performance, for valence

addpath('F:\Projects\DPPS\DefenseAgent\Scripts\rlsimplepps')
addpath('F:\Projects\DPPS\DefenseAgent\Scripts')
addpath('F:\Programs\Matlab\Utilities\')
addpath('F:\Programs\Matlab\Utilities\plotting\colormaps\')

% % load('F:\Projects\DPPS\DefenseAgent\Results\Performance\Valence\NetSizes\Mults_of_12_51Batch_plus2_plus4_rewards.mat')
% % load('F:\Projects\DPPS\DefenseAgent\Results\Performance\Valence\NetSizes\Mults_of_12_51Batch_plus2_minus4_rewards.mat')

% % load('F:\Projects\DPPS\DefenseAgent\Results\Performance\Valence\NetSizes\Mults_of_12_51Batch_plus2_minus4_rewards_NoHist.mat')

load('F:\Projects\DPPS\DefenseAgent\Results\Performance\Valence\NetSizes\Mults_of_9_51Batch_plus2_minus4_rewards_NoHist_DefRew-01_v2.mat')

% f = figure,
for iM = 1:length(ntRS)
    
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
    
    
    sFP.plt.plotThrFl = 0;
    
    sFP.plt.meanLimbCols = 1;
    sFP.plt.lmbCol=3:12;
%     
%     sFP.plt.meanLimbCols = 0;
%     sFP.plt.lmbCol= 9;
%     
    sFP.plt.plAct=2;
    
    sFP.plt.rowLag = 1;
    sFP.plt.colLag = 0;
    
    sFP.plt.stimRow=[3:size(w.world2D,1)-1];
    sFP.plt.stimCol=[2:size(w.world2D,2)-1];
    %     sFP.plt.pltType = 'Binned';
    
    sFP = DefaultSettings(sFP);
    
    
    % 	set(0, 'currentfigure', f);
    subplot(1 + ceil(length(ntRS)/2),2,1)
    tmp = s.prf.skipBatches;
    plot(ntRS(iM).perf.rewPerAct(1:tmp:end,1),'LineWidth',2);
    hold on;
    
    %     f2 = figure,
    subplot(1 + ceil(length(ntRS)/2),2,1+iM)
    [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
    DisplActValsFun(sFP,w,Q);
    
end
subplot(1 + ceil(length(ntRS)/2),2,1)
% set(0, 'currentfigure', f);
legend('1','2','3','4','5','6','7','8','9','10','11');
% hold off


%% Calculate Q values for different valence values (positive only)


% TEMPORARY:
% % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Plus4_minus01_movecost.mat')

% % % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Plus4_minus01_movecost_NoHist_2_2RandSpr_v3_MoreNetSizes.mat')

load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Plus4_minus01_movecost_NoHist_v2.mat')
% load('F:\Projects\DPPS\DefenseAgent\Results\Performance\FullModel\NetSizes\Body_50_130Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE_FewerSpeeds_SmallWorld_V3.mat')
% rS(end+1) = ntRS;

aTP = [0 1];

for iM = 1:length(rS)
    for iV = 1:length(aTP)
        vRS(iM,iV)=rS(iM);
    end
end

Q = nan([14 15 14 15 9 length(rS) 2]);
for iM = 1:length(rS)
    
    % loop through valence
    for iV = 1:length(aTP)
        
        sFP = DefaultSettings(rS(iM).s);
        w = rS(iM).w;
        net = rS(iM).net;
        Qtable = rS(iM).Qtable;
        
        sFP.plt.meanLimbCols = 1;
        sFP.plt.lmbCol = [2:(sFP.wrld.size(2)-1)];
        sFP.plt.stimCol= [2:(sFP.wrld.size(2)-1)] ;
        sFP.plt.stimRow = [3:size(w.world2D,1)-1];
            
        sFP.plt.plotThrFl = aTP(iV);
        
        [qTmp,tmpNA] = CalcNetOutput(sFP,w,net);
        % % % % %         allNeurAct(:,:,:,:,:,:,iM,iV) = tmpNA;
        
        
        [vRS(iM,iV).rDistN, vRS(iM,iV).pDistN, vRS(iM,iV).hProxN, aD, rD, cD] = CalcDistCorr(sFP,w,permute(tmpNA,[3 4 1 2 5 6]));
        [vRS(iM,iV).rDistQ, vRS(iM,iV).pDistQ, vRS(iM,iV).hProxQ, aDQ, rDQ, cDQ] = CalcDistCorr(sFP,w,qTmp);
        
        
        Q(size(Q,1) + 1 - size(qTmp,1) : size(Q,1),:,...
            size(Q,3) + 1 - size(qTmp,3) : size(Q,3),:,1:size(qTmp,5),iM,iV) = qTmp;
        
        if iM == 1 & iV == 1
            nAall = nan([size(tmpNA(:,:,:,:,1)) max(arrayfun(@(x) length(rS(x).s.lp.netS) , [1:length(rS)])) max(arrayfun(@(x) max(rS(x).s.lp.netS) , [1:length(rS)])) length(rS) 2]);
        end
        nAall(:,:,:,:,1:size(tmpNA,5),1:size(tmpNA,6),iM,iV) = permute(tmpNA,[3 4 1 2 5 6]);
        
        
    end
end

% $$$ Now do I plot the r value as a function of valence or something?

%%

[rMat_lVal] = DisplayProxStats(vRS(:,1));
[rMat_hVal] = DisplayProxStats(vRS(:,2));

%% Check whether higher valence gives higher Q-values (kind of trivial)


% $$$ NEXT FIGURE OUT WHY THE ZVALS FOR HIGHER Q VALUES ARE SO WEAK FOR MODEL 2
% --> DO I JUST DO STATS FOR ALL MODELS TOGETHER INSTEAD OR SOMETHING?...

clear qVec
for iM = 1:length(rS)
    for iAct = 1:rS(1).s.act.numA
        for iV = 1:length(aTP)
            
            tmp = Q(1,2:end-1,1:11,2:end-1,iAct,iM,iV);
            qVec(:,iV) = tmp(:);
        end
        
        % Perform stats to see which stimulus' q-value is higher
        opts.testType = 'signrank';
        [figString zVal(iAct,iM) pVal(iAct,iM)] = CorrelationTitle([],qVec,opts);
        
    end
end

pVal = CalcFDR(pVal);

[pMax pMaxInd] = max(pVal(:));

disp(['Higher reward gives bigger magnitude Q fields: ' ])
disp(['P <= ' num2str(pMax) '. |Z| <= |' num2str(zVal(pMaxInd)) '|' ])

disp('effect size')
zVal ./ sqrt(numel(qVec))

disp('mean Z')
nanmean(zVal(:))

disp('std Z')
nanstd(zVal(:))



% % % %% Use thresholds to see whether space expands
% % % 
% % % 
% % % overThresh = NaN([3 9 4]);
% % % 
% % % clear rV pV allV tmpV allOverThresh overThresh
% % % 
% % % for iM = 1:length(rS)
% % %     
% % %     for iAct = 1:1 %rS(iM).s.act.numA
% % %         
% % %         clear overThresh tmpRL
% % %         
% % %         % lC sR sC cL rL
% % % %         qTmp = squeeze(Qall(end,:,:,:,iAct,:,:,iM));
% % %         qTmp = squeeze(nanmean(Q(end,:,:,:,:,iM,:),5));
% % %         
% % % %         qThresh = prctile(qTmp(:),85);
% % %         
% % %         % Find the extent of the field for each limb column
% % %         lCs = 2:14;
% % %         for iLC = 1:length(lCs)
% % %             cLC = lCs(iLC);
% % %             Qline = squeeze(nanmean(qTmp(cLC,1:11,cLC,:,:),[1 3]));
% % %             
% % %             qThresh = prctile(Qline(:),90); % $$$ SEE IF THIS GIVES SIMILAR RESULTS CAUSE THAT WOULD BE BETTER
% % %             
% % %             figure,plot( (Qline - max(Qline))./ (max(Qline) - min(Qline)) )
% % %             
% % %             
% % %             
% % %             for iV = 1:size(Qline,2) % loop through ROW LAGS to find expansion size
% % %                 try
% % %                     overThresh(iV,iLC) = find(Qline(:,iV) > qThresh,1);
% % %                 catch
% % %                     overThresh(iV,iLC) = NaN;
% % %                 end
% % %                 tmpV(iV,iLC) = iV;
% % %             end
% % %         end
% % %         
% % %         allV(:,:,iAct,iM)           = tmpV;
% % %         allOverThresh(:,:,iAct,iM)  = overThresh;
% % %         
% % %         % Column Lag p-values
% % %         [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(tmpV(~isnan(overThresh)),overThresh(~isnan(overThresh)));
% % %         rV(iAct,iM) = rTmp(2);
% % %         pV(iAct,iM) = pTmp(2);
% % %         
% % %         
% % %     end
% % %     
% % % end
% % % 
% % % rV
% % % pV



%% Use blurred thresholds to plot action value funciton


overThresh = NaN([3 9 4]);

clear rV pV allV tmpV allOverThresh overThresh

for iM = 1:length(rS)
    
    for iAct = 1:1 %rS(iM).s.act.numA
        
        clear overThresh tmpRL
        
        % lC sR sC cL rL
%         qTmp = squeeze(Q(end,:,:,:,iAct,:,:,iM));
        qTmp = squeeze(nanmean(Q(end,:,:,:,:,iM,:),5));
        qTmp = qTmp - min(qTmp(:)); % makeit all positive so that it is easy to take ratios
        
        qThresh = prctile(qTmp(:),90);
        
        % Find the extent of the field for each limb column
        lCs = 2:14;
        for iLC = 1:length(lCs)
            cLC = lCs(iLC);
            Qline = squeeze(nanmean(qTmp(cLC,1:11,cLC,:,:),[1 3]));
            
%             figure,plot( (Qline - max(Qline))./ (max(Qline) - min(Qline)) )
%             figure,plot(Qline ./ (Qline + qThresh))
            
            % relative contact relevance
%             relContRel(:,:,iLC,iM) = (Qline ./ (Qline + qThresh));
%             relContRel(:,:,iLC,iM) = Qline ./ (max(qTmp(:)) + qThresh);
%             relContRel(:,:,iLC,iM) = Qline ./ qThresh;

            % relative contact relevance, with random noise
            relContRel(:,:,iLC,iM) = normcdf(Qline,qThresh,qThresh.*0.25);

            
            for iV = 1:size(relContRel,2)
                % Center of mass
                coM(iV,iLC,iM) = nansum(relContRel(:,iV,iLC,iM).*(1:size(relContRel,1))') ./ nansum(relContRel(:,iV,iLC,iM));
                % Point of first increase
                tmpFI = find(relContRel(:,iV,iLC,iM)>0.5,1);;
                if isempty(tmpFI)
                    pFI(iV,iLC,iM) = NaN;
                else
                    pFI(iV,iLC,iM) = find(relContRel(:,iV,iLC,iM)>0.5,1);
                end
            end



            
        end
        
    end
    
end

figure,plot(nanmean(relContRel(:,:,:),3))

%% Do stats to show that movement happens at a greater distance

disp('-----------------------')
disp('Higher absolute valence leads to more actions')
clear pActDist signRankActDist zActDist
for iM = 1:size(pFI,3)
%     [pTmp,hTmp,stats] = signrank(pFI(2,:,iM)-pFI(1,:,iM),'method','approximate');
    [pTmp,hTmp,stats] = signrank(pFI(2,:,iM),pFI(1,:,iM),'method','approximate');
    pActDist(iM) = pTmp;
    signRankActDist(iM) = stats.signedrank;
    zActDist(iM) = stats.zval;
end

[dmy pActCorr] = fdr(pActDist)
signRankActDist
zActDist



disp('effect size')
zActDist ./ sqrt(numel(pFI(:,:,1)))

disp('mean Z')
nanmean(zActDist(:))

disp('std Z')
nanstd(zActDist(:))

%% Show neural activation is higher when valence is higher?

size(nAall)
% % % testLC = 5:11;
testLC = 2:14;
clear tmpD

% limb row, col, stim row, col, layers, ???, models, valence <--- FIND OUT
% WHAT THE ??? DIMENSION IS --> is it neurons? --> then we would be testing
% for all stimulus positions, but neurons separately
tmpD1 = permute(nAall(w.lmb.row,testLC,2:11,testLC,:,:,:,:),[5 6 7 8  1 2 3 4]);
tmpD1 = tmpD1(:,:,:,:,:);
allSD = nanstd(tmpD1,[],5);




% for iL = 1:size(nAall,5)
%     for iN = 1:size(nAall,6)
%         
%         boxplot(allSD(:,:,1,1)')
%         
%     end
% end



A = (allSD(:,:,:,2) -  allSD(:,:,:,1)) ./ (allSD(:,:,:,2) +  allSD(:,:,:,1));

disp('--------------------')
disp('neural increase in response to higher valence stimuli')
[ppp hhh statss] = signrank(A(~isnan(A)))

disp('effect size')
statss.zval ./ sqrt(sum(~isnan(A),'all'))


% % % % DO it for separate networks?
% % % disp('--------------------')
% % % disp('And separately for each network:')
% % % 
% % % for iM = 1:size(nAall,7)
% % %     clear tmpD
% % %     tmpD1 = permute(nAall(w.lmb.row,testLC,2:11,testLC,:,:,iM,:),[5 6 7 8  1 2 3 4]);
% % %     tmpD1 = tmpD1(:,:,:,:,:);
% % %     allSD = nanstd(tmpD1,[],5);
% % %     A = (allSD(:,:,:,2) -  allSD(:,:,:,1)) ./ (allSD(:,:,:,2) +  allSD(:,:,:,1));
% % %     disp('neural increase in response to higher valence stimuli')
% % %     [ppp hhh statss] = signrank(A(~isnan(A)))
% % % end


% for iV = 1:2
% figure,
% tmpM  = nanmean(allSD(:,:,:,iV),[2 3]);
% tmpSD = nanstd(allSD(:,:,:,iV),[],[2 3]);
% bar(tmpM) % one value for each of the 5 layers
% hold on
%     
% arrayfun(@(n) plot([n n],tmpM(n) + [ - tmpSD(n), tmpSD(n)],'k','LineWidth',2) ,1:length(tmpM));
%         
% end
% 
% 
% figure,
% tmpM  = nanmean(A,[2 3]);
% tmpSD = nanstd(A,[],[2 3]);
% bar(tmpM) % one value for each of the 5 layers
% hold on
%     
% arrayfun(@(n) plot([n n],tmpM(n) + [ - tmpSD(n), tmpSD(n)],'k','LineWidth',2) ,1:length(tmpM));
% 
% A = allSD(:,:,:,1) -  allSD(:,:,:,2);

% figure,histogram(A(:));


%% Calculate proportion of neurons for which valence increases variability

rMat        = rMat_lVal;
for iEl = 1:numel(rMat.hProxN)
    if ~isnan(rMat.hProxN(iEl))
        rMat.hProxN(iEl) = rMat.hProxN(iEl)  & rMat_hVal.hProxN(iEl); 
    end
end

for iM = 1:length(rS)
    
    for iL = 1:length(rS(iM).s.lp.netS)
        for iN = 1:rS(iM).s.lp.netS(iL)
            
            rMat.propOfCorrVarNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(rMat.hProxN(1:length(rS(iM).s.lp.netS),:,iM) == 1 & ...
                A(1:length(rS(iM).s.lp.netS),:,iM) > 0,2) ./ nansum(rMat.hProxN(1:length(rS(iM).s.lp.netS),:,iM),2);
            
            rMat.propCorrAndVarNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(rMat.hProxN(1:length(rS(iM).s.lp.netS),:,iM) == 1 & ...
                A(1:length(rS(iM).s.lp.netS),:,iM) > 0,2) ./ rS(iM).s.lp.netS';
            
            rMat.propVarNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum( A(1:length(rS(iM).s.lp.netS),:,iM) > 0,2) ./ rS(iM).s.lp.netS';
            
        end
    end

end



% % % %%
% % % 
% % % 
% % % tmpD = tmpD(:,:);
% % % 
% % % nanmean(tmpD,2)
% % % nanmedian(tmpD,2)
% % % 
% % % nanstd(tmpD,[],2)
% % % 
% % % 
% % % 
% % % [p h] = signrank(tmpD(1,:),tmpD(2,:))
% % % 
% % % 
% % % % $$$ CHECK WHY THE SIGNRANK WITH AND WITHOUT different models doesn't
% % % % workj --> does it go in dofferent directions?...
% % % % --> SHOULD I USE ABSOLUTE? probably


%% Show neural activation is MORE VARIABLE/ affected by position when valence is higher?





% % %% Find distance at which staying still becomes more valuable than moving, and 
% % %  calculate sigmoid midpoints of best action choice
% % 
% % stayBest = Q(:,:,:,:,2,:,:) > Q(:,:,:,:,1,:,:) & Q(:,:,:,:,2,:,:) > Q(:,:,:,:,3,:,:);
% % 
% % size(stayBest(:,:,1:5,:,:,:,:)) = NaN;
% % 
% % sFP.plt.fitSigmoid = 1;
% % 
% % testLC = [2:sFP.wrld.size(2)-1];
% % 
% % sFP.plt.stimRow = 2:sFP.wrld.size(1)-2;
% % 
% % 
% % clear sigMidQ_Gl sigMidQ_Thr
% % tic
% % for iM = 1:length(rS)
% %     for iLC = 1:length(testLC)
% %         
% %         LC = testLC(iLC);
% %         sFP.plt.lmbCol = LC;
% %         
% %         
% %         [rDistQ, pDistQ, hProxQ, aDQ, rDQ, cDQ, ...
% %             sigmQ_Gl, sigMidQ_Gl(iLC,iM)] = ...
% %             CalcDistCorr(sFP,w,double(stayBest(:,:,:,:,1,iM,1)));
% %         
% %         [rDistQ, pDistQ, hProxQ, aDQ, rDQ, cDQ, ...
% %             sigmQ_Thr, sigMidQ_Thr(iLC,iM)] = ...
% %             CalcDistCorr(sFP,w,double(stayBest(:,:,:,:,1,iM,2)));
% %     end
% %     iM
% % end
% % toc
% % 
% % %% Sanity check plot stayBest
% % 
% % % $$$ MAYBE I SHOULD MAKE STAYRANKING, instead of staybest?
% % 
% % % $$$ OR MAYBE the difference is so small due to the random movement
% % % possibilities?... --> make the directions more causal?
% % 
% % % $$$ ALSO check how many different levels of valence there were and if
% % % there were even threats present in the currently loaded ntRS
% % 
% % iM = 1;
% % iV = 1;
% % 
% % testLC = [2:sFP.wrld.size(2)-1];
% % 
% % clear stayBestLines
% % for iLC = 1:length(testLC)
% %     
% %     cLimb = testLC(iLC);
% %     
% % %     squeeze(nanmean(Q(12,cLimb,:,cLimb,2,:,iV),2) - nanmean(Q(12,cLimb,:,cLimb,[1 3],:,iV),[2 5]) );
% % 
% %     stayBestLines(:,iLC,:,:) = squeeze(stayBest(12,cLimb,2:11,cLimb,1,:,:));
% % 
% % end
% % 
% % %% Plot staybestlines for the valences separately $$$
% % 
% % for iV = 1:2
% %     tmpD = stayBestLines(:,:,:,iV);
% %     tmpD = tmpD(:,:);
% %     
% %     figure,plot(nanmean(tmpD,2),'-x')
% % end
% % 
% % %%
% % 
% % 
% % for iLC = 1:length(testLC)
% %     for iM = 1:length(rS)
% %         for iV = 1:2
% % 
% %  % Fit a sigmoid as well
% %             X = [size(stayBestLines):-1:1]';
% %             Y = double(stayBestLines(:,iLC,iM,iV));
% %                        
% %             fo = fitoptions('Method','NonlinearLeastSquares',...
% %                 'Lower',s.plt.fitSigmLow,...
% %                 'Upper',s.plt.fitSigmUp, ...
% %                 'StartPoint',s.plt.fitSigmStart);
% %             
% %             ft = fittype('a/ (1+exp(-b*(x-c))) + d' , 'options', fo);
% %            
% %             [curve{iLC,iM,iV}, gof{iLC,iM,iV}] = fit(X,Y,ft);
% %             tmp = coeffvalues(curve{iLC,iM,iV});
% %             sigMid(iLC,iM,iV) = tmp(3);
% %             
% %         end
% %     end
% % end

%%

% subplot(2,1,1)
% imagesc(squeeze(Q(12,cLimb,:,:,3,iM,1))); colorbar
% 
% subplot(2,1,2)
% 
% figure,plot(squeeze(Q(12,cLimb,:,:,3,iM,iV)) );
% 
% figure,plot(squeeze(nanmean(Q(12,cLimb,:,cLimb,2,iM,iV),2) - nanmean(Q(12,cLimb,:,cLimb,3,iM,iV),2) ) );
% 
% hold on
% 
% plot(squeeze(nanmean(Q(12,cLimb,:,cLimb,2,iM,iV),2) - nanmean(Q(12,cLimb,:,cLimb,1,iM,iV),2) ) );
% 
% 
% % figure,plot(squeeze(Q(12,cLimb,:,:,2,iM,iV)));
%  
% figure,plot(squeeze(stayBest(12,cLimb,:,cLimb,1,iM,iV)))

% % %% Do stats on the best action choice sigmoid distance
% % 
% % for iM = 1:length(rS)
% %     for iLC = 1:length(testLC)
% %         sMG(iLC,iM) = sigMidQ_Gl(iLC,iM).aD;
% %         sMT(iLC,iM) = sigMidQ_Thr(iLC,iM).aD;
% %     end
% % end
% % 
% % valenceVals = [ones(size(sMT)).*sFP.act.GoalRew; ones(size(sMT)).*sFP.act.ThreatRew];
% % 
% % % opts.testType = 'lme';
% % % [figString zVal pVal] = CorrelationTitle(valenceVals,[sMG; sMT],opts)
% % 
% % figure,
% % plot(valenceVals,[sMG; sMT],'x')
% % 
% % opts.testType = 'signrank';
% % [figString zVal pVal] = CorrelationTitle([],[sMG(:), sMT(:)],opts)
% % 
% % hist([sMG(:), sMT(:)])
% % legend('2 rew','4 rew')
% % title(figString)
% % xlabel('distance at which staying becomes the preferred movement')

%% Check whether threats (negative reward stimuli) also give PPS fields

% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_2Batch_Plus2_Minus2_minus01_movecost_NoHist.mat')
% tmp = rS;


% $$$ THIS TO BE PUT BACK IN
% =========================================================================
% TEMPORARY EXTRA v2:
for iVV = 1 %15 %20
% load(['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_2Batch_Plus2_Minus2_minus01_movecost_NoHist_V' num2str(iVV) '.mat']);
% load(['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Minus2_minus01_movecost_NoHist_V' num2str(iVV) '.mat']);

% % % % % % $$$ THIS IS THE small networks
% % % load(['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Minus2_minus03_movecost_NoHist_V' num2str(iVV) '.mat'])

% % % $$$ THIS IS THE RIGHT ONE! But just check if the smaller networks
% work too

% % % % ====================================
% % % % MAIN TEXT OPTION (TANSIG)
% % % % ====================================
% % % load(['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_SuperCompRelearn_B_V' num2str(iVV) '.mat']);

% % % % ====================================
% % % % RELU (POSLIN)
% % % % ====================================
% % % load(['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_poslin_B_V' num2str(iVV) '.mat'])

% ====================================
% LogSig
% ====================================
load(['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_logsig_B_V' num2str(iVV) '.mat'])

% % % % ====================================
% % % % SoftMax
% % % % ====================================
% % % load(['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_softmax_B_V' num2str(iVV) '.mat'])


% % % % ====================================
% % % % TriBas
% % % % ====================================
% % % load(['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_L1_Regularization_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_tribas_B_V' num2str(iVV) '.mat'])


% % % % ====================================
% % % % RadBas
% % % % ====================================
% % % load(['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_L1_Regularization_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_radbas_B_V' num2str(iVV) '.mat'])


% % % % ====================================
% % % % TANSIG - L1 Regularised
% % % % ====================================
% % % load(['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_L1_Regularization_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_tansig_B_V' num2str(iVV)])

% load(['F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Minus2_minus01_movecost_NoHist_NEWLEARNINGPARAMS_SmallNetRelearn_B_V' num2str(iVV) '.mat']);
if iVV == 1
    tmp = rS;
end
tmp(end+1:end+length(rS)) = rS;
end
rS = tmp;
% =========================================================================

% $$$ REMEMBER TO PUT THIS BACK IN --> COMBO MODEL $$$$$$$$$$$$$$$$$$$$$$$$
% % % % 
% % % load('F:\Projects\DPPS\DefenseAgent\Results\Performance\FullModel\NetSizes\Body_50_130Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE_FewerSpeeds_SmallWorld_V3.mat')
% % % rS(end+1:end+length(ntRS)) = ntRS;



% =========================================================================
% TEMPORARY EXTRA v2b:
% rS = tmp;
% =========================================================================

% % % % =========================================================================
% % % % TEMPORARY EXTRA:
% % % 
% % % % -------------------------------------------------------------------------
% % % % extra model which includes all other stuff (except multple limbs)
% % % load('F:\Projects\DPPS\DefenseAgent\Results\Net_GoalAndThreat_MultSpeeds_Goal2_Thr-4_Stay-0_1_400Retrains');
% % % rS(length(rS)+1).s = s;
% % % rS(end).w = w;
% % % rS(end).net = net;
% % % 
% % % 
% % % 
% % % % -------------------------------------------------------------------------
% % % % extra FULL MODEL (but not completed learning)
% % % load('F:\Projects\DPPS\DefenseAgent\Results\Net_GoalAndThreat_MultSpeeds_Goal2_Thr-4_Stay-0_1_400Retrains_DIFFNETSTRUCT.mat');
% % % 
% % % rS(length(rS)+1).s = s;
% % % rS(end).w = w;
% % % rS(end).net = net;
% % % 
% % % 
% % % % -------------------------------------------------------------------------
% % % % extra almost full model (again)
% % % load('F:\Projects\DPPS\DefenseAgent\Results\Net_GoalAndThreat_FewerSpeeds_Goal2_Thr-2_Stay-0_1_400Retrains_DIFFNETSTRUCT.mat');
% % % rS(length(rS)+1).s = s;
% % % rS(end).w = w;
% % % rS(end).net = net;
% % % 
% % % % =========================================================================







sFP = rS(1).s;

% test Columns (columns at which to test the sigmoids
tC = 2:sFP.wrld.size(2)-1;

% all Threat Presences to test the extents with. 0 is just goal, 1 is just
% 'threat' (which could also have positive valence)
aTP = [0 1];

% Preallocate space for the extent of the value and neural fields
nLay = 1; maxN = 1;
for iM=1:length(rS)
    nLay = max([nLay, length(rS(iM).s.lp.netS)]);
    maxN = max([maxN, max(rS(iM).s.lp.netS)]);
end
qExtent = nan([length(aTP),length(tC),rS(1).s.act.numA,length(rS)]);
nExtent = nan([length(aTP),length(tC),nLay,maxN,length(rS)]);



% Loop through valence values (threat presence) and calculate 'extent of PPS'
for iV = 1:length(aTP)
    rS_sepV(iV).rS = rS;
    rS_sepV(iV).V = iV;
    
    for iM = 1:length(rS)
        
        sFP = DefaultSettings(rS(iM).s); 
        w = rS(iM).w;
        net = rS(iM).net;
        Qtable = rS(iM).Qtable;
        
        sFP.plt.plotThrFl = aTP(iV);
           
        sFP.plt.rowLims = [1.5 sFP.wrld.size(1)-0.5];
        sFP.plt.stimRow=[3:size(w.world2D,1)-3];
        
        % -------------------------------------------------------------------------
        % Calculate the correlati-ness of each velocity, for Q-values and neural
        % activity
        sFP.plt.meanLimbCols = 1;
        sFP.plt.fitSigmoid = 0;
        sFP.plt.lmbCol = [2:(sFP.wrld.size(2)-1)];
        sFP.plt.stimCol= [2:(sFP.wrld.size(2)-1)] ;
        [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
        
        if size(Q,5)>3
            warning('Multiple limbs are part of model, but only assessing actions 1-3');
        end
        if iM == 1 & iV == 1
            Qall(:,:,:,:,:,iM,iV) = Q(:,:,:,:,1:3);
        end
        
        Qall(size(Qall,1) + 1 - size(Q,1) : size(Qall,1),:,...
            size(Qall,3) + 1 - size(Q,3) : size(Qall,3),:,:,iM,iV) = Q(:,:,:,:,1:3);
        
        [rS_sepV(iV).rS(iM).rDistQ, rS_sepV(iV).rS(iM).pDistQ, ...
            rS_sepV(iV).rS(iM).hProxQ] = ...
            CalcDistCorr(sFP,w,Q);
        
        [rS_sepV(iV).rS(iM).rDistN, rS_sepV(iV).rS(iM).pDistN, ...
            rS_sepV(iV).rS(iM).hProxN] = ...
            CalcDistCorr(sFP,w,permute(allNeurAct,[3 4 1 2 5 6]));
        
        % Find best movement type for each limb position
        lC = 2:sFP.wrld.size(2)-1;
        sFP.plt.ON = 0;
        for iLC = 1:length(lC)
            sTmp=sFP;
            sTmp.plt.lmbCol = lC(iLC);
            sTmp.plt.pltType = 'Imagesc';
            [mT squishQs QbyType] = FindMoveTypes(net,sTmp,w);
            % Don't take the middle column, cause it is hard to classify
            moveType(size(Qall,1) + 1 - size(Q,1) : size(Qall,1),:,iLC,iM,iV) = mT(:,[2:7 9:end-1]);
        end
        
        iM
    end
    iV
end


%% Display proximity-dependence stats for goal and treat

% $$$ HERE HERE 
% For all the models
[rMat_Gl] = DisplayProxStats(rS_sepV(1).rS);
[rMat_Thr] = DisplayProxStats(rS_sepV(2).rS);

%% Check whether limbs move more towards or away from stimuli when goals or threats are present
% 
% for iMT = 1:length(rS)
% glMT  = moveType(:,:,:,iMT,1);
% thrMT = moveType(:,:,:,iMT,2);
% 
% opts.testType = 'signrank';
% [figString zVal(iMT) pVal(iMT)] = CorrelationTitle([],[glMT(:), thrMT(:)],opts);
% end
% 
% [pMax pMaxInd] = max(pVal);
% 
% disp('The limb moves more towards goals, and more away from threats:')
% disp(['P <= ' num2str(pMax) '. |Z| <= |' num2str(zVal(pMaxInd)) '|' ])
% 
% %%
% 
% rS(15).perf.rewPerAct(1:25:end,:)
% 
% 
% %%


for iM = 1:length(rS); % 1:size(moveType,4);

glMT2   = moveType(:,:,:,iM,1);
glMT2  = glMT2(:);
thrMT2  = moveType(:,:,:,iM,2);
thrMT2 = thrMT2(:);

% nanmean(glMT2)
% nanmedian(glMT2)
% 
% nanmean(thrMT2)
% nanmedian(thrMT2)
% [hh pp statss] = ttest(glMT2,thrMT2)

[ppp(iM) hhh(iM) statsss] = signrank(glMT2,thrMT2)
zzz(iM) = statsss.zval;

end


glMT2   = moveType(:,:,:,iM,1);
glMT2  = glMT2(:);
thrMT2  = moveType(:,:,:,iM,2);
thrMT2 = thrMT2(:);

[pppAll hhhAll(iM) statsss] = signrank(glMT2,thrMT2)
zzzAll = statsss.zval;


[dmy pppCorr] = fdr(ppp);

[maxP maxInd] = max(pppCorr(:));
disp(['Grab-ness difference pos neg. Maximum p for Q-values: P < =  ' num2str(maxP,3) '. |Z| > =  |' num2str(zzz(maxInd),3) '|' ])
disp(' ')

disp('effect size')
zzz(maxInd) ./ sqrt(numel(glMT2))

disp('mean Z')
nanmean(zzz(:))

disp('std Z')
nanstd(zzz(:))



% figure, 
histogram((glMT2 - thrMT2),linspace(-2,2,6))

title('goal movetype - threat movetype. -ve is towards stimulus')
%% Show artificial neuron preferences - which prefer goals more, and which prefer threats more

disp([''])
disp(['Proportion of neurons that correlate with goal proximity'])
disp(rMat_Gl.propCorrNeur)
disp(['Av: ' num2str(nanmean(rMat_Gl.propCorrNeur)) '+-' num2str(nanstd(rMat_Gl.propCorrNeur)) ])


disp([''])
disp(['Proportion of neurons that correlate with threat proximity'])
disp(rMat_Thr.propCorrNeur)
disp(['Av: ' num2str(nanmean(rMat_Thr.propCorrNeur)) '+-' num2str(nanstd(rMat_Thr.propCorrNeur)) ])

% Calcualte the proportion of neurons that correlates with both goals and threats
bothHN = (rMat_Gl.hProxN==1 & rMat_Thr.hProxN==1);
% Also the proportion of neurons that correlates with EITHER goals or threats
eitherHN = (rMat_Gl.hProxN==1 & rMat_Thr.hProxN==1) | (rMat_Gl.hProxN==1 | rMat_Thr.hProxN==1);
for iM = 1:length(rS)
% proportion of neurons which correlate with BOTH, by layer and across layers
propBothPerLay(1:length(rS(iM).s.lp.netS),iM) = ...
    sum(bothHN(1:length(rS(iM).s.lp.netS),:,iM),2)./(rS(iM).s.lp.netS)';
propBoth(iM) = sum(bothHN(:,:,iM),'all')./sum(rS(iM).s.lp.netS);

% proportion of neurons which correlate with EITHER, by layer and across layers
propEitherPerLay(1:length(rS(iM).s.lp.netS),iM) = ...
    sum(eitherHN(1:length(rS(iM).s.lp.netS),:,iM),2)./(rS(iM).s.lp.netS)';
propEither(iM) = sum(eitherHN(:,:,iM),'all')./sum(rS(iM).s.lp.netS); 
end

disp([''])
disp(['Proportion of neurons that correlate with BOTH goal AND threat proximity'])
disp(propBoth)
disp(['Av: ' num2str(nanmean(propBoth)) '+-' num2str(nanstd(propBoth)) ])

disp([''])
disp(['Proportion of neurons that correlate with EITHER of goal and threat proximity'])
disp(propEither)
disp(['Av: ' num2str(nanmean(propEither)) '+-' num2str(nanstd(propEither)) ])

warning('I SHOULD ADD IN LAYER DEPTH SO THAT ANY DEPTH CALCULATIONS WORK CORRECTLY!!!')
warning('I SHOULD ADD IN LAYER DEPTH SO THAT ANY DEPTH CALCULATIONS WORK CORRECTLY!!!')
warning('I SHOULD ADD IN LAYER DEPTH SO THAT ANY DEPTH CALCULATIONS WORK CORRECTLY!!!')


%% See whether networks are consistently more ordered in some way

% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\SuccessorState\NetworkAnalysisWorkSpace_15iVV.mat')

nPm = 100;

clear allBinNPFpm allBinNPF allSumBinNPFdist allSumBinNPFpmdist

structMetricToCalc = 'DistSimCorr';
structMetricToCalc = 'GraphColour';

for iM = 1:length(rS) %-1
    
    iM
    rSall(iM) = rS(iM);
    rSall(iM).s.plt.nPm = 100 ;
    [rSall(iM),allNetAR(iM),allSumBinNPFdist(iM),allSumBinNPFpmdist(:,iM),allBinNPF(:,:,iM),allBinNPFpm(:,:,:,iM)] = AssessStructure(rSall(iM));


    
    s=DefaultSettings(rS(iM).s);
    w=rS(iM).w;
    net=rS(iM).net;
    
    sFP = s;
    
    % sFP.nta.comparison = 'BodyPart';
    sFP.nta.comparison = 'Valence';
    
    % 'Absolute' 'Row' 'Column' 'AbsRow' 'AbsColumn'
    % sFP.plt.distanceType = 'AbsColumn';
    sFP.plt.distanceType = 'Absolute';
%     sFP.plt.distanceType = 'Row';
    % sFP.plt.distanceType = 'AbsRow';
    
    sFP.plt.startLayer = 1;
    sFP.plt.stopLayer  = length(s.lp.netS) ;
    
    sFP.plt.ON = 0;
    
    [netAR,s] = PlotConnections(net,sFP,w);
    
    allNetAR(iM) = netAR;

    

    
    
    % ===========================================================
    % Try to do a proximity 'colour' thing
    
    
    figure,
    p = plot(netAR.G,'Layout','force','WeightEffect','inverse','UseGravity','on');
    
    % Make an indicator of layer depth
    tmp = netAR.A_B_rat';
    layNum = [1:size(netAR.A_B_rat,1)] .*    ones(size(tmp));
    layNum = layNum(:);
    layNum(isnan(tmp(:))) = [];
    
    A_B_rat = NanRemFlatten(netAR.A_B_rat'); classType = 'rat';
%     A_B_rat = NanRemFlatten(netAR.A_B_diff'); classType = 'diff';
%     A_B_rat = NanRemFlatten(netAR.AorB'); classType = 'bool';
    G = netAR.G;
    
    
% % %     figure,
    hold on;
    % First sort nodes by distance to all other nodes (i.e. make a distance
    % matrix)
    clear nodeDists closestDists closestNodes
% % %     scatter(p.XData,p.YData,1000,[zeros(size(A_B_rat)) 1-A_B_rat A_B_rat], '.')
    for iNode=1:size(G.Nodes,1)
        
        if strcmp(s.nta.compScale,'Graded')
            %         scatter(p.XData(iNode),p.YData(iNode),1000,[0 1-A_B_rat(iNode) A_B_rat(iNode)], '.')
            %         scatter(p.XData(iNode),p.YData(iNode),1000,[0 .5-4.*A_B_rat(iNode) .5+4.*A_B_rat(iNode)], '.')
        end
        
        nodeDists(:,iNode)    = sqrt(sum(([p.XData(iNode),p.YData(iNode)] - [p.XData(:),p.YData(:)])'.^2));
        [closestDists(:,iNode) closestNodes(:,iNode)] = sort(nodeDists(:,iNode));
    end
    close all

    switch structMetricToCalc
        case 'DistSimCorr'
% % % 
% % %             stimPrefSim  = 1 - abs(A_B_rat - A_B_rat');
% % %             [rho(iM) pval(iM)] = corr(stimPrefSim(:),nodeDists(:));

        case 'GraphColour'
    
    
    tmpX = p.XData;
    tmpY = p.YData;
    
    % % % % [theta rho] = cart2pol(tmpX,tmpY);
    % % %
    % % % % theta = deg2rad(55);
    % % % theta = deg2rad(0);
    % % % % theta = deg2rad(-30);
    % % % % theta = deg2rad(-45);
    % % % % theta = deg2rad(-60);
    % % % R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    % % %
    % % % rotXY = R*[tmpX;tmpY];
    % % %
    % % % tmpX = rotXY(1,:);
    % % % tmpY = rotXY(2,:);
    
    
    
    
    % ==============================================================
    % Try to make it into a 2d plot and convolve with a 2D kernel
    
    % For non-rotated x and y
    nBinsX = 100;
    nBinsY = 100;
    sgm = 5;
    sz  = 15;
    
    
    myGausFilt = images.internal.createGaussianKernel([sgm sgm], [sz sz]);
    
    binXw = (max(tmpX) - min(tmpX))./nBinsX;
    binYw = (max(tmpY) - min(tmpY)) ./nBinsY;
    
    % Binned neural preferences
    binNP  = nan([nBinsX nBinsY]);
    binNPF = nan([nBinsX nBinsY]);
    binCsx = linspace(min(tmpX),max(tmpX),nBinsX);
    binCsy = linspace(min(tmpY),max(tmpY),nBinsY);
    
    for iBx = 1:nBinsX
        inX =   tmpX >= binCsx(iBx) - binXw ./2 & ...
            tmpX <= binCsx(iBx) + binXw ./2 ;
        for iBy = 1:nBinsY
            inY =   tmpY >= binCsy(iBy) - binYw ./2 & ...
                tmpY <= binCsy(iBy) + binYw ./2 ;
            
            %         if sum(inX & inY) > 0
            %             binNP(iBx,iBy) = nanmean(A_B_rat(inX & inY) - .5);
            if strcmp(classType,'diff')
                binNP(iBx,iBy) = nansum(A_B_rat(inX & inY));
            elseif strcmp(classType,'rat') | strcmp(classType,'bool')
                binNP(iBx,iBy) = nansum(A_B_rat(inX & inY) - .5);
            end
            %         end
        end
    end
    
    
    
    %  Pad binNP so that it can be convolved nicely
    tmpPad = nan(size(binNP) + [1 1] .*2 .* ceil(sz/2));
    tmpPad(ceil(sz/2) + 1 : end - ceil(sz/2) , ceil(sz/2) + 1 : end - ceil(sz/2) ) = binNP;
    
    
    % Make my own filtered version with hookers and cocaine
    basePoint = ceil(sz/2) + 1;
    for iBx = 1:nBinsX
        for iBy = 1:nBinsY
            binNPF(iBx,iBy) = nansum(tmpPad(basePoint + [iBx - floor(sz/2) : iBx + floor(sz/2)] , ...
                basePoint + [iBy - floor(sz/2) : iBy + floor(sz/2)]  ) .* ...
                myGausFilt, 'all');
        end
    end
    
    
    % Try finding classification distance to neighbours
    for iBx = 2:nBinsX-1
        for iBy = 2:nBinsY-1
            tmpNeighbours = binNPF(iBx + [-1 1],iBy + [-1 1]);
            binNPFdist(iBx,iBy) = sum((binNPF(iBx,iBy) - tmpNeighbours(:)).^2);
        end
    end
    figure,imagesc(binNPFdist); colorbar;
    sumBinNPFdist = sum(binNPFdist(:))
    
    
    figure,ImagescInvisNan(1,binCsx,binCsy,binNPF')
    axis xy
    colormap(redbluecmapRory(5,6))
    SymColAxis;
    colorbar
    
    
    
    figure,ImagescInvisNan(1,binCsx,binCsy,binNPF')
    axis xy
    colormap(redbluecmapRory(5,6))
    SymColAxis;
    colorbar
    
    colormap jet
    hold on
    
    if strcmp(classType,'diff')
                    scatter(tmpX ,tmpY ,30,[.5+A_B_rat zeros(size(A_B_rat)) .5-A_B_rat], ...
        'o','Filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.0,'MarkerEdgeColor','w')
    scatter(tmpX ,tmpY ,30,[A_B_rat zeros(size(A_B_rat)) 1-A_B_rat], ...
        'o','Filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.0,'MarkerEdgeColor','w')
    end
    
    colormap(redbluecmapRory(100,10))
    end

    switch structMetricToCalc
            case 'DistSimCorr'
            case 'GraphColour'

    % ==============================================================
    % Then perform permutations

    % Binned neural preferences
    binNPpm  = nan([nBinsX nBinsY nPm]);
    binNPFpm = nan([nBinsX nBinsY nPm]);


    for iPm = 1:nPm


        % % %         % $$$ permute across layers
        % % %         randOrd = randperm(numel(tmpX));
        % % %         tmpXpm = tmpX(randOrd);
        % % %         tmpYpm = tmpY(randOrd);

        % $$$ permute within layers $$$ HERE HERE5
        iN = 0;
        %for iL = 1:length(rS(iM).s.lp.netS)
        clear tmpXpm tmpYpm
        for iL = sFP.plt.startLayer: sFP.plt.stopLayer

            nNInLay = sum(~isnan(netAR.A_B_diff(iL,:)));

            iN2 = iN + nNInLay;

            randOrd = iN + randperm(nNInLay);
            tmpXpm(1 + iN:iN2) = tmpX(randOrd);
            tmpYpm(1 + iN:iN2) = tmpY(randOrd);

            iN = iN2;
        end


        for iBx = 1:nBinsX
            inX =   tmpXpm >= binCsx(iBx) - binXw ./2 & ...
                tmpXpm <= binCsx(iBx) + binXw ./2 ;
            for iBy = 1:nBinsY
                inY =   tmpYpm >= binCsy(iBy) - binYw ./2 & ...
                    tmpYpm <= binCsy(iBy) + binYw ./2 ;

                %         if sum(inX & inY) > 0
                if strcmp(classType,'diff')
                    binNPpm(iBx,iBy,iPm) = nansum(A_B_rat(inX & inY));
                elseif strcmp(classType,'rat') | strcmp(classType,'bool')
                    binNPpm(iBx,iBy,iPm) = nansum(A_B_rat(inX & inY) - .5);
                end
                %             binNPpm(iBx,iBy,iPm) = nansum(A_B_rat(inX & inY));
                %         end
            end
        end


        % Pad the data so I can have full binNPFpm
        %  Pad binNP so that it can be convolved nicely
        tmpPad = nan([size(binNP) + [1 1] .*2 .* ceil(sz/2)]);
        tmpPad(ceil(sz/2) + 1 : end - ceil(sz/2) , ceil(sz/2) + 1 : end - ceil(sz/2)) = binNPpm(:,:,iPm);

        % Make my own filtered version with hookers and cocaine
        basePoint = ceil(sz/2) + 1;
        for iBx = 1:nBinsX
            for iBy = 1:nBinsY

                binNPFpm(iBx,iBy,iPm) = nansum(tmpPad(basePoint + [iBx - floor(sz/2) : iBx + floor(sz/2)] , ...
                    basePoint + [iBy - floor(sz/2) : iBy + floor(sz/2)]  ) .* ...
                    myGausFilt, 'all');
            end
        end


        % Try finding classification distance to neighbours
        for iBx = 2:nBinsX-1
            for iBy = 2:nBinsY-1
                tmpNeighbours = binNPFpm(iBx + [-1 1],iBy + [-1 1],iPm);
                binNPFpmdist(iBx,iBy,iPm) = sum((binNPFpm(iBx,iBy,iPm) - tmpNeighbours(:)).^2);
            end
        end
        sumBinNPFpmdist(iPm) = sum(binNPFpmdist(:,:,iPm),'all');


    end

    allSumBinNPFdist(iM)        = sumBinNPFdist;
    allSumBinNPFpmdist(:,iM)    = sumBinNPFpmdist;

    sum(sumBinNPFpmdist < sumBinNPFdist)

    sum(squeeze(sum(abs(binNPFpm),[1 2])) > squeeze(sum(abs(binNPF),[1 2])))

    allBinNPFpm(:,:,:,iM) = binNPFpm;
    allBinNPF(:,:,iM)     = binNPF;


    % ===========================================================
    % See if the preferences change as a function of depth
    for iR = 1:size(binNPF,1)
        for iC = 1:size(binNPF,2)
            pPermSpace(iR,iC) = sum(abs(binNPFpm(iR,iC,:)) > abs(binNPF(iR,iC)) ) ./ size(binNPFpm,3);
        end
        tmpPm = squeeze(nanmean(abs(binNPFpm(iR,:,:)),2))  ;
        tmpOg = nanmean(abs(binNPF(iR,:)),2);
        if isnan(tmpOg)
            pPermLine(iR) = NaN;
        else
            pPermLine(iR) = sum( tmpPm > tmpOg ) ./ size(binNPFpm,3);
        end
    end


    for iC = 1:size(binNPF,2)
        tmpPm = squeeze(nanmean(abs(binNPFpm(:,iC,:)),1))  ;
        tmpOg = nanmean(abs(binNPF(:,iC)),1);
        if isnan(tmpOg)
            pPermLine2(iC) = NaN;
        else
            pPermLine2(iC) = sum( tmpPm > tmpOg ) ./ size(binNPFpm,3);
        end
    end



    figure,
    subplot(2,1,1)
    plot(nanmean(abs(binNPF),2))
    hold on, plot(nanmean(abs(binNPFpm),[2 3]))
    xlim([0 200])
    subplot(2,1,2)
    [pthr,pcor,padj] = fdr(pPermLine);
    plot(pPermLine); hold on
    plot(pcor,'x')
    ylim([0 0.1])
    xlim([0 200])
    title('pval row')


    figure,
    subplot(2,1,1)
    plot(nanmean(abs(binNPF),1))
    hold on, plot(nanmean(abs(binNPFpm),[1 3]))
    xlim([0 200])
    subplot(2,1,2)
    [pthr,pcor,padj] = fdr(pPermLine2);
    plot(pPermLine2); hold on
    plot(pcor,'x')
    ylim([0 0.1])
    xlim([0 200])
    title('pval column')

    figure,
    subplot(3,1,1)
    imagesc(abs(binNPF))
    subplot(3,1,2)
    imagesc(nanmean(abs(binNPFpm),[ 3]))
    subplot(3,1,3)
    imagesc(pPermSpace)
    title('p-val')
    caxis([0 0.05])

    end

end

%% Assess the network strength stats

% netAR = allNetAR(1)

netAR.AorB

% for small network
% tmpP = arrayfun(@(x) allNetAR(x).difdifP ,[1:length(rS)]); tmpP(1:3) = [];
% tmpP = arrayfun(@(x) allNetAR(x).difdifP ,[1:48]); tmpP(1:3) = [];

% tmpR = arrayfun(@(x) allNetAR(x).difdifR ,[1:48]); tmpR(1:3) = [];

% for large network
tmpP = arrayfun(@(x) allNetAR(x).difdifP ,[1:length(rS)-1])
tmpR = arrayfun(@(x) allNetAR(x).difdifR ,[1:48]); tmpR(1:3) = [];


% A = arrayfun(@(x) allNetAR(x).classdistP ,[1:length(rS)-1])
% A = arrayfun(@(x) allNetAR(x).difdistP ,[1:length(rS)-1])
% tmpP = arrayfun(@(x) allNetAR(x).ratratP ,[1:length(rS)-1])



[dmy corrP] = fdr(tmpP)

disp('mean rho')
nanmean(tmpR)
disp('std rho')
nanstd(tmpR)




disp('--------------')
sprintf('%i out of %i networks show structure using this metric', sum(corrP < .05), numel(corrP))


histogram(tmpP,20)


structMetric    = squeeze(sum(abs(allBinNPF(:,:,:)),[1 2]));
structMetricPm  = squeeze(sum(abs(allBinNPFpm(:,:,:,:)),[1 2]));

%% Test whether there is more structure than by chance

clear numSmaller
for iM = 1:size(allBinNPF,3)
    numSmaller(iM) = sum(structMetricPm(:,iM) > structMetric(iM))
end


structureOverChance = (structMetric - structMetricPm'   );

figure,histogram(structureOverChance(:))


[hh pp] = ttest(structureOverChance(:))

[hhh ppp ci statsss] =ttest2(structMetric,structMetricPm(:))
disp('struct metric mean')
nanmean(structMetric)

disp('permuted struct metric mean')
nanmean(structMetricPm(:))


figure,histogram(structMetric,10,'Normalization','probability'); hold on
histogram(structMetricPm(:),10,'Normalization','probability');

%% $$$ Do a proper permutation test

tmpD = structMetric - structMetricPm';

sum(tmpD >= 0, 'all') ./ numel(tmpD)


%% Test whether the structure increases as a function of network type

clear structMetricByNetType structMetricByNetTypePm

netShapes = unique(allNetW);
nNS = length(netShapes)
for iNS = 1:nNS
    
    inclM = allNetW == netShapes(iNS);
    
    tmpD    = structMetric(inclM);
    tmpDPm  = structMetricPm(:,inclM);

    structMetricByNetType(iNS,:)   = tmpD;
    structMetricByNetTypePm(iNS,:,:) = tmpDPm;

    disp('')
    disp(['Final Layer size: ' num2str(netShapes(iNS)) ]);
    [hhh ppp ci statsss] =ttest2(tmpD,tmpDPm(:));
    statsss

    sum(tmpD - tmpDPm' >= 0 ,'all') ./ numel(tmpD)
    
end


%% Try out a different metric of order: correlation between neural distance and neural response properties
load('Results\ForFigures\SuccessorState\NetworkAnalysisWorkSpace3.mat');

tic
for iM = 1:length(allNetAR)

    netAR = allNetAR(iM);

    % Make an indicator of layer depth
    tmp = netAR.A_B_rat';
    layNumTmp = [1:size(netAR.A_B_rat,1)] .*    ones(size(tmp));
    layNumTmp = layNumTmp(:);
    layNumTmp(isnan(tmp(:))) = [];
    layNum(:,iM) = layNumTmp;
    layNumDiff(:,:,iM)  = abs(bsxfun(@minus, layNumTmp', layNumTmp));

    A_B_rat = NanRemFlatten(netAR.A_B_rat'); classType = 'rat';
    G = netAR.G;

    p = plot(G,'Layout','force','WeightEffect','inverse','UseGravity','on');

    % First make a distance matrix
    for iNode=1:size(G.Nodes,1)
        nodeDists(:,iNode,iM)    = sqrt(sum(([p.XData(iNode),p.YData(iNode)] - [p.XData(:),p.YData(:)])'.^2));
    end
    % Also permute the distance matrix
    for iPm = 1:100
        permutedNodes = randperm(numel(p.XData));
        for iNode=1:size(G.Nodes,1)
            nodeDistsPm(:,iNode,iM,iPm)    = sqrt(sum(([p.XData(permutedNodes(iNode)),p.YData(permutedNodes(iNode))] - [p.XData(:),p.YData(:)])'.^2));
        end
    end

    % Then make a stimulus preference matrix
    stimPrefSimTmp = 1 - abs(A_B_rat - A_B_rat');
    stimPrefSim(:,:,iM)  = stimPrefSimTmp;

    % Calculate the correlation between the two
    nodeDistsTmp = nodeDists(:,:,iM);
    [rho(iM) pval(iM)] = corr(stimPrefSimTmp(nodeDistsTmp ~= 0),nodeDistsTmp(nodeDistsTmp ~= 0));
%     [rho(iM) pval(iM)] = corr(stimPrefSimTmp(:),nodeDistsTmp(:));
    % $$$ I need to account for stimulus proximity being 0

end
toc

[pValThr, pValcor, pValAdj] = fdr(pval);


fprintf('%i out of %i networks show structure using this metric\n', sum(pValAdj(4:48) < 0.05), 45);

disp('mean rho')
nanmean(rho)
disp('std rho')
nanstd(rho)

%% Perform LME for all models independently, accounting for layer differences

modelNum = repmat( permute((1:numel(allNetAR)), [1 3 2]), [size(stimPrefSim,[1 2]) 1]);
layNumOrg                       = repmat(permute(layNum,[1 3 2]) , [1 size(layNum,1) 1]);

allPerf     = arrayfun( @(im) sum(rS(im).perf.rewPerAct(end,:)) , [1:length(rS)-1] );
allPerfRep  = repmat(  permute(allPerf,[1 3 2]),[size(layNumOrg,[1 2])  1  ] );
allNetW     = arrayfun( @(im) rS(im).s.lp.netS(end) , [1:length(rS)-1] );
allNetWRep  = repmat(  permute(allNetW,[1 3 2]),[size(layNumOrg,[1 2])  1  ] );

allegedNoStructure      = pValAdj > 0.05;
allegedNoStructureRep   = repmat(  permute(allegedNoStructure,[1 3 2]),[size(layNumOrg,[1 2])  1  ] );

tbl                             = table(stimPrefSim(:));
tbl.Properties.VariableNames    = {'StimulusPreferenceSimilarity'};
% % % tbl.netOutputLayerWidth         = allNetW'; 
tbl.nodeDists                   = nodeDists(:);

% tbl.layNumDiff                  = abs(layNumDiff(:));
% tbl.layNumDiff                  = layNumDiff(:);
tbl.layNumDiff                  = categorical(abs(layNumDiff(:)));

tbl.layNumOrg                   = layNumOrg(:);
tbl.modelNum                    = modelNum(:);

tbl.allPerf                     = allPerfRep(:);
tbl.allNetW                     = categorical(allNetWRep(:));

tbl.allegedNoStructure          = allegedNoStructureRep(:);


% tmpTbl = tbl(tbl.nodeDists ~= 0, :);
% tmpTbl = tbl(tbl.nodeDists ~= 0 & tbl.layNumDiff == 0,:);
% tmpTbl = tbl(tbl.nodeDists ~= 0 & tbl.layNumDiff == 0 & tbl.layNumOrg  == 4,:);
% tmpTbl = tbl(tbl.nodeDists ~= 0 & tbl.layNumOrg  == 4,:);

% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists*layNumDiff + (1|modelNum)')
% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists + layNumDiff + (1|modelNum)')
% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists + (1|modelNum)')
% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists + allNetW + allPerf + (1|modelNum)')

% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists + nodeDists:allNetW  + (1|modelNum)')

% anova(lme)

% Run LME for each network
tic
for iM = 1:length(allNetAR)
    tmpTbl          = tbl(tbl.nodeDists ~= 0 & tbl.modelNum == iM,:);
%     tmpLme          = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists*layNumDiff + (1|modelNum)');
    tmpLme          = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ layNumDiff:nodeDists + nodeDists + (1|modelNum)');


% % %     tmpTbl          = tbl(tbl.nodeDists ~= 0 & tbl.modelNum == iM & tbl.layNumDiff == 0,:);
% % %     tmpLme          = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists + (1|modelNum)');
    tmpAnova                = anova(tmpLme);
    estimateDistEff(iM)     = tmpLme.Coefficients{2,2};
    estimateStandErr(iM)    = tmpLme.Coefficients{2,3};
    pNodeDists(iM)          = tmpLme.Coefficients{2,6};

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
toc

[pValThr, pValcor, pValLMEAdj] = fdr(pNodeDists(4:48));

fprintf('%i out of %i networks show structure using this metric\n', sum(pValLMEAdj < 0.05), 45);


disp('Weighted mean estimate of node distance')
sum(  (estimateDistEff ./ estimateStandErr.^2)) ./ sum(1 ./ estimateStandErr.^2)
disp('standard error of effect of node distance')
sqrt(1 ./ sum(1 ./ estimateStandErr.^2))


%%

tmpTbl = tbl(tbl.nodeDists ~= 0, :);

% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists*layNumDiff + nodeDists*allPerf +  (1|modelNum)')
% lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists*allPerf*layNumDiff +  (1|modelNum)')
lme = fitlme(tmpTbl,'StimulusPreferenceSimilarity ~ nodeDists*allPerf +  (1|modelNum)')

anova(lme)

%% See if the performance is related to numsmaller

% % for small network
% allPerf = arrayfun( @(im) sum(rS(im).perf.rewPerAct(end,:)) , [1:length(rS)] );
% allNetW = arrayfun( @(im) rS(im).s.lp.netS(end) , [1:length(rS)] );

% for large network
allPerf = arrayfun( @(im) sum(rS(im).perf.rewPerAct(end,:)) , [1:length(rS)-1] );
allNetW = arrayfun( @(im) rS(im).s.lp.netS(end) , [1:length(rS)-1] );
rS

figure,plot(allPerf,numSmaller,'x'); lsline
figure,histogram2(allPerf,numSmaller,'FaceColor','flat')
[rho pval] = corr(allPerf',numSmaller')



figure,plot(allNetW,numSmaller,'x'); lsline
figure,histogram2(allNetW,numSmaller,'FaceColor','flat')
[rho pval] = corr(allNetW',numSmaller')



disp('Structure over chance correlates with network architecture:')
allNetW2 = repmat(allNetW',[1 size(structureOverChance,2)]);
[rho pval] = corr(allNetW2(:),structureOverChance(:))
figure,plot(allNetW2,structureOverChance,'x'); lsline
figure,histogram2(allNetW2,structureOverChance,'FaceColor','flat')

% $$$ CHECK WHETHER THE BIG HISTOGRAM STILL WORKS WHEN USING THE MORE
% RECENTLY TRAINED NETWORKS SECTION ABOVE)

%%

netShapes = unique(allNetW);
nNS = length(netShapes)
f1 = figure,
f2 = figure,
f3 = figure,
for iNS = 1:nNS
    
    inclM = allNetW == netShapes(iNS);
    
    tmpPM = structMetricPm(:,inclM);
    
    figure(f1)
    subplot(nNS,1,iNS)
    histogram(structMetric(inclM),5,'Normalization','probability'); hold on
    histogram(tmpPM,5,'Normalization','probability');
    title(sprintf('Final layer size: %i' , netShapes(iNS)))
    
    figure(f2)
    subplot(nNS,1,iNS)
    nanmedian(structureOverChance(inclM,:),'all')
    tmpPM = structureOverChance(inclM,:);
%     histogram(tmpPM(:),linspace(-4,8,7),'Normalization','probability'); hold on
    histogram(tmpPM(:),'Normalization','probability'); hold on
    tmpPM2(1:numel(tmpPM(:)),iNS) = tmpPM(:);
    title(sprintf('Final layer size: %i' , netShapes(iNS)))    
    xlim([-4 8])
    
    figure(f3)
%     subplot(nNS,1,iNS)
    tmpPM = structureOverChance(inclM,:);
%     histogram(tmpPM(:),linspace(-4,8,7),'Normalization','probability'); hold on
    histogram(tmpPM(:),'Normalization','probability'); hold on
    tmpPM2(1:numel(tmpPM(:)),iNS) = tmpPM(:);
    title('Functional sub-networks naturally arise, and are more likely to arise when information is condensed')    
    xlabel('Increase in network order over chance')
    xlim([-4 8])

end

legend('Narrowing networks', 'constant width networks', 'Widening networks')

% figure,
% hist(tmpPM2,20)

%%

rS(6).s.lp.netS

















































%%
%%
%%







%% Check connectivity for all networks (threats and goal networks)

for iM = 1:length(rS)
    
      iM
    
    s=rS(iM).s;
    w=rS(iM).w;
    net=rS(iM).net;
    
    sFP = s;
    
    % sFP.nta.comparison = 'BodyPart';
    sFP.nta.comparison = 'Valence';
    
    % 'Absolute' 'Row' 'Column' 'AbsRow' 'AbsColumn'
    % sFP.plt.distanceType = 'AbsColumn';
    sFP.plt.distanceType = 'Absolute';
%     sFP.plt.distanceType = 'Row';
    % sFP.plt.distanceType = 'AbsRow';
    
    sFP.plt.startLayer = 1;
    sFP.plt.stopLayer  = length(s.lp.netS) ;
    
    sFP.plt.ON = 0;
    
    [netAR,s] = PlotConnections(net,sFP,w);
    
    allNetAR(iM) = netAR;
    
    close all
    
end

%% Show separability of threat and goal networks

for iM = 1:length(rS)
    pDistCorr(iM) = allNetAR(iM).difdifP;
    rDistCorr(iM) = allNetAR(iM).difdifR;
    
    pClassDist(iM) = allNetAR(iM).classdistP;
    zClassDist(iM) = allNetAR(iM).classdistZ;
end

[pDistCorr] = CalcFDR(pDistCorr);
[pClassDist] = CalcFDR(pClassDist);

[pDMax pDMaxInd] = max(pDistCorr);
[pCMax pCMaxInd] = max(pClassDist);

disp('Correlation between (covariance difference) of neuron and previous neurons')
disp(['P <= ' num2str(pDMax) ', rho >= ' num2str(rDistCorr(pDMaxInd))])

disp('Ranksum of distances to all other goal and threat neurons, between goal and threat neurons')
disp(['P <= ' num2str(pCMax) ', Z >= ' num2str(zClassDist(pCMaxInd))])


%% Try to alter the network (most of it from TestZone)

% Check number of 'towards/away' moves, for different bias amounts
% $$$ HERE HERE --> MAYBE try doing it for each neuron in the last layer or
% 2 independently?

% $$$ OR EVEN BETTER: change them in the direction of the correlation! $$$

% $$$ Maybe use only the neurons which are 'goal centric' from previous analyses?

% set to 1 to alter full net. set to 0 to alter single neurons
alterNetFL = 1;

for iM = 1:length(rS) % $$$ WATCH OUT HERE
    
    rS(iM).allMT_GL = [];
    rS(iM).allMT_THR = [];
    rS(iM).allQbyType_GL = [];
    rS(iM).allQbyType_THR = [];

s = DefaultSettings(rS(iM).s);
w = rS(iM).w;
net = rS(iM).net;
Qtable = rS(iM).Qtable;



aS.ListOrMat='List';
% aS.alterType='Bias';

% =========================================================================
% % % biasVal=linspace(-0.02,0.02,9);
% % % aS.alterType='DirectionalConnectionBias';


% biasVal=linspace(-0.01,0.01,9);
% biasVal=linspace(-0.05,0.05,21);
biasVal=linspace(-0.10,0.10,21);
% biasVal=linspace(-0.40,0.40,21);
% biasVal=linspace(-0.20,0.20,21);
% biasVal=[linspace(-0.15,-0.105,10) , 0 ,  linspace(0.105,0.15,10)] ;
% biasVal=linspace(-0.5,0.5,101);
% biasVal=linspace(-0.5,0.5,21);
% biasVal=linspace(-0.25,0.25,41);
% % biasVal=linspace(-0.2,0.2,51);
aS.alterType='MultiplicativeDirectionalConnectionBias';


% This is basically ablating (-1) and not ablating (0)
% biasVal = [-1 0];
% aS.alterType = 'MultiplicativeConnectionBias';
% $$$ LAST THING TO CHECK IS WHETHER ABLATING WORKS BETTER
% $$$ LAST THING TO CHECK IS WHETHER ABLATING WORKS BETTER
% $$$ LAST THING TO CHECK IS WHETHER ABLATING WORKS BETTER
% $$$ LAST THING TO CHECK IS WHETHER ABLATING WORKS BETTER
% $$$ LAST THING TO CHECK IS WHETHER ABLATING WORKS BETTER
% $$$ LAST THING TO CHECK IS WHETHER ABLATING WORKS BETTER
% $$$ LAST THING TO CHECK IS WHETHER ABLATING WORKS BETTER
% $$$ LAST THING TO CHECK IS WHETHER ABLATING WORKS BETTER
% =========================================================================

s.plt.ON=0;
% s.plt.ON=1;

s.plt.colLag=0;

s.plt.stimCol=2:s.wrld.size(2)-1;
s.plt.lmbCol=2:s.wrld.size(2)-1;
s.plt.meanLimbCols = 1;

% aS.SignedBias=1;
aS.SignedBias=0;

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% aS.ignoreLay=[1:floor(length(rS(iM).s.lp.netS)/2)];
% aS.ignoreLay=[1:(length(rS(iM).s.lp.netS)-2)];
% aS.ignoreLay=1:2;
aS.ignoreLay=[];

s.plt.OneActFl=0;
s.plt.pltType = 'Imagesc';


% settings for sigmoid fitting
testLC = [2:s.wrld.size(2)-1];
s2 = DefaultSettings(s);
s2.plt.stimRow = 2:s2.wrld.size(1)-2;
s2.plt.stimCol = [2:s2.wrld.size(2)-1];
s2.plt.fitSigmoid = 1;

clear allMT


% Goal and threat stimulus
stimTypes=[0 1];
% Goal and threat network (stupidly set the other way round...)
netTypes=[1 0];



% =============================================================
% [aS.lay aS.neur]=find(netAR.AorB==netTypes(iNT));
% aS.lay=9; aS.neur=5;

% [glL glN] = ... 
%     find(rMat_Gl.hProxN(:,:,iM) == 1);
% [thrL thrN] = ...
%     find(rMat_Thr.hProxN(:,:,iM) == 1);


% % % $$$ RECENT
% [glL glN] = ...
%     find(allNetAR(iM).AorB == 1 );

[thrL thrN] = ...
    find(allNetAR(iM).AorB == 0 );



% % [glL glN] = ...
% %     find(rMat_Gl.hProxN(:,:,iM) == 1 & ~(rMat_Thr.hProxN(:,:,iM) == 1) );
% % 
% % [thrL thrN] = ...
% %     find(rMat_Thr.hProxN(:,:,iM) == 1 & ~(rMat_Gl.hProxN(:,:,iM) == 1) );


% nN = size(allNetAR(iM).AorB,2);
% [glL glN] = ...
%     find(rMat_Gl.hProxN(:,1:nN,iM) == 1 & (allNetAR(iM).AorB==1) );
% [thrL thrN] = ...
%     find(rMat_Thr.hProxN(:,1:nN,iM) == 1 & (allNetAR(iM).AorB==0) );



nN = size(allNetAR(iM).AorB,2);
[glL glN] = ...
    find(rMat_Gl.hProxN(:,1:nN,iM) == 1 & ~(rMat_Thr.hProxN(:,1:nN,iM) == 1) & (allNetAR(iM).AorB==1) );
[thrL thrN] = ...
    find(rMat_Thr.hProxN(:,1:nN,iM) == 1 & ~(rMat_Gl.hProxN(:,1:nN,iM) == 1) & (allNetAR(iM).AorB==0) );

% =============================================================

for iB = 1:numel(biasVal) % bias
    
    % $$$$$ HERE --> maybe SET the biasvals' signs depending on the sign of the correlation?
    aS.biasVals=biasVal(iB);
    
    nStimLoops=length(stimTypes);
    
    nNetLoops=length(netTypes);
    
    
    for iST=1:nStimLoops
        
        % $$$ HERE HERE - bear in mind I have to probably make a separate varibale
        % for Thr and Gl Qs and MovementTypes. Also need to remove the '8' limbcol,
        % and make sure that the extreme limbcols aren't included. Then the idea is
        % to calculate those values for the alteration of every neuron in every net
        % type
        
        iNT = 1;
        if alterNetFL == 1
            nNL = 1;
        else
            nNL = length(glL);
        end
        for iNL = 1:nNL
            if alterNetFL == 1
                aS.lay = glL;
                aS.neur = glN;
            else
                aS.lay = glL(iNL);
                aS.neur = glN(iNL);
            end
            % Posibly change the sign of the correlation
            if aS.SignedBias==1
                clear aS.biasVals
                for iN=1:numel(aS.neur)
                    % netAR.rThr and rGl are (confusingly..) positive if the neuron's
                    % firing rate becomes more NEGATIVE with proximity
                    % (because it is the correlation with distance)
                    if netTypes(iNT)==1
                        aS.biasVals(iN)=-sign(netAR.rGl(aS.lay(iN),aS.neur(iN))).*biasVal(iB);
                    elseif netTypes(iNT)==0
                        aS.biasVals(iN)=-sign(netAR.rThr(aS.lay(iN),aS.neur(iN))).*biasVal(iB);
                    end
                    
                end
            end
            s.plt.plotThrFl=stimTypes(iST);
            if iM == 3
                disp('test')
            end
            [moveType_2 squishQs QbyType tmpQ]=FindMoveTypes(net,s,w,aS);
            rS(iM).allQbyType_GL(:,:,iB,iNL,iST,:)=QbyType(:,[2:7 9:s.wrld.size(2)-1],:);
            rS(iM).allMT_GL(:,:,iB,iNL,iST)=moveType_2(:,[2:7 9:s.wrld.size(2)-1]);
            
            
            % ============================================================
            % work in progress
            s2.plt.lmbCol = 8;
            moveType_2b = nan(size(tmpQ(:,:,:,:,1)));
            moveType_2b(size(w.world2D,1)-2,s2.plt.lmbCol,:,:) = permute(moveType_2,[3 4 1 2]);
            moveType_2c = moveType_2b;
            moveType_2c(moveType_2b == 0) = 1;
            moveType_2c(moveType_2b ~= 0) = 0;
            % Only look at zone with stimulus-relevant movetypes (i.e. goal-actions)
            if iST == 1
                moveType_2b(moveType_2b >= 0) = 0;
                moveType_2b(moveType_2b >= 0) = 1;
                moveType_2b(moveType_2b < 0) = 0;
            elseif iST == 2
                moveType_2b(moveType_2b <= 0) = -1;
                moveType_2b(moveType_2b > 0) = 0;
                moveType_2b(moveType_2b < 0) = 1;
            end
            
            s2.plt.fitSigmDistTypes = {'aD'};
            s2.plt.fitSigmStart = [1, .5, 0, 1]; 
            [~,~,~,~,~,~,sigmOut,sigMidQ] = ...
                CalcDistCorr(s2,w,double(moveType_2b));
            rS(iM).sigMidQ_Gl_aD(iB,iNL,iST) = sigMidQ.aD;
            rS(iM).sigmOut_Gl(iB,iNL,iST) = sigmOut;
            
            [~,~,~,~,~,~,sigmOut,sigMidQ] = ...
                CalcDistCorr(s2,w,double(moveType_2c));
            rS(iM).sigMidQ_Gl_aD_Anymove(iB,iNL,iST) = sigMidQ.aD;
            rS(iM).sigmOut_Gl_Anymove(iB,iNL,iST) = sigmOut;
            % ============================================================
            
        end
        
        
        
        iNT = 2;
        if alterNetFL == 1
            nNL = 1;
        else
            nNL = length(thrL);
        end
        for iNL = 1:nNL
            if alterNetFL == 1
                aS.lay = thrL;
                aS.neur = thrN;
            else
                aS.lay = thrL(iNL);
                aS.neur = thrN(iNL);
            end
            
            % Posibly change the sign of the correlation
            if aS.SignedBias==1
                clear aS.biasVals
                for iN=1:numel(aS.neur)
                    % netAR.rThr and rGl are (confusingly..) positive if the neuron's
                    % firing rate becomes more NEGATIVE with proximity
                    % (because it is the correlation with distance)
                    if netTypes(iNT)==1
                        aS.biasVals(iN)=-sign(netAR.rGl(aS.lay(iN),aS.neur(iN))).*biasVal(iB);
                    elseif netTypes(iNT)==0
                        aS.biasVals(iN)=-sign(netAR.rThr(aS.lay(iN),aS.neur(iN))).*biasVal(iB);
                    end
                    
                end
            end
            s.plt.plotThrFl=stimTypes(iST);
            [moveType_2 squishQs QbyType tmpQ]=FindMoveTypes(net,s,w,aS);
            rS(iM).allQbyType_THR(:,:,iB,iNL,iST,:)=QbyType(:,[2:7 9:s.wrld.size(2)-1],:);
            rS(iM).allMT_THR(:,:,iB,iNL,iST)=moveType_2(:,[2:7 9:s.wrld.size(2)-1]);
            
            
            % ============================================================
            % work in progress
            s2.plt.lmbCol = 8;
            moveType_2b = nan(size(tmpQ(:,:,:,:,1)));
            moveType_2b(size(w.world2D,1)-2,s2.plt.lmbCol,:,:) = permute(moveType_2,[3 4 1 2]);
            moveType_2c = moveType_2b;
            moveType_2c(moveType_2b == 0) = 1;
            moveType_2c(moveType_2b ~= 0) = 0;% Only look at zone with stimulus-relevant movetypes (i.e. goal-actions)
            if iST == 1
                moveType_2b(moveType_2b >= 0) = 0;
                moveType_2b(moveType_2b >= 0) = 1;
                moveType_2b(moveType_2b < 0) = 0;
            elseif iST == 2
                moveType_2b(moveType_2b <= 0) = -1;
                moveType_2b(moveType_2b > 0) = 0;
                moveType_2b(moveType_2b < 0) = 1;
            end
            
            s2.plt.fitSigmDistTypes = {'aD'};
            s2.plt.fitSigmStart = [1, .5, 0, 1]; 
            [~,~,~,~,~,~,sigmOut,sigMidQ] = ...
                CalcDistCorr(s2,w,double(moveType_2b));
            rS(iM).sigMidQ_Thr_aD(iB,iNL,iST) = sigMidQ.aD;
            rS(iM).sigmOut_Thr(iB,iNL,iST) = sigmOut;
            
            [~,~,~,~,~,~,sigmOut,sigMidQ] = ...
                CalcDistCorr(s2,w,double(moveType_2c));
            rS(iM).sigMidQ_Thr_aD_Anymove(iB,iNL,iST) = sigMidQ.aD;
            rS(iM).sigmOut_Thr_Anymove(iB,iNL,iST) = sigmOut;
            % ============================================================
            
            
        end
    end
end
iM
end

%% Try the expansiveness thing  HERE
clear sigMidQ_Gl_aD sigMidQ_Thr_aD
for iM = 1:length(rS)
    
%     sigMidQ_Gl_aD(:,:,:,iM) = rS(iM).sigMidQ_Gl_aD;
%     sigMidQ_Thr_aD(:,:,:,iM) = rS(iM).sigMidQ_Thr_aD;

    sigMidQ_Gl_aD(:,:,:,iM) = rS(iM).sigMidQ_Gl_aD_Anymove;
    sigMidQ_Thr_aD(:,:,:,iM) = rS(iM).sigMidQ_Thr_aD_Anymove;
    
end

figure,

inclM = pDistCorr<=1 ;
inclM = pDistCorr<=0.05 ;
% inclM = pDistCorr<=1 & (1:length(rS))<=44;
% inclM = pDistCorr<=1 & (1:length(rS))>=45;
% inclM = pDistCorr<=0.05 ;
% inclM = pDistCorr<=0.05 & (1:length(rS))>=39;

subplot(2,2,1)
plot(biasVal,squeeze(nanmedian(sigMidQ_Gl_aD(:,1,1,inclM ),4)),'-o')
title('goal net,goal stim')

subplot(2,2,2)
plot(biasVal,squeeze(nanmedian(sigMidQ_Gl_aD(:,1,2,inclM ),4)),'-o')
title('goal net,thr stim')

subplot(2,2,3)
plot(biasVal,squeeze(nanmedian(sigMidQ_Thr_aD(:,1,1,inclM ),4)),'-o')
title('threat net,goal stim')

subplot(2,2,4)
plot(biasVal,squeeze(nanmedian(sigMidQ_Thr_aD(:,1,2,inclM ),4)),'-o')
title('threat net,thr stim')


%% Do stats on the best action choice sigmoid distance
% $$$ taken from above

% % % for iM = 1:1%length(rS)
% % %     for iLC = 1:length(testLC)
% % %         sMG(iLC,iM) = sigMidQ_Gl(iLC,iM).aD;
% % %         sMT(iLC,iM) = sigMidQ_Thr(iLC,iM).aD;
% % %     end
% % % end
% % % 
% % % valenceVals = [ones(size(sMT)).*sFP.act.GoalRew; ones(size(sMT)).*sFP.act.ThreatRew];
% % % 
% % % % opts.testType = 'lme';
% % % % [figString zVal pVal] = CorrelationTitle(valenceVals,[sMG; sMT],opts)
% % % 
% % % figure,
% % % plot(valenceVals,[sMG; sMT],'x')
% % % 
% % % opts.testType = 'signrank';
% % % [figString zVal pVal] = CorrelationTitle([],[sMG(:), sMT(:)],opts)
% % % 
% % % hist([sMG(:), sMT(:)])
% % % legend('2 rew','4 rew')
% % % title(figString)
% % % xlabel('distance at which staying becomes the preferred movement')







%% Plot things for bias adjustment
clear bvG bvT
for iM = 1:length(rS)

% figure,
for iST=1:length(stimTypes)
    subplot(2,length(stimTypes),iST)
    
    for iN=1:size(rS(iM).allMT_GL,4)
        plot(biasVal,squeeze(sum(rS(iM).allMT_GL(:,:,1:length(biasVal),iN,iST),[1 2])),'-o','LineWidth',1); hold on
    end
    if iST == 1
        title('GOALneur change. Response to GOALS')
    elseif iST == 2
        title('GOALneur change. Response to THREATS')
    end
    
    plot(biasVal,squeeze(nanmean(sum(rS(iM).allMT_GL(:,:,1:length(biasVal),:,iST),[1 2]),4)),'-ok','LineWidth',1); hold on
%     plot(biasVal,squeeze(nanmedian(sum(rS(iM).allMT_GL(:,:,1:length(biasVal),:,iST),[1 2]),4)),'-ok','LineWidth',2); hold on
    xlabel('bias amount')
    ylabel('Defensiveness of actions')
    
    bvG(:,iST,iM) = squeeze(nanmean(sum(rS(iM).allMT_GL(:,:,1:length(biasVal),:,iST),[1 2]),4));
end

for iST=1:length(stimTypes)
    subplot(2,length(stimTypes),iST+2)
    
    for iN=1:size(rS(iM).allMT_THR,4)
        plot(biasVal,squeeze(sum(rS(iM).allMT_THR(:,:,1:length(biasVal),iN,iST),[1 2])),'-o','LineWidth',1); hold on
    end
    if iST == 1
        title('THREATneur change. Response to GOALS')
    elseif iST == 2
        title('THREATneur change. Response to THREATS')
    end
    
    plot(biasVal,squeeze(nanmean(sum(rS(iM).allMT_THR(:,:,1:length(biasVal),:,iST),[1 2]),4)),'-ok','LineWidth',1); hold on
%     plot(biasVal,squeeze(nanmedian(sum(rS(iM).allMT_THR(:,:,1:length(biasVal),:,iST),[1 2]),4)),'-ok','LineWidth',2); hold on
    xlabel('bias amount')
    ylabel('Defensiveness of actions')
    
    bvT(:,iST,iM) = squeeze(nanmean(sum(rS(iM).allMT_THR(:,:,1:length(biasVal),:,iST),[1 2]),4));
end

end


for iST=1:length(stimTypes)
    subplot(2,length(stimTypes),iST)
    plot(biasVal,nanmedian(bvG(:,iST,:),3),'-or','LineWidth',2);
    
    subplot(2,length(stimTypes),iST+2)
    plot(biasVal,nanmedian(bvT(:,iST,:),3),'-or','LineWidth',2);
end

%% $$$ MAYBE 1 more thing to try is use ALLLL the neurons from all models?..

allMT = [];
allNT = [];
allM = [];
allST = [];
allB = [];
allSEPpDC = [];
allSEPpCD = [];

for iM = 1:length(rS)

% figure,
for iST=1:length(stimTypes)
    subplot(1,length(stimTypes),iST)
    
    
   
% %         boxplot(squeeze( sum(rS(iM).allMT_GL(:,:,1,:,iST),[1 2]))) 
% %         boxplot(squeeze( sum(rS(iM).allMT_THR(:,:,1,:,iST),[1 2]))); hold on

    % Make a big table to store all the single neuron effects
    for iB =1:size(rS(iM).allMT_GL,3)
        tmpMT = squeeze( sum(rS(iM).allMT_GL(:,:,iB,:,iST),[1 2]));
        tmpNT = ones(size(tmpMT)); % net type
        tmpM = ones(size(tmpMT)).*iM;
        tmpST = ones(size(tmpMT)).*iST;
        tmpB = ones(size(tmpMT)).*biasVal(iB);
        tmpSEPpDC = ones(size(tmpMT)).*pDistCorr(iM);
        tmpSEPpCD = ones(size(tmpMT)).*pClassDist(iM);
        allMT = [allMT ; tmpMT]; 
        allNT= [allNT ; tmpNT];
        allM = [allM ; tmpM];
        allST = [allST; tmpST];
        allB = [allB; tmpB];
        allSEPpDC = [allSEPpDC; tmpSEPpDC];
        allSEPpCD = [allSEPpCD; tmpSEPpCD];

        tmpMT = squeeze( sum(rS(iM).allMT_THR(:,:,iB,:,iST),[1 2]));
        tmpNT = -ones(size(tmpMT));
        tmpM = ones(size(tmpMT)).*iM;
        tmpST = ones(size(tmpMT)).*iST;
        tmpB = ones(size(tmpMT)).*biasVal(iB);
        tmpSEPpDC = ones(size(tmpMT)).*pDistCorr(iM);
        tmpSEPpCD = ones(size(tmpMT)).*pClassDist(iM);
        allMT = [allMT ; tmpMT]; 
        allNT = [allNT ; tmpNT];
        allM = [allM ; tmpM];
        allST = [allST; tmpST];
        allB = [allB; tmpB];
        allSEPpDC = [allSEPpDC; tmpSEPpDC];
        allSEPpCD = [allSEPpCD; tmpSEPpCD];
        

    end
    
    
    plot(biasVal, ...
        squeeze(nanmean(sum(rS(iM).allMT_GL(:,:,1:length(biasVal),:,iST),[1 2]),4)) - ...
        squeeze(nanmean(sum(rS(iM).allMT_THR(:,:,1:length(biasVal),:,iST),[1 2]),4)),'-ok','LineWidth',2); hold on
%     plot(biasVal,squeeze(nanmedian(sum(rS(iM).allMT_GL(:,:,1:length(biasVal),:,iST),[1 2]),4)),'-ok','LineWidth',2); hold on
    xlabel('bias amount')
    ylabel('Goal - Threat Defensiveness of actions')
    
    
    if iST == 1
        title('GOAL - THR neur change. Response to GOALS')
    elseif iST == 2
        title('GOAL - THR neur change. Response to THREATS')
    end
    
    [h p] = ttest2(squeeze( sum(rS(iM).allMT_GL(:,:,1,:,iST),[1 2])),squeeze( sum(rS(iM).allMT_THR(:,:,1,:,iST),[1 2])));
%     [p h] = ranksum(squeeze( sum(rS(iM).allMT_GL(:,:,1,:,iST),[1 2])),squeeze( sum(rS(iM).allMT_THR(:,:,1,:,iST),[1 2])));
% hold off;
end
end



subplot(1,length(stimTypes),1)
plot(biasVal,nanmedian(bvG(:,1,:),3)- nanmedian(bvT(:,1,:),3),'-or','LineWidth',2);

subplot(1,length(stimTypes),2)
plot(biasVal,nanmedian(bvG(:,2,:),3)- nanmedian(bvT(:,2,:),3),'-or','LineWidth',2);

% transform stim type to +1 -1
% allST = -(allST-1.5).*2;
allST = categorical(allST);
% allNT = categorical(allNT);


%% Try to compare the effects of positive and negative bias by subtracting 
% them (so looking for asymmetries in the response pattern
tbl = table(allMT,allNT,allST,allM,allB,allSEPpDC,allSEPpCD);

diffTbl = tbl(tbl.allB>=0,:);

aNT = unique(diffTbl.allNT);
aST = unique(diffTbl.allST);
aM  = unique(diffTbl.allM);
aB = unique(diffTbl.allB);

for iNT = 1:length(aNT)
    for iST = 1:length(aST)
        for iM = 1:length(aM)
        for iB = 1:length(aB)
            
            nt = aNT(iNT);
            st = aST(iST);
            m = aM(iM);
            b = aB(iB);
            
            currDiffRow = find(diffTbl.allNT == nt & diffTbl.allST == st & ...
                diffTbl.allM == m & diffTbl.allB == b);
            currNegRow = find(tbl.allNT == nt & tbl.allST == st & ...
                tbl.allM == m & (tbl.allB > -b - 0.000001) & (tbl.allB < -b + 0.000001));
            if length(currDiffRow) > 1 || length(currNegRow) > 1
                warning(['too many elements!'])
            end
            if diffTbl.allB(currDiffRow) >= 0
                diffTbl.allMT(currDiffRow) = diffTbl.allMT(currDiffRow) - tbl.allMT(currNegRow);
            end
        end
        end

    end
end



%% Run stats on the diff table


A=unique(allST);
B=unique(allNT);

% non-zero bias
diffTbl.nZB = diffTbl.allB > 0;

% maxB = 0.4;
% maxB = 0.2;
maxB = 0.1;
% maxB = 0.05;
% maxB = 0.03;


sD = diffTbl.allMT < 10000;
% sD = diffTbl.allB ~= 0;
% sD = diffTbl.allMT < 10000 & diffTbl.allM<=44 ;
% sD = diffTbl.allMT < 10000 & diffTbl.allM<=44  & diffTbl.allSEPpDC <=0.05 ;
% sD = diffTbl.allMT < 10000 & diffTbl.allM<=44  & diffTbl.allSEPpDC <=0.05 & diffTbl.allB<=maxB;
% sD = diffTbl.allMT < 10000 & diffTbl.allM>=45;
% sD = diffTbl.allMT < 10000 & diffTbl.allM>=45  & diffTbl.allSEPpDC <=0.05 & diffTbl.allB<=maxB;
% sD = diffTbl.allMT < 10000 & diffTbl.allM>=40 & diffTbl.allB<=maxB & diffTbl.allSEPpDC <=0.05 ;
% sD = diffTbl.allMT < 10000 ;
% sD = diffTbl.allMT < 10000 & diffTbl.allB<=maxB;
% sD = diffTbl.allMT < 10000 & diffTbl.allM<=40 & diffTbl.allB<=maxB;
% sD = diffTbl.allMT < 10000 & diffTbl.allM>=41 & diffTbl.allM<=80 & diffTbl.allB<=maxB;
% sD = diffTbl.allMT < 10000 & diffTbl.allM>=91 & diffTbl.allB<=maxB;

% sD = diffTbl.allMT < 10000 & diffTbl.allM>20 & diffTbl.allM<59;

% sD = diffTbl.allMT < 10000 & diffTbl.allB<=maxB & diffTbl.allSEPpDC <=0.05 ;
% sD = diffTbl.allMT < 10000 & diffTbl.allSEPpDC <=0.05 ;

diffTbl2 = diffTbl(sD,:);
diffLme = fitlme(diffTbl2,'allMT ~ allB*allST*allNT - allNT + (1|allM)')
% 
% diffLmeNZB = fitlme(diffTbl2,'allMT ~ nZB*allST*allNT - allNT + (1|allM)')

% 
% sD = diffTbl.allST == A(2) & diffTbl.allM>=81 & diffTbl.allB <= maxB ;
% sD = diffTbl.allST == A(1) & diffTbl.allB <= maxB & diffTbl.allSEPpDC <=1;
% sD = diffTbl.allST == A(1) & diffTbl.allB<=maxB & diffTbl.allM>=40 & diffTbl.allSEPpDC <=0.05 ; 
diffTbl2 = diffTbl(sD,:);
diffLme = fitlme(diffTbl2,'allMT ~ allB*allNT - allNT + (1|allM)')


% sD = diffTbl.allNT == B(2) & diffTbl.allB <= maxB ;
% % sD = diffTbl.allNT == B(2) & diffTbl.allB <= maxB & diffTbl.allSEPpDC <=0.05;
% % sD = diffTbl.allNT == B(2) & diffTbl.allB <= maxB & diffTbl.allM<=40 & diffTbl.allSEPpDC <=0.05 ;
% diffTbl2 = diffTbl(sD,:);
% diffLme = fitlme(diffTbl2,'allMT ~ allB*allST + (1|allM)')

 
sD = diffTbl.allST == A(2) & diffTbl.allNT  == B(1);
diffTbl2 = diffTbl(sD,:);
diffLme = fitlme(diffTbl2,'allMT ~ allB + (1|allM)')



figure, scatter(diffTbl2.allB,diffTbl2.allMT); lsline

figure, imagesc(hist3([diffTbl2.allB,diffTbl2.allMT])'); axis xy

% $$$ I NEED TO test at which points/whether the MT difference is
% larger/smaller than zero!



%%

% selectedData
% sD = allB >= 0;
sD = allB <= 10000;


tbl = table(allMT,allNT,allST,allM,allB);
tbl = tbl(sD,:);

% lme1 = fitlme(tbl,'allMT ~ allNT*allST + (1|allM)')
% lme2 = fitlme(tbl,'allMT ~ allST + (1|allM)')
% lme3 = fitlme(tbl,'allMT ~ allNT + allST + (1|allM)')

% lme1 = fitlme(tbl,'allMT ~ allB + allNT:allB + allST + (1|allM)')

lme = fitlme(tbl,'allMT ~ allB*allST*allNT - allNT - allNT:allB+ (1|allM)')
altlme = fitlme(tbl,'allMT ~ allB*allST*allNT - allNT + (1|allM)')


%       H0: Observed response vector was generated by model LME. 
%       H1: Observed response vector was generated by model ALTLME.
compare(lme,altlme)

% $$$ INVESTIGATE THIS TRIPLE INTERACTION BY SPLITTING IT UP!!


% $$$ CHECK HOW THIS WORKS
% compare(lme2,lme3)

%%

A=unique(allST);
sD = allST == A(1);

% sD = allB <= 0 & allST == A(2);
% sD = allB <= 1000;


tbl = table(allMT,allNT,allST,allM,allB);
tbl = tbl(sD,:);

% lme1 = fitlme(tbl,'allMT ~ allNT*allST + (1|allM)')
% lme2 = fitlme(tbl,'allMT ~ allST + (1|allM)')
% lme3 = fitlme(tbl,'allMT ~ allNT + allST + (1|allM)')

% lme1 = fitlme(tbl,'allMT ~ allB + allNT:allB + allST + (1|allM)')

lme1 = fitlme(tbl,'allMT ~ allB*allNT - allNT + (1|allM)')

%%


sD = allNT == 1;

% sD = allB <= 0 & allST == A(2);
% sD = allB <= 1000;


tbl = table(allMT,allNT,allST,allM,allB);
tbl = tbl(sD,:);

% lme1 = fitlme(tbl,'allMT ~ allNT*allST + (1|allM)')
% lme2 = fitlme(tbl,'allMT ~ allST + (1|allM)')
% lme3 = fitlme(tbl,'allMT ~ allNT + allST + (1|allM)')

% lme1 = fitlme(tbl,'allMT ~ allB + allNT:allB + allST + (1|allM)')

lme1 = fitlme(tbl,'allMT ~ allB*allST + (1|allM)')


%%
% sD = allB <= 0 & allST == A(2) & allNT == 1;
sD = allST == A(2) & allNT == -1;

tbl = table(allMT,allNT,allST,allM,allB);
tbl = tbl(sD,:);

lme1 = fitlme(tbl,'allMT ~ allB + (1|allM)')

% - so it looks like (with 15 models, and covariance definition) the response to goals is affected by
% biasing the goal network consistently: (nly for negative bias) more
% NEGATIVE bias leads to MORE defensive actions --> at least that fits.. 
% - This is for bias -0.05 to o (although I tested to +0.05), with the
% whole network altered. And 15 Models.

% IT ALSO LOOKS LIke (with 26 models, and covariance definition), the
% GOALNET change makes it so that bias has an effect: more bias means
% response to goals is more defensive, while response to threats is less
% defensive. THIS IS TRUE FOR -0.05 to +0.05 at least, and I 
% NEED TO DOUBLE CHECK
% WHETHER IT IS more or less true for bigger ranges of biases -->
% --> SURPRISINGLY SEEMS THE SAME FOR -.25 to +.25!
% --> Retrying this now with 39(?)models --> with -0.1 to +0.1, goal net
% has an effect with bias, but threat net doesn't... what the hell?
% ---> BUT WHAT ABOUT THE bias difference? --> THere is a net type bi bias
% interaction, which is very neat. BUT it might be driven almost completely
% by the goal net? --> again, looks like effect of bias with goal net is
% one way for threats and the other way for goals, while bias with threat
% net doesn't really do anything.. Does that mean I have to check how many
% neruons are goal or threat? Or is it cause negative valence stimuli lead
% to smaller Q-fields? --> hmmm.
% $$$ ONE THING IT DOES SHOW IS THAT the 'less defensiveness' is bigger
% than the 'more defensiveness', so that on the whole, biasing goal neurons
% positively makes things a bit more appettitive than biasing them
% negatively.
% $$$ --> I think try again with a SMALLER RANGE OF BIASES AND WITH STRICT
% COMBO?
% $$$ --> RIGHT NOW, I'm trying for all layer biasing (not just last 2) -->
% then after that I've got to try it for different bias range probably?
% --> SO FOR -0.1 to+0.1 IT WORKS PRETTY WELL! defensiveness goes in the
% direction expected for the different network types, and there isn't too
% much difference between stimulus types (well, except for threat when
% biasing threat network, which isn't very affected) --> TRYING FOR SMALLER
% BIAS RANGE ( for -0.01 to +0.01 it doesn't work,and for -0.05 to +0.05
% it's not very strong, but mainly the goal network instead of the threat
% network doing work. I'll try [-0.15 to -0.1 plus 0.1 to 0.15] now.. -->
% that is Ok ish, but mainly for goal net again? Although when leaving out
% 0, the threatnet goes in a similar diretion to goal net but less
% strongly?..
% --> Anyway, when going from -0.2 to +0.2, it's pretty good as well
% $$$ --> I tried it with more models, and the last 40 or 50 are less
% good... Even go in the opposite direction. Population wise, the original
% results still hold, but they are weaker..
% --> TRYING TO CALCULATE MORE OF THEM (models) AGAIN NOW -- increase iM

% :( ). NOW WILL HAVE TO TRY FOR THE STRICT COV plus HPROXN --> THEN
% will need to run more models again
% --> THE STRICT ONE LOOKS PRetTY GOOD, although some minor exceptions,
% goalnet looks like it makes things more grabby, and threat net less
% grabby

% $$$ NOW TRYING JUST hPROXN --> doesn't work
% $$$ NEXT TRY EXCLUSIVE HPROXN?

% $$$ --> now running again, but with ONLY the hPROXN

% --> for the loose combo covariance and other one defintion, I only get an
% effect of bias when using goals as stimuli, and the bias doesn't depend
% on network type :(

% $$$ AND individual neurons?

% $$$ AND I NEED TO PLOT THE AVERAGE RESPONSE FIELD TO GET A FEEL (haha)
% FOR HOW EXPANDED IT is/isn't maybe?

% $$$ AND I should still try to calculate the ACTUAL EXTENT of the
% appettitive/defensive responses, to get a better feel for whether these
% things expand or not

% $$$ I CAN ALSO try to make the combo one less selective --> have only as
% an extra condition that the neurons like goals statistically, but not
% that they don't like threats
% --> doing this now --> Doesn;t seem to work very well, for 29 models at
% least

% $$$ FOR THE COMBO ONE, it only works with THREATENING STIMULI... Also, I
% think the bias effect is the 'wrong way round'?...

% $$$ AND ULTIMATELY I need to make more model instances then, to ensure
% that the results are actually robust - I've been fishing a bit so far

% $$$ ONE THING I'M TRYING is subtracting the movement types from eachother
% for each level of BIAS - positive and negative -- maybe it works, who
% knows..














































%% $$$$$ !!!!! Try with other network and try with smaller bias value range

% $$$$$$ ALS, this I can fix quickly maybe: staying still sould count as an
% appetitive movement if the stimulus is directly above the limb
% $$$$$ OR OR jsut set everything on the midline to 0 because it is hard to
% interpret? --> DECIDE


iM = 2;

qBTDiff=rS(iM).allQbyType(:,:,:,:,:,2)-rS(iM).allQbyType(:,:,:,:,:,1);

figure,
for iNT=1:length(netTypes)
    for iST=1:length(stimTypes)
        plot(biasVal,squeeze(sum(qBTDiff(:,:,:,iNT,iST),[1 2])),'-o','LineWidth',2); hold on
        %        plot(biasVal,squeeze(sum(totalDef(:,:,:,iNT,iST),[1 2])),'-o','LineWidth',2); hold on
        %         plot(biasVal,squeeze(sum(totalGrab(:,:,:,iNT,iST),[1 2])),'-o','LineWidth',2); hold on
    end
end

xlabel('bias amount')
ylabel('Q value diff (def-appet)')
legend('Glnet, GlStim','GlNet, ThrStim','ThrNet, GlStim','ThrNet, ThrStim')

totalDef=rS(iM).allMT>0;
totalGrab=rS(iM).allMT<0;
figure,
for iNT=1:length(netTypes)
    for iST=1:length(stimTypes)
        plot(biasVal,squeeze(sum(rS(iM).allMT(:,:,:,iNT,iST),[1 2])),'-o','LineWidth',2); hold on
        %        plot(biasVal,squeeze(sum(totalDef(:,:,:,iNT,iST),[1 2])),'-o','LineWidth',2); hold on
        %         plot(biasVal,squeeze(sum(totalGrab(:,:,:,iNT,iST),[1 2])),'-o','LineWidth',2); hold on
    end
end

xlabel('bias amount')
ylabel('Defensiveness of actions')
legend('Glnet, GlStim','GlNet, ThrStim','ThrNet, GlStim','ThrNet, ThrStim')


figure,
for iNT=1:length(netTypes)
    
    plot(biasVal,squeeze(sum(rS(iM).allMT(:,:,:,iNT,:),[1 2 5])),'-o','LineWidth',2); hold on
    %        plot(biasVal,squeeze(sum(totalDef(:,:,:,iNT,iST),[1 2])),'-o','LineWidth',2); hold on
    %         plot(biasVal,squeeze(sum(totalGrab(:,:,:,iNT,iST),[1 2])),'-o','LineWidth',2); hold on
    
end

xlabel('bias amount')
ylabel('Defensiveness of actions')
legend('Glnet change','ThrNet change')

clear diffMTSpl
figure,
for iNT=1:length(netTypes)
    
    
    subplot(1,2,iNT)
    
    tmpMT=squeeze(sum(rS(iM).allMT(:,:,:,iNT,:),[1 2 5]));
    diffMT=(tmpMT-tmpMT(ceil(end/2)));
    diffMTSpl(1,:)=diffMT(1:floor(end/2));
    try % in case the size is even
        diffMTSpl(2,:)=diffMT(ceil(end/2):end);
    catch % in case the size is odd
        diffMTSpl(2,:)=diffMT((ceil(end/2)+1):end);
    end
    sdMT=(nanstd(diffMTSpl,[],2)./2)./sqrt(size(diffMTSpl,1));
    
    bar( nanmean(diffMTSpl,2) ); hold on
    plot([1 1 ; 2 2]', [nanmean(diffMTSpl,2) + [-sdMT sdMT]]','k','LineWidth',2);
    if iNT==1, title('changing goal network');
    elseif iNT==2, title('changing threat network'); end
    ylabel('Difference in defensivenes (+ve) vs goalness (-ve)')
    xlabel('Decr bias, Incr bias')
    
end






figure,
for iNT=1:length(netTypes)
    
    subplot(1,2,iNT)
    
    tmpMT=squeeze(sum(totalGrab(:,:,:,iNT,:),[1 2 5]));
    diffMT=(tmpMT./tmpMT(ceil(end/2))).*100 -100;
    diffMTSpl(1,:)=diffMT(1:floor(end/2));
    try % in case the size is even
        diffMTSpl(2,:)=diffMT(ceil(end/2):end);
    catch % in case the size is odd
        diffMTSpl(2,:)=diffMT((ceil(end/2)+1):end);
    end
    sdMT=(nanstd(diffMTSpl,[],2)./2)./sqrt(size(diffMTSpl,1));
    
    bar( nanmean(diffMTSpl,2) ); hold on
    plot([1 1 ; 2 2]', [nanmean(diffMTSpl,2) + [-sdMT sdMT]]','k','LineWidth',2);
    if iNT==1, title('changing goal network');
    elseif iNT==2, title('changing threat network'); end
    ylabel('percent change of appetitive responses')
    xlabel('Decr bias, Incr bias')
    
end

figure,
for iNT=1:length(netTypes)
    
    subplot(1,2,iNT)
    
    tmpMT=squeeze(sum(totalDef(:,:,:,iNT,:),[1 2 5]));
    diffMT=(tmpMT./tmpMT(ceil(end/2))).*100 -100;
    diffMTSpl(1,:)=diffMT(1:floor(end/2));
    try % in case the size is even
        diffMTSpl(2,:)=diffMT(ceil(end/2):end);
    catch % in case the size is odd
        diffMTSpl(2,:)=diffMT((ceil(end/2)+1):end);
    end
    sdMT=(nanstd(diffMTSpl,[],2)./2)./sqrt(size(diffMTSpl,1));
    
    bar( nanmean(diffMTSpl,2) ); hold on
    plot([1 1 ; 2 2]', [nanmean(diffMTSpl,2) + [-sdMT sdMT]]','k','LineWidth',2);
    if iNT==1, title('changing goal network');
    elseif iNT==2, title('changing threat network'); end
    ylabel('percent change of defensive responses')
    xlabel('Decr bias, Incr bias')
    
end


f1 = figure,
f2 = figure,
for iNT=1:length(netTypes)
    
    figure(f1)
    subplot(1,2,iNT)
    
    tmpRatio=(sum(totalDef(:,:,:,iNT,:),[1 2 5])-sum(totalGrab(:,:,:,iNT,:),[1 2 5])) ./ ...
        (sum(totalDef(:,:,:,iNT,:),[1 2 5])+sum(totalGrab(:,:,:,iNT,:),[1 2 5]));
    tmpMT=squeeze(tmpRatio);
    diffMT=(tmpMT./tmpMT(ceil(end/2))).*100 -100;
    diffMTSpl(1,:)=diffMT(1:floor(end/2));
    try % in case the size is even
        diffMTSpl(2,:)=diffMT(ceil(end/2):end);
    catch % in case the size is odd
        diffMTSpl(2,:)=diffMT((ceil(end/2)+1):end);
    end
    sdMT=(nanstd(diffMTSpl,[],2)./2)./sqrt(size(diffMTSpl,1));
    
    bar( nanmean(diffMTSpl,2) ); hold on
    plot([1 1 ; 2 2]', [nanmean(diffMTSpl,2) + [-sdMT sdMT]]','k','LineWidth',2);
    if iNT==1, title('changing goal network');
    elseif iNT==2, title('changing threat network'); end
    ylabel('defense to grab ratio relative to baseline (%)')
    xlabel('Decr bias, Incr bias')
    
    
    figure(f2)
    subplot(1,2,iNT)
    
    plot(biasVal,diffMT,'-o','LineWidth',2);
    ylabel('defense to grab ratio relative to baseline (%)')
    if iNT==1
        title('Bias Goal net')
    else
        title('Bias Threat net')
    end
    xlabel('biasval')
    
    
end











































































%% Check whether the midpoint of a Q-value sigmoid correlates with stimulus VALENCE
a=tic

% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\Valence_51Batch_Plus2_Plus4');
    
sFP = rS(iM).s;

% test Columns (columns at which to test the sigmoids
tC = 2:sFP.wrld.size(2)-1;

% all Threat Presences to test the extents with. 0 is just goal, 1 is just
% 'threat' (which could also have positive valence)
aTP = [0 1];

% Preallocate space for the extent of the value and neural fields
nLay = 1; maxN = 1;
for iM=1:length(rS)
    nLay = max([nLay, length(rS(iM).s.lp.netS)]);
    maxN = max([maxN, max(rS(iM).s.lp.netS)]);
end
qExtent = nan([length(aTP),length(tC),rS(1).s.act.numA,length(rS)]);
nExtent = nan([length(aTP),length(tC),nLay,maxN,length(rS)]);



% Loop through valence values (threat presence) and calculate 'extent of PPS'
for iV = 1:length(aTP)
    rS_sepV(iV).rS = rS;
    rS_sepV(iV).V = iV;
    
    b=tic
    for iM = 1:length(rS)
        
        sFP = rS(iM).s; 
        w = rS(iM).w;
        net = rS(iM).net;
        Qtable = rS(iM).Qtable;
        
        sFP.plt.plotThrFl = aTP(iV);
           
        sFP.plt.rowLims = [1.5 sFP.wrld.size(1)-0.5];
        sFP.plt.stimRow=[3:size(w.world2D,1)-3];
        
        % -------------------------------------------------------------------------
        % Calculate the correlati-ness of each velocity, for Q-values and neural
        % activity
        sFP.plt.meanLimbCols = 1;
        sFP.plt.fitSigmoid = 0;
        sFP.plt.lmbCol = [2:(sFP.wrld.size(2)-1)];
        sFP.plt.stimCol= [2:(sFP.wrld.size(2)-1)] ;
        [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
        
        [rS_sepV(iV).rS(iM).rDistQ, rS_sepV(iV).rS(iM).pDistQ, ...
            rS_sepV(iV).rS(iM).hProxQ] = ...
            CalcDistCorr(sFP,w,Q);
        
        [rS_sepV(iV).rS(iM).rDistN, rS_sepV(iV).rS(iM).pDistN, ...
            rS_sepV(iV).rS(iM).hProxN] = ...
            CalcDistCorr(sFP,w,permute(allNeurAct,[3 4 1 2 5 6]));
        
        
        % -------------------------------------------------------------------------
        % Loop through columns and calculate 'extent of PPS', for Q and neural
        % activity
        sFP.plt.meanLimbCols = 0;
        sFP.plt.fitSigmoid = 1;
        for iC = 1:length(tC)
            
            % current column
            cC = tC(iC);
            
            sFP.plt.lmbCol = cC;
            sFP.plt.stimCol= cC ;
                       
            [Q,allNeurAct] = CalcNetOutput(sFP,w,net);

            [rDistQ, pDistQ, ...
                hProxQ, aDQ, rDQ, cDQ, ...
                sigmQ, sigMidQ] = ...
                CalcDistCorr(sFP,w,Q);
            
            qExtent(iV,iC,:,iM) = sigMidQ.rD;
            
            [rDistN, pDistN, ...
                hProxN, aDN, rDN, cDN, ...
                sigmN, sigMidN] = ...
                CalcDistCorr(sFP,w,permute(allNeurAct,[3 4 1 2 5 6]));
            
            nExtent(iV,iC,1:size(sigMidN.rD,1),1:size(sigMidN.rD,2),iM) = sigMidN.rD;

        end
        iM
    end
    iV
    toc(b)
end
toc(a)

save('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\ExtentProcessed','-v7.3')

%% DISPLAY RESULTS - Q-value EXTENT dependence on valence (only positive)

% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\Valence\ExtentProcessed')

allV = repmat([sFP.act.GoalRew sFP.act.ThreatRew]',[1 size(qExtent,2)]);

clear p r rL rU
for iM = 1:length(rS)
    for iAct = 1:sFP.act.numA
        qEtmp = qExtent(:,:,iAct,iM);
%         hist(qEtmp'); legend('1','2')
% %         [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(allV(:),qEtmp);
% %         r(iAct,iM) = rTmp(2);
% %         p(iAct,iM) = pTmp(2);
% %         rL(iAct,iM) = rLTmp(2);
% %         rU(iAct,iM) = rUTmp(2);

        [pTmp,hTmp,stats] = signrank(allV(:),qEtmp(:));
        p(iAct,iM) = pTmp;
        zVal(iAct,iM) = stats.zval;
    end
end

p = CalcFDR(p);
    
[pMax pMaxInd] = max(p(:));

disp(['Relationship between stimulus valence and Q-field extent:'])
disp('')
disp(['P <= ' num2str(pMax) '. zVal > ' num2str(zVal(pMaxInd))])


%% Movement dependence on valence (positive and negative)

% % % $$$ TEMPORARY:
% % load('F:\Projects\DPPS\DefenseAgent\Results\Performance\Valence\NetSizes\Mults_of_12_51Batch_plus2_minus2_rewards.mat')
% % rS = ntRS([4 6]);

% Goal and threat stimulus
stimTypes=[0 1];

% Load movement types for every model and stimulus type/valence
for iM = 1:length(rS)
    
    sFP = rS(iM).s;
    w = rS(iM).w;
    net = rS(iM).net;
    Qtable = rS(iM).Qtable;
    
    sFP.plt.ON = 0;
    sFP.plt.rowLims = [1.5 sFP.wrld.size(1)-0.5];
    sFP.plt.lmbCol = 3:sFP.wrld.size(2)-2;
    sFP.plt.meanLimbCols = 1;
    sFP.plt.stimCol = [2:size(w.world2D,2)-1];
    
    
    for iST = 1:length(stimTypes)
        
        sFP.plt.plotThrFl=stimTypes(iST);
        
        [moveType(:,:,iM,iST) squishQs(:,:,:,iM,iST) QbyType(:,:,:,iM,iST)]=FindMoveTypes(net,sFP,w);
        
    end
    
end

% Reshape the movement types into vectors for stats
for iST = 1:length(stimTypes)
    tmp = moveType(:,:,:,iST);
    vecMT(:,iST) = tmp(:);
end


% Perform stats
opts.testType = 'signrank';
[figString zVal pVal] = CorrelationTitle([],vecMT,opts);

figure, hist(vecMT,3)
xlabel('movement type: -1 towards, 0 stay, +1 away')
ylabel('count')
legend('goal','threat')
title(figString);


%% Look at movement release as a function of distance?

iAct = 1;
iM = 2;
iST = 2;
% iST = 1:2;

imagesc(squishQs(:,:,iAct,iM,iST)); colorbar

% plot(squeeze(nanmean(squishQs(2:end-1,2:end-1,iAct,iM,iST),2)));

%% $$$ FOR NETWORK CHECKING

s.plt.startLayer=1;
% s.plt.stopLayer=9;
[netAR,s] = PlotConnections(net,s,w);


%% $$$ TEMP - checking whether 2 reward results are consistent wiht the only goal stuff 
% (cause maybe the 4 reward makes the Q-propagation less good)

% figure,

iCC=14;
sFP.plt.stimCol=iCC;
sFP.plt.lmbCol=iCC;

sFP.plt.plotThrFl = 1;
[QThr,allNeurAct] = CalcNetOutput(sFP,w,net);
sFP.plt.plotThrFl = 0;
[QGl,allNeurAct] = CalcNetOutput(sFP,w,net);


sFP.plt.pltType = 'Binned'
sFP.plt.distanceType='Row'
sFP.plt.fS.nBins=20
sFP.plt.fS.pltType = 'bar'

DisplActValsFun(sFP,w, QThr);
DisplActValsFun(sFP,w, QGl);

hold off

%%


function [npks pks locs] = RowPeakFind(s,Q)
    Qr = RecenterQNA(s,Q);
    Qline = squeeze(nanmean(Qr(1,:,:,:,:,:),[4]));
    for iLC = 2:size(Qline,2)-1
        [pks{iLC-1} locs{iLC-1}] = findpeaks(Qline(iLC,:,3));
        npks(iLC-1) = length(locs{iLC-1});
    end
    
end