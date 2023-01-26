%% Assess the effect of network size on performance

addpath('F:\Projects\DPPS\DefenseAgent\Scripts\rlsimplepps')
addpath('F:\Projects\DPPS\DefenseAgent\Scripts')
addpath('D:\Old_D\Programs\Matlab\Utilities\')
addpath('D:\Old_D\Programs\Matlab\Utilities\plotting\')
addpath('D:\Old_D\Programs\Matlab\Utilities\plotting\colormaps\')

load('F:\Projects\DPPS\DefenseAgent\Results\Performance\VelocityDirection\NetSizes\Mults_of_12_51Batch_V3');
load('F:\Projects\DPPS\DefenseAgent\Results\Performance\FullModel\NetSizes\Body_50_130Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE_FewerSpeeds_SmallWorld_V3.mat')


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
    
    % %     sFP.plt.meanLimbCols = 1;
    % %     sFP.plt.lmbCol=3:12;
    
    sFP.plt.meanLimbCols = 0;
    sFP.plt.lmbCol=8;
    
    sFP.plt.bdyCol=4;
    
    sFP.plt.plAct=2;
    
    sFP.plt.rowLag = 2;
    sFP.plt.colLag = 1;
    
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


%% Check whether the midpoint of a Q-value sigmoid correlates with movement speed
a=tic
% % % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\Fig_Vel_Dir_Dependence_101Batches_OLDRUNRELEARN.mat');
% % % load('F:\Projects\DPPS\DefenseAgent\Results\Performance\FullModel\NetSizes\Body_50_130Batch_plus2_2_minus2_4_rew_YesHist_DefRew-01_MOREBATCHSIZE_FewerSpeeds_SmallWorld_V3.mat')

% $$$ $$$ NEXT DO IS RUN FROM HERE, and plot velocity direction figure

rS(end+1) = ntRS;

iM = 1;
sFP = DefaultSettings(rS(iM).s);

% test Columns (columns at which to test the sigmoids
tC = 2:sFP.wrld.size(2)-1;


% figure

nLay = 1; maxN = 1;
for iM=1:length(rS)
    nLay = max([nLay, length(rS(iM).s.lp.netS)]);
    maxN = max([maxN, max(rS(iM).s.lp.netS)]);
end

% $$$ qExtent = nan([length(sFP.gol.alSpR),length(tC),rS(1).s.act.numA,length(rS)]);
qExtent = nan([length(sFP.gol.alSpR),length(tC),max([arrayfun(@(i) rS(i).s.act.numA, 1:length(rS))]),length(rS)]);
nExtent = nan([length(sFP.gol.alSpR),length(tC),nLay,maxN,length(rS)]);

% Loop through velocities (row lags) and calculate 'extent of PPS'
for iV = 1:3 %sFP.gol.alSpR
    rS_sepV(iV).rS = rS;
    rS_sepV(iV).V = iV;
    
    b=tic
    for iM = 1:length(rS)
        
        
        sFP = DefaultSettings(rS(iM).s); % $$$ HERE AND BELOW, change so that it takes the right iM $$$
        w = rS(iM).w;
        net = rS(iM).net;
        Qtable = rS(iM).Qtable;
        
        sFP.plt.rowLag = iV;
        
        
        
        sFP.plt.rowLims = [1.5 sFP.wrld.size(1)-0.5];
        sFP.plt.stimRow=[3:size(w.world2D,1)-3];
        
        
        % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        % $$$ maybe HERE CALCULATE MEAN FIRING OF EACH NEURON for each velocity
        % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
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
        
        
% % %         % -------------------------------------------------------------------------
% % %         % Loop through columns and calculate 'extent of PPS', for Q and neural
% % %         % activity
% % %         sFP.plt.meanLimbCols = 0;
% % %         sFP.plt.fitSigmoid = 1;
% % %         for iC = 1:length(tC)
% % %             
% % %             % current column
% % %             cC = tC(iC);
% % %             
% % %             sFP.plt.lmbCol = cC;
% % %             sFP.plt.stimCol= cC ;
% % %             
% % %             
% % %             
% % %             [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
% % %             
% % %             % %     [rS(iM).bdy.rDistQ, rS(iM).bdy.pDistQ, rS(iM).bdy.hProxQ, aDQ, rDQ, cDQ, sigm] = ...
% % %             % %                 CalcDistCorr(sFP,w,Q);
% % %             
% % %             % %             [rS_sepV(iV).rS(iM).rDistQ, rS_sepV(iV).rS(iM).pDistQ, ...
% % %             % %                 rS_sepV(iV).rS(iM).hProxQ, aDQ, rDQ, cDQ, ...
% % %             % %                 rS_sepV(iV).rS(iM).sigmQ, rS_sepV(iV).rS(iM).sigMidQ] = ...
% % %             % %                 CalcDistCorr(sFP,w,Q);
% % %             % %
% % %             % %             [rS_sepV(iV).rS(iM).rDistN, rS_sepV(iV).rS(iM).pDistN, ...
% % %             % %                 rS_sepV(iV).rS(iM).hProxN, aDN, rDN, cDN, ...
% % %             % %                 rS_sepV(iV).rS(iM).sigmN, rS_sepV(iV).rS(iM).sigMidN] = ...
% % %             % %                 CalcDistCorr(sFP,w,permute(allNeurAct,[3 4 1 2 5 6]));
% % %             
% % %             [rDistQ, pDistQ, ...
% % %                 hProxQ, aDQ, rDQ, cDQ, ...
% % %                 sigmQ, sigMidQ] = ...
% % %                 CalcDistCorr(sFP,w,Q);
% % %             
% % %             qExtent(iV,iC,1:length(sigMidQ.rD),iM) = sigMidQ.rD;
% % %             
% % %             [rDistN, pDistN, ...
% % %                 hProxN, aDN, rDN, cDN, ...
% % %                 sigmN, sigMidN] = ...
% % %                 CalcDistCorr(sFP,w,permute(allNeurAct,[3 4 1 2 5 6]));
% % %             
% % %             nExtent(iV,iC,1:size(sigMidN.rD,1),1:size(sigMidN.rD,2),iM) = sigMidN.rD;
% % %             
% % %             %             for iAct = 1:sFP.act.numA
% % %             %
% % %             %                 % $$$ Figure out whether this works for neurons as well, or if I need
% % %             %                 % to change more things...
% % %             %                 tmp = coeffvalues(sigmN.rD(iAct).curve);
% % %             %                 qExtent(iV,iC,iAct,iM) = tmp(3);
% % %             %
% % %             %                 tmp = coeffvalues(sigmQ.rD(iAct).curve);
% % %             %                 nExtent(iV,iC,iAct,iM) = tmp(3);
% % %             %             end
% % %             
% % %             
% % %         end
        iM
    end
    iV
    toc(b)
end
toc(a)
save('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\ExtentProcessed_V3','-v7.3')

%% DISPLAY RESULTS - Q-value dependence on speed

% load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\ExtentProcessed_V2');
load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\ExtentProcessed_V3');

allV = repmat([1:3]',[1 size(qExtent,2)]);

clear p r rL rU
for iM = 1:length(rS)
    for iAct = 1:sFP.act.numA
        qEtmp = qExtent(:,:,iAct,iM);
        [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(allV(:),qEtmp);
        r(iAct,iM) = rTmp(2);
        p(iAct,iM) = pTmp(2);
        rL(iAct,iM) = rLTmp(2);
        rU(iAct,iM) = rUTmp(2);
    end
end

p = CalcFDR(p,0.05,0);
    
[pMax pMaxInd] = max(p(:));

disp(['Relationship between stimulus velocity Q-field extent:'])
disp('')
disp(['P <= ' num2str(pMax) '. Rho > ' num2str(r(pMaxInd))])



%% DISPLAY RESULTS - Q-value dependence on direction (CALCULATE)

clear qDiffVec pCL rCL Qall

allColLags = [-2 -1 0 1 2];
allRowLags = [1 2 3];

Qall = NaN([14 15 14 15 9 length(allColLags) length(allRowLags) length(rS)]);


% $$$ Find Q-values for all models, row lags and column lags
for iM = 1:length(rS)
for iRL = 1:length(allRowLags)
for iCL = 1:length(allColLags)
    
    sFP = DefaultSettings(rS(iM).s);
    w = rS(iM).w;
    net = rS(iM).net;
    
    
    sFP.plt.rowLag = allRowLags(iRL);
    
    
    sFP.plt.colLag = allColLags(iCL);
    
    
    sFP.plt.meanLimbCols = 1;
    sFP.plt.fitSigmoid = 0;
    sFP.plt.lmbCol = [2:(sFP.wrld.size(2)-1)];
    sFP.plt.stimCol= [2:(sFP.wrld.size(2)-1)] ;
    
    % $$$ PUT BODY UNDER LIMB IF BODY IS INCLUDED
    if sFP.act.numA == 3
        [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
    elseif sFP.act.numA == 9
        clear Q
        for iLC = [1:(sFP.wrld.size(2))]
            sFP.plt.lmbCol = iLC;
            sFP.plt.bdyCol = iLC;
            [Qtmp,allNeurAct] = CalcNetOutput(sFP,w,net);
            Q(:,iLC,:,:,:) = Qtmp(:,iLC,:,:,:);
        end
    end
    
    
    
%     if iM  > 1 & size(Q,5) > size(Qall,5)
%         Qall(:,:,:,:,size(Qall,5)+1:size(Q,5),:,:,:) = NaN;
%     end
    
%     Qall(:,:,:,:,:,iCL,iRL,iM) = Q;
    Qall(size(Qall,1) + 1 - size(Q,1) : size(Qall,1),:,...
         size(Qall,3) + 1 - size(Q,3) : size(Qall,3),:,1:size(Q,5),iCL,iRL,iM) = Q;
    
end
iRL
end
iM
end


% Shift Q-values to the center to make it all comparable
QallOrig = Qall;
[Qall]   = RecenterQNA(sFP,Qall);



%% DISPLAY RESULTS - Q-value dependence on speed VERSION 2

% MAYBE FIND THE FIRST DISTANCE AT WHICH Q is ABOVE SOME VALUE?

%      lR lC sR sC A cL rL l
% Qall(14 15 14 15 9 length(allColLags) length(allRowLags) length(rS)]);

overThresh = NaN([3 9 4]);

clear rRL pRL allRL allOverThresh

for iM = 1:length(rS)
    
    for iAct = 1:1 %rS(iM).s.act.numA
        
        clear overThresh tmpRL
        
        % lC sR sC cL rL
%         qTmp = squeeze(Qall(end,:,:,:,iAct,:,:,iM));
        qTmp = squeeze(nanmean(Qall(end,:,:,:,:,:,:,iM),5));
        
%         qThresh = prctile(qTmp(:),85);
        qThresh = prctile(qTmp(:),90);
        
        % Find the extent of the field for each limb column
        lCs = 2:14;
        for iLC = 1:length(lCs)
            cLC = lCs(iLC);
            Qline = squeeze(nanmean(qTmp(cLC,:,cLC,:,:),[1 3 4]));
            
                    qThresh = prctile(Qline(:),90);
            
            
            
            for iRL = 1:size(Qline,2) % loop through ROW LAGS to find expansion size
                try
                    overThresh(iRL,iLC) = find(Qline(:,iRL) > qThresh,1);
%                     % Center of mass
%                     overThresh(iRL,iLC) = nansum(Qline(:,iRL).*(1:size(Qline,1))') ./ nansum(Qline(:,iRL))
                catch
                    overThresh(iRL,iLC) = NaN;
                end
                tmpRL(iRL,iLC) = iRL;
            end
        end
        
        allRL(:,:,iAct,iM)         = tmpRL;
        allOverThresh(:,:,iAct,iM) = overThresh;
        
        % Column Lag p-values
        [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(tmpRL(~isnan(overThresh)),overThresh(~isnan(overThresh)));
        rRL(iAct,iM) = rTmp(2);
        pRL(iAct,iM) = pTmp(2);
        
        
    end
    
end

format long
rRL
[dmy pRL2] = fdr(pRL(1:4));

disp('max p:')
[maxP maxInd] = max(pRL2);
maxP
disp('weakest rho')
minRHo = rRL(maxInd)
disp('mean rho')
meanRHo = mean(rRL)
disp('std rho')
stdRHo = std(rRL)
format short

% figure,
% plot(allRL(:),allOverThresh(:))


%% DISPLAY RESULTS - Q-value dependence on * Direction* VERSION 2

% MAYBE FIND THE FIRST DISTANCE AT WHICH Q is ABOVE SOME VALUE?

%      lR lC sR sC A cL rL l
% Qall(14 15 14 15 9 length(allColLags) length(allRowLags) length(rS)]);

clear rCL pCL allCL allOverThreshC tmpCL alloverThreshC overThreshC

for iM = 1:length(rS)
    
    for iAct = 1:1 %rS(iM).s.act.numA
        
        clear overThresh tmpCL
        
        % lC sR sC cL rL
%         qTmp = squeeze(Qall(end,:,:,:,iAct,:,:,iM));
        qTmp = squeeze(nanmean(Qall(end,:,:,:,:,:,:,iM),5));
        
        qThresh = prctile(qTmp(:),50);
        
        % Find the extent of the field for each limb row
        sRs = 5:12;
        for iSR = 1:length(sRs)
            cSR = sRs(iSR);
            Qline = squeeze(nanmean(qTmp(2:14,cSR,2:end-1,:,:),[1 2 5]));
            % Make it absolute positive
            Qline = Qline - min(Qline,[],1);
            
            
            
            for iCL = 1:size(Qline,2) % loop through COLUMN LAGS to find extension size
                try
%                     overThreshC(iCL,iR) = find(Qline(:,iCL) > qThresh,1); % $$$ I NEED TO FIND A BETTER THRESHOLD!
%                     [dmy overThreshC(iCL,iSR)] = max(Qline(:,iCL)); % $$$ I NEED TO FIND A BETTER THRESHOLD!
                    % Calculate center of mass
                    overThreshC(iCL,iSR) = sum(Qline(:,iCL).*(1:size(Qline,1))') ./ sum(Qline(:,iCL));
                catch
                    overThreshC(iCL,iSR) = NaN;
                end
                tmpCL(iCL,iSR) = iCL;
            end
        end
        
        % shift to mean for each limb column
        overThreshC = overThreshC - nanmean(overThreshC,1);
        
        allCL(:,:,iAct,iM)         = tmpCL;
        alloverThreshC(:,:,iAct,iM) = overThreshC;
        
        % Column Lag p-values
        [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(tmpCL(~isnan(overThreshC)),overThreshC(~isnan(overThreshC)));
        rCL(iAct,iM) = rTmp(2);
        pCL(iAct,iM) = pTmp(2);
        
        
    end
    
end

rCL
pCL

format long
rCL
[dmy pCL2] = fdr(pCL(1:4));

disp('min p:')
[maxP maxInd] = max(pCL2);
maxP
disp('weakest rho')
minRHo = rCL(maxInd)
disp('mean rho:')
meanRho = mean(rCL)
disp('std rho:')
stdRho = std(rCL)
format short

figure,plot(allCL(:),alloverThreshC(:))

% %%
% 
% % figure,
% plot(squeeze(nanmean(qTmp(2:end-1,11,2:end-1,1:5,1:3),[1 5])),'LineWidth',3); legend


%% Create RESULTS - neural activity's dependence on speed and direction (CALCULATE)
% $$$ NOTE: I also need to figure out which direction the receptive fields
% SHOULD be expanding in, I think?


clear qDiffVec pCL rCL

allColLags = rS(1).s.gol.alSpC;
allRowLags = rS(1).s.gol.alSpR;

% -------------------------------------------------------------------------
% Find Neural activation values for all models, row lags and column lags
for iM = 1:length(rS)
for iRL = 1:length(allRowLags)
for iCL = 1:length(allColLags)
    
    sFP = DefaultSettings(rS(iM).s);
    w = rS(iM).w;
    net = rS(iM).net;
    
    
    sFP.plt.rowLag = allRowLags(iRL);
    
    
    sFP.plt.colLag = allColLags(iCL);
    
    
    sFP.plt.meanLimbCols = 1;
    sFP.plt.fitSigmoid = 0;
    sFP.plt.lmbCol = [2:(sFP.wrld.size(2)-1)];
    sFP.plt.stimCol= [2:(sFP.wrld.size(2)-1)] ;
    [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
    
    allNeurAct = permute(allNeurAct,[3 4 1 2 5 6]);
    
    if iM ==1 && iRL == 1 && iCL == 1
        nAall = nan([size(allNeurAct(:,:,:,:,1)) 11 50 length(allColLags) length(allRowLags) length(rS)]);
    end
    
%     nAall(:,:,:,:,1:size(allNeurAct,5),1:size(allNeurAct,6),iCL,iRL,iM) = allNeurAct;
    nAall(size(nAall,1) + 1 - size(allNeurAct,1) : size(nAall,1),:, ...
          size(nAall,3) + 1 - size(allNeurAct,3) : size(nAall,3),:, ...
          1:size(allNeurAct,5),1:size(allNeurAct,6),iCL,iRL,iM) = allNeurAct; 
    
end
iRL
end
iM
end


% Shift Q-values to the center to make it all comparable
[nAall] = RecenterQNA(sFP,nAall);

% Average over row lags (stimulus velocities)
nAmean = nanmean(nAall,8);
% Average over rows, to find left-right peak position
nAline = nanmean(nAmean(1,:,:,:,:,:,:,1,:),[3]);

%% DISPLAY RESULTS - NEURACT- dependence on speed VERSION 2

% MAYBE FIND THE FIRST DISTANCE AT WHICH Q is ABOVE SOME VALUE?

%      lR lC sR sC lay neur cL rL l
% Qall(14 15 14 15 11 50 length(allColLags) length(allRowLags) length(rS)]);

[rMat] = DisplayProxStats(rS_sepV(3).rS);
% [rMat] = DisplayProxStats(rS);

overThreshFBN = NaN([3 11 50 4]);

clear rRL pRL allRL alloverThreshFBN

rRLN                = nan([11,50,4]);
pRLN                = nan([11,50,4]);
allRLN              = nan([3,13,11,50,4]);
alloverThreshFBN    = nan([3,13,11,50,4]);

for iM = 1:4 %length(rS)
    
    for iL = 1:length(rS(iM).s.lp.netS)
        for iN = 1:rS(iM).s.lp.netS(iL)
        
        clear overThreshFBN tmpRL
        
        % lC sR sC cL rL
%         naTmp = squeeze(nAall(end,:,:,:,iL,iN,:,:,iM));
%         naTmp = squeeze(nanmean(nAall(end,:,:,:,:,:,:,iM),5));
        % z-score, then absolute
%         naTmp = squeeze(nAall(end,:,:,:,iL,iN,:,:,iM));
%         naTmp = abs( (naTmp - nanmean(naTmp(:)))./nanstd(naTmp(:)));
        
        % z-score, then flip depending on direction of correlation
        naTmp = sign(rMat.rDistN(iL,iN,iM)) .* squeeze(nAall(end,:,:,:,iL,iN,:,:,iM));
        naTmp = abs( (naTmp - nanmean(naTmp(:)))./nanstd(naTmp(:)) - nanmin(naTmp(:)) );
        
        
%         naThresh = prctile(naTmp(:),85);
%         naThresh = prctile(naTmp(:),90);
        
        % Find the extent of the field for each limb column
        lCs = 2:14;
        for iLC = 1:length(lCs)
            cLC = lCs(iLC);
            Qline = squeeze(nanmean(naTmp(cLC,:,cLC,:,:),[1 3 4]));
            
            naThresh = prctile(Qline(:),90);
            
            for iRL = 1:size(Qline,2) % loop through ROW LAGS to find expansion size
                try
%                     overThreshFBN(iRL,iLC) = find(Qline(:,iRL) > naThresh,1); % $$$ HERE HERE CHANGE THIS MAYBE
                    % $$$ TRY INTERPOLATING
                    qLTmp = Qline(~isnan(Qline(:,iRL)),iRL);
                    overThreshFBN(iRL,iLC) = find(interp1(1:size(qLTmp,1),qLTmp,1:1:size(qLTmp,1)) > naThresh,1);
%                     overThreshFBN(iRL,iLC) = sum(Qline(:,iRL).*(1:size(Qline,1))') ./ sum(Qline(:,iRL));
                catch
                    overThreshFBN(iRL,iLC) = NaN;
                end
                tmpRL(iRL,iLC) = iRL;
            end
        end
        
        allRLN(:,:,iL,iN,iM)           = tmpRL;
        alloverThreshFBN(:,:,iL,iN,iM) = overThreshFBN;
        
        % Column Lag p-values
        [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(tmpRL(~isnan(overThreshFBN)),overThreshFBN(~isnan(overThreshFBN)));
        if numel(rTmp) > 1
            rRLN(iL,iN,iM) = rTmp(2);
            pRLN(iL,iN,iM) = pTmp(2);
        else
            % don't assign p-values if only 1 value
            rRLN(iL,iN,iM) = NaN;
            pRLN(iL,iN,iM) = NaN;
        end
        
        end
    end
    
end

% rRLN 
% pRLN

disp(' ')
disp(' -------------------- ')
disp('Proportion of row lag neurons overall')
rMat.overallRLneurProp = ...
                 sum(rMat.pDistN<0.05 & ...
                     pRLN<0.05 ,[1 2]) ./ sum(rMat.pDistN<0.05,[1 2]);
nanmean(rMat.overallRLneurProp)
nanstd(rMat.overallRLneurProp)

%% $$$ trying to figure out if z score then absolute is the right thing to do above

% figure,
% plot(allRL(:),alloverThreshFBN(:))

% Must also ind neurons which do and don't correlate with distance
% [rMat] = DisplayProxStats(rS);
[rMat] = DisplayProxStats(rS_sepV(1).rS);


rMat.pRowLagLN      = pRLN;
rMat.rRowLagLN      = rRLN;

rMat.propCorrNeur = nan([11,4]); % $$$
rMat.propRLExpNeur = nan([11,4]); % neurons whose field Expands with distance
rMat.propCorrAndRLNeur = nan([11,4]); % neurons who 'have a field' and is affected by distance
rMat.propCorrAndRLExpNeur = nan([11,4]); % neurons whos field EXPANDS with distance
rMat.propOfCorrRLNeur = nan([11,4]); % porportion of neurons which have a field that expand
rMat.propOfCorrRLExpNeur = nan([11,4]); % porportion of field neurons which have a field that expand

for iM = 1:length(rS)
    
    for iL = 1:length(rS(iM).s.lp.netS)
            
            
            rMat.propCorrNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2) ./ (rS(iM).s.lp.netS)';
        
            rMat.propRLExpNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(pRLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2) ./ (rS(iM).s.lp.netS)';
        
            rMat.propRLExpNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(pRLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2) ./ (rS(iM).s.lp.netS)';
            
            rMat.propCorrAndRLNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                pRLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2) ./(rS(iM).s.lp.netS)';
            
            rMat.propCorrAndRLExpNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                pRLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                rRLN(1:length(rS(iM).s.lp.netS),:,iM)<0,2) ./(rS(iM).s.lp.netS)'
            
            
            rMat.propOfCorrRLNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                pRLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2) ./ sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2);
            
            rMat.propOfCorrRLExpNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                pRLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                rRLN(1:length(rS(iM).s.lp.netS),:,iM)<0,2) ./ sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2);
            
    end

end

% $$$ THEN after all this stuff, I can go on to make the pie
% chart and the COLUMN LAG effect :) hooray
            
figure,plot(rMat.relLayDep,rMat.propCorrAndRLNeur)
figure,plot(rMat.relLayDep,rMat.propOfCorrRLNeur)



%% DISPLAY RESULTS - NEURACT dependence on * Direction* VERSION 2

% MAYBE FIND THE FIRST DISTANCE AT WHICH Q is ABOVE SOME VALUE?

%      lR lC sR sC layer neuron cL rL l
% Qall(14 15 14 15 11 50 length(allColLags) length(allRowLags) length(rS)]);

clear rCLN pCLN allCLN allOverThreshCN tmpCL alloverThreshCN overThreshCN

overThreshCLN = NaN([3 11 50 4]);

rCLN                = nan([11,50,4]);
pCLN                = nan([11,50,4]);
allCLN              = nan([5,11,11,50,4]); % col lag, stim row, layer, neuron, model
alloverThreshCN     = nan([5,11,11,50,4]);

for iM = 1:length(rS)
    
    
    for iL = 1:length(rS(iM).s.lp.netS)
        for iN = 1:rS(iM).s.lp.netS(iL)
            
         % $$$ HERE REPLACE WITH LAYER AND NEURON NUMBER
        
        clear overThreshC tmpCL
        
        % lC sR sC cL rL
% % %         naTmp = squeeze(nAall(end,:,:,:,iL,iN,:,:,iM));
%         naTmp = squeeze(nanmean(nAall(end,:,:,:,:,:,:,:,iM),5));
        naTmp = sign(rMat.rDistN(iL,iN,iM)) .* squeeze(nAall(end,:,:,:,iL,iN,:,:,iM));
        naTmp = abs( (naTmp - nanmean(naTmp(:)))./nanstd(naTmp(:)) - nanmin(naTmp(:)) );
        
        naThresh = prctile(naTmp(:),50);
        
        % Find the extent of the field for each limb row
        sRs = 2:12;
        for iSR = 1:length(sRs)
            cSR = sRs(iSR);
            
            % Because the neural values aren't guaranteed to be positive,
            % first z-score and absolute the values, so that the abs is the
            % outlier, corresponding to the peak of the bump
            % (although dividing by std isn't necessary cause we're just
            % looking for the max anyway)
% % %             Qline = abs(squeeze(nanmean(naTmp(2:13,cSR,2:end-1,:,:),[1 2 5])) - ...
% % %             nanmean(squeeze(nanmean(naTmp(2:13,cSR,2:end-1,:,:),[1 2 5])),1));
            % try multiplying it by the direction of correlation instead,
            % as above
%             Qline = sign(rMat.rDistN(iL,iN,iM)) .* (squeeze(nanmean(naTmp(2:13,cSR,2:end-1,:,:),[1 2 5])) - ...
%             nanmean(squeeze(nanmean(naTmp(2:13,cSR,2:end-1,:,:),[1 2 5])),1));

            Qline = squeeze(nanmean(naTmp(2:13,cSR,2:end-1,:,:),[1 2 5]));
            
            qThresh = prctile(Qline(:),90);

            for iCL = 1:size(Qline,2) % loop through COLUMN LAGS to find extension size
                try
%                     overThreshC(iCL,iR) = find(Qline(:,iCL) > naThresh,1); % $$$ I NEED TO FIND A BETTER THRESHOLD!
%                     [dmy overThreshCN(iCL,iSR)] = max(Qline(:,iCL)); % $$$ I NEED TO FIND A BETTER THRESHOLD!
                    overThreshCN(iCL,iSR) = sum(Qline(:,iCL).*(1:size(Qline,1))') ./ sum(Qline(:,iCL));
                catch
                    overThreshCN(iCL,iSR) = NaN;
                end
                tmpCL(iCL,iSR) = iCL;
            end
        end
        
        % shift to mean for each limb column
        overThreshCN = overThreshCN - nanmean(overThreshCN,1);
        
        allCLN(:,:,iL,iN,iM)          = tmpCL;
        alloverThreshCN(:,:,iL,iN,iM) = overThreshCN;
        
        % Column Lag p-values
        [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(tmpCL(~isnan(overThreshCN)),overThreshCN(~isnan(overThreshCN)));
        rCLN(iL,iN,iM) = rTmp(2);
        pCLN(iL,iN,iM) = pTmp(2);
        
        end 
    end
    
end

rCLN
pCLN



% Then put the data into rMat

rMat.pColLagLN      = pCLN;
rMat.rColLagLN      = rCLN;

rMat.propCLExpNeur = nan([11,4]); % neurons whose field Rxpands with distance
rMat.propCorrAndCLNeur = nan([11,4]); % neurons who 'have a field' and it expands with distance
rMat.propOfCorrCLNeur = nan([11,4]); % porportion of neurons which have a field that expand

rMat.propOfCorrVELNeur = nan([11,4]); % porportion of neurons which have a field that depends on all velocity things

for iM = 1:length(rS)
    
    for iL = 1:length(rS(iM).s.lp.netS)
        for iN = 1:rS(iM).s.lp.netS(iL)
                        
            rMat.propCLExpNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(pCLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2) ./ (rS(iM).s.lp.netS)';
            
            rMat.propCorrAndCLNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                pCLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2) ./(rS(iM).s.lp.netS)';
            
            rMat.propOfCorrCLNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                pCLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2) ./ sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2);
            
            rMat.propOfCorrVELNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                 sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                     pCLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                     pRLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2) ./ sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2);
                 
             rMat.propCorrAndVELNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                 sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                     pCLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                     pRLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2) ./ (rS(iM).s.lp.netS)';
                 
             rMat.propCorrNotAnyVELNeur(1:length(rS(iM).s.lp.netS),iM) = ...
                 sum(rMat.pDistN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                     pCLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05 & ...
                     pRLN(1:length(rS(iM).s.lp.netS),:,iM)<0.05,2) ./ (rS(iM).s.lp.netS)';
            
        end
    end

end

% $$$ THEN after all this stuff, I can go on to make the pie
% chart and the COLUMN LAG effect :) hooray
            
% figure,plot(rMat.relLayDep,rMat.propCorrAndCLNeur)
figure,plot(rMat.relLayDep,rMat.propOfCorrCLNeur)

% Run an LME
figure,plot(rMat.relLayDep,rMat.propOfCorrVELNeur)

opts.testType = 'lme';
[figString rhoVal pVal] = CorrelationTitle(rMat.relLayDep,rMat.propOfCorrVELNeur,opts);

figure, plot (rMat.relLayDep,rMat.propOfCorrVELNeur,'-o','LineWidth',2);
title(figString);
xlabel('Relative layer depth')
ylabel('Proportion of neurons with PPS properties')

disp(['Relationship between neural layer PPS-ness: ' figString]);disp(' ')


disp(' ')
disp(' -------------------- ')
disp('Proportion of velocity neurons overall')
% % % rMat.overallVELneurProp = ...
% % %                  sum(rMat.pDistN<0.05 & ...
% % %                      pCLN<0.05 & ...
% % %                      pRLN<0.05,[1 2]) ./ sum(rMat.pDistN<0.05,[1 2]);
rMat.overallVELneurProp = ...
                 sum(rMat.pDistN<0.05 & ...
                     pCLN<0.05 & ...
                     pRLN<0.05,[1 2]) ./ sum(rMat.pDistN<100,[1 2]);
nanmean(rMat.overallVELneurProp)
nanstd(rMat.overallVELneurProp)

disp(' ')
disp(' -------------------- ')
disp('Proportion of row lag neurons overall')
% % % rMat.overallRLneurProp = ...
% % %                  sum(rMat.pDistN<0.05 & ...
% % %                      pRLN<0.05 ,[1 2]) ./ sum(rMat.pDistN<0.05,[1 2]);
rMat.overallRLneurProp = ...
                 sum(pRLN<0.05 ,[1 2]) ./ sum(rMat.pDistN<100,[1 2]);
nanmean(rMat.overallRLneurProp)
nanstd(rMat.overallRLneurProp)

disp(' ')
disp(' -------------------- ')
disp('Proportion of column lag neurons overall')
% % % rMat.overallCLneurProp = ...
% % %                  sum(rMat.pDistN<0.05 & ...
% % %                      pCLN<0.05 ,[1 2]) ./ sum(rMat.pDistN<0.05,[1 2]);
rMat.overallCLneurProp = ...
                 sum(pCLN<0.05 ,[1 2]) ./ sum(rMat.pDistN<100,[1 2]);
nanmean(rMat.overallCLneurProp)
nanstd(rMat.overallCLneurProp)

% 
% %% DISPLAY RESULTS - Q-value dependence on direction (DISPLAY)
% 
% % $$$ SO HERE I'M CHECKING WHICH ROWS TO INCLUDE
% 
% % Average over row lags (stimulus velocities)
% Qmean = nanmean(Qall,7);
% % Average over rows, to find left-right peak position
% Qline = nanmean(Qmean(end,:,2:13,:,:,:,1,:),[3]); % $$$ MAYTBE MAKE qline already average over limb cols
% 
% rCL = nan([9 length(rS)]);
% pCL = nan([9 length(rS)]);
% 
% for iM = 1:length(rS)
% 
%     % ---------------------------------------------------------------------
%     % Set up and run LEFT vs RIGHT difference LME
%     for iCL = 1:length(allColLags)
%         qDiffTmp = squeeze(Qmean(end,:,:,1:floor(end/2),:,iCL,1,iM) - Qmean(end,:,:,ceil(1+end/2):end,:,iCL,1,iM));
%         tmp = qDiffTmp(:);
%         qDiffVec(:,iCL) = tmp;
%     end
%     % column Lag vector
%     cLvec = repmat([-2:2],[size(qDiffVec,1) 1]);
%     opts.testType = 'lme';
%     [figStrings{iM} rhoVals(iM) pVals(iM)] = CorrelationTitle(cLvec',qDiffVec',opts);
%     
%     
%     % ---------------------------------------------------------------------
%     % Set up and run CORRELATION between PEAK POSITION and column lag    
%     % Shift the function, so that extreme deviations from the norm are easy to find
%     baseQ = nanmedian(Qmean(2:end-1,2:end-1,2:end-1,2:end-1,:,:,1,iM),[1 2 3 4]);
%     Qshift = squeeze(Qline - baseQ);
% 
%     % Find position at which Q-values are most extreme
%     [maxAbsQ maxQPos ] = nanmax(abs(Qshift([2:end-1],[2:end-1],:,:,:)),[],2);
%     maxQPos = squeeze(maxQPos);
% %     [maxAbsQ maxQPos ] = nanmax(nanmedian(abs(Qshift([2:end-1],[2:end-1],:,:,:)),1),[],2);
% %     maxQPos = permute(squeeze(maxQPos),[4 1 2 3]);
%     
%     % all column lag, for the format of the correlation tests
%     clear tmpCL
%     tmpCL(1,1,1:length(allColLags)) = [allColLags];
%     allCL = repmat(tmpCL,[size(maxQPos,1), 1]);
%     
%     % Column Lag p-values
%     for iAct = 1:rS(iM).s.act.numA
%         posTmp = maxQPos(:,iAct,:,iM);
%         [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(allCL(:),posTmp(:));
%         rCL(iAct,iM) = rTmp(2);
%         pCL(iAct,iM) = pTmp(2);
%     end
% 
% end
% 
% pVals = CalcFDR(pVals)
% pCL = CalcFDR(pCL);
% 
% disp(['Relationship between stimulus direction and Q-field direction (using left-right difference), at least:'])
% disp('')
% disp([figStrings{pVals==max(pVals)}])
% 
% 
% 
% [pClMax pClMaxInd] = max(pCL(:));
% 
% disp(['Relationship between direction and Q-field extent:'])
% disp('')
% disp(['P <= ' num2str(pClMax) '. Rho < ' num2str(rCL(pClMaxInd))])
% 
% % $$$ TO PLOT RESPONSE AROUND LIMB for 2nd action 
% %     in limb co0lumn 7 FOR the 4th model:
% % figure,imagesc(squeeze(Qmean(end,7,:,:,2,3,1,4)))
% 
% % $$$ NEXT: keep trying to figure out why the hand-centred stuff is bad - maybe try doing body-centred instead?
% % $$$ OR if it doesn't wokr, just use the body's action fields as an excuse
% 
% %% DISPLAY RESULTS - Neural activity's dependence on speed
% 
% % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\ExtentProcessed');
% 
% % -------------------------------------------------------------------------
% % First find neurons which have body-part centric coding
% for iV  = 1:length(rS_sepV)
%     [rMat(iV)] = DisplayProxStats(rS_sepV(iV).rS);
%     hProxN_all(:,:,:,iV) = rMat(iV).hProxN;
% end
% hProxN = min(hProxN_all,[],4);
% 
% 
% 
% allV = repmat([1:3]',[1 size(nExtent,2)]);
% 
% % -------------------------------------------------------------------------
% % Fins neurons whose response profile correlates with stimulus velocity
% s = rS(1).s;
% r = nan([length(s.lp.netS), max(s.lp.netS), length(rS)]);
% p = r; rL = r; rU = r; 
% for iM = 1:length(rS)
%     s = rS(iM).s;
%     for iL = 1:length(s.lp.netS)
%         for iN = 1:s.lp.netS(iL)
%             qEtmp = nExtent(:,:,iL,iN,iM);
%             [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(allV(:),qEtmp);
%             r(iL,iN,iM) = rTmp(2);
%             p(iL,iN,iM) = pTmp(2);
%             rL(iL,iN,iM) = rLTmp(2);
%             rU(iL,iN,iM) = rUTmp(2);
%         end
%     end
% end
% 
% % % % % hSpeedN = CalcFDR(p)<0.05 & hProxN == 1;
% % % % % hExpN = CalcFDR(p)<0.05 & hProxN == 1 & r>0;
% % % % % nNeur = sum(~isnan(p),2);
% % % % % propSpeedNPerLay = squeeze(sum(hSpeedN ,2)./nNeur);
% % % % % propExpNPerLay = squeeze(sum(hExpN,2)./nNeur);
% % % % % 
% % % % % propSpeedN = squeeze(sum( hSpeedN,[1 2])./sum(nNeur,1));
% % % % % propExpN = squeeze(sum(hExpN,[1 2])./sum(nNeur,1));
% 
% 
% % Calculate proportions of speed-sensitive and expanding neurons. All as a
% % function of the neurons which also are limb-centric
% hSpeedN = CalcFDR(p)<0.05 == 1 & hProxN == 1;
% hExpN = CalcFDR(p)<0.05 == 1 & r>0 & hProxN == 1;
% nNeur = sum(~isnan(p) & hProxN == 1,2);
% propSpeedNPerLay = squeeze(sum(hSpeedN ,2)./nNeur);
% propExpNPerLay = squeeze(sum(hExpN,2)./nNeur);
% 
% propSpeedN = squeeze(sum( hSpeedN,[1 2])./sum(nNeur,1));
% propExpN = squeeze(sum(hExpN,[1 2])./sum(nNeur,1));
% 
% 
% 
% relLayDep = nan(size(propExpNPerLay));
% for iM = 1:length(rS)
%     maxNeur = length(rS(iM).s.lp.netS);
%     relLayDep(1:maxNeur,iM) = linspace(0,1,maxNeur);
% end
% 
% figure, plot(relLayDep,propExpNPerLay./propSpeedNPerLay,'-o','LineWidth',2)
% ylabel('proportion of SPEED neurons with EXPANDING receptive fields');
% xlabel('layer depth')
% opts.testType = 'lme';
% [figString rhoVal pVal] = CorrelationTitle(relLayDep,propExpNPerLay./propSpeedNPerLay,opts);
% title(figString);
% 
% 
% figure, plot(relLayDep,propSpeedNPerLay,'-o','LineWidth',2)
% ylabel('proportion of neurons with SPEED-SENSITIVE receptive fields');
% xlabel('layer depth')
% opts.testType = 'lme';
% [figString rhoVal pVal] = CorrelationTitle(relLayDep,propSpeedNPerLay,opts);
% title(figString);
% 
% 
% figure, plot(relLayDep,propExpNPerLay,'-o','LineWidth',2)
% ylabel('proportion of neurons with EXPANDING receptive fields');
% xlabel('layer depth')
% opts.testType = 'lme';
% [figString rhoVal pVal] = CorrelationTitle(relLayDep,propExpNPerLay,opts);
% title(figString);
% 
% 
% disp(['Proportion of limb-centric neurons that CORRELATE with movement speed'])
% disp(['Av: ' num2str(nanmean(propSpeedN)) '+-' num2str(nanstd(propSpeedN)) ])
% 
% disp(['Proportion of limb-centric neurons that EXPAND with movement speed'])
% disp(['Av: ' num2str(nanmean(propExpN)) '+-' num2str(nanstd(propExpN)) ])
% 
% 
% disp(['Relationship layer depth and proportion of "expanding" limb-centric neurons:'])
% disp('')
% disp([figString])
% 
% 
% %% DISPLAY RESULTS - neural activity's dependence on direction
% % $$$ NOTE: I also need to figure out which direction the receptive fields
% % SHOULD be expanding in, I think?
% 
% 
% % -------------------------------------------------------------------------
% % First find neurons which have body-part centric coding
% for iV  = 1:length(rS_sepV)
%     [rMat(iV)] = DisplayProxStats(rS_sepV(iV).rS);
%     hProxN_all(:,:,:,iV) = rMat(iV).hProxN;
% end
% hProxN = min(hProxN_all,[],4);
% 
% 
% clear qDiffVec pCL rCL
% 
% 
% allColLags = rS(1).s.gol.alSpC;
% allRowLags = rS(1).s.gol.alSpR;
% 
% % -------------------------------------------------------------------------
% % Find Neural activation values for all models, row lags and column lags
% for iM = 1:length(rS)
% for iRL = 1:length(allRowLags)
% for iCL = 1:length(allColLags)
%     
%     sFP = DefaultSettings(rS(iM).s);
%     w = rS(iM).w;
%     net = rS(iM).net;
%     
%     
%     sFP.plt.rowLag = allRowLags(iRL);
%     
%     
%     sFP.plt.colLag = allColLags(iCL);
%     
%     
%     sFP.plt.meanLimbCols = 1;
%     sFP.plt.fitSigmoid = 0;
%     sFP.plt.lmbCol = [2:(sFP.wrld.size(2)-1)];
%     sFP.plt.stimCol= [2:(sFP.wrld.size(2)-1)] ;
%     [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
%     
%     allNeurAct = permute(allNeurAct,[3 4 1 2 5 6]);
%     
%     if iM ==1 && iRL == 1 && iCL == 1
%         nAall = nan([size(allNeurAct(:,:,:,:,1)) 11 50 length(allColLags) length(allRowLags) length(rS)]);
%     end
%     
% %     nAall(:,:,:,:,1:size(allNeurAct,5),1:size(allNeurAct,6),iCL,iRL,iM) = allNeurAct;
%     nAall(size(nAall,1) + 1 - size(allNeurAct,1) : size(nAall,1),:, ...
%           size(nAall,3) + 1 - size(allNeurAct,3) : size(nAall,3),:, ...
%           1:size(allNeurAct,5),1:size(allNeurAct,6),iCL,iRL,iM) = allNeurAct;
%     
%     
% end
% iRL
% end
% iM
% end
% 
% % Shift Q-values to the center to make it all comparable
% [nAall] = RecenterQNA(sFP,nAall);
% 
% % Average over row lags (stimulus velocities)
% nAmean = nanmean(nAall,8);
% % Average over rows, to find left-right peak position
% nAline = nanmean(nAmean(1,:,:,:,:,:,:,1,:),[3]);
% 
% clear nAall
% 
% 
% % Preset some values to NaN matrices so that 0s don't get counted wrong
% maxL = 1;
% maxN = 1;
% for iM=1:length(rS)
%     maxL = max([maxL length(rS(iM).s.lp.netS)]);
%     maxN = max([maxN rS(iM).s.lp.netS]);
% end
% pCL = nan([maxL, maxN, length(rS)]);
% rCL = pCL;
% 
% clear tmpCL
% for iM = 1:length(rS)
% 
% % %     % ---------------------------------------------------------------------
% % %     % Set up and run LEFT vs RIGHT difference LME
% % %     for iCL = 1:length(allColLags)
% % %         qDiffTmp = squeeze(nAmean(1,:,:,1:floor(end/2),:,:,iCL,1,iM) - nAmean(1,:,:,ceil(1+end/2):end,:,:,iCL,1,iM));
% % %         tmp = qDiffTmp(:);
% % %         qDiffVec(:,iCL) = tmp;
% % %     end
% % %     % column Lag vector
% % %     cLvec = repmat([-2:2],[size(qDiffVec,1) 1]);
% % %     opts.testType = 'lme';
% % %     [figStrings{iM} rhoVals(iM) pVals(iM)] = CorrelationTitle(cLvec',qDiffVec',opts);
%     
%     
%     % ---------------------------------------------------------------------
%     % Set up and run CORRELATION between PEAK POSITION and column lag    
%     % Shift the function, so that extreme deviations from the norm are easy to find
%     baseNA = nanmedian(nAmean(2:end-1,2:end-1,2:end-1,2:end-1,:,:,:,1,iM),[1 2 3 4]);
%     nAshift = squeeze(nAline(:,:,:,:,:,:,:,:,iM) - baseNA);
% 
%     % Find position at which Q-values are most extreme
%     [maxAbsnA maxnAPos ] = nanmax(abs(nAshift([2:end-1],[2:end-1],:,:,:,:)),[],2);
%     maxnAPos = squeeze(maxnAPos);
%     
%     % all column lag, for the format of the correlation tests
%     tmpCL(1,1,1,1:length(allColLags)) = [allColLags];
%     allCL = repmat(tmpCL,[size(maxnAPos,1), 1]);
%     
%     % Column Lag p-values
%     for iL = 1:length(rS(iM).s.lp.netS)
%         for iN = 1:rS(iM).s.lp.netS(iL)
%             posTmp = maxnAPos(:,iL,iN,:);
%             [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(allCL(:),posTmp(:));
%             rCL(iL,iN,iM) = rTmp(2);
%             pCL(iL,iN,iM) = pTmp(2);
%         end
%     end
% 
% end
% 
% pVals = CalcFDR(pVals);
% pCL = CalcFDR(pCL);
% 
% 
% 
% % -------------------------------------------------------------------------
% % Calculate proportions of direction sensitive and direciton-expanding neurons. All as a
% % function of the neurons which also are limb-centric
% hSpeedN = CalcFDR(pCL)<0.05 & hProxN == 1;
% hExpN = CalcFDR(pCL)<0.05 & rCL<0 & hProxN == 1;
% nNeur = sum(~isnan(pCL) & hProxN == 1,2);
% propSpeedNPerLay = squeeze(sum(hSpeedN ,2)./nNeur);
% propExpNPerLay = squeeze(sum(hExpN,2)./nNeur);
% 
% propSpeedN = squeeze(sum( hSpeedN,[1 2])./sum(nNeur,1));
% propExpN = squeeze(sum(hExpN,[1 2])./sum(nNeur,1));
% 
% 
% relLayDep = nan(size(propExpNPerLay));
% for iM = 1:length(rS)
%     maxNeur = length(rS(iM).s.lp.netS);
%     relLayDep(1:maxNeur,iM) = linspace(0,1,maxNeur);
% end
% 
% figure, plot(relLayDep,propExpNPerLay./propSpeedNPerLay,'-o','LineWidth',2)
% ylabel('proportion of DIRECTIONAL neurons with PROPERLY DIRECTIONS receptive fields');
% xlabel('layer depth')
% opts.testType = 'lme';
% [figString rhoVal pVal] = CorrelationTitle(relLayDep,propExpNPerLay./propSpeedNPerLay,opts);
% title(figString);
% 
% 
% figure, plot(relLayDep,propSpeedNPerLay,'-o','LineWidth',2)
% ylabel('proportion of neurons with DIRECTIONAL receptive fields');
% xlabel('layer depth')
% opts.testType = 'lme';
% [figString rhoVal pVal] = CorrelationTitle(relLayDep,propSpeedNPerLay,opts);
% title(figString);
% 
% 
% figure, plot(relLayDep,propExpNPerLay,'-o','LineWidth',2)
% ylabel('proportion of neurons with PROPERLY DIRECTIONAL receptive fields');
% xlabel('layer depth')
% opts.testType = 'lme';
% [figString rhoVal pVal] = CorrelationTitle(relLayDep,propExpNPerLay,opts);
% title(figString);
% 
% 
% disp(['Proportion of limb-centric neurons that change their left/right responsiveness'])
% disp(['Av: ' num2str(nanmean(propSpeedN)) '+-' num2str(nanstd(propSpeedN)) ])
% 
% disp(['Proportion of limb-centric neurons that EXPAND in the direction of incoming movement'])
% disp(['Av: ' num2str(nanmean(propExpN)) '+-' num2str(nanstd(propExpN)) ])
% 
% 
% disp(['Relationship layer depth and proportion of "expanding" limb-centric neurons:'])
% disp('')
% disp([figString])
% 
% 
% 
% %%
% % figure, 
% 
% iM = 3;
% 
% iL = 9;
% iN = 7;
% iLC = 11;
% % plot(Qshift(:,iLC,iL,iN,iCL)')
% 
% iSideVel = 3;
% 
% iSR =9;
% plot(squeeze(nAmean(1,iLC,iSR,:,iL,iN,iSideVel,1,iM)))
% 
% 
% % $$$ IT LOOKS LIKE stimrow doesn't matter for some reason??...
% %%
% 
% squeeze(nAmean(1,iLC,iSR+0,:,iL+0,iN,iSideVel,1,iM))
% 
% squeeze(nAall(1,iLC,iSR+0,:,iL+-1,iN,iSideVel,1,iM))
% 
% 
% %%
% 
% squeeze(allNeurAct(5,iLC,iSR+1,:,iL,iN))
% 
% %%
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% Make plots to get better idea of what stats should be
% 
% 
% % load('F:\Projects\DPPS\DefenseAgent\Results\ForFigures\VelocityDirection\Fig_Vel_Dir_Dependence_101Batches');
% 
% iM=2;
% 
% sFP = rS(iM).s; % $$$ HERE AND BELOW, change so that it takes the right iM $$$
% w = rS(iM).w;
% net = rS(iM).net;
% Qtable = rS(iM).Qtable;
% 
% sFP.plt.rowLims = [1.5 sFP.wrld.size(1)-0.5];
% sFP.plt.lmbCol = [2:size(w.world2D,2)-1];
% 
% % sFP.plt.stimRow=[3:size(w.world2D,1)-1];
% sFP.plt.stimRow=[3:size(w.world2D,1)-2];
% sFP.plt.stimCol=[2:size(w.world2D,2)-1];
% 
% % %
% % % sFP.plt.rowLims = [1.5 sFP.wrld.size(1)-0.5];
% % % sFP.plt.lmbCol = 8;
% % %
% % % sFP.plt.stimRow = [3:size(w.world2D,1)-1];
% % % sFP.plt.stimCol= 8 ;
% 
% sFP.plt.meanLimbCols = 1;
% 
% sFP.plt.fS.nBins = 5;
% 
% sFP.plt.colLag = 0;
% sFP.plt.rowLag = 3;
% 
% 
% % sFP.plt.distanceType = 'AbsColumn';
% sFP.plt.distanceType = 'Row';
% 
% % sFP.plt.pltType = 'Imagesc';
% sFP.plt.pltType = 'Binned';
% 
% % sFP.plt.fS.pltType = 'shadeplot';
% sFP.plt.fS.pltType = 'bar';
% 
% sFP.plt.varargs={'LineWidth',2};
% 
% 
% 
% 
% 
% 
% [Q,allNeurAct] = CalcNetOutput(sFP,w,net);
% DisplActValsFun(sFP,w,Q);
% 
% % Q3=squeeze(Q(12,:,:,:,:));

%% MOre for plot

% THIS IS PLOTTING - for figure 3
% % %     cAct = 1:3;
% % %     
% % %     qDiff = squeeze(Q(1,:,:,1:floor(end/2),cAct) - Q(1,:,:,ceil(1+end/2):end,cAct));
% % %     
% % %     qDifff(:,iCL) = qDiff(:);
% % %     %
% % %     %     hist(qDiff(:),50);
% % %     %     hold on
% % %     
% % %     % -------------------------------------------------------------------------
% % %     % Plot the column-dependence
% % %     subplot(length(allColLags),1,iCL)
% % %     sFP.plt.rowLims = [1.5 sFP.wrld.size(1)-0.5];
% % %     sFP.plt.colLims = [1.5 sFP.wrld.size(2)-0.5];
% % %     
% % %     sFP.plt.stimRow = [2:sFP.wrld.size(1)-2];
% % %     sFP.plt.stimCol = [2:sFP.wrld.size(2)-1];
% % %     
% % %     
% % %     sFP.plt.pltType = 'Binned';
% % %     sFP.plt.distanceType = 'Column';
% % %     sFP.plt.fS.nBins = 13;
% % % %     sFP.plt.fS.pltType = 'shadeplot';
% % %     sFP.plt.fS.pltType = 'bar';
% % % %     sFP.plt.fS.pltType = 'plot';
% % % 
% % % %     sFP.plt.fS.eBars = 'centile';
% % %     sFP.plt.fS.eBars = 'SD';
% % %     
% % %     DisplActValsFun(sFP,w, Q); hold on
% % %     vline(0,'-.k')
% % %     % -------------------------------------------------------------------------


%% More for plot; 

% 
% sFP.plt.rowLims = [1.5 sFP.wrld.size(1)-0.5];
% sFP.plt.colLims = [1.5 sFP.wrld.size(2)-0.5];
% 
% sFP.plt.stimRow = [2:sFP.wrld.size(1)-2];
% sFP.plt.stimCol = [2:sFP.wrld.size(2)-1];
% 
% 
% % sFP.plt.pltType = 'Imagesc'
% % sFP.plt.pltType = 'Binned'
% % sFP.plt.distanceType = 'Column'
% 
% DisplActValsFun(sFP,w, Qmean(:,:,:,:,:,5));

%% More for plot, with direction preferences
% 
% % figure,
% plot(squeeze(nanmean(Qshift(3,:,3,5),[1 3])))