function [netAR,s] = PlotConnections(net,s,w)
%[netAR,s] = PlotConnections(net,s,w)
%   Plots the network and checks whether there are subnetworks
%
%   Inputs:
%   - s.nta.comparison      The type of sub-networks looked for.
%                           Options: 'Valence', 'BodyPart'
%                           $$$ I SHOULD ADD IN 'None'


% % % load('D:\DPPS\DefenseAgent\Results\Net_GoalAndThreat_MultSpeeds_Goal2_Thr-4_Stay-0_1_400Retrains')
% % % load('D:\DPPS\DefenseAgent\Results\Net_Goal_3_2_Rew_BodyAndHand_BodyMoves_InfLR.mat')

% % % s.plt.stimRow=1:s.wrld.size(1)-1;
s.plt.stimRow=1:s.wrld.size(1)-3;
s.plt.stimCol=1:s.wrld.size(2)-1;
s.plt.meanLimbCols=1;

if s.plt.ON == 1
figure,
end
% Calculate the correlations
s=DefaultSettings(s); s.plt.pltType='Binned'; s.plt.sequentialLimbCols=0;

% $$$ I NEED TO CHECK WHETHER WEIGHTS ARE CALCULATED CORRECTLY AND IN THE
% RIGHT ORDER!! $$$

% $$ AND THEN make a weight thing which goes further back?

% $$$ ALSO ALSO: for the body/limb preference: check whether I've been
% using correlation or covariance or something else. Consider making the
% option for correlation, if necessary..

switch s.nta.comparison
    case 'Valence'
        s2=s;
        s2.plt.ON = 0;
        s2.plt.lmbCol=2:s.wrld.size(2)-1;
        s2.plt.plotThrFl=0;
        [Q,allNeurAct] = CalcNetOutput(s2,w,net); [rGl pGl covGl] = NetAnalysis(s2,w,allNeurAct); % Goal correlation
        glAct=squeeze(nanmean(allNeurAct(:,:,:,s.plt.lmbCol,s.plt.startLayer:s.plt.stopLayer,:),[1 2 3 4]));
        s2.plt.plotThrFl=1;
        [Q,allNeurAct] = CalcNetOutput(s2,w,net); [rThr pThr covThr] = NetAnalysis(s2,w,allNeurAct); % Threat correlation
        thrAct=squeeze(nanmean(allNeurAct(:,:,:,s.plt.lmbCol,s.plt.startLayer:s.plt.stopLayer,:),[1 2 3 4]));
    case 'BodyPart'
        testBodyCols=1:size(w.world2D,2);
        for iBdyCol=1:length(testBodyCols)
            s2=s;
            s2.plt.ON = 0;
            s2.plt.bdyCol=testBodyCols(iBdyCol);
            [Q,allNeurAct_PlusBody(:,:,:,:,:,:,iBdyCol)] = CalcNetOutput(s2,w,net);
        end
        
        allNeurAct_MeanBody = nanmean(allNeurAct_PlusBody,7);
        allNeurAct_MeanLimb = permute(nanmean(allNeurAct_PlusBody,4),[1 2 3 7 5 6 4]);
        
        s2.plt.lmbRow = size(w.world2D,1)-2;
        [rGl pGl covGl] = NetAnalysis(s2,w,allNeurAct_MeanBody ); % limb correlations
        glAct=squeeze(nanmean(allNeurAct_MeanBody (:,:,:,:,s.plt.startLayer:s.plt.stopLayer,:),[1 2 3 4]));
        
        s2.plt.lmbRow = size(w.world2D,1)-1;
        [rThr pThr covThr] = NetAnalysis(s2,w,allNeurAct_MeanLimb); % body correlation
        thrAct=squeeze(nanmean(allNeurAct_MeanLimb(:,:,:,:,s.plt.startLayer:s.plt.stopLayer,:),[1 2 3 4]));
end
% close all

%%

[GL_THR_diff GL_THR_rat GLorTHR]=DiffRat(rGl,rThr);
% $$$ covVal here
[GL_THR_diff_cov GL_THR_rat_cov GLorTHR_cov]=DiffRat(covGl,covThr);

% $$$ TEMP:
GL_THR_diff = GL_THR_diff_cov;
GL_THR_rat = GL_THR_rat_cov;
GLorTHR = GLorTHR_cov; 
rGl=covGl; rThr=covThr;

% Set everything before startlayer to NaN
if s.plt.startLayer > 1
    rGl(1 : s.plt.startLayer - 1,:)         = NaN;
    rThr(1 : s.plt.startLayer - 1,:)        = NaN;
    GLorTHR(1 : s.plt.startLayer - 1,:)     = NaN;
    GL_THR_diff(1 : s.plt.startLayer - 1,:) = NaN;
    GL_THR_rat(1 : s.plt.startLayer - 1,:)  = NaN;
end

% Set everything after stoplayer to NaN
rGl(s.plt.stopLayer + 1 : end,:)            = NaN;
rThr(s.plt.stopLayer + 1 : end,:)           = NaN;
GLorTHR(s.plt.stopLayer + 1 : end,:)        = NaN;
GL_THR_diff(s.plt.stopLayer + 1 : end,:)    = NaN;
GL_THR_rat(s.plt.stopLayer + 1 : end,:)     = NaN;


tpS=size(GL_THR_diff); % temp size
% this is a measure of how much the feeding in neurons 'prefer' goals or threats
wghtByGl=nan(tpS); wghtByThr=nan(tpS); wghtByDiff=nan(tpS); wghtByRat=nan(tpS); wghtByClass=nan(tpS);

x=[];y=[]; % Neuron plot positions
src=[];% source
tgt=[];% target
wht=[];% weights
prvS=0;
nNames=[]; % Neuron Names
for iL=s.plt.startLayer:s.plt.stopLayer
    x=[x, 1:net.layers{iL}.dimensions ];
    y=[y, iL*ones(1,net.layers{iL}.dimensions)];
    
    tmpNms= [repmat(iL,[1 net.layers{iL}.dimensions]); [1:net.layers{iL}.dimensions]];
    nNames=[nNames, sprintfc('%d_%d',tmpNms(:)' ) ];
    
    if iL<s.plt.stopLayer % No connections from the last layer forward
        % normalised weights by input to neuron
        wNorm=abs((net.LW{iL+1,iL}-mean(net.LW{iL+1,iL},2))./std(net.LW{iL+1,iL},[],2));
        % sorted weights by input to neuron
        [dmy,I]=sort(net.LW{iL+1,iL},2);
        [dmy,wSort]=sort(I,2);
        wSort=1+((wSort-min(wSort,[],2))./max(wSort,[],2));
        for iRow=1:net.layers{iL+1}.dimensions
            for iCol =1:net.layers{iL}.dimensions
                src=[src iCol+prvS];
                tgt=[tgt iRow+net.layers{iL}.dimensions+prvS];
                
                switch s.nta.useNetWeight
                    case 'Raw',     
                        wht=[wht abs(net.LW{iL+1,iL}(iRow,iCol))];
                        whtMat(iL,iRow,iCol) = abs(net.LW{iL+1,iL}(iRow,iCol));
                    case 'Norm',    
                        wht=[wht wNorm(iRow,iCol)];
                        whtMat(iL,iRow,iCol) = wNorm(iRow,iCol);
                    case 'Ranked',  
                        wht=[wht wSort(iRow,iCol)];
                        whtMat(iL,iRow,iCol) = wSort(iRow,iCol);
                end
            end
            
            
        end
        
        % This calculates the sum of the weights multiplied by the 'goalness'
        % For every neuron in layer iL+1 relative to the neurons in the
        % layer below, i.e. iL
        switch s.nta.useNetWeight
            case 'Raw',     useWght=abs(net.LW{iL+1,iL});
            case 'Norm',    useWght=wNorm;
            case 'Ranked',  useWght=wSort;
        end

        wghtByGl(iL+1,1:size(net.LW{iL+1,iL},1))=sum( ... 
            abs(useWght) .* ... Weights between layers
            GLorTHR(iL,1:size(useWght,2))  ... == 1 if neuron prefers goals
            ,2);
        wghtByThr(iL+1,1:size(useWght,1))=sum( ...
            abs(useWght) .* ...
            ~GLorTHR(iL,1:size(useWght,2)) ... == 1 if neuron prefers threats
            ,2);
        wghtByClass(iL+1,1:size(useWght,1))=sum( ...
            abs(useWght) .* ...
            (GLorTHR(iL,1:size(useWght,2))-.5).*2 ... == -1 neuron prefers threats, 1 if neuron prefers goals
            ,2);
        wghtByDiff(iL+1,1:size(useWght,1))=sum( ... 
            abs(useWght) .* ... Weights between layers
            GL_THR_diff(iL,1:size(useWght,2))  ... == +ve if neuron prefers goals
            ,2);
        wghtByRat(iL+1,1:size(useWght,1))=sum( ... 
            abs(useWght) .* ... Weights between layers
            GL_THR_rat(iL,1:size(useWght,2))  ... > 1 if neuron prefers goals
            ,2);
               
        prvS=max(src);
    end
end


% Calc and plot the correlation between neurons' preference and the previous layers'
% preferences
[difdifR difdifP]=corrcoef(GL_THR_diff(~isnan(wghtByDiff)),wghtByDiff(~isnan(wghtByDiff)));
[ratratR ratratP]=corrcoef(GL_THR_rat(~isnan(wghtByRat)),wghtByRat(~isnan(wghtByRat)));
% % % [GlGlR GlGlP]=corrcoef(GLorTHR(~isnan(GLorTHR)),wghtByClass(~isnan(GLorTHR)));
[difclassR difclassP]=corrcoef(GL_THR_diff(~isnan(wghtByDiff)),wghtByClass(~isnan(wghtByDiff)))

% % % figure,
% % % subplot(2,2,1)
% % % plot(GL_THR_diff(:),wghtByClass(:),'x');
% % % lsline; title(['Difference plot. r=' num2str(GlGlR(2)) '. p=' num2str(GlGlP(2))])
% % % xlabel('Neuron gl thr difference'); ylabel('weights*prev layer neuron difference');



% if s.plt.ON == 1
    figure,
    subplot(2,2,1)
    plot(GL_THR_diff(:),wghtByDiff(:),'x');
    lsline; title(['Difference plot. r=' num2str(difdifR(2)) '. p=' num2str(difdifP(2))])
    xlabel('Neuron gl thr difference'); ylabel('weights*prev layer neuron difference');
    
    subplot(2,2,2)
    plot(GL_THR_rat(:),wghtByRat(:),'x');
    
    lsline; title(['Ratio plot. r=' num2str(ratratR(2)) '. p=' num2str(ratratP(2))])
    xlabel('Neuron gl thr ratio'); ylabel('weights*prev layer neuron ratio');
    
    subplot(2,2,3)
%     boxplot([wghtByGl(GLorTHR==1),wghtByThr(GLorTHR==1)],'labels',{'WbGl','WbThr'}); hold on
%     plot([wghtByGl(GLorTHR==1),wghtByThr(GLorTHR==1)]')
%     title('Goal preferring neurons');
    plot(GL_THR_diff(:),wghtByClass(:),'x'); hold on
    lsline; title(['Difference class plot. r=' num2str(difclassR(2)) '. p=' num2str(difclassP(2))])
    xlabel('Neuron gl thr difference'); ylabel('weights*prev layer neuron class');
    
    
    subplot(2,2,4)
    boxplot([wghtByGl(GLorTHR==0),wghtByThr(GLorTHR==0)],'labels',{'WbGl','WbThr'}); hold on
    plot([wghtByGl(GLorTHR==0),wghtByThr(GLorTHR==0)]')
    title('Threat preferring neurons');
% end


G=digraph(src,tgt,wht);

% $$$ ALSO, can still try plotting with direct weighting, just to see what
% it looks like

% $$$$ NEED TO also make 2 subgraphs, which have weights 1/w, and which
% only include either goal or threat preferring neurons
GL_THR_diff_fl = NanRemFlatten((GL_THR_diff)');
GL_THR_rat_fl  =  NanRemFlatten((GL_THR_rat)');
GLorTHR_fl     =  NanRemFlatten((GLorTHR)');

% Make a graph with inverse weights for calulating weighted distances
Ginv=digraph(src,tgt,1./wht);
% Calculate weighted distances
dInv=distances(Ginv,1:max([src tgt]));
dInv(isinf(dInv))=NaN; dInv(1:1+size(dInv,1):end)=NaN;
% weighted distances FROM all goal and threat neurons. 
dGl=dInv; dThr=dInv;
dGl(~GLorTHR_fl,:)=NaN; 
dThr(~~GLorTHR_fl,:)=NaN;
dGlMn=nanmedian(dGl,1);
dThrMn=nanmedian(dThr,1);

% Correlation between median path length leading to neuron (from threat and
% goal neurons), and the amount of threat/goalness of a neuron
distDiff=dThrMn-dGlMn;
tmpBl=~isnan(dThrMn);
[difdistR difdistP]=corrcoef(GL_THR_diff_fl(tmpBl),distDiff(tmpBl));
[ddPl] = MakeEvenMatrix({distDiff(GLorTHR_fl==1),distDiff(GLorTHR_fl==0)});
[classdistP h stats] = ranksum(ddPl(1,:),ddPl(2,:));
classdistZ = stats.zval;
% if s.plt.ON == 1
    figure,
    subplot(1,2,1)
    plot(GL_THR_diff_fl(tmpBl),distDiff(tmpBl),'x')
    title(['Difference vs gl thr dist diff. r=' num2str(difdistR(2)) '. p=' num2str(difdistP(2))]);
    lsline
    xlabel('Gl Thr pref dif')
    ylabel('median goal dist - median thr dist');
    
    subplot(1,2,2)
    violin(ddPl','xlabel',{'GL neur ','THR neur'}); hold on; legend off
    ylabel('median goal dist - median thr dist');
    title(['RS p= ' num2str(classdistP)]);
% end


% Display the network
if s.plt.ON == 1
    p=DisplayNetwork(s,G,x,y,nNames,rGl,rThr);
    title('Goal/threat COVARIANCE colours');
end

% p=DisplayNetwork(G,x,y,nNames,glAct,thrAct);
% title('Goal/threat ACTIVATION colour');

% network Analysis Results
netAR.G         = G;
netAR.rGl       = rGl;          netAR.rThr   =rThr;
netAR.pGl       = pGl;          netAR.pThr   =pThr;
netAR.A_B_diff  = GL_THR_diff; 
netAR.A_B_rat   = GL_THR_rat;
netAR.AorB      = GLorTHR;
netAR.wht       = wht;
netAR.whtMat    = whtMat;

netAR.difdifR   = difdifR(2);
netAR.ratratR   = ratratR(2);
netAR.difdistR  = difdistR(2);
netAR.classdistZ= classdistZ;
netAR.difdifP   = difdifP(2);
netAR.ratratP   = ratratP(2);
netAR.difdistP  = difdistP(2);
netAR.classdistP= classdistP;



function [p] = DisplayNetwork(s,G,x,y,nNames,rA,rB)
    % rA is correlation/cov with A, rB is correlation/cov with B

% LWidths = abs(3*G.Edges.Weight/max(abs(G.Edges.Weight))).^2.5;
% LWidths = abs(5*G.Edges.Weight/max(abs(G.Edges.Weight)));
% LWidths = abs(7.5*G.Edges.Weight/max(abs(G.Edges.Weight)));
LWidths = abs(1.3*G.Edges.Weight/max(abs(G.Edges.Weight))).^6;
LWidths = abs(1.3*G.Edges.Weight/max(abs(G.Edges.Weight))).^9;

figure,
% p=plot(G,'XData',x,'YData',y,'NodeLabel',nNames,'LineWidth',LWidths,'EdgeColor','k');
% p=plot(G,'Layout','force','NodeLabel',nNames,'LineWidth',LWidths)
% p=plot(G,'Layout','force','WeightEffect','inverse','XStart',x,'YStart',y,'NodeLabel',nNames,'LineWidth',LWidths)

p=plot(G,'Layout','force','WeightEffect','inverse','UseGravity','on','NodeLabel',nNames,'LineWidth',LWidths);
p.MarkerSize=10e-10;

% Plot red dot for threat correlation larger than goal correlation, and
% blue dot for goal correlation larger than threat correlation
hold on
% First remove NaNs and make correlation values fit though


[A_B_diff A_B_rat AorB]=DiffRat(rA,rB);


% Set everything before startlayer to NaN
if s.plt.startLayer > 1
    rA(1 : s.plt.startLayer - 1,:)  = NaN;
    rB(1 : s.plt.startLayer - 1,:) = NaN;
    A_B_diff(1 : s.plt.startLayer - 1,:)   = NaN;
    A_B_rat(1 : s.plt.startLayer - 1,:)    = NaN;
    AorB(1 : s.plt.startLayer - 1,:)       = NaN;
end

% Set everything after stoplayer to NaN
rA(s.plt.stopLayer + 1 : end,:)         = NaN;
rB(s.plt.stopLayer + 1 : end,:)         = NaN;
A_B_diff(s.plt.stopLayer + 1 : end,:)   = NaN;
A_B_rat(s.plt.stopLayer + 1 : end,:)    = NaN;
AorB(s.plt.stopLayer + 1 : end,:)       = NaN;



A_B_diff = NanRemFlatten((A_B_diff)');
A_B_rat =  NanRemFlatten((A_B_rat)');
AorB =  NanRemFlatten((AorB)');

plot(p.XData(~AorB),p.YData(~AorB),'g.','MarkerSize',35);
plot(p.XData(~~AorB),p.YData(~~AorB),'b.','MarkerSize',35); % 'not not' to make it a boolean

for iNode=1:size(G.Nodes,1)
    %     plot(p.XData(iNode),p.YData(iNode),'.','markerfacecolor',...
    %         [0.5+A_B_diff(iNode) 0 0],'MarkerSize',25)
    
    % plot(p.XData(iNode),p.YData(iNode),'.','markerfacecolor',...
    %         [A_B_rat(iNode) A_B_rat(iNode) A_B_rat(iNode)],'MarkerSize',25)
    
    %     scatter(p.XData(iNode),p.YData(iNode),1500,[1-A_B_rat(iNode) 1-A_B_rat(iNode) 1-A_B_rat(iNode)], '.')
    
    %     scatter(p.XData(iNode),p.YData(iNode),1500,[1-A_B_rat(iNode) 0 A_B_rat(iNode)], '.')
    
    if strcmp(s.nta.compScale,'Graded')
        scatter(p.XData(iNode),p.YData(iNode),1500,[0 1-A_B_rat(iNode) A_B_rat(iNode)], '.')
    end
    
end
% scatter(p.XData,p.YData,150,'k', 'o')

end


function [AB_diff AB_rat AorB]=DiffRat(A,B)
%calculate the difference and or ratio between quantities A and B, remove
%NaNs and flatten

% make sure AorB has the right NaNs
AorB=nan(size(A));
tmpLog=~isnan(A) & ~isnan(B);
AorB(tmpLog)=abs(A(tmpLog))>abs(B(tmpLog));

AB_diff=abs(A)-abs(B);
AB_rat=(((abs(A)-abs(B))./(abs(A)+abs(B)))+1)./2;
end

function [outVar]= NanRemFlatten(inVar)
% Removes Nans and flattens
outVar=inVar(:);
outVar(isnan(outVar))=[];
end

function [neatMatrix] = MakeEvenMatrix(inputVectors)
    maxL=1;
    nVecs=length(inputVectors);
    for iV=1:nVecs
        tmpL=length(inputVectors{iV});
        if tmpL>maxL
            maxL=tmpL;
        end
    end
    
    neatMatrix=nan([nVecs,maxL]);
    
    for iV=1:nVecs
        vecL=length(inputVectors{iV});
        neatMatrix(iV,1:vecL)=inputVectors{iV};
    end
    
    
    
end


end



