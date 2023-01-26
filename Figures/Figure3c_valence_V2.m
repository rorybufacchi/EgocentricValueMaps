%% MAKE SURE TO RUN Valence.m FIRST, or load this:
load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\Valence\ValenceWorkspaceForFigure.mat')

%% Plot heatmaps

fS.gridXstart = -4.5;
fS.gridXstep = 1;
fS.gridYstart = 3.5;
fS.gridYstep = 1;

f.ValenceMagnMap.f  = figure('Position',[20 20 600 400]),
f.ValenceMagnLine.f = figure('Position',[20 20 300 400])
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


allRowLags = rS(cM).s.gol.alSpR;

condN = {'Low Valence' , 'High Valence'};
plCols = {[0 0 1] , [.8 0 0]};
opts.FaceAlpha = 0.3;

for iCond = 1:2
    
    figure(f.ValenceMagnMap.f);

        ax{iCond} = subplot(1,2,iCond);
                
%         sFP.plt.intrpFact = .25;
        sFP.plt.intrpFact = 1;
       
        qTemp = nanmean(Q(:,:,:,:,:,:,iCond),[6]);
        
        DisplActValsFun(sFP,w,qTemp); hold on
        
        GridOverImage(fS,ax{iCond});
        
        caxis([0 4])
        
        colorbar off
        
        title([condN{iCond}])
        
        box off
        
        figure(f.ValenceMagnLine.f);
        
        tmpM  = nanmean(relContRel(:,iCond,:),3);
        tmpSD = nanstd(relContRel(:,iCond,:),[],3);
        opts.c = plCols{iCond};
        ShadedPlot([size(relContRel,1):-1:1]',...
            tmpM,tmpM - tmpSD,tmpM + tmpSD, opts)
%         ShadedPlot([1:11]',...
%             tmpM,tmpM - tmpSD,tmpM + tmpSD)
        
        xlabel('distance from limb')
        ylabel('Relative probability of contact-actions')

        ylim([0 1])
        xlim([1 8])

end

for iCond = 1:2
        colormap(ax{iCond},whitetocol(100,plCols{iCond}));
%         colormap(redbluecmapRory(5,5))
end

legend(condN)

%% Then individual neuron responsiveness to valence?


% A = (allSD(:,:,:,2) -  allSD(:,:,:,1)) ;
A = (allSD(:,:,:,2) -  allSD(:,:,:,1)) ./ (allSD(:,:,:,2) +  allSD(:,:,:,1));

% figure,histogram(A(:)); hold on
f.VarIncrease.f = figure('Position',[20 20 300 400])
histogram(A(:),linspace(-1,1,5),'Normalization','probability'); hold on
yLims = ylim;
plot([0 0],yLims,'-.k')

title('Most neurons increase their variability when exposed to a more valuable stimulus')

xlabel('Change in neuron variability')
ylabel('probability')

box off

% B = permute(allSD ./ (allSD(:,:,:,2) +  allSD(:,:,:,1)),[4 1 2 3]);
% figure,
% bar(nanmean(B(:,:),2))
% violin(B(:,:))



% %%
% 
% figure,
% for iL = 1:4
%     subplot(1,4,iL)
%     
%     
% %     plot([-1.1 1.1],[.1 .1],'-.k'); hold on
% %     plot([-1.1 1.1],[.2 .2],'-.k')
%     
%     
%     A2 = A(iL,:,:); A2 = A2(:);
% %     histogram(A2,[-1:0.25:1]);
% %     histogram(A2,[-1:0.125:1]);
% % % %     histogram(A2,'Normalization','probability'); hold on
%     histogram(A2,linspace(-1,1,11),'Normalization','probability'); hold on
% %     xlim([-1 1])
% ylabel('probability')
% 
%     plot([0 0],[0 .3],'-.k')
%     
% end
% xlabel('Proportion change in neuron variability')
% sgtitle('Most neurons in later layers increase their variability when exposed to a more valuable stimulus')


%% Save figures


allFields = fields(f);
for iF = 1:length(allFields)
    cF = allFields{iF};
    
    set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\Valence\BitsAndPieces\FromMatlab\' cF '.eps'] , 'epsc')
    saveas(f.(cF).f,['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\Valence\BitsAndPieces\FromMatlab\' cF '.pdf'] , 'pdf')
end




