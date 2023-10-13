% $$$ NEXT is make model type 5 use the network for plotting :) $$$

load('Results\ForFigures\Fig1_Results_v3')
s=rS(end).s;
w=rS(end).w;
addpath(genpath('Scripts\EgocentricValueMaps'))


% Make 2 plots, one with original gammas, and one with gamma == 1
allGammas = [0.7 0.7 ; 1 1; 0.3 0.7]

for iFig = 1:3

% settings for plot
sFP=s;
sFP.plt.lmbRow = s.wrld.size(1)-2;
sFP.plt.rowLims=[6.5 s.wrld.size(1)-1.5];
sFP.plt.colLims=[3.5 s.wrld.size(2)-3.5];
% % % sFP.plt.rowLims=[6.5 s.wrld.size(1)-1.5];
% % % sFP.plt.colLims=[3.5 s.wrld.size(2)-3.5];
sFP.plt.cBarFl=0;
sFP.plt.meanLimbCols=1;
sFP=DefaultSettings(sFP);
sFP.plt.axesVis=0;
            
% figure settings
fS.w = 32; % Width (cm)
fS.h = 15; % Height
fS.x = 5; % x screen position
fS.y = 5; % y screen position


fS.fillW = .99; % fraction of width that's going to be filled
fS.fillH = .99; % fraction of height that's going to be filled

fS.sub.w = [0.09 0.18 0.18 0.18 0.18]; % fraction of space that each subplot takes (width)
fS.sub.h = [0.07 0.35 0.35]; % fraction of space that each subplot takes (height)

% % fS.sub.w = [0.12 0.17 0.17 0.17 0.17]; % fraction of space that each subplot takes (width)
% % fS.sub.h = [0.1 0.33 0.33]; % fraction of space that each subplot takes (height)

fS.gridXstart = -4.5;
fS.gridXstep = 1;
fS.gridYstart = 3.5;
fS.gridYstep = 1;


nWpl = length(fS.sub.w); % number of plots in width
nHpl = length(fS.sub.h); % number of plots in height
wSpace = (fS.fillW-sum(fS.sub.w))./(nWpl-1) ; % fraction of space between subplots (width)
hSpace = (fS.fillH-sum(fS.sub.h))./(nHpl-1) ; % space between subplots (height)

%subplot positions
spX = MakeSubplotPos(fS.fillW, fS.sub.w, wSpace);
spY = flip(MakeSubplotPos(fS.fillH, flip(fS.sub.h), hSpace));

fig{iFig} = figure('Units', 'centimeters', 'Position', [fS.x fS.y fS.w fS.h ]);

%  ------------------------------------------------------------------------
% Top 4 plots
spl = subplot(40,1000,1);
set(spl, 'Position', [spX(2) spY(2) fS.sub.w(2)+wSpace/4 fS.sub.h(2) ]);
sFP.plt.plAct=1;
sFP.plt.lmbCol=8;

s.clc.startSC = 12;
s.clc.gammaVal = allGammas(iFig,1);
s.clc.randSpread =     [0];
s.clc.spreadProb = 1 ./[1];
s.clc.actConsequence = [0];
[newQ optQ] = CalcQDirect(s);
DisplActValsFun(sFP,w,newQ); hold on
GridOverImage(fS,spl);
% % colorbar

spl = subplot(40,1000,1)
set(spl, 'Position', [spX(3)-wSpace/4 spY(2) fS.sub.w(3)+wSpace/4 fS.sub.h(2) ]);
s.clc.gammaVal = allGammas(iFig,2);
[newQ optQ] = CalcQDirect(s);
DisplActValsFun(sFP,w,newQ); hold on
GridOverImage(fS,spl);
% % colorbar

spl = subplot(40,1000,1)
set(spl, 'Position', [spX(4) spY(2) fS.sub.w(4)+wSpace/4 fS.sub.h(2) ]);
s.clc.gammaVal = allGammas(iFig,1);
s.clc.actConsequence = [-1 0 1];
[newQ optQ] = CalcQDirect(s);
% Create average value across actions
if iFig == 3
    tmpNewQ = repmat(newQ,[1 1 1 1 3]);
    tmpNewQ(:,:,:,2:end  ,1)   = newQ(:,:,:,1:end-1);
    tmpNewQ(:,:,:,1:end-1,3) = newQ(:,:,:,2:end);
    newQ = nanmean(tmpNewQ,5);
    caxis([0 0.7]);
end
DisplActValsFun(sFP,w,newQ); hold on
GridOverImage(fS,spl);
% % colorbar

spl = subplot(40,1000,1)
set(spl, 'Position', [spX(5)-wSpace/4 spY(2) fS.sub.w(5)+wSpace/4 fS.sub.h(2) ]);

s.clc.gammaVal = allGammas(iFig,2);
[newQ optQ] = CalcQDirect(s);
% Create average value across actions
if iFig == 3
    tmpNewQ = repmat(newQ,[1 1 1 1 3]);
    tmpNewQ(:,:,:,2:end  ,1)   = newQ(:,:,:,1:end-1);
    tmpNewQ(:,:,:,1:end-1,3) = newQ(:,:,:,2:end);
    newQ = nanmean(tmpNewQ,5);
    caxis([0 0.7]);
end
DisplActValsFun(sFP,w,newQ); hold on
GridOverImage(fS,spl);
% % colorbar

%%%  ------------------------------------------------------------------------



% Bottom 4 plots
spl = subplot(40,1000,1)
set(spl, 'Position', [spX(2) spY(3) fS.sub.w(2)+wSpace/4 fS.sub.h(3) ]);

s.clc.startSC = 12;
s.clc.gammaVal = allGammas(iFig,1);
s.clc.randSpread =     [-1 0 1];
s.clc.spreadProb = 1 ./[ 3 3 3];
s.clc.actConsequence = [0];
[newQ optQ] = CalcQDirect(s);
DisplActValsFun(sFP,w,newQ); hold on
GridOverImage(fS,spl);
% % colorbar


spl = subplot(40,1000,1)
set(spl, 'Position', [spX(3)-wSpace/4 spY(3) fS.sub.w(3)+wSpace/4 fS.sub.h(2) ]);
s.clc.gammaVal = allGammas(iFig,2);
[newQ optQ] = CalcQDirect(s);
DisplActValsFun(sFP,w,newQ); hold on
GridOverImage(fS,spl);
% % colorbar

spl = subplot(40,1000,1)
set(spl, 'Position', [spX(4) spY(3) fS.sub.w(4)+wSpace/4 fS.sub.h(2) ]);
s.clc.gammaVal = allGammas(iFig,1);
s.clc.actConsequence = [-1 0 1];    
[newQ optQ] = CalcQDirect(s);
% Create average value across actions
if iFig == 3
    tmpNewQ = repmat(newQ,[1 1 1 1 3]);
    tmpNewQ(:,:,:,2:end  ,1)   = newQ(:,:,:,1:end-1);
    tmpNewQ(:,:,:,1:end-1,3) = newQ(:,:,:,2:end);
    newQ = nanmean(tmpNewQ,5);
    caxis([0 0.7]);
end
DisplActValsFun(sFP,w,newQ); hold on
GridOverImage(fS,spl);
% % colorbar


spl = subplot(40,1000,1)
set(spl, 'Position', [spX(5)-wSpace/4 spY(3) fS.sub.w(5)+wSpace/4 fS.sub.h(2) ]);
s.clc.gammaVal = allGammas(iFig,2);
s.clc.actConsequence = [-1 0 1];
[newQ optQ] = CalcQDirect(s);
% Create average value across actions
if iFig == 3
    tmpNewQ = repmat(newQ,[1 1 1 1 3]);
    tmpNewQ(:,:,:,2:end  ,1)   = newQ(:,:,:,1:end-1);
    tmpNewQ(:,:,:,1:end-1,3) = newQ(:,:,:,2:end);
    newQ = nanmean(tmpNewQ,5);
    caxis([0 0.7]);
end
DisplActValsFun(sFP,w,newQ); hold on
GridOverImage(fS,spl);
% % colorbar


%  ------------------------------------------------------------------------
% gamma labels
glPosx = [3 2 5 4 3 2 5 4]; % gamma label x-positions
wST= (wSpace/4) .*[-1 0 -1 0 -1 0 -1 0]; % width shift temp
glPosy = [2 2 2 2 3 3 3 3]; % gamma label y-positions
for iGP = 1: length(glPosx)% gamma position
    
    gm=rS(iGP).s.lp.gamma;
    tmpH=(nHpl-2)./(nHpl-1); % extra space needed for the label
    spl = subplot(40,1000,1)
    set(spl, 'Position', [spX(glPosx(iGP))+wST(iGP) spY(glPosy(iGP))+fS.sub.h(glPosy(iGP))+hSpace.*.14 fS.sub.w(glPosx(iGP))+wSpace/4 fS.sub.h(1)/2 ]);
    fill([-1 1 1 -1],[-1 -1 1 1],0.5*[1 1 1],'EdgeColor','none')
    set(gca,'visible','off')
    FNT = {'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'Color', 1.*[1 1 1], };
    txt=text(spl,0, 0, ['\gamma = ' num2str(gm)],'VerticalAlignment', 'middle',  FNT{:});
end

%  ------------------------------------------------------------------------
% Action group 1 label
% extra space needed for the label
spl = subplot(40,1000,1)
set(spl, 'Position', [spX(2) spY(1)-hSpace/3 fS.sub.w(2)+fS.sub.w(3)+wSpace fS.sub.h(1)+hSpace/3 ]);
fill([-1 1 1 -1],[-1 -1 1 1],0.5*[1 1 1],'EdgeColor','none')
set(gca,'visible','off')
FNT = {'HorizontalAlignment', 'left', ...
    'FontSize', 14, 'FontWeight', 'bold','Color', 1.*[1 1 1], };
% txt=text(spl,-.9, 0, 'Available actions:','VerticalAlignment', 'bottom',  FNT{:});
txt=text(spl,-.9, 0, 'Stay','VerticalAlignment', 'middle',  FNT{:});


%  ------------------------------------------------------------------------
% Action group 2 label
spl = subplot(40,1000,1)
set(spl, 'Position', [spX(4) spY(1)-hSpace/3 fS.sub.w(3)+fS.sub.w(4)+wSpace fS.sub.h(1)+hSpace/3 ]);
fill([-1 1 1 -1],[-1 -1 1 1],0.5*[1 1 1],'EdgeColor','none')
set(gca,'visible','off')
FNT = {'HorizontalAlignment', 'left', ...
    'FontSize', 14, 'FontWeight', 'bold','Color', 1.*[1 1 1], };
% txt=text(spl,-.9, 0, 'Available actions:','VerticalAlignment', 'bottom',  FNT{:});
txt=text(spl,-.9, 0, 'Left, Stay, Right','VerticalAlignment', 'middle',  FNT{:});

%  ------------------------------------------------------------------------
% Dynamics group 1 label
spl = subplot(40,1000,1)
set(spl, 'Position', [spX(1) spY(2) fS.sub.w(1)/3 fS.sub.h(2) ]);
fill([-1 1 1 -1],[-1 -1 1 1],0.5*[1 1 1],'EdgeColor','none')
set(gca,'visible','off')
FNT = {'Rotation', 90, ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 14, 'FontWeight', 'bold','Color', 1.*[1 1 1], };
txt=text(spl,0, 0, 'Deterministic','VerticalAlignment', 'middle',  FNT{:});


% % % Dynamics group 2 image
% % bS = 0.2; % base Size
% % bDispx = [0 0 0];
% % bDispy = bS.*1.5.*[2 1 0];
% % xLIMS=[0 bS];
% % yLIMS=[0 max(bDispy) + bS];
% % xWarpFact=(diff(xLIMS)/(fS.sub.w(1)/2))
% % 
% % spl = subplot(40,1000,1)
% % xTmp=spX(1)+(fS.sub.w(1)/2)+wSpace/2;
% % set(spl, 'Position', [xTmp spY(2) fS.sub.w(1)/2 fS.sub.h(2) ]);
% % for iB=1:length(bDispx) % box index
% %     fill(bDispx(iB) + (bS.*[0 1 1 0]), bDispy(iB) + (bS.*[0 0 1 1]),1*[1 1 1],'EdgeColor','k'); hold on
% %     if iB<length(bDispx)
% %         annotation(gcf,'arrow', xTmp+fS.sub.w(1).*[1/4 1/4],spY(2)+fS.sub.h(2).*((3/8).*(iB-1)+[(3/8) (1/4)]));
% %         text(spl,(1/16), ((5/16).*(iB-1))+(1/4), 'P=1/2','VerticalAlignment', 'middle','HorizontalAlignment', 'right');
% %     end
% % end
% % set (gca,'xlim',xLIMS);
% % set (gca,'ylim',yLIMS);
% % set(gca,'visible','off')



%  ------------------------------------------------------------------------
% Dynamics group 2 label
spl = subplot(40,1000,1)
set(spl, 'Position', [spX(1) spY(3) fS.sub.w(1)/3 fS.sub.h(3) ]);
fill([-1 1 1 -1],[-1 -1 1 1],0.5*[1 1 1],'EdgeColor','none')
set(gca,'visible','off')
FNT = {'Rotation', 90, ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 14, 'FontWeight', 'bold','Color', 1.*[1 1 1], };
txt=text(spl,0, 0, 'Probabilistic','VerticalAlignment', 'middle',  FNT{:});



% % % Dynamics group 3 image $$$ HERE HERE
% % bS = 0.2; % base Size
% % bDispx = bS.*[0 0 -1 1];
% % bDispy = bS.*1.5.*[1 0  0 0];
% % 
% % spl = subplot(40,1000,1)
% % xTmp=spX(1)+(fS.sub.w(1)/2)+wSpace/2;
% % set(spl, 'Position', [xTmp spY(3) fS.sub.w(1)/2 fS.sub.h(3) ]);
% % for iB=1:length(bDispx) % box index
% %     fill(bDispx(iB) + (bS.*[0 1 1 0]), bDispy(iB) + (bS.*[0 0 1 1]),1*[1 1 1],'EdgeColor','k'); hold on
% % end
% % 
% % annotation(gcf,'arrow', xTmp+fS.sub.w(1).*[1/4 1/4],spY(3)+fS.sub.h(2).*[3/4 1/4]);
% % text(spl,(1/16), ((5/16).*(iB-1))+(1/4), 'P=1/2','VerticalAlignment', 'middle','HorizontalAlignment', 'right');
% % 
% % set (gca,'xlim',[-bS bS]);
% % set (gca,'ylim',[0 max(bDispy) + bS]);
% % set(gca,'visible','off')


%  ------------------------------------------------------------------------
% Rewards label
spl = subplot(40,1000,1)
set(spl, 'Position', [spX(1) spY(1)+(fS.sub.h(1)*(1/3))  fS.sub.w(1)+(wSpace*(3/4)) fS.sub.h(1)*(2/3) ]);
fill([-1 1 1 -1],[-1 -1 1 1],0.5*[1 1 1],'EdgeColor','none')
set(gca,'visible','off')
FNT = {'HorizontalAlignment', 'center', ...
    'FontSize', 12, 'FontWeight', 'bold','Color', 1.*[1 1 1], };
txt=text(spl,0, 0, 'Rewards','VerticalAlignment', 'middle',  FNT{:});

colormap(whitetocol(100,[0 0 .7]))

end

%% Save figures


for iF = 1:length(fig)
    set(fig{iF}, 'Renderer', 'painters'); % default, opengl
    saveas(fig{iF},['Results\ForFigures\TheoreticalFig\BitsAndPieces\FromMatlab\BasePlots_DirectCalc' num2str(iF) '.eps'] , 'epsc')
    saveas(fig{iF},['Results\ForFigures\TheoreticalFig\BitsAndPieces\FromMatlab\BasePlots_DirectCalc' num2str(iF) '.pdf'] , 'pdf')
end

%% FUNCTIONS

function plotPos = MakeSubplotPos(fillSize, subSize, subSpace)
plotPos=(1-fillSize)./2;
for iX=2:length(subSize)
    plotPos = [plotPos plotPos(end)+subSize(iX-1) + subSpace] ;
end
end


function [newQ optQ] = CalcQDirect(s)
newQ = zeros([s.wrld.size s.wrld.size]);
optQ = zeros([s.wrld.size s.wrld.size]);

% first initialise with the contact reward
for iLC = 1:length(newQ)
    newQ(:,iLC,s.clc.startSC,iLC) = 1;
    optQ(:,iLC,s.clc.startSC,iLC) = 1;
end

% next, loop through rows backwards and add sums
for iLC = 3:length(newQ)-2
    for iSR = s.clc.startSC - 1 : -1 : 1
        
        for iSC = 3:size(newQ,4) - 2
            colsBelow = iSC + s.clc.randSpread;
            
            % For each action, take the maximum expected value that it could lead
            % to (i.e. accounting for random spread after action is taken)
            newQ(:,iLC,iSR,iSC) = s.clc.gammaVal .* ...
                    sum(s.clc.spreadProb .* ...
                    squeeze(optQ(:,iLC,iSR + 1,colsBelow )),2);
                
            for iAct = 1:length(s.clc.actConsequence)
                optQ(:,iLC,iSR,iSC) = max(optQ(:,iLC,iSR,iSC), s.clc.gammaVal .* ...
                    sum(s.clc.spreadProb .* ...
                    squeeze(optQ(:,iLC,iSR + 1,colsBelow + s.clc.actConsequence(iAct) )),2));
            end
                
        end
        
    end
end

% Then set the value of the touch condition back to 0
for iLC = 1:length(newQ)
    newQ(:,iLC,s.clc.startSC,iLC) = 0;
end

end
    

