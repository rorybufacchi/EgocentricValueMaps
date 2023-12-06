% NOTE: first run ProximityPosition.m, or load this:
% load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\WorkSpaceForProxPos.mat')
% load('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\WorkSpaceForProxPos_3.mat')

% % % addpath('D:\Old_D\DPPS\DefenseAgent\Scripts\rlsimplepps')
% % % addpath('D:\Old_D\DPPS\DefenseAgent\Scripts')
% % % addpath('D:\Old_D\Programs\Matlab\Utilities\')
% % % addpath('D:\Old_D\Programs\Matlab\Utilities\plotting\colormaps\')


% % % %% Old HBR data
% % % 
% % % f.HBR.s.LineSpec = {'LineWidth',2};
% % % f.HBR.s.xLim     = [4.5 45];
% % % f.HBR.s.yLim     = [-25 25];
% % % f.HBR.s.cLim     = [0 .5];
% % % 
% % % f.HBR.f = figure('Position',[20 20 600 600]);
% % % 
% % % % -------------------------------------------------------------------------
% % % f.HBR.ax{1} = axes('Position',[.1 .1 .6 .6]);
% % % 
% % % addpath('D:\Old_D\DPPS\Scripts\DesktopBackup\Scripts\geometry_plots')
% % % addpath('D:\Old_D\DPPS')
% % % load('D:\Old_D\DPPS\TriNeur\hitProbFor2Dall_BothHalves_noFace_noSpace.mat')
% % % 
% % % xVals=-20:0.5:80;
% % % yVals=-25:0.5:25;
% % % circSize=350;
% % % hpFor2DAv=nanmean(hpFor2Dall,3);
% % % 
% % % imagesc(xVals,yVals,nanmean(hpFor2Dall,3))
% % % xlim(f.HBR.s.xLim)
% % % ylim(f.HBR.s.yLim)
% % % 
% % % xlabel('Distance from midline (cm)')
% % % ylabel('Distance from midline (cm)')
% % % 
% % % caxis(f.HBR.s.cLim);
% % % 
% % % CustomColourMap2;
% % % colorbar off
% % % 
% % % 
% % % % -------------------------------------------------------------------------
% % % f.HBR.ax{2} = axes('Position',[.1 .75 .6 .15]);
% % % inclX = xVals >= 9.5 + 3.5 & xVals < f.HBR.s.xLim(2);
% % % plot(xVals(inclX),movmean(hpFor2DAv(floor(end/2),inclX),5),f.HBR.s.LineSpec{:});
% % % xlim(f.HBR.s.xLim)
% % % f.HBR.ax{2}.XAxis.Visible = 'off';
% % % f.HBR.ax{2}.Box = 'off';
% % % ylabel('HBR magnitude')
% % % 
% % % % -------------------------------------------------------------------------
% % % f.HBR.ax{3} = axes('Position',[.75 .1 .15 .6]);
% % % plot(movmean(hpFor2DAv(:,xVals == 13.5),5), yVals,f.HBR.s.LineSpec{:});
% % % ylim(f.HBR.s.yLim)
% % % f.HBR.ax{3}.YAxis.Visible = 'off';
% % % f.HBR.ax{3}.Box = 'off';
% % % xlabel('HBR magnitude')
% % % 
% % % % -------------------------------------------------------------------------
% % % f.HBR.ax{4} = axes('Position',[.6 .15 .1 .5]);
% % % f.HBR.cb = colorbar('Color','w')
% % % f.HBR.ax{4}.Visible = 'off'
% % % f.HBR.cb.Position = f.HBR.cb.Position + [0 0 0.025 0];
% % % f.HBR.cb.Ticks = [0 1];
% % % f.HBR.cb.TickLabels = {'Min','Max'};
% % % 
% % % 
% % % % -------------------------------------------------------------------------
% % % f.HBR.ax{5} = axes('Position',[.025 .325 .15 .15]);
% % % tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\HeadTemplate.jpg');
% % % alphaChan = ~all(tmpFig > 250, 3);
% % % image(tmpFig,'AlphaData',alphaChan); axis off



%% Old HBR data ROTATED

f.HBR.s.LineSpec = {'LineWidth',2};
f.HBR.s.xLim     = [4.5 45];
f.HBR.s.yLim     = [-25 25];
f.HBR.s.cLim     = [0 .5];

f.HBR.f = figure('Position',[20 20 600 600]);

% -------------------------------------------------------------------------
f.HBR.ax{1} = axes('Position',[.1 .1 .6 .6]);

addpath('D:\Old_D\DPPS\Scripts\DesktopBackup\Scripts\geometry_plots')
addpath('D:\Old_D\DPPS')
load('D:\Old_D\DPPS\TriNeur\hitProbFor2Dall_BothHalves_noFace_noSpace.mat')

xVals=-20:0.5:80;
yVals=-25:0.5:25;
circSize=350;
hpFor2DAv=nanmean(hpFor2Dall,3);

imagesc(yVals,xVals,nanmean(hpFor2Dall,3)'); axis xy
xlim(f.HBR.s.yLim)
ylim(f.HBR.s.xLim)

ylabel('Distance from midline (cm)')
xlabel('Distance from midline (cm)')

caxis(f.HBR.s.cLim);

% CustomColourMap2;
colormap(whitetocol(100,[0 0 1]));
colorbar off


% -------------------------------------------------------------------------
f.HBR.ax{2} = axes('Position',[.1 .75 .6 .15]);
plot(yVals,movmean(hpFor2DAv(:,xVals == 13.5),5),f.HBR.s.LineSpec{:});
xlim(f.HBR.s.yLim)
f.HBR.ax{2}.XAxis.Visible = 'off';
f.HBR.ax{2}.Box = 'off';
ylabel('HBR magnitude')

% -------------------------------------------------------------------------
f.HBR.ax{3} = axes('Position',[.75 .1 .15 .6]);
inclX = xVals >= 9.5 + 3.5 & xVals < f.HBR.s.xLim(2);
plot(movmean(hpFor2DAv(floor(end/2),inclX),5),xVals(inclX),f.HBR.s.LineSpec{:});
ylim(f.HBR.s.xLim)
f.HBR.ax{3}.YAxis.Visible = 'off';
f.HBR.ax{3}.Box = 'off';
xlabel('HBR magnitude')

% -------------------------------------------------------------------------
f.HBR.ax{4} = axes('Position',[.6 .15 .1 .5]);
% f.HBR.cb = colorbar('Color','w')
f.HBR.cb = colorbar('Color','k')
f.HBR.ax{4}.Visible = 'off'
f.HBR.cb.Position = f.HBR.cb.Position + [0 0 0.025 0];
f.HBR.cb.Ticks = [0 1];
f.HBR.cb.TickLabels = {'Min','Max'};


% -------------------------------------------------------------------------
f.HBR.ax{5} = axes('Position',[.325 .025  .15 .15 ]);
tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\HeadTemplate.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off




%% HBR data equivalent

f.HBRQ.s.LineSpec = {'LineWidth',2};
f.HBRQ.f = figure('Position',[20 20 600 600]);


fS.gridXstart = -4.5;
fS.gridXstep = 1;
fS.gridYstart = 3.5;
fS.gridYstep = 1;


clear Qtmp
for iM = 1:3 %length(rS)
    s   = DefaultSettings(rS(iM).s);
    w   = rS(iM).w;
    net = rS(iM).net;
    
    
    s.plt.lmbCol=5;
    s.plt.bdyCol=10;
    s.plt.rowLims=[1.5 s.wrld.size(1)-0.5];
    [Qtmp(:,:,:,:,:,iM),allNeurAct] = CalcNetOutput(s,w,net);
end

Q = nanmean(Qtmp,6);
for iM = 1:length(rS)
    Q(:,:,:,:,iM) = nanmean(Q,5);
end

s.plt.intrpFact=1;%.2; 
s.plt.OneActFl=1;
s.plt.plAct = 2;
s.plt.meanLimbCols = 1;
s.plt.lmbCol=2:12;
s.plt.bdyCol=0;
s.plt.stimRow=[1:size(w.world2D,1)];
s.plt.stimCol=[2:size(w.world2D,1)];

% -------------------------------------------------------------------------
f.HBRQ.ax{1} = axes('Position',[.1 .1 .6 .6]);

DisplActValsFun(s,w,Q);
GridOverImage(fS,f.HBRQ.ax{1});
xlabel('Stimulus x-position'); ylabel('Stimulus y-position');
title('Network value output')

ylim([1 11])
xlim([-5 5])

% CustomColourMap2
colormap(whitetocol(100,[0 0 1]));
colorbar off

ylabel('Distance from midline (pixels)')
xlabel('Distance from midline (pixels)')


handRow = 12;
% % % handRow = ;
kAction = 1;

Qtable = rS(1).Qtable;


cM = 1;

s = rS(cM).s;
w = rS(cM).w;

iC=0; clear cD rD aD
for hC=2:15
    iC=iC+1;
    [cD(:,:,iC),rD(:,:,iC),aD(:,:,iC)] = CalcDist1D(s,w,hC,handRow);
    qTemp(:,:,iC)=squeeze(Qtable(handRow,hC,:,:,kAction));
    % % %     plot(aD(:,:,iC),qTemp(:,:,iC),'x');hold on
end


% -------------------------------------------------------------------------
f.HBRQ.ax{2} = axes('Position',[.1 .75 .6 .15]);
nBins=11; clear binQ
[Y,E] = discretize(cD(:),nBins);
for iBin=1:nBins
    binQ(iBin)=nanmean(qTemp(Y==iBin));
end
% plot((E(1:end-1)+E(2:end))./2,binQ); 
plot([-6:6],squeeze(nanmean(Q(w.lmb.row,8,10,2:14,:),5)),f.HBRQ.s.LineSpec{:})
xlabel('Distance from "limb"'); ylabel('Average action value');
box off
f.HBRQ.ax{2}.XAxis.Visible = 'off';
f.HBRQ.ax{2}.Box = 'off';



% -------------------------------------------------------------------------
f.HBRQ.ax{3} = axes('Position',[.75 .1 .15 .6]);
nBins=10; clear binQ
[Y,E] = discretize(rD(:),nBins);
for iBin=1:nBins
    binQ(iBin)=nanmean(qTemp(Y==iBin));
end
% plot(binQ,(E(1:end-1)+E(2:end))./2); 
plot(squeeze(nanmean(Q(w.lmb.row,8,1:11,8,:),5)),11:-1:1,f.HBRQ.s.LineSpec{:})
xlabel('Average action value');
box off
f.HBRQ.ax{2}.XAxis.Visible = 'off';
f.HBRQ.ax{2}.Box = 'off';


% -------------------------------------------------------------------------
f.HBRQ.ax{4} = axes('Position',[.6 .15 .1 .5]);
% f.HBRQ.cb = colorbar('Color','w')
f.HBRQ.cb = colorbar('Color','k')
f.HBRQ.ax{4}.Visible = 'off'
f.HBRQ.cb.Position = f.HBRQ.cb.Position + [0 0 0.025 0];
f.HBRQ.cb.Ticks = [0 1];
f.HBRQ.cb.TickLabels = {'Min','Max'};

% -------------------------------------------------------------------------
f.HBRQ.ax{5} = axes('Position',[.325 .025  .15 .15 ]);
tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\HeadTemplate.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off

%% Monkey Limb positions


% f.LPos.f = figure('Position',[20 20 600 400]);
f.LPos.f = figure('Position',[20 20 1200 600]);
f.LPos.s.xLim       = [3 15];
f.LPos.s.yLim       = [6 13];
f.LPos.s.limbPos    = [12 8];

% =========================================================================
% Real data

% -------------------------------------------------------------------------
f.LPos.ax{1} = axes('Position',[.1 .6 .32 .35 ]);
tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\MonkeyTopDown2Poses_V2_Narrow.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off


% -------------------------------------------------------------------------
f.LPos.ax{2} = axes('Position',[.1 .1 .32 .35 ]);
monkeyDat(:,:,1) = [0.916 1.049 1.309 1.080 0.858 ; ...
                    0.806 0.937 1.222 0.990 0.719 ; ...
                    0.701 0.816 1.135 0.896 0.576];
monkeyDat(:,:,2) = [0.736 0.740 1.007 1.562 1.125 ; ...
                    0.604 0.639 0.903 1.528 0.997 ; ...
                    0.469 0.535 0.795 1.497 0.868];
monkeyMax = 1.594;

opts.FaceAlpha = 0.3;
opts.c = [0 0 1];
ShadedPlot([1:5],monkeyDat(2,:,1),monkeyDat(3,:,1),monkeyDat(1,:,1),opts); hold on
opts.c = [1 0 0];
ShadedPlot([1:5],monkeyDat(2,:,2),monkeyDat(3,:,2),monkeyDat(1,:,2),opts); hold on



% =========================================================================
% Modelled data


fS2.gridXstart = f.LPos.s.xLim(1)+.5;
fS2.gridXstep = 1;
fS2.gridYstart = .5; %f.LPos.s.yLim(1);
fS2.gridYstep = 1;

% -------------------------------------------------------------------------
f.LPos.ax{3} = axes('Position',[.6 .8 .3 .15 ]);

s.plt.meanLimbCols = 0;
s.plt.lmbCol=8;

% DisplActValsFun(s,w,Q); hold on
% s.plt.lmbCol=8;
% DisplActValsFun(s,w,-Q);
mergeQPlot = zeros(size(Q));
mergeQPlot(:,f.LPos.s.limbPos(2),2:12,2:14,:) = Q(:,f.LPos.s.limbPos(1),2:12,2:14,:) -  Q(:,f.LPos.s.limbPos(2),4:14,2:14,:);
% mergeQPlot(:,8,3:14,2:14,:) = Q(:,11,2:13,2:14,:) -  Q(:,8,3:14,2:14,:);
% mergeQPlot(:,8,3:14,2:14,:) = Q(:,11,2:13,2:14,:) -  Q(:,6,2:13,2:14,:);
% imagesc(squeeze(nanmean( mergeQPlot(w.lmb.row,8,3:13,2:14,:) ,5)));
s.plt.lmbCol=f.LPos.s.limbPos(2);
s.plt.intrpFact = 1;%.2;
DisplActValsFun(s,w,Q); hold on
DisplActValsFun(s,w,-Q);
DisplActValsFun(s,w,mergeQPlot);

xlabel('Stimulus x-position'); ylabel('Stimulus y-position');
title('Network value output')

ylim(f.LPos.s.yLim)
xlim(f.LPos.s.xLim)

GridOverImage(fS2,f.LPos.ax{3});

colormap(redbluecmapRory(100,20))
caxis([-2 2])
% CustomColourMap2
colorbar off

f.LPos.ax{3}.Visible = 'off';


% -------------------------------------------------------------------------
f.LPos.ax{5} = axes('Position',[.6 .1 .32 .35 ]);



opts.FaceAlpha = 0.3;
opts.c = [0 0 1];
tmpM = squeeze(nanmean(Q(w.lmb.row,f.LPos.s.limbPos(2),10:11,2:14,:),[3 5]));
tmpL = tmpM - squeeze(nanstd(Q(w.lmb.row,f.LPos.s.limbPos(2),10:11,2:14,:),[],[3 5]));
tmpH = tmpM + squeeze(nanstd(Q(w.lmb.row,f.LPos.s.limbPos(2),10:11,2:14,:),[],[3 5]));
ShadedPlot([1:length(tmpM)],tmpM,tmpL,tmpH,opts); hold on
opts.c = [1 0 0];
tmpM = squeeze(nanmean(Q(w.lmb.row,f.LPos.s.limbPos(1),10:11,2:14,:),[3 5]));
tmpL = tmpM - squeeze(nanstd(Q(w.lmb.row,f.LPos.s.limbPos(1),10:11,2:14,:),[],[3 5]));
tmpH = tmpM + squeeze(nanstd(Q(w.lmb.row,f.LPos.s.limbPos(1),10:11,2:14,:),[],[3 5]));
ShadedPlot([1:length(tmpM)],tmpM,tmpL,tmpH,opts); hold on


% plot(squeeze(nanmean(Q(w.lmb.row,f.LPos.s.limbPos(2),2:12,2:14,:),[3 5])),'LineWidth',2); hold on 
% plot(squeeze(nanmean(Q(w.lmb.row,f.LPos.s.limbPos(1),2:12,2:14,:),[3 5])),'LineWidth',2); 
% % % legend('limb left','limb right')
xlabel('Stimulus x-position');
ylabel('Average action value')
xlim(f.LPos.s.xLim - [0 1])


% -------------------------------------------------------------------------
f.LPos.ax{4} = axes('Position',[.55 .6 .32 .35 ]);
tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\MonkeyTopDown2Poses_NoBubbles_Narrow.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off


%% Individual neural activity 

cM = 4;

f.NeurA.f = figure('Position',[20 20 1200 1400]);

% REAL DATA
% =========================================================================
% -------------------------------------------------------------------------
f.NeurA.ax{1} = axes('Position',[.1 .8 .2 .15 ]);
tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\BrainAreas.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off


% -------------------------------------------------------------------------
f.NeurA.ax{2} = axes('Position',[.1 .45 .25 .25 ]);
% % MonkeyDatLimb = [0.151 0.151 0.143 0.355 0.231 0.058 03420 0.145 0.145 0.051 0.145 0.145 0.232 0.117 0.056 0.559 0.432 0.138 0.237 0.537 0.329 0.548 0.435 0.644 0.851 0.770 0.340 0.224 0.546 0.434 0.645 0.544 1.269 1.591 1.374 1.072 1.480 1.374 0.967 1.061 1.164 0.967 1366 1.183 0.863 0.961 1.164 1.063 0.546 0.653 ; ...
% %                  0.151 0.242 0.132 0.132 0.049 0.118 0.058 0.243 0.131 0.058 0.338 0.121 0.121 0.050 0.232 0.346 0.240 0.047 0.863 0.135 0.749 0.329 1.067 0.445 1.156 0.340 0.641 0.953 0.432 0.227 0.522 0.329 0.741 0.660 1.158 0.743 1.374 1.072 1.48 0.746 1.48 0.967 1.061 1.686 1.379 1.076 1.270 0.937 0.648 1.171 1.062 0.950 0.539; ...
% %                  0.121 0.227 0.109 0.045 0.128 0.049 0.049 0.132 0.231 0.049 0.461 0.461 0.114 0.227 0.947 0.637 0.322 0.730 0.448 0.937 0.324 0.644 0.851 0.529 0.323 0.213 0.938 0.229 0.846 0.329 0.837 1.586 1.269 1.374 1.266 1.266 0.736 0.736 0.736 0.545 1.061 0.837 0.740 1.063 0.936 0.419]
monkeyDatLimb = [0.231 0.058 0.342 0.145 0.145 0.051 0.145 0.145 0.232 0.117 0.056 0.559 0.432 0.138 0.237 0.537 0.329 0.548 0.435 0.644 0.851 0.770 0.340 0.224 0.546 0.434 0.645 0.544 1.269 1.591 1.374 1.072 1.480 1.374 0.967 1.061 1.164 0.967 1.366 1.183 0.863 0.961 1.164 1.063 0.546 0.653 ; ...
                 0.243 0.131 0.058 0.338 0.121 0.121 0.050 0.232 0.346 0.240 0.047 0.863 0.135 0.749 0.329 1.067 0.445 1.156 0.340 0.641 0.953 0.432 0.227 0.522 0.329 0.741 0.660 1.158 0.743 1.374 1.072 1.480 0.746 1.480 0.967 1.061 1.686 1.379 1.076 1.270 0.937 0.648 1.171 1.062 0.950 0.539 ; ...
                 0.121 0.227 0.109 0.045 0.128 0.049 0.049 0.132 0.231 0.049 0.461 0.461 0.114 0.227 0.947 0.637 0.322 0.730 0.448 0.937 0.324 0.644 0.851 0.529 0.323 0.213 0.938 0.229 0.846 0.329 0.837 1.586 1.269 1.374 1.266 1.266 0.736 0.736 0.736 0.545 1.061 0.837 0.740 1.063 0.936 0.419];
monkeyDatLimb   = (monkeyDatLimb ./ 0.975) .*50;
monkeyTimes     = [1:size(monkeyDatLimb,2)] ./ 1.2 - .38;
tmpD     = nanmean(monkeyDatLimb);
tmpX     = [1:length(tmpD)] - length(tmpD)/2;
tmpGaus  = normpdf(tmpX,14,8);
imagesc([tmpGaus' * movmean(tmpD,10)]');
f.NeurA.ax{2}.XAxis.Visible = 'off';
f.NeurA.ax{2}.YAxis.Visible = 'off';
title('Real Limb-centred Neuron')
box off


% -------------------------------------------------------------------------
f.NeurA.ax{3} = axes('Position',[.375 .45 .075 .25 ]);
f.NeurA.RealBarsLimb = bar(monkeyTimes, nanmean(monkeyDatLimb),'FaceColor','flat'); hold on; box off
f.NeurA.RealBarsLimb.CData = zeros(size(f.NeurA.RealBarsLimb.CData));
plot([0 0], [0 1.8],'-.k')

view([90 90])

f.NeurA.ax{3}.XAxis.Visible = 'off';
ylabel('Activity (Hz)')



% -------------------------------------------------------------------------
% Use the monkeydatlimb, and smear it in the x-direction,
% using Imagesc
monkeyDatHead = [0.240 0.174 0.0588 0.091 0.129 0.092 0.092 0.053 0.025 0.053 0.136 0.061 0.0160 0.129 0.058 0.094 0.053 0.167 0.098 0.125 0.094 0.049 0.089 0.056 0.037 0.125 0.099 0.245 0.164 0.245 0.195 0.096 0.240 0.158 0.059 0.244 0.127 0.089 0.197 0.125 0.316 0.348 0.237 0.125 0.424 0.308 0.240 .0309 0.424 0.277 0.348 0.421 0.280 0.193 0.238 0.346 0.428 0.384 0.447 0.535 0.316 0.535 0.353 0.384 0.500 0.417 0.752 0.707 0.672 0.787 0.787 0.752 0.719 0.676 0.786 0.646 0.856 0.721 0.823 0.823 0.823 0.496 0.669 0.669 0.494 1.009 0.602 0.642 0.716 0.501 0.449 ; ...
                 0.450 0.421 0.235 0.079 0.039 0.1988 0.051 0.046 0.046 0.014 0.080 0.018 0.084 0.084 0.044 0.018 0.023 0.0880 0.047 0.124 0.013 0.084 0.084 0.051 0.086 0.047 0.047 0.047 0.084 0.047 0.047 0.084 0.047 0.120 0.080 0.125 0.088 0.125 0.059 0.117 0.082 0.112 0.086 0.160 0.087 0.113 0.164 0.087 0.120 0.311 0.193 0.235 0.120 0.271 0.165 0.125 0.369 0.304 0.408 0.374 0.412 0.417 0.231 0.458 0.369 0.523 0.523 0.627 0.598 0.667 0.530 0.744 0.528 0.738 0.698 0.959 0.823 1.112 0.856 0.954 0.924 1.042 0.967 0.608 0.667 0.601 0.648 0.860 0.431 0.745 0.415 ];
monkeyDatHead   = (monkeyDatHead ./ 0.73) .*100;
monkeyTimes     = [1:size(monkeyDatHead,2)] ./ (1/5) - .22;
f.NeurA.ax{4} = axes('Position',[.1 .1 .25 .25 ]);
tmpD     = nanmean(monkeyDatHead);
tmpX     = [1:length(tmpD)] - length(tmpD)/2;
tmpGaus  = normpdf(tmpX,-27,15);
imagesc([tmpGaus' * movmean(tmpD,10)]');
f.NeurA.ax{4}.YAxis.Visible = 'off';
title('Real Head-centred Neuron')

box off

% -------------------------------------------------------------------------
f.NeurA.ax{5} = axes('Position',[.375 .1 .075 .25 ]);
% monkeyDatHead = [0.606 0.461 0.382 0.275 0.207 0.240 0.174 0.0588 0.091 0.129 0.092 0.092 0.053 0.025 0.053 0.136 0.061 0.0160 0.129 0.058 0.094 0.053 0.167 0.098 0.125 0.094 0.049 0.089 0.056 0.037 0.125 0.099 0.245 0.164 0.245 0.195 0.096 0.240 0.158 0.059 0.244 0.127 0.089 0.197 0.125 0.316 0.348 0.237 0.125 0.424 0.308 0.240 .0309 0.424 0.277 0.348 0.421 0.280 0.193 0.238 0.346 0.428 0.384 0.447 0.535 0.316 0.535 0.353 0.384 0.500 0.417 0.752 0.707 0.672 0.787 0.787 0.752 0.719 0.676 0.786 0.646 0.856 0.721 0.823 0.823 0.823 0.496 0.669 0.669 0.494 1.009 0.602 0.642 0.716 0.501 0.449 ; ...
%                  0.450 0.421 0.235 0.079 0.039 0.1988 0.051 0.046 0.046 0.014 0.080 0.018 0.084 0.084 0.044 0.018 0.023 0.0880 0.047 0.124 0.013 0.084 0.084 0.051 0.086 0.047 0.047 0.047 0.084 0.047 0.047 0.084 0.047 0.120 0.080 0.125 0.088 0.125 0.059 0.117 0.082 0.112 0.086 0.160 0.087 0.113 0.164 0.087 0.120 0.311 0.193 0.235 0.120 0.271 0.165 0.125 0.369 0.304 0.408 0.374 0.412 0.417 0.231 0.458 0.369 0.523 0.523 0.627 0.598 0.667 0.530 0.744 0.528 0.738 0.698 0.959 0.823 1.112 0.856 0.954 0.924 1.042 0.967 0.608 0.667 0.601 0.648 0.860 0.431 0.745 0.415 ];
f.NeurA.RealBarsHead = bar(monkeyTimes, nanmean(monkeyDatHead),'FaceColor','flat'); hold on; box off
f.NeurA.RealBarsHead.CData = zeros(size(f.NeurA.RealBarsHead.CData));
plot([0 0], [0 1.8],'-.k')

view([90 90])

f.NeurA.ax{5}.XAxis.Visible = 'off';
ylabel('Activity (Hz)')







             

% MODELLED DATA
% =========================================================================

% -------------------------------------------------------------------------
f.NeurA.ax{6} = axes('Position',[.6 .8 .2 .15 ]);
tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\NeurNET.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off




% -------------------------------------------------------------------------
f.NeurA.ax{7} = axes('Position',[.6 .45 .25 .25 ]);
CM = 4;

s = rS(cM).s;
w = rS(cM).w;

% $$$ 5 1 is limb centred. 4 2
cL = 4;
cN = 2;

s.plt.intrpFact = 1;%.2;
s.plt.plAct = 1;
% % % s.plt.lmbCol=2:14;
% % % s.plt.meanLimbCols = 1;

s.plt.lmbCol= 11;
s.plt.meanLimbCols = 0;

s.plt.bdyCol=2:14;
s.plt.pltType = 'Imagesc';
DisplActValsFun(s,w,-neurActForPlot_MeanBody(:,:,:,:,cL,cN));
CustomColourMap2
colorbar off
box off

xlabel('Stimulus x-position'); ylabel('Stimulus y-position');
title('Artificial limb-centred neuron')

% caxis([-.6 .8])
caxis([-.7 .6])

ylim([1.5 9.5])
% % % xlim([-6.5 6.5])
xlim([3 13])

GridOverImage(fS,f.NeurA.ax{7});


% -------------------------------------------------------------------------
f.NeurA.ax{8} = axes('Position',[.875 .45 .075 .25 ]);
allLC = 2:14;
clear tmpPlt
for iLC = 1:length(allLC)
    cLC = allLC(iLC);
    tmpPlt(:,iLC) = squeeze( nanmean(-neurActForPlot_MeanBody(w.bdy.row,cLC,:,cLC,cL,cN)   ,[2 4]  )  );
end
tmpPlt(9:end,:) = NaN;
tmpPlt = abs(tmpPlt - nanmin(tmpPlt(:)));
f.NeurA.ModelBarsLimb = bar([8:-1:1], nanmean(tmpPlt(1:8,:),2),'FaceColor','flat'); hold on; box off
f.NeurA.ModelBarsLimb.CData = zeros(size(f.NeurA.ModelBarsLimb.CData));
view([90 -90])

xlim([-.5 7.5])
f.NeurA.ax{8}.XAxis.Visible = 'off';
ylabel('Activity')







% -------------------------------------------------------------------------
f.NeurA.ax{9} = axes('Position',[.6 .1 .25 .25 ]);
CM = 4;

s = rS(cM).s;
w = rS(cM).w;

cL = 4;
cN = 2;

s.plt.intrpFact = 1;%.2;
s.plt.plAct = 1;
% % % s.plt.lmbCol=2:14;
% % % s.plt.meanLimbCols = 1;

s.plt.lmbCol= 5;
s.plt.meanLimbCols = 0;

s.plt.bdyCol=2:14;
s.plt.pltType = 'Imagesc';
DisplActValsFun(s,w,-neurActForPlot_MeanLimb(:,:,:,:,cL,cN));
CustomColourMap2
colorbar off
box off

xlabel('Stimulus x-position'); ylabel('Stimulus y-position');
title('Artificial body-centred neuron')

caxis([-.6 .4])

ylim([1.5 9.5])
% % % xlim([-6.5 6.5])
xlim([3 13])

GridOverImage(fS,f.NeurA.ax{9});


% -------------------------------------------------------------------------
f.NeurA.ax{10} = axes('Position',[.875 .1 .075 .25 ]);
allLC = 2:14;
clear tmpPlt
for iLC = 1:length(allLC)
    cLC = allLC(iLC);
    tmpPlt(:,iLC) = squeeze( nanmean(-neurActForPlot_MeanLimb(w.bdy.row,cLC,:,cLC,cL,cN)   ,[2 4]  )  );
end
tmpPlt = tmpPlt - nanmin(tmpPlt(:));
f.NeurA.ModelBarsBody = bar([8:-1:1], nanmean(tmpPlt(1:8,:),2),'FaceColor','flat'); hold on; box off
f.NeurA.ModelBarsBody.CData = zeros(size(f.NeurA.ModelBarsBody.CData));
view([90 -90])

xlim([-.5 7.5])
f.NeurA.ax{10}.XAxis.Visible = 'off';
ylabel('Activity')





% -------------------------------------------------------------------------
f.LPos.ax{11} =  axes('Position',[.6 .43 .25 .1 ]);
tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\MonkeyTopDown1Pose.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off




% -------------------------------------------------------------------------
f.LPos.ax{12} =  axes('Position',[.6 .08 .25 .1 ]);
tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\MonkeyTopDown1Pose.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off



% -------------------------------------------------------------------------
f.NeurA.ax{13} = axes('Position',[.1 .43 .25 .1 ]);
% f.NeurA.ax{13} = axes('Position',[.1 .43 .25 .15 ]);
% tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\\MonkeyTopDown1Pose_WithLimbCentredAct.jpg');
tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\\MonkeyTopDown1Pose.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off

% -------------------------------------------------------------------------
f.NeurA.ax{14} = axes('Position',[.1 .08 .25 .1 ]);
% f.NeurA.ax{14} = axes('Position',[.1 .08 .25 .15 ]);
% tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\\MonkeyTopDown1Pose_WithHeadCentredAct.jpg');
tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\\MonkeyTopDown1Pose.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off




% -------------------------------------------------------------------------
yRange = .2;
%bar(monkeyTimes, nanmean(monkeyDatLimb),'FaceColor','flat'); hold on; box off
% f.NeurA.DepthProp.CData = zeros(size(f.NeurA.RealBarsLimb.CData));

fS.nBins = 4;
% fS.pltType = 'shadeplot';
fS.pltType = 'bar';


% create depth groups so that I can make pie charts
clear splGroups
for iM  = 1:length(rS)
    nLay = length(rS(iM).s.lp.netS);
    grSize = round(nLay./fS.nBins);
    cInd =1;
    for iSpl = 1:fS.nBins-1
        splGroups{iM,iSpl} = cInd : iSpl .* grSize;
        cInd = iSpl.* grSize + 1;
    end
    splGroups{iM,iSpl + 1} = cInd : nLay;
end


for iSpl = 1:fS.nBins
    
    f.NeurA.ax{15 + iSpl-1} = axes('Position',...
        [.85 (.75 + (iSpl-1) .* (yRange/fS.nBins)) .1 (yRange/fS.nBins) ]);

    for iM = 1:length(rS)
             
        cL = splGroups{iM,iSpl};
        propNeurs(iM,iSpl) = nanmean(rMat.propCorrNeurPerLay(cL,iM),[2 1]);
    end
   
%     p = pie([1 - nanmean(propNeurs(:,iSpl),1); nanmean(propNeurs(:,iSpl),1)],{'',''})
    p = pie([1 - nanmean(propNeurs(:,iSpl),1); nanmean(propNeurs(:,iSpl),1)])
    colormap(coltowhite(100,[0 0 1]));
    set(p,'EdgeColor','k')
    
    if iSpl == 1
        l = legend('% proximity-dependent units')
        l.Position = l.Position - [.1 .02 0 0]
    end
end



% f.NeurA.DepthProp = BinPlot(rMat_All.relLayDep,rMat_All.propCorrNeurPerLay,fS,'LineWidth',2)
xlabel('binned layer depth')
ylabel('binned proportion of neurons that correlate with proximity')

view([90 90]);



colormap(f.NeurA.ax{2}, whitetocol(100,[0 0 1]));

colormap(f.NeurA.ax{4}, whitetocol(100,[0 .6 0]));

colormap(f.NeurA.ax{7}, whitetocol(100,[0 0 1]));

colormap(f.NeurA.ax{9}, whitetocol(100,[0 .6 0]));


%% Effect on reaction time

% f.RT.f = figure('Position',[20 20 600 400]);
f.RT.f = figure('Position',[20 20 1200 800]);

f.RT.ax{1} = axes('Position',[.1 .1 .3 .5]);

humanMax30 = 2.132;

humanDat = ([1.577 -0.304] ./ humanMax30) .* 30;
humanSD  = ([1.577-1.209 , - .304+0.765] ./ humanMax30) .* 30;
f.RT.RealBars = bar(humanDat,'FaceColor','flat'); hold on
f.RT.RealBars.CData(2,:) = [1 0 0];
colrs = {'k','k'}
for iM = 1:2
    plot([iM iM] ,humanDat(iM) + [ - humanSD(iM) humanSD(iM)],colrs{iM},'LineWidth',2);
end
xlim([.5 2.5])
ylim([-12 40])

ylabel('Far - Near Reaction Time difference (s)')

f.RT.ax{1}.XTickLabel = {'Baseline', 'r-TMS'}

box off



% -------------------------------------------------------------------------
f.RT.ax{2} = axes('Position',[.55 .1 .3 .5]);

distances = [2:6 ; 8:12];

allLC = 2:14;
clear bpAvOverCols
for iM = 1:3
for iD = 1:size(distances,1)
for iLC = 2:length(allLC)
    cLC = allLC(iLC);
    bpAvOverCols(iLC,iD,iM) = nanmean(squeeze(rtRS(iM).avBP(w.lmb.row,cLC,distances(iD,:),cLC)));
end
end
end

% don't include the pure learned model
bpAvOverCols(:,:,1) = [];


tmpD    = squeeze(bpAvOverCols(:,2,:) - bpAvOverCols(:,1,:)   );
tmpM    = squeeze(nanmean( bpAvOverCols(:,2,:) - bpAvOverCols(:,1,:)   ));
tmpSD	= squeeze(nanstd( bpAvOverCols(:,2,:) - bpAvOverCols(:,1,:)  ))  ./sqrt(size(bpAvOverCols,1));
f.RT.ModelBars = bar(tmpM','FaceColor','flat'); box off
f.RT.ModelBars.CData(2,:) = [1 0 0];
hold on
boxplot(tmpD)
plotRowAveragesWithSpread_alt(tmpD');
hold on

for iM = 1:2
    for iD = 1:2
        plot([iM iM],tmpM(iM) + [ - tmpSD(iM) tmpSD(iM)],'k','LineWidth',2);
    end
end
xlim([.5 2.5])
ylim([-.2 1])

ylabel('Far - Near Reaction Time difference (a.u.)')

f.RT.ax{2}.XTickLabel = {'Baseline', sprintf('simulated \\newline TMS')}


% -------------------------------------------------------------------------
f.RT.ax{3} = axes('Position',[.1 .65 .3 .3]);
tmpFig = imread('Results\ForFigures\ProximityPosition\BitsAndPieces\TMS_v2.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off
title('TMS abolishes spatial RT effect')

% -------------------------------------------------------------------------
% f.RT.ax{4} = axes('Position',[.55 .65 .3 .3]);
f.RT.ax{4} = axes('Position',[.55 .725 .3 .225]);
tmpFig = imread('Results\ForFigures\ProximityPosition\BitsAndPieces\TMS_NeurNet_NoStimulator.jpg');
alphaChan = ~all(tmpFig > 250, 3);
image(tmpFig,'AlphaData',alphaChan); axis off
title('simulated TMS abolishes spatial RT effect')

tmpD = squeeze(bpAvOverCols(:,2,:) - bpAvOverCols(:,1,:) );


[pRT hRT statsRT] = signrank(tmpD(:,1),zeros(size(tmpD(:,1))),'method','approximate')
effSizeNoLesion = statsRT.zval ./ sqrt(numel(tmpD(:,1)))

[pRT hRT statsRT] = signrank(tmpD(:,2),zeros(size(tmpD(:,1))),'method','approximate')
effSizeLesion = statsRT.zval ./ sqrt(numel(tmpD(:,2)))

[pRTcomp hRTcomp statsRTcomp] = ranksum(tmpD(:,1),tmpD(:,2),'method','approximate')
effSizeComp = statsRTcomp.zval ./ sqrt(numel(tmpD))

%% Supplementary figure: individual neural activities

% NEXT ISSUE: the square networks give to wide plots... SO, should make the
% number of columns fixed across all the networks --> figure out how to do

clear v;

for iM = 2
    
    sFP     = DefaultSettings(rS(iM).s);
    w       = rS(iM).w;
    net     = rS(iM).net;
    
    sFP.plt.axesVis     = 0;
    sFP.plt.grid        = 0;
    sFP.plt.plAct       = 1;
    sFP.plt.startLayer  = 1;
    sFP.plt.colLims     = [1.5 14.5];
    sFP.plt.rowLims     = [1.5 14.5];
    sFP.plt.showLimb    = 1;
    
    sFP.plt.meanLimbCols = 0;
    
    % Ensure the same widths in the plot
    if iM == 2
        sFP.plt.nSubCols = 12;
    elseif iM == 5
        sFP.plt.nSubCols = 13;
    else
        sFP.plt.nSubCols = [];
    end
    
    if iM < 4
        sFP.plt.stopLayer   = 4;
        
% % %         % Single position plot
% % %         sFP.plt.vidFl       = 0;
% % %         figN                = ['MultNeur_M' num2str(iM)];
% % %         f.(figN).f          = figure('Position',[20 20 1800 600]);
% % %         sFP.plt.lmbCol      = 8;
% % %         NetAnalysis(sFP,w,permute(neurActForPlot{iM},[3 4 1 2 5 6]));
        
% % %             % Multi position video
% % %             sFP.plt.vidFl         = 1;
% % %             sFP.plt.vidFileName   = ['D:\Old_D\DPPS\DefenseAgent\Results\Videos\LimbRight_M' num2str(iM) '.avi'];
% % %             sFP.plt.lmbCol      = 2:14;
% % %             figure('Position',[20 20 1800 600]);
% % %             NetAnalysis(sFP,w,permute(neurActForPlot{iM},[3 4 1 2 5 6]))

            % Multi position video - forward and back
            sFP.plt.vidFl         = 1;
            sFP.plt.vidFileName   = ['D:\Old_D\DPPS\DefenseAgent\Results\Videos\LimbRight_M' num2str(iM) '_ForwardBack.avi'];
            sFP.plt.lmbCol      = [2:14 13:-1:3];
            sFP.plt.vidFR       = 2.5;
            figure('Position',[20 20 1800 600]);
            NetAnalysis(sFP,w,permute(neurActForPlot{iM},[3 4 1 2 5 6]))
    else
        
        sFP.plt.stopLayer   = 5;
        
        % Single position plot
        sFP.plt.vidFl       = 0;
        figN                = ['MultNeur_M' num2str(iM)];
        f.(figN).f          = figure('Position',[20 20 1800 900]);
        sFP.plt.lmbCol      = 5;
        sFP.plt.showBodyCol = 11;
        NetAnalysis(sFP,w,neurActForPlot{iM}(:,:,:,:,:,:,sFP.plt.showBodyCol));
        
        
% % %             % Multiple position video
% % %             clear v
% % %             sFP.plt.vidFl         = 0;
% % %             sFP.plt.vidFileName   = ['D:\Old_D\DPPS\DefenseAgent\Results\Videos\BodyLeftLimbRight_M' num2str(iM) '.avi'];
% % %             v                     = VideoWriter(sFP.plt.vidFileName);
% % %             v.FrameRate           = sFP.plt.vidFR;
% % %             open(v)
% % %             
% % %             vidFileName2          = ['D:\Old_D\DPPS\DefenseAgent\Results\Videos\BodyStayLimbRight_M' num2str(iM) '.avi'];
% % %             v2                    = VideoWriter(vidFileName2);
% % %             v2.FrameRate          = sFP.plt.vidFR;
% % %             open(v2)
% % %             
% % %             vidFileName3          = ['D:\Old_D\DPPS\DefenseAgent\Results\Videos\BodyLeftLimbStay_M' num2str(iM) '.avi'];
% % %             v3                    = VideoWriter(vidFileName3);
% % %             v3.FrameRate          = sFP.plt.vidFR;
% % %             open(v3)
% % %             
% % %             vidFileName4          = ['D:\Old_D\DPPS\DefenseAgent\Results\Videos\BodyAvLimbRight_M' num2str(iM) '.avi'];
% % %             v4                    = VideoWriter(vidFileName4);
% % %             v4.FrameRate          = sFP.plt.vidFR;
% % %             open(v4)
% % %             
% % %             allC = 2:14;
% % %             figure('Position',[20 20 1800 900]);
% % %             for iC = 1:length(allC)
% % %                 
% % %                 % Move opposite directions
% % %                 sFP.plt.lmbCol      = allC(iC);
% % %                 sFP.plt.showBodyCol = 15 - sFP.plt.lmbCol;
% % %                 NetAnalysis(sFP,w,neurActForPlot{cM}(:,:,:,:,:,:,sFP.plt.showBodyCol));
% % %                 frame = getframe(gcf);
% % %                 writeVideo(v,frame);
% % %                 hold off;
% % %                 
% % %                 % Move only limb
% % %                 sFP.plt.showBodyCol = 11;
% % %                 NetAnalysis(sFP,w,neurActForPlot{cM}(:,:,:,:,:,:,sFP.plt.showBodyCol));
% % %                 frame = getframe(gcf);
% % %                 writeVideo(v2,frame);
% % %                 hold off;
% % %                 
% % %                 % Move only body
% % %                 sFP.plt.lmbCol = 3;
% % %                 sFP.plt.showBodyCol = 15 - allC(iC);
% % %                 NetAnalysis(sFP,w,neurActForPlot{cM}(:,:,:,:,:,:,sFP.plt.showBodyCol));
% % %                 frame = getframe(gcf);
% % %                 writeVideo(v3,frame);
% % %                 hold off;
% % %                 
% % %                 % Move only limb, av over body
% % %                 sFP.plt.lmbCol = allC(iC);
% % %                 sFP.plt.showBodyCol = [];
% % %                 NetAnalysis(sFP,w,nanmean(neurActForPlot{cM}(:,:,:,:,:,:,2:end-1),7));
% % %                 frame = getframe(gcf);
% % %                 writeVideo(v4,frame);
% % %                 hold off;
% % %                 
% % %             end
% % %             close(v)
% % %             close(v2)
% % %             close(v3)
% % %             close(v4)
        
    end
    
end

%% Save figures

allFields = fields(f);
for iF = 1:length(allFields)
    cF = allFields{iF};
    
%     cFig = figure(f.(cF).f)
%     print(cFig,['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\FromMatlab\' cF ,'-dpfd'])
    set(f.(cF).f, 'Renderer', 'painters'); % default, opengl
    saveas(f.(cF).f,['Results\ForFigures\ProximityPosition\BitsAndPieces\FromMatlab\' cF '_V2.eps'] , 'epsc')
    saveas(f.(cF).f,['Results\ForFigures\ProximityPosition\BitsAndPieces\FromMatlab\' cF '_V2.pdf'] , 'pdf')
%     saveas(f.(cF).f,['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\FromMatlab\' cF '.svg'] , 'svg')
end

%%

% % % addpath('D:\Old_D\Programs\Matlab\Toolboxes\export_fig')
% % % % $$$$ USE THIS!
% % % 
% % % 
% % % allFields = fields(f);
% % % for iF = 1:length(allFields)
% % %     cF = allFields{iF};
% % %     
% % %     cFig = figure(f.(cF).f)
% % % 
% % %     cName = ['D:\Old_D\DPPS\DefenseAgent\Results\ForFigures\ProximityPosition\BitsAndPieces\FromMatlab\' cF '.pdf'];
% % %     
% % %     switch cF
% % %         case 'HBR'
% % %             export_fig HBR.pdf
% % %     end
% % % end





% % % 
% % % %% Combine figures
% % % 
% % % 
% % % 
% % % % Prepare subplots
% % % figure
% % % h(1)=subplot(1,2,1);
% % % h(2)=subplot(1,2,2);
% % % % Paste figures on the subplots
% % % copyobj(allchild(get(f.NeurAx.f,'CurrentAxes')),h(1));
% % % copyobj(allchild(get(k,'CurrentAxes')),h(2));


% % % % % 
% % % % % 
% % % % % %%
% % % % % s.plt.plAct = 2;
% % % % % s.plt.lmbCol = 2:14;
% % % % % s.plt.meanLimbCols = 1;
% % % % % s.plt.bdyCol = 2;
% % % % % s.plt.stimCol = 2:14;
% % % % % s.plt.stimRow=[1:size(w.world2D,1)];
% % % % % 
% % % % % s.plt.intrpFact = 1;
% % % % % 
% % % % % Q_MeanBody2 = Q_MeanBody;
% % % % % Q_MeanLimb2 = Q_MeanLimb;
% % % % % for iAct = 1:size(Q_MeanBody,5)
% % % % %     Q_MeanBody2(:,:,:,:,iAct) = nanmean(Q_MeanBody,5);
% % % % %     Q_MeanLimb2(:,:,:,:,iAct) = nanmean(Q_MeanLimb,5);
% % % % % end
% % % % % 
% % % % % hold off
% % % % % DisplActValsFun(s,w,Q_MeanLimb2);
% % % % % % DisplActValsFun(s,w,Q_MeanBody2);
% % % % % ylim([2 10.5])
% % % % % xlim([-7.5 7.5])
% % % % % hold off
% % % % % 
% % % % % % $$$ DO I set collag and rowlag to a specific value?
% % % % % % $$$ DO I PLOT NEURACT FOR SPECIFIC LIMBS?
% % % % % 
% % % % % %%
% % % % % 
% % % % % cM = 4; %length(rS)-1;
% % % % % 
% % % % % % [r c] = find(rS(cM).bdy.pDistN.max<0.05);
% % % % % [r c] = find(rS(cM).lmb.pDistN.max<0.05);
% % % % % 
% % % % % s.plt.meanLimbCols = 1;
% % % % % s.plt.lmbCol = [2:14];
% % % % % s.plt.stimRow = [2:10];
% % % % % 
% % % % % iN = 6; % 6 IS GOOD
% % % % % 
% % % % % s.plt.plAct = 1;
% % % % % s.plt.pltType='Imagesc';
% % % % % 
% % % % % subplot(1,2,1)
% % % % % DisplActValsFun(s,w,neurActForPlot_MeanBody(:,:,:,:,r(iN),c(iN)))
% % % % % ylim([2 10.5])
% % % % % xlim([-7.5 7.5])
% % % % % title('MeanBody')
% % % % % 
% % % % % rS(cM).lmb.pDistN.max(r(iN),c(iN))
% % % % % 
% % % % % 
% % % % % subplot(1,2,2)
% % % % % DisplActValsFun(s,w,neurActForPlot_MeanLimb(:,:,:,:,r(iN),c(iN)))
% % % % % ylim([2 10.5])
% % % % % xlim([-7.5 7.5])
% % % % % title('MeanLimb')
% % % % % 
% % % % % rS(cM).bdy.pDistN.max(r(iN),c(iN))
% % % % % 
% % % % % % $$$ LEt's have a look how good the single neuron plots are for this
% % % % % % network
% % % % % 
% % % % % 
% % % % % s = rS(cM).s;
% % % % % w = rS(cM).w;
% % % % % 
% % % % % 
% % % % % s.plt.intrpFact = 0.1;
% % % % % s.plt.stimRow=[1:size(w.world2D,1)];
% % % % % s.plt.stimCol=[2:size(w.world2D,1)];
% % % % % 
% % % % % 
% % % % % figure('Position',[0 0 800 400])
% % % % % s.plt.lmbCol=5;
% % % % % s.plt.bdyCol=2:14;
% % % % % s.plt.pltType='Imagesc';
% % % % % s.plt.stopLayer = 5;
% % % % % NetAnalysis(s,w,permute(neurActForPlot_MeanLimb,[3 4 1 2 5 6]));
% % % % % sgtitle('meanLimb')
% % % % % CustomColourMap2
% % % % % 
% % % % % 
% % % % % figure('Position',[0 0 800 400])
% % % % % s.plt.pltType='Imagesc';
% % % % % s.plt.stopLayer = 5;
% % % % % s.plt.lmbCol = 8;
% % % % % NetAnalysis(s,w,permute(neurActForPlot_MeanBody,[3 4 1 2 5 6]));
% % % % % sgtitle('meanBody')
% % % % % CustomColourMap2
% % % % % 
% % % % % 
% % % % % % BODY PLOT
% % % % % figure
% % % % % cL = 5;
% % % % % cN = 2;
% % % % % s.plt.plAct = 1;
% % % % % s.plt.lmbCol=5;
% % % % % s.plt.bdyCol=2:14;
% % % % % s.plt.pltType = 'Imagesc';
% % % % % DisplActValsFun(s,w,neurActForPlot_MeanLimb(:,:,:,:,cL,cN));
% % % % % figure,
% % % % % % s.plt.pltType = 'Binned';
% % % % % % s.plt.fS.pltType = 'bar'; %'shadeplot';
% % % % % % s.plt.bdyCol=5;
% % % % % % DisplActValsFun(s,w,neurActForPlot_MeanLimb(:,:,:,:,cL,cN));
% % % % % 
% % % % % allLC = 2:14;
% % % % % clear tmpPlt
% % % % % for iLC = 1:length(allLC)
% % % % %     cLC = allLC(iLC);
% % % % %     tmpPlt(:,iLC) = squeeze( nanmean(neurActForPlot_MeanLimb(w.bdy.row,cLC,:,cLC,cL,cN)   ,[2 4]  )  );
% % % % % end
% % % % % tmpPlt = tmpPlt - nanmin(tmpPlt(:));
% % % % % bar([8:-1:1], nanmean(tmpPlt(1:8,:),2));
% % % % % view([90 90])
% % % % % 
% % % % % 
% % % % % 
% % % % % %%
% % % % % % % % % -------------------------------------------------------------------------
% % % % % % % % f.RT.ax{1} = axes('Position',[.1 .1 .3 .6]);
% % % % % % % % 
% % % % % % % % humanDat = [600 583 ; 603 607]';
% % % % % % % % humanSD  = [ 34  32 ; 33 31]';
% % % % % % % % bar(humanDat'); hold on
% % % % % % % % for iM = 1:2
% % % % % % % %     for iD = 1:2
% % % % % % % %         if iD == 1
% % % % % % % %             barOffset = 1 - .8555;
% % % % % % % %         else
% % % % % % % %             barOffset = -1 + .8555;
% % % % % % % %         end
% % % % % % % %         plot([iM iM] - barOffset,humanDat(iD,iM) + [ - humanSD(iD,iM) humanSD(iD,iM)],'k','LineWidth',2);
% % % % % % % %         plot([iM] - barOffset,humanDat(iD,iM) ,'xk','LineWidth',2);
% % % % % % % %     end
% % % % % % % % end
% % % % % % % % 
% % % % % % % % ylim([540 660])
% % % % % % -------------------------------------------------------------------------
% % % % % 
% % % % % %% Q-value correlation with distance
% % % % % 
% % % % % handRow = 12;
% % % % % kAction = 1;
% % % % % 
% % % % % Qtable = rS(1).Qtable;
% % % % % 
% % % % % 
% % % % % cM = 1;
% % % % % 
% % % % % s = rS(cM).s;
% % % % % w = rS(cM).w;
% % % % % 
% % % % % iC=0; clear cD rD aD
% % % % % for hC=2:15
% % % % %     iC=iC+1;
% % % % %     [cD(:,:,iC),rD(:,:,iC),aD(:,:,iC)] = CalcDist1D(s,w,hC,handRow);
% % % % %     qTemp(:,:,iC)=squeeze(Qtable(handRow,hC,:,:,kAction));
% % % % %     % % %     plot(aD(:,:,iC),qTemp(:,:,iC),'x');hold on
% % % % % end
% % % % % % % % title('Action value falls of as a function of distance')
% % % % % % % % xlabel('Distance to Limb'); ylabel('Value of action')
% % % % % % % % allAD=unique(aD);
% % % % % % % % for iAD=1:length(allAD)
% % % % % % % %     meanQ(iAD)=nanmean(qTemp(aD(:)==allAD(iAD)));
% % % % % % % %     sdQ(iAD)=nanstd(qTemp(aD(:)==allAD(iAD)));
% % % % % % % % end
% % % % % % % % plot(allAD(2:end),meanQ(2:end),'k','LineWidth',2); hold on
% % % % % % % % plot(allAD(2:end),meanQ(2:end)+sdQ(2:end)/2,'r');
% % % % % % % % plot(allAD(2:end),meanQ(2:end)-sdQ(2:end)/2,'r')
% % % % % % % % hold off
% % % % % 
% % % % % % $$$ Maybe I should change this into a plot showing the response near the
% % % % % % hand? $$$ Either way, needs to be made prettier
% % % % % figure,
% % % % % tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Documentation\Figures\FromOtherPapers\MonkeyForwardOnly.png');
% % % % % image(tmpFig); axis off
% % % % % title('"Peripersonal space" neuron firing rate falls of as a function of distance')
% % % % % % set(gcf,
% % % % % 
% % % % % % figure('Position',[0 0 600 500]),
% % % % % % tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Documentation\Figures\FromOtherPapers\BlinkExample_Bufacchi2018.png');
% % % % % % image(tmpFig); axis off
% % % % % % title('Behavioural "Peripersonal space" responses also fall of as a function of distance')
% % % % % 
% % % % % figure,
% % % % % nBins=10; clear binQ
% % % % % [Y,E] = discretize(aD(:),nBins);
% % % % % for iBin=1:nBins
% % % % %     binQ(iBin)=nanmean(qTemp(Y==iBin));
% % % % % end
% % % % % bar((E(1:end-1)+E(2:end))./2,binQ); xlabel('Distance from "limb"'); ylabel('Average action value');
% % % % % title('Action value also falls of as a function of distance to a limb: absolute distance')
% % % % % box off
% % % % % 
% % % % % 
% % % % % 
% % % % % figure,
% % % % % nBins=10; clear binQ
% % % % % [Y,E] = discretize(cD(:),nBins);
% % % % % for iBin=1:nBins
% % % % %     binQ(iBin)=nanmean(qTemp(Y==iBin));
% % % % % end
% % % % % bar((E(1:end-1)+E(2:end))./2,binQ); xlabel('Distance from "limb"'); ylabel('Average action value');
% % % % % title('Action value also falls of as a function of distance to a limb: left-right distance')
% % % % % box off
% % % % % 
% % % % % 
% % % % % 
% % % % % figure,
% % % % % nBins=10; clear binQ
% % % % % [Y,E] = discretize(rD(:),nBins);
% % % % % for iBin=1:nBins
% % % % %     binQ(iBin)=nanmean(qTemp(Y==iBin));
% % % % % end
% % % % % bar((E(1:end-1)+E(2:end))./2,binQ); xlabel('Distance from "limb"'); ylabel('Average action value');
% % % % % title('Action value also falls of as a function of distance to a limb: forward-back distance')
% % % % % box off
% % % % % 
% % % % % %% Response moving with limb
% % % % % 
% % % % % handRow=12 ;
% % % % % kAction=2;
% % % % % 
% % % % % cM = 1;
% % % % % 
% % % % % Qtable = rS(cM).Qtable;
% % % % % 
% % % % % s = rS(cM).s;
% % % % % w = rS(cM).w;
% % % % % 
% % % % % figure('Position',0.6.*[0 0 1000 1000])
% % % % % subplot(10,1,[1:4])
% % % % % tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Documentation\Figures\FromOtherPapers\MonkeyMoveHand_GrazianoGross1993_2.png');
% % % % % image(tmpFig); axis off
% % % % % title('Top-down sketch of neuron response field around hand')
% % % % % subplot(10,1,6:10)
% % % % % tmpFig = imread('D:\Old_D\DPPS\DefenseAgent\Documentation\Figures\FromOtherPapers\MonkeyMoveHand_GrazianoGross1997.png');
% % % % % image(tmpFig); axis off
% % % % % title('Neuron response while visual stimulus moves towards body')
% % % % % sgtitle('"Peripersonal space" neuron visual response field moves with limb position')
% % % % % 
% % % % % figure,
% % % % % subplot(2,2,1)
% % % % % handCol=8;
% % % % % % imagesc(squeeze(Qtable(handRow,handCol,:,:,kAction))); colorbar; axis square;
% % % % % ImagescInterp(squeeze(Qtable(handRow,handCol,:,:,kAction)),1); colorbar; axis square;
% % % % % A=gca;
% % % % % box off
% % % % % xlabel('Stimulus x-position');
% % % % % ylabel('Stimulus y-position'); title('"Limb" on the left')
% % % % % 
% % % % % subplot(2,2,3)
% % % % % handCol=8; colorbar
% % % % % h=plot(squeeze(nanmean(Qtable(handRow,handCol,:,:,kAction),3)),'-o'); axis square; %colorbar
% % % % % A=gca; set(A,'Position',[0.1300    0.1100    0.2215    0.3412]);
% % % % % xlabel('Stimulus x-position'); ylabel('Average value')
% % % % % 
% % % % % subplot(2,2,2)
% % % % % handCol=11;
% % % % % imagesc(squeeze(Qtable(handRow,handCol,:,:,kAction))); colorbar; axis square;
% % % % % box off
% % % % % xlabel('Stimulus x-position');
% % % % % ylabel('Stimulus y-position'); title('"Limb" on the right')
% % % % % 
% % % % % subplot(2,2,4)
% % % % % handCol=11; colorbar
% % % % % h=plot(squeeze(nanmean(Qtable(handRow,handCol,:,:,kAction),3)),'-o'); axis square; %colorbar
% % % % % A=gca; set(A,'Position',[0.5703    0.1100    0.2215    0.3412]);
% % % % % xlabel('Stimulus x-position'); ylabel('Average value')
% % % % % 
% % % % % colormap(whitetocol(100,[0 0 0.7]));
% % % % % 
% % % % % sgtitle('Value fields also move with limb position')
% % % % % 
% % % % % figure, 
% % % % % handCol=8; colorbar
% % % % % h=plot(squeeze(nanmean(Qtable(handRow,handCol,:,:,kAction),3)),'-o','LineWidth',2); axis square; %colorbar
% % % % % hold on 
% % % % % handCol=11;
% % % % % h=plot(squeeze(nanmean(Qtable(handRow,handCol,:,:,kAction),3)),'-o','LineWidth',2); axis square; %colorbar
% % % % % legend('limb left','limb right')
% % % % % xlabel('Stimulus x-position'); ylabel('Average value')
% % % % % title('Value fields also move with limb position')
% % % % % 
% % % % % 
% % % % % %% 2D plots
% % % % % cM = 1;
% % % % % 
% % % % % s   = rS(cM).s;
% % % % % w   = rS(cM).w;
% % % % % net = rS(cM).net;
% % % % % 
% % % % % s.plt.intrpFact = 0.1
% % % % % s.plt.lmbCol=5;
% % % % % s.plt.bdyCol=10;
% % % % % s.plt.rowLims=[1.5 s.wrld.size(1)-0.5];
% % % % % [Q,allNeurAct] = CalcNetOutput(s,w,net);
% % % % % 
% % % % % s.plt.intrpFact=.1;
% % % % % s.plt.OneActFl=1;
% % % % % 
% % % % % figure('Position',[0 0 600 400])
% % % % % 
% % % % % s.plt.meanLimbCols = 1;
% % % % % 
% % % % % s.plt.lmbCol=2:12;
% % % % % s.plt.bdyCol=5;
% % % % % s.plt.stimRow=[1:size(w.world2D,1)];
% % % % % s.plt.stimCol=[2:size(w.world2D,1)];
% % % % % 
% % % % % DisplActValsFun(s,w,Q);
% % % % % xlabel('Stimulus x-position'); ylabel('Stimulus y-position');
% % % % % title('Network value output')
% % % % % 
% % % % % colormap(whitetocol(100,[0 0 0.7])); grid off; colorbar
% % % % % 
% % % % % 
% % % % % figure,
% % % % % clear rowAct rowDist
% % % % % for iR = 1:size(allNeurAct,1)
% % % % %     for iC = 1:size(allNeurAct,1)
% % % % %         tmpAct = abs(allNeurAct(iR,:,1,:,:,:) - nanmedian(allNeurAct(iR,:,1,:,:,:),2) );
% % % % %         rowAct(iR,:) = tmpAct(:);
% % % % %         rowDist(iR,:) = repmat(iR,[1 size(rowAct,2)]);
% % % % %     end
% % % % % end
% % % % % % % bar([11:-1:1],nanmedian(rowAct(1:11,:),2))
% % % % % bar([10:-1:1],nanmedian(rowAct(1:10,:),2))
% % % % % xlabel('Row distance from hand')
% % % % % ylabel('Absolute neural activation')
% % % % % title('Median activation across all neurons in network')
% % % % % 
% % % % % 
% % % % % s.plt.lmbCol=8;
% % % % % s.plt.bdyCol=8;
% % % % % 
% % % % % figure('Position',[0 0 800 400])
% % % % % s.plt.pltType='Imagesc';
% % % % % NetAnalysis(s,w,allNeurAct);
% % % % % s.plt.stopLayer = 5
% % % % % NetAnalysis(s,w,permute(neurActForPlot_MeanLimb,[3 4 1 2 5 6]));
% % % % % 
% % % % % sgtitle('Later layers are more consistently affected by stimulus distance and limb position: response fields')
% % % % % 
% % % % % 
% % % % % % colormap(whitetocol(100,[0 0 0.7])); grid off; %colorbar
% % % % % CustomColourMap2
% % % % % 

% %% Neurons' correlation with distance
% 
% % fS.nDepthBins = 5;
% % fS.npropBins = 3;
% 
% fS.nDepthBins = 3;
% fS.npropBins = 3;
% figure,
% 
% [freqs XY] =hist3([rMat_All.propCorrNeurPerLay(:), rMat_All.relLayDep(:)],[fS.npropBins fS.nDepthBins]);
% imagesc(XY{2},XY{1},freqs); axis xy
% % % ImagescInterp(freqs,0.3,XY{2},XY{1}); axis xy
% 
% % % imagesc([0:0.2:1],[0:0.5:1],hist3([rMat_All.propCorrNeurPerLay(:), rMat_All.relLayDep(:)],[fS.npropBins fS.nDepthBins])); axis xy
% colormap(whitetocol(100,[0 0 0.7]));
% hold on
% plot([rMat_All.relLayDep],[rMat_All.propCorrNeurPerLay],'-ok','LineWidth',2)
% 
% xlabel('layer depth')
% ylabel('proportion of neurons that correlate with proximity')
% 
% fS.nBins = 4;
% % fS.pltType = 'shadeplot';
% fS.pltType = 'bar';
% figure,
% BinPlot(rMat_All.relLayDep,rMat_All.propCorrNeurPerLay,fS,'LineWidth',2)
% xlabel('binned layer depth')
% ylabel('binned proportion of neurons that correlate with proximity')

% % % % % 
% % % % % %% 'Psychophysical' responses
% % % % % 
% % % % % % figure('Position',[0 0 1200 500])
% % % % % 
% % % % % for iM = 1:2
% % % % % 
% % % % % rtRS(iM).s.fl.extraInput
% % % % % 
% % % % % sFP = DefaultSettings(sFP);
% % % % % 
% % % % % % sFP.plt.pltType = 'Imagesc';
% % % % % sFP.plt.pltType = 'Binned';
% % % % % sFP.plt.ON = 1;
% % % % % sFP.plt.meanLimbCols = 1;
% % % % % 
% % % % % sFP.plt.bdyCol = 2:12;
% % % % % 
% % % % % 
% % % % % % sFP.plt.stimRow = [3:size(w.world2D,1)-1];
% % % % % sFP.plt.stimRow = [1:size(w.world2D,1)-4];
% % % % % 
% % % % % sFP.plt.stimCol = [2:size(w.world2D,2)-1];
% % % % % sFP.plt.lmbCol = [2:size(w.world2D,2)-1];
% % % % % % sFP.plt.stimCol = 8;
% % % % % % sFP.plt.lmbCol = 8;
% % % % % % sFP.plt.stimCol = [6:10];
% % % % % % sFP.plt.lmbCol = 8;
% % % % % 
% % % % % 
% % % % % % sFP.plt.fS.pltType = 'bar';
% % % % % sFP.plt.fS.pltType = 'shadeplot';
% % % % % sFP.plt.fS.nBins = 10;
% % % % % 
% % % % % sFP.plt.varargs = {'LineWidth',2};
% % % % % 
% % % % % % sFP.plt.distanceType = 'Absolute';
% % % % % % sFP.plt.distanceType = 'Column';
% % % % % % sFP.plt.distanceType = 'AbsColumn';
% % % % % sFP.plt.distanceType = 'AbsRow';
% % % % % % sFP.plt.distanceType = 'Row';
% % % % % 
% % % % % 
% % % % % sFP.plt.plAct = 1;
% % % % % % subplot(1,2,iM);
% % % % % figure,
% % % % % DisplActValsFun(sFP,w, rtRS(iM).avBP);
% % % % % ylim([0 1]);
% % % % % xlabel('Distance to limb');
% % % % % ylabel('Response strength to "tactile" stimulus')
% % % % % 
% % % % % if iM == 1
% % % % %     title(sprintf('No PPS subnetwork: dist correlation rho = %.3f, P = %.3f',...
% % % % %         rtRS(iM).rDistBP.rDabs, rtRS(iM).pDistBP.rDabs))
% % % % % else
% % % % %     title(sprintf('With PPS subnetwork: dist correlation rho = %.3f, P = %.3f',...
% % % % %         rtRS(iM).rDistBP.rDabs, rtRS(iM).pDistBP.rDabs))
% % % % % end
% % % % % 
% % % % % 
% % % % % end
% % % % % 
% % % % % % sgtitle('Drive to "press button" ')
% % % % % 




























