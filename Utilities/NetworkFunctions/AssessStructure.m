function [rS,netAR,sumBinNPFdist,sumBinNPFpmdist,binNPF,binNPFpm] = AssessStructure(rS)
% [rS,netAR,allNetAR,sumBinNPFdist,sumBinNPFpmdist,binNPF,binNPFpm] = AssessStructure(rS)
%   Calculates various metrics of network structure




s=DefaultSettings(rS.s);
w=rS.w;
net=rS.net;
nPm = rS.s.plt.nPm;

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



figure,
hold on;
% First sort nodes by distance to all other nodes (i.e. make a distance
% matrix)
clear nodeDists closestDists closestNodes
scatter(p.XData,p.YData,1000,[zeros(size(A_B_rat)) 1-A_B_rat A_B_rat], '.')
for iNode=1:size(G.Nodes,1)

    if strcmp(s.nta.compScale,'Graded')
        %         scatter(p.XData(iNode),p.YData(iNode),1000,[0 1-A_B_rat(iNode) A_B_rat(iNode)], '.')
        %         scatter(p.XData(iNode),p.YData(iNode),1000,[0 .5-4.*A_B_rat(iNode) .5+4.*A_B_rat(iNode)], '.')
    end

    nodeDists(:,iNode)    = sqrt(sum(([p.XData(iNode),p.YData(iNode)] - [p.XData(:),p.YData(:)])'.^2));
    [closestDists(:,iNode) closestNodes(:,iNode)] = sort(nodeDists(:,iNode));
end


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

close all


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
    %for iL = 1:length(rS.s.lp.netS)
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


% % % % ===========================================================
% % % % See if the preferences change as a function of depth
% % % for iR = 1:size(binNPF,1)
% % %     for iC = 1:size(binNPF,2)
% % %         pPermSpace(iR,iC) = sum(abs(binNPFpm(iR,iC,:)) > abs(binNPF(iR,iC)) ) ./ size(binNPFpm,3);
% % %     end
% % %     tmpPm = squeeze(nanmean(abs(binNPFpm(iR,:,:)),2))  ;
% % %     tmpOg = nanmean(abs(binNPF(iR,:)),2);
% % %     if isnan(tmpOg)
% % %         pPermLine(iR) = NaN;
% % %     else
% % %         pPermLine(iR) = sum( tmpPm > tmpOg ) ./ size(binNPFpm,3);
% % %     end
% % % end
% % % 
% % % 
% % % for iC = 1:size(binNPF,2)
% % %     tmpPm = squeeze(nanmean(abs(binNPFpm(:,iC,:)),1))  ;
% % %     tmpOg = nanmean(abs(binNPF(:,iC)),1);
% % %     if isnan(tmpOg)
% % %         pPermLine2(iC) = NaN;
% % %     else
% % %         pPermLine2(iC) = sum( tmpPm > tmpOg ) ./ size(binNPFpm,3);
% % %     end
% % % end


end