function [rr pp covVal] = NetAnalysis(s,w,allNeurAct)
% Plot neural activation as a function of distance from the hand

if isempty(s.plt.lmbRow)
    handRow=size(w.world2D,1)-2;
else
    handRow=s.plt.lmbRow;
end

if s.plt.vidFl==1
    v = VideoWriter(s.plt.vidFileName);
    v.FrameRate=s.plt.vidFR;
    open(v)
end

if s.plt.sequentialLimbCols == 1
    % current Hand Column
    for iiHC = 1:length(s.plt.lmbCol)

        cHC = s.plt.lmbCol(iiHC);

        [r p covVal] = MultNeurPlotFun(s,w,allNeurAct,handRow,cHC);


        rr(:,:,iiHC) = r;
        pp(:,:,iiHC) = p;


        if s.plt.vidFl == 1
            frame=getframe(gcf)
            writeVideo(v,frame)
        end
    end

    if s.plt.vidFl == 1
        close(v)
    end

else


    [rr pp covVal n] = MultNeurPlotFun(s,w,allNeurAct(s.plt.stimRow,s.plt.stimCol,   :,:,:,:),handRow,s.plt.lmbCol);

end

% replace meaningless zeros with NaNs
for iL=1:length(s.lp.netS)
    if size(rr,1) < iL
        rr(iL,:) = NaN;
        pp(iL,:) = NaN;
        rr(iL,:) = NaN;
    end
    rr(iL,s.lp.netS(iL)+1:end)=NaN;
    pp(iL,s.lp.netS(iL)+1:end)=NaN;
    covVal(iL,s.lp.netS(iL)+1:end)=NaN;
end



    function [ r2 p2 covVal nA cD rD aD ] = MultNeurPlotFun(s,w,allNeurAct,handRow,cHC)

        [cD,rD,aD] = CalcDist1D(s,w,cHC,handRow);

        switch s.plt.distanceType
            case 'Absolute'
                distanceVals=aD;
            case 'Row'
                distanceVals=rD;
            case 'Column'
                distanceVals=cD;
            case 'AbsRow'
                distanceVals=abs(rD);
            case 'AbsColumn'
                distanceVals=abs(cD);
        end
        distanceVals=distanceVals(s.plt.stimRow,s.plt.stimCol,:);

        for cL=s.plt.startLayer:s.plt.stopLayer

            % plLoc = (cL-s.plt.startLayer)*max(s.lp.netS(s.plt.startLayer:s.plt.stopLayer)) + 1;
            plLoc = (cL-s.plt.startLayer)*max(s.lp.netS(s.plt.startLayer:s.plt.stopLayer)) + 1 + ...
                floor(( max(s.lp.netS) - s.lp.netS(cL) )./2);

            % Adjust current subplot position if making the subplots a
            % specific width
            if ~isempty(s.plt.nSubCols)
                plLoc =  plLoc + floor(  ( s.plt.nSubCols - max(s.lp.netS(s.plt.startLayer:s.plt.stopLayer))  )./2     .* ...
                    ( (cL-1).*2 + 1) );
            end

            for cN=1:s.lp.netS(cL)

                if s.plt.ON == 1
                    if isempty(s.plt.nSubCols)
                        subplot(s.plt.stopLayer-(s.plt.startLayer-1),...
                            max(s.lp.netS(s.plt.startLayer:s.plt.stopLayer)),plLoc);
                    else
                        subplot(s.plt.stopLayer-(s.plt.startLayer-1),...
                            s.plt.nSubCols,plLoc);
                    end
                    plLoc=plLoc+1;
                end

                % Neural Activation
                if strcmp(s.plt.pltType,'Binned')
                    nA=squeeze(allNeurAct(:,:,handRow,cHC,cL,cN)); % $$$ Neuron 1 is giving no activity for some reason

                    [Y,E] = discretize(distanceVals(:),s.plt.nBins);
                    for iBin=1:s.plt.nBins
                        binA(iBin,1:sum(Y==iBin))=nA(Y==iBin);
                        binAMean(iBin)=nanmean(nA(Y==iBin));
                        %                         binAsd(iBin) =  nanstd(nA(Y==iBin));
                        binPct(iBin,1) = prctile(nA(Y==iBin),10);
                        binPct(iBin,2) = prctile(nA(Y==iBin),90);
                    end
                    binA(binA==0)=NaN;
                    xVals=(E(1:end-1)+E(2:end))./2;

                    if s.plt.separateLimbCols==1 && s.plt.ON == 1
                        plot(repmat(xVals(:),[size(binA,2) 1]),binA(:),'.','LineWidth',0.5); hold off
                    end

                    x = [xVals  flip(xVals)];
                    y = [binPct(:,1) ; flip(binPct(:,2))];
                    c = [0.7 0.7 0.7];
                    if s.plt.meanLimbCols == 1 && s.plt.ON == 1
                        patch(x,y,c,'EdgeColor','None','FaceAlpha',.5); hold on
                        plot(xVals,binAMean,'LineWidth',2); %hold off
                    end
                elseif strcmp(s.plt.pltType,'Imagesc')

                    nA=squeeze(nanmean(allNeurAct(:,:,handRow,cHC,cL,cN),4));

                    if s.plt.ON == 1
                        ImagescInterp(nA,s.plt.intrpFact); colormap jet
                        ylim(s.plt.rowLims); xlim(s.plt.colLims);
                        hold on;

                        if s.plt.grid == 1
                            fS.gridXstart = s.plt.colLims(1) ;
                            fS.gridXstep = 1;
                            fS.gridYstart = s.plt.rowLims(1) ;
                            fS.gridYstep = 1;
                            GridOverImage(fS,gca);
                        end
                        if s.plt.axesVis == 0
                            set(gca,'visible','off');
                        end
                        if s.plt.showLimb == 1
                            % scatter(cHC,handRow,'filled','k');
                            scatter(cHC,handRow,'filled','MarkerEdgeColor',[0 0 0],...
                                'MarkerFaceColor',[.8 .8 .8],...
                                'LineWidth',2);
                        end
                        if ~isempty(s.plt.showBodyCol)
                            scatter(s.plt.showBodyCol,handRow + 1,'filled','MarkerEdgeColor',[0 0 0],...
                                'MarkerFaceColor',[.5 .5 .5],...
                                'LineWidth',2);
                        end

                        hold off;
                    end

                elseif strcmp(s.plt.pltType,'OnlyCalc')

                end

                % [TBD: Sort this out]
                try
                    nAtemp=nA(s.plt.stimRow,s.plt.stimCol,:);
                catch
                    nAtemp=nA;
                end
                try
                    [r1,p1,rL,rU] = corrcoef(distanceVals(:),nAtemp(:));
                    cov1 = nancov(distanceVals(:),nAtemp(:));
                    covVal(cL,cN)=cov1(2);
                    r2(cL,cN)=r1(2);
                    p2(cL,cN)=p1(2);
                catch
                    r2=NaN;
                    p2=NaN;
                end



            end
        end

        if s.plt.ON == 1
            colormap(redbluecmapRory(100,10));
        end

    end



end
