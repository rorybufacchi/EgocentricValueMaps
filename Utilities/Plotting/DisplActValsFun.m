function [bestQ squishQs rDist pDist] = DisplActValsFun(s,w,Q)

s = DefaultSettings(s);

flipValFlag=0;

if isempty(s.plt.lmbRow)
    handRow=size(w.world2D,1)-2;
else
    handRow=s.plt.lmbRow;
end


% % % -------------------------------------------------------------------------
% % % Calulculate distance correlations
% % [cD,rD,aD] = CalcDist1D(s,w,s.plt.lmbCol,handRow);
% %
% % switch s.plt.distanceType
% %     case 'Absolute'
% %         distanceVals=aD;
% %     case 'Row'
% %         distanceVals=rD;
% %     case 'Column'
% %         distanceVals=cD;
% %     case 'AbsRow'
% %         distanceVals=abs(rD);
% %     case 'AbsColumn'
% %         distanceVals=abs(cD);
% % end
% % distanceVals=distanceVals(s.plt.stimRow,s.plt.stimCol,:);
% %
% % for iAct = 1:s.act.numA
% %     qTmp = squeeze(Q(handRow,s.plt.lmbCol,s.plt.stimRow,s.plt.stimCol,iAct));
% %
% %     [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(distanceVals(:),qTmp(:));
% %
% %     rDist(iAct)=rTmp(2);
% %     pDist(iAct)=pTmp(2);
% % end




% -------------------------------------------------------------------------
% Determine which limb positions to plot, and average accordingly
if s.plt.meanLimbCols==1
    centerVal=ceil(size(w.world2D,2)/2);
    for lmbCol=s.plt.lmbCol
        shiftAmount=centerVal - lmbCol;
        Q(:,lmbCol,:,2:end-1,:)=circshift(Q(:,lmbCol,:,2:end-1,:),shiftAmount,4);
    end
    Q=nanmean(Q(:,s.plt.lmbCol,:,:,:),2);
    sepCols=1;
else
    sepCols=s.plt.lmbCol;
end


% Plotting starts here ---------------------------------------------------

kAction=s.plt.plAct;

switch s.plt.pltType
    case 'Imagesc'
        for handCol=sepCols
            
            
            
            % Plot only one action value
            if s.plt.OneActFl==1
                
                rowVal=[1:size(w.world2D,1)];
                colVal=[1:size(w.world2D,2)];
                if s.plt.meanLimbCols==1
                    colVal=colVal-centerVal;
                end
                
                
                
                Qtemp=squeeze(Q(handRow,handCol,:,:,kAction));
                squishQs=Qtemp;
                if s.plt.ON==1
                    if flipValFlag==1
                        ImagescInterp(-squeeze(Q(handRow,handCol,:,:,kAction)),s.plt.intrpFact,colVal,rowVal); if s.plt.cBarFl==1, colorbar; end
                    else
                        ImagescInterp(squeeze(Q(handRow,handCol,:,:,kAction)),s.plt.intrpFact,colVal,rowVal); if s.plt.cBarFl==1, colorbar; end
                    end
                    if s.plt.axesVis==0
                        set(gca,'visible','off')
                    end
                    % caxis([-1 1])
                    %     title(['Q(' s.act.Name{kAction} '), NN'])
                    title(['Q(' s.act.Name{kAction} ')'])
                    ylim(s.plt.rowLims);
                    if s.plt.meanLimbCols==1
                        xlim(colVal(1)+s.plt.colLims+[0 -1]);
                    else
                        xlim(colVal(1)+s.plt.colLims);
                    end

                    % Plot contour lines if requested
                    if s.plt.contours == 1
                        tmpQ = squeeze(Q(handRow,handCol,:,:,kAction));
                        % Create x- and y- axes
                        [X, Y] = meshgrid(colVal, rowVal);
 
                        % Up-sample data using interp2
                        factor = 2; % Up-sample factor
                        [Xq, Yq] = meshgrid(linspace(min(X(:)), max(X(:)), factor * size(tmpQ, 2)), linspace(min(Y(:)), max(Y(:)), factor * size(tmpQ, 1)));
                        tmpQ2 = interp2(X,Y,tmpQ, Xq, Yq, 'cubic');

                        % Overlay contour plot
                        hold on;
                        [C, h] = contour(Xq, Yq, tmpQ2, s.plt.contourVals, 'k','LineWidth',2); % '10' is the number of contour levels, adjust as needed
                        % clabel(C, h); % Optional: Add labels to the contour levels
                        hold off;
                    end

                    drawnow
                    
                    colormap jet
                end
                squishQs=squeeze(Q(handRow,handCol,:,:,:));
                bestQ=[];
            else
                
                x=1:size(w.world2D,1);
                y=1:size(w.world2D,2);
                [X,Y]=meshgrid(x,y);
                X=X; Y=Y;
                
                squishQs=squeeze(Q(handRow,handCol,:,:,:));
                
                bestQ=squeeze(Q(handRow,handCol,:,:,3))-squeeze(Q(handRow,handCol,:,:,1));
                QmaxTemp=max(squeeze(Q(handRow,handCol,:,:,[1 3])),[],3);
                Qstay=squeeze(Q(handRow,handCol,:,:,2))>QmaxTemp;
                bestQ(Qstay)=0;
                
                % % %         [maxVal bestQ]=max(Q(handRow,handCol,:,:,:),[],5);
                % % %         bestQ=squeeze(bestQ);
                % % % % %
                % % % % %         axis xy
                
                if s.plt.ON==1
                    imagesc(bestQ); axis xy; hold on
                    
                    quiver(Y,X,bestQ',zeros(size(bestQ))',0.8,'w','LineWidth',1.2);  hold on;
                    
                    [r c]=find(Qstay);
                    scatter(c,r,'o','filled'); hold off
                    
                    if s.plt.axesVis==0
                        set(gca,'visible','off')
                    end
                    
                    view(180,90);
                end
                
                
                
            end
            
            
        end
        
        
    case 'Binned'
        
        Qtemp=squeeze(Q(handRow,sepCols,s.plt.stimRow,s.plt.stimCol,kAction));
        if s.plt.meanLimbCols == 1
            [cD,rD,aD] = CalcDist1D(s,w,centerVal,handRow);
        else
            [cD,rD,aD] = CalcDist1D(s,w,sepCols,handRow);
        end
        
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
        
        BinPlot(distanceVals,Qtemp,s.plt.fS,s.plt.varargs{:})
        % %             BINPLOT HERE
        % % %     % THIS IS OLD, keep FOR LEGACY PURPOSES. Remove later
        % % %     else
        % % %
        % % %         %     set(0,'CurrentFigure',hActVal)
        % % %         % figure,
        % % %         for kAction=1:s.act.numA
        % % %
        % % %             try
        % % %                 subplot(2,s.act.numA*2,kAction)
        % % %                 ImagescInterp(squeeze(Qtable(handRow,handCol,:,:,kAction)),s.plt.intrpFact); if s.plt.cBarFl==1, colorbar; end
        % % %                 % caxis([-1 1])
        % % %                 title(['Q(' s.act.Name{kAction} '), table'])
        % % %                 ylim(rowLims);
        % % %                 drawnow
        % % %             end
        % % %
        % % %             subplot(2,s.act.numA*2,s.act.numA + kAction )
        % % %             ImagescInterp(squeeze(Q(handRow,handCol,:,:,kAction)),s.plt.intrpFact); if s.plt.cBarFl==1, colorbar; end
        % % %             % caxis([-1 1])
        % % %             title(['Q(' s.act.Name{kAction} '), NN'])
        % % %             ylim(rowLims);
        % % %             drawnow
        % % %
        % % %         end
        % % %
        % % %         try
        % % %             subplot(2,s.act.numA*2,s.act.numA*2 + 2)
        % % %             pickRight=Qtable(:,:,:,:,rVal)-Qtable(:,:,:,:,lVal);
        % % %             ImagescInterp(squeeze(pickRight(handRow,handCol,:,:,1)),s.plt.intrpFact); if s.plt.cBarFl==1, colorbar; end
        % % %             caxis([-1 1])
        % % %             title('TABLE: RIGHT = 1, LEFT = 0.');
        % % %             colLims=max(abs(caxis));
        % % %             caxis([-colLims colLims]);
        % % %             ylim(rowLims);
        % % %             drawnow
        % % %         end
        % % %
        % % %         subplot(2,s.act.numA*2,s.act.numA*4 - 2)
        % % %         pickRight=Q(:,:,:,:,rVal)-Q(:,:,:,:,lVal);
        % % %         ImagescInterp(squeeze(pickRight(handRow,handCol,:,:,1)),s.plt.intrpFact); if s.plt.cBarFl==1, colorbar; end
        % % %         colLims=max(abs(caxis));
        % % %         caxis([-colLims colLims]);
        % % %         ylim(rowLims);
        % % %         title('NN: RIGHT +ve, LEFT -ve ');
        % % %         colormap jet
        % % %         drawnow
        % % %
        % % %         subplot(2,s.act.numA*2,s.act.numA*4-1)
        % % %         pickRight=Q(:,:,:,:,rVal)>Q(:,:,:,:,lVal);
        % % %         ImagescInterp(squeeze(pickRight(handRow,handCol,:,:,1)),s.plt.intrpFact); if s.plt.cBarFl==1, colorbar; end
        % % %         colLims=max(abs(caxis));
        % % %         caxis([-colLims colLims]);
        % % %         ylim(rowLims);
        % % %         title('NN: RIGHT = 1, LEFT = 0.');
        % % %         colormap jet
        % % %         drawnow
        % % %
        % % %         subplot(2,s.act.numA*2,s.act.numA*4)
        % % %         [maxValue maxAction]=max(Q(:,:,:,:,:),[],5);
        % % %         ImagescInterp(squeeze(maxAction(handRow,handCol,:,:,1)),s.plt.intrpFact); if s.plt.cBarFl==1, colorbar; end
        % % %         caxis([1 s.act.numA]);
        % % %         ylim(rowLims);
        % % %         title(['NN: LEFT = 1, STAY = 2, RIGHT =3 ']);
        % % %         colormap jet
        % % %         drawnow
        % % %
        % % %         subplot(2,s.act.numA*2,s.act.numA*2 + 1)
        % % %         handPoss=zeros(size(w.world2D));
        % % %         handPoss(handRow,handCol)=1;
        % % %         ImagescInterp(handPoss,s.plt.intrpFact);
        % % %         caxis([-1 1]);
        % % %         ylim(rowLims);
        % % %         %     title(['Hand Position, actionNum=' num2str(countActions)]);
        % % %         title(['Hand Position']);
        % % %         if s.plt.cBarFl==1, colorbar; end
        % % %         drawnow
        % % %
        % % %     end
end


%%

% $$$ OK THIS IS A PROBLEM that its' been putting things onto handrow 19
% even though that's wrong --> find out why
% ImagescInterp(squeeze(nanmean(Qtable(14,9,:,:,1),[5 2])),s.plt.intrpFact); if s.plt.cBarFl==1, colorbar; end;
