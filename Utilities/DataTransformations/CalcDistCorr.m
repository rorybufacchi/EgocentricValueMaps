function [rDist, pDist, hProx, aD, rD, cD, sigmOut, sigMid, rUDist, rLDist] = CalcDistCorr(s,w,Q)
%CalcDistCOrr Calculates whether a set of Q values are correlated with
%various measures of distance
%
%   [rDist, pDist, hProx, aD, rD, cD] = CalcDistCorr(s,w,Q)
%
%   outputs:
%   - rDist:    rho values
%   - pDist:    p values
%   - hProx:    rho values


if isempty(s.plt.lmbRow)
    handRow=size(w.world2D,1)-2;
else
    handRow=s.plt.lmbRow;
end
if s.plt.DistFromTool == 1
    handRow = handRow - max(s.lmb.ToolRows);
end

if ~isfield(s.plt,'fitSigmDistTypes')
    distTypes={'aD','rD','cD','rDabs','cDabs'};
else
    distTypes=s.plt.fitSigmDistTypes;
end

[cD,rD,aD] = CalcDist1D(s,w,s.plt.lmbCol,handRow);

% loop through all distance types
for iDT = 1:length(distTypes)
    
    dT = distTypes{iDT};
    
    switch dT
        case 'aD'
            distanceVals=aD;
        case 'rD'
            distanceVals=rD;
        case 'cD'
            distanceVals=cD;
        case 'rDabs'
            distanceVals=abs(rD);
        case 'cDabs'
            distanceVals=abs(cD);
    end
    
    distanceVals=distanceVals(s.plt.stimRow,s.plt.stimCol,:);
    
    rDist.(dT)  = nan([size(Q,5) size(Q,6)]);
    rLDist.(dT) = nan([size(Q,5) size(Q,6)]);
    rUDist.(dT) = nan([size(Q,5) size(Q,6)]);
    pDist.(dT)  = nan([size(Q,5) size(Q,6)]);
    sigMid.(dT) = nan([size(Q,5) size(Q,6)]);
    
    for iAct = 1:size(Q(:,:,:,:,:),5) % so that it can be done with neural activation as well
        % % %         qTmp = squeeze(Q(handRow,s.plt.lmbCol,s.plt.stimRow,s.plt.stimCol,iAct));
        % % %
        % % %         %match the dimensions in the right order
        % % %         qTmp = permute (qTmp,[2 3 1]);
        
        % -------------------------------------------------------------------------
        % NOT SURE IF THIS IS RIGHT YET --> I think it works for the valence setup,
        % but not sure about the others
        qTmp = Q(handRow,s.plt.lmbCol,s.plt.stimRow,s.plt.stimCol,iAct);
        
        %match the dimensions in the right order
        qTmp = permute (qTmp,[3 4 2 1 5]);
        
        
        
        
        [rTmp,pTmp,rLTmp,rUTmp] = corrcoef(distanceVals(:),double(qTmp(:)));
        
        rDist.(dT)(iAct)    = rTmp(2);
        pDist.(dT)(iAct)    = pTmp(2);
        rLDist.(dT)(iAct)   = rLTmp(2);
        rUDist.(dT)(iAct)   = rUTmp(2);
        
        
        if s.plt.fitSigmoid == 1 
            %------------------------------------------------------------------
            % Fit a sigmoid as well
            X = distanceVals(:);
            Y = qTmp(:);
            
            % If the absolute distance is taken, the zero point has got to
            % be ignored (cause you can't have a stimulus there, and there
            % is only one repetition of that point, making it a bias in the
            % wrong direction
            switch dT
                case 'aD'
                    X(distanceVals(:) == 0) = [];
                    Y(distanceVals(:) == 0) = [];
            end
            
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',s.plt.fitSigmLow,...
                'Upper',s.plt.fitSigmUp, ...
                'StartPoint',s.plt.fitSigmStart);
            
            ft = fittype('a/ (1+exp(-b*(x-c))) + d' , 'options', fo);
            if ~min(isnan(Y))
                [sigmOut.(dT)(iAct).curve,sigmOut.(dT)(iAct).gof] = fit(X,Y,ft);
                tmp = coeffvalues(sigmOut.(dT)(iAct).curve);
                sigMid.(dT)(iAct) = tmp(3);
            else
                sigmOut.(dT)(iAct).curve = [];
                sigmOut.(dT)(iAct).gof = [];
                sigMid.(dT)(iAct) = NaN;
            end
            
% %             
% %                             figure,
% % %                     if iAct == 2 && strcmp(dT, 'rD')
% %                         plot(sigmOut.(dT)(iAct).curve); hold on;
% %                         scatter(X(:),Y(:))
% % %                     end
        else
            sigmOut = [];
            sigMid  = [];
        end
    end
    
    
end

% Just catch errors in case not all the necessary distance values are
% searched for
try

% % hProx = CalcFDR(pDist.cDabs)<0.05 & CalcFDR(pDist.rD)<0.05 & CalcFDR(pDist.aD)<0.05;
% % % find the higest p-value for reporting
% % allP = cat(3,pDist.cDabs, pDist.rD, pDist.aD);
% % and the corresponding rho/r value
% allR = cat(3,rDist.cDabs, rDist.rD, rDist.aD);

hProx = CalcFDR(pDist.cDabs,0.05,0)<0.05 & CalcFDR(pDist.rDabs,0.05,0)<0.05 & CalcFDR(pDist.aD,0.05,0)<0.05;
% find the higest p-value for reporting
allP = cat(3,pDist.cDabs, pDist.rDabs, pDist.aD);
% and the corresponding rho/r value
allR = cat(3,rDist.cDabs, rDist.rDabs, rDist.aD);

[pDist.max maxInd]= max(allP,[],3);


iI=1;
for iR = 1:size(allR,1)
    for iC = 1:size(allR,2)
        rDist.maxP(iR,iC) = allR(iR,iC,maxInd(iI));
        iI=iI+1;
    end
end

catch
    warning('Problem in assigning "hProx" and similar variables. If necessary,')
    warning('ensure that "s.plt.fitSigmDistTypes" includes cDabs, rDabs and aD')
    
    rDist = []; pDist = []; hProx = [];
end



end

