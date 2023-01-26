function [rMat] = DisplayProxStats(rS)
%DisplayProxStats Displays and stores stats for all input model results
%
%   Inputs:
%   rS:     Structure of results from different models
%
%   rMat:   Strucure of all the most extreme stats, for reporting

for iM = 1:length(rS)
    
    % expand the matrices if necessary
    if iM >1 && size(pDistQ,1)<size(rS(iM).pDistQ.max,1)
        pDistQ = cat(1,pDistQ,nan([size(rS(iM).pDistQ.max,1) - size(pDistQ,1) , size(pDistQ,2)]));
        rDistQ = cat(1,rDistQ,nan([size(rS(iM).rDistQ.maxP,1) - size(rDistQ,1) , size(rDistQ,2)]));
    end
    
    pDistQ(1:length(rS(iM).pDistQ.max),iM) = rS(iM).pDistQ.max;
    rDistQ(1:length(rS(iM).rDistQ.maxP),iM) = rS(iM).rDistQ.maxP;
    % Put Nans in if necessary
    pDistQ(length(rS(iM).pDistQ.max)+1:end,iM) = nan;
    rDistQ(length(rS(iM).rDistQ.maxP)+1:end,iM) = nan;
    
    % Only do neural stats if they have been calculated
    if isfield(rS(iM),'pDistN')
        
    
    nLay = length(rS(iM).s.lp.netS);
    maxN = max(rS(iM).s.lp.netS) ; % maximum neurons in layer
    
    if iM > 1
        % Reshape matrices if necessary
        if nLay > size(pDistN,1)
            pDistN = cat(1,pDistN,nan([nLay - size(pDistN,1), size(pDistN,2), size(pDistN,3)])  );
            rDistN = cat(1,rDistN,nan([nLay - size(rDistN,1), size(rDistN,2), size(rDistN,3)])  );
            hProxN = cat(1,hProxN,nan([nLay - size(hProxN,1), size(hProxN,2), size(hProxN,3)])  );
            propCorrNeurPerLay = cat(1,propCorrNeurPerLay,nan([nLay - size(propCorrNeurPerLay,1), size(propCorrNeurPerLay,2)])  );
        end
        if maxN > size(pDistN,2)
            pDistN = cat(2,pDistN,nan([max([nLay size(pDistN,1)]), maxN - size(pDistN,2), size(pDistN,3)]));
            rDistN= cat(2,rDistN,nan([max([nLay size(rDistN,1)]), maxN - size(rDistN,2), size(rDistN,3)]));
            hProxN = cat(2,hProxN,nan([max([nLay size(hProxN,1)]), maxN - size(hProxN,2), size(hProxN,3)]));
        end
    end
    
       pDistN(1:nLay ,1:maxN,iM) = rS(iM).pDistN.max(:,1:max(rS(iM).s.lp.netS));
        rDistN(1:nLay ,1:maxN,iM) = rS(iM).rDistN.maxP(:,1:max(rS(iM).s.lp.netS));
        hProxN(1:nLay ,1:maxN,iM) = double(rS(iM).hProxN(:,1:max(rS(iM).s.lp.netS)));
        
        % proportion of neurons which correlate, by layer
        propCorrNeurPerLay(1:length(rS(iM).s.lp.netS),iM) = sum(rS(iM).hProxN,2)./(rS(iM).s.lp.netS)';
        
        % proportion of neurons which correlate, across all layers
        propCorrNeur(iM) = sum(rS(iM).hProxN,'all')./sum(rS(iM).s.lp.netS);
        
        % Replace zeros with NaNs where there are no neurons
        for iL = 1:length(rS(iM).s.lp.netS)
            pDistN(iL,(rS(iM).s.lp.netS(iL)+1):end,iM) = NaN;
            hProxN(iL,(rS(iM).s.lp.netS(iL)+1):end,iM) = NaN;
        end
    end
end

% Only do neural stats if they have been calculated
if isfield(rS(iM),'pDistN')
    
    % Relative depth of a layer - used for calculating the correlation between
    % neuron layer depth and correlation with limb position
    relLayDep = nan(size(propCorrNeurPerLay));
    for iM = 1:length(rS)
        nNeur = length(rS(iM).s.lp.netS);
        relLayDep(1:nNeur,iM) = linspace(0,1,nNeur);
    end
    
    opts.testType = 'lme';
    % [figString rhoVal pVal] = CorrelationTitle(relLayDep(:),propCorrNeurPerLay(:),opts);
    [figString rhoVal pVal] = CorrelationTitle(relLayDep,propCorrNeurPerLay,opts);
    
    figure, plot (relLayDep, propCorrNeurPerLay,'-o','LineWidth',2);
    title(figString);
    xlabel('Relative layer depth')
    ylabel('Proportion of neurons with PPS properties')
    
    disp(['Relationship between neural layer PPS-ness: ' figString]);disp(' ')
    
    
    %results matrix
    rMat.pDistN = pDistN;
    rMat.rDistN = rDistN;
    rMat.hProxN = hProxN;
    rMat.relLayDep = relLayDep;
    rMat.propCorrNeurPerLay = propCorrNeurPerLay;
    rMat.propCorrNeur = propCorrNeur;
    

end

pDistQOrg    = pDistQ(:); % store non-corrected values as well
[dmy pDistQ] = fdr(pDistQ(:));

[maxP maxInd] = max(pDistQ(:));
disp(['PPS-ness of action value. Maximum p for Q-values: P < =  ' num2str(maxP,3) '. |rho| > =  |' num2str(rDistQ(maxInd),3) '|' ])

% extract mean rho vals and confidence intervals
allR = [];
for iM = 1:length(rS)
    allR = [allR; rS(iM).rDistQ.cDabs(:); rS(iM).rDistQ.rDabs(:); rS(iM).rDistQ.aD(:)];
end
meanR = nanmean(abs(allR));
ciR   = prctile(abs(allR),[5 95]); % confidence intervals
stdR  = nanstd(allR); % standard dev

disp(['MEAN |rho| =  ' num2str(meanR,3) ' STD |rho| =  ' num2str(stdR,3) '. 95% CI |rho| =  ' num2str(ciR(1),3) ' - ' num2str(ciR(2),3) ])
disp(' ')
disp(['Uncorrected P < =  ' num2str(max(pDistQOrg(:)),3) '' ])
disp(' ')

    

%results matrix
rMat.pDistQOrg  = pDistQOrg;
rMat.pDistQ     = pDistQ;
rMat.rDistQ     = rDistQ;
rMat.maxP       = maxP;

end

