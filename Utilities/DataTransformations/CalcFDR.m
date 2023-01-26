function [pValuesFDR fdrThresh]  =  CalcFDR(pValues,alpha, plotFlag)
%CALCFDR Calculates False-detection-rate corrected p-values
%
%   [pValuesFDR fdrThresh]  =  CalcFDR(pValues,alpha, plotFlag)
%
%   inputs:
%   - pValues:     original p-values. Can be vector, array or table
%   - alpha:       'significance' cutoff. default 0.05
%   - plotFlag:    1 if the corrected p-values are to be plotted. default 0

if istable(pValues)
    pVTbl = pValues;
    pValues = table2array(pValues);
end

if ~exist('alpha','var')
    alpha = 0.05;
end

if ~exist('plotFlag','var')
    plotFlag = 0;
end

numPvals = sum(~isnan(pValues(:)));

pValIncrement = alpha.*(1./numPvals);

pValComp = pValIncrement:pValIncrement:alpha;

[sortedpVals sortPos] = sort(pValues(:));

pValuesFDR = pValues;
pValuesFDR(sortPos( ...
    sortedpVals(~isnan(sortedpVals))>pValComp(:))  )  =  ...
    1;

% change matrix back into table if input was a table
if exist('pVTbl','var')
    pValuesFDR = array2table(pValuesFDR,...
        'VariableNames',pVTbl.Properties.VariableNames,...
        'RowNames',pVTbl.Properties.RowNames);
end

fdrThresh = max(pValuesFDR(pValuesFDR ~= 1));

if plotFlag == 1
    plot(sortedpVals,'x'); hold on; plot(pValComp); hold off
    xlabel('Index of sorted p-value')
    ylabel('p-value')
    legend('uncorrected p-values','correction cutoff','Location','northwest')
    title(['FDR Corrected p-cutoff = ' num2str(fdrThresh) ' (given \alpha = ' num2str(alpha) ')']) 
end

end

