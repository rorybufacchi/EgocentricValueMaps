rafunction [figString, rhoVal, pVal, fullStats] = CorrelationTitle(xVals,yVals,opts)
%CORRELATIONTITLE Returns the correlation between two matrices, and also a
%string for use in a figure title
%
%[figString, rhoVal, pVal, fullStats] = CorrelationTitle(xVals,yVals,opts)
%
%   Inputs:
%   opts.testType:  determines what type of p-value is shown
%                   - standard is Pearsons's correlation
%                   - If lme is chosen, yVals ~ xVals + (1|colNum) is run,
%                     where colNum is the column of the x- and y-matrixes.
%                     So columns are the separate observations/subjects
%                   - if 'ttest', a t-test is performed. note that xVals
%                     should only take 2 values for a ttest
%                   - if 'signrank', a signed rank test is performed. note that xVals
%                     should only take 2 values for a ttest


if ~exist('opts','var')
    opts.testType = 'corr';
end

if isempty(xVals)
    xVals = repmat(1:size(yVals,2),[size(yVals,1) 1]);
end

switch opts.testType
    case 'corr'
        
        % Flatten and calculate correlation
        xVals = xVals(:);
        yVals = yVals(:);
        
        useVals = ~isnan(xVals) & ~isnan(yVals);
        
        [rhoVal, pVal]=corr(xVals(useVals), yVals(useVals));
        
        figString = ['RHO = ' num2str(rhoVal,3) ', P = ' num2str(pVal,3)];
        
    case 'lme'
        
        colNum = repmat([1:size(xVals,2)], [size(xVals,1) 1]);
%         colNum = repmat([1:size(xVals,2)]', [size(xVals,1) 1]);
%         colNum = repmat([1:size(xVals,2)], [size(xVals,1) 1]);
        
        % create table for LME
%         tbl = table(xVals(:),yVals(:),colNum(:));

        xVals = xVals(:); yVals = yVals(:); colNum = colNum(:);
        tbl = table(xVals,yVals,colNum);

        
        fullStats = fitlme(tbl,'yVals ~ xVals + (1|colNum)');

        % partial variance
        fixEffs     = fixedEffects(fullStats);
        partialY    = tbl.xVals .* fixEffs(2);
        varp        = nanvar(partialY);
        % partial eta square [effect size]
        pEta2       =  varp ./ (varp + nanvar(tbl.yVals));
        
        pVal = fullStats.coefTest;
        rhoVal = fullStats.Coefficients{2,2};
        
        figString = ['Coef = ' num2str(rhoVal,3) ' [' num2str(fullStats.Coefficients{2,7}) '-' num2str(fullStats.Coefficients{2,8}) '], P = ' num2str(pVal,3) '. partialEta2 = ' num2str(pEta2,3)];
        
    case 'signrank'
        
        % Flatten data
        xVals = xVals(:);
        yVals = yVals(:);
        
        xTypes = unique(xVals);
        
        if numel(xTypes) > 2 
            warning(['Signrank test requested, but more than 2 data ' ... 
                'regression points were given']);
        else
            [pVal hVal stats] = signrank(yVals(xVals == xTypes(1)), yVals(xVals == xTypes(2)));
            if isfield(stats,'zval')
                rhoVal      = stats.zval;
                effSizR1    = stats.zval ./ sqrt(numel(xVals.*2));
            else
                rhoVal = Inf;
            end
        end
        
        figString = ['Z-val = ' num2str(rhoVal,3) ', P = ' num2str(pVal,3) '. EffSize r = ' num2str(effSizR1,3)];
        
    case 'ttest'
        
        % Flatten data
        xVals = xVals(:);
        yVals = yVals(:);
        
        xTypes = unique(xVals);
        
        if numel(xTypes) > 2 
            warning(['t-test requested, but more than 2 data ' ... 
                'regression points were given']);
        else
            [hVal pVal ci stats] = ttest(yVals(xVals == xTypes(1)), yVals(xVals == xTypes(2)));
            rhoVal = stats.tstat;
        end
        
        figString = ['T-val = ' num2str(rhoVal,3) ', P = ' num2str(pVal,3)];
end


end

