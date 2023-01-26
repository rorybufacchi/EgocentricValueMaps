function [] = SubplotSigTitles(pVal1,pVal2,names)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



for iPlot=1:length(pVal1(:))
    if pVal1(iPlot)<=0.05
        if pVal2(iPlot)<=0.05 
        subplot(size(pVal1,2),size(pVal1,1),iPlot);
        title([ names{1}  '& ' names{2}])
        else
        subplot(size(pVal1,2),size(pVal1,1),iPlot);
        title([ names{1}])
        end
    elseif pVal2(iPlot)<=0.05 
        subplot(size(pVal1,2),size(pVal1,1),iPlot);
        title([names{2}])
    elseif pVal2(iPlot)>0.05 
        subplot(size(pVal1,2),size(pVal1,1),iPlot);
        title('None')        
    end
end



end

