function [outVar]= NanRemFlatten(inVar)
% Removes Nans and flattens
outVar=inVar(:);
outVar(isnan(outVar))=[];
end
