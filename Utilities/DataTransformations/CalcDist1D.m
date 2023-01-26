function [cD,rD,aD]  = CalcDist1D(s,w,cHC,handRow)
% [cD,rD,aD]  = CalcDist1D(s,w,cHC,handRow)
% CalcResp1DCalculates the 1D distances
%   cHC is current hand column

% object columns
objCs = repmat([1:size(w.world2D,2)],[size(w.world2D,1) 1]);
objRs = repmat([1:size(w.world2D,1)]',[1 size(w.world2D,2)]);

 for iHC = 1:length(cHC)     

% col-distance
if strcmp(s.wrld.resetType,'BottomTop_InfLR')
    % If making infinite left-right world, the extreme columns should be
    % adjacent
    xtraObjCs(:,:,1) = objCs + objCs(1,end-2); % +2 because we skip the world edges
    xtraObjCs(:,:,2) = objCs - objCs(1,end-2);
    allDists = cHC(iHC)-cat(3,objCs,xtraObjCs);
    [dmycD minDistPos] = min(abs(allDists),[],3);
%     cD = max(cHC-cat(3,objCs,flip(objCs,2)),[],3);
    cD(:,:,iHC) = allDists(:,:,1);
    for iR = 1:size(cD,1)
        for iC = 1:size(cD,2)
            cD(iR,iC,iHC) = allDists(iR,iC,minDistPos(iR,iC));
        end
    end
else
    cD(:,:,iHC) = (cHC(iHC)-objCs);
end
% row-distance
rD(:,:,iHC) = (handRow-objRs);
% absolute distance
aD(:,:,iHC) = sqrt(cD(:,:,iHC).^2 + rD(:,:,iHC).^2);

 end

end

