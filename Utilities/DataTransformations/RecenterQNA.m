function [outVal] = RecenterQNA(s,inVal)
%[outVal] = RecenterQNA(inputArg1,inputArg2)
%   Shifts all the Q or Neural Activation values, so that the central-most
%   limb Column is the position of the actual limb. e.g. if normal Q has
%   limb Col == 2, for Q(:,2,:,:,:), then the outVal of this function has
%   limb Col == 8, by using circshift to move the limb column round
%
%   NOTE: remember to set s.plt.lmb


    centerVal=ceil(s.wrld.size(2)/2);
    for lmbCol=2:s.wrld.size(2)-1
        shiftAmount=centerVal - lmbCol;
        inVal(:,lmbCol,:,2:end-1,:)=circshift(inVal(:,lmbCol,:,2:end-1,:),shiftAmount,4);
    end
    outVal = inVal;


end

