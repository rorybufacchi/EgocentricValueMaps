function [w] = StoreStateAttributes(storedEps,w)
%StoreStateAttributes stores things like the spread of state values in w, so
%that neural activation functions can be correctly computed
%   $$$$ Detailed explanation goes here

w.stateAttr.minS=nanmin(cell2mat(storedEps.S),[],1);
w.stateAttr.maxS=nanmax(cell2mat(storedEps.S),[],1);

end

