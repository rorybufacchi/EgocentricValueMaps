function [net] = ChangeNeuronType(net,s)
%[net] = ChangeNeuronType(net,s)
%
%   Changes the type of neurons in the network

% Loop through layers and change transfer function for each
for iL = 1:length(net.layers)-1
    net.layers{iL}.transferFcn = s.lp.neurTypes;
end

end