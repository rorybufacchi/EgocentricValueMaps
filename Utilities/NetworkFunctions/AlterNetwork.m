function net = AlterNetwork(net,aS)
% aS is the 'aterSettings'




% For 'lesion'/'ablation', need to set the output weights from a specific
% (or multiple) neurons to 0

if ~isfield(aS,'ListOrMat')
    aS.ListOrMat='Mat';
end

% either do as a matrix,
if strcmp(aS.ListOrMat,'Mat')
    
    % Ignore specific layers
    aS.lay(ismember(aS.lay,aS.ignoreLay))=[];
    
    switch aS.alterType
        case 'Ablate'
            net.lw{aS.lay+1,aS.lay}(:,aS.neur)=0;
        case 'Bias'
            net.b{aS.lay}(aS.neur)=net.b{aS.lay}(aS.neur)+aS.biasVals;
        otherwise
            error('Non-existent alteration type entered')
    end
    
    % or as a list
elseif strcmp(aS.ListOrMat,'List') & numel(aS.lay)==numel(aS.neur)
    
    % Ignore specific layers
    aS.neur(ismember(aS.lay,aS.ignoreLay))=[];
    aS.lay(ismember(aS.lay,aS.ignoreLay))=[];

    
    for iN=1:numel(aS.lay)
        tmpL=aS.lay(iN);
        tmpN=aS.neur(iN);
        
        if length(aS.biasVals)>1
            bV=aS.biasVals(iN);
        else
            bV=aS.biasVals;
        end
        
        switch aS.alterType
            case 'Ablate'
                net.lw{tmpL+1,tmpL}(:,tmpN)=0;
            case 'Bias'
                net.b{tmpL}(tmpN)=net.b{tmpL}(tmpN)+bV;
            % For this and the next, adapt the weights in the direction that they exist
            % already, so that the effect of the neuron - whatever
            % direction it is in - becomes stronger
            case 'DirectionalConnectionBias'
                for iNN = 1:size(net.lw{tmpL+1,tmpL},1) % loop through 'NextNeuron'
                    tmpSign = sign(net.lw{tmpL+1,tmpL}(iNN,tmpN));
                    net.lw{tmpL+1,tmpL}(iNN,tmpN) = ...
                        net.lw{tmpL+1,tmpL}(iNN,tmpN) + bV.*tmpSign;
                end
            case 'MultiplicativeDirectionalConnectionBias'
                for iNN = 1:size(net.lw{tmpL+1,tmpL},1) % loop through 'NextNeuron'
                    tmpSign = sign(net.lw{tmpL+1,tmpL}(iNN,tmpN));
                    net.lw{tmpL+1,tmpL}(iNN,tmpN) = ...
                        net.lw{tmpL+1,tmpL}(iNN,tmpN).* ( 1+ (bV.*tmpSign) );
                end
            case 'MultiplicativeConnectionBias'
                for iNN = 1:size(net.lw{tmpL+1,tmpL},1) % loop through 'NextNeuron'
                    net.lw{tmpL+1,tmpL}(iNN,tmpN) = ...
                        net.lw{tmpL+1,tmpL}(iNN,tmpN).* ( 1+ (bV) );
                end
            case 'MultiplicativeConnectionBias_WithinNet'
                
                if tmpL < length(net.iw)-1
                    % if it is not a last layer neuron, increase connection
                	% strengths to all other neurons WITHIN THE SUBNETWORK
                    nextNeurList = aS.neur(aS.lay == tmpL + 1);
                else
                    % If it is a last layer neuron, increase connection
                    % strength to the output directly
                    nextNeurList = 1:size(net.lw{tmpL+1,tmpL},1);
                end
                for iNN = nextNeurList
                    net.lw{tmpL+1,tmpL}(iNN,tmpN) = ...
                        net.lw{tmpL+1,tmpL}(iNN,tmpN).* ( 1+ (bV) );
                end

            otherwise
                error('Non-existent alteration type entered')
        end
        
    end
    
else
    error('wrong number of layers and neurons')
end


end




