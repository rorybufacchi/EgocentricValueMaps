function [Q,allNeurAct] = CalcNetOutput(s,w,net,extraInput)
%CalcNetOutput Find value and neural output of a trained network
%
%   [Q,allNeurAct] = CalcNetOutput(s,w,net,extraInput)
%   
%   Output:
%   - Q             Action values. Dimensions:   Lmb Row, Lmb Col, Obj Row, Obj Col, layer, neuron
%   - AllNeurAct    Neural activation Dimensions:Obj Row, Obj Col, LMb Row, Lmb Col, layer, neuron

addpath('D:\Old_D\DPPS\DefenseAgent\Scripts')

clear tempQ2

% $$ Maybe also try making thrat a negative visibility value
% handVisVal=60;
% goalVisVal=200;
% thrVisVal=-200;

% Number of other state variables - 1 if not taking the mean
if s.plt.meanOSVfl==0
    nOSV=[1 1];
else
    nOSV=s.wrld.size-2;
end

if isempty(s.plt.lmbRow);
    s.plt.lmbRow=size(w.world2D,1)-2;
end


try
    if s.fl.tch==1
        touchVal=zeros(size(touchV));
    end
catch
    warning('s.fl.tch might not exist. Setting it to 0')
    s.fl.tch=0;
end

for kRow=1:nOSV(1)
    for kCol=1:nOSV(2)
        
        if s.plt.meanOSVfl==0
            % other state variables - only used if s.plt.meanOSVfl is 0
            inpOSV=[s.plt.lmbRow-6 5];
            % inpOSV=[4 4];
        else
            inpOSV=[kRow kCol]+1
        end
        
        % -----------------------------------------------------------------
        % Create fake input
        if ~exist('extraInput','var'), extraInput=[]; end
        [inpAll] = CreateFullInput(s,inpOSV,extraInput);
        
        
        % If not plotting neurons, create q-value output
        if s.fl.mo==1
            if s.fl.deepNet==1
                tempQ=predict(net,reshape(inpAll,[size(inpAll,1) 1 1 size(inpAll,2)]))';
            else
                tempQ=net(inpAll);
            end
            tempQ2(:,:,:,:)=permute(reshape(tempQ,[size(tempQ,1) s.wrld.size(1) s.wrld.size(2) s.wrld.size(2)]),[2 3 4 1]); % $$$ HERE - this is a guess right now
        else
            tempQ=net(inp);
            tempQ2(:,:,:,kAct)=reshape(tempQ,[s.wrld.size(1) s.wrld.size(2) s.wrld.size(2)]);
        end
        
        
        
        [neurAct2] = CalcNeurAct(s,net,inpAll,w);
        
        for iLayer=1:size(neurAct2,4);
            for iNeur=1:size(neurAct2,5);
                
                tempNeurAct(:,:,:,:)=repmat(squeeze(neurAct2(:,:,:,iLayer,iNeur,:)),[1 1 1 size(tempQ2,4)]);
                
                % Not sure why this is correct but apparently % it is...
                tempNeurAct2=permute(tempNeurAct,[5 2 1 3 4]);
                if s.plt.plotThrFl==1 % $$$$ NOT SURE ABOUT THIS
                    tempNeurAct2=permute(tempNeurAct2,[1 2 3 4 5]);
                end
                tempNeurAct2=repmat(tempNeurAct2,[size(w.world2D,1) 1 1 1 1]);
                
                for hRow=1:size(neurAct2,1);
                    for hCol=1:size(neurAct2,2);
                        allNeurAct(:,:,hRow,hCol,iLayer,iNeur)=squeeze(tempNeurAct2(hRow,hCol,:,:,s.plt.plAct));
                    end
                end
                
            end
        end
        
        % Not sure why this is correct but apparently % it is...
        QQ=permute(tempQ2,[5 2 1 3 4]);
        if s.plt.plotThrFl==1 % $$$$ NOT SURE ABOUT THIS
            QQ=permute(QQ,[1 2 3 4 5]);
        end
        QQ=repmat(QQ,[s.wrld.size(1) 1 1 1 1]);
        
        % Obj Row, Obj Column, Hand Row, Hand Column, layer, neuron
        % Replace all zero values with NaN
        allNeurAct(allNeurAct==0)=NaN;
        
        
        
        if s.plt.meanOSVfl==0
            Q=QQ;
        else
            Qall(:,:,:,:,:,kRow,kCol)=QQ;
        end
        
    end
end

if s.plt.meanOSVfl==1
    Q=nanmean(nanmean(Qall,6),7);
end

end

