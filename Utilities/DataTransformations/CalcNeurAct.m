function [neurAct2] = CalcNeurAct(s,net,inpAll,w)
%CalcNeurAct Calculates the activity of neurons in a network
%   [neurAct2] = CalcNeurAct(s,net,inpAll,w)

if s.fl.deepNet==1
    
    for kLay=1:size(s.lp.netS,2)
        
    tempNeurAct = squeeze(activations(net,reshape(inpAll,[size(inpAll,1) 1 1 size(inpAll,2)]),kLay+1));
    
    for kNeur=1:s.lp.netS(kLay)
        %                     neurAct2(:,:,:,1,kNeur,:)=permute(reshape(tempNeurAct(kNeur,:,1),[s.wrld.size(1) s.wrld.size(2) s.wrld.size(2) ]),[2 3 4 1]);
        neurAct2(:,:,:,kLay,kNeur,:)=reshape(tempNeurAct(kNeur,:),[s.wrld.size(1) s.wrld.size(2) s.wrld.size(2) ]);
    end
    
    end
    
else
    
    % Dimensions
    x1=inpAll;
    % IF there is only 0 input for x1, then add in a 1 somewhere so
    % that there are actual values going through instead of Infs
    % and resulting NaNs
    x1(min(x1==0,[],2),end)=1;
    % ^ Check above
    
    % x1_step1_xoffset = min(x1,[],2);
    x1_step1_xoffset=w.stateAttr.minS';
    
    % x1_step1_xmax = max(x1,[],2);
    x1_step1_xmax=w.stateAttr.maxS';
    
    x1_step1_gain = 2./(x1_step1_xmax -x1_step1_xoffset); tempinvGain=(1./x1_step1_xmax );
    x1_step1_gain(x1_step1_gain==Inf)=tempinvGain(x1_step1_gain==Inf);
    x1_step1_ymin = -1;
    % % %     y1_step1_ymin = -1;
    % % %     y1_step1_gain = 2/(max(y)-min(y));
    % % %     y1_step1_xoffset = min(y);
    % Input 1
    xp1All = MapminmaxApply(x1,x1_step1_gain,x1_step1_xoffset,x1_step1_ymin); %figure,plot(xp1(2,:),x(2,:),'x')
    
    % In case some of the variables had no variance in the training data, exclude them altogether
    if size(net.iw{1},2) ~= size(xp1All,1)
        warning('It looks like some variables had no variance in the training data.')
        xp1All([net.inputs{1}.processSettings{1}.remove],:)=[];
    end
    tempNeurAct(1:s.lp.netS(1),:,1) = ApplyTransferFcn(s,net,1,xp1All,[]);
% % %     tempNeurAct(1:s.lp.netS(1),:,1)=TanSigApply(net.iw{1}*xp1All+net.b{1});
    for kNeur=1:s.lp.netS(1)
        %                     neurAct2(:,:,:,1,kNeur,:)=permute(reshape(tempNeurAct(kNeur,:,1),[s.wrld.size(1) s.wrld.size(2) s.wrld.size(2) ]),[2 3 4 1]);
        neurAct2(:,:,:,1,kNeur,:)=reshape(tempNeurAct(kNeur,:,1),[s.wrld.size(1) s.wrld.size(2) s.wrld.size(2) ]);
    end
    for kLay=2:size(s.lp.netS,2)
        tempNeurAct(1:s.lp.netS(kLay),:,kLay) = ApplyTransferFcn(s,net,kLay,xp1All,tempNeurAct);
% % %         tempNeurAct(1:s.lp.netS(kLay),:,kLay)=TanSigApply(net.lw{kLay,kLay-1}*tempNeurAct(1:s.lp.netS(kLay-1),:,kLay-1)+net.b{kLay});
        for kNeur=1:s.lp.netS(kLay)
            %                         neurAct2(:,:,:,kLay,kNeur,:)=permute(reshape(tempNeurAct(kNeur,:,kLay),[s.wrld.size(1) s.wrld.size(2) s.wrld.size(2) ]),[2 3 4 1]);
            neurAct2(:,:,:,kLay,kNeur,:)=reshape(tempNeurAct(kNeur,:,kLay),[s.wrld.size(1) s.wrld.size(2) s.wrld.size(2) ]);
        end
    end
end



end


function [actOut] = ApplyTransferFcn(s,net,kLay,xp1All,tempNeurAct)

if kLay == 1
    actIn = net.iw{1}*xp1All+net.b{1};
else
    actIn = net.lw{kLay,kLay-1}*tempNeurAct(1:s.lp.netS(kLay-1),:,kLay-1)+net.b{kLay};
end

actOut = feval(net.layers{kLay}.transferFcn,actIn);


end

