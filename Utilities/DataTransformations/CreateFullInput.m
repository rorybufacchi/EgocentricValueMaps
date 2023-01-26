function [inpAll] = CreateFullInput(s,inpOSV,extraInput)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[allCol,allGoalR,allGoalC]=meshgrid(1:s.wrld.size(2),1:s.wrld.size(1),1:s.wrld.size(2));

% all other state variables are summarised here - $$$ TO BE DONE $$$
if s.fl.thr==1
    oSV=inpOSV;
    
    aoS=oSV.*ones([size(allCol(:),1) 2]);
    
    if s.plt.plotThrFl==1
        allPosTemp=[allGoalR(:) allGoalC(:) aoS];
        aoS(:,1)=allPosTemp(:,1);
        aoS(:,2)=allPosTemp(:,2);
        allGoalR=allPosTemp(:,3);
        allGoalC=allPosTemp(:,4);
    end
else
    aoS=[];
end

% -------------------------------------------------------------------------
% Make 'visual' input
if s.fl.vis==1
    if s.fl.thr==0
        inpworlds=PosToVis(baseworld,repmat(handRow,size(allCol(:))),allCol(:),...
            allGoalR(:),allGoalC(:),1,1,s.plt.lmbVisVal,s.plt.golVisVal,s.plt.golVisVal,...
            s.fl.thr);
    else
        inpworlds=PosToVis(baseworld,repmat(handRow,size(allCol(:))),allCol(:),...
            allGoalR(:),allGoalC(:),aoS(:,1),aoS(:,2),s.plt.lmbVisVal,s.plt.golVisVal,s.plt.golVisVal,...
            s.fl.thr);
    end
    if s.fl.vpr==1
        inpAll=[allCol(:)' ; reshape(inpworlds,[size(inpworlds,1)*size(inpworlds,2) size(inpworlds,3) ]) ];
    else
        inpAll=[ reshape(inpworlds,[size(inpworlds,1)*size(inpworlds,2) size(inpworlds,3) ]) ];
    end
elseif s.fl.tch == 1
    inpAll=[allCol(:),allGoalR(:),allGoalC(:), aoS, touchVal.*ones(size(allGoalC(:)))]';
else
    inpAll=[allCol(:),allGoalR(:),allGoalC(:), aoS]';
end

% -------------------------------------------------------------------------
% Add body position in if necessary
if s.fl.bdyMov == 1
    inpAll=[inpAll; s.plt.bdyCol.*ones(size(allGoalC(:)))'];
end

% -------------------------------------------------------------------------
% Add tool informaiton if necessary
if s.fl.ToolChange==1
    inpAll=[inpAll; s.plt.ToolPresent.*ones(size(allGoalC(:)))'];
end

% -------------------------------------------------------------------------
% Add reaction time information if necessary
if s.fl.rtTask== 1
    inpAll=[inpAll; s.plt.rtTouchVal.*ones(size(allGoalC(:)))'];
end

% -------------------------------------------------------------------------
% Add extra input if necessary $$$$$ STILL WORK IN PROGRESS
if s.fl.extraInput==1
    inpAll=[inpAll; nan([size(allGoalC(:),1) size(extraInput,5)] )'];
    for iInp=1:size(inpAll,2)
        
        if s.fl.vis==1
            error('Extra input doesnt work with pixel input yet');
        else
            if s.xtra.useThr==1
                inpInp=inpAll([1 4 5], iInp);
            else
                inpInp=inpAll([1 2 3], iInp);
            end
        end
        
        inpAll(end-2:end,iInp) = ...
            squeeze(extraInput(1,inpInp(1),inpInp(2),inpInp(3),:));
    end
end
% $$$$$ HERE HERE squeeze(extraInput(1,inpAll(1,1),inpAll(2,1),inpAll(3,1),:  ))

% % %         for kAct=2:s.act.numA
% % %             if s.fl.vis==1
% % %                 inpAll=[inpAll   [inpAll(1:end-1,1:mel(allCol(:))) ] ];
% % %             else
% % %                 inpAll=[inpAll, [allCol(:),allGoalR(:),allGoalC(:),aoS]'];
% % %             end
% % %         end

% -------------------------------------------------------------------------
% Add history information if necessary (i.e. previous timepoint)
if s.fl.hist==1
% % %     if size(inpAll,1) < size(S,2)
        if s.fl.vis~=1
            if s.fl.mo==1
                inpAll=repmat(inpAll,[2 1]) ;
            else
                inpAll=[repmat(inpAll(1:end-1,:),[2 1]) ; inpAll(end,:)];
            end
            if s.plt.plotThrFl==0
                inpAll([2 ],:)=inpAll([2 ],:)-s.plt.rowLag;
                inpAll([3 ],:)=inpAll([3 ],:)-s.plt.colLag;
                if s.fl.thr==1
                    inpAll([4 ],:)=inpAll([4 ],:)-s.plt.nvRl; % $$$ NEED TO FIND OUT WHAT TO DO HERE <-- use harddrive to find out where nvrl comes from
                    inpAll([5 ],:)=inpAll([5 ],:)-s.plt.nvCl;
                end
            else
                inpAll([2 ],:)=inpAll([2 ],:)-s.plt.nvRl;
                inpAll([3 ],:)=inpAll([3 ],:)-s.plt.nvCl;
                inpAll([4],:)=inpAll([4],:)-s.plt.rowLag;
                inpAll([5],:)=inpAll([5],:)-s.plt.colLag;
            end
        else
            % If using history and vision,  first make new 'history' images
            if s.plt.plotThrFl==0
                inpworlds=PosToVis(baseworld,repmat(handRow,size(allCol(:))),allCol(:),...
                    allGoalR(:)-s.plt.rowLag,allGoalC(:)-s.plt.colLag,aoS(:,1)-s.plt.nvRl,aoS(:,2)-s.plt.nvCl,s.plt.lmbVisVal,s.plt.golVisVal,s.plt.golVisVal,...
                    s.fl.thr);
            else
                inpworlds=PosToVis(baseworld,repmat(handRow,size(allCol(:))),allCol(:),...
                    allGoalR(:)-s.plt.nvRl,allGoalC(:)-s.plt.nvCl,aoS(:,1)-s.plt.rowLag,aoS(:,2)-s.plt.colLag,s.plt.lmbVisVal,s.plt.golVisVal,s.plt.golVisVal,...
                    s.fl.thr);
            end 
            if s.fl.vpr==1
                prvInpAll=[allCol(:)' ; reshape(inpworlds,[size(inpworlds,1)*size(inpworlds,2) size(inpworlds,3) ]) ];
            else
                prvInpAll=[ reshape(inpworlds,[size(inpworlds,1)*size(inpworlds,2) size(inpworlds,3) ]) ];
            end
            
            inpAll=[prvInpAll ; inpAll];
            
        end
% % %     end
end


end

