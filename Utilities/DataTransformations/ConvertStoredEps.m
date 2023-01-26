%% Convert to vision from position

handVisVal=60;
goalVisVal=100;
thrVisVal=10;

baseMaze=maze2D;
baseMaze( baseMaze~=0 ) = 50;


% ------------------------------------------------------------------------
matS=cell2mat(storedEps.S);
newSTmp=PosToVis(baseMaze,repmat(handRow,[size(matS,1) 1]),matS(:,1),...
    matS(:,2),matS(:,3),matS(:,4),matS(:,5),handVisVal,goalVisVal,thrVisVal,thrFlag);
clear newS
for kTrial=1:size(matS,1)
    tmpTmp=newSTmp(:,:,kTrial);
    newS(:,kTrial)=tmpTmp(:);
end
newS=mat2cell(newS',ones(size(newS',1),1));
storedEps.S(1:size(newS,1))=newS;
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Also for previous S
matS=cell2mat(storedEps.prvS);
newSTmp=PosToVis(baseMaze,repmat(handRow,[size(matS,1) 1]),matS(:,1),...
    matS(:,2),matS(:,3),matS(:,4),matS(:,5),handVisVal,goalVisVal,thrVisVal,thrFlag);

clear newS
for kTrial=1:size(matS,1)
    tmpTmp=newSTmp(:,:,kTrial);
    newS(:,kTrial)=tmpTmp(:);
end
newS=mat2cell(newS',ones(size(newS',1),1));
storedEps.prvS(1:size(newS,1))=newS;
% ------------------------------------------------------------------------

%% Convert to vis+propr from vision


handVisVal=60;
goalVisVal=100;
thrVisVal=10;
% thrVisVal=50;

clear newS
% ------------------------------------------------------------------------
matS=cell2mat(storedEps.S);
matS=reshape(matS',[mazeSize(1) mazeSize(2) size(matS,1)]);
[hR hC] =VisToPos(matS,handVisVal,goalVisVal,thrVisVal);
for kTrial=1:size(matS,3)
    tmpTmp=matS(:,:,kTrial);
    newS(:,kTrial)=tmpTmp(:);
end
newS=[hC; newS];
newS=mat2cell(newS',ones(size(newS',1),1));
storedEps.S(1:size(newS,1))=newS;
% ------------------------------------------------------------------------

clear newS
% ------------------------------------------------------------------------
matS=cell2mat(storedEps.prvS);
matS=reshape(matS',[mazeSize(1) mazeSize(2) size(matS,1)]);
[hR hC] =VisToPos(matS,handVisVal,goalVisVal,thrVisVal);
for kTrial=1:size(matS,3)
    tmpTmp=matS(:,:,kTrial);
    newS(:,kTrial)=tmpTmp(:);
end
newS=[hC; newS];
newS=mat2cell(newS',ones(size(newS',1),1));
storedEps.prvS(1:size(newS,1))=newS;
% ------------------------------------------------------------------------
