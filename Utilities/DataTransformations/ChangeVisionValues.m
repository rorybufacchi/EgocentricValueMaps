%% Change reward visibility value

matS=cell2mat(storedEps.S);
matS(matS==100)=200;
tempS=mat2cell(matS,ones(size(matS,1),1));
storedEps.S(1:size(tempS,1))=tempS;


%% Change hand visibility value

matS=cell2mat(storedEps.S);
matS(matS==60)=100;
tempS=mat2cell(matS,ones(size(matS,1),1));
storedEps.S(1:size(tempS,1))=tempS;


%% Change threat visibility value

matS=cell2mat(storedEps.S);
matS(matS==10)=-200;
tempS=mat2cell(matS,ones(size(matS,1),1));
storedEps.S(1:size(tempS,1))=tempS;



%%
storedEps.R(storedEps.R==-2)=0;
