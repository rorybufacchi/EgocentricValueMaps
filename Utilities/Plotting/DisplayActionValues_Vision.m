addpath('C:\Users\Rory\Documents\DPPS\DefenseAgent\Scripts')

lVal=1; mVal=2; rVal=3; 

handVisVal=60;
goalVisVal=100;
thrVisVal=10;
% thrVisVal=50;

% $$ Maybe also try making thrat a negative visibility value
% handVisVal=60;
% goalVisVal=200;
% thrVisVal=-200;

plotNeurFlag=1;

plotThrFlag=0;

meanOSVFlag=0;

plot1ActFlag=1;
plAction=1;

% cAxVals=[0 1];
cAxVals=[];

% currNeur=12;
% currLayer=1;

% currNeur=4;
% currLayer=2;

currNeur=2;
currLayer=3;

% currNeur=5;
% currLayer=8;

handRow=size(maze2D,1)-2;
handCol=4;

if meanOSVFlag==0
    nOSV=[1 1];
else
    nOSV=mazeSize-2;
end

clear neurAct2 tempNeurAct

baseMaze=maze2D;
baseMaze( baseMaze~=0 ) = 50;

[allCol,allGoalR,allGoalC]=meshgrid(1:mazeSize(2),1:mazeSize(1),1:mazeSize(2));

rowLims=[1 mazeSize(1)-1];

rowLag=1;
colLag=0;

% non-visible object lag (i.e. thr if goal is vis and vice-versa)
nvRl=1;
nvCl=0;

for kRow=1:nOSV(1)
    for kCol=1:nOSV(2)
if meanOSVFlag==0
% other state variables - only used if meanOsvFlag is 0
inpOSV=[handRow-6 5];
% inpOSV=[4 4];
else
    inpOSV=[kRow kCol]+1
end

% all other state variables are summarised here - $$$ TO BE DONE $$$
if thrFlag==1
    oSV=inpOSV; 
    
    aoS=oSV.*ones([size(allCol(:),1) 2]);
    
    if plotThrFlag==1
        allPosTemp=[allGoalR(:) allGoalC(:) aoS];
        aoS(:,1)=allPosTemp(:,1);
        aoS(:,2)=allPosTemp(:,2);
        allGoalR=allPosTemp(:,3);
        allGoalC=allPosTemp(:,4);
    end
else
    aoS=[];
end

% Make 'visual' input 
if visFlag==1
    
    inpMazes=PosToVis(baseMaze,repmat(handRow,size(allCol(:))),allCol(:),...
        allGoalR(:),allGoalC(:),aoS(:,1),aoS(:,2),handVisVal,goalVisVal,thrVisVal,...
        thrFlag);
    if visPropFlag==1
        inpAll=[allCol(:)' ; reshape(inpMazes,[size(inpMazes,1)*size(inpMazes,2) size(inpMazes,3) ]) ; ones(size(allCol(:)))'];
    else
        inpAll=[ reshape(inpMazes,[size(inpMazes,1)*size(inpMazes,2) size(inpMazes,3) ]) ; ones(size(allCol(:)))'];
    end
else
    inpAll=[allCol(:),allGoalR(:),allGoalC(:), aoS, ones(size(allCol(:)))]';
end
for kAct=2:possActions
    if visFlag==1
        inpAll=[inpAll   [inpAll(1:end-1,1:numel(allCol(:))) ; kAct*ones(size(allCol(:)))' ] ];
    else
        inpAll=[inpAll, [allCol(:),allGoalR(:),allGoalC(:),aoS, kAct*ones(size(allCol(:)))]'];
    end
end

if histFlag==1
    if size(inpAll,1) < size(S,2)
        if visFlag~=1
            inpAll=[repmat(inpAll(1:end-1,:),[2 1]) ; inpAll(end,:)];
            if plotThrFlag==0
                inpAll([2 ],:)=inpAll([2 ],:)-rowLag;
                inpAll([3 ],:)=inpAll([3 ],:)-colLag;
                inpAll([4 ],:)=inpAll([4 ],:)-nvRl;
                inpAll([5 ],:)=inpAll([5 ],:)-nvCl;
            else
                inpAll([2 ],:)=inpAll([2 ],:)-nvRl;
                inpAll([3 ],:)=inpAll([3 ],:)-nvCl;
                inpAll([4],:)=inpAll([4],:)-rowLag;
                inpAll([5],:)=inpAll([5],:)-colLag;
            end
        else
            % If using history and vision,  first make new 'history' images
            if plotThrFlag==0
                inpMazes=PosToVis(baseMaze,repmat(handRow,size(allCol(:))),allCol(:),...
                    allGoalR(:)-rowLag,allGoalC(:)-colLag,aoS(:,1)-nvRl,aoS(:,2)-nvCl,handVisVal,goalVisVal,thrVisVal,...
                    thrFlag);
            else 
                inpMazes=PosToVis(baseMaze,repmat(handRow,size(allCol(:))),allCol(:),...
                    allGoalR(:)-nvRl,allGoalC(:)-nvCl,aoS(:,1)-rowLag,aoS(:,2)-colLag,handVisVal,goalVisVal,thrVisVal,...
                    thrFlag);
            end %%%%
            if visPropFlag==1
                prvInpAll=[allCol(:)' ; reshape(inpMazes,[size(inpMazes,1)*size(inpMazes,2) size(inpMazes,3) ]) ];
            else
                prvInpAll=[ reshape(inpMazes,[size(inpMazes,1)*size(inpMazes,2) size(inpMazes,3) ]) ];
            end
            
            inpAll=[repmat(prvInpAll,[1 possActions]) ; inpAll];
            
        end
    end
end


if plotNeurFlag==1
    % Dimensions
    x1=inpAll;
    x1_step1_xoffset = min(x1,[],2);
    x1_step1_gain = 2./(max(x1,[],2)-min(x1,[],2)); tempinvGain=(1./max(x1,[],2));
    x1_step1_gain(x1_step1_gain==Inf)=tempinvGain(x1_step1_gain==Inf)
    x1_step1_ymin = -1;
% % %     y1_step1_ymin = -1;
% % %     y1_step1_gain = 2/(max(y)-min(y));
% % %     y1_step1_xoffset = min(y);
    % Input 1
    xp1All = MapminmaxApply(x1,x1_step1_gain,x1_step1_xoffset,x1_step1_ymin); %figure,plot(xp1(2,:),x(2,:),'x')
end

for kAct=1:possActions
    
    % select only the input for one action
    inp=inpAll(:,  [1:numel(allCol(:))]  + ((kAct-1)*numel(allCol(:))));
    %     inp=[allCol(:),allGoalR(:),allGoalC(:), aoS, kAct*ones(size(allCol(:)))]';
    
    tempQ=net(inp);
    tempQ2(:,:,:,kAct)=reshape(tempQ,[mazeSize(1) mazeSize(2) mazeSize(2)]);
    
    % Calculate activity of each neuron
    if plotNeurFlag==1
        xp1=xp1All(:,[1:floor(size(xp1All,2)/possActions)] + ((kAct-1)*(floor(size(xp1All,2)/possActions))) );
        tempNeurAct(1:netSize(1),:,1)=TanSigApply(net.iw{1}*xp1+net.b{1});
        for kNeur=1:netSize(1)
            neurAct2(:,:,:,1,kNeur,kAct)=reshape(tempNeurAct(kNeur,:,1),[mazeSize(1) mazeSize(2) mazeSize(2) ]);
        end
        for kLay=2:size(netSize,2)
            tempNeurAct(1:netSize(kLay),:,kLay)=TanSigApply(net.lw{kLay,kLay-1}*tempNeurAct(1:netSize(kLay-1),:,kLay-1)+net.b{kLay});
            for kNeur=1:netSize(kLay)
                neurAct2(:,:,:,kLay,kNeur,kAct)=reshape(tempNeurAct(kNeur,:,kLay),[mazeSize(1) mazeSize(2) mazeSize(2) ]);
            end
        end
        tempQ2(:,:,:,kAct)=squeeze(neurAct2(:,:,:,currLayer,currNeur,kAct));
    end
    
end


% Not sure why this is correct but apparently % it is...
QQ=permute(tempQ2,[5 2 1 3 4]);
if plotThrFlag==1 % $$$$ NOT SURE ABOUT THIS
    QQ=permute(QQ,[1 2 3 4 5]);
end
QQ=repmat(QQ,[mazeSize(1) 1 1 1 1]);

if meanOSVFlag==0
    Q=QQ;
else
    Qall(:,:,:,:,:,kRow,kCol)=QQ;
end

    end
end

if meanOSVFlag==1
    Q=nanmean(nanmean(Qall,6),7);
end


% Plotting starts here ---------------------------------------------------

if plot1ActFlag==1
kAction=plAction;
    imagesc(squeeze(Q(handRow,handCol,:,:,kAction))); colorbar
    % caxis([-1 1])
%     title(['Q(' actNames{kAction} '), NN'])
if plotNeurFlag ==1
    title(['Activation(' actNames{kAction} ')'])
else
    title(['Q(' actNames{kAction} ')'])
end
    ylim(rowLims)
    if ~isempty(cAxVals)
        caxis(cAxVals);
    end
    drawnow
    colormap jet
else
    
    set(0,'CurrentFigure',hActVal)
% figure,
for kAction=1:possActions
    
    subplot(2,possActions*2,kAction)
    imagesc(squeeze(Qtable(handRow,handCol,:,:,kAction))); colorbar
    % caxis([-1 1])
    title(['Q(' actNames{kAction} '), table'])
    ylim(rowLims)
    drawnow
    
    subplot(2,possActions*2,possActions + kAction )
    imagesc(squeeze(Q(handRow,handCol,:,:,kAction))); colorbar
    % caxis([-1 1])
    title(['Q(' actNames{kAction} '), NN'])
    ylim(rowLims)
    drawnow
    
end

subplot(2,possActions*2,possActions*2 + 2)
pickRight=Qtable(:,:,:,:,rVal)-Qtable(:,:,:,:,lVal);
imagesc(squeeze(pickRight(handRow,handCol,:,:,1))); colorbar
caxis([-1 1])
title('TABLE: RIGHT = 1, LEFT = 0.');
colLims=max(abs(caxis));
caxis([-colLims colLims])
ylim(rowLims)
drawnow

subplot(2,possActions*2,possActions*4 - 2)
pickRight=Q(:,:,:,:,rVal)-Q(:,:,:,:,lVal);
imagesc(squeeze(pickRight(handRow,handCol,:,:,1))); colorbar
colLims=max(abs(caxis));
caxis([-colLims colLims])
ylim(rowLims)
title('NN: RIGHT +ve, LEFT -ve ');
colormap jet
drawnow

subplot(2,possActions*2,possActions*4-1)
pickRight=Q(:,:,:,:,rVal)>Q(:,:,:,:,lVal);
imagesc(squeeze(pickRight(handRow,handCol,:,:,1))); colorbar
colLims=max(abs(caxis));
caxis([-colLims colLims])
ylim(rowLims)
title('NN: RIGHT = 1, LEFT = 0.');
colormap jet
drawnow

subplot(2,possActions*2,possActions*4)
[maxValue maxAction]=max(Q(:,:,:,:,:),[],5);
imagesc(squeeze(maxAction(handRow,handCol,:,:,1))); colorbar
caxis([1 possActions])
ylim(rowLims)
title(['NN: LEFT = 1, STAY = 2, RIGHT =3 ']);
colormap jet
drawnow

subplot(2,possActions*2,possActions*2 + 1)
handPoss=zeros(size(maze2D));
handPoss(handRow,handCol)=1;
imagesc(handPoss)
caxis([-1 1])
ylim(rowLims)
title(['Hand Position, actionNum=' num2str(countActions)]);
colorbar
drawnow

end