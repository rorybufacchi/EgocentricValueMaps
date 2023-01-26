addpath('C:\Users\Rory\Documents\DPPS\DefenseAgent\Scripts')

set(0,'CurrentFigure',hActVal)

plotNeurFlag=1;

plotThrFlag=0;

currNeur=1;
currLayer=5;

handRow=size(maze2D,1)-2;
handCol=5;

clear neurAct2 tempNeurAct

% Mean osv flag
mOSVFL=0;

% other state variables - only used if meanOsvFlag is 0
% inpOSV=[handRow-6 5];
inpOSV=[11 4];


[allCol,allgoalR,allgoalC]=meshgrid(1:mazeSize(2),1:mazeSize(1),1:mazeSize(2));

% all other state variables are summarised here - $$$ TO BE DONE $$$
if size(batchEps.S{1},2)>3
    if mOSVFL==1
        % other state values - average value
        oSV=cell2mat(batchEps.S);
        oSV=nanmean(oSV(:,4:end),1);
    else
        oSV=inpOSV;
    end
    
    aoS=oSV.*ones([size(allCol(:),1) size(batchEps.S{1},2)-3]);
    
    if plotThrFlag==1 
        allPosTemp=[allgoalR(:) allgoalC(:) aoS];
        aoS(:,1)=allPosTemp(:,1);
        aoS(:,2)=allPosTemp(:,2);
        allgoalR=allPosTemp(:,3);
        allgoalC=allPosTemp(:,4);
    end
else 
    aoS=[];
end

if plotNeurFlag==1
    inpAll=[allCol(:),allgoalR(:),allgoalC(:), aoS, ones(size(allCol(:)))]';
    for kAct=2:possActions
        inpAll=[inpAll, [allCol(:),allgoalR(:),allgoalC(:),aoS, kAct*ones(size(allCol(:)))]'];
    end
    % Dimensions
    x1=inpAll;
    x1_step1_xoffset = min(x1,[],2);
    x1_step1_gain = 2./(max(x1,[],2)-min(x1,[],2)); tempinvGain=(1./max(x1,[],2));
    x1_step1_gain(x1_step1_gain==Inf)=tempinvGain(x1_step1_gain==Inf)
    x1_step1_ymin = -1;
    y1_step1_ymin = -1;
    y1_step1_gain = 2/(max(y)-min(y));
    y1_step1_xoffset = min(y);
    % Input 1
    xp1All = MapminmaxApply(x1,x1_step1_gain,x1_step1_xoffset,x1_step1_ymin); %figure,plot(xp1(2,:),x(2,:),'x')
end

for kAct=1:possActions
    inp=[allCol(:),allgoalR(:),allgoalC(:), aoS, kAct*ones(size(allCol(:)))]';
    
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
Q=permute(tempQ2,[5 2 1 3 4]);
if plotThrFlag==1 % $$$$ NOT SURE ABOUT THIS
    Q=permute(Q,[1 2 3 4 5]);
end
Q=repmat(Q,[mazeSize(1) 1 1 1 1]);

% figure,
for kAction=1:possActions
    
    subplot(2,possActions*2,kAction)
    imagesc(squeeze(Qtable(handRow,handCol,:,:,kAction))); colorbar
    % caxis([-1 1])
    title(['Q(' actNames{kAction} '), table'])
    drawnow
    
    subplot(2,possActions*2,possActions + kAction )
    imagesc(squeeze(Q(handRow,handCol,:,:,kAction))); colorbar
    % caxis([-1 1])
    title(['Q(' actNames{kAction} '), NN'])
    drawnow
    
end

subplot(2,possActions*2,possActions*2 + 2)
pickRight=Qtable(:,:,:,:,2)-Qtable(:,:,:,:,1);
imagesc(squeeze(pickRight(handRow,handCol,:,:,1))); colorbar
caxis([-1 1])
title('TABLE: RIGHT = 1, LEFT = 0.');
colLims=max(abs(caxis));
caxis([-colLims colLims])
drawnow

subplot(2,possActions*2,possActions*4 - 2)
pickRight=Q(:,:,:,:,2)-Q(:,:,:,:,1);
imagesc(squeeze(pickRight(handRow,handCol,:,:,1))); colorbar
colLims=max(abs(caxis));
caxis([-colLims colLims])
title('NN: RIGHT +ve, LEFT -ve ');
colormap jet
drawnow

subplot(2,possActions*2,possActions*4-1)
pickRight=Q(:,:,:,:,2)>Q(:,:,:,:,1);
imagesc(squeeze(pickRight(handRow,handCol,:,:,1))); colorbar
colLims=max(abs(caxis));
caxis([-colLims colLims])
title('NN: RIGHT = 1, LEFT = 0.');
colormap jet
drawnow

subplot(2,possActions*2,possActions*4)
[maxValue maxAction]=max(Q(:,:,:,:,:),[],5);
imagesc(squeeze(maxAction(handRow,handCol,:,:,1))); colorbar
caxis([1 possActions])
title(['NN: LEFT = 1, RIGHT = 2, Stay=3 ']);
colormap jet
drawnow

subplot(2,possActions*2,possActions*2 + 1)
handPoss=zeros(size(maze2D));
handPoss(handRow,handCol)=1;
imagesc(handPoss)
caxis([-1 1])
title(['Hand Position, actionNum=' num2str(countActions)]);
colorbar
drawnow