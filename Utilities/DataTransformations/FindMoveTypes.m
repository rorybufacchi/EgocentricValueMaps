function [moveType squishQs QbyType Q] = FindMoveTypes(net,s,w,aS)
%[moveType squishQs QbyType] = FindMoveTypes(net,s,w,aS)
%   aS is alteration settings. If empty, it just runs on the normal network

if s.plt.meanLimbCols ~= 1
    warning('Not averaging over limb cols when culculating best movement types')
    warning('setting it to 1 now')
    s.plt.meanLimbCols = 1;
end

if exist('aS')
    ablNet = AlterNetwork(net,aS);
else
    ablNet = net;
end
s.plt.OneActFl=0;
[Q,allNeurAct] = CalcNetOutput(s,w,ablNet);
if s.plt.ON==1, figure; end
[bestQ squishQs] = DisplActValsFun(s,w,Q);


% Define type of movement as 1 for defensive and -1 for goal oriented
% figure,plot(bestQ')

% This is -1 for left of center and +1 for right of center
nC=size(bestQ,2);
Pos=ones(size(bestQ));
Pos(:,floor(nC/2)+1:end)=-1;
% if mod(nC,2)==1, Pos(:,ceil(nC/2))=0; end

% Defensive moves are -1, and appettitive moves are +1
moveType=nan(size(bestQ));
moveType(bestQ>0)=1;
moveType(bestQ<0)=-1;
moveType(bestQ==0)=0;

% If the stimulus is directly above the limb, things are different
if mod(nC,2)==1
    moveType(abs(moveType(:,ceil(nC/2)))==1,ceil(nC/2))=-1; % any movement is defensive
    moveType(abs(moveType(:,ceil(nC/2)))==0,ceil(nC/2))=1; % and staying is appetitive
end

moveType=moveType.*Pos;

if s.plt.ON==1
    figure, imagesc(moveType); axis xy; view(180,90);
end

% Define Appetitive Qs
grabQs=squishQs(:,:,1);
grabQs(:,floor(nC/2)+1:end)=squishQs(:,floor(nC/2)+1:end,3);
% Define Defensive Qs
dodgeQs=squishQs(:,:,3);
dodgeQs(:,floor(nC/2)+1:end)=squishQs(:,floor(nC/2)+1:end,1);
% If the stimulus is directly above the limb, things are again different
if mod(nC,2)==1
    grabQs(:,ceil(nC/2))=squishQs(:,ceil(nC/2),2);
    dodgeQs(:,ceil(nC/2))=max(squishQs(:,ceil(nC/2),[1 3]),[],3);
end

QbyType(:,:,1)=grabQs;
QbyType(:,:,2)=dodgeQs;


end

