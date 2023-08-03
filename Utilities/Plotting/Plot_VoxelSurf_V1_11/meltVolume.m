function out = meltVolume(data,friction,drain)
%"melts"solid described in voxel clouds in data. This script has to be
%called iteratively until the desired amount of melting is archieved.
%The values of data have to be between 0 and 1.
%min friction is sum of mf =17.8065
%1 means solid
%0 means fully melted
%add option to keep single pixels!!
friction=max(friction,17.8065);
 nb=[-1  1  0  0  0  0 -1 -1 -1 -1  0  0  0  0  1  1  1  1  -1 -1 -1 -1 +1  1  1  1;...
      0  0 -1  1  0  0 -1 +1  0  0 -1 -1  1  1 -1  1  0  0  -1 -1 +1 +1 -1 -1  1  1;...
      0  0  0  0 -1  1  0  0 -1 +1 -1  1 -1  1  0  0 -1  1  -1 +1 -1 +1 -1  1 -1  1];
mf=zeros(1,26);
mf(1:6)=1;
mf(7:18)=1/sqrt(2);
mf(9:26)=1/sqrt(3);
     [dataX,dataY,dataZ]=size(data);
b=zeros(27,dataX,dataY,dataZ);
for xx=1:26
    dx=nb(1,xx);
    dy=nb(2,xx);
    dz=nb(3,xx);
    temp=data(2+dx:end-1+dx,2+dy:end-1+dy,2+dz:end-1+dz);
    b(xx,2:end-1,2:end-1,2:end-1)=mf(xx)*(squeeze(data(2:end-1,2:end-1,2:end-1))-temp);
end
out=min(data,data-squeeze(sum(b,1))./friction);
out(out<drain)=0;

% survivors=(data>=drain);
% b2=zeros(6,dataX,dataY,dataZ);
% for xx=1:6
%     dx=nb(1,xx);
%     dy=nb(2,xx);
%     dz=nb(3,xx);
%     b2(xx,2:end-1,2:end-1,2:end-1)=survivors(2+dx:end-1+dx,2+dy:end-1+dy,2+dz:end-1+dz);
% end
% survivors=survivors.*(1-squeeze(b2(1,:,:,:))).*...
%     (1-squeeze(b2(2,:,:,:))).*...
%     (1-squeeze(b2(3,:,:,:))).*...
%     (1-squeeze(b2(4,:,:,:))).*...
%     (1-squeeze(b2(5,:,:,:))).*...
%     (1-squeeze(b2(6,:,:,:)));
% 
% out=out+drain*survivors;