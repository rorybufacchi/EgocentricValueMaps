
%%
rng=5;
testx=rand(8,8,8)>0.85;
%change this to the 3d Matrix you want to display:
input=testx;
alpha = 1 ; % transparency
%------------------------
disp('calculating...');

%cropping the matrix ( reduces lag for large empty matrices)
[yy,xx]=find(sum(input,3));
xmin = min(xx(:));
xmax = max(xx(:));
ymin = min(yy(:));
ymax = max(yy(:));
[xx,zz] = find(squeeze(sum(input,1)));
zmin=min(zz(:));
zmax = max(zz(:));
testcrop=input(ymin:ymax,xmin:xmax,zmin:zmax);


%uncomment this to shrink the volume of the clusters
%testcrop=meltVolume(testcrop,18,0.7);

%color different clusters differently:
testcrop=colorVoxelGroups(double(testcrop));

figure(2);
voxelSurf(testcrop,true,[0.1 0.8 0.1 0.8 0.1 0.8],alpha );

%use this if you want to have a separate color for each voxelgroup:
colormap(lines(max(testcrop(:))));

%colormap(lines(1));

disp('done');
