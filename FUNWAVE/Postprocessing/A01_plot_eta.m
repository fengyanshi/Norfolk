
clear all
close all

froot='/Volumes/BigSur_2022/Norfolk_FUNWAVE/Results/';
fcase='R1';

fdir=[froot fcase '/'];
fdir_input='../Grid/';
fdir_data=['../Data/' fcase '/'];

files=[25];
dt=10.0;

Mindep=0.1;
uv_deplimit = 0.2;

% read info



m=12544;
n=5888;

DimsX={[m n]};
% dimensions
dx=0.75;
dy=0.75;
x=[0:m-1]*dx;
y=[0:n-1]*dy;
sk=4;

fname=[fdir 'dep.out'];
fileID=fopen(fname);
dep=fread(fileID,DimsX{1},'*single');
fclose(fileID);
dep=dep';

fig=figure(1);
colormap jet

wid=12;
len=8;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

for k=1:length(files) 

numb=files(k);

fnum=sprintf('%.5d',numb);

% read files -----------------------

fname=[fdir 'eta_' fnum];
fileID=fopen(fname);
eta=fread(fileID,DimsX{1},'*single');
fclose(fileID);
eta=eta';


fname=[fdir 'mask_' fnum];
fileID=fopen(fname);
mask=fread(fileID,DimsX{1},'*single');
fclose(fileID);
mask=mask';


eta(mask<1)=NaN;

clf

pcolor(x(1:sk:end),y(1:sk:end),eta(1:sk:end,1:sk:end)),shading flat
time1=['surface at ' num2str(numb*dt) ' s'];
title(time1)
caxis([-1.5 4])
xlabel('x (m)')
ylabel('y (m)')
axis equal

eval(['mkdir ' 'plots/' fcase])

pname=['plots/' fcase '/' 'eta_' num2str(numb) '.jpg'];
print('-djpeg100',pname);

end



