
clear all
close all

froot='/Volumes/BigSur_2022/Norfolk_FUNWAVE/Results/';
fcase='R1_smdm';

fdir=[froot fcase '/'];
fdir_input='../Grid/';
fdir_data=['../Data/' fcase '/'];

files=[28];
dt=10.0;

Mindep=0.1;
uv_deplimit = 0.2;

% read info



m=7744;
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


% define movie file and parameters
myVideo = VideoWriter(['plots/eta_' fcase '.mp4'],'MPEG-4');
myVideo.FrameRate = 10;  
myVideo.Quality = 100;
%vidHeight = 576; %this is the value in which it should reproduce
%vidWidth = 1024; %this is the value in which it should reproduce
open(myVideo);

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

pause(0.1)
% save image
F = print('-RGBImage','-r300');
%J = imresize(F,[vidHeight vidWidth]);
mov(k).cdata = F;

writeVideo(myVideo,mov(k).cdata);

end
close(myVideo)


