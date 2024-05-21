
clear all
close all

froot='/Volumes/BigSur_2022/Norfolk_FUNWAVE/Results/';
fcase='R2';

fdir=[froot fcase '/'];
fdir_input='../Grid/';
fdir_data=['../Data/' fcase '/'];

files=[18];
dt=100.0;

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

fname=[fdir 'dep.out'];
fileID=fopen(fname);
dep=fread(fileID,DimsX{1},'*single');
fclose(fileID);
dep=dep';


fig=figure(1);
colormap jet

wid=6;
len=12;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

for k=1:length(files) 

numb=files(k);

fnum=sprintf('%.5d',numb);

% read files -----------------------

fname=[fdir 'etamean_' fnum];
fileID=fopen(fname);
etamean=fread(fileID,DimsX{1},'*single');
fclose(fileID);
etamean=etamean';

fname=[fdir 'umean_' fnum];
fileID=fopen(fname);
umean=fread(fileID,DimsX{1},'*single');
fclose(fileID);
umean=umean';

fname=[fdir 'vmean_' fnum];
fileID=fopen(fname);
vmean=fread(fileID,DimsX{1},'*single');
fclose(fileID);
vmean=vmean';

[vort vort1]=curl(x,y,umean,vmean);

fname=[fdir 'Hsig_' fnum];
fileID=fopen(fname);
hs=fread(fileID,DimsX{1},'*single');
fclose(fileID);
hs=hs';


etamean(dep<0.0)=NaN;
hs(dep<0.0)=NaN;


clf

subplot(311)
pcolor(x,y,etamean),shading flat
time1=['mean surface at ' num2str(1000+numb*dt) ' s'];
title(time1)
caxis([-0.2 0.2 ])

xlabel('x (m)')
ylabel('y (m)')
axis equal

subplot(312)
pcolor(x,y,hs),shading flat
time1=['Hsig at ' num2str(1000+numb*dt) ' s'];
title(time1)
caxis([0 2 ])

xlabel('x (m)')
ylabel('y (m)')
axis equal

subplot(313)
pcolor(x,y,vort),shading flat
time1=['Vorticity at ' num2str(1000+numb*dt) ' 1/s'];
title(time1)
caxis([-0.1 0.1 ])

xlabel('x (m)')
ylabel('y (m)')
axis equal

eval(['mkdir ' 'plots/' fcase])

pname=['plots/' fcase '/' 'hs_' num2str(numb) '.jpg'];
print('-djpeg100',pname);

end



