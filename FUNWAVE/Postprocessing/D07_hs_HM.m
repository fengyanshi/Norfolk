
clear all
close all

froot='/Volumes/BigSur_2022/Norfolk_FUNWAVE/Results/';
%fcase={'R1_smdm','R2_smdm','R3_smdm','R4_smdm','R5_smdm','R6_smdm'};
case_name={'ERA5','HM','Acc GN 1m','RMW F1.25','SLR 1.3m','WSF 1.225'};
fcase={'R2_smdm'};
case_name={'HM'};


files=[25];
dt=100.0;

ax=[-76.299 -76.246 36.94 36.98];
cax=[0 2.5];
Mindep=0.1;
uv_deplimit = 0.2;

% read info

m=7744;
n=5888;

DimsX={[m n]};
% dimensions
% load x and y 
load([froot 'lonlat.mat']);

x0=x(1,1);
y0=y(1,1);
%dx=0.75;
%dy=0.75;
%x=[0:m-1]*dx;
%y=[0:n-1]*dy;
sk=2;


fig=figure(1);
colormap jet
clf

wid=10;
len=8;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);


numb=files(1);

fnum=sprintf('%.5d',numb);

for k=1:length(fcase)
% read files -----------------------

fdir=[froot fcase{k} '/'];

fname=[fdir 'dep.out'];
fileID=fopen(fname);
dep=fread(fileID,DimsX{1},'*single');
fclose(fileID);
dep=dep';

fname=[fdir 'Hsig_' fnum];
fileID=fopen(fname);
hs=fread(fileID,DimsX{1},'*single');
fclose(fileID);
hs=hs';


hs(dep<0.0)=NaN;


plot([x0 x0+0.03],[y0 y0+0.04],'.r','MarkerSize',1)
plot_google_map('maptype','satellite','APIKey','AIzaSyBeu2oRBtLClpcm4i2VDIXltuzMAOY5yX4')
hold on


b=pcolor(x(1:sk:end,1:sk:end),y(1:sk:end,1:sk:end),hs(1:sk:end,1:sk:end)),shading flat
set(b,'FaceAlpha',0.8);
set(b,'AlphaData',~isnan(hs))
time1=[case_name{k}];
title(time1)
caxis(cax)
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' Hs (m) ');

xlabel('Lon (deg)')
ylabel('Lat (deg)')

daspect([1 1 0.8])
%axis equal
axis(ax)

end

eval(['mkdir ' 'plots/snap'])

pname=['plots/snap/' case_name{1} '_' 'hs.jpg'];
print('-djpeg100',pname);



