
clear all
close all

froot='/Volumes/BigSur_2022/Norfolk_FUNWAVE/Results/';
fcase={'R1_smdm','R2_smdm','R3_smdm','R4_smdm','R5_smdm','R6_smdm'};
case_name={'ERA5','HM','Acc GN 1m','RMW F1.25','SLR 1.3m','WSF 1.225'};
%fcase={'R1_smdm'};
%case_name={'ERA5'};


files=[25];
dt=100.0;

%ax=[-76.3 -76.246 36.94 36.975];
ax=[-76.299 -76.275 36.96 36.972];
cax=[-0.1 0.1];
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
dx=0.75;
dy=0.75;
xc=[0:m-1]*dx;
yc=[0:n-1]*dy;
sk=4;
skx=128;
sky=128;
fac=0.005;
mindepuv=0.5;


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

[vort vort1]=curl(xc,yc,umean,vmean);

vort(dep<mindepuv)=NaN;
umean(dep<mindepuv)=NaN;
vmean(dep<mindepuv)=NaN;

subplot(3,2,k)


plot([x0 x0+0.03],[y0 y0+0.04],'.r','MarkerSize',1)
plot_google_map('maptype','satellite','APIKey','AIzaSyBeu2oRBtLClpcm4i2VDIXltuzMAOY5yX4')
hold on


b=pcolor(x(1:sk:end,1:sk:end),y(1:sk:end,1:sk:end),vort(1:sk:end,1:sk:end)),shading flat
set(b,'FaceAlpha',0.8);
set(b,'AlphaData',~isnan(vort))

time1=[case_name{k}];
title(time1)
%hold on
%quiver(x(1:skx:end,1:sky:end),y(1:sky:end,1:sky:end),umean(1:sky:end,1:skx:end)*fac,vmean(1:sky:end,1:skx:end)*fac,0,'Color','k');
caxis(cax)

cbar=colorbar;
set(get(cbar,'ylabel'),'String',' Vorticity (1/s) ');

if(k==5 || k==6)
xlabel('Lon (deg)')
end
if(k==1 || k==3 ||k==5)
ylabel('Lat (deg)')
end

daspect([1 1 0.8])
axis(ax)

end

eval(['mkdir ' 'plots/six_cases'])

pname=['plots/six_cases' '/' 'vort.jpg'];
print('-djpeg100',pname);




