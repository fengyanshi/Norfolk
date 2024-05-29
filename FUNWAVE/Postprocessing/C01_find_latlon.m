
clear all
close all

froot='/Volumes/BigSur_2022/Norfolk_FUNWAVE/Results/';
fcase={'R1_smdm','R2_smdm','R3_smdm','R4_smdm','R5_smdm','R6_smdm'};
case_name={'ERA5','HM','Acc GN 1m','RMW F1.25','SLR 1.3m','WSF 1.225'};

fcase={'R1_smdm'};
case_name={'ERA5','HM','Acc GN 1m','RMW F1.25','SLR 1.3m','WSF 1.225'};


%fdir_input='../Grid/';
%fdir_data=['../Data/' fcase_whole{1} '/'];

files=[25];
dt=100.0;

ax=[-76.31 -76.24 36.94 36.98];
cax=[0 2.5];
Mindep=0.1;
uv_deplimit = 0.2;

% read info


m=7744;
n=5888;

DimsX={[m n]};
% dimensions


x0=-76.305+0.00525;
y0=36.9425-0.0025;
angle=11.0*pi/180.
dx1=0.75;
dy1=0.75;
R=6371000.0;
Sx=R*pi/180.0*cos(pi/180.0*36.96);
Sy=R*pi/180.0;
rat=Sx/Sy;

x1=x0+[0:m-1]*dx1/Sx*1.00;
y1=y0+[0:n-1]*dy1/Sy*0.975;
sk=4;

[x2 y2]=meshgrid(x1,y1);
x=x0+(x2-x0).*cos(angle)-(y2-y0).*sin(angle);
y=y0+(x2-x0).*sin(angle)+(y2-y0).*cos(angle);

save([froot 'lonlat.mat'], 'x', 'y');

fig=figure(1);
colormap parula
clf

wid=10;
len=10;
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


%subplot(3,2,k)


plot([x0 x0+0.07],[y0 y0+0.09],'.r','MarkerSize',1)
plot_google_map('maptype','satellite','APIKey','AIzaSyBeu2oRBtLClpcm4i2VDIXltuzMAOY5yX4')
hold on

b=pcolor(x(1:sk:end,1:sk:end),y(1:sk:end,1:sk:end),hs(1:sk:end,1:sk:end)),shading flat
set(b,'FaceAlpha',0.6);
set(b,'AlphaData',~isnan(hs))

time1=[case_name{k}];
title(time1)
caxis(cax)
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' Hs (m) ');

if(k==5 || k==6)
xlabel('x (m, model coordinate)')
end
if(k==1 || k==3 ||k==5)
ylabel('y (m)')
end
%axis equal
daspect([1 1 rat])
axis(ax)

end

eval(['mkdir ' 'plots/six_cases'])

pname=['plots/six_cases' '/' 'hs.jpg'];
%print('-djpeg100',pname);




