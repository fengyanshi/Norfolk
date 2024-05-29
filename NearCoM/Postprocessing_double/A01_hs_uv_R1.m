
clear all
close all

froot='/Volumes/BigSur_2022/Norfolk_nearcom/Results/';
fcase='D_IR_Base_ERA5';
casename='ERA5';

fdir=[froot fcase '/'];
fdir_input='../Grid_double/';
fdir_data=['../Data_double/' fcase '/'];

files=[91];
DT=600.0;

Mindep=0.1;
uv_deplimit = 0.2;

% read info

filename=[fdir_data 'info4plot.txt'];
T=readtable(filename);
dt=datetime(T.Var1);
daynum=datenum(dt);


m=768;
n=304;

X1=load([fdir_input 'lon_curv.txt']);
Y1=load([fdir_input 'lat_curv.txt']);
dep=load([fdir_input 'dep_circ.txt']);

X=X1(1:n,1:m);
Y=Y1(1:n,1:m);

fac=0.004;
%xyaxis=[-76.333 -76.26 36.94 36.975];
xyaxis=[-76.31 -76.255 36.94 36.975];
%xyaxis=[-76.305 -76.265 36.96 36.974]; % small
skx=4;
sky=4;
skxW=8;
skyW=8;



%set(0,'DefaultFigureColormap',feval('jet'));
figure(1)
wid=10;
len=6;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[2 2 wid len],'paperposition',[0 0 wid len]);

colormap parula

icount=0;
for num=1:length(files)
icount=icount+1;

fname=sprintf('%.4d',files(num));
eta=load([fdir 'eta_' fname]);
mask=load([fdir 'mask_' fname]);
u=load([fdir 'u_' fname]);
v=load([fdir 'v_' fname]);
hs=load([fdir 'hs_' fname]);
wdir=load([fdir 'wdir_' fname]);

wx=cos(wdir*pi/180)*0.25;
wy=sin(wdir*pi/180)*0.25;

uu=sqrt(u.^2+v.^2);


[ux,uy]=gradient(u,1,1);
[vx,vy]=gradient(v,1,1);
ww=(uy-vx);

[n,m]=size(eta);

Dep=dep(1:n,1:m);
MASK=ones(size(mask));
MASK(eta+Dep<Mindep)=0;
hs(MASK<1)=NaN;
eta(MASK<1)=NaN;
u(MASK<1)=NaN;
v(MASK<1)=NaN;
ww(MASK<1)=NaN;
wx(MASK<1)=NaN;
wy(MASK<1)=NaN;

u(Dep<uv_deplimit)=NaN;
v(Dep<uv_deplimit)=NaN;

Dtime=daynum+(files(num)-1)*DT/3600/24;
Dtime_datetime=datetime(Dtime,'ConvertFrom','datenum');
titledate=datestr(Dtime_datetime);

% elevation scale
%maxele=max(max(eta));
%minele=maxele-0.1;
cax_hs=[0 2.5];
cax_eta=[-0.5 2.5];

maxele=1.70;
minele=-0.5;
maxhs=2.0;
minhs=0.0;

figure(1)
clf

x0=X(1,1);
y0=Y(1,1);

plot([x0 x0+0.07],[y0 y0+0.09],'.r','MarkerSize',1)
plot_google_map('maptype','satellite','APIKey','AIzaSyBeu2oRBtLClpcm4i2VDIXltuzMAOY5yX4')
hold on

b=pcolor(X,Y,hs),shading interp
set(b,'FaceAlpha',0.8);
set(b,'AlphaData',~isnan(hs))

cbar=colorbar;
caxis(cax_hs)
set(get(cbar,'ylabel'),'String',' Hs (m) ');

hold on

quiver(X(1:skyW:end,1:skxW:end),Y(1:skyW:end,1:skxW:end),wx(1:skyW:end,1:skxW:end)*fac,wy(1:skyW:end,1:skxW:end)*fac,0,'Color','k');

xlabel('Lon(deg) ','fontsize',12,'fontweight','bold');
ylabel('Lat(deg) ','fontsize',12,'fontweight','bold');
axis(xyaxis);
title(['CASE: ' casename  ' ( ' titledate ')'])
caxis([minhs maxhs])

eval(['mkdir ' 'plots/report/']);
fname=['plots/report/hs_' fcase '_' num2str(files(num)) '.jpg'];
print('-djpeg100', fname)

figure(2)
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[2 2 wid len],'paperposition',[0 0 wid len]);
clf

plot([x0 x0+0.07],[y0 y0+0.09],'.r','MarkerSize',1)
plot_google_map('maptype','satellite','APIKey','AIzaSyBeu2oRBtLClpcm4i2VDIXltuzMAOY5yX4')
hold on

b=pcolor(X,Y,eta),shading interp
set(b,'FaceAlpha',0.8);
set(b,'AlphaData',~isnan(hs))

caxis(cax_eta)
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' \eta (m) ');
hold on
quiver(X(1:sky:end,1:skx:end),Y(1:sky:end,1:skx:end),u(1:sky:end,1:skx:end)*fac,v(1:sky:end,1:skx:end)*fac,0,'Color','k');

xlabel('Lon(deg) ','fontsize',12,'fontweight','bold');
ylabel('Lat(deg) ','fontsize',12,'fontweight','bold');
axis(xyaxis);
title(['CASE: ' casename  ' ( ' titledate ')'])

eval(['mkdir ' 'plots/report/']);
fname=['plots/report/ele_uv_' fcase '_' num2str(files(num)) '.jpg'];
print('-djpeg100', fname)
end



