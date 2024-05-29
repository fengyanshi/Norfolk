
clear all
close all

froot='/Volumes/BigSur_2022/Norfolk_nearcom/Results/';
fcase={'D_IR_Base_ERA5','D_IR_Base_HM','D_IR_Bathy_Acc_GN_1m','D_IR_RMW_F1.25','D_IR_SLR_1.3M','D_IR_WSF_1.225'};
%fcase={'D_IR_Base_ERA5'};
case_name={'ERA5','HM','Acc GN 1m','RMW F1.25','SLR 1.3m','WSF 1.225'};

fdir_input='../Grid_double/';

files=[91 98 98 97 98 91];
DT=600.0;

ax=[-76.3 -76.246 36.94 36.975];

cax_hs=[0 2.5];
cax_eta=[0 3];

Mindep=[0.1 0.4 0.4 0.4 0.4 0.2];
uv_deplimit = 0.2;

% read info

m=768;
n=304;

X1=load([fdir_input 'lon_curv.txt']);
Y1=load([fdir_input 'lat_curv.txt']);

Xx=load([fdir_input 'x_circ.txt']);
Yy=load([fdir_input 'y_circ.txt']);

dxx=diff(Xx,1,2);
dyx=diff(Yy,1,2);

dxy=diff(Xx);
dyy=diff(Yy);


dx=sqrt(dxx.^2+dyx.^2);
dy=sqrt(dxy.^2+dyy.^2);

m1=768;
n1=303;
dx1=dx(1:n1,1:m1);
dy1=dy(1:n1,1:m1);

area=dx1.*dy1;

X=X1(1:n,1:m);
Y=Y1(1:n,1:m);

fac=0.004;
xyaxis=[-76.31 -76.255 36.94 36.975];
skx=4;
sky=4;
skxW=8;
skyW=8;

fig=figure(1);
colormap jet
clf

wid=10;
len=10;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);


for k=1:length(fcase)
%for k=5:5
% read files -----------------------

if k==5 
dep=load([fdir_input 'dep_circ_SLR.txt']);
else
dep=load([fdir_input 'dep_circ.txt']);
end

fdir=[froot fcase{k} '/'];


if k==5 
fdir0=[froot fcase{1} '/'];
eta0=load([fdir0 'eta_0001']);
mask0=load([fdir0 'mask_0001']);
else
eta0=load([fdir 'eta_0001']);
mask0=load([fdir 'mask_0001']);
end

[n,m]=size(eta0);

Dep=dep(1:n,1:m);
if k==5
MASK0=mask0;
else
MASK0=ones(size(mask0));
MASK0(eta0+Dep<Mindep(k))=0;
end

%eta(MASK<1)=NaN;


fname=sprintf('%.4d',files(k));
eta=load([fdir 'eta_' fname]);
mask=load([fdir 'mask_' fname]);

MASK=ones(size(mask));
MASK(eta+Dep<Mindep(k))=0;


floodmask=MASK-MASK0;
floodmask_c=floodmask;

floodmask(floodmask<1)=NaN;

% ---------------
subplot(3,2,k)

x0=X(1,1);
y0=Y(1,1);

plot([x0 x0+0.07],[y0 y0+0.09],'.r','MarkerSize',1)
plot_google_map('maptype','satellite','APIKey','AIzaSyBeu2oRBtLClpcm4i2VDIXltuzMAOY5yX4')
hold on

b=pcolor(X,Y,floodmask),shading flat
set(b,'FaceAlpha',0.8);
set(b,'AlphaData',~isnan(floodmask))

caxis(cax_eta)
%cbar=colorbar;
%set(get(cbar,'ylabel'),'String',' \eta (m) ');

axis(xyaxis);
title([case_name{k}])


if(k==5 || k==6)
xlabel('Lon (deg)')
end
if(k==1 || k==3 ||k==5)
ylabel('Lat (deg)')
end

daspect([1 1 0.8])
%axis equal
axis(xyaxis)

% calculate area
floodmask1=floodmask_c(1:n1,1:m1);
floodarea1=floodmask1.*area;
floodarea(k)=sum(sum(floodarea1));


end

eval(['mkdir ' 'plots/six_cases'])

pname=['plots/six_cases' '/' 'flood.jpg'];
print('-djpeg100',pname);

save('-ASCII','plots/six_cases/flood_area.txt','floodarea');


