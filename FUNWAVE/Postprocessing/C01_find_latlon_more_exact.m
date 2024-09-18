
clear all
close all

froot='/Volumes/BigSur_2022/Norfolk_FUNWAVE/Results/';


% linear coordinate transformation ------

x2=[-76.291282 -76.270701];
y2=[36.947564 36.963229];

x1=[4984 7896];
y1=[830 2516];

a=diff(x2./y1)/diff(x1./y1);
b=diff(x2./x1)/diff(y1./x1);
c=diff(y2./y1)/diff(x1./y1);
d=diff(y2./x1)/diff(y1./x1);

% ----------------------------------------

ax=[-76.31 -76.24 36.94 36.98];
cax=[0 2.5];
Mindep=0.1;


m=7744;
n=5888;

%X1d=[1:m];
%Y1d=[1:n];

X1=100 ;
Y1=100;

%[X1,Y1]=meshgrid(X1d,Y1d);

X2=a*X1+b*Y1;
Y2=c*X1+d*Y1;


sk=4;



%save([froot 'lonlat.mat'], 'x', 'y');

fig=figure(1);
colormap parula
clf

wid=10;
len=10;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);


dep=load('/Volumes/BigSur_2022/Norfolk_FUNWAVE/dep_075m_12544x5888.txt');

DEP=dep(1:2:end,1:2:end);

clf
pcolor(-DEP),shading flat
caxis([-0.1 0.1])

