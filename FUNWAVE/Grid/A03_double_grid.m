clear all
fdir='/Users/fengyanshi/OUTSIDE_Google/Users/FUNWAVE_NORFOLK/Chris/Chris_Norfolk/';
dep=load([fdir 'model_1p5m_7871x6098.txt']);
Zold=dep(1:3000,1472:end);

Znew=double_grid(Zold);

[n m]=size(Znew);

dx=0.75;
dy=0.75;
x=[0:m-1]*dx;
y=[0:n-1]*dy;

Znew(Znew<-5)=-5;

fig=figure(1)
clf
wid=10;
len=5;
colormap jet
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
pcolor(x(1:4:end,1:4:end),y(1:4:end,1:4:end),-Znew(1:4:end,1:4:end)),shading flat
xlabel('x (m)')
ylabel('y (m)')
demcmap(-Znew)
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' bathymetry (m) ');

print -djpeg100 funwave_grid_double.jpg

dep=Znew;

% filtering
edge=5; % set edge points yo dont want filter

K=0.11111*ones(3);
tmp = conv2(dep,K,'same');

dep(1+edge:end-edge,1+edge:end-edge)=tmp(1+edge:end-edge,1+edge:end-edge);

% heavy filtering at the building
%m1=9200;
%m2=9900;
%n1=1000;
%n2=1450;

m1=8500;
m2=10500;
n1=500;
n2=1800;

clear tmp
tmp=dep(n1:n2,m1:m2);

for kk=1:10
tmp1=tmp;
K=0.11111*ones(3);
% 5 points
% K=ones(5)/25.;
tmp2 = conv2(tmp1,K,'same');

tmp(1+edge:end-edge,1+edge:end-edge)=tmp2(1+edge:end-edge,1+edge:end-edge);

end

dep(n1+wid:n2-wid,m1+wid:m2-wid)=tmp(1+wid:end-wid,1+wid:end-wid);

Znew=dep;
Znew(Znew>8.0)=8.0;

dep1=Znew(1:5888,1:12544);


save('-ASCII', [fdir 'dep_075m_12544x5888_filter.txt'], 'dep1');