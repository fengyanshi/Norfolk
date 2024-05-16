clear all
dep=load('model_1p5m_7871x6098.txt');
dep1=dep(1:3000,1472:end);
[n m]=size(dep1);
dx=1.5;
dy=1.5;
x=[0:m-1]*dx;
y=[0:n-1]*dy;

dep1(dep1<-5)=-5;

fig=figure(1)
clf
wid=10;
len=5;
colormap jet
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
pcolor(x,y,-dep1),shading flat
xlabel('x (m)')
ylabel('y (m)')
demcmap(-dep1)
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' bathymetry (m) ');

print -djpeg100 funwave_grid_2d.jpg

dep1(dep1>8.0)=8.0;

save -ASCII dep_1p5m_6400x3000.txt dep1