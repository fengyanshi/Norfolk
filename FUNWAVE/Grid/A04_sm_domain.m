clear all
fdir='/Users/fengyanshi/OUTSIDE_Google/Users/FUNWAVE_NORFOLK/Chris/Chris_Norfolk/';
dep=load([fdir 'dep_075m_12544x5888_filter.txt']);
Znew=dep(1:5888,3800:end-1000);

[n m]=size(Znew);
dx=0.75;
dy=0.75;
x=[0:m-1]*dx;
y=[0:n-1]*dy;

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

print -djpeg100 funwave_grid_sm.jpg

Znew(Znew>6.0)=6.0;
Znew(:,7300:end)=6.0;
dep=Znew(:,1:7744);

save('-ASCII', [fdir 'dep_075m_7744x5888_filter.txt'], 'dep');