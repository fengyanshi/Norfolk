clear all
dep=load('model_1p5m_7871x6098.txt');
dep1=dep(1:3000,1500:end);
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

% profs
x0_1=[4000 4400 4800 5200 5600 6000];
y0_1=[2800 2660 2570 2450 2200 1850];
angle1=[60.0 62.5 62.5 50.0 45 40]*pi/180.;
dx1=1.0;
m1=1500;

hold on
for np=1:length(x0_1)
x1(np,:)=x0_1(np)+dx1*[0:m1-1]*cos(angle1(np));
y1(np,:)=y0_1(np)+dx1*[0:m1-1]*sin(angle1(np));
plot(x1(np,:),y1(np,:),'k-','LineWidth',1);
name1=['P_' num2str(np)];
text(x1(np,end),y1(np,end),name1)
end
print -djpeg100 funwave_grid.jpg

dep_1=griddata(x,y,dep1,x1,y1);

fig=figure(2)
clf
wid=8;
len=10;
colormap jet
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

for np=1:length(x0_1)
subplot(length(x0_1),1,np)
plot(-dep_1(np,:),'LineWidth',1)
grid
ylabel('z (m)')
axis([0 1500 -7 5])
hold on
plot([0 1500],[0 0],'k-')
name1=['P_' num2str(np)];
text(1400,3,name1)
end
xlabel('x1 (m)')

print -djpeg100 funwave_prof.jpg



