
clear all

fdir='/Users/fengyanshi/TMP/tmp2/';
fdir_input='../grid_new/';
X1=load([fdir_input 'lon_curv.txt']);
Y1=load([fdir_input 'lat_curv.txt']);
dep=load([fdir_input 'dep_curv.txt']);

X=X1(1:end-1,1:end-1);
Y=Y1(1:end-1,1:end-1);

lon_m=min(min(X));
lat_m=min(min(Y));

R=6371000.0;

x=(X-lon_m)*R*pi/180.0*cos(pi/180.0*36.96);
y=(Y-lat_m)*R*pi/180.0;

fac=0.002;
%xyaxis=[-76.336 -76.25 36.94 36.98];
xyaxis=[-76.305 -76.27 36.94 36.97]; % small
xyaxis1=[3200 6000 1000 4000]; % x y
skx=2;
sky=2;

Dep=dep(1:end-1,1:end-1);

m1=188;
n1=60;
m2=331;
n2=35;
m3=338;
n3=70;
m4=194;
n4=109;

Xse=X(n1,m1);
Yse=Y(n1,m1);
Xsw=X(n2,m2);
Ysw=Y(n2,m2);

Xne=X(n3,m3);
Yne=Y(n3,m3);
Xnw=X(n4,m4);
Ynw=Y(n4,m4);




figure(1)
pcolor(X,Y,Dep),shading interp
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' Depth (m) ');
%demcmap(B);
hold on
%pcolor(X,Y,cc),shading interp
%colormap(jet)
axis(xyaxis)
caxis([-4 4])
plot([X(n1,m1) X(n2,m2) X(n3,m3) X(n4,m4) X(n1,m1)],[Y(n1,m1) Y(n2,m2) Y(n3,m3) Y(n4,m4) Y(n1,m1)],'wo','LineWidth',2)
Xp=[num2str(X(n1,m1)) ',  ' num2str(Y(n1,m1))];
text(X(n1,m1),Y(n1,m1),Xp,'Color','w')
Xp=[num2str(X(n2,m2)) ',  ' num2str(Y(n2,m2))];
text(X(n2,m2),Y(n2,m2),Xp,'Color','k')
Xp=[num2str(X(n3,m3)) ',  ' num2str(Y(n3,m3))];
text(X(n3,m3),Y(n3,m3),Xp,'Color','k')
Xp=[num2str(X(n4,m4)) ',  ' num2str(Y(n4,m4))];
text(X(n4,m4),Y(n4,m4),Xp,'Color','w')

xlabel('Lon(deg) ','fontsize',12,'fontweight','bold');
ylabel('Lat(deg) ','fontsize',12,'fontweight','bold');


% mod
%for kx=1:7
%xp=
%
%end

%axis(xyaxis);
print -djpeg100 harbor_ll.jpg



figure(2)
pcolor(x,y,Dep),shading interp
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' Depth (m) ');
%demcmap(B);
hold on
plot([x(n1,m1) x(n2,m2) x(n3,m3) x(n4,m4) x(n1,m1)],[y(n1,m1) y(n2,m2) y(n3,m3) y(n4,m4) y(n1,m1)],'wo','LineWidth',2)
%pcolor(X,Y,cc),shading interp
%colormap(jet)
axis(xyaxis1)
caxis([-4 4])
xlabel('East(m) ','fontsize',12,'fontweight','bold');
ylabel('North(m) ','fontsize',12,'fontweight','bold');
print -djpeg100 harbor_xy.jpg



