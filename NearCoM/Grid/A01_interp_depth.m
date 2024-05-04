clear all
fdir='/Volumes/DISK_2020_5/Norfolk_DEM/model/';
corr=load('../test1/new_393x153_xy.txt');
m=393;  
n=153;  
     for i=1:m;                                 
      for j=1:n;                                          
        x(j,i)=corr((j-1)*m+i,1)*1.0;               
        y(j,i)=corr((j-1)*m+i,2)*1.0;               
      end                                       
   end  

bathy=load([fdir 'dep_20m_larger.txt']);
lon=load([fdir 'lon_20m_larger.txt']);
lat=load([fdir 'lat_20m_larger.txt']);
[X Y]=meshgrid(lon,lat);
[Lon Lat]=meshgrid(lon,lat);


% interp
bathy_curv=griddata(X,Y,bathy,x,y);

figure(1)
wid=8;
len=10;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid+1 len+1],'paperposition',[0 0 wid len]);
clf
pcolor(x,y,bathy_curv),shading interp;
hold on
sk=1;
xx=x;
yy=y;

xx(bathy_curv>1.0)=NaN;

line(xx(1:sk:end,1:sk:end),yy(1:sk:end,1:sk:end),'Color','k')
line(xx(1:sk:end,1:sk:end)',yy(1:sk:end,1:sk:end)','Color','k')
xlabel('x (m)')
ylabel('y (m)')
grid
demcmap(bathy_curv);
print -djpeg100 grid_bathy_xy.jpg
print -depsc2 grid_bathy_xy.eps

save -ASCII dep_curv.txt bathy_curv
save -ASCII x_curv.txt x
save -ASCII y_curv.txt y

% latlon

lat_curv=griddata(X,Y,Lat,x,y);
lon_curv=griddata(X,Y,Lon,x,y);

figure(2)
wid=8;
len=7.3;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid+1 len+1],'paperposition',[0 0 wid len]);
clf
pcolor(lon_curv,lat_curv,bathy_curv),shading interp;
hold on
sk=1;
lonc=lon_curv;
latc=lat_curv;

lonc(bathy_curv>1.0)=NaN;

line(lonc(1:sk:end,1:sk:end),latc(1:sk:end,1:sk:end),'Color','k')
line(lonc(1:sk:end,1:sk:end)',latc(1:sk:end,1:sk:end)','Color','k')
xlabel('lon (deg)')
ylabel('lat (deg)')
grid
demcmap(bathy_curv);
print -djpeg100 grid_bathy_latlon.jpg
print -depsc2 grid_bathy_latlon.eps

save -ASCII lon_curv.txt lon_curv
save -ASCII lat_curv.txt lat_curv


