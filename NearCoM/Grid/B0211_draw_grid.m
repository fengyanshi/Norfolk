clear all

lon_curv=load('lon_curv.txt');
lat_curv=load('lat_curv.txt');
bathy_curv=load('dep_curv.txt');

ax1=[-76.345 -76.25 36.90 36.98];

figure(1)
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
plot([ax1(1) ax1(2) ax1(2) ax1(1) ax1(1)],[ax1(3) ax1(3) ax1(4) ax1(4) ax1(3)],'w-','LineWidth',1)
demcmap(bathy_curv);
print -djpeg100 grid_bathy_latlon_B.jpg
print -depsc2 grid_bathy_latlon_B.eps

% zoom
wid=8;
len=4.3;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid+1 len+1],'paperposition',[0 0 wid len]);
axis([-76.335 -76.25 36.938 36.975])
print -djpeg100 grid_bathy_latlon_zoom.jpg
print -depsc2 grid_bathy_latlon_zoom.eps

% write boundary points
clear west_b east_b north_b south_b
sk=30;
% west
west_b(:,1)=lat_curv(1:sk:end,end);
west_b(:,2)=lon_curv(1:sk:end,end);
% east
east_b(:,1)=lat_curv(1:sk:end,1);
east_b(:,2)=lon_curv(1:sk:end,1);
% north
north_b(:,1)=lat_curv(end,1:sk:end);
north_b(:,2)=lon_curv(end,1:sk:end);
% south
south_b(:,1)=lat_curv(1,1:sk:end);
south_b(:,2)=lon_curv(1,1:sk:end);

figure(2)

wid=8;
len=7.3;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid+1 len+1],'paperposition',[0 0 wid len]);
clf
plot(west_b(:,2),west_b(:,1),'ro-')
hold on
plot(east_b(:,2),east_b(:,1),'ro-')
plot(south_b(:,2),south_b(:,1),'ro-')
plot(north_b(:,2),north_b(:,1),'ro-')
contour(lon_curv,lat_curv,bathy_curv,[0 0],'Color','k')
xlabel('lon (deg)')
ylabel('lat (deg)')
grid

print -djpeg100 boundary.jpg
print -depsc2 boundary.eps

save -ASCII west.txt west_b
save -ASCII east.txt east_b
save -ASCII south.txt south_b
save -ASCII north.txt north_b


