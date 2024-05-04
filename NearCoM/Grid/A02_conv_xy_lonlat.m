clear all

lon=load('x_curv.txt');
lat=load('y_curv.txt');

lon_m=min(min(lon));
lat_m=min(min(lat));

R=6371000.0;

x=(lon-lon_m)*R*pi/180.0*cos(pi/180.0*36.6);
y=(lat-lat_m)*R*pi/180.0;

save -ASCII xx_curv.txt x
save -ASCII yy_curv.txt y