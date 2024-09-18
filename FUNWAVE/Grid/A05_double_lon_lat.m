clear all
fdir='/Users/fengyanshi/WORK/work/Chris_Norfolk/old/';
f_drive='/Volumes/BigSur_2022/Norfolk_FUNWAVE/';

% used 7871x6098 grid and dep1=dep(1:3000,1472:end);

load([fdir 'lon_lat_6400x3000']);

lon_new=double_grid(lon);
lat_new=double_grid(lat);

% small domain
lon_sm1=lon_new(1:5888,1:12544);
lat_sm1=lat_new(1:5888,1:12544);

lon_sm=lon_sm1(1:5888,3800:end-1001);
lat_sm=lat_sm1(1:5888,3800:end-1001);

save([f_drive 'lon_lat_7744x5888.mat'], 'lon_sm','lat_sm');