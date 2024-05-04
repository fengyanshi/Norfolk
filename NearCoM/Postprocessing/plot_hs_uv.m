
clear all

fdir='/Users/fengyanshi/TMP/tmp2/';
fdir_input='../grid_new/';
X1=load([fdir_input 'lon_curv.txt']);
Y1=load([fdir_input 'lat_curv.txt']);
dep=load([fdir_input 'dep_curv.txt']);

X=X1(1:end-1,1:end-1);
Y=Y1(1:end-1,1:end-1);

fac=0.002;
%xyaxis=[-76.336 -76.25 36.94 36.98];
xyaxis=[-76.305 -76.265 36.96 36.974]; % small
skx=2;
sky=2;

nstart=59;%input('nstart=');
nend=59;%input('nend=');

wid=10;
len=12;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[2 2 wid len],'paperposition',[0 0 wid len]);

icount=0;
for num=nstart:1:nend
icount=icount+1;

fname=sprintf('%.4d',num);
ele=load([fdir 'eta_' fname]);
mask=load([fdir 'mask_' fname]);
u=load([fdir 'u_' fname]);
v=load([fdir 'v_' fname]);
hs=load([fdir 'hs_' fname]);

uu=sqrt(u.^2+v.^2);
%ww=curl(X,Y,u,v);

[ux,uy]=gradient(u,1,1);
[vx,vy]=gradient(v,1,1);
ww=(uy-vx);

[n,m]=size(ele);

cc=ele;
cc(mask<1)=NaN;
Dep=dep(1:end-1,1:end-1);
B=-Dep;
B(mask==1)=hs(mask==1);
B(mask<1)=NaN;
u(mask<1)=NaN;
v(mask<1)=NaN;
ww(mask<1)=NaN;

u(Dep>-0.5)=NaN;
v(Dep>-0.5)=NaN;

clf
subplot(211)
pcolor(X,Y,B),shading interp
%caxis([0 2.0])
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' Hs (m) ');
%demcmap(B);
hold on
%pcolor(X,Y,cc),shading interp
%colormap(jet)
quiver(X(1:sky:end,1:skx:end),Y(1:sky:end,1:skx:end),u(1:sky:end,1:skx:end)*fac,v(1:sky:end,1:skx:end)*fac,0,'Color','y');

xlabel('Lon(deg) ','fontsize',12,'fontweight','bold');
ylabel('Lat(deg) ','fontsize',12,'fontweight','bold');
axis(xyaxis);
title(['color: wave height' '   (t=' num2str((num-1)*5) ' min)'])

subplot(212)
pcolor(X,Y,cc),shading interp
caxis([-0.1 0.1])
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' \eta (m) ');
hold on
quiver(X(1:sky:end,1:skx:end),Y(1:sky:end,1:skx:end),u(1:sky:end,1:skx:end)*fac,v(1:sky:end,1:skx:end)*fac,0,'Color','k');

xlabel('Lon(deg) ','fontsize',12,'fontweight','bold');
ylabel('Lat(deg) ','fontsize',12,'fontweight','bold');
axis(xyaxis);
title(['color: surface elevation' '  (t=' num2str((num-1)*5) ' min)'])

%title(['Time = '  num2str((num(1))*5) ' min'])
%M(:,icount)=getframe(gcf);
pause(0.1)
end
print -djpeg100 hs_eta_uv_sm.jpg


