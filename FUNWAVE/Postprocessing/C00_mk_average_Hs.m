
clear all
close all

froot='/Volumes/BigSur_2022/Norfolk_FUNWAVE/Results/';
fcase={'R1_smdm','R2_smdm','R3_smdm','R4_smdm','R5_smdm','R6_smdm'};
case_name={'ERA5','HM','Acc GN 1m','RMW F1.25','SLR 1.3m','WSF 1.225'};


%fdir_input='../Grid/';
%fdir_data=['../Data/' fcase_whole{1} '/'];

files=[25];
dt=100.0;

ax=[100 5560 100 4300];
cax=[0 2.5];
Mindep=0.1;
uv_deplimit = 0.2;

% read info

m=7744;
n=5888;

DimsX={[m n]};
% dimensions
dx=0.75;
dy=0.75;
x=[0:m-1]*dx;
y=[0:n-1]*dy;
sk=4;


fig=figure(1);
colormap parula
clf

wid=10;
len=10;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);


numb=files(1);

fnum=sprintf('%.5d',numb);

for k=1:length(fcase)
% read files -----------------------

fdir=[froot fcase{k} '/'];

fname=[fdir 'dep.out'];
fileID=fopen(fname);
dep=fread(fileID,DimsX{1},'*single');
fclose(fileID);
dep=dep';

fname=[fdir 'Hsig_' fnum];
fileID=fopen(fname);
hs=fread(fileID,DimsX{1},'*single');
fclose(fileID);
hs=hs';


hs(dep<0.0)=NaN;


subplot(3,2,k)
pcolor(x(1:sk:end),y(1:sk:end),hs(1:sk:end,1:sk:end)),shading flat
time1=[case_name{k}];
title(time1)
caxis(cax)
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' Hs (m) ');

if(k==5 || k==6)
xlabel('x (m, model coordinate)')
end
if(k==1 || k==3 ||k==5)
ylabel('y (m)')
end
axis equal
axis(ax)

end

eval(['mkdir ' 'plots/six_cases'])

%pname=['plots/six_cases' '/' 'hs.jpg'];
%print('-djpeg100',pname);




