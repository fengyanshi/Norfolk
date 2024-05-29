clear all
close all

fa=load('plots/six_cases/flood_area.txt');
flooda=fa';
case_name={'ERA5','HM','Acc GN 1m','RMW F1.25','SLR 1.3m','WSF 1.225'};

fig=figure(1);
clf
wid=8;
len=5;
set(fig,'units','inches','paperunits','inches','papersize',...
    [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
y=[flooda];
bar(y,0.5)
ylabel('flood area (m^2)');
xlabel('scenarios')
grid
%text(1,9.6,['Relative RMSE = ' num2str(RRMSE*100,3) '%' ])
set(gca,'XtickLabel',case_name)

print('-djpeg100', 'plots/six_cases/nearcom_bars.jpg');