clear all
close all
froot='/Users/fengyanshi/OUTSIDE_Google/Users/NearCoM_Norfolk/';
fmodel='Flow/';
fcase='IR_Base_ERA5';
fdir=[froot fmodel 'NearCoM_flow_' fcase '/'];

% east, north, west
fpoints={'Flow_W_36.993116_-76.229796.csv','Flow_N_36.989496_-76.317953.csv', 'Flow_E_36.956904_-76.336848.csv'};

figure(1)
hold on
for k=1:length(fpoints)
fname=[fdir fpoints{k}];
Tbl = readtable(fname);
eta = Tbl.wl;
Wx = Tbl.windx;
Wy = Tbl.windy;

figure(1)
hold on
plot(eta)
text(170,eta(100),fpoints{k})
figure(2)
hold on
plot(Wx)
plot(Wy)
text(170,Wx(100),fpoints{k})

end
%axis([0 180 0 7])

eval(['mkdir ' 'plots/']);
figure(1)
print('-djpeg100', ['plots/' 'etaWNW' fcase '.jpg'])

figure(2)
print('-djpeg100', ['plots/' 'Wind_WNW' fcase '.jpg'])





