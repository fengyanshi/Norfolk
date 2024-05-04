clear all
close all
froot='/Users/fengyanshi/OUTSIDE_Google/Users/NearCoM_Norfolk/';
fmodel='Flow/';
fcase='IR_Base_ERA5';
fdir=[froot fmodel 'NearCoM_flow_' fcase '/'];

fpoints={'Flow_E_36.933015_-76.336177.csv', 'Flow_E_36.956904_-76.336848.csv', 'Flow_N_36.983655_-76.325546.csv', 'Flow_N_36.995337_-76.310359.csv', 'Flow_N_37.009834_-76.299035.csv', 'Flow_W_36.941266_-76.224995.csv', 'Flow_W_37.014133_-76.23519.csv', 'Flow_E_36.942476_-76.336425.csv',	'Flow_E_36.972628_-76.338057.csv', 'Flow_N_36.986576_-76.321749.csv', 'Flow_N_36.998402_-76.306887.csv', 'Flow_N_37.016391_-76.295918.csv', 'Flow_W_36.961676_-76.22595.csv',	'Flow_W_37.078983_-76.266537.csv', 'Flow_E_36.947238_-76.336557.csv', 'Flow_N_36.973921_-76.338202.csv', 'Flow_N_36.989496_-76.317953.csv', 'Flow_N_37.001662_-76.30385.csv',	'Flow_N_37.023493_-76.292604.csv', 'Flow_W_36.97425_-76.227078.csv', 'Flow_E_36.951646_-76.336693.csv', 'Flow_N_36.979647_-76.330292.csv', 'Flow_N_36.992416_-76.314156.csv', 'Flow_N_37.005188_-76.301607.csv', 'Flow_N_37.075656_-76.272944.csv', 'Flow_W_36.993116_-76.229796.csv'};

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
print('-djpeg100', ['plots/' fcase '.jpg'])







