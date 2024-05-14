clear all
close all
froot='/Users/fengyanshi/OUTSIDE_Google/Users/NearCoM_Norfolk/';
fmodel='Waves/';
fcase='IR_Base_ERA5_NBC';
fdir=[froot fmodel 'NearCoM_waves_' fcase '/'];

fpoints={'waves_IR_36.603_-74.837.csv',	'waves_IR_36.9895_-76.318.csv','waves_IR_36.621_-75.8509.csv',	'waves_IR_36.9924_-76.3142.csv','waves_IR_36.933_-76.3362.csv',	'waves_IR_36.9931_-76.2298.csv','waves_IR_36.9413_-76.225.csv',	'waves_IR_36.9953_-76.3104.csv','waves_IR_36.9425_-76.3364.csv',	'waves_IR_36.9984_-76.3069.csv','waves_IR_36.9472_-76.3366.csv',	'waves_IR_37.0017_-76.3038.csv','waves_IR_36.9516_-76.3367.csv',	'waves_IR_37.0052_-76.3016.csv','waves_IR_36.9569_-76.3368.csv',	'waves_IR_37.0098_-76.299.csv','waves_IR_36.9617_-76.2259.csv',	'waves_IR_37.0141_-76.2352.csv','waves_IR_36.9726_-76.3381.csv',	'waves_IR_37.0164_-76.2959.csv','waves_IR_36.9739_-76.3382.csv',	'waves_IR_37.0235_-76.2926.csv','waves_IR_36.9742_-76.2271.csv',	'waves_IR_37.0757_-76.2729.csv','waves_IR_36.9796_-76.3303.csv',	'waves_IR_37.079_-76.2665.csv','waves_IR_36.9837_-76.3255.csv',	'waves_IR_38.46_-74.692.csv','waves_IR_36.9866_-76.3217.csv'};

figure(1)
hold on
for k=1:length(fpoints)
fname=[fdir fpoints{k}];
Tbl = readtable(fname);
Hs = Tbl.Hs;
Hs(Hs<-0.)=NaN;
Tp= Tbl.tp;
Tp(Tp<0.)=NaN;
Dm=Tbl.mwd; 
Dm(Dm<0.)=NaN;
Ds=Tbl.Dspr; 
Ds(Ds<0.)=NaN;

figure(1)
hold on
plot(Hs)
text(170,Hs(100),fpoints{k})
figure(2)
hold on
plot(Tp)
text(170,Tp(100),fpoints{k})
figure(3)
hold on
plot(Dm)
text(170,Dm(100),fpoints{k})

end
%axis([0 180 0 7])

eval(['mkdir ' 'plots/']);
%print('-djpeg100', ['plots/' fcase '.jpg'])







