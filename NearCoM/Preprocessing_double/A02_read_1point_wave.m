clear all
froot='/Users/fengyanshi/OUTSIDE_Google/Users/NearCoM_Norfolk/';
fmodel='Waves/';
fcase='IR_Base_ERA5_NBC';
fdir=[froot fmodel 'NearCoM_waves_' fcase '/'];

fpoints={'waves_IR_36.9742_-76.2271.csv'};

figure(1)
hold on
for k=1:length(fpoints)
fname=[fdir fpoints{k}];
Tbl = readtable(fname);
Hs = Tbl.Hs;
Tp= Tbl.tp;
Dm=Tbl.mwd; 
Ds=Tbl.Dspr; 

figure(1)
plot(Hs)
text(170,Hs(100),fpoints{k})
figure(2)
plot(Tp)
text(170,Tp(100),fpoints{k})
figure(3)
plot(Dm)
text(170,Dm(100),fpoints{k})

end
%axis([0 180 0 7])

eval(['mkdir ' 'plots/']);
%print('-djpeg100', ['plots/' fcase '.jpg'])







