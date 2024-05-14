clear all
froot='/Users/fengyanshi/OUTSIDE_Google/Users/NearCoM_Norfolk/';
fmodel='Waves/';
fcase='IR_Base_ERA5/NearCom/';
%fpoint=''waves_model_36.916_-76.4271.csv'';
fdir=[froot fmodel fcase];

fpoints={'waves_model_37.2018_-76.1837.csv'};

figure(1)
hold on
for k=1:length(fpoints)
fname=[fdir fpoints{k}];
Tbl = readtable(fname);
Hs = Tbl.Hs;
Tm= Tbl.tm;
Dm=Tbl.mwd;  
plot(Hs)
text(170,Hs(169),fpoints{k})
end
axis([0 180 0 7])

eval(['mkdir ' 'plots/' fcase]);
print('-djpeg100', ['plots/' fcase 'tmp_1point.jpg'])






