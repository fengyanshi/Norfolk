clear all
froot='/Users/fengyanshi/OUTSIDE_Google/Users/NearCoM_Norfolk/';
fmodel='Waves/';
fcase='IR_Base_ERA5_NBC';
fdir=[froot fmodel 'NearCoM_waves_' fcase '/'];

fpoints={'waves_IR_36.9742_-76.2271.csv'};

nstart=80;
nend=105;

figure(1)
hold on
for k=1:length(fpoints)
fname=[fdir fpoints{k}];
Tbl = readtable(fname);
Time1=Tbl.Date;
Hs1 = Tbl.Hs;
Tp1= Tbl.tp;
Dm1=Tbl.mwd;  
Ds1=Tbl.Dspr; 
plot(Time1,Hs1)
text(170,Hs1(169),fpoints{k})

Time=Time1(nstart:nend);
Hs = Hs1(nstart:nend) ;
Tp= Tp1(nstart:nend);
Dm=Dm1(nstart:nend);  
Ds=Ds1(nstart:nend); 

plot(Time,Hs,'r','LineWidth',2)

end
%axis([0 180 0 7])

eval(['mkdir ' 'plots/']);
print('-djpeg100', ['plots/' fcase '_1point.jpg'])

% write
for k=1:length(Time)
   nyear=num2str(year(Time(k)),'%.4d');
   nmonth=num2str(month(Time(k)),'%.2d');
   ndate=num2str(day(Time(k)),'%.2d');
   nhour=num2str(hour(Time(k)),'%.2d');
   nminute=num2str(minute(Time(k)),'%.2d');
   yeartime=[nyear nmonth ndate '.' nhour nminute '00'];
   hh=num2str(Hs(k),'%4.1f');
   wp=num2str(Tp(k),'%4.1f');
   wa=num2str(Dm(k),'%4.1f');
   ws=num2str(Ds(k),'%4.1f');
   tot=[yeartime ' ' hh ' ' wp ' ' wa ' ' ws];
%   wavenew(k,1:length(tot))=tot;
wavenew{k}=tot;
end

eval(['mkdir ' 'data/']);
eval(['mkdir ' 'data/' fcase]);

fid=fopen(['data/' fcase '/' 'wave.txt'],'w','n');
fprintf(fid,'%s\n','TPAR');
for k=1:length(wavenew)
%fprintf(fid, '%s\n', strtrim(wavenew(k,:)));
fprintf(fid, '%s\n', wavenew{k});
end
fclose(fid);





