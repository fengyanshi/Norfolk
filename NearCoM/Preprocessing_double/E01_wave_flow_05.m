clear all
% (1) data root
froot='/Users/fengyanshi/OUTSIDE_Google/Users/NearCoM_Norfolk/';

% (2) case --------------------
fcase='IR_SLR_1.3M';
%         ---------------------

% dont need change
fdir_wave=[froot 'Waves/' 'NearCoM_waves_' fcase '/'];
fdir_flow=[froot 'Flow/' 'NearCoM_flow_' fcase '/'];

% (3) wave point
fpoints_wave={'waves_IR_36.9742_-76.2271.csv'};
% (4) flow point
fpoints_flow={'Flow_W_36.993116_-76.229796.csv'};

% (4) output location
foutput=['/Users/fengyanshi/OUTSIDE_Google/GITHUB/Norfolk/NearCoM/' 'Data_double/D_' fcase '/'];


% (5) start and end
% -------------
nstart=80;
nend=105;
% ------------

% (6) SLR
SLR=1.3;

eval(['mkdir ' foutput]);


% wave

fname=[fdir_wave fpoints_wave{1}];
Tbl = readtable(fname);
Time1=Tbl.Date;
Hs1 = Tbl.Hs;
Tp1= Tbl.tp;
Dm1=Tbl.mwd;  
Ds1=Tbl.Dspr; 


Time=Time1(nstart:nend);
Hs = Hs1(nstart:nend) ;
Tp= Tp1(nstart:nend);
Dm=Dm1(nstart:nend);  
Ds=Ds1(nstart:nend); 


set(0,'DefaultFigureColormap',feval('jet'));
figure(1)
subplot(211)
plot(Time1,Hs1)
%text(170,Hs1(169),fpoints{k})
hold on
plot(Time,Hs,'r','LineWidth',2)

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
wavenew{k}=tot;
end

%eval(['mkdir ' 'data/']);
%eval(['mkdir ' 'data/' fcase]);

fid=fopen([foutput 'wave.txt'],'w','n');
fprintf(fid,'%s\n','TPAR');
for k=1:length(wavenew)
%fprintf(fid, '%s\n', strtrim(wavenew(k,:)));
fprintf(fid, '%s\n', wavenew{k});
end
fclose(fid);

% flow ------
fname=[fdir_flow fpoints_flow{1}];
Tbl = readtable(fname);
eta = Tbl.wl;
Wx = Tbl.windx;
Wy = Tbl.windy;
time=Tbl.Date;

time2=time(nstart:nend);
eta2=eta(nstart:nend);

subplot(212)
plot(time,eta)
hold on
plot(time2,eta2,'r','LineWidth',2)

print('-djpeg100', [foutput 'wave_flow.jpg'])

coupling_filename={'coupling.txt'}
logfile='log.txt';

time=([nstart:nend]-nstart)*3600.0;
eastpoints=[1:304];
westpoints=[1:304];

for sp=1:length(eastpoints)
for ti=1:length(time)
eastfine(sp,1,ti)=0.0;  % u
eastfine(sp,2,ti)=0.0;  % v
eastfine(sp,3,ti)=eta(ti+nstart-1)-SLR;  %eta
end
end

southfine=[];
northfine=[];
westfine=[];


filename=[foutput coupling_filename{1}];
FIN = fopen(filename,'w');           

Finfo = fopen(logfile,'w');  
% log file
fprintf(Finfo,'coupling data\nboundary info: num of points, start point');

if isempty(eastfine)==0
fprintf(Finfo,'\nEAST\n\t%d\t\t%d',size(eastfine,1),1);
else
fprintf(Finfo,'\nEAST\n\t%d\t\t%d',-1,1);
end

if isempty(westfine)==0
fprintf(Finfo,'\nWEST\n\t%d\t\t%d',size(westfine,1),1);
else
fprintf(Finfo,'\nWEST\n\t%d\t\t%d',-1,1);
end

if isempty(southfine)==0
fprintf(Finfo,'\nSOUTH\n\t%d\t\t%d',size(southfine,1),1);
else
fprintf(Finfo,'\nSOUTH\n\t%d\t\t%d',-1,1);
end

if isempty(northfine)==0
fprintf(Finfo,'\nNORTH\n\t%d\t\t%d',size(northfine,1),1);
else
fprintf(Finfo,'\nNORTH\n\t%d\t\t%d',-1,1);
end

% end log file

%%
fprintf(FIN,'coupling data\nboundary info: num of points, start point');

if isempty(eastfine)==0
fprintf(FIN,'\nEAST\n\t%d\t\t%d',size(eastfine,1),1);
else
fprintf(FIN,'\nEAST\n\t%d\t\t%d',-1,1);
end

if isempty(westfine)==0
fprintf(FIN,'\nWEST\n\t%d\t\t%d',size(westfine,1),1);
else
fprintf(FIN,'\nWEST\n\t%d\t\t%d',-1,1);
end

if isempty(southfine)==0
fprintf(FIN,'\nSOUTH\n\t%d\t\t%d',size(southfine,1),1);
else
fprintf(FIN,'\nSOUTH\n\t%d\t\t%d',-1,1);
end

if isempty(northfine)==0
fprintf(FIN,'\nNORTH\n\t%d\t\t%d',size(northfine,1),1);
else
fprintf(FIN,'\nNORTH\n\t%d\t\t%d',-1,1);
end

%%


%fprintf(FIN,'coupling data\nboundary info: num of points, start point');
%fprintf(FIN,'\nEAST\n\t%d\t\t%d',size(eastfine,1),1);
%fprintf(FIN,'\nWEST\n\t%d\t\t%d',size(westfine,1),1);
%fprintf(FIN,'\nSOUTH\n\t%d\t\t%d',-1,1);
%fprintf(FIN,'\nNORTH\n\t%d\t\t%d',-1,1);

fprintf(FIN,'\nTIME SERIES');
for t = 1:length(time)
    disp(sprintf('Writing Time Step No. %d    of   %d',t,length(time) ))
    fprintf(FIN,'\n\t%f',time(t));
    printside(FIN,'EAST',eastfine,t)
    printside(FIN,'WEST',westfine,t)
    printside(FIN,'SOUTH',southfine,t)
    printside(FIN,'NORTH',northfine,t)
end
fclose(FIN);
disp('Finished!')

fprintf(Finfo,['\n' 'Total time steps: ' num2str(length(time))]);
fprintf(Finfo,['\n' 'Time interval is AROUND: ' num2str(time(2)-time(1))]);

%clear eta_fine u_fine v_fine eastfine westfine southfine northfine

fclose(Finfo);

% info
fid=fopen([foutput '/info.txt'],'w','n');
fprintf(fid, '%s\n', wavenew{1});
fprintf(fid, '%s\n', wavenew{end});

fclose(fid);

% info
fid=fopen([foutput '/info4plot.txt'],'w','n');
fprintf(fid,'%s\n',Time(1));
fclose(fid);


