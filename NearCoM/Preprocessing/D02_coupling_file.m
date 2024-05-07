clear all
close all

froot='/Users/fengyanshi/OUTSIDE_Google/Users/NearCoM_Norfolk/';
fmodel='Flow/';
fcase='IR_Base_ERA5';
fdir=[froot fmodel 'NearCoM_flow_' fcase '/'];
foutput='/Users/fengyanshi/OUTSIDE_Google/GITHUB/Norfolk/NearCoM/Preprocessing/data/IR_Base_ERA5_NBC/';

% east, north, west
fpoints={'Flow_W_36.993116_-76.229796.csv'};

nstart=80;
nend=105;

eval(['mkdir ' foutput])
% input

fname=[fdir fpoints{1}];
Tbl = readtable(fname);
eta = Tbl.wl;
Wx = Tbl.windx;
Wy = Tbl.windy;
time=Tbl.Date;

time2=time(nstart:nend);
eta2=eta(nstart:nend);

figure(1)
hold on
plot(time,eta)
%text(170,eta(100),fpoints{k})
hold on
plot(time2,eta2,'r','LineWidth',2)

figure(1)
print('-djpeg100', ['plots/' 'eta' fcase '.jpg'])

%figure(2)
%hold on
%plot(Wx)
%plot(Wy)
%text(170,Wx(100),fpoints{k})


coupling_filename={'coupling.txt'}
logfile='log.txt';


time=([nstart:nend]-nstart)*3600.0;
eastpoints=[1:152];
westpoints=[1:152];

for sp=1:length(eastpoints)
for ti=1:length(time)
eastfine(sp,1,ti)=0.0;  % u
eastfine(sp,2,ti)=0.0;  % v
eastfine(sp,3,ti)=eta(ti+nstart-1);  %eta
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
%disp(sprintf('NOTE: This coupling file starts at time = %d sec',sample(start,1)))
%disp(sprintf('      ends at time = %d sec',sample(stop,1)))

%fprintf(Finfo,['\n' sprintf('NOTE: This coupling file starts at time = %d sec',sample(start,1))]);
%fprintf(Finfo,['\n' sprintf('      ends at time = %d sec',sample(stop,1))]);
fprintf(Finfo,['\n' 'Total time steps: ' num2str(length(time))]);
fprintf(Finfo,['\n' 'Time interval is AROUND: ' num2str(time(2)-time(1))]);

%clear eta_fine u_fine v_fine eastfine westfine southfine northfine

fclose(Finfo);



