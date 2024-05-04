clear all
%wave=load('wavecon.n1f');
% time
dt=5.0; % every 10 sec
% Hs
Hs=2.0;
Ig=1.8;
%
Tp=8.0;
%
Wdir=300.0;
%
Wspr=15.0;
% envelope period 
TT = 60.0;
% time
hours=3.0;
nmin=hours*3600.0/dt; %  every 10 sec 3hr 3*3600/10
wave=zeros(nmin,5);
wave(:,1)=[0:nmin-1];
tmp_time=wave(:,1);
wave_time=datenum('June 1, 2012 00:00:00.000 AM')+tmp_time*dt/24/3600.0;
wave_date=datevec(wave_time);
for k=1:length(wave_date)
   nyear=int2str(wave_date(k,1));
     if wave_date(k,2)<10
    nmonth=['0' int2str(wave_date(k,2))];
   else
    nmonth=int2str(wave_date(k,2));      
     end
      if wave_date(k,3)<10
    ndate=['0' int2str(wave_date(k,3))];
   else
    ndate=int2str(wave_date(k,3))      ;
      end
     if wave_date(k,4)<10
    nhour=['0' int2str(wave_date(k,4))];
   else
    nhour=int2str(wave_date(k,4))      ;
     end
     if wave_date(k,5)<10
    nminute=['0' int2str(wave_date(k,5))];
   else
    nminute=int2str(wave_date(k,5))      ;
     end
     if wave_date(k,6)<10
    nsecond=['0' int2str(wave_date(k,6))];
   else
    nsecond=int2str(wave_date(k,6))      ;
     end
     yeartime(k,:)=[nyear nmonth ndate '.' nhour nminute nsecond];
end


%plot(wind_date,wind(:,2));
% envelope period TT
for kk=1:length(wave)
tt=(kk-1)*dt;
wave_height(kk)=Hs+Ig*sin(2.0*pi/TT*tt); % 
end

wave_period(1:length(wave))=Tp;
wave_angle(1:length(wave))=Wdir;
wave_spreading(1:length(wave))=Wspr;

formatSpec = '%.2f';

for k=1:length(wave)
   yt=yeartime(k,1:15);
   hh=num2str(wave_height(k),formatSpec);
   wp=num2str(wave_period(k),formatSpec);
   wa=num2str(wave_angle(k),formatSpec);
   ws=num2str(wave_spreading(k),formatSpec);
   tot=[yt ' ' hh ' ' wp ' ' wa ' ' ws];
   wavenew(k,1:length(tot))=tot;
end

fid=fopen('wave_west.txt','w','n');
fprintf(fid,'%s\n','TPAR')
for k=1:length(wavenew)
fprintf(fid, '%s\n', strtrim(wavenew(k,:)));
end
fclose(fid);

fid=fopen('wave_north.txt','w','n');
fprintf(fid,'%s\n','TPAR')
for k=1:length(wavenew)
fprintf(fid, '%s\n', strtrim(wavenew(k,:)));
end
fclose(fid);

fid=fopen('wave_south.txt','w','n');
fprintf(fid,'%s\n','TPAR')
for k=1:length(wavenew)
fprintf(fid, '%s\n', strtrim(wavenew(k,:)));
end
fclose(fid);

