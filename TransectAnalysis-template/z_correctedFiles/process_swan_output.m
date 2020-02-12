%SWAN OUTPUT
%
%Author: Scott Hayward
%Company: Ransom Consulting, inc. 
%Project: 2018 FEMA appeal, York and Cumberland Counties 
%
%this script processes the results of each SWAN model located in the
%transectdata.xls spreadsheet. Results are plotted in comaprison to the
%ADCIRC model, and the wave height, period, and setup are extracted at the
%toe of the runup profile. Additonally, the maximum wave setup is saved
%into the excel spreadsheet.

% chk nld 20190912

clc; clear all;close all;fclose all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datafile='../data/transectdata.xls'; 
tDIR='../ADCIRC_returns/'; %location of transects
logpre='logfiles/';
swanpre='swanfiles/'; %prefix for swan file output
L_append=2;
savefile='../data/transectdata.xls'; imgfile='tpng'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,txt,raw]=xlsread(datafile);
for i=2:size(raw,1)
fnames{i-1}=raw{i,1};
end
wha_start=num(:,1);
wha_twl=num(:,2);
r_start=num(:,6);

for i=1:length(fnames)

logfile=[logpre fnames{i} '_log.txt'];
fname=[tDIR fnames{i} 'XYZSTA_RETURNS.csv'];
file=xlsread(fname);
lon{i}=file(:,1);lat{i}=file(:,2);sta=file(:,4);ele=file(:,3);adc_wl=file(:,7);adc_hs=file(:,6);

%file id's
fid=fopen(logfile,'a');
fid1=fopen([swanpre fnames{i} '.dat'],'r');

%read SWAN .dat file output
li=fgetl(fid1);li=fgetl(fid1);li=fgetl(fid1);li=fgetl(fid1);
li=fgetl(fid1);li=fgetl(fid1);li=fgetl(fid1);li=fgetl(fid1);
n=0;
clear('xswan','tp','hsswan','setup');
while li ~=-1
    n=n+1;
    A=str2num(li);
    setup(n)=A(10)*3.280833333;
    hsswan(n)=A(3);
    tp(n)=A(4);
    xswan(n)=(n-1)*3.280833333+wha_start(i);
    li=fgetl(fid1);
end

%plot adcirc, and SWAN outputs
figure
adc_wl(adc_wl<0)=NaN;
adc_hs(adc_hs<0)=NaN;
plot(sta,ele,'k');hold on
plot(sta,adc_wl,'b');
plot(xswan,setup+wha_twl(i),'g');
plot(sta,adc_hs,'b--')
plot(xswan,hsswan*3.28083,'g--');
title({'{\color{red}REVISED SEP-05-2019}','2-D ADCIRC+SWAN and SWAN 1-D results, Transect: ' fnames{i} })
ylim([0 20]);xlabel('Station (ft)');ylabel('Elevation, TWL, HS (ft)');
legend('Transect Profile', 'TWL (ADC 2D)', 'TWL (SWN 1D)', 'Hs (ADC 2D)', 'Hs (SWN 1D)','location','northwest')
% xlim([xswan(1) xswan(end)])
grid minor
set(gcf,'position',[100 100 900 600],'paperorientation','landscape');
print('-r450','-dpdf',[logpre fnames{i}]);
pause(.2)

r_hs(i)=interp1(xswan,hsswan,r_start(i),'nearest')*3.280833333;
r_per(i)=interp1(xswan,tp,r_start(i),'nearest');
r_setup(i)=interp1(xswan,setup,r_start(i),'nearest');
max_setup(i)=max(setup);

fprintf(fid,'%s\n',  ['SWAN maximum additional wave setup: ' num2str(max_setup(i)) ' feet']);


fprintf(fid,'%s\n',  ['SWAN output at toe:']);
fprintf(fid,'%s\n',  ['    SETUP- ' num2str(r_setup(i)) ' feet' ]);
fprintf(fid,'%s\n',  ['    HS-    ' num2str(r_hs(i)) ' feet' ]);
fprintf(fid,'%s\n\n',['    PER-   ' num2str(r_per(i)) ' seconds']);

fprintf(fid,'%s\n\n\n','PART 2 COMPLETE________________________________________');

%input filename
%swan filename
%other

% interpolate, and write the transect input file (1-foot increments)

fclose all;
end

%add output to excel file

for i=1:length(fnames)
    
    %write the transect data into excel file 
    d1={max_setup(i)};
    d2={r_setup(i)};
    d3={r_hs(i)};
    d4={r_per(i)};
    xlswrite(savefile,d1,'Sheet1',['F' num2str(L_append+i-1)]);
    xlswrite(savefile,d2,'Sheet1',['L' num2str(L_append+i-1)]);
    xlswrite(savefile,d3,'Sheet1',['M' num2str(L_append+i-1)]);
    xlswrite(savefile,d4,'Sheet1',['N' num2str(L_append+i-1)]);
    
end
