%SWAN MODEL 
%
%Author: Scott Hayward
%Company: Ransom Consulting, inc. 
%Project: 2018 FEMA appeal, York and Cumberland Counties 
%
%A SWAN model is created for each transect specified in transectdata.xls
%spreadsheet. The transect is interpolated onto a 1-meter grid. Imperial
%units are converted to metric, and written into SWAN files. In Addition, A
%windows .bat file is written to automatically execute all SWAN
%computations. 
%
%When complete, run: process_swan_output.m to plot all results, and write
%into the transectdata.xls spreadsheet. 
%
%chk nld 20181003

clc;clear all;close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datafile='../data/transectdata.xls'; 
tempfile='FEMAtemplate.txt';
tDIR='../ADCIRC_returns/'; %location of transects
logpre='logfiles/';
elevpre='gridfiles/'; %prefix for grid file
swanpre='swanfiles/'; %prefix for swan file
wvel='25.1'; %wind speed, meters per second. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,txt,raw]=xlsread(datafile);
for i=2:size(raw,1)
fnames{i-1}=raw{i,1};
end
wha_start=num(:,1);
wl=num(:,2)/3.280833333;
hs=num(:,3)/3.280833333;
per=num(:,4);



%read full template file. 
fid=fopen(tempfile);
line=fgetl(fid);
i=1;templatefile={};
while line ~= -1 
templatefile{end+1}=line;
line=fgetl(fid);
end

for i=1:length(fnames)

logfile=[logpre fnames{i} '_log.txt']
fid=fopen(logfile,'w');
fname=[tDIR fnames{i} 'XYZSTA_RETURNS.csv'];
file=xlsread(fname);
lon{i}=file(:,1);lat{i}=file(:,2);sta=file(:,4);ele=file(:,3);
% convert ele to meters
ele=ele./3.280833333;

fprintf(fid,'%s\n','_______________________________________________________');
fprintf(fid,'%s\n','PART 2: SWAN 1-D');
fprintf(fid,'%s\n','');

% interpolate, and write the transect input file (1-m increments)

%new transect from whafis start to shoreline
trim=ele;trim(find(sta<wha_start(i)))=-999;
x=wha_start(i):3.280833333:sta(min(find(trim>wl(i)))); %go a little past water line, in meters
z{i}=interp1(sta,ele,x);
elevfile=[elevpre fnames{i} 'zmeters_xmeters.grd'];
fid2=fopen(elevfile,'w');
dep=z{i}-wl(i); % water column depth in feet
for zi=1:length(z{i})
    fprintf(fid2,'%6.2f\n',dep(zi)); %elevation in meters
end
fid2=fopen([swanpre fnames{i} '.swn'],'w');

% read the template file, and replace each variable
for li=1:length(templatefile)
    swnfile{li}=templatefile{li}; %read each line of template
    
%     %$TWL$    (total water elevation (m))
%     swnfile{li}=strrep(swnfile{li},'$TWL$',num2str(wl(i)/3.28083));
    %$XLENC$  (length of transect (m))
    swnfile{li}=strrep(swnfile{li},'$XLENC$',num2str(length(x)-1));
    %$MXC$    (number of mesh cells (one less than number of points))
    swnfile{li}=strrep(swnfile{li},'$MXC$',num2str(length(x)-1));
    %$MXINP$  (number of mesh points (length of input file))
    swnfile{li}=strrep(swnfile{li},'$MXINP$',num2str(length(x)-1));
    %$FNAME1$ (name of 1-d transect grid)
    swnfile{li}=strrep(swnfile{li},'$FNAME1$',['../' elevfile]);
    %$VEL$    (Wind Velocity (m/s))
    swnfile{li}=strrep(swnfile{li},'$VEL$',wvel);
    %$DIR$    (Wind Direction)
    swnfile{li}=strrep(swnfile{li},'$DIR$','0'); %double-check this later. 
    %$HSIG$   (Significant Wave Height (m) )
    swnfile{li}=strrep(swnfile{li},'$HSIG$',num2str(hs(i))); 
    %$TPEAK$  (Peak wave Period (s))
    swnfile{li}=strrep(swnfile{li},'$TPEAK$',num2str(per(i)));
    %$FNAME$  (output file name)
    swnfile{li}=strrep(swnfile{li},'$FNAME$',[fnames{i} '.dat']);
    fprintf(fid2,'%s\n',swnfile{li});
end

%batch file

fprintf(fid,'%s\n',  ['swan input grid name: 2_swan/' elevfile]);
fprintf(fid,'%s\n',  ['swan file name:       2_swan/' swanpre fnames{i} '.swn']);
fprintf(fid,'%s\n\n',['swan output name:     2_swan/' swanpre fnames{i} '.dat']);

fprintf(fid,'%s\n',  ['Boundary Conditions:']);
fprintf(fid,'%s\n',  ['    TWL- ' num2str(wl(i)) ' meters' ]);
fprintf(fid,'%s\n',  ['    HS-  ' num2str(hs(i)) ' meters' ]);
fprintf(fid,'%s\n\n',['    PER- ' num2str(per(i)) ' seconds']);

fprintf(fid,'%s\n\n',['Batch File: 2_swan/' swanpre 'runswan.dat']);
% fprintf(fid,'%s\n\n\n','PART 2 COMPLETE________________________________________');



fclose all;
end

bfname='swanfiles/runswan.bat';
fid3=fopen(bfname,'w');
fprintf(fid3,'%s\n',['@ECHO OFF']);
for i=1:length(fnames)
    fnamein=[fnames{i}];
    fprintf(fid3,'%s\n',['call swanrun ' fnamein]);
end
fprintf(fid3,'%s\n',['@END']);

fclose all
    



