clc;clear all;close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datafile='../data/transectdata.xls'; 
datafile='../data/transectdata.xls';
tDIR='../ADCIRC_returns/'; %location of transects
imgfile='tpng'; runupname='CM-runup';
picktoetop=0; %to pick toetop for all stations, set equal to 1
L_append=2;
engineer='sjh';
jobstr='job 2';
% filename='YK_R2_T.in';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,txt,raw]=xlsread(datafile);
for i=2:size(raw,1)
Tnames{i-1}=raw{i,1};
client{i-1}='FEMA';
eng{i-1}=engineer;
job{i-1}=jobstr;
run{i-1}=1;
desc{i-1}=['RUNUP2 transect: ' Tnames{i-1}];
filename{i-1}=[Tnames{i-1} '.in'];
rough{i-1}=num(i-1,15);
toe{i-1}=num(i-1,1); %using WHAFIS start rather than TAW
top{i-1}=num(i-1,8);
wave{i-1}=num(i-1,3);%*0.626; % Multiply Hs by 0.626 to get Hmean
twl{i-1}=num(i-1,2); %This is the twl 
per{i-1}=num(i-1,4);%*0.85; % convert significant wave period from SWAN to mean period
end



%read data from each transect
for i=1:size(Tnames,2)
%     fname=[tDIR fnames{i} 'XYZSTA_RETURNS.csv'];
    fname=[tDIR Tnames{i} 'XYZSTA_RETURNS.csv'];
    file=xlsread(fname);
    lon{i}=file(:,1);lat{i}=file(:,2);sta{i}=file(:,4);ele{i}=file(:,3);
    sta{i}(isnan(ele{i}))=[];
    ele{i}(isnan(ele{i}))=[];
    lon{i}(isnan(ele{i}))=[];
    lat{i}(isnan(ele{i}))=[];
end
%Create waterlevel input variable
for i=1:length(wave)
    wlvl{i}=[twl{i} wave{i} per{i}];
end


fbat=fopen('RUNUP2_files/run.bat','w');

for i=1:size(Tnames,2)
    [output{i},ele_fill(i),tslope{i},xsct(i),indices{i},offsets(i),Hm_deshoal(i),D_shallow(i),strt(i),sta_shift(i),HuntStr{i}]=runup2_format(client(i),...
        eng(i),job(i),desc(i),run(i),sta(i),ele(i),rough(i),wlvl(i),...
        filename{i},toe(i),top(i),Tnames(i));

    outfile=strrep(filename{i},'in','out')
    fprintf(fbat,'RUNUP2 %s %s\n',filename{i},outfile);
end

fclose(fbat)




%% Open DOSBOX, mount drive
    % for example, this current directory may be E:/Projects/FEMA_appeal/Kittery/Runup2
    % In Dosbox, type: 
    % mount E E:/Projects/FEMA_appeal/Kittery/Runup2
    % will mount this folder to the E drive in Dosbox. 
    % e:
    % cd to working directory
    % run run.bat
        
% Next, Run the RUNUP2 Format script for each transect, create a different file for each one

