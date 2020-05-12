%use this to shorten the TAW diary
%
clc;clear all;close all
%%%%%%%%%%%%%%%  CONFIG  %%%%%%%%%%%%%%%%%%%%%
datafile='../data/transectdata.xls'; 
lenTAWfile=401;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,txt,raw]=xlsread(datafile);
for i=2:size(raw,1)
fnames{i-1}=raw{i,1};
end

line=0;

for i=1:length(fnames)
    clear line
    name=['logfiles/' fnames{i} '-DIARY.txt']
    fid=fopen(name,'r');
    n=1;
    line{n}=fgetl(fid);
    while ischar(line{n})
        n=n+1;
        line{n}=fgetl(fid);
    end
    fclose all
    
    %print
    fid=fopen(name,'w');
    for n=1:length(line)
        if length(line{n}) == 0 && n >= lenTAWfile
        else
           fprintf(fid,'%s\n',line{n});
        end
    end
        
end

fclose all 
