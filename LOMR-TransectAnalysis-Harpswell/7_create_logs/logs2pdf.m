
clc;clear all;close all
%%%%%%%%%%%%%%%  CONFIG  %%%%%%%%%%%%%%%%%%%%%
exefile='text2pdf.win95.exe';
%see text2pdf documentation: https://github.com/philips/text2pdf 

datafile='../data/transectdata.xls'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,txt,raw]=xlsread(datafile);
for i=2:size(raw,1)
fnames{i-1}=raw{i,1};
end

fid=fopen('convert_logs.bat','w');
fprintf(fid,'%s\n','@ECHO OFF')

for i=1:length(fnames)
    %create text file that executes text2pdf
    fprintf(fid,'%s\n',[exefile ' -f"Courier" -s10 -c216 -v10 -x<612> -y<792> '  '../1_input/logfiles/'  fnames{i} '_log.txt > '   '../1_input/logfiles/'  fnames{i} '_log.pdf']);
    fprintf(fid,'%s\n',[exefile ' -f"Courier" -s10 -c216 -v10 -x<612> -y<792> '  '../2_swan/logfiles/'   fnames{i} '_log.txt > '   '../2_swan/logfiles/'   fnames{i} '_log.pdf']);
    fprintf(fid,'%s\n',[exefile ' -f"Courier" -s8 -c216 -v8 -x<612> -y<792> ' '../2_swan/swanfiles/' fnames{i} '.prt > '        '../2_swan/logfiles/'   fnames{i} '_1swanprint.pdf']);
    fprintf(fid,'%s\n',[exefile ' -f"Courier" -s8 -c216 -v8 -x<792> -y<612>  ' '../2_swan/swanfiles/' fnames{i} '.dat > '        '../2_swan/logfiles/'   fnames{i} '_3swanout.pdf']);
    fprintf(fid,'%s\n',[exefile ' -f"Courier" -s10 -c216 -v10 -x<612> -y<792> '  '../3_whafis/logfiles/' fnames{i} '_log.txt > '   '../3_whafis/logfiles/' fnames{i} '_log.pdf']);
    fprintf(fid,'%s\n',[exefile ' -f"Courier" -s7 -c216 -v7 -x<612> -y<792> ' '../3_whafis/logfiles/' fnames{i} '_short.out > ' '../3_whafis/logfiles/' fnames{i} '_output.pdf']);
    fprintf(fid,'%s\n',[exefile ' -f"Courier" -s10 -c216 -v10 -x<612> -y<792> '  '../4_taw/logfiles/'    fnames{i} '_log.txt > '   '../4_taw/logfiles/'    fnames{i} '_log.pdf']);
    fprintf(fid,'%s\n',[exefile ' -f"Courier" -s8 -c216 -v8 -x<612> -y<792> '    '../4_taw/logfiles/'      fnames{i} '-DIARY.txt > ' '../4_taw/logfiles/'    fnames{i} '-DIARY.pdf']);
    fprintf(fid,'%s\n',[exefile ' -f"Courier" -s10 -c216 -v10 -x<612> -y<792> '  '../5_runup2/logs/'     fnames{i} '_log.txt > '   '../5_runup2/logs/' fnames{i} '_log.pdf']);
    fprintf(fid,'%s\n',[exefile ' -f"Courier" -s8 -c216 -v8 -x<612> -y<792> '    '../5_runup2/RUNUP2_files/' fnames{i} '.in > '    '../5_runup2/logs/' fnames{i} '_in.pdf']);
    fprintf(fid,'%s\n',[exefile ' -f"Courier" -s8 -c216 -v8 -x<792> -y<612> '  '../5_runup2/RUNUP2_files/' fnames{i} '.OUT > '  '../5_runup2/logs/' fnames{i} '_out.pdf']);
      
end

fclose all 
