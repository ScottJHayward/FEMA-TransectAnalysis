
clc;clear all;close all
%%%%%%%%%%%%%%%  CONFIG  %%%%%%%%%%%%%%%%%%%%%
datafile='../data/transectdata.xls'; 
%%%%%%%%%%%%% end config %%%%%%%%%%%%%%%%%%%%%

[num,txt,raw]=xlsread(datafile);
for i=2:size(raw,1)
    fnames{i-1}=raw{i,1};
    TAWvalid(i-1)=raw{i,21};
end

for i=1:length(fnames)
    %delete file directory, re-create it
    if exist(fnames{i})~=0
        delete([fnames{i} '/*']);
        rmdir(fnames{i});
    end
    mkdir(fnames{i})
    
    %now copy all files into the directory we created.
    copyfile(['../1_input/logfiles/' fnames{i}  '_title.pdf'],     [fnames{i},'/section 0.0- title page.pdf']);
    copyfile(['../1_input/logfiles/' fnames{i}  '_log.pdf'],       [fnames{i},'/section 1.0- user input.pdf']);
    copyfile(['../1_input/logfiles/' fnames{i}  '_data.pdf'],      [fnames{i},'/section 1.1- user input figure.pdf']);
    copyfile(['../2_swan/logfiles/' fnames{i}   '_log.pdf'],       [fnames{i},'/section 2.0- swan.pdf']);
    copyfile(['../2_swan/logfiles/' fnames{i}   '.pdf'],           [fnames{i},'/section 2.1- swan figure.pdf']);
    copyfile(['../2_swan/logfiles/' fnames{i}   '_1swanprint.pdf'],[fnames{i},'/section 2.2- swan print.pdf']);
    copyfile(['../2_swan/logfiles/' fnames{i}   '_3swanout.pdf'],  [fnames{i},'/section 2.3- swan output.pdf']);
    copyfile(['../3_whafis/logfiles/' fnames{i} '_log.pdf'],       [fnames{i},'/section 3.0- whafis.pdf']);
    copyfile(['../3_whafis/' fnames{i} '.pdf'],           [fnames{i},'/section 3.1- whafis figure.pdf']);
    copyfile(['../3_whafis/logfiles/' fnames{i} '_output.pdf'],    [fnames{i},'/section 3.3- whafis output.pdf']);
    copyfile(['../6_plot/' fnames{i} 'CrestElev1.pdf'],            [fnames{i},'/section 3.4- whafis crest elevation.pdf']);
    % copyfile(['../4_taw/logfiles/' fnames{i}    '_log.pdf'],       [fnames{i},'/section 4.0- TAW log.pdf']);
    copyfile(['../4_taw/logfiles/' fnames{i}    '-runup.pdf'],     [fnames{i},'/section 4.1- TAW figure.pdf']);
    copyfile(['../4_taw/logfiles/' fnames{i}    '-DIARY.pdf'],     [fnames{i},'/section 4.2- TAW Diary.pdf']);
    copyfile(['../5_runup2/logs/' fnames{i}     '_log.pdf'],       [fnames{i},'/section 5.0- Runup2 log.pdf']);
    copyfile(['../5_runup2/logs/' fnames{i}     '_in.pdf'],        [fnames{i},'/section 5.1- Runup2 input.pdf']);
    copyfile(['../5_runup2/logs/' fnames{i}     '_out.pdf'],       [fnames{i},'/section 5.2- Runup2 output.pdf']);
    copyfile(['../5_runup2/logs/' fnames{i}     '.pdf'],           [fnames{i},'/section 5.3- Runup2 figure.pdf']);
    
end


