%write log file, determine whether or not to add information from the
datafile='../data/transectdata.xls';

[num,txt,raw]=xlsread(datafile);

for i=2:size(raw,1)
fnames{i-1}=raw{i,1};
end

setup=num(:,11);
maxsetup=num(:,5);
twl=num(:,2);
hs=num(:,12);
per=num(:,13);
toex=num(:,6);
toez=num(:,7);
topx=num(:,8);
topz=num(:,9);
gberm=num(:,14);
grough=num(:,15);
gbeta=num(:,16);
gperm=num(:,17);
runup=num(:,18);
valid=num(:,20);

for i=4:length(fnames)
    %create log file
    fid=fopen(['logfiles/' fnames{i} '_log.txt'] , 'w')
    
    %write header for file, 
    fprintf(fid,'%s\n','_______________________________________________________');
    fprintf(fid,'%s\n','PART 4: TAW');
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n',  ['Input Paramters:']);
    fprintf(fid,'%s\n',  ['    TWL- ' num2str(twl(i)) ' feet' ]);
    fprintf(fid,'%s\n',  ['    HS-  ' num2str(hs(i))  ' feet' ]);
    fprintf(fid,'%s\n',  ['    PER- ' num2str(per(i)) ' seconds']);
    fprintf(fid,'%s\n',  ['    TOE- x: ' num2str(toex(i)) ' , z: ' num2str(toez(i)) ' feet']);
    fprintf(fid,'%s\n',  ['    TOP- x: ' num2str(topx(i)) ' , z: ' num2str(topz(i)) ' feet']);
    fprintf(fid,'%s\n',  ['    GBERM-   ' num2str(gberm(i))]);
    fprintf(fid,'%s\n',  ['    GGROUGH- ' num2str(grough(i))]);
    fprintf(fid,'%s\n',  ['    GBETA-   ' num2str(gbeta(i))]);
    fprintf(fid,'%s\n',  ['    GPERM-   ' num2str(gperm(i))]);
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n','RUNNING TAW:');
    fprintf(fid,'%s\n','...');
    fprintf(fid,'%s\n',['MATLAB DIARY: /4_taw/logfiles/' fnames{i} '-DIARY.txt']);
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n','CHECKING VALIDITY: ');
    fprintf(fid,'%s\n','...');
    if valid(i)==1
        fprintf(fid,'%s\n','TAW method is valid!');
        fprintf(fid,'%s\n','Using TAW runup to detemine runup elevation');
        fprintf(fid,'%s\n',['TAW 2% runup: ' num2str(runup(i)) ' feet']);
    else
        fprintf(fid,'%s\n','TAW method is not valid!');
        fprintf(fid,'%s\n','Runup elevation to be calculated using another method');
    end
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n','PART 4 COMPLETE________________________________________');
    






  
end
    
