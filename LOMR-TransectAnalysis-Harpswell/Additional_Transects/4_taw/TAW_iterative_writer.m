%chk nld 20200220

clc;clear all;close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datafile='../data/transectdata.xls'; 
tDIR='../ADCIRC_returns/'; %location of transects
imgfile='tpng'; runupname='CM-runup';
L_append=2;
csvoutpre='inpfiles/';
templatefile='TAW_template.txt';
templatelines=401;

%config
CITYNAME='The Town of Harpswell';
COUNTY='Cumberland';
ENGINEER='SJH';
DATE=date;
tawfilename='TAW_iterative.m'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,txt,raw]=xlsread(datafile);
for i=2:size(raw,1)
fnames{i-1}=raw{i,1};
end
setup=num(:,11);
setup(isnan(setup))=0;
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

% for i=[1:17 19:32 34:length(fnames)]
for i=1:size(raw,1)-1

fname=[tDIR fnames{i} 'XYZSTA_RETURNS.csv'];
file=xlsread(fname);
lon{i}=file(:,1);lat{i}=file(:,2);sta{i}=file(:,4);ele{i}=file(:,3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%add output to excel file
% for i=[1:17 19:32 34:length(fnames)]
for i=1:size(raw,1)-1
    
    %write a transect file for each run
    toesta=find(abs(sta{i}-toex(i))==min(abs(sta{i}-toex(i))));
    topsta=find(abs(sta{i}-topx(i))==min(abs(sta{i}-topx(i))));
    fid=fopen([csvoutpre,fnames{i},'sta_ele_include.csv'],'w')
    sta2=floor(sta{i}(toesta)):ceil(sta{i}(topsta));
    ele2=interp1(sta{i},ele{i},sta2);

    figure('position', [100 100 1000 600])
    plot(sta2,ele2,'kd');hold on
    plot(sta{i},ele{i},'color',[0.5 0.2 0.1])
    
    wl=twl(i)+setup(i);
    dd=wl-ele{i};
    kk=find((dd < 0), 1) ;
    shoreX=interp1([ele{i}(kk-1) ele{i}(kk)],[sta{i}(kk-1) sta{i}(kk)],wl);
    plot([min(sta{i}) shoreX],[twl(i)+setup(i) twl(i)+setup(i)],'b-','linewidth',2);


    grid on
    xlabel ('Station (ft)','fontsize',16);
    ylabel ('Elevation (ft-NAVD88)','fontsize',16);
    title ([fnames{i},'TAW Profile Selection'],'fontsize',16);
    xlim([min(sta2)-100 max(sta2)+100])
    ylim([0 40]);
    grid minor

    for s=1:length(sta2)
        fprintf(fid,'%.6f,%.6f,%d\n',[sta2(s) ele2(s) 1]);
    end
    fclose all
end

%now read the template file
fid=fopen(templatefile,'r');
for i=1:templatelines
    %WRITE INFO    
    template{i}=fgetl(fid);
end

fid2=fopen(tawfilename,'w')
% for i=[1:17 19:32 34:length(fnames)]
for i=1:size(raw,1)-1
for li=1:length(template)
    tawfile{li}=template{li}; %read each line of template
    tawfile{li}=strrep(tawfile{li},'$DIARYFILE$',['logfiles/' fnames{i} '-DIARY.txt']);
    tawfile{li}=strrep(tawfile{li},'$CITYNAME$',CITYNAME);
    tawfile{li}=strrep(tawfile{li},'$COUNTY$',COUNTY);
    tawfile{li}=strrep(tawfile{li},'$DATE$',DATE);
    tawfile{li}=strrep(tawfile{li},'$TID$',fnames{i});
    tawfile{li}=strrep(tawfile{li},'$PNGNAME$',['logfiles/' fnames{i} '-runup']);
    tawfile{li}=strrep(tawfile{li},'$ENGINEER$',ENGINEER);
    tawfile{li}=strrep(tawfile{li},'$CSVNAME$',[csvoutpre,fnames{i},'sta_ele_include.csv']);
    tawfile{li}=strrep(tawfile{li},'$CSVNAME2$',[tDIR fnames{i} 'XYZSTA_RETURNS.csv']);
    tawfile{li}=strrep(tawfile{li},'$SWEL$',num2str(twl(i)));
    tawfile{li}=strrep(tawfile{li},'$HSIG$',num2str(hs(i)));
    tawfile{li}=strrep(tawfile{li},'$TPER$',num2str(per(i)));
    tawfile{li}=strrep(tawfile{li},'$GBERM$',num2str(gberm(i)));
    tawfile{li}=strrep(tawfile{li},'$GROUGH$',num2str(grough(i)));
    tawfile{li}=strrep(tawfile{li},'$GBETA$',num2str(gbeta(i)));
    tawfile{li}=strrep(tawfile{li},'$GPERM$',num2str(gperm(i)));
    tawfile{li}=strrep(tawfile{li},'$TOESETUP$',num2str(setup(i)));
    tawfile{li}=strrep(tawfile{li},'$MAXSETUP$',num2str(maxsetup(i)));
    tawfile{li}=strrep(tawfile{li},'$PLOTTITLE$',['Iterative TAW for ' fnames{i}]);
    tawfile{li}=strrep(tawfile{li},'$TNAME$',fnames{i});
    
    
    %print to file
    fprintf(fid2,'%s\n',tawfile{li});
    
end

    fprintf(fid2,'%s','xlswrite(''../data/transectdata.xls'',Z2,''Sheet1'',''S');
    fprintf(fid2,'%s\n',[num2str(i+1) ''');']);
    
    fprintf(fid2,'%s','xlswrite(''../data/transectdata.xls'',TAW_ALWAYS_VALID,''Sheet1'',''U');
    fprintf(fid2,'%s\n',[num2str(i+1) ''');']);
    
    fprintf(fid2,'%s','xlswrite(''../data/transectdata.xls'',gamma_berm,''Sheet1'',''O');
    fprintf(fid2,'%s\n',[num2str(i+1) ''');']);

end
fclose all





