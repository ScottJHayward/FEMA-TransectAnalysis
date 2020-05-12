% chk nld 20181016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datafile='../data/transectdata.xls'; 
%datafile='../data/transectdata.xls';
%tDIR='../ADCIRC_returns/'; %location of transects
%imgfile='tpng'; runupname='YK-runup';
%picktoetop=0; %to pick toetop for all stations, set equal to 1
%L_append=2;
%engineer='sjh';
%jobstr='job 2';
% filename='YK_R2_T.in';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[num,txt,raw]=xlsread(datafile);
%for i=2:size(raw,1)
%Tnames{i-1}=raw{i,1};
%client{i-1}='FEMA';
%eng{i-1}=engineer;
%job{i-1}=jobstr;
%run{i-1}=1;
%desc{i-1}=['RUNUP2 transect: ' Tnames{i-1}];
%filename{i-1}=[Tnames{i-1} '.in'];
%rough{i-1}=num(i-1,15);
%toe{i-1}=num(i-1,1); %using WHAFIS start rather than TAW
%top{i-1}=num(i-1,8);
%wave{i-1}=num(i-1,3);%*0.626; % Multiply Hs by 0.626 to get Hmean
%twl{i-1}=num(i-1,2); %This is the twl 
%per{i-1}=num(i-1,4);%*0.85; % convert significant wave period from SWAN to mean period
%end

%read data from each transect
%for i=1:size(Tnames,2)
%     fname=[tDIR fnames{i} 'XYZSTA_RETURNS.csv'];
%    fname=[tDIR Tnames{i} 'XYZSTA_RETURNS.csv'];
%    file=xlsread(fname);
%    lon{i}=file(:,1);lat{i}=file(:,2);sta{i}=file(:,4);ele{i}=file(:,3);
%    sta{i}(isnan(ele{i}))=[];
%    ele{i}(isnan(ele{i}))=[];
%    lon{i}(isnan(ele{i}))=[];
%    lat{i}(isnan(ele{i}))=[];
%end
%Create waterlevel input variable
%for i=1:length(wave)
%    wlvl{i}=[twl{i} wave{i} per{i}];
%end

%% Open DOSBOX, mount drive
    % for example, this current directory may be E:/Projects/FEMA_appeal/Kittery/Runup2
    % In Dosbox, type: 
    % mount E E:/Projects/FEMA_appeal/Kittery/Runup2
    % will mount this folder to the E drive in Dosbox. 
        
% Next, Run the RUNUP2 Format script for each transect, create a different file for each one

%for i=1:size(Tnames,2)
%#    [output{i},ele_fill(i),tslope{i},xsct(i),indices{i},offsets(i),Hm_deshoal(i),strt(i),sta_shift(i)]=runup2_format(client(i),...
%#        eng(i),job(i),desc(i),run(i),sta(i),ele(i),rough(i),wlvl(i),...
%#        filename{i},toe(i),top(i),Tnames(i));
%#end

% [output,ele_fill,tslope,xsct,indices,offsets]=runup2_format(client,...
%         eng,job,desc,run,sta,ele,rough,wlvl,...
%         filename,toe,top,Tnames);

%% Read runup results and calculate 2% runup, write logs and calculate the ACES runup

% Create filenames
fdir='RUNUP2_files\';
for i=2:size(raw,1)
    filename2{i-1}=[fdir Tnames{i-1} '.out'];
end

% Open output files, read the data, save it to the csv
for i=1:length(filename2)
    
    fid=fopen(filename2{i},'r');
    
    % Register li_count to the first line of the output table
    pre_li=0;
    while ~feof(fid)
        li=fgetl(fid);
        
        if length(li)==56 && strcmp(li(45:47),'OUT')
            li_count=1;
            break
        end
        pre_li=pre_li+1;
    end
    
    % Run through the output table and collect results
    r_li=0;
    messages=0;
    while ~feof(fid)
        
        li=fgetl(fid);
        A=str2num(li);
        li_count=li_count+1;
        if  isempty(A) && ~isempty(li) && li_count>11
            messages=1;
            continue
        end
        if isempty(A)
%         if length(li) < 2
            continue
        end
        % Add in if statement for when the solution does not converge-
        % to flag it as problematic
        r_li=r_li+1;
        swel(r_li,1)=A(1);
        Hm_output(r_li,1)=A(2);
        per2(r_li,1)=A(3);
%         break_slope(r_li)=A(4);
%         runup_slope(r_li)=A(5);
        runup(r_li,1)=A(6);
        
    end
    fclose(fid);
    
    % Calculate stats on runup values
    exist runup 'var';
    if ans
        runupflag=1;
    else
        runupflag=0;
    end
    
    if runupflag %Runup is valid and computed properly
        swel_raw{i,1}=swel;
        Hm_raw{i,1}=Hm_output;
        per_raw{i,1}=per2;
        R_raw{i,1}=runup;
        R_m(i,1)=mean(runup);
        R_2per(i,1)=R_m(i)*2.2; % calculate 2% runup value
        swel_toe(i,1)=swel(1);
        R_total(i,1)=swel(1)+R_2per(i); %calculate total runup elevation
        R_std(i,1)=std(runup);
        Hs_toe(i,1)=wave{i}; %Significant wave height at toe
        Tp_toe(i,1)=per{i}; %Period of significant waveheight at toe
        Hm_toe(i,1)=wave{i}*0.626; %Mean waveheight at toe
        Tpm_toe(i,1)=per{i}*0.85; %Mean period at toe
%         Hm_deshoal(i,1)=Hm_output(5); %deshoaled waveheight output from RUNUP2        
    elseif ~runupflag %Runup program did not complete calculations, not valid
        swel_raw{i,1}=twl{i};
        Hm_raw{i,1}=-9999;
        per_raw{i,1}=-9999;
        R_raw{i,1}=-9999;
        R_2per(i,1)=-9999;
        R_m(i,1)=-9999;
        swel_toe(i,1)=-9999;
        R_total(i,1)=-9999;
        Hs_toe(i,1)=-9999;
        Hs_toe(i,1)=wave{i}; %Significant wave height at toe
        Tp_toe(i,1)=per{i}; %Period of significant waveheight at toe
        Hm_toe(i,1)=wave{i}*0.626; %Mean waveheight at toe
        Tpm_toe(i,1)=per{i}*0.85; %Mean period at toe
%         Hm_deshoal(i,1)=Hm_output(5); %deshoaled waveheight output from RUNUP2
    end
    
    % Check for the presence of messages, if they exist, flag them
    if messages == 1 && R_2per(i)~=-9999
        % These are transects that have valid numbers, but have warning or
        % error messages
        messageflag(i,1)=-1;
        messagelist{i}='Nonfatal Error, Check Output';
    elseif R_2per(i)==-9999 
        messageflag(i,1)=1;
        messagelist{i}='RUNUP2 Failed';
    else
        messageflag(i,1)=0;
        messagelist{i}='No Messages';
    end
    
     
    clear fid messages swel Hm_output per2 break_slope runup_slope runup
    
    %Write to Diary File
    diaryfile=sprintf('logs/%s_log.txt',Tnames{i});
    fid=fopen(diaryfile,'w');
    fprintf(fid,'\n_______________________________________________________\n');
    fprintf(fid,'PART 5: RUNUP2\n');
    fprintf(fid,'\n');
    fprintf(fid,'        for transect: %s\n',Tnames{i});
    fprintf(fid,'\nStation locations shifted by: %.2f feet from their\noriginal location to set the shoreline to\nelevation 0 for RUNUP2 input\n',offsets{i});
    fprintf(fid,'\n');
    fprintf(fid,'\n______________RUNUP2 INPUT CONVERSIONS_________________\n');
    fprintf(fid,'        for transect: %s\n',Tnames{i});
    fprintf(fid,'\nIncident significant wave height: %.2f feet\n',Hs_toe(i));
    fprintf(fid,'\nPeak wave period: %.2f seconds\n',Tp_toe(i));
    fprintf(fid,'\nMean wave height: %.2f feet\n',Hm_toe(i));
    fprintf(fid,'\nLocal Depth below SWEL: %.2f feet\n',D_shallow(i));
    fprintf(fid,'\nMean wave height deshoaled using Hunt approximation for\ncelerity assuming constant wave energy flux.\n');
    fprintf(fid,' References: R.G. Dean and R.A. Dalrymple. 2000.  Water\n\n');
    fprintf(fid,'             Wave Mechanics for Engineers and Scientists. World\n');
    fprintf(fid,'             Scientific Publishing Company, River Edge New Jersy\n');
    fprintf(fid,'\n');    
    fprintf(fid,'             USACE (1985), Direct Methods for Calculating Wavelength, CETN-1-17\n');
    fprintf(fid,'             US Army Engineer Waterways Experiment Station Coastel Engineering\n');
    fprintf(fid,'             Research Center, Vicksburg, MS\n');
    fprintf(fid,'\n');      
    fprintf(fid,'             also see Coastal Engineering Manual Part II-3\n'); 
    fprintf(fid,'             for discussion of shoaling coefficient\n');
    fprintf(fid,'\n');      

    for ih=1:length(HuntStr{i})
            fprintf(fid,HuntStr{i}{ih});
    end
    fprintf(fid,'\nDeepwater mean wave height: %.2f feet\n',Hm_deshoal(i));
    fprintf(fid,'\n______________END RUNUP2 CONVERSIONS___________________\n');
    fprintf(fid,'\n______________RUNUP2 RESULTS___________________________\n');
    fprintf(fid,'        for transect: %s\n',Tnames{i});
    fprintf(fid,'\nRUNUP2 SWEL:\n');
    fprintf(fid,'%.2f\n',swel_raw{i});
    fprintf(fid,'\nRUNUP2 deepwater mean wave heights:\n');
    fprintf(fid,'%.2f\n',Hm_raw{i});
    fprintf(fid,'\nRUNUP2 mean wave periods:\n');
    fprintf(fid,'%.2f\n',per_raw{i});
    fprintf(fid,'\nRUNUP2 runup above SWEL:\n');
    fprintf(fid,'%.2f\n',R_raw{i});
    fprintf(fid,'\nRUNUP2 Mean runup height above SWEL: %.2f feet\n',R_m(i));
    fprintf(fid,'\nRUNUP2 2-percent runup height above SWEL: %.2f feet\n',R_2per(i));    
    fprintf(fid,'\nRUNUP2 2-percent runup elevation: %.2f feet-NAVD88\n',R_total(i));    
    fprintf(fid,'\nRUNUP2 Messages:\n');
    fprintf(fid,'%s\n',messagelist{i});
    fprintf(fid,'\n______________END RUNUP2 RESULTS_______________________\n\n\n');



    fprintf(fid,'\n________________ACES BEACH RUNUP_______________________\n');
    fprintf(fid,'\nIncident significant wave height: %.2f feet\n',Hs_toe(i));
    fprintf(fid,'\nSignificant wave height is mean wave height divided by 0.626\n');
    fprintf(fid,'Reference: D.2.8.1.2.1 Atlanic and Gulf of Mexico G&S Feb. 2007\n');
    % calculate deepwater Hs


    fprintf(fid,'\nDeepwater significant wave height: %.2f feet\n',Hm_deshoal(i)./0.626); %multiply mean deshoaled by 1/0.626 to get Hs deshoaled
    fprintf(fid,'\nPeak wave period: %.2f seconds\n',Tp_toe(i));
    ACES_slope=(ele_fill{i}(indices{i}(1))-ele_fill{i}(indices{i}(end)))/(sta_shift{i}(indices{i}(1))-sta_shift{i}(indices{i}(end))); %this is the rise/run
    fprintf(fid,'\nAverage beach Slope: 1:%.2f (H:V) \n',1/ACES_slope);
    %calculate runup using aces method
    [AR,ARlog]=Aces_Beach_Runup(Hm_deshoal(i)./0.626,Tp_toe(i),ACES_slope); %use Significant Wave height for ACES beach runup
    fprintf(fid,'\n\n%s',ARlog);
    fprintf(fid,'\nACES RUNUP CALCULATED USING ''Aces_Beach_Runup.m''\n');
    fprintf(fid,'\nACES Beach 2-percent runup height above SWEL: %.2f feet\n',AR(2));    
    fprintf(fid,'\nACES Beach 2-percent runup elevation: %.2f feet-NAVD88\n',AR(2)+swel_raw{i}(1));    
    if ACES_slope <= 1/5
            fprintf(fid,'\nACES BEACH RUNUP is valid\n');
    else 
            fprintf(fid,'\n!!!ACES BEACH RUNUP is NOT valid\n');
    end
    fprintf(fid,'\n_____________END ACES BEACH RESULTS____________________\n');
    fprintf(fid,'\nPART 5 COMPLETE________________________________________\n');
    fclose(fid);
    
    
    
    
    % Plot profiles with Runup values
    
    figure('position', [100 100 900 900])
    hold on
    plot(sta_shift{i}(indices{i}),ele_fill{i}(indices{i}),'linewidth',3,'color',[0.5 0.2 0.1]);
    if swel_toe(i)>0
        title(['Runup2 2% runup elevation for Transect: ' Tnames{i}])
        plot([min(sta_shift{i}(indices{i})) max(sta_shift{i}(indices{i}))],[swel_toe(i) swel_toe(i)],'color','b','linewidth',2)
        plot([min(sta_shift{i}(indices{i})) max(sta_shift{i}(indices{i}))],[R_total(i) R_total(i)],'color','r','linewidth',2)
    else
        title('Runup2 error, see log sheet')
    end
    %plot(sta,dep,'linewidth',2,'color',[0.2 0.2 0.4]);
    grid on
    yylm=get (gca,'ylim');
    ylim([yylm(1) yylm(2)+5])
    dztxt=(yylm(end)-yylm(1))/20;
    str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',R_total(i));    
    text(mean([min(sta_shift{i}(indices{i})) max(sta_shift{i}(indices{i}))]),R_total(i)+dztxt,str);
    str=sprintf('Run-up Height is %5.1f Feet',R_2per(i));    
    text(mean([min(sta_shift{i}(indices{i})) max(sta_shift{i}(indices{i}))]),R_total(i)-dztxt,str);
    xlabel('Station (ft)');ylabel('Elevation (ft)')
    imgname=['logs/' Tnames{i}];
    legend('Ground','Swel','2% Runup','location','best')
    print('-r450','-dpdf','-fillpage',imgname);

%     %plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%     %scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
%     plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
%     scatter(topStaAll,R2_all+SWEL,30,'r*');
%     %scatter(toe_sta,Db,70,'k*');
%     scatter(toe_sta,Ztoe,70,'k*');
%     % plot berm segs
%     for nn=1:length(Berm_Segs)
%         seg=Berm_Segs(nn);
%         h=Berm_Heights(nn);
%         plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
%     end
%     xlim([toe_sta-50 max(topStaAll)+50])
%     grid on
%     hleg=legend('Ground','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
%     set (hleg,'location','southeast');
%     xlabel('Distance from zero NAVD88 shoreline (feet)');
%     ylabel('Elevation (Feet-NAVD88)');
%     str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
%     yylm=get (gca,'ylim');
%     dztxt=(yylm(end)-yylm(1))/15;
%     text(50,Z2+dztxt,str);
%     str=sprintf('Run-up Height is %5.1f Feet',R2_new);
%     text(60,Z2-dztxt,str);
%     title (plotTitle);
%     set(gcf,'color','w');
%     print('-r450','-dpdf',imgname);
end
% sprintf('%.2f\t%.2f\t%.2f\t%.2f\n',swel_raw{i}, Hm_raw{i}, per_raw{i}, R_raw{i})
% write the results to the existing transect data file
xlswrite('../data/transectdata.xls',R_total,1,'T2');
xlswrite('../data/transectdata.xls',messageflag,1,'V2');
xlswrite('../data/transectdata.xls',{'RUNUP2 Message Flag'},1,'V1');
    
