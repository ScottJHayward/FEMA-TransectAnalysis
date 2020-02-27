%USER INPUT
%
%Author: Scott Hayward
%Company: Ransom Consulting, inc.
%Project: 2018 FEMA appeal/2020 FEMA LOMR's, York and Cumberland Counties
%
%when this script is executed, MATLAB will prompt the user to choose the
%starting locations of the 1-D wave model, followed by the toe/top of the
%transect to be used in either iterative TAW or Runup2.
%
%If there is a defined bluff or dune
%crest, the top should be choosen as close to the peak as possible. If the
%slope levels off, the top should be the first break of the slope.
%
%Results from this script are written into a excel spreadsheet, and the
%rest of the transect-based calculations will use this sheet to read
%inputs, and publish results.
%
% chk nld 20181003


clc;clear all;close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tDIR='../ADCIRC_returns/'; %location of transects
tname='CM-';
tnums = {'144-1','158-1','158-2'}
pickSTA=ones(size(tnums)); %by default. can be array of 0's and 1's [0 1 1]
savefile='../data/transectdata.xls'; imgfile='tpng';
logpre=['logfiles/' tname];
L_append=2; %line to append to... this will overwrite previous data
map_matlab_path='C:\Users\shayward\Documents\MATLAB\Map-Matlab-master';
% map_matlab_path='C:\Users\nathan.dill\Documents\MatLab\Map-Matlab';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end config %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(map_matlab_path)
%%
for i=1:length(tnums)
    close all
    figure(1);hold off
    %open log file. Next scripts will read this file, and add data to it
    logfile=[logpre tnums{i} '_log.txt'];
    fid=fopen(logfile,'w');
    fprintf(fid,'%s\n','_______________________________________________________');
    fprintf(fid,'%s\n',['DATA LOG FOR TRANSECT ID: ' tname tnums{i}]);
    fprintf(fid,'%s\n','_______________________________________________________');
    fprintf(fid,'%s\n','');fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n','_______________________________________________________');
    fprintf(fid,'%s\n','PART 1: USER INPUT');
    fprintf(fid,'%s\n','');
    
    %read csv file and break up variables
    fname=[tDIR tname tnums{i} 'XYZSTA_RETURNS.csv'];
    file=xlsread(fname);
    lon{i}=file(:,1);lat{i}=file(:,2);sta=file(:,4);z_fema=file(:,3);z_adc=file(:,5);hs=file(:,6);wlev=file(:,7);per=file(:,8);heading=file(:,9),perheading=file(:,10);
    hs(hs<0)=NaN;per(per<0)=NaN;wlev(wlev<0)=NaN;
    %loop through the transects, extend the wave height and water
    %elevations, this removes the '-999' filler for missing data
    %     for j=2:length(sta)
    %             if wlev(j)<=0;
    %                 wlev(j)=wlev(j-1);
    %             end
    %             if hs(j)<=0;
    %                 hs(j)=hs(j-1);
    %             end
    %         end
    
    %replace -999 with NaN
    
    %only for each station specified, select toe/top
        %plot data, get user feedback
        plot(sta,z_fema,'k-'); hold on;
        plot(sta,z_adc,'k--');
        plot(sta,wlev,'b-');
        plot(sta,hs,'g');
        plot(sta,per,'r');
        ylim([-10 50])
        %plot the elevation of 2*Hs below the TWL.
        cp=min(find(max(wlev)<=z_fema)); %cross point, or shoreline
        first_cross=sta(cp);
        %         dep=max(wlev)-z_fema; dep(dep<0)=NaN; bwh=dep*.78;
        grid minor
        title(['transect' tname tnums{i} ': Pick SWAN-1D/WHAFIS start']);
        xlabel('Cross-shore Distance (feet)');
        ylabel('Elevation (ft), Significant Wave Height (ft), Period (sec)');
        legend('Transect Profile','ADCIRC Profile','TWL','HS','TP','location','best');
        if pickSTA(i)==1
        [x,z]=ginput(1);
        else 
        data=xlsread(savefile)
        x=data(i,1);
        end
        loc=find(abs(sta-x)==min(abs(sta-x)));loc=loc(1)
        wha_start(i)=sta(loc);
        
        
        fprintf(fid,'%s\n','');
        wha_wl(i)=wlev(loc);
        wha_hs(i)=hs(loc);
        wha_per(i)=per(loc);
        
        
        
        title(['transect' tname tnums{i} ': Pick Toe/Top Stations for TAW and Runup2'])
        axis([first_cross-500 first_cross+80 -20 40])
        
        if pickSTA(i)==1
        [x,z]=ginput(2);
        else 
        data=xlsread(savefile)
        x=[data(i,6) data(i,8)];
        end
        
       
        %assign the user input to toe/top variables
        name{i}=[tname tnums{i}];
        toex(i)=x(1);
        topx(i)=x(2);
        
        %interpolate the water elevation and wave height
        toex(i)=interp1(sta,sta,x(1),'nearest');
        toez(i)=interp1(sta,z_fema,x(1),'nearest');
        topx(i)=interp1(sta,sta,x(2),'nearest');
        topz(i)=interp1(sta,z_fema,x(2),'nearest');
        wl_toe(i)=interp1(sta,wlev,toex(i),'nearest');
        hs_toe(i)=interp1(sta,hs,toex(i),'nearest');
        tp_toe(i)=wha_per(i); %use SWAN/WHAFIS input
        
        
        
        fprintf(fid,'%s\n','SWAN 1-D / WHAFIS input');
        fprintf(fid,'%s\n','_______________________');
        fprintf(fid,'%s\n',['station:    ' num2str(sta(loc)) ' ft']);
        fprintf(fid,'%s\n',['LON:        ' num2str(lon{i}(loc),6) ' deg E']);
        fprintf(fid,'%s\n',['LAT:        ' num2str(lat{i}(loc),6) ' deg N']);
        fprintf(fid,'%s\n',['Bottom ELEV:       ' num2str(z_fema(loc)) ' ft-NAVD88']);
        fprintf(fid,'%s\n',['TWL:        ' num2str(wlev(loc)) ' ft-NAVD88']);
        fprintf(fid,'%s\n',['HS:         ' num2str(hs(loc)) ' ft']);
        fprintf(fid,'%s\n',['TP:         ' num2str(per(loc)) ' sec']);
        fprintf(fid,'%s\n',['Wave Direction bin: ' num2str(perheading(1)) ' deg CCW from East (90 deg sector)']);
        fprintf(fid,'%s\n\n\n',['Transect Direction:   ' num2str(heading(1)) ' deg CCW from East']);
        
        %re-set axes, plot all selected points, save a fig
        title(['transect' tname tnums{i} ': User-Selected inputs'])
        axis([wha_start(i)-100 topx(i)+500 -20 40])
        plot([wha_start(i) wha_start(i)],[-100 100],'c')
        scatter([toex(i) topx(i)],[toez(i) topz(i)],'md')
        legend('Transect Profile','ADCIRC Profile','TWL','HS','TP','1-D Wave start','TAW toe/top','location','best')
        set(gcf,'position',[100 100 900 600],'paperorientation','landscape')
        print('-r450','-dpdf',[logpre tnums{i} '_data']);
        pause(.2)
        
        
        figure(2)
        titlename=['Transect Number: ' tname tnums{i}];
        subplot(2,1,1)
        hold off
        plot([lon{i}(1) lon{i}(end)],[lat{i}(1) lat{i}(end)],'r','linewidth',1.25);hold on
        title(titlename)
        xlabel('Longitude');ylabel('Latitude');
        set(gca,'tickdir','out');
        set(gcf,'position',[100 100 800 1000])
        Map(gca);
        axis([(lon{i}(1)+lon{i}(end))/2-.03 (lon{i}(1)+lon{i}(end))/2+.03  (lat{i}(1)+lat{i}(end))/2-.02 (lat{i}(1)+lat{i}(end))/2+.02]);pause(.1)
        subplot(2,1,2)
        hold off
        plot(sta,z_fema,'k','linewidth',1.25);
        xlabel('Station (ft)')
        ylabel('Elevation (ft)')
        grid minor;
        title('Elevation profile')
        set(gca,'tickdir','out');
        %print the title page for the transect
        print('-r450','-dpdf',[logpre tnums{i} '_title']);
        
        station{i}=sta;
        elevation{i}=z_fema;
        adcircelevation{i}=z_adc;
        
        fprintf(fid,'%s\n','TAW/RUNUP input');
        fprintf(fid,'%s\n','_______________');
        fprintf(fid,'%s\n',['toe sta:     ' num2str(toex(i)) ' ft']);
        fprintf(fid,'%s\n',['toe elev:    ' num2str(toez(i)) ' ft-NAVD88']);
        % fprintf(fid,'%s\n',['toe TWL:     ' num2str(wl_toe(i))]);
        % fprintf(fid,'%s\n',['toe HS:      ' num2str(hs_toe(i))]);
        % fprintf(fid,'%s\n',['toe TP:      ' num2str(tp_toe(i))]);
        fprintf(fid,'%s\n',['top sta:     ' num2str(topx(i)) ' ft']);
        fprintf(fid,'%s\n',['top elev:    ' num2str(topz(i)) ' ft-NAVD88']);
        fprintf(fid,'%s\n','*Wave and water level conditions at toe to be calculated in SWAN 1-D*');
        
        fprintf(fid,'%s\n','');
        
        fprintf(fid,'%s\n','PART 1 COMPLETE________________________________________')
        fprintf(fid,'%s\n','')
        
        hold off
        fclose all
end
% Now, populate the XLS data summary file. this will be used by other scripts.
%%
slope=(topz-toez)./(topx-toex);

%header for transect data csv
d={'Transect Name','WHAFIS/SWAN Start station (ft)','WHAFIS/SWAN TWL (ft)','WHAFIS/SWAN Hs (ft)','WHAFIS/SWAN Period (sec)','WHAFIS/SWAN Max Setup (ft)','Runup toe station (ft)','Runup toe elevation (ft)','Runup top station (ft)','Runup top elevation (ft)','Slope','Setup at toe (ft)','Hs at toe (ft)','Per at toe (ft)','GBERM','GROUGH','GBETA','GPERM','TAW elev (ft)','Runup2 elev (ft)','TAW valid?'};
xlswrite(savefile,d,'Sheet1','A1');

for i=1:length(tnums)
    if pickSTA(i)==1
        %write the transect data into excel file
        d={name{i},wha_start(i),wha_wl(i),wha_hs(i),wha_per(i),'TBD',toex(i),toez(i),topx(i),topz(i),slope(i),'TBD','TBD','TBD',1,1,1,1,'TBD','TBD',1,1};%re-read gamma if pickstta=0
        xlswrite(savefile,d,'Sheet1',['A' num2str(L_append+i-1)]);
    else    %write the transect data into excel file
        data=xlsread(savefile);
        gberm=data(i,14);
        grough=data(i,15);
        gbeta=data(i,16);
        gperm=data(i,17);
        d={name{i},wha_start(i),wha_wl(i),wha_hs(i),wha_per(i),'TBD',toex(i),toez(i),topx(i),topz(i),slope(i),'TBD','TBD','TBD',gberm,grough,gbeta,gperm};%re-read gamma if pickstta=0
        xlswrite(savefile,d,'Sheet1',['A' num2str(L_append+i-1)]);
    end
end

% data.sta=station;
% data.ele=elevation;
% data.adc_ele=adcircelevation;
% data.name=name;
% data.toex=toex;
% data.toez=toez;
% data.topx=topx;
% data.topz=topz;
% data.wl=wl_toe;
% data.hs=hs_toe;
% data.per=tp_toe;
% data.lon=lon;
% data.lat=lat;
% data.wha_start=wha_start;
% save(savefile,'data')
