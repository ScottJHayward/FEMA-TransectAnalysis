clc;clear all;close all
format long
%%%%%%%%%%%%%%%%%%   config   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datafile='../data/transectdata.xls'; 
tDIR='../ADCIRC_returns/'; %location of transects
fid=fopen('runup.kml','w'); %kml file name
prefix='CM-';
%%%%%%%%%%%%%%%%%%   end config   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,txt,raw]=xlsread(datafile);
for i=2:size(raw,1)
fnames{i-1}=raw{i,1};
end
setup=num(:,11);
twl=num(:,2)+setup;
hs=num(:,12);
per=num(:,13);
toex=num(:,6);
toez=num(:,7);
topx=num(:,8);
topz=num(:,9);
runup=num(:,18);
r2=num(:,19);
valid=num(:,20);%Look for TAW valid logic in csv, select RUNUP2 if not valid

for i=1:length(fnames)

fname=[tDIR fnames{i} 'XYZSTA_RETURNS.csv'];
file=xlsread(fname);
lon{i}=file(:,1);lat{i}=file(:,2);sta{i}=file(:,4);ele{i}=file(:,3);

end

kmlHEAD(fid)

for i=1:length(fnames)
        close all
    figure
    %determine runup elevation. 
    %find closest station to top of structure/slope
    f=find(abs(sta{i}-ceil(topx(i)))==min(abs(sta{i}-ceil(topx(i)))));
    top=sta{i}(f);
    f=f(end);
    
    %plot elevation profile
    plot(sta{i},ele{i});hold on
    plot([toex(i) topx(i)],[toez(i) topz(i)],'k')
    scatter(sta{i}(f),ele{i}(f),'g');
    
    %is taw valid? if so, use TAW. If not, use the next column over.
    %(sometimes runup2, sometimes another method)
    if valid(i)==0;
        RUP(i)=r2(i);
    else
        RUP(i)=runup(i);
    end
    
    %now move ve zone 30 ft landward of structure, collect location and
    %Lat/Lon of point, for plotting later.
    if RUP(i)>=topz(i)
        RUP_sta(i)=topx(i)+30;
        RUP_lon(i)=interp1(sta{i},lon{i},RUP_sta(i));
        RUP_lat(i)=interp1(sta{i},lat{i},RUP_sta(i));        
    else
        x=find(abs(sta{i}-ceil(toex(i)))==min(abs(sta{i}-ceil(toex(i)))));
        x=x(end);
        while RUP(i) > ele{i}(x)
            x=x+1;
        end
        RUP_sta(i)=interp1([ele{i}(x-1) ele{i}(x)],[sta{i}(x-1) sta{i}(x)],RUP(i));
        RUP_lon(i)=interp1(sta{i},lon{i},RUP_sta(i));
        RUP_lat(i)=interp1(sta{i},lat{i},RUP_sta(i));
    end
    
    %write placemark for runup
    kmlPMARKred(fid,RUP_lon(i),RUP_lat(i),['RUP: ' num2str(RUP(i))]);
    %write placemark for the selected top
    CREST_lon(i)=interp1(sta{i},lon{i},topx(i));
    CREST_lat(i)=interp1(sta{i},lat{i},topx(i)); 
    kmlPMARKred(fid,CREST_lon(i),CREST_lat(i),['TOP: ' num2str(topz(i))]);

    
    
    %plot the runup
    scatter(RUP_sta(i),RUP(i),'r');
    pause(.5)
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WHAFIS DATA %%%%%%%%%%%%%%%%%%%%%%%%%
% %   whafisCSV=[WhafisDir 'HW-' num2str(i) '_parsed.csv'];
%     whafisCSV=[prefix num2str(i) '_parsed.csv'];
%         
%     %read csv, plot wave heights 
%     data=xlsread(whafisCSV);
%     WHA.sta=data(:,1);
%     WHA.wl=data(:,6);
%     WHA.crest=data(:,9);
%     WHA.wave=data(:,8);
%     WHA.zone=zeros(size(WHA.sta));
%     
%     
%     plot(WHA.sta,WHA.crest,'b');
%     
%     
%     %limits for graph
%     xlim([toex(i)-50 topx(i)+50])
%     ylim([toez(i)-10 topz(i)+10])
%     %determine cards. 
%     WHA.zone(find(WHA.wave>0))=1;
%     WHA.zone(find(WHA.wave>=3))=2;
%     ele=round(WHA.crest);
%    
%     WHA.cindex={};
%     WHA.cnames={};
%     
%     for j=1:length(WHA.zone)-1
% %         if ele(j) ~= ele(j-1) || ele(j) ~= ele(j+1) || WHA.zone(j) ~= WHA.zone(j-1) || WHA.zone(j) ~= WHA.zone(j+1)
%         if ele(j) ~= ele(j+1) || WHA.zone(j) ~= WHA.zone(j+1)
%             %this might need to be changed. do we use ceiling, or round
%             %to nearest foot? i am assuming round is correct. 
%             ele2=num2str(ele(j));
%             WHA.cindex{end+1}=j;
%             if WHA.zone(j) == 2
%                 WHA.cnames{end+1}=['VE' ele2];
%             elseif WHA.zone(j) == 1
%                 WHA.cnames{end+1}=['AE' ele2];
%             else
%                 WHA.cnames{end+1}=['X'];
%             end
%         end
%     end
%     WHA.cnames{end+1}=['X'];
%     if length(WHA.cindex) >= 1
%         WHA.cindex{end+1}=WHA.cindex{end}+1;
%     else
%         WHA.cindex{1}=1;
%     end
%     
%     %now create a kml pin and label for each WHAFIS card
%     WHstart=find(sta{i}==WHA.sta(1))-1;
%     for n=1:length(WHA.cindex)
%         %go back and fix long, and lat. they should be cells of arrays. 
%         index(n)=find(WHA.sta(n)==sta{i});
%         kmlPMARK(fid,lon{i}(index(n)),lat{i}(index(n)),WHA.cnames{n});
%     end
    pause(.5)
end


kmlFOOT(fid);
fclose all

function kmlHEAD(fid)
fprintf(fid,'%s\n','<?xml version="1.0" encoding="utf-8"?> ');
fprintf(fid,'%s\n',' <kml>  ');
fprintf(fid,'%s\n','  <Document>  ');
fprintf(fid,'%s\n','   <name>WHAFIS and RUNUP calculations</name>  ');
end

function kmlPMARK(fid,x,y,name)
fprintf(fid,'%s\n','<Placemark>  ');
fprintf(fid,'%s\n',['    <name>' name '</name>  ']);
fprintf(fid,'%s\n','    <Point>  ');
fprintf(fid,'%s\n',['     <coordinates>' num2str(x,'%3.8f') ',' num2str(y,'%3.8f') '</coordinates>  ']);
fprintf(fid,'%s\n','    </Point>  ');
fprintf(fid,'%s\n','   </Placemark>  ');
end

function kmlPMARKred(fid,x,y,name)
fprintf(fid,'%s\n','<Placemark>  ');
fprintf(fid,'%s\n',['    <name>' name '</name>  ']);
fprintf(fid,'%s\n','    <styleUrl>#msn_R</styleUrl>  ');
fprintf(fid,'%s\n','    <Point>  ');
fprintf(fid,'%s\n',['     <coordinates>' num2str(x,'%3.8f') ',' num2str(y,'%3.8f') '</coordinates>  ']);
fprintf(fid,'%s\n','    </Point>  ');
fprintf(fid,'%s\n','   </Placemark>  ');
end

function kmlFOOT(fid)
fprintf(fid,'%s\n','  </Document>');  
fprintf(fid,'%s\n',' </kml>  ');

end
