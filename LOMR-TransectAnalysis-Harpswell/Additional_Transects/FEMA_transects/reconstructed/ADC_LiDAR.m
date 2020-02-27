%this script is meant to combine the ADCIRC and LIDAR interpolated
%transects that were extracted along FEMA's transect. 
clc; clear all; close all
lidfile='./lidar/';
adcfile='./bathy/';

num={'144-1','158-1','158-2'};


for i=1:length(num)
    source={};
    %input and output file names
    lidname=[lidfile 'CM-' num{i} 'XYZSTA.csv'];
    adcname=[adcfile 'CM-' num{i} 'XYZSTA.csv'];
    csvname=['../CM-' num{i} 'XYZSTA.csv'];
       
    %read and assign lidar to variables
    [DATA1,text,raw]=xlsread(lidname);
    sta=DATA1(:,4);
    ele=DATA1(:,3);
    lon=DATA1(:,1);
    lat=DATA1(:,2);
    for n=1:length(lat)
        source{n}=raw{n,5};
    end
    
    %read Bathy data and assign to variable
    DATA2=xlsread(adcname);
    sta_adc=DATA2(:,4);
    ele_adc=DATA2(:,3);
%     hs=DATA2(:,5);
%     twl=DATA2(:,6);
       
    %combine data
    for i2=1:length(ele)
        if ele_adc(i2)<=0
            ele(i2)=ele_adc(i2);
            source{i2}='Maine ADCIRC model'; %use adcirc if below zero
        end
    end
        ele(ele<=-999)=NaN;

    %downsample data
    n1=1;
    sta2=[];
    ele2=[];
    lon2=[];
    lat2=[];
    hs2=[];
    twl2=[];
    source2={};
    while n1 <= length(sta)-1
        if isnan(ele(n1))
        else
        n2=n1;
        j=1;
        while ele(n1+1)==ele(n1) && n1<=length(sta)-2
            n2(j)=n1;
            n1=n1+1;
            j=j+1;
        end
        %now combine indices
            sta2(end+1)=mean(sta(n2));
            ele2(end+1)=mean(ele(n2));
            lon2(end+1)=mean(lon(n2));
            lat2(end+1)=mean(lat(n2));
            source2{end+1}=source{n2(1)};
%             hs2(end+1)=mean(hs(n2));
%             twl2(end+1)=mean(twl(n2));
        end
                    n1=n1+1;
    end
    
   %shift stations to be zero at the shoreline, and round station to
   %nearest foot
   x=1;
   while ele2(x)<0
       x=x+1;
   end
   sta2=sta2-sta2(x);
   
   
    figure
    plot(sta2,ele2);pause(.1)
    
    fid=fopen(csvname,'w');
    DATA=[lon2;lat2;ele2;twl2;sta2]';
    fprintf(fid,'%s\n','X,Y,Z,STA,source');
    for l=1:length(lon2)
    fprintf(fid,'%6.10f,%6.10f,%6.10f,%6.2f,',DATA(l,:));
    fprintf(fid,'%s\n',source2{l});
    end
fclose all
end
