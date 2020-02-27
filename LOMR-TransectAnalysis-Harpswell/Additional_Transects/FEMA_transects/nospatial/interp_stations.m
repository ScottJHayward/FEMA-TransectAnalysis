%read start and end locations of each transect, and interpolate 1-ft
clear all;close all;fclose all
format long

DATA=xlsread('start_end.csv');
fpre='../reconstructed/rc_';
transects={'CM-144-1','CM-158-1','CM-158-2'};

for i=1:length(transects)
    
    %read lon/lat data
    start_x=DATA(i,1);
    start_y=DATA(i,2);
    end_x=DATA(i,3);
    end_y=DATA(i,4);
    
    lon=[start_x,end_x];
    lat=[start_y,end_y];
    
    [xutm,yutm,utmzone]=deg2utm(lat,lon);
    
    len=sqrt((xutm(1)-xutm(2))^2+(yutm(1)-yutm(2))^2);
    pts=(1:1/3.28083:len);
    x=interp1([0;len],xutm,pts);
    y=interp1([0;len],yutm,pts);    
    
%     %extend offshore by 1000 
%     dx=x(2)-x(1);
%     dy=y(2)-y(1);
%     for n=1:1000
%         x=[x(1)-dx x];
%         y=[y(1)-dy y];
%     end
    %loop through each index and convert to UTM
    for j=1:length(x)
        [y2(j),x2(j)]=utm2deg(x(j),y(j),utmzone(1,:));
    end
    %
    tx{i}=x2;
    ty{i}=y2;
    sta{i}=1:length(tx{i});

    fname=[fpre transects{i} 'XY.csv'];
    d2=[tx{i}' ty{i}' sta{i}']; 
    dlmwrite(fname,d2,'delimiter',',','precision',10);
end