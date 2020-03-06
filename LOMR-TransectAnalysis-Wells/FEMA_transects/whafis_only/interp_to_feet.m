
%interpolate the station/elevation of the WHAFIS transects to one foot
%increments
format long
for i=[102 103 105 106 107 108 109 110 120 121 124 125 127]
% for i=['07'];
    

    %file names
    f_in  = ['YK-' num2str(i) 'XY.csv'];
    f_out = ['../reconstructed/rc_YK-' num2str(i) 'XY.csv'];
%     f_in  = ['YK-07XY.csv'];
%     f_out = ['../reconstructed/rc_YK-07XY.csv'];
    
    %read the input file
    DATA=xlsread(f_in);
    x=DATA(:,1);
    y=DATA(:,2);
    ele=DATA(:,3);
    sta=DATA(:,4);
    
    %sort the stations
    [sta,index]=sort(sta);
    ele=ele(index);
    y=y(index);
    x=x(index);
    
    
    s2=sta(1):1:sta(end);
    x2=interp1(sta,x,s2);
    y2=interp1(sta,y,s2);
    ele2=interp1(sta,ele,s2);
    
    DATA2=[x2;y2;s2;ele2]';
    
    dlmwrite(f_out,DATA2,'delimiter',',','precision',10);
    
end

