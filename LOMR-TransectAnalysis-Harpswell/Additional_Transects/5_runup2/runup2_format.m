function [output,ele_fill,tslope,xsct,indices,offsets,Hm_deshoal,D_shallow,strt,sta_shift,HuntStr]=runup2_format(client,eng,job,desc,run,sta,ele,rough,wlvl,filename,toe,top,Tnames)

% This function is designed to take the input variables and create an
% output file called transectname.in to be used as an uput file for the
% FEMA RUNUP 2.0 program to calculate wave runup. The output file is a
% fixed-width text file that follows the format outlined in the FEMA RUNUP
% 2.0 Readme File.

% Inputs
%   Client Name:    client  Cell, Maximum of 26 characters each
%   Engineers Name: eng     Cell, Maximum of 10 characters each
%   Job Number:     job     Cell, Maximum of 10 characters each
%   Transect Name:  desc    Cell, Maximum of 74 characters each
%   Run Number:     run     Cell, Maximum of 4 digits each
%   Last Slope:     lslope  Cell, Slope of the beach profile landward of
%                           the most landward point in the transect. If you
%                           would like to have the script calculate the
%                           last slope this value should be "0".
%   Station:        sta     Cell of the x-location along the transect in
%                           feet. The station should be an integer, and
%                           negative values should be seaward of the point
%                           at which elevation=0, while positive values
%                           should be landward of that point.
%   Elevation:      elev    Cell of the z-location (elevation) at each
%                           station location along the transect. negative
%                           values are elevations below the datum, positive
%                           values are above the datum.
%   Roughness:      rough   Cell of the surface roughness at each
%                           location along the transect line OR If you would
%                           like all the surface roughness values to be
%                           equal to the same value for a given transect, 
%                           just enter that value as a single value. The
%                           roughness value represents the roughness
%                           between the station point is corresponds to and
%                           the next landward station point.
%   Water Level:   wlvl     Cell, each cell represents an individual
%                           transect and contains a 1-d array that has the
%                           following columns: 
%                           column 1: total (still) water level (ft
%                           NAVD88) 
%                           column 2: mean
%                           waveheight at the toe of the slope (ft) 
%                           column 3: mean wave period in seconds
%                           The function takes the single input line for
%                           each transect and computes +/- 5% of the mean
%                           wave heigh and wave period for 9 total lines.
%   filename:               A string containing the name for the output
%                           file including the extension (.in)

% %% Simple Test Data
% clc
% clear all
% client={'Test1' 'Test_2'};
% eng={'Ransom' 'Ransom'};
% job={'21' '5126B'};
% desc={'French' 'purple'};
% run={1 1};
% lslope={10 10};
% xsct1=[-30 -200 1;...
%     0 0 1;...
%     8.4 544 .9;...
%     15.6 699.69 .9];
% xsct2=[-20 -100 1;...
%     0 0 1;...
%     10 50 1;...
%     10.75 56.75 1;...
%     20 85 1];
% xsct={xsct1,xsct2};
% wlvl1=[13.8 5.46 5.13;13.8 6.04 5.67];
% wlvl2=[13.8 5.46 5.13;13.8 6.04 5.67];
% wlvl={wlvl1,wlvl2};

%% Function

n=length(client); %number of transects
g=32.17404856;%gravitational constant

for i=1:n
    % Check for roughness. If a single value is given, create a vector of
    % all the same value for th length of the transect
    if length(rough{i})==1
        value=rough{i};
        rough_temp(1:length(sta{i}),1)=value;
        rough{i}=rough_temp;
        rough_temp=[];
    else
    end
    % Start Diary Files
    Tname=Tnames{i};
%     diaryfile{i}=sprintf('logs/%s_log.txt',Tname);
end

% Convert significant wave height and period to mean wave height and period
for i=1:n
wlvl{i}(2)=wlvl{i}(2)*0.626;
wlvl{i}(3)=wlvl{i}(3)*0.85;
end

% Shift the station locations to reference the point of elevation 0' as
% station 0'
for i=1:n
    ele_temp=ele{i};
    sta_temp=sta{i};
    for j=2:length(ele_temp)-1
        if ele_temp(j)<=0
%             zero_neg=ele_temp(j-1);
%             zero_pos=ele_temp(j);
        else
            break
        end
    end
    slope=(ele_temp(j)-ele_temp(j-1))/(sta_temp(j)-sta_temp(j-1));
    zero_dist=ele_temp(j-1)*slope^(-1);
    offset=sta_temp(j-1)-zero_dist;
    offsets{i}=offset;
    sta_shift{i}=sta_temp-offset;
    ele_temp=[];
    sta_temp=[];
    
%     % Write to Diary File
%     fid=fopen(diaryfile{i},'a');
%     fprintf(fid,'*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*\n');
%     fprintf(fid,'*--------- RUNUP2 Runup Calculations ---------*\n');
%     fprintf(fid,'*_____________________________________________*\n');
%     fprintf(fid,'        for transect: %s\n',Tname);
%     fprintf(fid,' Station locations shifted by: %.2f feet from their\noriginal location to set the shoreline to\nelevation 0 for RUNUP2 input\n',offset);
%     fclose(fid);
end
    
% Reduce the resolution of the transect inputs
num=20;
% disp('Select the seaward starting point for your RUNUP transect first, then select the landward ending point. Press ENTER to cotinue')
for i=1:n
    m=length(sta_shift{i});
    strt=find(sta_shift{i}>=(toe{i}-offsets{i}));strt(2:end)=[];
    ele_temp=ele{i};
    for j=strt:m-1
        if ele_temp(j+1)<ele_temp(j)
            ele_temp(j+1)=ele_temp(j);
        else
        end
    end
    ele_fill(i)={ele_temp};
    ele_temp=[];
end

for i=1:n
    st_temp=sta_shift{i};
    ele_temp=ele_fill{i};
    rough_temp=rough{i};
    swel=wlvl{i}(1);
%     plot(st{i},ele{i},'r',st_temp(st_temp<=0),ele_temp(st_temp<=0),'b',st_temp(st_temp>=0),ele_temp(st_temp>=0),'k')
%     hold on
%     [x,~]=ginput;
%     begin_pt=floor(x(1));
%     end_pt=ceil(x(2));
    begin_pt=floor(toe{i}-offsets{i});
    end_pt=ceil(top{i}-offsets{i});
    [indices]=res_reduce_tbn(st_temp,ele_temp,begin_pt,end_pt,num);
    end_ele=ele_temp(indices(end)); % Elevation of the final point in the transect
    end_ele_2=ele_temp(indices(end-1)); % Elevation of the second to last point in the transect
    end_pt_2=st_temp(indices(end-1)); % Station point of the second to last point in the transect
    lslope{i}=round((end_pt-end_pt_2)/(end_ele-end_ele_2));
    xsct{i}=([ele_temp(indices) st_temp(indices) rough_temp(indices)]);
%     ax=gca;
%     hold on
%     ax=plot(st{i},ele{i},'k',st_temp(indices),ele_temp(indices),'gd');
%     xlim([min(st_temp(indices(1))-20) max(st_temp(indices(end))+20)])
%     hold off
%     uiwait(msgbox('If the reduced resolution profile is ok, click to continue'))
%     close

    % De-shoal the wave height
    D=swel-ele_temp(indices(1));
    D_shallow=D;    

    [~,H0_H,HuntStr]=deshoal(wlvl{i}(2),wlvl{i}(3),D,g);
    wlvl{i}(2)=H0_H;
    Hm_deshoal(i)=H0_H;
    
    % Clear variable for next loop
    st_temp=[];
    ele_temp=[];
    rough_temp=[];
    begin_pt=[];
    end_pt=[];
    x=[];
    H0_H=[];
end

% Calculate the wave heights and period +/- 5%
[wlvl_out]=wlvl_range(wlvl,n);

% Begin creating fixed-width file
holder='';
% ouput=repmat(' ',[
counter=0;
startline=1;
% i=1;
for i=1:n
    % Name Line
    header(1,:)=sprintf('%2s%-26s%32s%-10s%+10s',holder,client{i},holder,eng{i},job{i});
    % Job Description Line
    header(2,:)=sprintf('%2s%-74s%4d',holder,desc{i},run{i});
    % Last Slope Line
    % If statement for when the slope is 4 characters, then make it no
    % decimal, otherwise add a decimal
    if lslope{i}>999
        tslope{i}='999.';
    elseif lslope{i}>=100
        tslope{i}=[num2str(lslope{i}) '.'];
    elseif lslope{i}<100&&lslope{i}>=10
        tslope{i}=[num2str(lslope{i}) '.0'];
    elseif lslope{i}<10
        tslope{i}=[num2str(lslope{i}) '.00'];
    else
    end        
    header(3,:)=sprintf('%s%76s',tslope{i},holder);
    counter=counter+3;
    % Profile Lines
    tempx=xsct{i};
    m=size(tempx,1);
    profile=repmat(' ',[m 80]);
    for j=1:m
        counter=counter+1;
        if j==m
            lpflag='1';
        else
            lpflag=[];
        end
        temp=sprintf('%1s %-5.2f %-6.1f %-5.1f%60s',lpflag,tempx(j,1),tempx(j,2),tempx(j,3),holder);
        profile(j,:)=temp(1:80);
%         ele_str=num2str(round(tempx(j,1),2));
%         st_str=num2str(round(tempx(j,2),2));
%         ele_num=round(tempx(j,1),1);
%         st_num=round(tempx(j,1),1);
%         profile(j,:)=sprintf('%1s %-5f %-6f %-5.1f%60s',lpflag,ele_num,st_num,tempx(j,3),holder);
%         ele_str=[];
%         st_str=[];
        ele_num=[];
        st_num=[];
    end
    % Wave Level Lines
    tempw=wlvl_out{i};
    m=size(tempw,1);
    wave=repmat(' ',[m 80]);
    for j=1:m
        counter=counter+1;
        if j==m&&i~=n
            ntflag='1';
        else
            ntflag='';
        end
        wave_pre=sprintf('%1s%-5.1f %-5.2f %-5.2f%62s',ntflag,tempw(j,1),tempw(j,2),tempw(j,3),holder);
%         fid=fopen(diaryfile{i},'a');
%         fprintf(fid,'Input Total Water Level= %5.1f\nInput Mean Wave Height= %5.2f\nInput Mean Wave Period= %5.2f\n',...
%             tempw(j,1),tempw(j,2),tempw(j,3));
%         fclose(fid);
        wave(j,:)=wave_pre(1:80);
    end
    output_temp=[header;profile;wave];
    output(startline:counter,:)=output_temp;
    startline=counter+1;
    tempx=[];
    tempw=[];
    profile=[];
    wave=[];
    output_temp=[];
end

% Print RUNUP2 Input to Diary File
% for i=1:n
%     fid=fopen(diaryfile{i},'a');
%     fprintf(fid,'*_____________________________________________*\n');
%     fprintf(fid,'*_____________RUNUP2 INPUT FILE_______________*\n');
%     for j=1:size(output,1)
%         fprintf(fid,'%s\r\n',output(j,:));
%     end
%     fprintf(fid,'*___________END RUNUP2 INPUT FILE_____________*\n');
%     fclose(fid);
% end
fdirname=['RUNUP2_files/' filename];
fid=fopen(fdirname,'w');
for i=1:size(output,1)
fprintf(fid,'%s\r\n',output(i,:));
end
fclose(fid);


%% Function to reduce the number of transect stations to 20
    function [indices] = res_reduce_tbn(st,ele2,begin_pt,end_pt,num)
        %UNTITLED Summary of this function goes here
        %   Detailed explanation goes here
        % Modified by TBN on 8/16/18
        indices=find(st>=begin_pt & st<=end_pt);
        for f=2:length(indices)
            if ele2(indices(f)) < ele2(indices(f-1))
                ele2(indices(f))=ele2(indices(f-1));
            end
        end
        %set up algorithm
        % indices=begin_pt:end_pt;
        x2=st(indices);
        z2=ele2(indices);
        err=[];
        while length(x2)>num
            x1=x2;z1=z2; %reset x1 to equal x2, etc.
            for f_i=1:length(indices) %calculate error from removing each point
                x2=x1;z2=z1;
                x2(f_i)=[];z2(f_i)=[];
                a1=trapz(x1,z1);
                a2=trapz(x2,z2);
                err(f_i)=abs(a1-a2);
            end
            pt=find(err==min(err));
            indices(pt)=[];%remove point with least error
            x2=st(indices);z2=ele2(indices);
            err=[];
            plot(st,ele2,'k',st(indices),ele2(indices),'b');
            disp(['length: ' num2str(length(indices))])
        end
        
    end

%% Function to calculate +/- 5% wave height and wave period
    function [wlvl_out]=wlvl_range(wlvl,n)
    fac=0.05; % the factor to subtract and add from wave height and wave period
        for v=1:n
            wlvl_temp=wlvl{v};
            w_height(1)=wlvl_temp(2)-(fac*wlvl_temp(2));
            w_height(2)=wlvl_temp(2);
            w_height(3)=wlvl_temp(2)+(fac*wlvl_temp(2));
            w_per(1)=wlvl_temp(3)-(fac*wlvl_temp(3));
            w_per(2)=wlvl_temp(3);
            w_per(3)=wlvl_temp(3)+(0.05*wlvl_temp(3));
            count=0;
            for y=1:3
                for z=1:3
                    count=count+1;
                    wlvl_temp_out(count,:)=[wlvl_temp(1) w_height(y) w_per(z)];
                end
            end
            wlvl_out(v)={wlvl_temp_out};
            wlvl_temp_out=[];
            wlvl_temp=[];
        end
    end
                    
        
%% Function to de-shoal wave height
    function [H0_E,H0_H,HuntStr]=deshoal(H,T,D,g)
        % function [H0_E,H0_H]=deepwaterWaveheight(H,T,D,g);
        %
        %  H=wave height
        %  T=wave period
        %  D=depth where wave height was recorded
        %  g=gravitational acceleration
        %
        %  H0_E = deepwater wave height Eckart eqn
        %  H0_H = deepwater wave height Hunt eqn
        %
        %  computes deep water wave height from linear theory
        %  and conservation of energy assumptions
        %
        %  uses both the Eckart and Hunt equations
        %  for approximate solution of dispersion relation.
        %
        %  Hunt is probably more accurate.
        %
        %  Reference:
        %
        %  R.G. Dean and R.A. Dalrymple. 2000.  Water
        %  Wave Mechanics for Engineers and Scientists. World
        %  Scientific Publishing Company, River Edge New Jersy
        %
        %
        %  USACE (1985), Direct Methods for Calculating Wavelength, CETN-1-17
        %  US Army Engineer Waterways Experiment Station Coastel Engineering
        %  Research Center, Vicksburg, MS
        %
        %  also see CEM Part II-3 for discussion of shoaling coefficient
        %
        
        twopi=2*3.141592653589793;
        si=1;

        % get deep water celerity
        L0=g.*T.*T./twopi;
        C0=L0./T;
        HuntStr{si}=sprintf('    Depth, D = %.2f\n',D);si=si+1;
        HuntStr{si}=sprintf('    Period, T = %.2f\n',T);si=si+1;
        HuntStr{si}=sprintf('    Waveheight, H = %.2f\n',H);si=si+1;
        HuntStr{si}='Deep water wavelength, L0 (ft)\n';si=si+1;
        HuntStr{si}='    L0 = g*T*T/twopi\n';si=si+1;
        HuntStr{si}=sprintf('    L0 = %.2f*%.2f*%.2f/%.2f = %.2f\n',g,T,T,twopi,L0);si=si+1;
        HuntStr{si}='Deep water wave celerity, C0 (ft/s)\n';si=si+1;
        HuntStr{si}='    C0 = L0/T\n';si=si+1;
        HuntStr{si}=sprintf('    C0 = %.2f/%.2f = %.2f \n',L0,T,C0);si=si+1;
        
        
        % angular frequency
        sigma=twopi ./ T;
        sigmasqOg=sigma.*sigma./g;
        HuntStr{si}=sprintf('Angular frequency, sigma (rad/s)\n');si=si+1;
        HuntStr{si}='    sigma = twopi/T \n';si=si+1;
        HuntStr{si}=sprintf('    sigma = %.2f/%.2f = %.2f\n',twopi,T,sigma);si=si+1;
                
        % local wave length  (Eckart, 1951)
        k=sigmasqOg./sqrt(tanh(sigmasqOg.*D));
        L=twopi./k;
        
        
        % get local celerity at depth D
        C1=L./T;
        
        % or use Hunt's (1979) approximation for Celerity
        y=sigma.*sigma.*D./g;
        C1H=sqrt( g.*D./(y+1./(1 + 0.6522.*y + 0.4622.*y.^2 + 0.0864.*y.^4 +0.0675.*y.^5)) );
        HuntStr{si}='Hunts (1979) approximation for Celerity C1H (ft/s) at Depth D (ft)\n';si=si+1;
        HuntStr{si}='    y = sigma.*sigma.*D./g\n';si=si+1;
        HuntStr{si}=sprintf('    y = %.2f*%.2f*%.2f/%.2f = %.2f\n',sigma,sigma,D,g,y);si=si+1;
        HuntStr{si}='    C1H = sqrt( g.*D./(y+1./(1 + 0.6522.*y + 0.4622.*y.^2 + 0.0864.*y.^4 +0.0675.*y.^5)) )\n';si=si+1;
        HuntStr{si}=sprintf('    C1H = %.2f\n',C1H);si=si+1;
                
        % shoaling coefficient
        Ks=sqrt(C0./C1);
        KsH=sqrt(C0./C1H);
        HuntStr{si}='Shoaling Coefficient KsH \n';si=si+1;
        HuntStr{si}='    KsH = sqrt(C0/C1H)\n';si=si+1;
        HuntStr{si}=sprintf('    KsH = sqrt(%.2f/%.2f) = %.2f\n',C0,C1H,KsH);si=si+1;
        
        % get deep water wave height
        H0_E=H./Ks;
        H0_H=H./KsH;
        HuntStr{si}='Deepwater Wave Height H0_H (ft) \n';si=si+1;
        HuntStr{si}='    H0_H = H/KsH\n';si=si+1;
        HuntStr{si}=sprintf('    H0_H = %.2f/%.2f = %.2f\n ',H,KsH,H0_H);si=si+1;
    end

end

