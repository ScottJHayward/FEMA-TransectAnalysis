%WHAFIS INPUT
%
%Author: Scott Hayward
%Company: Ransom Consulting, inc. 
%Project: 2018 FEMA appeal, York and Cumberland Counties 
%
%This script is used to create WHAFIS input files for each transect.
%Additionally, a .bat file is written in order to execute WHAFIS$ for each
%transect. Also plots the elevation profiles used as input for WHAFIS. 
%

% chk nld 20190912

%%%%%%%%%%%%%%%%%%%%%%%%%%  config  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all
format long
fpre='whafis4/';                 %prefix
bfilename='whafis4/runWHAFIS.bat';
datafile='../data/transectdata.xls'; 
tDIR='../ADCIRC_returns/'; %location of transects
prefix='YK-';
swanpre='../2_swan/swanfiles/';
logpre='logfiles/'
%%%%%%%%%%%%%%%%%%%%%%%%  end config  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,txt,raw]=xlsread(datafile);
for i=2:size(raw,1)
   fnames{i-1}=raw{i,1};
end

twl=num(:,2);
hs=num(:,3);
per=num(:,4);
startx=num(:,1);

for i=1:length(fnames)

   fname=[tDIR fnames{i} 'XYZSTA_RETURNS.csv'];
   file=xlsread(fname);
   lon{i}=file(:,1);lat{i}=file(:,2);sta{i}=file(:,4);ele{i}=file(:,3);wl_adc{i}=file(:,7);

end

for i=1:length(fnames)
    
   %read swan returns
   fid1=fopen([swanpre fnames{i} '.dat'],'r');

   %read SWAN .dat file output
   li=fgetl(fid1);li=fgetl(fid1);li=fgetl(fid1);li=fgetl(fid1);
   li=fgetl(fid1);li=fgetl(fid1);li=fgetl(fid1);li=fgetl(fid1);
   n=0;
   while li ~=-1
      n=n+1;
      A=str2num(li);
      setup{i}(n)=A(10)*3.280833333;%setup in feet
      xswan{i}(n)=(n-1)*3.280833333+startx(i);
      li=fgetl(fid1);
   end

   %remove points where setup might reported as -9
   NaNpts=find(setup{i}<-8.9);
   setup{i}(NaNpts)=[];
   xswan{i}(NaNpts)=[];

   setup{i}(setup{i}<0)=0; % not considering setdown 

   %Save maximum adcirc water level for use in backwater areas
   adcwl_max{i}=max(wl_adc{i});

   % find points where SWAN 1-d shows setup
   lll=find(setup{i} > 0 );
   newStas{i}=xswan{i}(lll);
   newWls=twl(i)+setup{i}(lll);
   newEles=interp1(sta{i},ele{i},newStas{i});
   newLons=interp1(sta{i},lon{i},newStas{i});
   newLats=interp1(sta{i},lat{i},newStas{i});
   

   % remove original points that overlap the swan points
   % go 1 extra foot in each direction to ensure no duplicate WHAFIS cards
   lll=find(sta{i} >= min(newStas{i})-1 & sta{i} <= max(newStas{i})+1);
   sta{i}(lll)=[];
   ele{i}(lll)=[];
   wl_adc{i}(lll)=[];
   lon{i}(lll)=[];
   lat{i}(lll)=[];
  

   % add in the swan 1-d points 
   sta{i}=[sta{i}; newStas{i}']
   ele{i}=[ele{i}; newEles'];
   wl_adc{i}=[wl_adc{i}; newWls'];
   lon{i}=[lon{i}; newLons'];
   lat{i}=[lat{i}; newLats'];

   % sort them
   [sta{i},II]=sort(sta{i},'ascend');
   ele{i}=ele{i}(II);
   wl_adc{i}=wl_adc{i}(II);
   lon{i}=lon{i}(II);
   lat{i}=lat{i}(II);

   %interpolate SWAN 1d TWL onto the transect
  % wl_swn{i}=interp1(xswan,setup+twl(i),sta{i});

   %plot some things as a check
  % figure(i)
  % hold on
  % plot(sta{i},ele{i},'k-sq');
  % plot(sta{i},wl_adc{i},'b:');
  % plot(xswan{i},setup{i}+twl(i),'y*');
  %%% plot(sta{i},wl_swn{i},'g-o');
  % ylim([-20 20]);


end % loop over transects





%replace adcirc water levels with swan 1D TWL
%swanKeepers={};
%for n=1:length(fnames)
%   swanKeepers{n}=[];
%   for i=1:length(wl_adc{n})
%      if isnan(wl_swn{n}(i)) == 0 
%         wl_adc{n}(i)=wl_swn{n}(i)
%         swanKeepers{n}=[swanKeepers{n} i];
%      end
%   end
%end

for i=1:length(wl_adc)
    %now fill any remaining sections that did not contain adcirc elevations
    badjs=find(wl_adc{i} < 0 & ele{i} <= max(wl_adc{i}));%adcwl_max); note, wl_adc now includes swan 1-D wl
    if ~isempty (badjs)
        j=badjs(1); badjs(1)=[];
        while length(badjs) >= 0
            %look left for wl > -999 turn right if you hit land
            j_wl=j-1;
            wl_left=wl_adc{i}(j_wl);
            el_left=ele{i}(j_wl);
            while wl_left < -99 %|| el_left > max(wl_adc{i})
                if j_wl == 1 || el_left > max(wl_adc{i})
                    break
                end
                j_wl=j_wl-1;
                wl_left=wl_adc{i}(j_wl);
                el_left=ele{i}(j_wl);
            end
            if wl_left > -99
                wl_adc{i}(j)=wl_adc{i}(j_wl);
                if length(badjs)==0
                   break
                end
                j=badjs(1); badjs(1)=[];
                continue;
            end
            % we didn't find it left, go right
            j_wl=j+1;
            wl_rt=wl_adc{i}(j_wl);
            el_rt=ele{i}(j_wl);
            while wl_rt < -99
                if j_wl == length(ele{i}) || el_rt > max(wl_adc{i})
                    break
                end
                j_wl=j_wl+1;
                wl_rt=wl_adc{i}(j_wl);
                el_rt=ele{i}(j_wl);
            end
            if wl_rt > -99
                wl_adc{i}(j)=wl_adc{i}(j_wl);
                if length(badjs)==0
                   break
                end
                j=badjs(1); badjs(1)=[];
                %           continue;
            else
               wl_adc{i}(j)=adcwl_max{i};
               if length(badjs)==0
                   break
               end
               j=badjs(1); badjs(1)=[];
            end
            
        end
    end
end





batchfile=fopen(bfilename,'w');
fprintf(batchfile,'%s\n',['@ECHO OFF']);

%initial variables
% these are just place holder strings that replaced below
TITLE='header';
IE='IE  0.00 3.55695       1       1 11.5168     4.7    12.4   56.14   56.14';
IF='IF     4    3.50    0.00    0.00    0.00    0.00    0.00    0.00    0.00';
AS='AS    47   11.52    0.00    0.00    0.00    0.00    0.00    0.00    0.00';
ET='ET';
for i=1:length(fnames)
    
    fname=[fpre fnames{i} '.dat'];
    fnamein=[fnames{i} '.dat'];
    fnamein2=[fnames{i} '_ck.dat'];
    
    fnameout=[fnames{i} '.out'];
    fid=fopen(fname,'w');
    wl=num2str(twl(i));wl=[wl '       '];       %add spaces
    hhs=num2str(1.6*hs(i));hhs=[hhs '       '];   %Critical Wave Height
    pper=num2str(per(i));pper=[pper '       '];   
    
    %find the location of the toe, which will be used for WHAFIS
    j=find(abs(sta{i}-startx(i))==min(abs(sta{i}-startx(i))));
    shift=sta{i}(j);
    %shift station
    sta{i}=sta{i}-sta{i}(j);
    lonj1=lon{i}(j);latj1=lat{i}(j);
    
    % before reducing resolution interpolate points at the TWL and make a list of indices to keep
    xx=sta{i};  zz=ele{i}; ww=wl_adc{i}; lo=lon{i}; la=lat{i};
    k=0;
    xnew=[];
    znew=[];
    wnew=[];
    lonew=[];
    lanew=[];
    keepers=[j];% j+1];
    while ( length(xx) > 1)
        % pop the first value
        x0=xx(1);  xx(1)=[];  
        z0=zz(1);  zz(1)=[];
        w0=ww(1);  ww(1)=[];
        lo0=lo(1); lo(1)=[];
        la0=la(1); la(1)=[];
        
        % check for upcrossing
        if (w0 >= z0) & (ww(1) <= zz(1))
            xn=interp1([z0 zz(1)],[x0 xx(1)],w0);
            loN=interp1([z0 zz(1)],[lo0 lo(1)],w0);
            laN=interp1([z0 zz(1)],[la0 la(1)],w0);
            zn=w0;
            wn=w0
            % push two points
            xnew=[xnew; x0; xn];
            znew=[znew; z0; zn];
            wnew=[wnew; w0; wn];
            lonew=[lonew; lo0; loN];
            lanew=[lanew; la0; laN];
            k=k+2;
            keepers=[keepers k];
         % check for downcrossing
         elseif (w0 < z0) & (ww(1) > zz(1))
            xn=interp1([z0 zz(1)],[x0 xx(1)],ww(1));
            loN=interp1([z0 zz(1)],[lo0 lo(1)],ww(1));
            laN=interp1([z0 zz(1)],[la0 la(1)],ww(1));
            zn=ww(1);
            wn=ww(1)
            % push two points
            xnew=[xnew; x0; xn];
            znew=[znew; z0; zn];
            wnew=[wnew; w0; wn];
            lonew=[lonew; lo0; loN];
            lanew=[lanew; la0; laN];
            k=k+2;
            keepers=[keepers k];
         else
            % push poped point only
            xnew=[xnew; x0];
            znew=[znew; z0];
            wnew=[wnew; w0];
            lonew=[lonew; lo0];
            lanew=[lanew; la0];
            k=k+1;
         end
        % check for downcrossing          
    end   
    sta{i}=xnew;
    ele{i}=znew;
    wl_adc{i}=wnew;
    lon{i}=lonew;
    lat{i}=lanew;
    
    %make sure we dont delete locations with Setup
    swanKeepers{i}=find(sta{i} <= newStas{i}(end)-shift & sta{i} >= newStas{i}(1)-shift);
    indices_=res_reduce(sta{i},ele{i},find(sta{i}==0),length(sta{i}),500);

    % add the keepers back
    % but first make sure swanKeepers doesn't put us over 1000 cards
   nCanKeep=999-length(indices_)-length(keepers);
   if length(swanKeepers{i}) > nCanKeep
       swanKeepers{i}(1:length(swanKeepers{i})-nCanKeep)=[];
   end
 
   indices=unique([indices_ keepers swanKeepers{i}']);
%     
    
    %populate until end of transect
    j=indices(1);
    %write IE card
    elev=[num2str(ele{i}(j)) '      '];
    %IE card
    line=IE;
    %add trimmed variables to the file. 
    %field 0: IE 
    %field 1: station at beginning: 
    station=['0.0      '];
    IE(3:8)=station(1:6);
    %field 2: elevation at station 1: 
    IE(10:16)=elev(1:7);
    %field 3: fetch length (miles): 
    %field 4: 10% SWEL
    %field 5: 1%  SWEL
    wl= [num2str(wl_adc{i}(j)) '       ']; 
    IE(34:40)=wl(1:7);
    %field 6: initial Wave Height
    IE(42:48)=hhs(1:7);
    %field 7: wave period
    IE(50:56)=pper(1:7);
    fprintf(fid,'%s\n',TITLE);
    fprintf(fid,'%s\n',IE);
    
    n=2;

    
    while n < length(indices)
        j=indices(n);
        %%%%%%%%%%%%%%%%%% IF card %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ele{i}(j) < wl_adc{i}(j)
    	   elev= [num2str(ele{i}(j)) '      '];
           %station=  [num2str(sta{i}(j)) '      '];
           station=  [ sprintf('%0.1f',sta{i}(j))    '      '];
           wl=   [num2str(wl_adc{i}(j)) '       ']; 
           %field 1: station 
           IF(3:8)=station(1:6);
           %field 2: elevation
           IF(10:16)=elev(1:7);
           %field 3: water level
           IF(26:32)=wl(1:7);
           if ele{i}(j) <= 0
             IF(1:1)='O';
           else
             IF(1:1)='I';
           end
           fprintf(fid,'%s\n',IF);
           n=n+1;
        elseif  ele{i}(j) == wl_adc{i}(j)
           % write an IF card right at the intersection of TWL with ground, then an AS card at the last point above surge
           elev= [num2str(ele{i}(j)) '      '];
           station=  [ sprintf('%0.1f',sta{i}(j))    '      '];
           wl=   [num2str(wl_adc{i}(j)) '       ']; 
           IF(3:8)=station(1:6);
           IF(10:16)=elev(1:7);
           IF(26:32)=wl(1:7);
           if ele{i}(j) <= 0
             IF(1:1)='O';
           else
             IF(1:1)='I';
           end
           fprintf(fid,'%s\n',IF);
           % now loop through points until you find the next one below the surge, write the AS card, then continue
           while (1==1)
               if n==length(indices)
                  break;
               end
               if ele{i}(indices(n+1)) < wl_adc{i}(indices(n+1)) %|| n+1==length(indices)
                  %%%%%%%%%%%%%%%%%% AS card %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  station=  [ sprintf('%0.1f',sta{i}(indices(n)))    '      '];
                  %wl=   [num2str(wl_adc{i}(indices(n))) '       '];
                  %field 1: station 
                  AS(3:8)=station(1:6);
                  %field 2: elevation
                  aselev= [num2str(wl_adc{i}(indices(n+1))) '       ']; %Get AS elevation from next water point
                  AS(9:16)=[' ',aselev(1:7)];  %print elevation to AS card
                  AS(25:32)=[' ',aselev(1:7)];  %print water elevation to AS card
                  fprintf(fid,'%s\n',AS);
                  n=n+1;
                  break
               else
                  n=n+1;
               end
           end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%% PS card %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    st=find(abs(sta{i}-startx(i))==min(abs(sta{1,i}-startx(i))));
    st=1;

    en=j;
    lonj2=lon{i}(indices(1));latj2=lat{i}(indices(1));
    [sx,sy,utmz]=deg2utm(lat{i}(indices(1)),lon{i}(indices(1)));
    [ex,ey,utmz]=deg2utm(lat{i}(indices(end)),lon{i}(indices(end)));
    PS1=['PS START(' num2str(sx) ',' num2str(sy) ')'];
    PS2=['PS END(' num2str(ex) ',' num2str(ey) ')'];
    fprintf(fid,'%s\n',PS1);
    fprintf(fid,'%s\n',PS2);


    %%%%%%%%%%%%%%%%%%% ET card %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fid,'%s\n',ET);
    fprintf(fid,'%s\n',' ');
    fprintf(fid,'%s\n',' ');
    
    
    
    fprintf(fid,'%s\n',' ');

  
    fprintf(batchfile,'%s\n',['WHAFIS4.exe ' fnamein ' ' fnameout]);
    
    
    
    %%%%%%%%%%%%%%%%% PRINT TRANSECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    tstr=['WHAFIS Elevation, and Total stillWater Level, transect: ' fnames{i}];
    wl_adc{i}(wl_adc{i} < -99)=NaN;
    plot(sta{i},ele{i},'k');hold on
    plot(xswan{i}-shift,setup{i}+twl(i),'o','color',[0.8500, 0.3250, 0.0980],'markersize',4);
    plot(sta{i}(indices),wl_adc{i}(indices),'b-d','markersize',2);
    legend('Transect Profile','SWAN 1-D points', 'TSWL','location','best')
    title({'{\color{red}REVISED SEP-05-2019}',tstr})
    xlabel('Station (ft)');ylabel('Elevation (ft)')
    ylim([0 20])
    set(gcf,'position',[100 100 900 600],'paperorientation','Landscape')
    grid on 
    set(gca,'xminorgrid','on','yminorgrid','on');
    print('-r450','-dpdf',fnames{i});
    pause(.2)
    
    %%%%%%%%%%%%%%%%% LOG FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    logfile=[logpre fnames{i} '_log.txt']
    fid=fopen(logfile,'w');
    
    fprintf(fid,'%s\n','_______________________________________________________');
    fprintf(fid,'%s\n','PART 3: WHAFIS');
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n',['WHAFIS input:  ' fnamein] );
    fprintf(fid,'%s\n',['WHAFIS output: ' fnameout] );
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n','PART 3 COMPLETE________________________________________');
    
    
end
fclose all
