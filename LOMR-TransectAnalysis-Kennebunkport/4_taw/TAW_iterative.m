clear all
close all
format long g

diary logfiles/YK-92-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Kennebunkport, York county, Maine
% TRANSECT ID: YK-92
% calculation by SJH, Ransom Consulting, Inc. 02-Apr-2020
% 100-year wave runup using TAW methodology
% including berm and weighted average with foreshore if necessary
%
% chk nld 20200220
%
% This script assumes that the incident wave conditions provided
% as input in the configuration section below are the 
% appropriate values located at the end of the foreshore
% or toe of the slope on which the run-up is being calculated
% the script does not attempt to apply a depth limit or any other
% transformation to the incident wave conditions other than 
% conversion of the peak wave period to the spectral mean wave 
% as recommended in the references below
% 
% references:
%
% Van der Meer, J.W., 2002. Technical Report Wave Run-up and 
% Wave Overtopping at Dikes. TAW Technical Advisory Committee on
% Flood Defence, The Netherlands.
%
% FEMA. 2007,  Atlantic Ocean and Gulf of Mexico Coastal Guidelines Update
%
% 
%

%--------------------------------------------------------------------------
% CONFIG
%--------------------------------------------------------------------------
fname='inpfiles/YK-92sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-92-runup';                            
SWEL=8.8306;  % 100-yr still water level including wave setup. 
H0=5.1811;    % significant wave height at toe of structure 
Tp=14.019;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.80989;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.9778; 
maxSetup=1.4808;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-92'


% END CONFIG
%-------------------------------

SWEL=SWEL+setupAtToe  
SWEL_fore=SWEL+maxSetup


% FIND WAVELENGTH USING DEEPWATER DISPERSION RELATION
% using English units
L0=32.15/(2*pi)*T0^2

% Find Hb (Munk, 1949)
%Hb=H0/(3.3*(H0/L0)^(1/3)) 
%Db=-Hb/.78+SWEL;  % depth at breaking

% The toe elevation here is only used to determine the average 
% structure slope, it is not used to depth limit the wave height.
% Any depth limiting or other modification of the wave height 
% to make it consitent with TAW guidance should be performed 
% prior to the input of the significant wave height given above. 
Ztoe=SWEL-1.5*H0

% read the transect
[sta,dep,inc] = textread(fname,'%n%n%n%*[^\n]','delimiter',',','headerlines',0);

% remove unselected points
k=find(inc==0);
sta(k)=[];
dep(k)=[];

sta_org=sta;  % used for plotting purposes
dep_org=dep;


% initial guess at maximum run-up elevation to estimate slope
Z2=SWEL+1.5*H0

% determine station at the max runup and -1.5*H0 (i.e. the toe)
top_sta=-999;
toe_sta=-999;
for kk=1:length(sta)-1
    if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
       top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
    end
    if  ((Ztoe > dep(kk)) & (Ztoe <= dep(kk+1)))    % here is the intersection of Ztoe with profile
       toe_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Ztoe)
    end
end
% check to make sure we got them, if not extend the end slopes outward
S=diff(dep)./diff(sta);
if toe_sta==-999
   dy=dep(1)-Ztoe;
   toe_sta=sta(1)-dy/S(1)
end
if top_sta==-999
   dy=Z2-dep(end);
   top_sta=sta(end)+dy/S(end)
end

% just so the reader can tell the values aren't -999 anymore
top_sta
toe_sta



% check for case where the toe of slope is below SWL-1.5*H0
% in this case interpolate setup from the setupAtToe(really setup as first station), and the max setup
% also un-include points seaward of SWL-1.5*H0
if Ztoe > dep(1)
   dd=SWEL_fore-dep;
   k=find(dd<0,1); % k is index of first land point
   staAtSWL=interp1(dep(k-1:k),sta(k-1:k),SWEL_fore);   
   dsta=staAtSWL-sta(1);
   dsetup=maxSetup-setupAtToe;
   dsetdsta=dsetup/dsta;
   setup=setupAtToe+dsetdsta*(toe_sta-sta(1));
   sprintf('-!!- Location of SWEL-1.5*H0 is %4.1f ft landward of toe of slope',dsta)
   sprintf('-!!- Setup is interpolated between setup at toe of slope and max setup')
   sprintf('-!!-       setup is adjusted to %4.2f feet',setup)
   SWEL=SWEL-setupAtToe+setup;
   sprintf('-!!-       SWEL is adjusted to %4.2f feet',SWEL)
   k=find(dep < SWEL-1.5*H0)
   sta(k)=[];
   dep(k)=[];
else
   sprintf('-!!- The User has selected a starting point that is %4.2f feet above the elevation of SWEL-1.5H0\n',dep(1)-Ztoe)
   sprintf('-!!- This may be reasonable for some cases.  However the user may want to consider:\n')
   sprintf('-!!-   1) Selecting a starting point that is at or below %4.2f feet elevation, or\n', Ztoe)
   sprintf('-!!-   2) Reducing the incident wave height to a depth limited condition.\n')
end



% now iterate converge on a runup elevation
tol=0.01;  % convergence criteria
R2del=999;
R2_new=3*H0; %initial guess
R2=R2_new;
iter=0;
R2_all=[];
topStaAll=[];
Berm_Segs=[];
TAW_ALWAYS_VALID=1;
while(abs(R2del) > tol && iter <= 25)
    iter=iter+1;
    sprintf ('!----------------- STARTING ITERATION %d -----------------!',iter)
    % elevation of toe of slope
    Ztoe
    % station of toe slope (relative to 0-NAVD88 shoreline
    toe_sta 
    % station of top of slope/extent of 2% run-up
    top_sta
    % elevation of top of slope/extent of 2% run-up
    Z2
    % incident significant wave height
    H0
    % incident spectral peak wave period
    Tp
    % incident spectral mean wave period
    T0
        
    R2=R2_new  
    Z2=R2+SWEL
    % determine slope for this iteration
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end)
    end
    
    % get the length of the slope (not accounting for berm)
    Lslope=top_sta-toe_sta

    % loop over profile segments to determine berm factor 
    % re-calculate influence of depth of berm based on this run-up elevation
    % check for berm, berm width, berm height
    berm_width=0;
    rdh_sum=0;
    Berm_Segs=[];
    Berm_Heights=[];
    for kk=1:length(sta)-1
       ddep=dep(kk+1)-dep(kk);
       dsta=sta(kk+1)-sta(kk);
       s=ddep/dsta;
       if  (s < 1/15)     % count it as a berm if slope is flatter than 1:15 (see TAW manual)
          sprintf ('Berm Factor Calculation: Iteration %d, Profile Segment: %d',iter,kk) 
          berm_width=berm_width+dsta;   % tally the width of all berm segments
          % compute the rdh for this segment and weight it by the segment length
          dh=SWEL-(dep(kk)+dep(kk+1))/2  
          if  dh < 0
              chi=R2;
          else 
              chi=2* H0;
          end
          if (dh <= R2 & dh >=-2*H0)
             rdh=(0.5-0.5*cos(3.14159*dh/chi)) ;
          else
             rdh=1;
          end
          rdh_sum=rdh_sum + rdh * dsta
          Berm_Segs=[Berm_Segs, kk];
          Berm_Heights=[Berm_Heights, (dep(kk)+dep(kk+1))/2];
       end
       if dep(kk) >= Z2  % jump out of loop if we reached limit of run-up for this iteration 
          break
       end
    end
    sprintf ('!------- End Berm Factor Calculation, Iter: %d ---------!',iter)
    berm_width
    rB=berm_width/Lslope 
    if (berm_width > 0)
       rdh_mean=rdh_sum/berm_width
    else
       rdh_mean=1
    end
    gamma_berm=1- rB * (1-rdh_mean)
    if gamma_berm > 1
       gamma_berm=1
    end
    if gamma_berm < 0.6
       gamma_berm =0.6
    end
    % Iribarren number
    slope=(Z2-Ztoe)/(Lslope-berm_width)
    Irb=(slope/(sqrt(H0/L0)))
    % runup height
    gamma_berm
    gamma_perm
    gamma_beta
    gamma_rough
    gamma=gamma_berm*gamma_perm*gamma_beta*gamma_rough

    % check validity
    TAW_VALID=1;
    if (Irb*gamma_berm < 0.5 | Irb*gamma_berm > 10 )
       sprintf('!!! - - Iribaren number: %6.2f is outside the valid range (0.5-10), TAW NOT VALID - - !!!\n', Irb*gamma_berm)
       TAW_VALID=0;
    else
       sprintf('!!! - - Iribaren number: %6.2f is in the valid range (0.5-10), TAW RECOMMENDED - - !!!\n', Irb*gamma_berm)
    end
    islope=1/slope;  
    if (slope < 1/8 | slope > 1)
       sprintf('!!! - - slope: 1:%3.1f V:H is outside the valid range (1:8 - 1:1), TAW NOT VALID - - !!!\n', islope)
       TAW_VALID=0;
    else
       sprintf('!!! - - slope: 1:%3.1f V:H is in the valid range (1:8 - 1:1), TAW RECOMMENDED - - !!!\n', islope)
    end
    if TAW_VALID == 0
       TAW_ALWAYS_VALID=0;
    end
         
    if (Irb*gamma_berm < 1.8)
       R2_new=gamma*H0*1.77*Irb
    else
       R2_new=gamma*H0*(4.3-(1.6/sqrt(Irb))) 
    end


    % check to see if we need to evaluate a shallow foreshore
    if berm_width > 0.25 * L0;
       disp ('!   Berm_width is greater than 1/4 wave length')
       disp ('!   Runup will be weighted average with foreshore calculation assuming depth limited wave height on berm');
       % do the foreshore calculation 
       fore_H0=0.78*(SWEL_fore-min(Berm_Heights))
       % get upper slope
       fore_toe_sta=-999;
       fore_toe_dep=-999;
       for kk=length(dep)-1:-1:1
          ddep=dep(kk+1)-dep(kk);
          dsta=sta(kk+1)-sta(kk);
          s=ddep/dsta;
          if s < 1/15
             break
          end
          fore_toe_sta=sta(kk);
          fore_toe_dep=dep(kk);
          upper_slope=(Z2-fore_toe_dep)/(top_sta-fore_toe_sta)
       end
       fore_Irb=upper_slope/(sqrt(fore_H0/L0));
       fore_gamma=gamma_perm*gamma_beta*gamma_rough;
       if (fore_Irb < 1.8)
          fore_R2=fore_gamma*fore_H0*1.77*fore_Irb;
       else
          fore_R2=fore_gamma*fore_H0*(4.3-(1.6/sqrt(fore_Irb))); 
       end
       if berm_width >= L0
          R2_new=fore_R2
          disp ('berm is wider than one wavelength, use full shallow foreshore solution');
       else
          w2=(berm_width-0.25*L0)/(0.75*L0)
          w1=1-w2
          R2_new=w2*fore_R2 + w1*R2_new
       end
    end % end berm width check

    % convergence criterion
    R2del=abs(R2-R2_new)
    R2_all(iter)=R2_new;

    % get the new top station (for plot purposes)
    Z2=R2_new+SWEL
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end);
    end
    topStaAll(iter)=top_sta;

end

% final 2% runup elevation
Z2=R2_new+SWEL

diary off


% make some plots


% get the full transect



file=xlsread('../ADCIRC_returns/YK-92XYZSTA_RETURNS.csv');
sta_full=file(:,4);
ele_full=file(:,3);


%determine SWL station
SWL_sta=-999;
for kk=1:length(sta)-1
    if ((SWEL > dep(kk)) & (SWEL <= dep(kk+1)))    % here is the intersection of z2 with profile
       SWL_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),SWEL)
       break;
    end
end

figure('position', [100 100 900 900])
hold on
plot(sta_full,ele_full,'linewidth',2,'color',[0.5 0.2 0.1]);
plot(sta,dep,'md');
grid on
plot([min(sta),SWL_sta],[SWEL SWEL],'color','b','linewidth',2)
%plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
scatter(topStaAll,R2_all+SWEL,30,'r*');
%scatter(toe_sta,Db,70,'k*');
scatter(toe_sta,Ztoe,70,'k*');
% plot berm segs
for nn=1:length(Berm_Segs)
   seg=Berm_Segs(nn);
   h=Berm_Heights(nn);
   plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
end
xlim([toe_sta-50 max(topStaAll)+50])
grid minor
hleg=legend('Ground','TAW Segments','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
set (hleg,'location','southeast');
xlabel('Distance from zero NAVD88 shoreline (feet)');
ylabel('Elevation (Feet-NAVD88)');
str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
yylm=get (gca,'ylim');
dztxt=(yylm(end)-yylm(1))/15;
text(toe_sta,Z2+dztxt/4,str);
str=sprintf('Run-up Height is %5.1f Feet',R2_new);
text(toe_sta,Z2-dztxt/4,str);
title (plotTitle);
set(gcf,'color','w');
print('-r450','-dpdf',imgname);

-1.000000e+00
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S2');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U2');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O2');
clear all
close all
format long g

diary logfiles/YK-99-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Kennebunkport, York county, Maine
% TRANSECT ID: YK-99-1
% calculation by SJH, Ransom Consulting, Inc. 02-Apr-2020
% 100-year wave runup using TAW methodology
% including berm and weighted average with foreshore if necessary
%
% chk nld 20200220
%
% This script assumes that the incident wave conditions provided
% as input in the configuration section below are the 
% appropriate values located at the end of the foreshore
% or toe of the slope on which the run-up is being calculated
% the script does not attempt to apply a depth limit or any other
% transformation to the incident wave conditions other than 
% conversion of the peak wave period to the spectral mean wave 
% as recommended in the references below
% 
% references:
%
% Van der Meer, J.W., 2002. Technical Report Wave Run-up and 
% Wave Overtopping at Dikes. TAW Technical Advisory Committee on
% Flood Defence, The Netherlands.
%
% FEMA. 2007,  Atlantic Ocean and Gulf of Mexico Coastal Guidelines Update
%
% 
%

%--------------------------------------------------------------------------
% CONFIG
%--------------------------------------------------------------------------
fname='inpfiles/YK-99-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-99-1-runup';                            
SWEL=9.4014;  % 100-yr still water level including wave setup. 
H0=2.142;    % significant wave height at toe of structure 
Tp=7.9359;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.0047506; 
maxSetup=0.30665;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-99-1'


% END CONFIG
%-------------------------------

SWEL=SWEL+setupAtToe  
SWEL_fore=SWEL+maxSetup


% FIND WAVELENGTH USING DEEPWATER DISPERSION RELATION
% using English units
L0=32.15/(2*pi)*T0^2

% Find Hb (Munk, 1949)
%Hb=H0/(3.3*(H0/L0)^(1/3)) 
%Db=-Hb/.78+SWEL;  % depth at breaking

% The toe elevation here is only used to determine the average 
% structure slope, it is not used to depth limit the wave height.
% Any depth limiting or other modification of the wave height 
% to make it consitent with TAW guidance should be performed 
% prior to the input of the significant wave height given above. 
Ztoe=SWEL-1.5*H0

% read the transect
[sta,dep,inc] = textread(fname,'%n%n%n%*[^\n]','delimiter',',','headerlines',0);

% remove unselected points
k=find(inc==0);
sta(k)=[];
dep(k)=[];

sta_org=sta;  % used for plotting purposes
dep_org=dep;


% initial guess at maximum run-up elevation to estimate slope
Z2=SWEL+1.5*H0

% determine station at the max runup and -1.5*H0 (i.e. the toe)
top_sta=-999;
toe_sta=-999;
for kk=1:length(sta)-1
    if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
       top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
    end
    if  ((Ztoe > dep(kk)) & (Ztoe <= dep(kk+1)))    % here is the intersection of Ztoe with profile
       toe_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Ztoe)
    end
end
% check to make sure we got them, if not extend the end slopes outward
S=diff(dep)./diff(sta);
if toe_sta==-999
   dy=dep(1)-Ztoe;
   toe_sta=sta(1)-dy/S(1)
end
if top_sta==-999
   dy=Z2-dep(end);
   top_sta=sta(end)+dy/S(end)
end

% just so the reader can tell the values aren't -999 anymore
top_sta
toe_sta



% check for case where the toe of slope is below SWL-1.5*H0
% in this case interpolate setup from the setupAtToe(really setup as first station), and the max setup
% also un-include points seaward of SWL-1.5*H0
if Ztoe > dep(1)
   dd=SWEL_fore-dep;
   k=find(dd<0,1); % k is index of first land point
   staAtSWL=interp1(dep(k-1:k),sta(k-1:k),SWEL_fore);   
   dsta=staAtSWL-sta(1);
   dsetup=maxSetup-setupAtToe;
   dsetdsta=dsetup/dsta;
   setup=setupAtToe+dsetdsta*(toe_sta-sta(1));
   sprintf('-!!- Location of SWEL-1.5*H0 is %4.1f ft landward of toe of slope',dsta)
   sprintf('-!!- Setup is interpolated between setup at toe of slope and max setup')
   sprintf('-!!-       setup is adjusted to %4.2f feet',setup)
   SWEL=SWEL-setupAtToe+setup;
   sprintf('-!!-       SWEL is adjusted to %4.2f feet',SWEL)
   k=find(dep < SWEL-1.5*H0)
   sta(k)=[];
   dep(k)=[];
else
   sprintf('-!!- The User has selected a starting point that is %4.2f feet above the elevation of SWEL-1.5H0\n',dep(1)-Ztoe)
   sprintf('-!!- This may be reasonable for some cases.  However the user may want to consider:\n')
   sprintf('-!!-   1) Selecting a starting point that is at or below %4.2f feet elevation, or\n', Ztoe)
   sprintf('-!!-   2) Reducing the incident wave height to a depth limited condition.\n')
end



% now iterate converge on a runup elevation
tol=0.01;  % convergence criteria
R2del=999;
R2_new=3*H0; %initial guess
R2=R2_new;
iter=0;
R2_all=[];
topStaAll=[];
Berm_Segs=[];
TAW_ALWAYS_VALID=1;
while(abs(R2del) > tol && iter <= 25)
    iter=iter+1;
    sprintf ('!----------------- STARTING ITERATION %d -----------------!',iter)
    % elevation of toe of slope
    Ztoe
    % station of toe slope (relative to 0-NAVD88 shoreline
    toe_sta 
    % station of top of slope/extent of 2% run-up
    top_sta
    % elevation of top of slope/extent of 2% run-up
    Z2
    % incident significant wave height
    H0
    % incident spectral peak wave period
    Tp
    % incident spectral mean wave period
    T0
        
    R2=R2_new  
    Z2=R2+SWEL
    % determine slope for this iteration
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end)
    end
    
    % get the length of the slope (not accounting for berm)
    Lslope=top_sta-toe_sta

    % loop over profile segments to determine berm factor 
    % re-calculate influence of depth of berm based on this run-up elevation
    % check for berm, berm width, berm height
    berm_width=0;
    rdh_sum=0;
    Berm_Segs=[];
    Berm_Heights=[];
    for kk=1:length(sta)-1
       ddep=dep(kk+1)-dep(kk);
       dsta=sta(kk+1)-sta(kk);
       s=ddep/dsta;
       if  (s < 1/15)     % count it as a berm if slope is flatter than 1:15 (see TAW manual)
          sprintf ('Berm Factor Calculation: Iteration %d, Profile Segment: %d',iter,kk) 
          berm_width=berm_width+dsta;   % tally the width of all berm segments
          % compute the rdh for this segment and weight it by the segment length
          dh=SWEL-(dep(kk)+dep(kk+1))/2  
          if  dh < 0
              chi=R2;
          else 
              chi=2* H0;
          end
          if (dh <= R2 & dh >=-2*H0)
             rdh=(0.5-0.5*cos(3.14159*dh/chi)) ;
          else
             rdh=1;
          end
          rdh_sum=rdh_sum + rdh * dsta
          Berm_Segs=[Berm_Segs, kk];
          Berm_Heights=[Berm_Heights, (dep(kk)+dep(kk+1))/2];
       end
       if dep(kk) >= Z2  % jump out of loop if we reached limit of run-up for this iteration 
          break
       end
    end
    sprintf ('!------- End Berm Factor Calculation, Iter: %d ---------!',iter)
    berm_width
    rB=berm_width/Lslope 
    if (berm_width > 0)
       rdh_mean=rdh_sum/berm_width
    else
       rdh_mean=1
    end
    gamma_berm=1- rB * (1-rdh_mean)
    if gamma_berm > 1
       gamma_berm=1
    end
    if gamma_berm < 0.6
       gamma_berm =0.6
    end
    % Iribarren number
    slope=(Z2-Ztoe)/(Lslope-berm_width)
    Irb=(slope/(sqrt(H0/L0)))
    % runup height
    gamma_berm
    gamma_perm
    gamma_beta
    gamma_rough
    gamma=gamma_berm*gamma_perm*gamma_beta*gamma_rough

    % check validity
    TAW_VALID=1;
    if (Irb*gamma_berm < 0.5 | Irb*gamma_berm > 10 )
       sprintf('!!! - - Iribaren number: %6.2f is outside the valid range (0.5-10), TAW NOT VALID - - !!!\n', Irb*gamma_berm)
       TAW_VALID=0;
    else
       sprintf('!!! - - Iribaren number: %6.2f is in the valid range (0.5-10), TAW RECOMMENDED - - !!!\n', Irb*gamma_berm)
    end
    islope=1/slope;  
    if (slope < 1/8 | slope > 1)
       sprintf('!!! - - slope: 1:%3.1f V:H is outside the valid range (1:8 - 1:1), TAW NOT VALID - - !!!\n', islope)
       TAW_VALID=0;
    else
       sprintf('!!! - - slope: 1:%3.1f V:H is in the valid range (1:8 - 1:1), TAW RECOMMENDED - - !!!\n', islope)
    end
    if TAW_VALID == 0
       TAW_ALWAYS_VALID=0;
    end
         
    if (Irb*gamma_berm < 1.8)
       R2_new=gamma*H0*1.77*Irb
    else
       R2_new=gamma*H0*(4.3-(1.6/sqrt(Irb))) 
    end


    % check to see if we need to evaluate a shallow foreshore
    if berm_width > 0.25 * L0;
       disp ('!   Berm_width is greater than 1/4 wave length')
       disp ('!   Runup will be weighted average with foreshore calculation assuming depth limited wave height on berm');
       % do the foreshore calculation 
       fore_H0=0.78*(SWEL_fore-min(Berm_Heights))
       % get upper slope
       fore_toe_sta=-999;
       fore_toe_dep=-999;
       for kk=length(dep)-1:-1:1
          ddep=dep(kk+1)-dep(kk);
          dsta=sta(kk+1)-sta(kk);
          s=ddep/dsta;
          if s < 1/15
             break
          end
          fore_toe_sta=sta(kk);
          fore_toe_dep=dep(kk);
          upper_slope=(Z2-fore_toe_dep)/(top_sta-fore_toe_sta)
       end
       fore_Irb=upper_slope/(sqrt(fore_H0/L0));
       fore_gamma=gamma_perm*gamma_beta*gamma_rough;
       if (fore_Irb < 1.8)
          fore_R2=fore_gamma*fore_H0*1.77*fore_Irb;
       else
          fore_R2=fore_gamma*fore_H0*(4.3-(1.6/sqrt(fore_Irb))); 
       end
       if berm_width >= L0
          R2_new=fore_R2
          disp ('berm is wider than one wavelength, use full shallow foreshore solution');
       else
          w2=(berm_width-0.25*L0)/(0.75*L0)
          w1=1-w2
          R2_new=w2*fore_R2 + w1*R2_new
       end
    end % end berm width check

    % convergence criterion
    R2del=abs(R2-R2_new)
    R2_all(iter)=R2_new;

    % get the new top station (for plot purposes)
    Z2=R2_new+SWEL
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end);
    end
    topStaAll(iter)=top_sta;

end

% final 2% runup elevation
Z2=R2_new+SWEL

diary off


% make some plots


% get the full transect



file=xlsread('../ADCIRC_returns/YK-99-1XYZSTA_RETURNS.csv');
sta_full=file(:,4);
ele_full=file(:,3);


%determine SWL station
SWL_sta=-999;
for kk=1:length(sta)-1
    if ((SWEL > dep(kk)) & (SWEL <= dep(kk+1)))    % here is the intersection of z2 with profile
       SWL_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),SWEL)
       break;
    end
end

figure('position', [100 100 900 900])
hold on
plot(sta_full,ele_full,'linewidth',2,'color',[0.5 0.2 0.1]);
plot(sta,dep,'md');
grid on
plot([min(sta),SWL_sta],[SWEL SWEL],'color','b','linewidth',2)
%plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
scatter(topStaAll,R2_all+SWEL,30,'r*');
%scatter(toe_sta,Db,70,'k*');
scatter(toe_sta,Ztoe,70,'k*');
% plot berm segs
for nn=1:length(Berm_Segs)
   seg=Berm_Segs(nn);
   h=Berm_Heights(nn);
   plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
end
xlim([toe_sta-50 max(topStaAll)+50])
grid minor
hleg=legend('Ground','TAW Segments','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
set (hleg,'location','southeast');
xlabel('Distance from zero NAVD88 shoreline (feet)');
ylabel('Elevation (Feet-NAVD88)');
str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
yylm=get (gca,'ylim');
dztxt=(yylm(end)-yylm(1))/15;
text(toe_sta,Z2+dztxt/4,str);
str=sprintf('Run-up Height is %5.1f Feet',R2_new);
text(toe_sta,Z2-dztxt/4,str);
title (plotTitle);
set(gcf,'color','w');
print('-r450','-dpdf',imgname);

-1.000000e+00
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S4');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U4');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O4');
clear all
close all
format long g

diary logfiles/YK-99-2-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Kennebunkport, York county, Maine
% TRANSECT ID: YK-99-2
% calculation by SJH, Ransom Consulting, Inc. 02-Apr-2020
% 100-year wave runup using TAW methodology
% including berm and weighted average with foreshore if necessary
%
% chk nld 20200220
%
% This script assumes that the incident wave conditions provided
% as input in the configuration section below are the 
% appropriate values located at the end of the foreshore
% or toe of the slope on which the run-up is being calculated
% the script does not attempt to apply a depth limit or any other
% transformation to the incident wave conditions other than 
% conversion of the peak wave period to the spectral mean wave 
% as recommended in the references below
% 
% references:
%
% Van der Meer, J.W., 2002. Technical Report Wave Run-up and 
% Wave Overtopping at Dikes. TAW Technical Advisory Committee on
% Flood Defence, The Netherlands.
%
% FEMA. 2007,  Atlantic Ocean and Gulf of Mexico Coastal Guidelines Update
%
% 
%

%--------------------------------------------------------------------------
% CONFIG
%--------------------------------------------------------------------------
fname='inpfiles/YK-99-2sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-99-2-runup';                            
SWEL=9.3261;  % 100-yr still water level including wave setup. 
H0=5.0379;    % significant wave height at toe of structure 
Tp=11.3516;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.94117;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.039442; 
maxSetup=0.93903;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-99-2'


% END CONFIG
%-------------------------------

SWEL=SWEL+setupAtToe  
SWEL_fore=SWEL+maxSetup


% FIND WAVELENGTH USING DEEPWATER DISPERSION RELATION
% using English units
L0=32.15/(2*pi)*T0^2

% Find Hb (Munk, 1949)
%Hb=H0/(3.3*(H0/L0)^(1/3)) 
%Db=-Hb/.78+SWEL;  % depth at breaking

% The toe elevation here is only used to determine the average 
% structure slope, it is not used to depth limit the wave height.
% Any depth limiting or other modification of the wave height 
% to make it consitent with TAW guidance should be performed 
% prior to the input of the significant wave height given above. 
Ztoe=SWEL-1.5*H0

% read the transect
[sta,dep,inc] = textread(fname,'%n%n%n%*[^\n]','delimiter',',','headerlines',0);

% remove unselected points
k=find(inc==0);
sta(k)=[];
dep(k)=[];

sta_org=sta;  % used for plotting purposes
dep_org=dep;


% initial guess at maximum run-up elevation to estimate slope
Z2=SWEL+1.5*H0

% determine station at the max runup and -1.5*H0 (i.e. the toe)
top_sta=-999;
toe_sta=-999;
for kk=1:length(sta)-1
    if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
       top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
    end
    if  ((Ztoe > dep(kk)) & (Ztoe <= dep(kk+1)))    % here is the intersection of Ztoe with profile
       toe_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Ztoe)
    end
end
% check to make sure we got them, if not extend the end slopes outward
S=diff(dep)./diff(sta);
if toe_sta==-999
   dy=dep(1)-Ztoe;
   toe_sta=sta(1)-dy/S(1)
end
if top_sta==-999
   dy=Z2-dep(end);
   top_sta=sta(end)+dy/S(end)
end

% just so the reader can tell the values aren't -999 anymore
top_sta
toe_sta



% check for case where the toe of slope is below SWL-1.5*H0
% in this case interpolate setup from the setupAtToe(really setup as first station), and the max setup
% also un-include points seaward of SWL-1.5*H0
if Ztoe > dep(1)
   dd=SWEL_fore-dep;
   k=find(dd<0,1); % k is index of first land point
   staAtSWL=interp1(dep(k-1:k),sta(k-1:k),SWEL_fore);   
   dsta=staAtSWL-sta(1);
   dsetup=maxSetup-setupAtToe;
   dsetdsta=dsetup/dsta;
   setup=setupAtToe+dsetdsta*(toe_sta-sta(1));
   sprintf('-!!- Location of SWEL-1.5*H0 is %4.1f ft landward of toe of slope',dsta)
   sprintf('-!!- Setup is interpolated between setup at toe of slope and max setup')
   sprintf('-!!-       setup is adjusted to %4.2f feet',setup)
   SWEL=SWEL-setupAtToe+setup;
   sprintf('-!!-       SWEL is adjusted to %4.2f feet',SWEL)
   k=find(dep < SWEL-1.5*H0)
   sta(k)=[];
   dep(k)=[];
else
   sprintf('-!!- The User has selected a starting point that is %4.2f feet above the elevation of SWEL-1.5H0\n',dep(1)-Ztoe)
   sprintf('-!!- This may be reasonable for some cases.  However the user may want to consider:\n')
   sprintf('-!!-   1) Selecting a starting point that is at or below %4.2f feet elevation, or\n', Ztoe)
   sprintf('-!!-   2) Reducing the incident wave height to a depth limited condition.\n')
end



% now iterate converge on a runup elevation
tol=0.01;  % convergence criteria
R2del=999;
R2_new=3*H0; %initial guess
R2=R2_new;
iter=0;
R2_all=[];
topStaAll=[];
Berm_Segs=[];
TAW_ALWAYS_VALID=1;
while(abs(R2del) > tol && iter <= 25)
    iter=iter+1;
    sprintf ('!----------------- STARTING ITERATION %d -----------------!',iter)
    % elevation of toe of slope
    Ztoe
    % station of toe slope (relative to 0-NAVD88 shoreline
    toe_sta 
    % station of top of slope/extent of 2% run-up
    top_sta
    % elevation of top of slope/extent of 2% run-up
    Z2
    % incident significant wave height
    H0
    % incident spectral peak wave period
    Tp
    % incident spectral mean wave period
    T0
        
    R2=R2_new  
    Z2=R2+SWEL
    % determine slope for this iteration
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end)
    end
    
    % get the length of the slope (not accounting for berm)
    Lslope=top_sta-toe_sta

    % loop over profile segments to determine berm factor 
    % re-calculate influence of depth of berm based on this run-up elevation
    % check for berm, berm width, berm height
    berm_width=0;
    rdh_sum=0;
    Berm_Segs=[];
    Berm_Heights=[];
    for kk=1:length(sta)-1
       ddep=dep(kk+1)-dep(kk);
       dsta=sta(kk+1)-sta(kk);
       s=ddep/dsta;
       if  (s < 1/15)     % count it as a berm if slope is flatter than 1:15 (see TAW manual)
          sprintf ('Berm Factor Calculation: Iteration %d, Profile Segment: %d',iter,kk) 
          berm_width=berm_width+dsta;   % tally the width of all berm segments
          % compute the rdh for this segment and weight it by the segment length
          dh=SWEL-(dep(kk)+dep(kk+1))/2  
          if  dh < 0
              chi=R2;
          else 
              chi=2* H0;
          end
          if (dh <= R2 & dh >=-2*H0)
             rdh=(0.5-0.5*cos(3.14159*dh/chi)) ;
          else
             rdh=1;
          end
          rdh_sum=rdh_sum + rdh * dsta
          Berm_Segs=[Berm_Segs, kk];
          Berm_Heights=[Berm_Heights, (dep(kk)+dep(kk+1))/2];
       end
       if dep(kk) >= Z2  % jump out of loop if we reached limit of run-up for this iteration 
          break
       end
    end
    sprintf ('!------- End Berm Factor Calculation, Iter: %d ---------!',iter)
    berm_width
    rB=berm_width/Lslope 
    if (berm_width > 0)
       rdh_mean=rdh_sum/berm_width
    else
       rdh_mean=1
    end
    gamma_berm=1- rB * (1-rdh_mean)
    if gamma_berm > 1
       gamma_berm=1
    end
    if gamma_berm < 0.6
       gamma_berm =0.6
    end
    % Iribarren number
    slope=(Z2-Ztoe)/(Lslope-berm_width)
    Irb=(slope/(sqrt(H0/L0)))
    % runup height
    gamma_berm
    gamma_perm
    gamma_beta
    gamma_rough
    gamma=gamma_berm*gamma_perm*gamma_beta*gamma_rough

    % check validity
    TAW_VALID=1;
    if (Irb*gamma_berm < 0.5 | Irb*gamma_berm > 10 )
       sprintf('!!! - - Iribaren number: %6.2f is outside the valid range (0.5-10), TAW NOT VALID - - !!!\n', Irb*gamma_berm)
       TAW_VALID=0;
    else
       sprintf('!!! - - Iribaren number: %6.2f is in the valid range (0.5-10), TAW RECOMMENDED - - !!!\n', Irb*gamma_berm)
    end
    islope=1/slope;  
    if (slope < 1/8 | slope > 1)
       sprintf('!!! - - slope: 1:%3.1f V:H is outside the valid range (1:8 - 1:1), TAW NOT VALID - - !!!\n', islope)
       TAW_VALID=0;
    else
       sprintf('!!! - - slope: 1:%3.1f V:H is in the valid range (1:8 - 1:1), TAW RECOMMENDED - - !!!\n', islope)
    end
    if TAW_VALID == 0
       TAW_ALWAYS_VALID=0;
    end
         
    if (Irb*gamma_berm < 1.8)
       R2_new=gamma*H0*1.77*Irb
    else
       R2_new=gamma*H0*(4.3-(1.6/sqrt(Irb))) 
    end


    % check to see if we need to evaluate a shallow foreshore
    if berm_width > 0.25 * L0;
       disp ('!   Berm_width is greater than 1/4 wave length')
       disp ('!   Runup will be weighted average with foreshore calculation assuming depth limited wave height on berm');
       % do the foreshore calculation 
       fore_H0=0.78*(SWEL_fore-min(Berm_Heights))
       % get upper slope
       fore_toe_sta=-999;
       fore_toe_dep=-999;
       for kk=length(dep)-1:-1:1
          ddep=dep(kk+1)-dep(kk);
          dsta=sta(kk+1)-sta(kk);
          s=ddep/dsta;
          if s < 1/15
             break
          end
          fore_toe_sta=sta(kk);
          fore_toe_dep=dep(kk);
          upper_slope=(Z2-fore_toe_dep)/(top_sta-fore_toe_sta)
       end
       fore_Irb=upper_slope/(sqrt(fore_H0/L0));
       fore_gamma=gamma_perm*gamma_beta*gamma_rough;
       if (fore_Irb < 1.8)
          fore_R2=fore_gamma*fore_H0*1.77*fore_Irb;
       else
          fore_R2=fore_gamma*fore_H0*(4.3-(1.6/sqrt(fore_Irb))); 
       end
       if berm_width >= L0
          R2_new=fore_R2
          disp ('berm is wider than one wavelength, use full shallow foreshore solution');
       else
          w2=(berm_width-0.25*L0)/(0.75*L0)
          w1=1-w2
          R2_new=w2*fore_R2 + w1*R2_new
       end
    end % end berm width check

    % convergence criterion
    R2del=abs(R2-R2_new)
    R2_all(iter)=R2_new;

    % get the new top station (for plot purposes)
    Z2=R2_new+SWEL
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end);
    end
    topStaAll(iter)=top_sta;

end

% final 2% runup elevation
Z2=R2_new+SWEL

diary off


% make some plots


% get the full transect



file=xlsread('../ADCIRC_returns/YK-99-2XYZSTA_RETURNS.csv');
sta_full=file(:,4);
ele_full=file(:,3);


%determine SWL station
SWL_sta=-999;
for kk=1:length(sta)-1
    if ((SWEL > dep(kk)) & (SWEL <= dep(kk+1)))    % here is the intersection of z2 with profile
       SWL_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),SWEL)
       break;
    end
end

figure('position', [100 100 900 900])
hold on
plot(sta_full,ele_full,'linewidth',2,'color',[0.5 0.2 0.1]);
plot(sta,dep,'md');
grid on
plot([min(sta),SWL_sta],[SWEL SWEL],'color','b','linewidth',2)
%plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
scatter(topStaAll,R2_all+SWEL,30,'r*');
%scatter(toe_sta,Db,70,'k*');
scatter(toe_sta,Ztoe,70,'k*');
% plot berm segs
for nn=1:length(Berm_Segs)
   seg=Berm_Segs(nn);
   h=Berm_Heights(nn);
   plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
end
xlim([toe_sta-50 max(topStaAll)+50])
grid minor
hleg=legend('Ground','TAW Segments','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
set (hleg,'location','southeast');
xlabel('Distance from zero NAVD88 shoreline (feet)');
ylabel('Elevation (Feet-NAVD88)');
str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
yylm=get (gca,'ylim');
dztxt=(yylm(end)-yylm(1))/15;
text(toe_sta,Z2+dztxt/4,str);
str=sprintf('Run-up Height is %5.1f Feet',R2_new);
text(toe_sta,Z2-dztxt/4,str);
title (plotTitle);
set(gcf,'color','w');
print('-r450','-dpdf',imgname);

-1.000000e+00
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S5');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U5');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O5');
clear all
close all
format long g

diary logfiles/YK-100-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Kennebunkport, York county, Maine
% TRANSECT ID: YK-100
% calculation by SJH, Ransom Consulting, Inc. 02-Apr-2020
% 100-year wave runup using TAW methodology
% including berm and weighted average with foreshore if necessary
%
% chk nld 20200220
%
% This script assumes that the incident wave conditions provided
% as input in the configuration section below are the 
% appropriate values located at the end of the foreshore
% or toe of the slope on which the run-up is being calculated
% the script does not attempt to apply a depth limit or any other
% transformation to the incident wave conditions other than 
% conversion of the peak wave period to the spectral mean wave 
% as recommended in the references below
% 
% references:
%
% Van der Meer, J.W., 2002. Technical Report Wave Run-up and 
% Wave Overtopping at Dikes. TAW Technical Advisory Committee on
% Flood Defence, The Netherlands.
%
% FEMA. 2007,  Atlantic Ocean and Gulf of Mexico Coastal Guidelines Update
%
% 
%

%--------------------------------------------------------------------------
% CONFIG
%--------------------------------------------------------------------------
fname='inpfiles/YK-100sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-100-runup';                            
SWEL=9.32;  % 100-yr still water level including wave setup. 
H0=4.7158;    % significant wave height at toe of structure 
Tp=11.1241;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.75;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.033077; 
maxSetup=0.77566;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-100'


% END CONFIG
%-------------------------------

SWEL=SWEL+setupAtToe  
SWEL_fore=SWEL+maxSetup


% FIND WAVELENGTH USING DEEPWATER DISPERSION RELATION
% using English units
L0=32.15/(2*pi)*T0^2

% Find Hb (Munk, 1949)
%Hb=H0/(3.3*(H0/L0)^(1/3)) 
%Db=-Hb/.78+SWEL;  % depth at breaking

% The toe elevation here is only used to determine the average 
% structure slope, it is not used to depth limit the wave height.
% Any depth limiting or other modification of the wave height 
% to make it consitent with TAW guidance should be performed 
% prior to the input of the significant wave height given above. 
Ztoe=SWEL-1.5*H0

% read the transect
[sta,dep,inc] = textread(fname,'%n%n%n%*[^\n]','delimiter',',','headerlines',0);

% remove unselected points
k=find(inc==0);
sta(k)=[];
dep(k)=[];

sta_org=sta;  % used for plotting purposes
dep_org=dep;


% initial guess at maximum run-up elevation to estimate slope
Z2=SWEL+1.5*H0

% determine station at the max runup and -1.5*H0 (i.e. the toe)
top_sta=-999;
toe_sta=-999;
for kk=1:length(sta)-1
    if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
       top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
    end
    if  ((Ztoe > dep(kk)) & (Ztoe <= dep(kk+1)))    % here is the intersection of Ztoe with profile
       toe_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Ztoe)
    end
end
% check to make sure we got them, if not extend the end slopes outward
S=diff(dep)./diff(sta);
if toe_sta==-999
   dy=dep(1)-Ztoe;
   toe_sta=sta(1)-dy/S(1)
end
if top_sta==-999
   dy=Z2-dep(end);
   top_sta=sta(end)+dy/S(end)
end

% just so the reader can tell the values aren't -999 anymore
top_sta
toe_sta



% check for case where the toe of slope is below SWL-1.5*H0
% in this case interpolate setup from the setupAtToe(really setup as first station), and the max setup
% also un-include points seaward of SWL-1.5*H0
if Ztoe > dep(1)
   dd=SWEL_fore-dep;
   k=find(dd<0,1); % k is index of first land point
   staAtSWL=interp1(dep(k-1:k),sta(k-1:k),SWEL_fore);   
   dsta=staAtSWL-sta(1);
   dsetup=maxSetup-setupAtToe;
   dsetdsta=dsetup/dsta;
   setup=setupAtToe+dsetdsta*(toe_sta-sta(1));
   sprintf('-!!- Location of SWEL-1.5*H0 is %4.1f ft landward of toe of slope',dsta)
   sprintf('-!!- Setup is interpolated between setup at toe of slope and max setup')
   sprintf('-!!-       setup is adjusted to %4.2f feet',setup)
   SWEL=SWEL-setupAtToe+setup;
   sprintf('-!!-       SWEL is adjusted to %4.2f feet',SWEL)
   k=find(dep < SWEL-1.5*H0)
   sta(k)=[];
   dep(k)=[];
else
   sprintf('-!!- The User has selected a starting point that is %4.2f feet above the elevation of SWEL-1.5H0\n',dep(1)-Ztoe)
   sprintf('-!!- This may be reasonable for some cases.  However the user may want to consider:\n')
   sprintf('-!!-   1) Selecting a starting point that is at or below %4.2f feet elevation, or\n', Ztoe)
   sprintf('-!!-   2) Reducing the incident wave height to a depth limited condition.\n')
end



% now iterate converge on a runup elevation
tol=0.01;  % convergence criteria
R2del=999;
R2_new=3*H0; %initial guess
R2=R2_new;
iter=0;
R2_all=[];
topStaAll=[];
Berm_Segs=[];
TAW_ALWAYS_VALID=1;
while(abs(R2del) > tol && iter <= 25)
    iter=iter+1;
    sprintf ('!----------------- STARTING ITERATION %d -----------------!',iter)
    % elevation of toe of slope
    Ztoe
    % station of toe slope (relative to 0-NAVD88 shoreline
    toe_sta 
    % station of top of slope/extent of 2% run-up
    top_sta
    % elevation of top of slope/extent of 2% run-up
    Z2
    % incident significant wave height
    H0
    % incident spectral peak wave period
    Tp
    % incident spectral mean wave period
    T0
        
    R2=R2_new  
    Z2=R2+SWEL
    % determine slope for this iteration
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end)
    end
    
    % get the length of the slope (not accounting for berm)
    Lslope=top_sta-toe_sta

    % loop over profile segments to determine berm factor 
    % re-calculate influence of depth of berm based on this run-up elevation
    % check for berm, berm width, berm height
    berm_width=0;
    rdh_sum=0;
    Berm_Segs=[];
    Berm_Heights=[];
    for kk=1:length(sta)-1
       ddep=dep(kk+1)-dep(kk);
       dsta=sta(kk+1)-sta(kk);
       s=ddep/dsta;
       if  (s < 1/15)     % count it as a berm if slope is flatter than 1:15 (see TAW manual)
          sprintf ('Berm Factor Calculation: Iteration %d, Profile Segment: %d',iter,kk) 
          berm_width=berm_width+dsta;   % tally the width of all berm segments
          % compute the rdh for this segment and weight it by the segment length
          dh=SWEL-(dep(kk)+dep(kk+1))/2  
          if  dh < 0
              chi=R2;
          else 
              chi=2* H0;
          end
          if (dh <= R2 & dh >=-2*H0)
             rdh=(0.5-0.5*cos(3.14159*dh/chi)) ;
          else
             rdh=1;
          end
          rdh_sum=rdh_sum + rdh * dsta
          Berm_Segs=[Berm_Segs, kk];
          Berm_Heights=[Berm_Heights, (dep(kk)+dep(kk+1))/2];
       end
       if dep(kk) >= Z2  % jump out of loop if we reached limit of run-up for this iteration 
          break
       end
    end
    sprintf ('!------- End Berm Factor Calculation, Iter: %d ---------!',iter)
    berm_width
    rB=berm_width/Lslope 
    if (berm_width > 0)
       rdh_mean=rdh_sum/berm_width
    else
       rdh_mean=1
    end
    gamma_berm=1- rB * (1-rdh_mean)
    if gamma_berm > 1
       gamma_berm=1
    end
    if gamma_berm < 0.6
       gamma_berm =0.6
    end
    % Iribarren number
    slope=(Z2-Ztoe)/(Lslope-berm_width)
    Irb=(slope/(sqrt(H0/L0)))
    % runup height
    gamma_berm
    gamma_perm
    gamma_beta
    gamma_rough
    gamma=gamma_berm*gamma_perm*gamma_beta*gamma_rough

    % check validity
    TAW_VALID=1;
    if (Irb*gamma_berm < 0.5 | Irb*gamma_berm > 10 )
       sprintf('!!! - - Iribaren number: %6.2f is outside the valid range (0.5-10), TAW NOT VALID - - !!!\n', Irb*gamma_berm)
       TAW_VALID=0;
    else
       sprintf('!!! - - Iribaren number: %6.2f is in the valid range (0.5-10), TAW RECOMMENDED - - !!!\n', Irb*gamma_berm)
    end
    islope=1/slope;  
    if (slope < 1/8 | slope > 1)
       sprintf('!!! - - slope: 1:%3.1f V:H is outside the valid range (1:8 - 1:1), TAW NOT VALID - - !!!\n', islope)
       TAW_VALID=0;
    else
       sprintf('!!! - - slope: 1:%3.1f V:H is in the valid range (1:8 - 1:1), TAW RECOMMENDED - - !!!\n', islope)
    end
    if TAW_VALID == 0
       TAW_ALWAYS_VALID=0;
    end
         
    if (Irb*gamma_berm < 1.8)
       R2_new=gamma*H0*1.77*Irb
    else
       R2_new=gamma*H0*(4.3-(1.6/sqrt(Irb))) 
    end


    % check to see if we need to evaluate a shallow foreshore
    if berm_width > 0.25 * L0;
       disp ('!   Berm_width is greater than 1/4 wave length')
       disp ('!   Runup will be weighted average with foreshore calculation assuming depth limited wave height on berm');
       % do the foreshore calculation 
       fore_H0=0.78*(SWEL_fore-min(Berm_Heights))
       % get upper slope
       fore_toe_sta=-999;
       fore_toe_dep=-999;
       for kk=length(dep)-1:-1:1
          ddep=dep(kk+1)-dep(kk);
          dsta=sta(kk+1)-sta(kk);
          s=ddep/dsta;
          if s < 1/15
             break
          end
          fore_toe_sta=sta(kk);
          fore_toe_dep=dep(kk);
          upper_slope=(Z2-fore_toe_dep)/(top_sta-fore_toe_sta)
       end
       fore_Irb=upper_slope/(sqrt(fore_H0/L0));
       fore_gamma=gamma_perm*gamma_beta*gamma_rough;
       if (fore_Irb < 1.8)
          fore_R2=fore_gamma*fore_H0*1.77*fore_Irb;
       else
          fore_R2=fore_gamma*fore_H0*(4.3-(1.6/sqrt(fore_Irb))); 
       end
       if berm_width >= L0
          R2_new=fore_R2
          disp ('berm is wider than one wavelength, use full shallow foreshore solution');
       else
          w2=(berm_width-0.25*L0)/(0.75*L0)
          w1=1-w2
          R2_new=w2*fore_R2 + w1*R2_new
       end
    end % end berm width check

    % convergence criterion
    R2del=abs(R2-R2_new)
    R2_all(iter)=R2_new;

    % get the new top station (for plot purposes)
    Z2=R2_new+SWEL
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end);
    end
    topStaAll(iter)=top_sta;

end

% final 2% runup elevation
Z2=R2_new+SWEL

diary off


% make some plots


% get the full transect



file=xlsread('../ADCIRC_returns/YK-100XYZSTA_RETURNS.csv');
sta_full=file(:,4);
ele_full=file(:,3);


%determine SWL station
SWL_sta=-999;
for kk=1:length(sta)-1
    if ((SWEL > dep(kk)) & (SWEL <= dep(kk+1)))    % here is the intersection of z2 with profile
       SWL_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),SWEL)
       break;
    end
end

figure('position', [100 100 900 900])
hold on
plot(sta_full,ele_full,'linewidth',2,'color',[0.5 0.2 0.1]);
plot(sta,dep,'md');
grid on
plot([min(sta),SWL_sta],[SWEL SWEL],'color','b','linewidth',2)
%plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
scatter(topStaAll,R2_all+SWEL,30,'r*');
%scatter(toe_sta,Db,70,'k*');
scatter(toe_sta,Ztoe,70,'k*');
% plot berm segs
for nn=1:length(Berm_Segs)
   seg=Berm_Segs(nn);
   h=Berm_Heights(nn);
   plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
end
xlim([toe_sta-50 max(topStaAll)+50])
grid minor
hleg=legend('Ground','TAW Segments','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
set (hleg,'location','southeast');
xlabel('Distance from zero NAVD88 shoreline (feet)');
ylabel('Elevation (Feet-NAVD88)');
str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
yylm=get (gca,'ylim');
dztxt=(yylm(end)-yylm(1))/15;
text(toe_sta,Z2+dztxt/4,str);
str=sprintf('Run-up Height is %5.1f Feet',R2_new);
text(toe_sta,Z2-dztxt/4,str);
title (plotTitle);
set(gcf,'color','w');
print('-r450','-dpdf',imgname);

-1.000000e+00
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S6');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U6');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O6');
clear all
close all
format long g

diary logfiles/YK-103-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Kennebunkport, York county, Maine
% TRANSECT ID: YK-103
% calculation by SJH, Ransom Consulting, Inc. 02-Apr-2020
% 100-year wave runup using TAW methodology
% including berm and weighted average with foreshore if necessary
%
% chk nld 20200220
%
% This script assumes that the incident wave conditions provided
% as input in the configuration section below are the 
% appropriate values located at the end of the foreshore
% or toe of the slope on which the run-up is being calculated
% the script does not attempt to apply a depth limit or any other
% transformation to the incident wave conditions other than 
% conversion of the peak wave period to the spectral mean wave 
% as recommended in the references below
% 
% references:
%
% Van der Meer, J.W., 2002. Technical Report Wave Run-up and 
% Wave Overtopping at Dikes. TAW Technical Advisory Committee on
% Flood Defence, The Netherlands.
%
% FEMA. 2007,  Atlantic Ocean and Gulf of Mexico Coastal Guidelines Update
%
% 
%

%--------------------------------------------------------------------------
% CONFIG
%--------------------------------------------------------------------------
fname='inpfiles/YK-103sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-103-runup';                            
SWEL=9.0423;  % 100-yr still water level including wave setup. 
H0=6.6846;    % significant wave height at toe of structure 
Tp=13.8007;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.85269;   % this may get changed automatically below
gamma_rough=0.75;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.96939; 
maxSetup=1.5373;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-103'


% END CONFIG
%-------------------------------

SWEL=SWEL+setupAtToe  
SWEL_fore=SWEL+maxSetup


% FIND WAVELENGTH USING DEEPWATER DISPERSION RELATION
% using English units
L0=32.15/(2*pi)*T0^2

% Find Hb (Munk, 1949)
%Hb=H0/(3.3*(H0/L0)^(1/3)) 
%Db=-Hb/.78+SWEL;  % depth at breaking

% The toe elevation here is only used to determine the average 
% structure slope, it is not used to depth limit the wave height.
% Any depth limiting or other modification of the wave height 
% to make it consitent with TAW guidance should be performed 
% prior to the input of the significant wave height given above. 
Ztoe=SWEL-1.5*H0

% read the transect
[sta,dep,inc] = textread(fname,'%n%n%n%*[^\n]','delimiter',',','headerlines',0);

% remove unselected points
k=find(inc==0);
sta(k)=[];
dep(k)=[];

sta_org=sta;  % used for plotting purposes
dep_org=dep;


% initial guess at maximum run-up elevation to estimate slope
Z2=SWEL+1.5*H0

% determine station at the max runup and -1.5*H0 (i.e. the toe)
top_sta=-999;
toe_sta=-999;
for kk=1:length(sta)-1
    if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
       top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
    end
    if  ((Ztoe > dep(kk)) & (Ztoe <= dep(kk+1)))    % here is the intersection of Ztoe with profile
       toe_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Ztoe)
    end
end
% check to make sure we got them, if not extend the end slopes outward
S=diff(dep)./diff(sta);
if toe_sta==-999
   dy=dep(1)-Ztoe;
   toe_sta=sta(1)-dy/S(1)
end
if top_sta==-999
   dy=Z2-dep(end);
   top_sta=sta(end)+dy/S(end)
end

% just so the reader can tell the values aren't -999 anymore
top_sta
toe_sta



% check for case where the toe of slope is below SWL-1.5*H0
% in this case interpolate setup from the setupAtToe(really setup as first station), and the max setup
% also un-include points seaward of SWL-1.5*H0
if Ztoe > dep(1)
   dd=SWEL_fore-dep;
   k=find(dd<0,1); % k is index of first land point
   staAtSWL=interp1(dep(k-1:k),sta(k-1:k),SWEL_fore);   
   dsta=staAtSWL-sta(1);
   dsetup=maxSetup-setupAtToe;
   dsetdsta=dsetup/dsta;
   setup=setupAtToe+dsetdsta*(toe_sta-sta(1));
   sprintf('-!!- Location of SWEL-1.5*H0 is %4.1f ft landward of toe of slope',dsta)
   sprintf('-!!- Setup is interpolated between setup at toe of slope and max setup')
   sprintf('-!!-       setup is adjusted to %4.2f feet',setup)
   SWEL=SWEL-setupAtToe+setup;
   sprintf('-!!-       SWEL is adjusted to %4.2f feet',SWEL)
   k=find(dep < SWEL-1.5*H0)
   sta(k)=[];
   dep(k)=[];
else
   sprintf('-!!- The User has selected a starting point that is %4.2f feet above the elevation of SWEL-1.5H0\n',dep(1)-Ztoe)
   sprintf('-!!- This may be reasonable for some cases.  However the user may want to consider:\n')
   sprintf('-!!-   1) Selecting a starting point that is at or below %4.2f feet elevation, or\n', Ztoe)
   sprintf('-!!-   2) Reducing the incident wave height to a depth limited condition.\n')
end



% now iterate converge on a runup elevation
tol=0.01;  % convergence criteria
R2del=999;
R2_new=3*H0; %initial guess
R2=R2_new;
iter=0;
R2_all=[];
topStaAll=[];
Berm_Segs=[];
TAW_ALWAYS_VALID=1;
while(abs(R2del) > tol && iter <= 25)
    iter=iter+1;
    sprintf ('!----------------- STARTING ITERATION %d -----------------!',iter)
    % elevation of toe of slope
    Ztoe
    % station of toe slope (relative to 0-NAVD88 shoreline
    toe_sta 
    % station of top of slope/extent of 2% run-up
    top_sta
    % elevation of top of slope/extent of 2% run-up
    Z2
    % incident significant wave height
    H0
    % incident spectral peak wave period
    Tp
    % incident spectral mean wave period
    T0
        
    R2=R2_new  
    Z2=R2+SWEL
    % determine slope for this iteration
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end)
    end
    
    % get the length of the slope (not accounting for berm)
    Lslope=top_sta-toe_sta

    % loop over profile segments to determine berm factor 
    % re-calculate influence of depth of berm based on this run-up elevation
    % check for berm, berm width, berm height
    berm_width=0;
    rdh_sum=0;
    Berm_Segs=[];
    Berm_Heights=[];
    for kk=1:length(sta)-1
       ddep=dep(kk+1)-dep(kk);
       dsta=sta(kk+1)-sta(kk);
       s=ddep/dsta;
       if  (s < 1/15)     % count it as a berm if slope is flatter than 1:15 (see TAW manual)
          sprintf ('Berm Factor Calculation: Iteration %d, Profile Segment: %d',iter,kk) 
          berm_width=berm_width+dsta;   % tally the width of all berm segments
          % compute the rdh for this segment and weight it by the segment length
          dh=SWEL-(dep(kk)+dep(kk+1))/2  
          if  dh < 0
              chi=R2;
          else 
              chi=2* H0;
          end
          if (dh <= R2 & dh >=-2*H0)
             rdh=(0.5-0.5*cos(3.14159*dh/chi)) ;
          else
             rdh=1;
          end
          rdh_sum=rdh_sum + rdh * dsta
          Berm_Segs=[Berm_Segs, kk];
          Berm_Heights=[Berm_Heights, (dep(kk)+dep(kk+1))/2];
       end
       if dep(kk) >= Z2  % jump out of loop if we reached limit of run-up for this iteration 
          break
       end
    end
    sprintf ('!------- End Berm Factor Calculation, Iter: %d ---------!',iter)
    berm_width
    rB=berm_width/Lslope 
    if (berm_width > 0)
       rdh_mean=rdh_sum/berm_width
    else
       rdh_mean=1
    end
    gamma_berm=1- rB * (1-rdh_mean)
    if gamma_berm > 1
       gamma_berm=1
    end
    if gamma_berm < 0.6
       gamma_berm =0.6
    end
    % Iribarren number
    slope=(Z2-Ztoe)/(Lslope-berm_width)
    Irb=(slope/(sqrt(H0/L0)))
    % runup height
    gamma_berm
    gamma_perm
    gamma_beta
    gamma_rough
    gamma=gamma_berm*gamma_perm*gamma_beta*gamma_rough

    % check validity
    TAW_VALID=1;
    if (Irb*gamma_berm < 0.5 | Irb*gamma_berm > 10 )
       sprintf('!!! - - Iribaren number: %6.2f is outside the valid range (0.5-10), TAW NOT VALID - - !!!\n', Irb*gamma_berm)
       TAW_VALID=0;
    else
       sprintf('!!! - - Iribaren number: %6.2f is in the valid range (0.5-10), TAW RECOMMENDED - - !!!\n', Irb*gamma_berm)
    end
    islope=1/slope;  
    if (slope < 1/8 | slope > 1)
       sprintf('!!! - - slope: 1:%3.1f V:H is outside the valid range (1:8 - 1:1), TAW NOT VALID - - !!!\n', islope)
       TAW_VALID=0;
    else
       sprintf('!!! - - slope: 1:%3.1f V:H is in the valid range (1:8 - 1:1), TAW RECOMMENDED - - !!!\n', islope)
    end
    if TAW_VALID == 0
       TAW_ALWAYS_VALID=0;
    end
         
    if (Irb*gamma_berm < 1.8)
       R2_new=gamma*H0*1.77*Irb
    else
       R2_new=gamma*H0*(4.3-(1.6/sqrt(Irb))) 
    end


    % check to see if we need to evaluate a shallow foreshore
    if berm_width > 0.25 * L0;
       disp ('!   Berm_width is greater than 1/4 wave length')
       disp ('!   Runup will be weighted average with foreshore calculation assuming depth limited wave height on berm');
       % do the foreshore calculation 
       fore_H0=0.78*(SWEL_fore-min(Berm_Heights))
       % get upper slope
       fore_toe_sta=-999;
       fore_toe_dep=-999;
       for kk=length(dep)-1:-1:1
          ddep=dep(kk+1)-dep(kk);
          dsta=sta(kk+1)-sta(kk);
          s=ddep/dsta;
          if s < 1/15
             break
          end
          fore_toe_sta=sta(kk);
          fore_toe_dep=dep(kk);
          upper_slope=(Z2-fore_toe_dep)/(top_sta-fore_toe_sta)
       end
       fore_Irb=upper_slope/(sqrt(fore_H0/L0));
       fore_gamma=gamma_perm*gamma_beta*gamma_rough;
       if (fore_Irb < 1.8)
          fore_R2=fore_gamma*fore_H0*1.77*fore_Irb;
       else
          fore_R2=fore_gamma*fore_H0*(4.3-(1.6/sqrt(fore_Irb))); 
       end
       if berm_width >= L0
          R2_new=fore_R2
          disp ('berm is wider than one wavelength, use full shallow foreshore solution');
       else
          w2=(berm_width-0.25*L0)/(0.75*L0)
          w1=1-w2
          R2_new=w2*fore_R2 + w1*R2_new
       end
    end % end berm width check

    % convergence criterion
    R2del=abs(R2-R2_new)
    R2_all(iter)=R2_new;

    % get the new top station (for plot purposes)
    Z2=R2_new+SWEL
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end);
    end
    topStaAll(iter)=top_sta;

end

% final 2% runup elevation
Z2=R2_new+SWEL

diary off


% make some plots


% get the full transect



file=xlsread('../ADCIRC_returns/YK-103XYZSTA_RETURNS.csv');
sta_full=file(:,4);
ele_full=file(:,3);


%determine SWL station
SWL_sta=-999;
for kk=1:length(sta)-1
    if ((SWEL > dep(kk)) & (SWEL <= dep(kk+1)))    % here is the intersection of z2 with profile
       SWL_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),SWEL)
       break;
    end
end

figure('position', [100 100 900 900])
hold on
plot(sta_full,ele_full,'linewidth',2,'color',[0.5 0.2 0.1]);
plot(sta,dep,'md');
grid on
plot([min(sta),SWL_sta],[SWEL SWEL],'color','b','linewidth',2)
%plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
scatter(topStaAll,R2_all+SWEL,30,'r*');
%scatter(toe_sta,Db,70,'k*');
scatter(toe_sta,Ztoe,70,'k*');
% plot berm segs
for nn=1:length(Berm_Segs)
   seg=Berm_Segs(nn);
   h=Berm_Heights(nn);
   plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
end
xlim([toe_sta-50 max(topStaAll)+50])
grid minor
hleg=legend('Ground','TAW Segments','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
set (hleg,'location','southeast');
xlabel('Distance from zero NAVD88 shoreline (feet)');
ylabel('Elevation (Feet-NAVD88)');
str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
yylm=get (gca,'ylim');
dztxt=(yylm(end)-yylm(1))/15;
text(toe_sta,Z2+dztxt/4,str);
str=sprintf('Run-up Height is %5.1f Feet',R2_new);
text(toe_sta,Z2-dztxt/4,str);
title (plotTitle);
set(gcf,'color','w');
print('-r450','-dpdf',imgname);

-1.000000e+00
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S7');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U7');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O7');
clear all
close all
format long g

diary logfiles/YK-105-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Kennebunkport, York county, Maine
% TRANSECT ID: YK-105
% calculation by SJH, Ransom Consulting, Inc. 02-Apr-2020
% 100-year wave runup using TAW methodology
% including berm and weighted average with foreshore if necessary
%
% chk nld 20200220
%
% This script assumes that the incident wave conditions provided
% as input in the configuration section below are the 
% appropriate values located at the end of the foreshore
% or toe of the slope on which the run-up is being calculated
% the script does not attempt to apply a depth limit or any other
% transformation to the incident wave conditions other than 
% conversion of the peak wave period to the spectral mean wave 
% as recommended in the references below
% 
% references:
%
% Van der Meer, J.W., 2002. Technical Report Wave Run-up and 
% Wave Overtopping at Dikes. TAW Technical Advisory Committee on
% Flood Defence, The Netherlands.
%
% FEMA. 2007,  Atlantic Ocean and Gulf of Mexico Coastal Guidelines Update
%
% 
%

%--------------------------------------------------------------------------
% CONFIG
%--------------------------------------------------------------------------
fname='inpfiles/YK-105sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-105-runup';                            
SWEL=9.4073;  % 100-yr still water level including wave setup. 
H0=3.7367;    % significant wave height at toe of structure 
Tp=12.3575;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.8768;   % this may get changed automatically below
gamma_rough=0.6;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.10581; 
maxSetup=0.56123;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-105'


% END CONFIG
%-------------------------------

SWEL=SWEL+setupAtToe  
SWEL_fore=SWEL+maxSetup


% FIND WAVELENGTH USING DEEPWATER DISPERSION RELATION
% using English units
L0=32.15/(2*pi)*T0^2

% Find Hb (Munk, 1949)
%Hb=H0/(3.3*(H0/L0)^(1/3)) 
%Db=-Hb/.78+SWEL;  % depth at breaking

% The toe elevation here is only used to determine the average 
% structure slope, it is not used to depth limit the wave height.
% Any depth limiting or other modification of the wave height 
% to make it consitent with TAW guidance should be performed 
% prior to the input of the significant wave height given above. 
Ztoe=SWEL-1.5*H0

% read the transect
[sta,dep,inc] = textread(fname,'%n%n%n%*[^\n]','delimiter',',','headerlines',0);

% remove unselected points
k=find(inc==0);
sta(k)=[];
dep(k)=[];

sta_org=sta;  % used for plotting purposes
dep_org=dep;


% initial guess at maximum run-up elevation to estimate slope
Z2=SWEL+1.5*H0

% determine station at the max runup and -1.5*H0 (i.e. the toe)
top_sta=-999;
toe_sta=-999;
for kk=1:length(sta)-1
    if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
       top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
    end
    if  ((Ztoe > dep(kk)) & (Ztoe <= dep(kk+1)))    % here is the intersection of Ztoe with profile
       toe_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Ztoe)
    end
end
% check to make sure we got them, if not extend the end slopes outward
S=diff(dep)./diff(sta);
if toe_sta==-999
   dy=dep(1)-Ztoe;
   toe_sta=sta(1)-dy/S(1)
end
if top_sta==-999
   dy=Z2-dep(end);
   top_sta=sta(end)+dy/S(end)
end

% just so the reader can tell the values aren't -999 anymore
top_sta
toe_sta



% check for case where the toe of slope is below SWL-1.5*H0
% in this case interpolate setup from the setupAtToe(really setup as first station), and the max setup
% also un-include points seaward of SWL-1.5*H0
if Ztoe > dep(1)
   dd=SWEL_fore-dep;
   k=find(dd<0,1); % k is index of first land point
   staAtSWL=interp1(dep(k-1:k),sta(k-1:k),SWEL_fore);   
   dsta=staAtSWL-sta(1);
   dsetup=maxSetup-setupAtToe;
   dsetdsta=dsetup/dsta;
   setup=setupAtToe+dsetdsta*(toe_sta-sta(1));
   sprintf('-!!- Location of SWEL-1.5*H0 is %4.1f ft landward of toe of slope',dsta)
   sprintf('-!!- Setup is interpolated between setup at toe of slope and max setup')
   sprintf('-!!-       setup is adjusted to %4.2f feet',setup)
   SWEL=SWEL-setupAtToe+setup;
   sprintf('-!!-       SWEL is adjusted to %4.2f feet',SWEL)
   k=find(dep < SWEL-1.5*H0)
   sta(k)=[];
   dep(k)=[];
else
   sprintf('-!!- The User has selected a starting point that is %4.2f feet above the elevation of SWEL-1.5H0\n',dep(1)-Ztoe)
   sprintf('-!!- This may be reasonable for some cases.  However the user may want to consider:\n')
   sprintf('-!!-   1) Selecting a starting point that is at or below %4.2f feet elevation, or\n', Ztoe)
   sprintf('-!!-   2) Reducing the incident wave height to a depth limited condition.\n')
end



% now iterate converge on a runup elevation
tol=0.01;  % convergence criteria
R2del=999;
R2_new=3*H0; %initial guess
R2=R2_new;
iter=0;
R2_all=[];
topStaAll=[];
Berm_Segs=[];
TAW_ALWAYS_VALID=1;
while(abs(R2del) > tol && iter <= 25)
    iter=iter+1;
    sprintf ('!----------------- STARTING ITERATION %d -----------------!',iter)
    % elevation of toe of slope
    Ztoe
    % station of toe slope (relative to 0-NAVD88 shoreline
    toe_sta 
    % station of top of slope/extent of 2% run-up
    top_sta
    % elevation of top of slope/extent of 2% run-up
    Z2
    % incident significant wave height
    H0
    % incident spectral peak wave period
    Tp
    % incident spectral mean wave period
    T0
        
    R2=R2_new  
    Z2=R2+SWEL
    % determine slope for this iteration
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end)
    end
    
    % get the length of the slope (not accounting for berm)
    Lslope=top_sta-toe_sta

    % loop over profile segments to determine berm factor 
    % re-calculate influence of depth of berm based on this run-up elevation
    % check for berm, berm width, berm height
    berm_width=0;
    rdh_sum=0;
    Berm_Segs=[];
    Berm_Heights=[];
    for kk=1:length(sta)-1
       ddep=dep(kk+1)-dep(kk);
       dsta=sta(kk+1)-sta(kk);
       s=ddep/dsta;
       if  (s < 1/15)     % count it as a berm if slope is flatter than 1:15 (see TAW manual)
          sprintf ('Berm Factor Calculation: Iteration %d, Profile Segment: %d',iter,kk) 
          berm_width=berm_width+dsta;   % tally the width of all berm segments
          % compute the rdh for this segment and weight it by the segment length
          dh=SWEL-(dep(kk)+dep(kk+1))/2  
          if  dh < 0
              chi=R2;
          else 
              chi=2* H0;
          end
          if (dh <= R2 & dh >=-2*H0)
             rdh=(0.5-0.5*cos(3.14159*dh/chi)) ;
          else
             rdh=1;
          end
          rdh_sum=rdh_sum + rdh * dsta
          Berm_Segs=[Berm_Segs, kk];
          Berm_Heights=[Berm_Heights, (dep(kk)+dep(kk+1))/2];
       end
       if dep(kk) >= Z2  % jump out of loop if we reached limit of run-up for this iteration 
          break
       end
    end
    sprintf ('!------- End Berm Factor Calculation, Iter: %d ---------!',iter)
    berm_width
    rB=berm_width/Lslope 
    if (berm_width > 0)
       rdh_mean=rdh_sum/berm_width
    else
       rdh_mean=1
    end
    gamma_berm=1- rB * (1-rdh_mean)
    if gamma_berm > 1
       gamma_berm=1
    end
    if gamma_berm < 0.6
       gamma_berm =0.6
    end
    % Iribarren number
    slope=(Z2-Ztoe)/(Lslope-berm_width)
    Irb=(slope/(sqrt(H0/L0)))
    % runup height
    gamma_berm
    gamma_perm
    gamma_beta
    gamma_rough
    gamma=gamma_berm*gamma_perm*gamma_beta*gamma_rough

    % check validity
    TAW_VALID=1;
    if (Irb*gamma_berm < 0.5 | Irb*gamma_berm > 10 )
       sprintf('!!! - - Iribaren number: %6.2f is outside the valid range (0.5-10), TAW NOT VALID - - !!!\n', Irb*gamma_berm)
       TAW_VALID=0;
    else
       sprintf('!!! - - Iribaren number: %6.2f is in the valid range (0.5-10), TAW RECOMMENDED - - !!!\n', Irb*gamma_berm)
    end
    islope=1/slope;  
    if (slope < 1/8 | slope > 1)
       sprintf('!!! - - slope: 1:%3.1f V:H is outside the valid range (1:8 - 1:1), TAW NOT VALID - - !!!\n', islope)
       TAW_VALID=0;
    else
       sprintf('!!! - - slope: 1:%3.1f V:H is in the valid range (1:8 - 1:1), TAW RECOMMENDED - - !!!\n', islope)
    end
    if TAW_VALID == 0
       TAW_ALWAYS_VALID=0;
    end
         
    if (Irb*gamma_berm < 1.8)
       R2_new=gamma*H0*1.77*Irb
    else
       R2_new=gamma*H0*(4.3-(1.6/sqrt(Irb))) 
    end


    % check to see if we need to evaluate a shallow foreshore
    if berm_width > 0.25 * L0;
       disp ('!   Berm_width is greater than 1/4 wave length')
       disp ('!   Runup will be weighted average with foreshore calculation assuming depth limited wave height on berm');
       % do the foreshore calculation 
       fore_H0=0.78*(SWEL_fore-min(Berm_Heights))
       % get upper slope
       fore_toe_sta=-999;
       fore_toe_dep=-999;
       for kk=length(dep)-1:-1:1
          ddep=dep(kk+1)-dep(kk);
          dsta=sta(kk+1)-sta(kk);
          s=ddep/dsta;
          if s < 1/15
             break
          end
          fore_toe_sta=sta(kk);
          fore_toe_dep=dep(kk);
          upper_slope=(Z2-fore_toe_dep)/(top_sta-fore_toe_sta)
       end
       fore_Irb=upper_slope/(sqrt(fore_H0/L0));
       fore_gamma=gamma_perm*gamma_beta*gamma_rough;
       if (fore_Irb < 1.8)
          fore_R2=fore_gamma*fore_H0*1.77*fore_Irb;
       else
          fore_R2=fore_gamma*fore_H0*(4.3-(1.6/sqrt(fore_Irb))); 
       end
       if berm_width >= L0
          R2_new=fore_R2
          disp ('berm is wider than one wavelength, use full shallow foreshore solution');
       else
          w2=(berm_width-0.25*L0)/(0.75*L0)
          w1=1-w2
          R2_new=w2*fore_R2 + w1*R2_new
       end
    end % end berm width check

    % convergence criterion
    R2del=abs(R2-R2_new)
    R2_all(iter)=R2_new;

    % get the new top station (for plot purposes)
    Z2=R2_new+SWEL
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end);
    end
    topStaAll(iter)=top_sta;

end

% final 2% runup elevation
Z2=R2_new+SWEL

diary off


% make some plots


% get the full transect



file=xlsread('../ADCIRC_returns/YK-105XYZSTA_RETURNS.csv');
sta_full=file(:,4);
ele_full=file(:,3);


%determine SWL station
SWL_sta=-999;
for kk=1:length(sta)-1
    if ((SWEL > dep(kk)) & (SWEL <= dep(kk+1)))    % here is the intersection of z2 with profile
       SWL_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),SWEL)
       break;
    end
end

figure('position', [100 100 900 900])
hold on
plot(sta_full,ele_full,'linewidth',2,'color',[0.5 0.2 0.1]);
plot(sta,dep,'md');
grid on
plot([min(sta),SWL_sta],[SWEL SWEL],'color','b','linewidth',2)
%plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
scatter(topStaAll,R2_all+SWEL,30,'r*');
%scatter(toe_sta,Db,70,'k*');
scatter(toe_sta,Ztoe,70,'k*');
% plot berm segs
for nn=1:length(Berm_Segs)
   seg=Berm_Segs(nn);
   h=Berm_Heights(nn);
   plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
end
xlim([toe_sta-50 max(topStaAll)+50])
grid minor
hleg=legend('Ground','TAW Segments','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
set (hleg,'location','southeast');
xlabel('Distance from zero NAVD88 shoreline (feet)');
ylabel('Elevation (Feet-NAVD88)');
str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
yylm=get (gca,'ylim');
dztxt=(yylm(end)-yylm(1))/15;
text(toe_sta,Z2+dztxt/4,str);
str=sprintf('Run-up Height is %5.1f Feet',R2_new);
text(toe_sta,Z2-dztxt/4,str);
title (plotTitle);
set(gcf,'color','w');
print('-r450','-dpdf',imgname);

-1.000000e+00
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S8');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U8');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O8');
clear all
close all
format long g

diary logfiles/YK-106-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Kennebunkport, York county, Maine
% TRANSECT ID: YK-106
% calculation by SJH, Ransom Consulting, Inc. 02-Apr-2020
% 100-year wave runup using TAW methodology
% including berm and weighted average with foreshore if necessary
%
% chk nld 20200220
%
% This script assumes that the incident wave conditions provided
% as input in the configuration section below are the 
% appropriate values located at the end of the foreshore
% or toe of the slope on which the run-up is being calculated
% the script does not attempt to apply a depth limit or any other
% transformation to the incident wave conditions other than 
% conversion of the peak wave period to the spectral mean wave 
% as recommended in the references below
% 
% references:
%
% Van der Meer, J.W., 2002. Technical Report Wave Run-up and 
% Wave Overtopping at Dikes. TAW Technical Advisory Committee on
% Flood Defence, The Netherlands.
%
% FEMA. 2007,  Atlantic Ocean and Gulf of Mexico Coastal Guidelines Update
%
% 
%

%--------------------------------------------------------------------------
% CONFIG
%--------------------------------------------------------------------------
fname='inpfiles/YK-106sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-106-runup';                            
SWEL=9.3612;  % 100-yr still water level including wave setup. 
H0=2.2239;    % significant wave height at toe of structure 
Tp=12.5841;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.94278;   % this may get changed automatically below
gamma_rough=0.6;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.66269; 
maxSetup=0.86394;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-106'


% END CONFIG
%-------------------------------

SWEL=SWEL+setupAtToe  
SWEL_fore=SWEL+maxSetup


% FIND WAVELENGTH USING DEEPWATER DISPERSION RELATION
% using English units
L0=32.15/(2*pi)*T0^2

% Find Hb (Munk, 1949)
%Hb=H0/(3.3*(H0/L0)^(1/3)) 
%Db=-Hb/.78+SWEL;  % depth at breaking

% The toe elevation here is only used to determine the average 
% structure slope, it is not used to depth limit the wave height.
% Any depth limiting or other modification of the wave height 
% to make it consitent with TAW guidance should be performed 
% prior to the input of the significant wave height given above. 
Ztoe=SWEL-1.5*H0

% read the transect
[sta,dep,inc] = textread(fname,'%n%n%n%*[^\n]','delimiter',',','headerlines',0);

% remove unselected points
k=find(inc==0);
sta(k)=[];
dep(k)=[];

sta_org=sta;  % used for plotting purposes
dep_org=dep;


% initial guess at maximum run-up elevation to estimate slope
Z2=SWEL+1.5*H0

% determine station at the max runup and -1.5*H0 (i.e. the toe)
top_sta=-999;
toe_sta=-999;
for kk=1:length(sta)-1
    if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
       top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
    end
    if  ((Ztoe > dep(kk)) & (Ztoe <= dep(kk+1)))    % here is the intersection of Ztoe with profile
       toe_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Ztoe)
    end
end
% check to make sure we got them, if not extend the end slopes outward
S=diff(dep)./diff(sta);
if toe_sta==-999
   dy=dep(1)-Ztoe;
   toe_sta=sta(1)-dy/S(1)
end
if top_sta==-999
   dy=Z2-dep(end);
   top_sta=sta(end)+dy/S(end)
end

% just so the reader can tell the values aren't -999 anymore
top_sta
toe_sta



% check for case where the toe of slope is below SWL-1.5*H0
% in this case interpolate setup from the setupAtToe(really setup as first station), and the max setup
% also un-include points seaward of SWL-1.5*H0
if Ztoe > dep(1)
   dd=SWEL_fore-dep;
   k=find(dd<0,1); % k is index of first land point
   staAtSWL=interp1(dep(k-1:k),sta(k-1:k),SWEL_fore);   
   dsta=staAtSWL-sta(1);
   dsetup=maxSetup-setupAtToe;
   dsetdsta=dsetup/dsta;
   setup=setupAtToe+dsetdsta*(toe_sta-sta(1));
   sprintf('-!!- Location of SWEL-1.5*H0 is %4.1f ft landward of toe of slope',dsta)
   sprintf('-!!- Setup is interpolated between setup at toe of slope and max setup')
   sprintf('-!!-       setup is adjusted to %4.2f feet',setup)
   SWEL=SWEL-setupAtToe+setup;
   sprintf('-!!-       SWEL is adjusted to %4.2f feet',SWEL)
   k=find(dep < SWEL-1.5*H0)
   sta(k)=[];
   dep(k)=[];
else
   sprintf('-!!- The User has selected a starting point that is %4.2f feet above the elevation of SWEL-1.5H0\n',dep(1)-Ztoe)
   sprintf('-!!- This may be reasonable for some cases.  However the user may want to consider:\n')
   sprintf('-!!-   1) Selecting a starting point that is at or below %4.2f feet elevation, or\n', Ztoe)
   sprintf('-!!-   2) Reducing the incident wave height to a depth limited condition.\n')
end



% now iterate converge on a runup elevation
tol=0.01;  % convergence criteria
R2del=999;
R2_new=3*H0; %initial guess
R2=R2_new;
iter=0;
R2_all=[];
topStaAll=[];
Berm_Segs=[];
TAW_ALWAYS_VALID=1;
while(abs(R2del) > tol && iter <= 25)
    iter=iter+1;
    sprintf ('!----------------- STARTING ITERATION %d -----------------!',iter)
    % elevation of toe of slope
    Ztoe
    % station of toe slope (relative to 0-NAVD88 shoreline
    toe_sta 
    % station of top of slope/extent of 2% run-up
    top_sta
    % elevation of top of slope/extent of 2% run-up
    Z2
    % incident significant wave height
    H0
    % incident spectral peak wave period
    Tp
    % incident spectral mean wave period
    T0
        
    R2=R2_new  
    Z2=R2+SWEL
    % determine slope for this iteration
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end)
    end
    
    % get the length of the slope (not accounting for berm)
    Lslope=top_sta-toe_sta

    % loop over profile segments to determine berm factor 
    % re-calculate influence of depth of berm based on this run-up elevation
    % check for berm, berm width, berm height
    berm_width=0;
    rdh_sum=0;
    Berm_Segs=[];
    Berm_Heights=[];
    for kk=1:length(sta)-1
       ddep=dep(kk+1)-dep(kk);
       dsta=sta(kk+1)-sta(kk);
       s=ddep/dsta;
       if  (s < 1/15)     % count it as a berm if slope is flatter than 1:15 (see TAW manual)
          sprintf ('Berm Factor Calculation: Iteration %d, Profile Segment: %d',iter,kk) 
          berm_width=berm_width+dsta;   % tally the width of all berm segments
          % compute the rdh for this segment and weight it by the segment length
          dh=SWEL-(dep(kk)+dep(kk+1))/2  
          if  dh < 0
              chi=R2;
          else 
              chi=2* H0;
          end
          if (dh <= R2 & dh >=-2*H0)
             rdh=(0.5-0.5*cos(3.14159*dh/chi)) ;
          else
             rdh=1;
          end
          rdh_sum=rdh_sum + rdh * dsta
          Berm_Segs=[Berm_Segs, kk];
          Berm_Heights=[Berm_Heights, (dep(kk)+dep(kk+1))/2];
       end
       if dep(kk) >= Z2  % jump out of loop if we reached limit of run-up for this iteration 
          break
       end
    end
    sprintf ('!------- End Berm Factor Calculation, Iter: %d ---------!',iter)
    berm_width
    rB=berm_width/Lslope 
    if (berm_width > 0)
       rdh_mean=rdh_sum/berm_width
    else
       rdh_mean=1
    end
    gamma_berm=1- rB * (1-rdh_mean)
    if gamma_berm > 1
       gamma_berm=1
    end
    if gamma_berm < 0.6
       gamma_berm =0.6
    end
    % Iribarren number
    slope=(Z2-Ztoe)/(Lslope-berm_width)
    Irb=(slope/(sqrt(H0/L0)))
    % runup height
    gamma_berm
    gamma_perm
    gamma_beta
    gamma_rough
    gamma=gamma_berm*gamma_perm*gamma_beta*gamma_rough

    % check validity
    TAW_VALID=1;
    if (Irb*gamma_berm < 0.5 | Irb*gamma_berm > 10 )
       sprintf('!!! - - Iribaren number: %6.2f is outside the valid range (0.5-10), TAW NOT VALID - - !!!\n', Irb*gamma_berm)
       TAW_VALID=0;
    else
       sprintf('!!! - - Iribaren number: %6.2f is in the valid range (0.5-10), TAW RECOMMENDED - - !!!\n', Irb*gamma_berm)
    end
    islope=1/slope;  
    if (slope < 1/8 | slope > 1)
       sprintf('!!! - - slope: 1:%3.1f V:H is outside the valid range (1:8 - 1:1), TAW NOT VALID - - !!!\n', islope)
       TAW_VALID=0;
    else
       sprintf('!!! - - slope: 1:%3.1f V:H is in the valid range (1:8 - 1:1), TAW RECOMMENDED - - !!!\n', islope)
    end
    if TAW_VALID == 0
       TAW_ALWAYS_VALID=0;
    end
         
    if (Irb*gamma_berm < 1.8)
       R2_new=gamma*H0*1.77*Irb
    else
       R2_new=gamma*H0*(4.3-(1.6/sqrt(Irb))) 
    end


    % check to see if we need to evaluate a shallow foreshore
    if berm_width > 0.25 * L0;
       disp ('!   Berm_width is greater than 1/4 wave length')
       disp ('!   Runup will be weighted average with foreshore calculation assuming depth limited wave height on berm');
       % do the foreshore calculation 
       fore_H0=0.78*(SWEL_fore-min(Berm_Heights))
       % get upper slope
       fore_toe_sta=-999;
       fore_toe_dep=-999;
       for kk=length(dep)-1:-1:1
          ddep=dep(kk+1)-dep(kk);
          dsta=sta(kk+1)-sta(kk);
          s=ddep/dsta;
          if s < 1/15
             break
          end
          fore_toe_sta=sta(kk);
          fore_toe_dep=dep(kk);
          upper_slope=(Z2-fore_toe_dep)/(top_sta-fore_toe_sta)
       end
       fore_Irb=upper_slope/(sqrt(fore_H0/L0));
       fore_gamma=gamma_perm*gamma_beta*gamma_rough;
       if (fore_Irb < 1.8)
          fore_R2=fore_gamma*fore_H0*1.77*fore_Irb;
       else
          fore_R2=fore_gamma*fore_H0*(4.3-(1.6/sqrt(fore_Irb))); 
       end
       if berm_width >= L0
          R2_new=fore_R2
          disp ('berm is wider than one wavelength, use full shallow foreshore solution');
       else
          w2=(berm_width-0.25*L0)/(0.75*L0)
          w1=1-w2
          R2_new=w2*fore_R2 + w1*R2_new
       end
    end % end berm width check

    % convergence criterion
    R2del=abs(R2-R2_new)
    R2_all(iter)=R2_new;

    % get the new top station (for plot purposes)
    Z2=R2_new+SWEL
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end);
    end
    topStaAll(iter)=top_sta;

end

% final 2% runup elevation
Z2=R2_new+SWEL

diary off


% make some plots


% get the full transect



file=xlsread('../ADCIRC_returns/YK-106XYZSTA_RETURNS.csv');
sta_full=file(:,4);
ele_full=file(:,3);


%determine SWL station
SWL_sta=-999;
for kk=1:length(sta)-1
    if ((SWEL > dep(kk)) & (SWEL <= dep(kk+1)))    % here is the intersection of z2 with profile
       SWL_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),SWEL)
       break;
    end
end

figure('position', [100 100 900 900])
hold on
plot(sta_full,ele_full,'linewidth',2,'color',[0.5 0.2 0.1]);
plot(sta,dep,'md');
grid on
plot([min(sta),SWL_sta],[SWEL SWEL],'color','b','linewidth',2)
%plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
scatter(topStaAll,R2_all+SWEL,30,'r*');
%scatter(toe_sta,Db,70,'k*');
scatter(toe_sta,Ztoe,70,'k*');
% plot berm segs
for nn=1:length(Berm_Segs)
   seg=Berm_Segs(nn);
   h=Berm_Heights(nn);
   plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
end
xlim([toe_sta-50 max(topStaAll)+50])
grid minor
hleg=legend('Ground','TAW Segments','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
set (hleg,'location','southeast');
xlabel('Distance from zero NAVD88 shoreline (feet)');
ylabel('Elevation (Feet-NAVD88)');
str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
yylm=get (gca,'ylim');
dztxt=(yylm(end)-yylm(1))/15;
text(toe_sta,Z2+dztxt/4,str);
str=sprintf('Run-up Height is %5.1f Feet',R2_new);
text(toe_sta,Z2-dztxt/4,str);
title (plotTitle);
set(gcf,'color','w');
print('-r450','-dpdf',imgname);

-1.000000e+00
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S9');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U9');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O9');
clear all
close all
format long g

diary logfiles/YK-107-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Kennebunkport, York county, Maine
% TRANSECT ID: YK-107
% calculation by SJH, Ransom Consulting, Inc. 02-Apr-2020
% 100-year wave runup using TAW methodology
% including berm and weighted average with foreshore if necessary
%
% chk nld 20200220
%
% This script assumes that the incident wave conditions provided
% as input in the configuration section below are the 
% appropriate values located at the end of the foreshore
% or toe of the slope on which the run-up is being calculated
% the script does not attempt to apply a depth limit or any other
% transformation to the incident wave conditions other than 
% conversion of the peak wave period to the spectral mean wave 
% as recommended in the references below
% 
% references:
%
% Van der Meer, J.W., 2002. Technical Report Wave Run-up and 
% Wave Overtopping at Dikes. TAW Technical Advisory Committee on
% Flood Defence, The Netherlands.
%
% FEMA. 2007,  Atlantic Ocean and Gulf of Mexico Coastal Guidelines Update
%
% 
%

%--------------------------------------------------------------------------
% CONFIG
%--------------------------------------------------------------------------
fname='inpfiles/YK-107sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-107-runup';                            
SWEL=9.3596;  % 100-yr still water level including wave setup. 
H0=5.7654;    % significant wave height at toe of structure 
Tp=13.2478;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.70525;   % this may get changed automatically below
gamma_rough=0.6;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.0031758; 
maxSetup=0.61291;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-107'


% END CONFIG
%-------------------------------

SWEL=SWEL+setupAtToe  
SWEL_fore=SWEL+maxSetup


% FIND WAVELENGTH USING DEEPWATER DISPERSION RELATION
% using English units
L0=32.15/(2*pi)*T0^2

% Find Hb (Munk, 1949)
%Hb=H0/(3.3*(H0/L0)^(1/3)) 
%Db=-Hb/.78+SWEL;  % depth at breaking

% The toe elevation here is only used to determine the average 
% structure slope, it is not used to depth limit the wave height.
% Any depth limiting or other modification of the wave height 
% to make it consitent with TAW guidance should be performed 
% prior to the input of the significant wave height given above. 
Ztoe=SWEL-1.5*H0

% read the transect
[sta,dep,inc] = textread(fname,'%n%n%n%*[^\n]','delimiter',',','headerlines',0);

% remove unselected points
k=find(inc==0);
sta(k)=[];
dep(k)=[];

sta_org=sta;  % used for plotting purposes
dep_org=dep;


% initial guess at maximum run-up elevation to estimate slope
Z2=SWEL+1.5*H0

% determine station at the max runup and -1.5*H0 (i.e. the toe)
top_sta=-999;
toe_sta=-999;
for kk=1:length(sta)-1
    if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
       top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
    end
    if  ((Ztoe > dep(kk)) & (Ztoe <= dep(kk+1)))    % here is the intersection of Ztoe with profile
       toe_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Ztoe)
    end
end
% check to make sure we got them, if not extend the end slopes outward
S=diff(dep)./diff(sta);
if toe_sta==-999
   dy=dep(1)-Ztoe;
   toe_sta=sta(1)-dy/S(1)
end
if top_sta==-999
   dy=Z2-dep(end);
   top_sta=sta(end)+dy/S(end)
end

% just so the reader can tell the values aren't -999 anymore
top_sta
toe_sta



% check for case where the toe of slope is below SWL-1.5*H0
% in this case interpolate setup from the setupAtToe(really setup as first station), and the max setup
% also un-include points seaward of SWL-1.5*H0
if Ztoe > dep(1)
   dd=SWEL_fore-dep;
   k=find(dd<0,1); % k is index of first land point
   staAtSWL=interp1(dep(k-1:k),sta(k-1:k),SWEL_fore);   
   dsta=staAtSWL-sta(1);
   dsetup=maxSetup-setupAtToe;
   dsetdsta=dsetup/dsta;
   setup=setupAtToe+dsetdsta*(toe_sta-sta(1));
   sprintf('-!!- Location of SWEL-1.5*H0 is %4.1f ft landward of toe of slope',dsta)
   sprintf('-!!- Setup is interpolated between setup at toe of slope and max setup')
   sprintf('-!!-       setup is adjusted to %4.2f feet',setup)
   SWEL=SWEL-setupAtToe+setup;
   sprintf('-!!-       SWEL is adjusted to %4.2f feet',SWEL)
   k=find(dep < SWEL-1.5*H0)
   sta(k)=[];
   dep(k)=[];
else
   sprintf('-!!- The User has selected a starting point that is %4.2f feet above the elevation of SWEL-1.5H0\n',dep(1)-Ztoe)
   sprintf('-!!- This may be reasonable for some cases.  However the user may want to consider:\n')
   sprintf('-!!-   1) Selecting a starting point that is at or below %4.2f feet elevation, or\n', Ztoe)
   sprintf('-!!-   2) Reducing the incident wave height to a depth limited condition.\n')
end



% now iterate converge on a runup elevation
tol=0.01;  % convergence criteria
R2del=999;
R2_new=3*H0; %initial guess
R2=R2_new;
iter=0;
R2_all=[];
topStaAll=[];
Berm_Segs=[];
TAW_ALWAYS_VALID=1;
while(abs(R2del) > tol && iter <= 25)
    iter=iter+1;
    sprintf ('!----------------- STARTING ITERATION %d -----------------!',iter)
    % elevation of toe of slope
    Ztoe
    % station of toe slope (relative to 0-NAVD88 shoreline
    toe_sta 
    % station of top of slope/extent of 2% run-up
    top_sta
    % elevation of top of slope/extent of 2% run-up
    Z2
    % incident significant wave height
    H0
    % incident spectral peak wave period
    Tp
    % incident spectral mean wave period
    T0
        
    R2=R2_new  
    Z2=R2+SWEL
    % determine slope for this iteration
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end)
    end
    
    % get the length of the slope (not accounting for berm)
    Lslope=top_sta-toe_sta

    % loop over profile segments to determine berm factor 
    % re-calculate influence of depth of berm based on this run-up elevation
    % check for berm, berm width, berm height
    berm_width=0;
    rdh_sum=0;
    Berm_Segs=[];
    Berm_Heights=[];
    for kk=1:length(sta)-1
       ddep=dep(kk+1)-dep(kk);
       dsta=sta(kk+1)-sta(kk);
       s=ddep/dsta;
       if  (s < 1/15)     % count it as a berm if slope is flatter than 1:15 (see TAW manual)
          sprintf ('Berm Factor Calculation: Iteration %d, Profile Segment: %d',iter,kk) 
          berm_width=berm_width+dsta;   % tally the width of all berm segments
          % compute the rdh for this segment and weight it by the segment length
          dh=SWEL-(dep(kk)+dep(kk+1))/2  
          if  dh < 0
              chi=R2;
          else 
              chi=2* H0;
          end
          if (dh <= R2 & dh >=-2*H0)
             rdh=(0.5-0.5*cos(3.14159*dh/chi)) ;
          else
             rdh=1;
          end
          rdh_sum=rdh_sum + rdh * dsta
          Berm_Segs=[Berm_Segs, kk];
          Berm_Heights=[Berm_Heights, (dep(kk)+dep(kk+1))/2];
       end
       if dep(kk) >= Z2  % jump out of loop if we reached limit of run-up for this iteration 
          break
       end
    end
    sprintf ('!------- End Berm Factor Calculation, Iter: %d ---------!',iter)
    berm_width
    rB=berm_width/Lslope 
    if (berm_width > 0)
       rdh_mean=rdh_sum/berm_width
    else
       rdh_mean=1
    end
    gamma_berm=1- rB * (1-rdh_mean)
    if gamma_berm > 1
       gamma_berm=1
    end
    if gamma_berm < 0.6
       gamma_berm =0.6
    end
    % Iribarren number
    slope=(Z2-Ztoe)/(Lslope-berm_width)
    Irb=(slope/(sqrt(H0/L0)))
    % runup height
    gamma_berm
    gamma_perm
    gamma_beta
    gamma_rough
    gamma=gamma_berm*gamma_perm*gamma_beta*gamma_rough

    % check validity
    TAW_VALID=1;
    if (Irb*gamma_berm < 0.5 | Irb*gamma_berm > 10 )
       sprintf('!!! - - Iribaren number: %6.2f is outside the valid range (0.5-10), TAW NOT VALID - - !!!\n', Irb*gamma_berm)
       TAW_VALID=0;
    else
       sprintf('!!! - - Iribaren number: %6.2f is in the valid range (0.5-10), TAW RECOMMENDED - - !!!\n', Irb*gamma_berm)
    end
    islope=1/slope;  
    if (slope < 1/8 | slope > 1)
       sprintf('!!! - - slope: 1:%3.1f V:H is outside the valid range (1:8 - 1:1), TAW NOT VALID - - !!!\n', islope)
       TAW_VALID=0;
    else
       sprintf('!!! - - slope: 1:%3.1f V:H is in the valid range (1:8 - 1:1), TAW RECOMMENDED - - !!!\n', islope)
    end
    if TAW_VALID == 0
       TAW_ALWAYS_VALID=0;
    end
         
    if (Irb*gamma_berm < 1.8)
       R2_new=gamma*H0*1.77*Irb
    else
       R2_new=gamma*H0*(4.3-(1.6/sqrt(Irb))) 
    end


    % check to see if we need to evaluate a shallow foreshore
    if berm_width > 0.25 * L0;
       disp ('!   Berm_width is greater than 1/4 wave length')
       disp ('!   Runup will be weighted average with foreshore calculation assuming depth limited wave height on berm');
       % do the foreshore calculation 
       fore_H0=0.78*(SWEL_fore-min(Berm_Heights))
       % get upper slope
       fore_toe_sta=-999;
       fore_toe_dep=-999;
       for kk=length(dep)-1:-1:1
          ddep=dep(kk+1)-dep(kk);
          dsta=sta(kk+1)-sta(kk);
          s=ddep/dsta;
          if s < 1/15
             break
          end
          fore_toe_sta=sta(kk);
          fore_toe_dep=dep(kk);
          upper_slope=(Z2-fore_toe_dep)/(top_sta-fore_toe_sta)
       end
       fore_Irb=upper_slope/(sqrt(fore_H0/L0));
       fore_gamma=gamma_perm*gamma_beta*gamma_rough;
       if (fore_Irb < 1.8)
          fore_R2=fore_gamma*fore_H0*1.77*fore_Irb;
       else
          fore_R2=fore_gamma*fore_H0*(4.3-(1.6/sqrt(fore_Irb))); 
       end
       if berm_width >= L0
          R2_new=fore_R2
          disp ('berm is wider than one wavelength, use full shallow foreshore solution');
       else
          w2=(berm_width-0.25*L0)/(0.75*L0)
          w1=1-w2
          R2_new=w2*fore_R2 + w1*R2_new
       end
    end % end berm width check

    % convergence criterion
    R2del=abs(R2-R2_new)
    R2_all(iter)=R2_new;

    % get the new top station (for plot purposes)
    Z2=R2_new+SWEL
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end);
    end
    topStaAll(iter)=top_sta;

end

% final 2% runup elevation
Z2=R2_new+SWEL

diary off


% make some plots


% get the full transect



file=xlsread('../ADCIRC_returns/YK-107XYZSTA_RETURNS.csv');
sta_full=file(:,4);
ele_full=file(:,3);


%determine SWL station
SWL_sta=-999;
for kk=1:length(sta)-1
    if ((SWEL > dep(kk)) & (SWEL <= dep(kk+1)))    % here is the intersection of z2 with profile
       SWL_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),SWEL)
       break;
    end
end

figure('position', [100 100 900 900])
hold on
plot(sta_full,ele_full,'linewidth',2,'color',[0.5 0.2 0.1]);
plot(sta,dep,'md');
grid on
plot([min(sta),SWL_sta],[SWEL SWEL],'color','b','linewidth',2)
%plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
scatter(topStaAll,R2_all+SWEL,30,'r*');
%scatter(toe_sta,Db,70,'k*');
scatter(toe_sta,Ztoe,70,'k*');
% plot berm segs
for nn=1:length(Berm_Segs)
   seg=Berm_Segs(nn);
   h=Berm_Heights(nn);
   plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
end
xlim([toe_sta-50 max(topStaAll)+50])
grid minor
hleg=legend('Ground','TAW Segments','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
set (hleg,'location','southeast');
xlabel('Distance from zero NAVD88 shoreline (feet)');
ylabel('Elevation (Feet-NAVD88)');
str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
yylm=get (gca,'ylim');
dztxt=(yylm(end)-yylm(1))/15;
text(toe_sta,Z2+dztxt/4,str);
str=sprintf('Run-up Height is %5.1f Feet',R2_new);
text(toe_sta,Z2-dztxt/4,str);
title (plotTitle);
set(gcf,'color','w');
print('-r450','-dpdf',imgname);

-1.000000e+00
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S10');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U10');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O10');
clear all
close all
format long g

diary logfiles/YK-108-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Kennebunkport, York county, Maine
% TRANSECT ID: YK-108
% calculation by SJH, Ransom Consulting, Inc. 02-Apr-2020
% 100-year wave runup using TAW methodology
% including berm and weighted average with foreshore if necessary
%
% chk nld 20200220
%
% This script assumes that the incident wave conditions provided
% as input in the configuration section below are the 
% appropriate values located at the end of the foreshore
% or toe of the slope on which the run-up is being calculated
% the script does not attempt to apply a depth limit or any other
% transformation to the incident wave conditions other than 
% conversion of the peak wave period to the spectral mean wave 
% as recommended in the references below
% 
% references:
%
% Van der Meer, J.W., 2002. Technical Report Wave Run-up and 
% Wave Overtopping at Dikes. TAW Technical Advisory Committee on
% Flood Defence, The Netherlands.
%
% FEMA. 2007,  Atlantic Ocean and Gulf of Mexico Coastal Guidelines Update
%
% 
%

%--------------------------------------------------------------------------
% CONFIG
%--------------------------------------------------------------------------
fname='inpfiles/YK-108sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-108-runup';                            
SWEL=9.0347;  % 100-yr still water level including wave setup. 
H0=2.4638;    % significant wave height at toe of structure 
Tp=12.522;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.99585;   % this may get changed automatically below
gamma_rough=0.6;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.41366; 
maxSetup=0.68007;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-108'


% END CONFIG
%-------------------------------

SWEL=SWEL+setupAtToe  
SWEL_fore=SWEL+maxSetup


% FIND WAVELENGTH USING DEEPWATER DISPERSION RELATION
% using English units
L0=32.15/(2*pi)*T0^2

% Find Hb (Munk, 1949)
%Hb=H0/(3.3*(H0/L0)^(1/3)) 
%Db=-Hb/.78+SWEL;  % depth at breaking

% The toe elevation here is only used to determine the average 
% structure slope, it is not used to depth limit the wave height.
% Any depth limiting or other modification of the wave height 
% to make it consitent with TAW guidance should be performed 
% prior to the input of the significant wave height given above. 
Ztoe=SWEL-1.5*H0

% read the transect
[sta,dep,inc] = textread(fname,'%n%n%n%*[^\n]','delimiter',',','headerlines',0);

% remove unselected points
k=find(inc==0);
sta(k)=[];
dep(k)=[];

sta_org=sta;  % used for plotting purposes
dep_org=dep;


% initial guess at maximum run-up elevation to estimate slope
Z2=SWEL+1.5*H0

% determine station at the max runup and -1.5*H0 (i.e. the toe)
top_sta=-999;
toe_sta=-999;
for kk=1:length(sta)-1
    if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
       top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
    end
    if  ((Ztoe > dep(kk)) & (Ztoe <= dep(kk+1)))    % here is the intersection of Ztoe with profile
       toe_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Ztoe)
    end
end
% check to make sure we got them, if not extend the end slopes outward
S=diff(dep)./diff(sta);
if toe_sta==-999
   dy=dep(1)-Ztoe;
   toe_sta=sta(1)-dy/S(1)
end
if top_sta==-999
   dy=Z2-dep(end);
   top_sta=sta(end)+dy/S(end)
end

% just so the reader can tell the values aren't -999 anymore
top_sta
toe_sta



% check for case where the toe of slope is below SWL-1.5*H0
% in this case interpolate setup from the setupAtToe(really setup as first station), and the max setup
% also un-include points seaward of SWL-1.5*H0
if Ztoe > dep(1)
   dd=SWEL_fore-dep;
   k=find(dd<0,1); % k is index of first land point
   staAtSWL=interp1(dep(k-1:k),sta(k-1:k),SWEL_fore);   
   dsta=staAtSWL-sta(1);
   dsetup=maxSetup-setupAtToe;
   dsetdsta=dsetup/dsta;
   setup=setupAtToe+dsetdsta*(toe_sta-sta(1));
   sprintf('-!!- Location of SWEL-1.5*H0 is %4.1f ft landward of toe of slope',dsta)
   sprintf('-!!- Setup is interpolated between setup at toe of slope and max setup')
   sprintf('-!!-       setup is adjusted to %4.2f feet',setup)
   SWEL=SWEL-setupAtToe+setup;
   sprintf('-!!-       SWEL is adjusted to %4.2f feet',SWEL)
   k=find(dep < SWEL-1.5*H0)
   sta(k)=[];
   dep(k)=[];
else
   sprintf('-!!- The User has selected a starting point that is %4.2f feet above the elevation of SWEL-1.5H0\n',dep(1)-Ztoe)
   sprintf('-!!- This may be reasonable for some cases.  However the user may want to consider:\n')
   sprintf('-!!-   1) Selecting a starting point that is at or below %4.2f feet elevation, or\n', Ztoe)
   sprintf('-!!-   2) Reducing the incident wave height to a depth limited condition.\n')
end



% now iterate converge on a runup elevation
tol=0.01;  % convergence criteria
R2del=999;
R2_new=3*H0; %initial guess
R2=R2_new;
iter=0;
R2_all=[];
topStaAll=[];
Berm_Segs=[];
TAW_ALWAYS_VALID=1;
while(abs(R2del) > tol && iter <= 25)
    iter=iter+1;
    sprintf ('!----------------- STARTING ITERATION %d -----------------!',iter)
    % elevation of toe of slope
    Ztoe
    % station of toe slope (relative to 0-NAVD88 shoreline
    toe_sta 
    % station of top of slope/extent of 2% run-up
    top_sta
    % elevation of top of slope/extent of 2% run-up
    Z2
    % incident significant wave height
    H0
    % incident spectral peak wave period
    Tp
    % incident spectral mean wave period
    T0
        
    R2=R2_new  
    Z2=R2+SWEL
    % determine slope for this iteration
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end)
    end
    
    % get the length of the slope (not accounting for berm)
    Lslope=top_sta-toe_sta

    % loop over profile segments to determine berm factor 
    % re-calculate influence of depth of berm based on this run-up elevation
    % check for berm, berm width, berm height
    berm_width=0;
    rdh_sum=0;
    Berm_Segs=[];
    Berm_Heights=[];
    for kk=1:length(sta)-1
       ddep=dep(kk+1)-dep(kk);
       dsta=sta(kk+1)-sta(kk);
       s=ddep/dsta;
       if  (s < 1/15)     % count it as a berm if slope is flatter than 1:15 (see TAW manual)
          sprintf ('Berm Factor Calculation: Iteration %d, Profile Segment: %d',iter,kk) 
          berm_width=berm_width+dsta;   % tally the width of all berm segments
          % compute the rdh for this segment and weight it by the segment length
          dh=SWEL-(dep(kk)+dep(kk+1))/2  
          if  dh < 0
              chi=R2;
          else 
              chi=2* H0;
          end
          if (dh <= R2 & dh >=-2*H0)
             rdh=(0.5-0.5*cos(3.14159*dh/chi)) ;
          else
             rdh=1;
          end
          rdh_sum=rdh_sum + rdh * dsta
          Berm_Segs=[Berm_Segs, kk];
          Berm_Heights=[Berm_Heights, (dep(kk)+dep(kk+1))/2];
       end
       if dep(kk) >= Z2  % jump out of loop if we reached limit of run-up for this iteration 
          break
       end
    end
    sprintf ('!------- End Berm Factor Calculation, Iter: %d ---------!',iter)
    berm_width
    rB=berm_width/Lslope 
    if (berm_width > 0)
       rdh_mean=rdh_sum/berm_width
    else
       rdh_mean=1
    end
    gamma_berm=1- rB * (1-rdh_mean)
    if gamma_berm > 1
       gamma_berm=1
    end
    if gamma_berm < 0.6
       gamma_berm =0.6
    end
    % Iribarren number
    slope=(Z2-Ztoe)/(Lslope-berm_width)
    Irb=(slope/(sqrt(H0/L0)))
    % runup height
    gamma_berm
    gamma_perm
    gamma_beta
    gamma_rough
    gamma=gamma_berm*gamma_perm*gamma_beta*gamma_rough

    % check validity
    TAW_VALID=1;
    if (Irb*gamma_berm < 0.5 | Irb*gamma_berm > 10 )
       sprintf('!!! - - Iribaren number: %6.2f is outside the valid range (0.5-10), TAW NOT VALID - - !!!\n', Irb*gamma_berm)
       TAW_VALID=0;
    else
       sprintf('!!! - - Iribaren number: %6.2f is in the valid range (0.5-10), TAW RECOMMENDED - - !!!\n', Irb*gamma_berm)
    end
    islope=1/slope;  
    if (slope < 1/8 | slope > 1)
       sprintf('!!! - - slope: 1:%3.1f V:H is outside the valid range (1:8 - 1:1), TAW NOT VALID - - !!!\n', islope)
       TAW_VALID=0;
    else
       sprintf('!!! - - slope: 1:%3.1f V:H is in the valid range (1:8 - 1:1), TAW RECOMMENDED - - !!!\n', islope)
    end
    if TAW_VALID == 0
       TAW_ALWAYS_VALID=0;
    end
         
    if (Irb*gamma_berm < 1.8)
       R2_new=gamma*H0*1.77*Irb
    else
       R2_new=gamma*H0*(4.3-(1.6/sqrt(Irb))) 
    end


    % check to see if we need to evaluate a shallow foreshore
    if berm_width > 0.25 * L0;
       disp ('!   Berm_width is greater than 1/4 wave length')
       disp ('!   Runup will be weighted average with foreshore calculation assuming depth limited wave height on berm');
       % do the foreshore calculation 
       fore_H0=0.78*(SWEL_fore-min(Berm_Heights))
       % get upper slope
       fore_toe_sta=-999;
       fore_toe_dep=-999;
       for kk=length(dep)-1:-1:1
          ddep=dep(kk+1)-dep(kk);
          dsta=sta(kk+1)-sta(kk);
          s=ddep/dsta;
          if s < 1/15
             break
          end
          fore_toe_sta=sta(kk);
          fore_toe_dep=dep(kk);
          upper_slope=(Z2-fore_toe_dep)/(top_sta-fore_toe_sta)
       end
       fore_Irb=upper_slope/(sqrt(fore_H0/L0));
       fore_gamma=gamma_perm*gamma_beta*gamma_rough;
       if (fore_Irb < 1.8)
          fore_R2=fore_gamma*fore_H0*1.77*fore_Irb;
       else
          fore_R2=fore_gamma*fore_H0*(4.3-(1.6/sqrt(fore_Irb))); 
       end
       if berm_width >= L0
          R2_new=fore_R2
          disp ('berm is wider than one wavelength, use full shallow foreshore solution');
       else
          w2=(berm_width-0.25*L0)/(0.75*L0)
          w1=1-w2
          R2_new=w2*fore_R2 + w1*R2_new
       end
    end % end berm width check

    % convergence criterion
    R2del=abs(R2-R2_new)
    R2_all(iter)=R2_new;

    % get the new top station (for plot purposes)
    Z2=R2_new+SWEL
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end);
    end
    topStaAll(iter)=top_sta;

end

% final 2% runup elevation
Z2=R2_new+SWEL

diary off


% make some plots


% get the full transect



file=xlsread('../ADCIRC_returns/YK-108XYZSTA_RETURNS.csv');
sta_full=file(:,4);
ele_full=file(:,3);


%determine SWL station
SWL_sta=-999;
for kk=1:length(sta)-1
    if ((SWEL > dep(kk)) & (SWEL <= dep(kk+1)))    % here is the intersection of z2 with profile
       SWL_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),SWEL)
       break;
    end
end

figure('position', [100 100 900 900])
hold on
plot(sta_full,ele_full,'linewidth',2,'color',[0.5 0.2 0.1]);
plot(sta,dep,'md');
grid on
plot([min(sta),SWL_sta],[SWEL SWEL],'color','b','linewidth',2)
%plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
scatter(topStaAll,R2_all+SWEL,30,'r*');
%scatter(toe_sta,Db,70,'k*');
scatter(toe_sta,Ztoe,70,'k*');
% plot berm segs
for nn=1:length(Berm_Segs)
   seg=Berm_Segs(nn);
   h=Berm_Heights(nn);
   plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
end
xlim([toe_sta-50 max(topStaAll)+50])
grid minor
hleg=legend('Ground','TAW Segments','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
set (hleg,'location','southeast');
xlabel('Distance from zero NAVD88 shoreline (feet)');
ylabel('Elevation (Feet-NAVD88)');
str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
yylm=get (gca,'ylim');
dztxt=(yylm(end)-yylm(1))/15;
text(toe_sta,Z2+dztxt/4,str);
str=sprintf('Run-up Height is %5.1f Feet',R2_new);
text(toe_sta,Z2-dztxt/4,str);
title (plotTitle);
set(gcf,'color','w');
print('-r450','-dpdf',imgname);

-1.000000e+00
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S11');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U11');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O11');
clear all
close all
format long g

diary logfiles/YK-109-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Kennebunkport, York county, Maine
% TRANSECT ID: YK-109
% calculation by SJH, Ransom Consulting, Inc. 02-Apr-2020
% 100-year wave runup using TAW methodology
% including berm and weighted average with foreshore if necessary
%
% chk nld 20200220
%
% This script assumes that the incident wave conditions provided
% as input in the configuration section below are the 
% appropriate values located at the end of the foreshore
% or toe of the slope on which the run-up is being calculated
% the script does not attempt to apply a depth limit or any other
% transformation to the incident wave conditions other than 
% conversion of the peak wave period to the spectral mean wave 
% as recommended in the references below
% 
% references:
%
% Van der Meer, J.W., 2002. Technical Report Wave Run-up and 
% Wave Overtopping at Dikes. TAW Technical Advisory Committee on
% Flood Defence, The Netherlands.
%
% FEMA. 2007,  Atlantic Ocean and Gulf of Mexico Coastal Guidelines Update
%
% 
%

%--------------------------------------------------------------------------
% CONFIG
%--------------------------------------------------------------------------
fname='inpfiles/YK-109sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-109-runup';                            
SWEL=9.0102;  % 100-yr still water level including wave setup. 
H0=5.2103;    % significant wave height at toe of structure 
Tp=13.7882;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.61088;   % this may get changed automatically below
gamma_rough=0.6;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.19397; 
maxSetup=0.8081;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-109'


% END CONFIG
%-------------------------------

SWEL=SWEL+setupAtToe  
SWEL_fore=SWEL+maxSetup


% FIND WAVELENGTH USING DEEPWATER DISPERSION RELATION
% using English units
L0=32.15/(2*pi)*T0^2

% Find Hb (Munk, 1949)
%Hb=H0/(3.3*(H0/L0)^(1/3)) 
%Db=-Hb/.78+SWEL;  % depth at breaking

% The toe elevation here is only used to determine the average 
% structure slope, it is not used to depth limit the wave height.
% Any depth limiting or other modification of the wave height 
% to make it consitent with TAW guidance should be performed 
% prior to the input of the significant wave height given above. 
Ztoe=SWEL-1.5*H0

% read the transect
[sta,dep,inc] = textread(fname,'%n%n%n%*[^\n]','delimiter',',','headerlines',0);

% remove unselected points
k=find(inc==0);
sta(k)=[];
dep(k)=[];

sta_org=sta;  % used for plotting purposes
dep_org=dep;


% initial guess at maximum run-up elevation to estimate slope
Z2=SWEL+1.5*H0

% determine station at the max runup and -1.5*H0 (i.e. the toe)
top_sta=-999;
toe_sta=-999;
for kk=1:length(sta)-1
    if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
       top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
    end
    if  ((Ztoe > dep(kk)) & (Ztoe <= dep(kk+1)))    % here is the intersection of Ztoe with profile
       toe_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Ztoe)
    end
end
% check to make sure we got them, if not extend the end slopes outward
S=diff(dep)./diff(sta);
if toe_sta==-999
   dy=dep(1)-Ztoe;
   toe_sta=sta(1)-dy/S(1)
end
if top_sta==-999
   dy=Z2-dep(end);
   top_sta=sta(end)+dy/S(end)
end

% just so the reader can tell the values aren't -999 anymore
top_sta
toe_sta



% check for case where the toe of slope is below SWL-1.5*H0
% in this case interpolate setup from the setupAtToe(really setup as first station), and the max setup
% also un-include points seaward of SWL-1.5*H0
if Ztoe > dep(1)
   dd=SWEL_fore-dep;
   k=find(dd<0,1); % k is index of first land point
   staAtSWL=interp1(dep(k-1:k),sta(k-1:k),SWEL_fore);   
   dsta=staAtSWL-sta(1);
   dsetup=maxSetup-setupAtToe;
   dsetdsta=dsetup/dsta;
   setup=setupAtToe+dsetdsta*(toe_sta-sta(1));
   sprintf('-!!- Location of SWEL-1.5*H0 is %4.1f ft landward of toe of slope',dsta)
   sprintf('-!!- Setup is interpolated between setup at toe of slope and max setup')
   sprintf('-!!-       setup is adjusted to %4.2f feet',setup)
   SWEL=SWEL-setupAtToe+setup;
   sprintf('-!!-       SWEL is adjusted to %4.2f feet',SWEL)
   k=find(dep < SWEL-1.5*H0)
   sta(k)=[];
   dep(k)=[];
else
   sprintf('-!!- The User has selected a starting point that is %4.2f feet above the elevation of SWEL-1.5H0\n',dep(1)-Ztoe)
   sprintf('-!!- This may be reasonable for some cases.  However the user may want to consider:\n')
   sprintf('-!!-   1) Selecting a starting point that is at or below %4.2f feet elevation, or\n', Ztoe)
   sprintf('-!!-   2) Reducing the incident wave height to a depth limited condition.\n')
end



% now iterate converge on a runup elevation
tol=0.01;  % convergence criteria
R2del=999;
R2_new=3*H0; %initial guess
R2=R2_new;
iter=0;
R2_all=[];
topStaAll=[];
Berm_Segs=[];
TAW_ALWAYS_VALID=1;
while(abs(R2del) > tol && iter <= 25)
    iter=iter+1;
    sprintf ('!----------------- STARTING ITERATION %d -----------------!',iter)
    % elevation of toe of slope
    Ztoe
    % station of toe slope (relative to 0-NAVD88 shoreline
    toe_sta 
    % station of top of slope/extent of 2% run-up
    top_sta
    % elevation of top of slope/extent of 2% run-up
    Z2
    % incident significant wave height
    H0
    % incident spectral peak wave period
    Tp
    % incident spectral mean wave period
    T0
        
    R2=R2_new  
    Z2=R2+SWEL
    % determine slope for this iteration
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end)
    end
    
    % get the length of the slope (not accounting for berm)
    Lslope=top_sta-toe_sta

    % loop over profile segments to determine berm factor 
    % re-calculate influence of depth of berm based on this run-up elevation
    % check for berm, berm width, berm height
    berm_width=0;
    rdh_sum=0;
    Berm_Segs=[];
    Berm_Heights=[];
    for kk=1:length(sta)-1
       ddep=dep(kk+1)-dep(kk);
       dsta=sta(kk+1)-sta(kk);
       s=ddep/dsta;
       if  (s < 1/15)     % count it as a berm if slope is flatter than 1:15 (see TAW manual)
          sprintf ('Berm Factor Calculation: Iteration %d, Profile Segment: %d',iter,kk) 
          berm_width=berm_width+dsta;   % tally the width of all berm segments
          % compute the rdh for this segment and weight it by the segment length
          dh=SWEL-(dep(kk)+dep(kk+1))/2  
          if  dh < 0
              chi=R2;
          else 
              chi=2* H0;
          end
          if (dh <= R2 & dh >=-2*H0)
             rdh=(0.5-0.5*cos(3.14159*dh/chi)) ;
          else
             rdh=1;
          end
          rdh_sum=rdh_sum + rdh * dsta
          Berm_Segs=[Berm_Segs, kk];
          Berm_Heights=[Berm_Heights, (dep(kk)+dep(kk+1))/2];
       end
       if dep(kk) >= Z2  % jump out of loop if we reached limit of run-up for this iteration 
          break
       end
    end
    sprintf ('!------- End Berm Factor Calculation, Iter: %d ---------!',iter)
    berm_width
    rB=berm_width/Lslope 
    if (berm_width > 0)
       rdh_mean=rdh_sum/berm_width
    else
       rdh_mean=1
    end
    gamma_berm=1- rB * (1-rdh_mean)
    if gamma_berm > 1
       gamma_berm=1
    end
    if gamma_berm < 0.6
       gamma_berm =0.6
    end
    % Iribarren number
    slope=(Z2-Ztoe)/(Lslope-berm_width)
    Irb=(slope/(sqrt(H0/L0)))
    % runup height
    gamma_berm
    gamma_perm
    gamma_beta
    gamma_rough
    gamma=gamma_berm*gamma_perm*gamma_beta*gamma_rough

    % check validity
    TAW_VALID=1;
    if (Irb*gamma_berm < 0.5 | Irb*gamma_berm > 10 )
       sprintf('!!! - - Iribaren number: %6.2f is outside the valid range (0.5-10), TAW NOT VALID - - !!!\n', Irb*gamma_berm)
       TAW_VALID=0;
    else
       sprintf('!!! - - Iribaren number: %6.2f is in the valid range (0.5-10), TAW RECOMMENDED - - !!!\n', Irb*gamma_berm)
    end
    islope=1/slope;  
    if (slope < 1/8 | slope > 1)
       sprintf('!!! - - slope: 1:%3.1f V:H is outside the valid range (1:8 - 1:1), TAW NOT VALID - - !!!\n', islope)
       TAW_VALID=0;
    else
       sprintf('!!! - - slope: 1:%3.1f V:H is in the valid range (1:8 - 1:1), TAW RECOMMENDED - - !!!\n', islope)
    end
    if TAW_VALID == 0
       TAW_ALWAYS_VALID=0;
    end
         
    if (Irb*gamma_berm < 1.8)
       R2_new=gamma*H0*1.77*Irb
    else
       R2_new=gamma*H0*(4.3-(1.6/sqrt(Irb))) 
    end


    % check to see if we need to evaluate a shallow foreshore
    if berm_width > 0.25 * L0;
       disp ('!   Berm_width is greater than 1/4 wave length')
       disp ('!   Runup will be weighted average with foreshore calculation assuming depth limited wave height on berm');
       % do the foreshore calculation 
       fore_H0=0.78*(SWEL_fore-min(Berm_Heights))
       % get upper slope
       fore_toe_sta=-999;
       fore_toe_dep=-999;
       for kk=length(dep)-1:-1:1
          ddep=dep(kk+1)-dep(kk);
          dsta=sta(kk+1)-sta(kk);
          s=ddep/dsta;
          if s < 1/15
             break
          end
          fore_toe_sta=sta(kk);
          fore_toe_dep=dep(kk);
          upper_slope=(Z2-fore_toe_dep)/(top_sta-fore_toe_sta)
       end
       fore_Irb=upper_slope/(sqrt(fore_H0/L0));
       fore_gamma=gamma_perm*gamma_beta*gamma_rough;
       if (fore_Irb < 1.8)
          fore_R2=fore_gamma*fore_H0*1.77*fore_Irb;
       else
          fore_R2=fore_gamma*fore_H0*(4.3-(1.6/sqrt(fore_Irb))); 
       end
       if berm_width >= L0
          R2_new=fore_R2
          disp ('berm is wider than one wavelength, use full shallow foreshore solution');
       else
          w2=(berm_width-0.25*L0)/(0.75*L0)
          w1=1-w2
          R2_new=w2*fore_R2 + w1*R2_new
       end
    end % end berm width check

    % convergence criterion
    R2del=abs(R2-R2_new)
    R2_all(iter)=R2_new;

    % get the new top station (for plot purposes)
    Z2=R2_new+SWEL
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end);
    end
    topStaAll(iter)=top_sta;

end

% final 2% runup elevation
Z2=R2_new+SWEL

diary off


% make some plots


% get the full transect



file=xlsread('../ADCIRC_returns/YK-109XYZSTA_RETURNS.csv');
sta_full=file(:,4);
ele_full=file(:,3);


%determine SWL station
SWL_sta=-999;
for kk=1:length(sta)-1
    if ((SWEL > dep(kk)) & (SWEL <= dep(kk+1)))    % here is the intersection of z2 with profile
       SWL_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),SWEL)
       break;
    end
end

figure('position', [100 100 900 900])
hold on
plot(sta_full,ele_full,'linewidth',2,'color',[0.5 0.2 0.1]);
plot(sta,dep,'md');
grid on
plot([min(sta),SWL_sta],[SWEL SWEL],'color','b','linewidth',2)
%plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
scatter(topStaAll,R2_all+SWEL,30,'r*');
%scatter(toe_sta,Db,70,'k*');
scatter(toe_sta,Ztoe,70,'k*');
% plot berm segs
for nn=1:length(Berm_Segs)
   seg=Berm_Segs(nn);
   h=Berm_Heights(nn);
   plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
end
xlim([toe_sta-50 max(topStaAll)+50])
grid minor
hleg=legend('Ground','TAW Segments','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
set (hleg,'location','southeast');
xlabel('Distance from zero NAVD88 shoreline (feet)');
ylabel('Elevation (Feet-NAVD88)');
str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
yylm=get (gca,'ylim');
dztxt=(yylm(end)-yylm(1))/15;
text(toe_sta,Z2+dztxt/4,str);
str=sprintf('Run-up Height is %5.1f Feet',R2_new);
text(toe_sta,Z2-dztxt/4,str);
title (plotTitle);
set(gcf,'color','w');
print('-r450','-dpdf',imgname);

-1.000000e+00
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S12');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U12');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O12');
clear all
close all
format long g

diary logfiles/YK-110-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Kennebunkport, York county, Maine
% TRANSECT ID: YK-110
% calculation by SJH, Ransom Consulting, Inc. 02-Apr-2020
% 100-year wave runup using TAW methodology
% including berm and weighted average with foreshore if necessary
%
% chk nld 20200220
%
% This script assumes that the incident wave conditions provided
% as input in the configuration section below are the 
% appropriate values located at the end of the foreshore
% or toe of the slope on which the run-up is being calculated
% the script does not attempt to apply a depth limit or any other
% transformation to the incident wave conditions other than 
% conversion of the peak wave period to the spectral mean wave 
% as recommended in the references below
% 
% references:
%
% Van der Meer, J.W., 2002. Technical Report Wave Run-up and 
% Wave Overtopping at Dikes. TAW Technical Advisory Committee on
% Flood Defence, The Netherlands.
%
% FEMA. 2007,  Atlantic Ocean and Gulf of Mexico Coastal Guidelines Update
%
% 
%

%--------------------------------------------------------------------------
% CONFIG
%--------------------------------------------------------------------------
fname='inpfiles/YK-110sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-110-runup';                            
SWEL=9.0222;  % 100-yr still water level including wave setup. 
H0=3.8292;    % significant wave height at toe of structure 
Tp=15.3192;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.98222;   % this may get changed automatically below
gamma_rough=0.6;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.28621; 
maxSetup=0.69584;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-110'


% END CONFIG
%-------------------------------

SWEL=SWEL+setupAtToe  
SWEL_fore=SWEL+maxSetup


% FIND WAVELENGTH USING DEEPWATER DISPERSION RELATION
% using English units
L0=32.15/(2*pi)*T0^2

% Find Hb (Munk, 1949)
%Hb=H0/(3.3*(H0/L0)^(1/3)) 
%Db=-Hb/.78+SWEL;  % depth at breaking

% The toe elevation here is only used to determine the average 
% structure slope, it is not used to depth limit the wave height.
% Any depth limiting or other modification of the wave height 
% to make it consitent with TAW guidance should be performed 
% prior to the input of the significant wave height given above. 
Ztoe=SWEL-1.5*H0

% read the transect
[sta,dep,inc] = textread(fname,'%n%n%n%*[^\n]','delimiter',',','headerlines',0);

% remove unselected points
k=find(inc==0);
sta(k)=[];
dep(k)=[];

sta_org=sta;  % used for plotting purposes
dep_org=dep;


% initial guess at maximum run-up elevation to estimate slope
Z2=SWEL+1.5*H0

% determine station at the max runup and -1.5*H0 (i.e. the toe)
top_sta=-999;
toe_sta=-999;
for kk=1:length(sta)-1
    if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
       top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
    end
    if  ((Ztoe > dep(kk)) & (Ztoe <= dep(kk+1)))    % here is the intersection of Ztoe with profile
       toe_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Ztoe)
    end
end
% check to make sure we got them, if not extend the end slopes outward
S=diff(dep)./diff(sta);
if toe_sta==-999
   dy=dep(1)-Ztoe;
   toe_sta=sta(1)-dy/S(1)
end
if top_sta==-999
   dy=Z2-dep(end);
   top_sta=sta(end)+dy/S(end)
end

% just so the reader can tell the values aren't -999 anymore
top_sta
toe_sta



% check for case where the toe of slope is below SWL-1.5*H0
% in this case interpolate setup from the setupAtToe(really setup as first station), and the max setup
% also un-include points seaward of SWL-1.5*H0
if Ztoe > dep(1)
   dd=SWEL_fore-dep;
   k=find(dd<0,1); % k is index of first land point
   staAtSWL=interp1(dep(k-1:k),sta(k-1:k),SWEL_fore);   
   dsta=staAtSWL-sta(1);
   dsetup=maxSetup-setupAtToe;
   dsetdsta=dsetup/dsta;
   setup=setupAtToe+dsetdsta*(toe_sta-sta(1));
   sprintf('-!!- Location of SWEL-1.5*H0 is %4.1f ft landward of toe of slope',dsta)
   sprintf('-!!- Setup is interpolated between setup at toe of slope and max setup')
   sprintf('-!!-       setup is adjusted to %4.2f feet',setup)
   SWEL=SWEL-setupAtToe+setup;
   sprintf('-!!-       SWEL is adjusted to %4.2f feet',SWEL)
   k=find(dep < SWEL-1.5*H0)
   sta(k)=[];
   dep(k)=[];
else
   sprintf('-!!- The User has selected a starting point that is %4.2f feet above the elevation of SWEL-1.5H0\n',dep(1)-Ztoe)
   sprintf('-!!- This may be reasonable for some cases.  However the user may want to consider:\n')
   sprintf('-!!-   1) Selecting a starting point that is at or below %4.2f feet elevation, or\n', Ztoe)
   sprintf('-!!-   2) Reducing the incident wave height to a depth limited condition.\n')
end



% now iterate converge on a runup elevation
tol=0.01;  % convergence criteria
R2del=999;
R2_new=3*H0; %initial guess
R2=R2_new;
iter=0;
R2_all=[];
topStaAll=[];
Berm_Segs=[];
TAW_ALWAYS_VALID=1;
while(abs(R2del) > tol && iter <= 25)
    iter=iter+1;
    sprintf ('!----------------- STARTING ITERATION %d -----------------!',iter)
    % elevation of toe of slope
    Ztoe
    % station of toe slope (relative to 0-NAVD88 shoreline
    toe_sta 
    % station of top of slope/extent of 2% run-up
    top_sta
    % elevation of top of slope/extent of 2% run-up
    Z2
    % incident significant wave height
    H0
    % incident spectral peak wave period
    Tp
    % incident spectral mean wave period
    T0
        
    R2=R2_new  
    Z2=R2+SWEL
    % determine slope for this iteration
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end)
    end
    
    % get the length of the slope (not accounting for berm)
    Lslope=top_sta-toe_sta

    % loop over profile segments to determine berm factor 
    % re-calculate influence of depth of berm based on this run-up elevation
    % check for berm, berm width, berm height
    berm_width=0;
    rdh_sum=0;
    Berm_Segs=[];
    Berm_Heights=[];
    for kk=1:length(sta)-1
       ddep=dep(kk+1)-dep(kk);
       dsta=sta(kk+1)-sta(kk);
       s=ddep/dsta;
       if  (s < 1/15)     % count it as a berm if slope is flatter than 1:15 (see TAW manual)
          sprintf ('Berm Factor Calculation: Iteration %d, Profile Segment: %d',iter,kk) 
          berm_width=berm_width+dsta;   % tally the width of all berm segments
          % compute the rdh for this segment and weight it by the segment length
          dh=SWEL-(dep(kk)+dep(kk+1))/2  
          if  dh < 0
              chi=R2;
          else 
              chi=2* H0;
          end
          if (dh <= R2 & dh >=-2*H0)
             rdh=(0.5-0.5*cos(3.14159*dh/chi)) ;
          else
             rdh=1;
          end
          rdh_sum=rdh_sum + rdh * dsta
          Berm_Segs=[Berm_Segs, kk];
          Berm_Heights=[Berm_Heights, (dep(kk)+dep(kk+1))/2];
       end
       if dep(kk) >= Z2  % jump out of loop if we reached limit of run-up for this iteration 
          break
       end
    end
    sprintf ('!------- End Berm Factor Calculation, Iter: %d ---------!',iter)
    berm_width
    rB=berm_width/Lslope 
    if (berm_width > 0)
       rdh_mean=rdh_sum/berm_width
    else
       rdh_mean=1
    end
    gamma_berm=1- rB * (1-rdh_mean)
    if gamma_berm > 1
       gamma_berm=1
    end
    if gamma_berm < 0.6
       gamma_berm =0.6
    end
    % Iribarren number
    slope=(Z2-Ztoe)/(Lslope-berm_width)
    Irb=(slope/(sqrt(H0/L0)))
    % runup height
    gamma_berm
    gamma_perm
    gamma_beta
    gamma_rough
    gamma=gamma_berm*gamma_perm*gamma_beta*gamma_rough

    % check validity
    TAW_VALID=1;
    if (Irb*gamma_berm < 0.5 | Irb*gamma_berm > 10 )
       sprintf('!!! - - Iribaren number: %6.2f is outside the valid range (0.5-10), TAW NOT VALID - - !!!\n', Irb*gamma_berm)
       TAW_VALID=0;
    else
       sprintf('!!! - - Iribaren number: %6.2f is in the valid range (0.5-10), TAW RECOMMENDED - - !!!\n', Irb*gamma_berm)
    end
    islope=1/slope;  
    if (slope < 1/8 | slope > 1)
       sprintf('!!! - - slope: 1:%3.1f V:H is outside the valid range (1:8 - 1:1), TAW NOT VALID - - !!!\n', islope)
       TAW_VALID=0;
    else
       sprintf('!!! - - slope: 1:%3.1f V:H is in the valid range (1:8 - 1:1), TAW RECOMMENDED - - !!!\n', islope)
    end
    if TAW_VALID == 0
       TAW_ALWAYS_VALID=0;
    end
         
    if (Irb*gamma_berm < 1.8)
       R2_new=gamma*H0*1.77*Irb
    else
       R2_new=gamma*H0*(4.3-(1.6/sqrt(Irb))) 
    end


    % check to see if we need to evaluate a shallow foreshore
    if berm_width > 0.25 * L0;
       disp ('!   Berm_width is greater than 1/4 wave length')
       disp ('!   Runup will be weighted average with foreshore calculation assuming depth limited wave height on berm');
       % do the foreshore calculation 
       fore_H0=0.78*(SWEL_fore-min(Berm_Heights))
       % get upper slope
       fore_toe_sta=-999;
       fore_toe_dep=-999;
       for kk=length(dep)-1:-1:1
          ddep=dep(kk+1)-dep(kk);
          dsta=sta(kk+1)-sta(kk);
          s=ddep/dsta;
          if s < 1/15
             break
          end
          fore_toe_sta=sta(kk);
          fore_toe_dep=dep(kk);
          upper_slope=(Z2-fore_toe_dep)/(top_sta-fore_toe_sta)
       end
       fore_Irb=upper_slope/(sqrt(fore_H0/L0));
       fore_gamma=gamma_perm*gamma_beta*gamma_rough;
       if (fore_Irb < 1.8)
          fore_R2=fore_gamma*fore_H0*1.77*fore_Irb;
       else
          fore_R2=fore_gamma*fore_H0*(4.3-(1.6/sqrt(fore_Irb))); 
       end
       if berm_width >= L0
          R2_new=fore_R2
          disp ('berm is wider than one wavelength, use full shallow foreshore solution');
       else
          w2=(berm_width-0.25*L0)/(0.75*L0)
          w1=1-w2
          R2_new=w2*fore_R2 + w1*R2_new
       end
    end % end berm width check

    % convergence criterion
    R2del=abs(R2-R2_new)
    R2_all(iter)=R2_new;

    % get the new top station (for plot purposes)
    Z2=R2_new+SWEL
    top_sta=-999;
    for kk=1:length(sta)-1
       if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
          top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
          break;
       end
    end
    if top_sta==-999
       dy=Z2-dep(end);
       top_sta=sta(end)+dy/S(end);
    end
    topStaAll(iter)=top_sta;

end

% final 2% runup elevation
Z2=R2_new+SWEL

diary off


% make some plots


% get the full transect



file=xlsread('../ADCIRC_returns/YK-110XYZSTA_RETURNS.csv');
sta_full=file(:,4);
ele_full=file(:,3);


%determine SWL station
SWL_sta=-999;
for kk=1:length(sta)-1
    if ((SWEL > dep(kk)) & (SWEL <= dep(kk+1)))    % here is the intersection of z2 with profile
       SWL_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),SWEL)
       break;
    end
end

figure('position', [100 100 900 900])
hold on
plot(sta_full,ele_full,'linewidth',2,'color',[0.5 0.2 0.1]);
plot(sta,dep,'md');
grid on
plot([min(sta),SWL_sta],[SWEL SWEL],'color','b','linewidth',2)
%plot([berm_sta0,berm_sta0+berm_width],[berm_height berm_height],'k:','linewidth',3);
%scatter(sta_org(Berm_Segs),dep_org(Berm_Segs),30,'b*')
plot([min(sta),topStaAll(end)],[Z2, Z2],'g','linewidth',2);
scatter(topStaAll,R2_all+SWEL,30,'r*');
%scatter(toe_sta,Db,70,'k*');
scatter(toe_sta,Ztoe,70,'k*');
% plot berm segs
for nn=1:length(Berm_Segs)
   seg=Berm_Segs(nn);
   h=Berm_Heights(nn);
   plot([sta(seg) sta(seg+1)],[h h],'y:','linewidth',3);
end
xlim([toe_sta-50 max(topStaAll)+50])
grid minor
hleg=legend('Ground','TAW Segments','SWEL','Final 2% Run-up','R2 iterations','Toe|SWEL-1.5*H0','Berm Segments');
set (hleg,'location','southeast');
xlabel('Distance from zero NAVD88 shoreline (feet)');
ylabel('Elevation (Feet-NAVD88)');
str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
yylm=get (gca,'ylim');
dztxt=(yylm(end)-yylm(1))/15;
text(toe_sta,Z2+dztxt/4,str);
str=sprintf('Run-up Height is %5.1f Feet',R2_new);
text(toe_sta,Z2-dztxt/4,str);
title (plotTitle);
set(gcf,'color','w');
print('-r450','-dpdf',imgname);

-1.000000e+00
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S13');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U13');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O13');
