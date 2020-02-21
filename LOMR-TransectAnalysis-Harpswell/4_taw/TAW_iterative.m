clear all
close all
format long g

diary logfiles/CM-122-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-122-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-122-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-122-1-runup';                            
SWEL=9.0674;  % 100-yr still water level including wave setup. 
H0=2.1741;    % significant wave height at toe of structure 
Tp=5.0338;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.00035761; 
maxSetup=0.050741;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-122-1'


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



file=xlsread('../ADCIRC_returns/CM-122-1XYZSTA_RETURNS.csv');
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

diary logfiles/CM-123-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-123
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-123sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-123-runup';                            
SWEL=9.0414;  % 100-yr still water level including wave setup. 
H0=3.9604;    % significant wave height at toe of structure 
Tp=5.1353;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.9877;   % this may get changed automatically below
gamma_rough=0.9;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.013638; 
maxSetup=0.5482;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-123'


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



file=xlsread('../ADCIRC_returns/CM-123XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S3');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U3');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O3');
clear all
close all
format long g

diary logfiles/CM-123-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-123-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-123-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-123-1-runup';                            
SWEL=9.023;  % 100-yr still water level including wave setup. 
H0=1.0126;    % significant wave height at toe of structure 
Tp=5.1017;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.17342; 
maxSetup=0.17342;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-123-1'


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



file=xlsread('../ADCIRC_returns/CM-123-1XYZSTA_RETURNS.csv');
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

diary logfiles/CM-124-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-124
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-124sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-124-runup';                            
SWEL=9.0068;  % 100-yr still water level including wave setup. 
H0=4.0078;    % significant wave height at toe of structure 
Tp=5.0362;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.93753;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.01648; 
maxSetup=0.27738;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-124'


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



file=xlsread('../ADCIRC_returns/CM-124XYZSTA_RETURNS.csv');
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

diary logfiles/CM-124-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-124-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-124-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-124-1-runup';                            
SWEL=8.9825;  % 100-yr still water level including wave setup. 
H0=1.7682;    % significant wave height at toe of structure 
Tp=3.343;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.026568; 
maxSetup=0.089504;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-124-1'


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



file=xlsread('../ADCIRC_returns/CM-124-1XYZSTA_RETURNS.csv');
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

diary logfiles/CM-124-2-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-124-2
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-124-2sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-124-2-runup';                            
SWEL=8.9775;  % 100-yr still water level including wave setup. 
H0=2.6451;    % significant wave height at toe of structure 
Tp=3.3905;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.6;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.022201; 
maxSetup=0.29117;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-124-2'


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



file=xlsread('../ADCIRC_returns/CM-124-2XYZSTA_RETURNS.csv');
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

diary logfiles/CM-126-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-126-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-126-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-126-1-runup';                            
SWEL=8.8944;  % 100-yr still water level including wave setup. 
H0=3.2859;    % significant wave height at toe of structure 
Tp=11.1501;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.99698;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.084521; 
maxSetup=0.32708;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-126-1'


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



file=xlsread('../ADCIRC_returns/CM-126-1XYZSTA_RETURNS.csv');
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

diary logfiles/CM-127-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-127
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-127sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-127-runup';                            
SWEL=8.8651;  % 100-yr still water level including wave setup. 
H0=6.5053;    % significant wave height at toe of structure 
Tp=9.8292;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.83689;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.037592; 
maxSetup=0.81296;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-127'


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



file=xlsread('../ADCIRC_returns/CM-127XYZSTA_RETURNS.csv');
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

diary logfiles/CM-127-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-127-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-127-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-127-1-runup';                            
SWEL=8.8666;  % 100-yr still water level including wave setup. 
H0=3.3682;    % significant wave height at toe of structure 
Tp=6.1803;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.82539;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.0097211; 
maxSetup=0.46254;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-127-1'


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



file=xlsread('../ADCIRC_returns/CM-127-1XYZSTA_RETURNS.csv');
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

diary logfiles/CM-129-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-129-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-129-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-129-1-runup';                            
SWEL=8.8189;  % 100-yr still water level including wave setup. 
H0=5.7154;    % significant wave height at toe of structure 
Tp=12.3322;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.6;   % this may get changed automatically below
gamma_rough=0.6;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.018871; 
maxSetup=0.45369;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-129-1'


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



file=xlsread('../ADCIRC_returns/CM-129-1XYZSTA_RETURNS.csv');
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

diary logfiles/CM-130-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-130
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-130sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-130-runup';                            
SWEL=8.8141;  % 100-yr still water level including wave setup. 
H0=7.5866;    % significant wave height at toe of structure 
Tp=13.6829;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.804;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.01352; 
maxSetup=0.36725;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-130'


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



file=xlsread('../ADCIRC_returns/CM-130XYZSTA_RETURNS.csv');
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

diary logfiles/CM-131-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-131-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-131-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-131-1-runup';                            
SWEL=8.8473;  % 100-yr still water level including wave setup. 
H0=1.4239;    % significant wave height at toe of structure 
Tp=3.7042;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.99896;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.1644; 
maxSetup=0.2308;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-131-1'


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



file=xlsread('../ADCIRC_returns/CM-131-1XYZSTA_RETURNS.csv');
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
clear all
close all
format long g

diary logfiles/CM-133-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-133-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-133-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-133-1-runup';                            
SWEL=8.9057;  % 100-yr still water level including wave setup. 
H0=1.8629;    % significant wave height at toe of structure 
Tp=2.655;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.0085433; 
maxSetup=0.017969;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-133-1'


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



file=xlsread('../ADCIRC_returns/CM-133-1XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S14');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U14');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O14');
clear all
close all
format long g

diary logfiles/CM-132-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-132
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-132sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-132-runup';                            
SWEL=8.8527;  % 100-yr still water level including wave setup. 
H0=3.4596;    % significant wave height at toe of structure 
Tp=4.6732;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.017657; 
maxSetup=0.45563;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-132'


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



file=xlsread('../ADCIRC_returns/CM-132XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S15');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U15');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O15');
clear all
close all
format long g

diary logfiles/CM-133-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-133
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-133sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-133-runup';                            
SWEL=8.8742;  % 100-yr still water level including wave setup. 
H0=5.3885;    % significant wave height at toe of structure 
Tp=7.957;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.97614;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.011821; 
maxSetup=0.63595;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-133'


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



file=xlsread('../ADCIRC_returns/CM-133XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S16');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U16');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O16');
clear all
close all
format long g

diary logfiles/CM-134-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-134
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-134sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-134-runup';                            
SWEL=8.9062;  % 100-yr still water level including wave setup. 
H0=3.4461;    % significant wave height at toe of structure 
Tp=4.9744;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.96945;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.028737; 
maxSetup=0.28085;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-134'


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



file=xlsread('../ADCIRC_returns/CM-134XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S17');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U17');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O17');
clear all
close all
format long g

diary logfiles/CM-134-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-134-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-134-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-134-1-runup';                            
SWEL=8.8995;  % 100-yr still water level including wave setup. 
H0=2.5462;    % significant wave height at toe of structure 
Tp=3.2843;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.0065354; 
maxSetup=0.21697;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-134-1'


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



file=xlsread('../ADCIRC_returns/CM-134-1XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S18');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U18');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O18');
clear all
close all
format long g

diary logfiles/CM-135-2-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-135-2
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-135-2sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-135-2-runup';                            
SWEL=8.8436;  % 100-yr still water level including wave setup. 
H0=0.65495;    % significant wave height at toe of structure 
Tp=2.0989;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.26319; 
maxSetup=0.28398;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-135-2'


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



file=xlsread('../ADCIRC_returns/CM-135-2XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S20');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U20');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O20');
clear all
close all
format long g

diary logfiles/CM-136-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-136
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-136sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-136-runup';                            
SWEL=8.8177;  % 100-yr still water level including wave setup. 
H0=3.393;    % significant wave height at toe of structure 
Tp=9.8875;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.020909; 
maxSetup=0.53172;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-136'


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



file=xlsread('../ADCIRC_returns/CM-136XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S21');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U21');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O21');
clear all
close all
format long g

diary logfiles/CM-137-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-137
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-137sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-137-runup';                            
SWEL=8.8077;  % 100-yr still water level including wave setup. 
H0=4.9302;    % significant wave height at toe of structure 
Tp=9.6953;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.94202;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.0016634; 
maxSetup=0.73897;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-137'


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



file=xlsread('../ADCIRC_returns/CM-137XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S22');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U22');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O22');
clear all
close all
format long g

diary logfiles/CM-138-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-138
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-138sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-138-runup';                            
SWEL=8.7974;  % 100-yr still water level including wave setup. 
H0=8.8182;    % significant wave height at toe of structure 
Tp=11.5883;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.7418;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.14678; 
maxSetup=1.2538;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-138'


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



file=xlsread('../ADCIRC_returns/CM-138XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S23');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U23');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O23');
clear all
close all
format long g

diary logfiles/CM-139-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-139
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-139sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-139-runup';                            
SWEL=8.804;  % 100-yr still water level including wave setup. 
H0=6.2267;    % significant wave height at toe of structure 
Tp=10.2103;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.7549;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.023839; 
maxSetup=0.87295;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-139'


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



file=xlsread('../ADCIRC_returns/CM-139XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S24');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U24');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O24');
clear all
close all
format long g

diary logfiles/CM-139-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-139-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-139-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-139-1-runup';                            
SWEL=8.7974;  % 100-yr still water level including wave setup. 
H0=9.1162;    % significant wave height at toe of structure 
Tp=9.9055;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.071548; 
maxSetup=1.804;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-139-1'


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



file=xlsread('../ADCIRC_returns/CM-139-1XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S25');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U25');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O25');
clear all
close all
format long g

diary logfiles/CM-139-2-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-139-2
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-139-2sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-139-2-runup';                            
SWEL=8.7973;  % 100-yr still water level including wave setup. 
H0=7.609;    % significant wave height at toe of structure 
Tp=9.0889;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.94269;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.052129; 
maxSetup=1.14;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-139-2'


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



file=xlsread('../ADCIRC_returns/CM-139-2XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S26');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U26');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O26');
clear all
close all
format long g

diary logfiles/CM-140-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-140
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-140sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-140-runup';                            
SWEL=8.7522;  % 100-yr still water level including wave setup. 
H0=10.9296;    % significant wave height at toe of structure 
Tp=14.4675;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.93231;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=1.4979; 
maxSetup=2.5278;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-140'


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



file=xlsread('../ADCIRC_returns/CM-140XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S27');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U27');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O27');
clear all
close all
format long g

diary logfiles/CM-141-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-141
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-141sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-141-runup';                            
SWEL=8.7636;  % 100-yr still water level including wave setup. 
H0=6.4243;    % significant wave height at toe of structure 
Tp=15.2662;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.95946;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=1.1599; 
maxSetup=1.8064;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-141'


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



file=xlsread('../ADCIRC_returns/CM-141XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S28');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U28');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O28');
clear all
close all
format long g

diary logfiles/CM-142-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-142
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-142sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-142-runup';                            
SWEL=8.7833;  % 100-yr still water level including wave setup. 
H0=8.0408;    % significant wave height at toe of structure 
Tp=13.8504;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.90215;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.58856; 
maxSetup=1.3427;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-142'


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



file=xlsread('../ADCIRC_returns/CM-142XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S29');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U29');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O29');
clear all
close all
format long g

diary logfiles/CM-145-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-145
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-145sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-145-runup';                            
SWEL=8.8099;  % 100-yr still water level including wave setup. 
H0=8.5251;    % significant wave height at toe of structure 
Tp=11.4911;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.95254;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.03793; 
maxSetup=1.1002;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-145'


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



file=xlsread('../ADCIRC_returns/CM-145XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S30');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U30');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O30');
clear all
close all
format long g

diary logfiles/CM-149-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-149
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-149sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-149-runup';                            
SWEL=8.8429;  % 100-yr still water level including wave setup. 
H0=9.9977;    % significant wave height at toe of structure 
Tp=13.6834;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.91997;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.087333; 
maxSetup=1.6597;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-149'


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



file=xlsread('../ADCIRC_returns/CM-149XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S31');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U31');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O31');
clear all
close all
format long g

diary logfiles/CM-149-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-149-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-149-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-149-1-runup';                            
SWEL=8.8764;  % 100-yr still water level including wave setup. 
H0=5.21;    % significant wave height at toe of structure 
Tp=10.9082;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.062926; 
maxSetup=0.33328;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-149-1'


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



file=xlsread('../ADCIRC_returns/CM-149-1XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S32');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U32');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O32');
clear all
close all
format long g

diary logfiles/CM-150-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-150
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-150sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-150-runup';                            
SWEL=8.8664;  % 100-yr still water level including wave setup. 
H0=7.6185;    % significant wave height at toe of structure 
Tp=13.6394;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.98769;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.0071522; 
maxSetup=0;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-150'


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



file=xlsread('../ADCIRC_returns/CM-150XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S33');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U33');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O33');
clear all
close all
format long g

diary logfiles/CM-150-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-150-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-150-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-150-1-runup';                            
SWEL=8.9195;  % 100-yr still water level including wave setup. 
H0=1.6621;    % significant wave height at toe of structure 
Tp=2.2768;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.0022605; 
maxSetup=0.082149;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-150-1'


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



file=xlsread('../ADCIRC_returns/CM-150-1XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S35');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U35');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O35');
clear all
close all
format long g

diary logfiles/CM-150-2-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-150-2
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-150-2sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-150-2-runup';                            
SWEL=8.9231;  % 100-yr still water level including wave setup. 
H0=3.1328;    % significant wave height at toe of structure 
Tp=8.6753;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.025892; 
maxSetup=0.27652;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-150-2'


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



file=xlsread('../ADCIRC_returns/CM-150-2XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S36');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U36');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O36');
clear all
close all
format long g

diary logfiles/CM-151-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-151-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-151-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-151-1-runup';                            
SWEL=8.9177;  % 100-yr still water level including wave setup. 
H0=1.8096;    % significant wave height at toe of structure 
Tp=2.2995;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.00097113; 
maxSetup=0.082529;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-151-1'


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



file=xlsread('../ADCIRC_returns/CM-151-1XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S37');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U37');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O37');
clear all
close all
format long g

diary logfiles/CM-158-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-158-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-158-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-158-1-runup';                            
SWEL=8.8157;  % 100-yr still water level including wave setup. 
H0=1.7049;    % significant wave height at toe of structure 
Tp=2.6702;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.6;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.010217; 
maxSetup=8.2021e-05;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-158-1'


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



file=xlsread('../ADCIRC_returns/CM-158-1XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S38');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U38');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O38');
clear all
close all
format long g

diary logfiles/CM-159-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-159
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-159sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-159-runup';                            
SWEL=8.8217;  % 100-yr still water level including wave setup. 
H0=6.5832;    % significant wave height at toe of structure 
Tp=14.6161;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.99728;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.036975; 
maxSetup=1.0255;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-159'


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



file=xlsread('../ADCIRC_returns/CM-159XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S39');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U39');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O39');
clear all
close all
format long g

diary logfiles/CM-159-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-159-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-159-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-159-1-runup';                            
SWEL=8.811;  % 100-yr still water level including wave setup. 
H0=5.2342;    % significant wave height at toe of structure 
Tp=13.5919;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.6;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.01229; 
maxSetup=0.75187;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-159-1'


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



file=xlsread('../ADCIRC_returns/CM-159-1XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S40');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U40');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O40');
clear all
close all
format long g

diary logfiles/CM-161-1-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: CM-161-1
% calculation by SJH, Ransom Consulting, Inc. 20-Feb-2020
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
fname='inpfiles/CM-161-1sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/CM-161-1-runup';                            
SWEL=8.874;  % 100-yr still water level including wave setup. 
H0=1.7143;    % significant wave height at toe of structure 
Tp=2.1869;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.0060892; 
maxSetup=0;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for CM-161-1'


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



file=xlsread('../ADCIRC_returns/CM-161-1XYZSTA_RETURNS.csv');
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
xlswrite('../data/transectdata.xls',Z2,'Sheet1','S41');
xlswrite('../data/transectdata.xls',TAW_ALWAYS_VALID,'Sheet1','U41');
xlswrite('../data/transectdata.xls',gamma_berm,'Sheet1','O41');
