clear all
close all
format long g

diary logfiles/YK-05-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: YK-05
% calculation by SJH, Ransom Consulting, Inc. 19-Feb-2020
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
fname='inpfiles/YK-05sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-05-runup';                            
SWEL=9.0268;  % 100-yr still water level including wave setup. 
H0=3.5425;    % significant wave height at toe of structure 
Tp=6.337;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.98563;   % this may get changed automatically below
gamma_rough=0.75;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.040102; 
maxSetup=0.0028839;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-05'


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



file=xlsread('../ADCIRC_returns/YK-05XYZSTA_RETURNS.csv');
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

diary logfiles/YK-06F-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: YK-06F
% calculation by SJH, Ransom Consulting, Inc. 19-Feb-2020
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
fname='inpfiles/YK-06Fsta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-06F-runup';                            
SWEL=9.0235;  % 100-yr still water level including wave setup. 
H0=5.4882;    % significant wave height at toe of structure 
Tp=9.7138;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.96447;   % this may get changed automatically below
gamma_rough=0.6;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.02834; 
maxSetup=0.62428;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-06F'


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



file=xlsread('../ADCIRC_returns/YK-06FXYZSTA_RETURNS.csv');
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

diary logfiles/YK-06-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: YK-06
% calculation by SJH, Ransom Consulting, Inc. 19-Feb-2020
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
fname='inpfiles/YK-06sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-06-runup';                            
SWEL=9.0235;  % 100-yr still water level including wave setup. 
H0=5.4588;    % significant wave height at toe of structure 
Tp=9.7161;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.6;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.028035; 
maxSetup=0.73082;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-06'


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



file=xlsread('../ADCIRC_returns/YK-06XYZSTA_RETURNS.csv');
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

diary logfiles/YK-07-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: YK-07
% calculation by SJH, Ransom Consulting, Inc. 19-Feb-2020
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
fname='inpfiles/YK-07sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-07-runup';                            
SWEL=9.0273;  % 100-yr still water level including wave setup. 
H0=3.4318;    % significant wave height at toe of structure 
Tp=6.9867;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=1;   % this may get changed automatically below
gamma_rough=0.85;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=-0.02211; 
maxSetup=0;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-07'


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



file=xlsread('../ADCIRC_returns/YK-07XYZSTA_RETURNS.csv');
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

diary logfiles/YK-14-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: YK-14
% calculation by SJH, Ransom Consulting, Inc. 19-Feb-2020
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
fname='inpfiles/YK-14sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-14-runup';                            
SWEL=9.19;  % 100-yr still water level including wave setup. 
H0=4.9688;    % significant wave height at toe of structure 
Tp=13.8709;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.93233;   % this may get changed automatically below
gamma_rough=0.8;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.5381; 
maxSetup=1.1359;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-14'


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



file=xlsread('../ADCIRC_returns/YK-14XYZSTA_RETURNS.csv');
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

diary logfiles/YK-15-DIARY.txt  % open a diary file to record calculations
diary on         % begin recording

% FEMA appeal for The Town of Harpswell, Cumberland county, Maine
% TRANSECT ID: YK-15
% calculation by SJH, Ransom Consulting, Inc. 19-Feb-2020
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
fname='inpfiles/YK-15sta_ele_include.csv';   % file with station, elevation, include
                                      % third colunm is 0 for excluded points
imgname='logfiles/YK-15-runup';                            
SWEL=9.2819;  % 100-yr still water level including wave setup. 
H0=6.2414;    % significant wave height at toe of structure 
Tp=12.7769;    % peak period, 1/fma, 
T0=Tp/1.1;   

gamma_berm=0.97227;   % this may get changed automatically below
gamma_rough=1;  
gamma_beta=1;  
gamma_perm=1;

setupAtToe=0.18011; 
maxSetup=1.0731;   % only used in case of berm/shallow foreshore weighted average 

plotTitle='Iterative TAW for YK-15'


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



file=xlsread('../ADCIRC_returns/YK-15XYZSTA_RETURNS.csv');
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
