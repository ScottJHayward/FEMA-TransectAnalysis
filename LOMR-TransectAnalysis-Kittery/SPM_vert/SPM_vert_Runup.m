close all
clear all;

% runup for PNSY transect CM-135-1 on vertical wall
% following Atlantic and Gulf of Mexico Coastal Guidelines, 2007, D.2.8.1.4

% name of png file with SPM nomograph
pngname='SPM-Vert-Runup-NOMO.png';

% name for modified png with lines superimposed
outPng='T2-SPM-Runup_calc.png';

g=32.174; % ft/s/s

Hs=2 % significant wave height (ft) at toe
Tp=3.47  % peak period (s) at SWAN modelinput... no output at toe- assume constant
twltoe = 8.658598813  % setup + SWL
ztoe = 5
ds = twltoe-ztoe   % depth at toe (ft)


% determine mean wave height and period (D.2.8.1.2.1)
Ho=0.626*Hs;
T=0.85*Tp;

% determine depth ration
dsOH= ds/Ho


% value for horizontal axis
x=Ho/g/T/T

% now the next few statements will make a plot open up
% first click on the plot corresponding the x value
% computed above, then click again on the intersection 
% with the line that is drawn and the appropriate line on 
% the plot for the value of dsOh given above. dsOh=3 if 
% dsOh > 3
%
%
im=imread(pngname);
ypix=size(im,1);
xpix=size(im,2);
figure('position',[100 100 100+xpix 100+ypix]);
image(im);
hold on;

% write the value of x on the plot
str=sprintf('Sig. Waveht, Hs = %6.2f ft;  Mean WaveHt, Ho = 0.626*Hs = %6.2f ft; d_s/Ho = %6.2f',Hs,Ho,dsOH);
text(25,ypix-70,str,'color','r')
str=sprintf('Peak Period, Tp = %6.2f ft;  Mean Period T = 0.85*Tp = %6.2f ',Tp,T);
text(25,ypix-50,str,'color','r')
str=sprintf('x= Ho/g/T/T = %7.4f',x);
text(25,ypix-30,str,'color','r')
set(gcf,'color','w');

[x1,y1]=ginput(1);
plot([x1 x1],[0 ypix],'b-','linewidth',2)
 
[x2,y2]=ginput(1);
plot([0 xpix],[y2 y2],'r-','linewidth',2)



%
%
%
%
%   FIRST STOP HERE 
%   
%   enter value read from Y-axis below 
%   then continue
%
%%

y=2.1

Ro=y*Ho

R2=2.2*Ro

str=sprintf('y= R/Ho = %4.1f',y);
text(650,ypix-70,str,'color','r')
str=sprintf('R = y*Ho = %6.3f',Ro);
text(650,ypix-50,str,'color','r')
str=sprintf('R2 = 2.2*R = %6.3f',R2);
text(650,ypix-30,str,'color','r')

f=getframe(gcf);
im2=frame2im(f);
imwrite(im2,outPng,'png');







%
%
%
%  last stop here
%
%
%
