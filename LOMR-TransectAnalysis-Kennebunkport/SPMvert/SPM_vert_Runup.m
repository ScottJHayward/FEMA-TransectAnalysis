close all
clear all;

% runup for PNSY transect T2 on vertical wall
% following Atlantic and Gulf of Mexico Coastal Guidelines, 2007, D.2.8.1.4

% directory to work in
cd 'C:\0_PROJECTS\171.06013-PNSY-Kittery-FEMA\PostAppeal\SPM_vertical_wall';

% name of png file with SPM nomograph
pngname='SPM-Vert-Runup-NOMO.PNG';

% name for modified png with lines superimposed
outPng='T2-SPM-Runup_calc.png';

g=32.174; % ft/s/s

Hs=1.32 % significant wave height (ft)
Tp=1.9  % peak period (s)
ds=48.6   % depth at toe (ft)


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
str=sprintf('Sig. Waveht, H_s = %6.2f ft;      ==>     Mean WaveHt, H_o = 0.626*Hs = %6.2f ft',Hs,Ho);
text(25,ypix-125,str,'color','r')
str=sprintf('Peak Period, T_p = %6.2f sec;   ==>    Mean Period, T_m = 0.85*Tp = %6.2f sec ',Tp,T);
text(25,ypix-105,str,'color','r')
str=sprintf('x= H_o/g/T_m/T_m = %7.4f',x);
text(25,ypix-85,str,'color','r')
set(gcf,'color','w');

[x1,y1]=ginput(1);
plot([x1 x1],[0 ypix],'b-','linewidth',2)
 
[x2,y2]=ginput(1);
plot([0 xpix],[y2 y2],'r-','linewidth',2)





%calculate y from the pixel value 

pix=[439 57];     
dpix=pix(2)-pix(1)
Y=log10([0.1 6.0]);  
dY=Y(2)-Y(1);
l10y=Y(1)+(y2-pix(1))*dY/dpix;
y=10.^l10y;

% limit to 2 sig fig
y=0.1*round(10*y);
 

Ro=y*Ho

R2=2.2*Ro

str=sprintf('y= R/H_o = %4.1f',y);
text(25,ypix-60,str,'color','r')
str=sprintf('R = y*H_o = %6.3f ft',Ro);
text(25,ypix-40,str,'color','r')
str=sprintf('R_2_%% = 2.2*R = %6.3f ft',R2);
text(25,ypix-20,str,'color','r')

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
