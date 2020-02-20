
close all;
clear all;

cd 'C:\0_PROJECTS\171.06013-PNSY-Kittery-FEMA\Modeling\Transects\SPM_vert';
 
tfile='C:\0_PROJECTS\171.06013-PNSY-Kittery-FEMA\Modeling\Transects\T2\T2.csv'

TWL=8.58;
toex=-5.588;
topx=1.223;
toez=-40
topz=8.7545;

R2=3.818;  % from SPM vertical wall method
Z2=R2+TWL;

Tname='T2';

slope=(topz-toez)/(topx-toex);

[sta,dep] = textread(tfile,'%n%n%*[^\n]','delimiter',',','headerlines',1);

 SWL_sta=(TWL-toez)/slope + toex;


 %  top_sta=-999;
 %  for kk=1:length(sta)-1
 %      if ((Z2 > dep(kk)) & (Z2 <= dep(kk+1)))    % here is the intersection of z2 with profile
 %         top_sta=interp1(dep(kk:kk+1),sta(kk:kk+1),Z2)
 %         break;
 %      end
 %  end
   top_sta=(Z2-toez)/slope + toex;
   if (Z2-topz >=3)  % splashoverzone
     top_sta=top_sta+30;
   end

   figure('position',[100 100 1800 700]);
   hold on
   plot(sta,dep,'linewidth',2,'color',[0.5 0.2 0.1]);
   grid on
   
   plot([min(sta),SWL_sta],[TWL TWL],'color','b','linewidth',2)
   plot([min(sta),top_sta],[Z2, Z2],'g','linewidth',2);
   plot([min(sta),top_sta],[round(Z2), round(Z2)],'r','linewidth',2);
   plot([toex topx],[toez topz],'y:','linewidth',2);
   grid on
   hleg=legend('Ground','TWL','2% Run-up','BFE','Approx. Slope');
   set (hleg,'fontsize',16,'location','southeast');
   xlabel('Distance from zero NAVD88 shoreline (feet)','fontsize',16);
   ylabel('Elevation (Feet-NAVD88)','fontsize',16);
   str=sprintf('Two Percent Run-up Elevation is %5.1f Feet-NAVD88',Z2);
   ylim([-5 15]);
   xlim([-50 150]);
   yylm=get (gca,'ylim');
   dztxt=(yylm(end)-yylm(1))/15;
   text(40,8-dztxt,str,'fontsize',16);
   str=sprintf('Run-up Height is %5.1f Feet',R2);
   text(40,8-2*dztxt,str,'fontsize',16);
   plotTitle=sprintf('%s 100-year SPM Vertical Slope, 2%% RUNUP',Tname); 
   title (plotTitle,'fontsize',16);
   set(gcf,'color','w');
   set(gca,'XMinorGrid','on') 
   set(gca,'YMinorGrid','on') 
   f=getframe(gcf);
   im=frame2im(f);
   pngname=sprintf('%s-SPM.png',Tname);
   imwrite(im,pngname,'png');
