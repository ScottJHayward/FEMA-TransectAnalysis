close all;
clear all;
% 20190912
%%%%%%%%%%%%%%%%%%   config   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datafile='../data/transectdata.xls'; 
tDIR='../ADCIRC_returns/'; %location of transects
prefix='CM-';
%%%%%%%%%%%%%%%%%%   end config   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,txt,raw]=xlsread(datafile);
for i=2:size(raw,1)
fnames{i-1}=raw{i,1};
end
setup=num(:,11);
twl=num(:,2)+setup;
hs=num(:,12);
per=num(:,13);
toex=num(:,6);
toez=num(:,7);
topx=num(:,8);
topz=num(:,9);
startx=num(:,1)

for i=1:length(fnames)
  fname=[tDIR fnames{i} 'XYZSTA_RETURNS.csv'];
  file=xlsread(fname);
  sta_full{i}=file(:,4);ele_full{i}=file(:,3)
end


titles = {'100-year WHAFIS Output'}

for i=1:length(fnames)

   infile=[fnames{i} '_parsed.csv'];
   imname=[fnames{i} 'CrestElev1.pdf'];

  [card,sta,easting,northing,lon,lat,elev,swl,hc,crest]=textread(infile,'%s%n%n%n%n%n%n%n%n%n%*n','headerlines',1,'delimiter',',');

  sta3=sta;
  k=find (strcmp('AS',card));
  swl(k)=nan;
  crest(k)=nan;
  sta3(k)=nan;


  sta2=[min(sta):.1:max(sta)];
  crest2=interp1(sta,crest,sta2);
  hc2=interp1(sta,hc,sta2);
  bfe=round(crest2);
  vzone=NaN(size(bfe));
  vzone(find( hc2 >=3 ))=bfe(find( hc2 >=3 ));
  
  j=find(abs(sta_full{i}-startx(i))==min(abs(sta_full{1,i}-startx(i))));
  shift=sta_full{i}(j);   % this is shift for full transect profile to start of WHAFIS transect
  
  % find the zero station in WHAFIS
  for k=1:length(sta)-1
      if (elev(k) < 0) & (elev(k+1)) > 0
         sta0=interp1(elev(k:k+1),sta(k:k+1),0);
         lon0=interp1(elev(k:k+1),lon(k:k+1),0);
         lat0=interp1(elev(k:k+1),lat(k:k+1),0);
         break
      end
  end   

  % calc heading and write orientation string
  theta=atan2d((northing(end)-northing(1)),(easting(end)-easting(1)))
  orientStr1=sprintf('Zero Station: %13.8f, %13.8f',lon0,lat0);
  orientStr2=sprintf('Onshore Dir: %5.1f deg CCW from E',theta);

  figure('position', [100 100 1000 600])  
  hold on 
  plot (sta_full{i}-shift-sta0,ele_full{i},'color',[0.5 0.2 0.1],'linewidth',1);
  %plot (sta,elev,'color',[0.5 0.2 0.1],'linewidth',2);
  plot (sta2-sta0,vzone,'y','linewidth',5);
  plot (sta3-sta0,swl,'b:','linewidth',2)
  plot (sta3-sta0,crest,'g-.','linewidth',2);
  plot (sta2-sta0,bfe,'r-.','linewidth',2);

  
 % k1=find(sta<=VG_bnd1(i))
 % k1=max(k1);
 % k2=find(sta<=VG_bnd2(i))
 % k2=max(k2);
 % plot(sta(k1:k2),elev(k1:k2),'g*'); 
  
  grid on
  xlabel ('Station (ft)','fontsize',16);
  ylabel ('Elevation (ft-NAVD88)','fontsize',16);
  tstr=[fnames{i},titles,orientStr1,orientStr2]
  title(tstr,'fontsize',16)
  ylim([-5 20]);
  hleg=legend('Ground','V-zone','SWEL + Setup','Wave Crest Envelope','BFE','location','best');
  set (hleg,'fontsize',16);
  set(gca,'XMinorGrid','on') 
   %set(gca,'YMinorGrid','on') 
  % draw minor grid lines with 1 ft spacing
  xx=[0 max(sta)*2]-sta0;
  for yy=-30:30 
    plot(xx,[yy yy],'k:');
  end
  xlim([0 max(sta)]-sta0);  

  set (gcf,'color','w');
  set (gcf,'paperOrientation','landscape');
  print('-r450','-dpdf',imname);
  
end 
 

