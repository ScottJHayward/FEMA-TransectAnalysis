function [ierr]=plotTransectSources(transfile);
% HELP!!!



%transfile='YK-107XYZSTA.csv'

[x,y,z,sta,source]=textread(transfile,'%n%n%n%n%s','headerlines',1,'delimiter',',');

sources = unique(source)

figure
hold on
% plot(sta,z)

cmap=jet(length(sources))

for i=1:length(sources)
  plot(sta(strcmp(source,sources{i})),z(strcmp(source,sources{i})),'-o')
end


hleg=legend('line',sources')


