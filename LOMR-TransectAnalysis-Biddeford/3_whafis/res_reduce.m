function [indices] = res_reduce(st,ele,begin_pt,end_pt,num)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%set up algorithm
indices=begin_pt:end_pt;
indices(find(ele(indices)>=25))=[];
i=isnan(ele(indices));
indices(i==1)=[];
x2=st(indices);
z2=ele(indices);
disp(['length: ' num2str(length(indices))])
err=[];
while length(x2)>num
    x1=x2;z1=z2; %reset x1 to equal x2, etc. 
    %a1=trapz(x1,z1);
    xc=x1(2:end-1);
    xl=x1(1:end-2);
    xr=x1(3:end);

    zc=z1(2:end-1);
    zl=z1(1:end-2);
    zr=z1(3:end);
    
    a1=(zl+zc).*(xc-xl) + (zc+zr).*(xr-xc);
    a2=(zl+zr).*(zr-zc);

    err=abs(a1-a2);

    [minerr,ff]=min(err);

    indices(ff)=[];%remove point with least error
    x2=st(indices);z2=ele(indices);    

%     plot(st,ele,'k',st(indices),ele(indices),'b');pause(.00001);
    disp(['length: ' num2str(length(indices))])
end
    
end

