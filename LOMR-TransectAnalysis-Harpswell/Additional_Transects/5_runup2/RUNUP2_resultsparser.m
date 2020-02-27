
fdir='RUNUP2_files\';
fname='YK-5.OUT';
filename=[fdir fname];
fid=fopen(filename,'r');

% Register li_count to the first line of the output table
pre_li=0;
while ~feof(fid)
    li=fgetl(fid);
    
    if length(li)==56 && strcmp(li(45:47),'OUT')
        li_count=1;
        break
    end
    pre_li=pre_li+1;
end

% Run through the output table and collect results
r_li=0;
while ~feof(fid)

    li=fgetl(fid);
    A=str2num(li);
    li_count=li_count+1;
    if  isempty(A) && ~isempty(li) && li_count>11
        messages{r_li}=li;
        continue
    end
    if isempty(A)
        continue
    end
    r_li=r_li+1;
    swel(r_li)=A(1);
    H_m(r_li)=A(2);
    per(r_li)=A(3);
    break_slope(r_li)=A(4);
    runup_slope(r_li)=A(5);
    runup(r_li)=A(6);
        
end
fclose(fid)

% Calculate stats on runup values
R_m=mean(runup);
R_std=std(runup);

% Check for the presence of messages, if they exist, flag them
exist messages 'var'
if ans
    messageflag=1;
else
    messageflag=0;
end



