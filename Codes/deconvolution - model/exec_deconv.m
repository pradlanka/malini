function [op HRF adjust_global PARA]=exec_deconv(ip)
TR=1;
% event = zeros(1,91);
HRF = zeros(30,91);
adjust_global = zeros(1,91);
PARA = zeros(3,91);

y=0;

na=zscore(ip);
z1=1;
for i=1:size(ip,2)
    if(na(:,i)~=0)
        y=1;
        subip(:,z1)=na(:,i);
        nzero1(z1,1)=i;z1=z1+1;
    end
end

clear i z

if(y==1)
[data_deconv , ~, HRF1, adjust_global1, PARA1] = wgr_deconv_canonhrf_par(subip,1,8/TR,TR); %8/TR=8/0.6=13.33
% [m n]=size(PAR);
% PARA=PAR(1:min(m,3),1:min(n,91));
else
    op=zeros(size(ip,1),size(ip,2));
    return
end

op=zeros(size(ip,1),size(ip,2));
for i=1:size(subip,2)
    op(:,nzero1(i))=data_deconv(:,i);
    HRF(:,nzero1(i))=HRF1(:,i);
    adjust_global(nzero1(i))=adjust_global1(i);
    PARA(:,nzero1(i))=PARA1(:,i);
end

clear i

end