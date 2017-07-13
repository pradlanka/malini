% Dynamic EC
% data is extracted time series from regions, size(data)=M x N x S, where M=number of time points, N is number of regions and S is number of subjects
 % if order (bic_smooth) and forgetting factor (ff_actual) are not
 % specified, we want the fastest variation in connectivity with least
 % memory as we want to know connectivity at every time instant. So use
 % order=1 and forgetting factor =1.
if ~exist('bic_smooth','var'),bic_smooth=1;end
if ~exist('ff_actual','var'),ff_actual=1;end
%-----------------------

% Case 1: Use the module below if you want connectivity for each individual
% subject and not a group value. 
% conn_individual is individuals dynamic EC matrix of size; M x N x N x S,in the N x N section, direction of connectivity is from rows to columns
if group==0
cd tsa

for sub=1:size(data,3)
    sub
    [x1,~]=mvaar1(squeeze(data(:,:,sub)),bic_smooth,ff_actual);
    
    temp1=reshape(x1,size(data,1),size(data,2),size(data,2),bic_smooth);

    conn_individual(:,:,:,sub)=squeeze(sum(temp1,4)); 
    clear x1  temp1;
end

cd ..
end

%------------------------

% Case 2: Use the module below if you DO NOT want connectivity for each individual
% subject and want a group value. 
% conn_group is group dynamic EC matrix of size; M x N x N,in the N x N section, direction of connectivity is from rows to columns
if gourp==1
cd tsa

[x2,e,Kalman2,Q2]=mvaar(squeeze(data(:,:,sub)),bic_smooth,ff_actual);

for sub=2:size(data,3)
    sub
    [x2,e,Kalman2,Q2]=mvaar(squeeze(data(:,:,sub)),bic_smooth,ff_actual,0,Kalman2);
end
    
temp2=reshape(x2,size(data,1),size(data,2),size(data,2),bic_smooth);
conn_group=squeeze(sum(temp2,4));


cd ..
end
