
function bic_smooth=find_order(data) 
Filepath=pwd;
%use this function to find order bic_smooth for given data and forgetting factor ff_actual.
%order is found using Bayesion information criterion, forgetting factor is
%optimized by minimizing variance of error' power 
% data is extracted time series from regions, size(data)=M x N x S, where M=number of time points, N is number of regions and S is number of subjects
addpath(genpath(Filepath))
count=1;
for i=1:size(data,2)
    for j=1:size(data,3)
        if(data(:,i,j)~=zeros(size(data,1), 1));
        [bic(count),aic(count),bc,ac] = aks_find_model_order(squeeze(data(:,i,j))',2,20);
        count=count+1;
        end
    end
end

bic_smooth=floor(mean(bic));
    
