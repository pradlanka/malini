
function [bic_smooth,ff_actual]=find_order_ff(data)
%use this function to find order bic_smooth for given data and forgetting factor ff_actual.
%order is found using Bayesion information criterion, forgetting factor is
%optimized by minimizing variance of error' power 
% data is extracted time series from regions, size(data)=M x N x S, where M=number of time points, N is number of regions and S is number of subjects
 addpath(genpath('C:\Users\szg0071\Documents\MATLAB\gui'))

%   data = load('Sub_001.mat');
%   data = struct2cell(data);
%   data = cell2mat(data);
%   
count=1;
for k=1:size(data,2) %upto 190.
    for l=1:size(data,3)%size of this is 1.
        if (data(:,k,l)) ~= zeros(size(data,1),1)%the zeros part is 200 by 1 matrix.
        [bic(count),aic(count),bc,ac] = aks_find_model_order(squeeze(data(:,k,l))',2,20);
        count=count+1;
        end
    end
end

cd ('Codes\FC_EC_calculation code\DynamicEC\tsa')
        

bic_smooth=floor(mean(bic));

    for sub=1
        count2=1;  sub
        for UC=[10^-8 10^-7 10^-6 10^-5 10^-4 10^-3 10^-2 10^-1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]
            [~,e]=mvaar1(squeeze(data(:,:,sub)),bic_smooth,UC);    %mvaar1 is modified version of orginal mvaar.m and this version reduces the requirement of memory by half
            param(count2,sub,:)=var(e)./var(squeeze(data(:,:,sub)));
            count2=count2+1
        end
    end


cd ..


param_mean=squeeze(mean(param,2));


for k=1:size(param_mean,2)
    a(k)=find(squeeze(param_mean(:,k))==min(squeeze(param_mean(:,k))));
    ff(k)=ind2sub(size(squeeze(param_mean(:,k))),a(k));
end

ff_smooth=ff;



for k=1:length(ff_smooth)
    if ff_smooth(k)==1
        ff_actual(k)=10^-8;
    end
    if ff_smooth(k)==2
        ff_actual(k)=10^-7;
    end
    if ff_smooth(k)==3
        ff_actual(k)=10^-6;
    end
    if ff_smooth(k)==4
        ff_actual(k)=10^-5;
    end
    if ff_smooth(k)==5
        ff_actual(k)=10^-4;
    end
    if ff_smooth(k)==6
        ff_actual(k)=10^-3;
    end
    if ff_smooth(k)==7
        ff_actual(k)=10^-2;
    end
    if ff_smooth(k)==8
        ff_actual(k)=10^-1;
    end
    if ff_smooth(k)==9
        ff_actual(k)=0.2;
    end
    if ff_smooth(k)==10
        ff_actual(k)=0.3;
    end
    if ff_smooth(k)==11
        ff_actual(k)=0.4;
    end
    if ff_smooth(k)==12
        ff_actual(k)=0.5;
    end
    if ff_smooth(k)==13
        ff_actual(k)=0.6;
    end
    if ff_smooth(k)==14
        ff_actual(k)=0.7;
    end
    if ff_smooth(k)==15
        ff_actual(k)=0.8;
    end
    if ff_smooth(k)==16
        ff_actual(k)=0.9;
    end
end

save parameters bic_smooth ff_actual
    cd ..
    cd ..
    cd ..
    cd ..
