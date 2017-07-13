function [Decision_Surf]=classifier_train(data,group_no,varargin)
echo off
if nargin > 2
    kernel_width= varargin{1};
else
    kernel_width = 25;
end
P = length(unique(group_no))-1; N=size(data,1);
t	= zeros(N,P);
for p=1:P
t(find(group_no==p),p)=1;
end
kernel_	= '+gauss';
PHI=sbl_kernelFunction(data,data,kernel_,kernel_width);
[probs, weights,alphas, used, kernels_]= mc_rvm(PHI,t,kernel_,3);
for p=1:P
    kernel_data{p}=data(used{p},:);
end
Decision_Surf = struct('probs',{probs},'weights', {weights},'alphas',{alphas},'used',{used},'kernels',{kernels_},'P',P,'kernel_width', {kernel_width},'kernel_data',{kernel_data});