function [Decision_Surf]=classifier_train(data, group_no, varargin)
train_data = horzcat(group_no,data);
if nargin == 2
  Elm_Type =1;
  C = 1;
  kernel_type  = 'RBF_kernel';
  kernel_para = 100;
elseif  nargin >2   
  Hyperparameters =varargin{1};
  Optimizable_Paramters = length(Hyperparameters);
  Elm_Type =1;
  C = Hyperparameters(1);
  kernel_type  = 'RBF_kernel';
  kernel_para = Hyperparameters(2);
else 
    error('Please check the optimizable parameters in your classifier');
end
Decision_Surf =  elm_kernel_train(train_data, Elm_Type, C, kernel_type,kernel_para);