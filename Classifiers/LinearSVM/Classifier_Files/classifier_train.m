function [Decision_Surf]=classifier_train(data, group_no, varargin)
if nargin == 2
    option =sprintf('-s 0 -t 0 -h 0 -q -c 1');
elseif  nargin >2
Optimizable_Paramters = length(varargin);
option =sprintf('-s 0 -t 0 -h 0 -q -c %d', varargin{1}); % C value in SVM
else 
    error('Please check the optimizable parameters in your classifier');
end
Decision_Surf = svmtrain(group_no, data,option);