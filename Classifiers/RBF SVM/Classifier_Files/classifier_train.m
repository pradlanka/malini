function [Decision_Surf]=classifier_train(data, group_no, varargin)
if nargin == 2
    option =sprintf('-s 0 -t 3 -h 0 -g %d -c %d -q',1/size(data,2),1);
elseif  nargin >2
parameter = varargin{1};
option =sprintf('-s 0 -t 3 -h 0 -g %d -c %d -q', parameter(1),parameter(2)); % C value in SVM
else 
    error('Please check the optimizable parameters in your classifier');
end
Decision_Surf = svmtrain(group_no, data,option);
