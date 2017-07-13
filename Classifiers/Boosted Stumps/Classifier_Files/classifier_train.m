function [Decision_Surf]=classifier_train(data,group_no,varargin)
if nargin > 2
    nIter = varargin{1};
else
    nIter = 300;
end
Decision_Surf = trainadaboost(data,group_no,nIter);