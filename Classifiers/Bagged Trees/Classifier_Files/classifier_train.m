function [Decision_Surf]=classifier_train(data,group_no,varargin)
if nargin > 2
    nTrees = varargin{1};
else
    nTrees = 300;
end
EnsembleClass = fitensemble(data, group_no,'Bag',nTrees,'Tree','type','classification');
Decision_Surf = compact(EnsembleClass);