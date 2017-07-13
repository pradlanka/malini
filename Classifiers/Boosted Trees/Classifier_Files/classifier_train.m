function [Decision_Surf]=classifier_train(data,group_no,varargin)
if nargin > 2
    No_classifiers = varargin{1};
else
    No_classifiers = 300;
end
EnsembleClass = fitensemble(data, group_no,'LPBoost',No_classifiers,'Tree');
Decision_Surf = compact(EnsembleClass);