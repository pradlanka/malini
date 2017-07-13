function [Decision_Surf]=classifier_train(data,group_no,varargin)
if nargin > 2
    No_classifiers = varargin{1};
else
    No_classifiers = 300;
end
EnsembleClass = TreeBagger(No_classifiers, data, group_no,'FBoot',2/3,'SampleWithReplacement','On','OOBPred','On','OOBVarImp','On');
Decision_Surf= EnsembleClass;
%Decision_Surf =  compact(EnsembleClass);