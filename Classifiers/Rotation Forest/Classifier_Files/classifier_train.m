
function [Decision_Surf]=classifier_train(data,group_no,varargin)
if nargin > 2
    knn = varargin{1};
else
    knn = 1;
end
Decision_Surf = rotation_forest_train(data, group_no, floor(sqrt(size(data,2))), 0.75,300,knn);