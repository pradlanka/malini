function [Decision_Surf]=classifier_train(data,group_no,varargin)
if nargin > 2
   No_neighbors = varargin{1};
else
   No_neighbors = 1;
end
Decision_Surf = ClassificationKNN.fit(data, group_no,'NumNeighbors',No_neighbors);