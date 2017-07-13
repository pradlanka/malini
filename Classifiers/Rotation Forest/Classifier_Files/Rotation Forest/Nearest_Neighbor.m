function test_targets = Nearest_Neighbor(train_patterns, train_targets, test_patterns, Knn)

% Classify using the Nearest neighbor algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	Knn		        - Number of nearest neighbors 
%
% Outputs
%	test_targets	- Predicted targets

L			= size(train_targets,1);
Uc          = unique(train_targets);

if (L < Knn),
   error('You specified more neighbors than there are points.')
end

N               = size(test_patterns, 1);
test_targets    = zeros(N,1); 
for i = 1:N,
    dist            = sum((train_patterns -ones(L,1)* test_patterns(i,:)).^2,2);
    [m, indices]    = sort(dist);
    
    n               = hist(train_targets(indices(1:Knn),:), Uc);
    
    [m, best]       = max(n);
    
    test_targets(i,1) = Uc(best);
end
