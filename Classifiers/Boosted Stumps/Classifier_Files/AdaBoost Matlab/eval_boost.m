
function [h,sum_alpha] = eval_boost(model,X);

h = zeros(size(X,1),1);

sum_alpha = 0;
for i=1:length(model)
  h = h + model{i}.alpha*eval_stump(model{i},X);
  sum_alpha = sum_alpha + model{i}.alpha;
end;
