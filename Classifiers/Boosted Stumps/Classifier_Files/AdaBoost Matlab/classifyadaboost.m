function [label posprob] = classifyadaboost(xtest,multiclass_model)

% xtest - the data that needs to be classified 1xC vector.
% multiclass_model - the multi class trained adaboost model
% label - The classified label.
% This function classifies the test data using the trained adaboost model
% and returns the label.
phi = 10;
wsum =zeros(size(multiclass_model,2),1);
for j = 1 : size(multiclass_model,2)    
    h_est(j) = eval_boost(multiclass_model{j},xtest);     
    for i_wl = 1 : min(250,size(multiclass_model{j},1))
        wsum(j) = wsum(j) + multiclass_model{j}{i_wl}.alpha;
    end
end    
marginals = h_est./(ones(size(h_est,1),1)*wsum');

posprob = exp(phi.*marginals)./(exp(phi.*marginals)+1);

[t label] = max(posprob);

return;