function multiclass_model = trainadaboost(xtrain,ytrain,nIter)

% xtrain - The training data. RxC Rows - number of samples, Columns - the
% dimension of the samples.
% ytrain - The labels associated with the training data - Rx1 vector.
% nLabel - The number of classes.

nLabel = length(unique(ytrain));
trlabel = [];
cnt = 1;
% fp = {'b.--', 'rx-', 'k*-.', 'md--', 'go:'}; 
% figure;
for act = 1 : nLabel
    trlabel(find(ytrain ~= act)) = -1;
    trlabel(find(ytrain == act)) = 1;        
     [trerr{cnt} multiclass_model{cnt}] = call_boosting(xtrain,trlabel',nIter);        
%     plot(trerr{cnt}./size(ytrain,1),fp{cnt}), hold on
    cnt = cnt + 1;
end
% ylabel('Error');
% xlabel('Number of Boosting Iterations'), hold off;
